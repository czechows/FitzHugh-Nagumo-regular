/* -----------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing nonrigorous Newton method 
 * to find a good guess for the periodic orbit.
 * ------------------------------------------ */



/* ------------------------------------------------------------------------------------------------------------------ */
/* ------------------------------------------ Numerical guesses class ----------------------------------------------- */
/* ------------------------------------------------------------------------------------------------------------------ */


class FhnFindPeriodicOrbit
{
public:
  int pm_count;
  DMap *vectorField;
  DMap *vectorFieldRev;
  DTaylor solver;
  DTaylor solverRev;
  std::vector<DMatrix> P_list;
  std::vector<DAffineSection> section;
  int dim;
  int fast_dim;
  std::vector<DVector> correctedGuess;
  DVector integrationTime;

  FhnFindPeriodicOrbit( std::vector<DVector> initialGuess )
    : pm_count ( initialGuess.size() ),
      vectorField( Fhn_vf ),
      vectorFieldRev( Fhn_vf_rev ),
      solver(*vectorField, order),
      solverRev(*vectorFieldRev, order),
      P_list( pm_count ),
      section( pm_count, DAffineSection( initialGuess[0], initialGuess[1]-initialGuess[0]) ),
      dim(3),
      fast_dim(dim - 1),
      correctedGuess( initialGuess ),
      integrationTime( pm_count )
  {
    for( int i = 0; i < pm_count; i++ )
    {
      (P_list[ i % pm_count ]).resize( dim, dim );
      (P_list[ i % pm_count ]).setToIdentity();

      for( int j = 1; j <= dim; j++ )
       (P_list[ i % pm_count ])( j, dim ) = ( initialGuess[ (i+1) % pm_count ]-initialGuess[ i % pm_count ] )(j);

      orthogonalizeRelativeColumn( P_list[ i % pm_count ], dim - 1 );

      section[ i % pm_count ] = DAffineSection( initialGuess[ i % pm_count ], initialGuess[ (i+1) % pm_count ]-initialGuess[ i % pm_count ] );
    }
  }

  DVector convertToVector( std::vector<DVector> list )
  {
    DVector result( fast_dim*list.size() );

    for( unsigned int i = 0; i < list.size(); i++)
    {
      for( int j = 0; j < fast_dim; j++)
       result[ i*fast_dim+j ] = (list[i])[j];
    }

    return result;
  };

  std::vector<DVector> convertToList( DVector vect )
  {
    std::vector<DVector> result( vect.dimension()/fast_dim );

    for( unsigned int i = 0; i < result.size(); i++)
    {
      for( int j = 0; j < fast_dim; j++)
       result[i][j] = vect[ i*fast_dim+j ];
    }

    return result;
  };


  std::vector<DVector> oneNewtonStep( std::vector<DVector> xList, double& error )   // one Newton step for function xi - P_{i}(x_{i-1}), i in mathbbZ
  {
    DVector x( convertToVector( xList ) );
 
    DVector pm_result(3);
    DVector setToIntegrate(3);

    DMatrix monodromyMatrix(3,3);
    DMatrix poincareDer(3,3);

    DVector x_eval( x.dimension() ); // evaluates function xi - P_{i}(x_{i-1})
    DMatrix *derMatrix;
    derMatrix = new DMatrix( x.dimension(), x.dimension() );
    (*derMatrix).setToIdentity();

    for( int i = 0; i < pm_count; i++ )
    {
      monodromyMatrix.setToIdentity(); // probably obsolete since it is done automatically (?) in the solver code
  
      setToIntegrate = (section[i]).getOrigin() + (P_list[i])*DVector({ x[2*i], x[2*i + 1], 0 });

      //    cout << "Integration to section " << i+1 << " in progress, setToIntegrate = " << setToIntegrate << " \n";
      
      DPoincareMap pm( solver, section[ (i+1) % pm_count ], poincare::MinusPlus );

      double time(0.);
      pm_result = pm(setToIntegrate, monodromyMatrix, time);

      integrationTime[i] = time;
 //     cout << time << " ";

      poincareDer = pm.computeDP(pm_result, monodromyMatrix, time);
      
      //     cout << pm_result << "\n";

      pm_result = inverseMatrix( P_list[ (i+1) % pm_count ] )*( pm_result - (section[ (i+1)%pm_count ]).getOrigin() );
      poincareDer = inverseMatrix( P_list[ (i+1) % pm_count ] )*poincareDer*(P_list[i]);

      x_eval[ (2*i + 2) % (fast_dim*pm_count) ] = x[ (2*i + 2) % (fast_dim*pm_count) ] - pm_result[0];
      x_eval[ (2*i + 3) % (fast_dim*pm_count) ] = x[ (2*i + 3) % (fast_dim*pm_count) ] - pm_result[1];

      (*derMatrix)[ (2*i+2) % (fast_dim*pm_count) ][2*i] = -poincareDer[0][0];
      (*derMatrix)[ (2*i+2) % (fast_dim*pm_count) ][2*i+1] = -poincareDer[0][1];
      (*derMatrix)[ (2*i+3) % (fast_dim*pm_count) ][2*i] = -poincareDer[1][0]; 
      (*derMatrix)[ (2*i+3) % (fast_dim*pm_count) ][2*i+1] = -poincareDer[1][1];
    }

    error = sqrt( scalarProduct( x_eval, x_eval ) );

    std::vector<DVector> result( xList.size() );
 
    for( unsigned int i = 0; i < xList.size(); i++ )
      (result[i]).resize( fast_dim );

    result = convertToList( x - capd::matrixAlgorithms::gauss( *derMatrix, x_eval ) );

    delete derMatrix;
    
    return result;
  };


  std::vector<DVector> newtonAlgorithm( double _tolerance )
  {
    double error = 1.;
  
    std::vector<DVector> x1( pm_count );
 
    for( int i = 0; i < pm_count; i++ )
    {
      (x1[i]).resize( fast_dim );
      
      for( int j = 1; j <= fast_dim; j++ )
        (x1[i])(j) = 0.;
    }


    while( error > _tolerance )
    {
      x1 = oneNewtonStep( x1, error );
      //cout << error << "\n";
    }
    
    for( int i = 0; i < pm_count; i++ )
      correctedGuess[i] = section[i].getOrigin() + P_list[i] * DVector( (x1[i])(1), (x1[i])(2), 0. );

    return x1;
  };


  std::vector<DVector> getCorrectedGuess( interval integrationTimeBound )
  {
    std::vector<DVector> tempCorrectedGuess( 2*pm_count );
    std::vector<DVector> resultCorrectedGuess;

    for( int i = 0; i < 2*pm_count; i++ )
      tempCorrectedGuess[i].resize( dim );

    for( int i = 0; i < pm_count; i++ )
    {
      if( integrationTime[ i % pm_count ] > integrationTimeBound.rightBound() )
      {
        tempCorrectedGuess[ 2*i ] = correctedGuess[i];

        DTaylor tempSolver( *Fhn_vf_ext, order );
        DCoordinateSection tempSection( dim + 1, dim, (integrationTime[ i%pm_count ])/2. );
        DPoincareMap tempPM( tempSolver, tempSection );

        DVector extVector( dim + 1 );
        for( int j = 1; j <= dim; j++ )
          extVector(j) = (correctedGuess[i])(j);

        extVector( dim + 1 ) = 0.;

        tempCorrectedGuess[ 2*i + 1 ] = DVector( dim, tempPM( extVector ).begin() );
        
        cout << i << pm_count << "OK " << correctedGuess[i] << " " << tempCorrectedGuess[ 2*i + 1 ] <<  " " << correctedGuess[i+1] << "\n";
      }
      else if( integrationTime[ (i+1) % pm_count ] < integrationTimeBound.leftBound() )
      {
        tempCorrectedGuess[ 2*i ] = DVector({ 0., 0., 0. });
        tempCorrectedGuess[ 2*i + 1 ] = DVector({ 0., 0., 0. });
        cout << "WARNING! \n";
      }
      else
      { 
        tempCorrectedGuess[ 2*i ] = correctedGuess[i];
        tempCorrectedGuess[ 2*i + 1 ] = DVector({0.,0.,0.});
      }
    }

    for( unsigned int i = 0; i < tempCorrectedGuess.size(); i++ )
    {
      if( tempCorrectedGuess[i] != DVector({ 0., 0., 0. }) )
          resultCorrectedGuess.push_back( tempCorrectedGuess[i] );
    }

    cout << pm_count << resultCorrectedGuess.size() << "\n";

    if( resultCorrectedGuess.size() < 2 )
      throw "ORBIT TOO SHORT!";

    return resultCorrectedGuess;
  }
/*
  void integrateOneTime( double _tolerance ) // TO CORRECT OR REMOVE! (SAVE COMPUTATION TIME)
  {
    std::vector<DVector> vect( returnCorrectedOrbit( _tolerance ) );
    for( int i = 0; i < pm_count-1; i++ )
    {

    DPoincareMap pm( solver, section[i+1], poincare::MinusPlus );
    double temp_time = 0.;

    DVector x = pm( vect[i], temp_time );
    cout << x  - vect[i+1]  << "\n";
    cout << temp_time << "\n";
    }
  }
*/
};



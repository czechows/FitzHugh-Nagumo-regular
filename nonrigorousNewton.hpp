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
  DMap vectorField;
  DMap vectorFieldRev;
  DTaylor solver;
  DTaylor solverRev;
  std::vector<DVector> xStartlist;
  std::vector<DMatrix> P_list;
  std::vector<DAffineSection> section;
  int dim;
  int fast_dim;

  FhnFindPeriodicOrbit( const int _pm_count )
    : pm_count ( _pm_count ),
      vectorField( Fhn_vf ),
      vectorFieldRev( Fhn_vf_rev ),
      solver(vectorField, order),
      solverRev(vectorFieldRev, order),
      xStartlist( pm_count ),
      P_list( pm_count ),
      section( pm_count, DAffineSection( xPrecomputed[0], xPrecomputed[1]-xPrecomputed[0]) ),
      dim(3),
      fast_dim(dim - 1)
  {
    for( unsigned int i = 0; i < xStartlist.size(); i++ )
    {
      (xStartlist[i]).resize( fast_dim );
      
      for( int j = 1; j <= fast_dim; j++ )
        (xStartlist[i])(j) = 0.;
    }

    for( int i = 0; i < pm_count - 1; i++ )
    {
      (P_list[i]).resize( dim, dim );
      (P_list[i]).setToIdentity();

      for( int j = 1; j <= dim; j++ )
       (P_list[i])( j, dim ) = ( xPrecomputed[i+1]-xPrecomputed[i] )(j);

      orthogonalizeRelativeColumn( P_list[i], dim - 1 );

      section[i] = DAffineSection( xPrecomputed[i], xPrecomputed[i+1]-xPrecomputed[i] );
    }
 
    (P_list[ pm_count - 1 ]).resize( dim, dim );
    (P_list[ pm_count - 1 ]).setToIdentity();

    for( int j = 1; j <= dim; j++ )
      (P_list[ pm_count - 1 ])( j, dim ) = ( xPrecomputed[0]-xPrecomputed[pm_count-1] )(j);
 
    section[ pm_count - 1 ] = DAffineSection( xPrecomputed[pm_count - 1], xPrecomputed[0]-xPrecomputed[ pm_count - 1 ] );

    orthogonalizeRelativeColumn( P_list[pm_count-1], dim - 1 );

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
    double time(0.);

    DMatrix monodromyMatrix(3,3);
    DMatrix poincareDer(3,3);

    DVector x_eval( x.dimension() ); // evaluates function xi - P_{i}(x_{i-1})
    DMatrix derMatrix( x.dimension(), x.dimension() );
    derMatrix.setToIdentity();

    for( int i = 0; i < pm_count-1; i++ )
    {
      monodromyMatrix.setToIdentity(); // probably obsolete since it is done automatically (?) in the solver code
  
      setToIntegrate = (section[i]).getOrigin() + (P_list[i])*DVector({ x[2*i], x[2*i + 1], 0 });

 //     cout << "Integration to section " << i+1 << " in progress, setToIntegrate = " << setToIntegrate << " \n";
      
      DPoincareMap pm( solver, section[i+1], poincare::MinusPlus );

      pm_result = pm(setToIntegrate, monodromyMatrix, time);
      poincareDer = pm.computeDP(pm_result, monodromyMatrix, time);
      
 //     cout << pm_result << "\n";

      pm_result = inverseMatrix( P_list[i+1] )*( pm_result - (section[i+1]).getOrigin() );
      poincareDer = inverseMatrix( P_list[i+1] )*poincareDer*(P_list[i]);

      x_eval[2*i + 2] = x[2*i + 2] - pm_result[0];
      x_eval[2*i + 3] = x[2*i + 3] - pm_result[1];

      derMatrix[2*i+2][2*i] = -poincareDer[0][0];
      derMatrix[2*i+2][2*i+1] = -poincareDer[0][1];
      derMatrix[2*i+3][2*i] = -poincareDer[1][0]; 
      derMatrix[2*i+3][2*i+1] = -poincareDer[1][1];
    }


    monodromyMatrix.setToIdentity();
  
    setToIntegrate = (section[pm_count - 1]).getOrigin() + (P_list[pm_count - 1])*DVector({ x[2*pm_count - 2], x[2*pm_count - 1], 0 });

    DPoincareMap pm( solver, section[0], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, monodromyMatrix, time );
    poincareDer = pm.computeDP( pm_result, monodromyMatrix, time);
 
    pm_result = inverseMatrix( P_list[0] )*( pm_result - (section[0]).getOrigin() );
    poincareDer = inverseMatrix( P_list[0] )*poincareDer*(P_list[pm_count - 1]);

    x_eval[0] = x[0] - pm_result[0];
    x_eval[1] = x[1] - pm_result[1];

    derMatrix[0][2*pm_count - 2] = -poincareDer[0][0];
    derMatrix[0][2*pm_count - 1] = -poincareDer[0][1];
    derMatrix[1][2*pm_count - 2] = -poincareDer[1][0];
    derMatrix[1][2*pm_count - 1] = -poincareDer[1][1];

    error = sqrt( scalarProduct( x_eval, x_eval ) );

    std::vector<DVector> result( xList.size() );
 
    for( unsigned int i = 0; i < xList.size(); i++ )
      (result[i]).resize( fast_dim );

 //   cout << "AFTER INTEGRATION, BEFORE GAUSS \n";
    result = convertToList( x - capd::matrixAlgorithms::gauss( derMatrix, x_eval ) );
 //   cout << "AFTER GAUSS \n";
    
    return result;
  };


  std::vector<DVector> newtonAlgorithm( double _tolerance )
  {
    double error = 1.;
  
    DVector x1( convertToVector( xStartlist ) );

    while( error > _tolerance )
    {
      x1 = convertToVector( oneNewtonStep( convertToList(x1), error ) );
      //cout << error << "\n";
    }
    
    return convertToList( x1 );
  };


  std::vector<DVector> returnCorrectedOrbit( double _tolerance )
  {
    std::vector<DVector> xList = newtonAlgorithm( _tolerance );

    std::vector<DVector> result;

    for( int i = 0; i < pm_count; i++ )
      result.push_back( (section[i]).getOrigin() + (P_list[i]) * DVector({ (xList[i])[0], (xList[i])[1], 0. }) );

    return result;
  }

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

};



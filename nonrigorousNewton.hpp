/* -----------------------------------------
 * This is a header file for continuation_fhn.cpp
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
  std::vector<DMatrix> P_list;
  std::vector<DAffineSection> section;
  int dim;
  int fast_dim;
  std::vector<DVector> correctedGuess;
  DVector integrationTime;
  DEuclNorm vectorNorm;

  FhnFindPeriodicOrbit( std::vector<DVector> initialGuess )
    : pm_count ( initialGuess.size() ),
      vectorField( Fhn_vf ),
      vectorFieldRev( Fhn_vf_rev ),
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

 //#include "nonrigorousNewtonBetterCoordinates.hpp"  // diagonalizing coordinates for nonrigorous Newton helps a bit but is not worth the effort

  DVector convertToVector( std::vector<DVector> list )
  {
    DVector result( (dim-1)*list.size() );

    for( unsigned int i = 0; i < list.size(); i++)
    {
      for( int j = 0; j < dim-1; j++)
       result[ i*(dim-1)+j ] = (list[i])[j];
    }

    return result;
  };

  std::vector<DVector> convertToList( DVector vect )
  {
    std::vector<DVector> result( vect.dimension()/(dim-1) );

    for( unsigned int i = 0; i < result.size(); i++)
    {
      for( int j = 0; j < dim-1; j++)
       result[i][j] = vect[ i*(dim-1)+j ];
    }

    return result;
  };


  std::vector<DVector> oneNewtonStep( std::vector<DVector> xList, double& error )   // one Newton step for function xi - P_{i}(x_{i-1})
  {
    
    DVector *x;
    x = new DVector( convertToVector( xList ) );
 
    DVector pm_result(3);
    DVector setToIntegrate(3);

    DMatrix monodromyMatrix(3,3);
    DMatrix poincareDer(3,3);

    DVector *x_eval;
    x_eval = new DVector( (*x).dimension() ); // evaluates function xi - P_{i}(x_{i-1})
    DMatrix *derMatrix;
    derMatrix = new DMatrix( (*x).dimension(), (*x).dimension() );
    (*derMatrix).setToIdentity();



    for( int i = 0; i < pm_count; i++ )
    {
      monodromyMatrix.setToIdentity(); // probably obsolete since it is done automatically in the solver code
  
      setToIntegrate = (section[i]).getOrigin() + (P_list[i])*DVector({ (*x)[2*i], (*x)[2*i + 1], 0. }); // actually point to integrate, we use interval integration nomenclature here


      int local_order( order );
      bool isDomainError(0);
 /*     double maxStep(4.);

    while(true)
      {
        try{*/
      {
        DTaylor solver(*vectorField, local_order);
       // solver.setMaxStep(maxStep);
        DPoincareMap pm( solver, section[ (i+1) % pm_count ], poincare::MinusPlus );

        double time(0.);

        pm_result = pm(setToIntegrate, monodromyMatrix, time);

        if( !isDomainError )
          integrationTime[i] = time;

        poincareDer = pm.computeDP(pm_result, monodromyMatrix, time);
     /* }
          break;
        }
        catch(std::domain_error)  // this catch would decrease the order or maxStep so nonrigorous integration would not cross a section in one step 
                                    // (then integrator throws domain_error); for the parameter range in the paper not necessary and did not work that well
                                    // so we disable it, however we leave it to uncomment for future purposes
        {
          cout << "DOMAIN ERROR! \n";
         // local_order = local_order - 1;
          maxStep = maxStep/2.;
          integrationTime[i] = 0.;
          isDomainError = 1;
          cout << "MAX STEP DECREASED TO " << maxStep << "! \n";
         // cout << "ORDER DECREASED TO " << local_order << " ! \n"; *?
        } */
      }

      pm_result = inverseMatrix( P_list[ (i+1) % pm_count ] )*( pm_result - (section[ (i+1)%pm_count ]).getOrigin() );
      poincareDer = inverseMatrix( P_list[ (i+1) % pm_count ] )*poincareDer*(P_list[i]);

      (*x_eval)[ (2*i + 2) % (*x_eval).dimension() ] = (*x)[ (2*i + 2) % (*x_eval).dimension() ] - pm_result[0];
      (*x_eval)[ (2*i + 3) % (*x_eval).dimension() ] = (*x)[ (2*i + 3) % (*x_eval).dimension() ] - pm_result[1];

      (*derMatrix)[ (2*i+2) % (*derMatrix).numberOfRows() ][2*i] = -poincareDer[0][0];
      (*derMatrix)[ (2*i+2) % (*derMatrix).numberOfRows() ][2*i+1] = -poincareDer[0][1];
      (*derMatrix)[ (2*i+3) % (*derMatrix).numberOfRows() ][2*i] = -poincareDer[1][0]; 
      (*derMatrix)[ (2*i+3) % (*derMatrix).numberOfRows() ][2*i+1] = -poincareDer[1][1];
    }

    error = vectorNorm( *x_eval );

    std::vector<DVector> result( xList.size() );
 
    for( unsigned int i = 0; i < xList.size(); i++ )
      (result[i]).resize( dim-1 );

    cout << "Newton algorithm in progress...";
    result = convertToList( *x - capd::matrixAlgorithms::gauss( *derMatrix, *x_eval ) );
    cout << " done! ";

    delete derMatrix;
    delete x_eval;
    delete x;
    
    return result;
  };


  std::vector<DVector> newtonAlgorithm( double _tolerance )
  {
    //initD( _tolerance );  // these three lines apply setting diagonalized coordinates for nonrigorous Newton - commented because they don't help much
    //for( int i = 0; i < pm_count; i++ )
    // orthogonalizeRelativeColumn( P_list[ i % pm_count ], dim - 1 );

    double error = 1.;
  
    std::vector<DVector> x1( pm_count );
 
    for( int i = 0; i < pm_count; i++ )
    {
      (x1[i]).resize( dim-1 );
      
      for( int j = 1; j <= dim-1; j++ )
        (x1[i])(j) = 0.;
    }


    while( error > _tolerance )
    {
      x1 = oneNewtonStep( x1, error );
      cout << "Error size: " << error << "\n";

      if( error > 2. )
        throw "NEWTON ALGORITHM DIVERGENT! \n";
    }

    
    for( int i = 0; i < pm_count; i++ )
      correctedGuess[i] = section[i].getOrigin() + P_list[i] * DVector( (x1[i])(1), (x1[i])(2), 0. );

    return x1;
  };


  void saveIntegrationTime() // we save approximate integration times for future reference
  {
    std::ofstream savedIntegrationTime;
    savedIntegrationTime.open("savedIntegrationTime.log", std::ios::trunc);
    savedIntegrationTime.precision(6);

    double totalIntegrationTime(0.);

    for( unsigned int i = 0; i < integrationTime.dimension(); i++ )
    {
     savedIntegrationTime << "Approximate integration time till section given by: " << section[i].getOrigin() << " is: " << integrationTime[i] << "\n";
     totalIntegrationTime = totalIntegrationTime + integrationTime[i];
    }

    savedIntegrationTime << "Approximate total integration time: " << totalIntegrationTime;
  }

  void saveOrbit( std::vector<DVector> orbit )     // we store the orbit so one can resume the program later from the saved numerical guess 
                                                   // (we do not have to store the exact copy of the guess!"
  {
    std::ofstream savedOrbit;
    savedOrbit.open("savedOrbit.hpp", std::ios::trunc);
    savedOrbit.precision(10);
    savedOrbit << "std::vector<DVector> savedOrbit = \n { \n";

    for( unsigned int i = 0; i < orbit.size() - 1; i++ )
      savedOrbit << "DVector(" << orbit[i] << "), \n";

    savedOrbit << "DVector(" << orbit[orbit.size() - 1] << ") \n };";
  }


  std::vector<DVector> getCorrectedGuess( interval integrationTimeBound )
  {
    std::vector<DVector> resultCorrectedGuess;


    for( int i = 0; i < pm_count; i++ )
    {
      if( integrationTime[ i ] > integrationTimeBound.rightBound() )
      {
        resultCorrectedGuess.push_back( correctedGuess[i] );

        DTaylor tempSolver( *vectorField, order );
        DTimeMap timeMap( tempSolver );
        
        DVector pointToIntegrate( correctedGuess[i] ); // for safety as timeMap changes the argument too
        resultCorrectedGuess.push_back( timeMap( (integrationTime[ i ])/2., pointToIntegrate ) );

        cout << "Long integration time. Section added. \n";
      }
      else if( integrationTime[ (i+1)%pm_count ] < integrationTimeBound.leftBound() )
        cout << "Short integration time. Section removed. \n";
      else
        resultCorrectedGuess.push_back( correctedGuess[i] );
    }

    cout << "Number of sections changed from: " << pm_count << " to: " << resultCorrectedGuess.size() << ". \n";

    if( resultCorrectedGuess.size() < 2 )
      throw "ORBIT TOO SHORT!";

    saveOrbit( resultCorrectedGuess );
    saveIntegrationTime();

    return resultCorrectedGuess;
  }

  void setDParameters( interval theta, interval eps )
  {
    double thetaD = ( theta.leftBound() + theta.rightBound() )/2. ;
    double epsD = ( eps.leftBound() + eps.rightBound() )/2. ;

    (*vectorField).setParameter( "theta", thetaD );
    (*vectorField).setParameter( "eps", epsD );
 
    (*vectorFieldRev).setParameter( "theta", thetaD );
    (*vectorFieldRev).setParameter( "eps", epsD );
  }

};



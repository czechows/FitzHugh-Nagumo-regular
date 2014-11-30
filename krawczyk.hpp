//


class FhnKrawczyk : public FhnCovering
{
  public:
  FhnKrawczyk( std::vector<DVector> initialGuess )
    : FhnCovering( initialGuess )
  {
  }

  IVector convertToVector( std::vector<IVector> list )
  {
    IVector result( (dim-1)*list.size() );

    for( unsigned int i = 0; i < list.size(); i++)
    {
      for( int j = 0; j < dim-1; j++)
       result[ i*(dim-1)+j ] = (list[i])[j];
    }

    return result;
  };

  std::vector<IVector> convertToList( IVector vect )
  {
    std::vector<IVector> result( vect.dimension()/(dim-1) );

    for( unsigned int i = 0; i < result.size(); i++)
    {
      for( int j = 0; j < dim-1; j++)
       result[i][j] = vect[ i*(dim-1)+j ];
    }

    return result;
  };
 
  IMatrix fiDer_withParams( int i, IVector _x0, interval eps, interval &iterationTime, int disc_j=5, int disc_k=5, int disc_l=1 )
  {
    iterationTime = 0.;
    IVector pm_resultExt( dim + 1 );
    ITaylor solverExt( *IVectorFieldExt, rig_order );

    IMatrix monodromyMatrix(4,4);
    IMatrix poincareDer(4,4);

    for( int j = 1; j <= disc_j; j++ )
    {
      for( int k = 1; k <= disc_k; k++ )
      {
        for( int l = 1; l <= disc_l; l++ )
        {
          interval tj = interval(j-1,j)/disc_j;
          interval tk = interval(k-1,k)/disc_k;
          interval tl = interval(l-1,l)/disc_l;
 
          interval epsDiam = right(eps) - left(eps);
          interval epsRange( (-epsDiam/2.).leftBound(), (epsDiam/2.).rightBound() );

          interval x0_1j = _x0(1).leftBound() + tj * ( _x0(1).rightBound() - _x0(1).leftBound() );
          interval x0_2k = _x0(2).leftBound() + tk * ( _x0(2).rightBound() - _x0(2).leftBound() );
          interval epsRange_l = epsRange.leftBound() + tl * ( epsRange.rightBound() - epsRange.leftBound() );


          C1Rect2Set setToIntegrate( extendAffineSection( ISection[ i % pm_count ], eps ).getOrigin(), 
              extendP( IP_list[ i % pm_count ] ), IVector({ x0_1j, x0_2k, 0., epsRange_l }) );
          IAffineSection sectionExt = extendAffineSection( ISection[ (i+1)% pm_count ], eps );

          IPoincareMap pm( solverExt, sectionExt, poincare::MinusPlus );
 

          interval time(0.);
          monodromyMatrix.setToIdentity();

          pm_resultExt = pm( setToIntegrate, monodromyMatrix, time );

          if( j + k + l == 3 )
          {
            poincareDer = pm.computeDP( pm_resultExt, monodromyMatrix, time );
            poincareDer = inverseMatrix( extendP( IP_list[ (i+1)%pm_count] ) )*poincareDer*( extendP( IP_list[i] ) );          
            iterationTime = time;
          }
          else
          {
            poincareDer = intervalHull(  inverseMatrix( extendP( IP_list[(i+1)%pm_count] ) )*(pm.computeDP( pm_resultExt, monodromyMatrix, time ))*( extendP( IP_list[i] ) ), poincareDer );
            iterationTime = intervalHull( iterationTime, time );
          }
        }
      }  
    }

    IMatrix fiDer(2,2);

    for( int j = 1; j <= 2; j++ )
    {
      for( int k = 1; k <= 2; k++ )
      {
        fiDer(j,k)=poincareDer(j,k);
      }
    }
    return fiDer;
  }


  IMatrix computeDerMatrix( IVector x, interval eps )   
  {
    std::vector<IVector> x_list = convertToList( x );
 
    IMatrix derMatrix( x.dimension(), x.dimension() );
    derMatrix.setToIdentity();

    for( int i = 0; i < pm_count-1; i++ )
    {
      interval iterationTime(0.);
      IMatrix pDer = fiDer_withParams( i, x_list[i], eps, iterationTime, 5, 5, 1 );
      totalPeriod = totalPeriod + iterationTime;

      derMatrix[2*i+2][2*i] = -pDer[0][0];
      derMatrix[2*i+2][2*i+1] = -pDer[0][1];
      derMatrix[2*i+3][2*i] = -pDer[1][0]; 
      derMatrix[2*i+3][2*i+1] = -pDer[1][1];
    }
 
    interval iterationTime(0.);
    IMatrix pDer = fiDer_withParams( pm_count-1, x_list[pm_count-1], eps, iterationTime, 5, 5, 1 );
    totalPeriod = totalPeriod + iterationTime;

    derMatrix[0][2*pm_count - 2] = -pDer[0][0];
    derMatrix[0][2*pm_count - 1] = -pDer[0][1];
    derMatrix[1][2*pm_count - 2] = -pDer[1][0];
    derMatrix[1][2*pm_count - 1] = -pDer[1][1];

    return derMatrix;
  }

  IMatrix approxInverse( IMatrix derMatrix )
  {
    IMatrix tempDerMatrix( midVector( derMatrix ) );

    DMatrix DderMatrix( pm_count*(dim-1), pm_count*(dim-1) );

    for( int i=0; i< pm_count*(dim-1); i++ ) 
    {
      for( int j=0; j < pm_count*(dim-1); j++ )
        DderMatrix[i][j] = tempDerMatrix[i][j].leftBound();
    }

    DMatrix DapproxInverse = capd::matrixAlgorithms::gaussInverseMatrix( DderMatrix ); 

    return IMatrix( DapproxInverse );
  };

  IVector computeF( IVector x, interval eps )   // evaluates function xi - P_{i}(x_{i-1})
  { 
    std::vector<IVector> x_list = convertToList( x );

    IVector x_eval( x.dimension() ); 
    for( int i = 0; i < pm_count-1; i++ )
    {
      interval tempIterationTime(0.);
      IVector pm_result = fi_withParams( i, x_list[i], eps, tempIterationTime, 1, 1, 1 ); 
      
      x_eval[ (2*i + 2) % (2*pm_count) ] = x[ (2*i + 2) % (2*pm_count) ] - pm_result[0];
      x_eval[ (2*i + 3) % (2*pm_count) ] = x[ (2*i + 3) % (2*pm_count) ] - pm_result[1]; 
    }

    return x_eval;
  };



  void proveExistenceOfOrbitWithKrawczyk( interval theta, interval eps, double _tolerance, double _radius )
  {
    setDParameters( theta, eps );

    (*IVectorField).setParameter( "theta", theta );
    (*IVectorField).setParameter( "eps", eps );

    (*IVectorFieldExt).setParameter( "theta", theta );

    init( _tolerance, _radius );

    IVector X_vector = convertToVector( X );
    IVector x0_vector( midVector( X_vector ) );

    totalPeriod = 0.;
    IMatrix DF( computeDerMatrix( X_vector, eps ) ); 

    IMaxNorm vectorNorm;

    IVector Fx0 = computeF( x0_vector, eps );

    IVector N( 2*pm_count );
    IVector K( 2*pm_count );

    try{
      IMatrix DFinverse = capd::matrixAlgorithms::gaussInverseMatrix( DF );
      N = -DFinverse * Fx0 + x0_vector;

      if( subsetInterior( N, X_vector ) )
      {
        cout << "\nInterval Newton method succeeds! The periodic orbit is locally unique! \n";
        cout << "norm(N) = " << vectorNorm(N) << " norm(X) = " << vectorNorm(X_vector) << " \n \n"; 
      }
      else
      {
        cout << "\nExistence of periodic orbit using interval Newton not proven\n"; 
        cout << "norm(N) = " << vectorNorm(N) << " norm(X) = " << vectorNorm(X_vector) << " \n"; 
        cout << "norm(invDFXFx0) = " << vectorNorm( DFinverse * Fx0 ) << " \n \n"; 
        throw std::runtime_error("INTERVAL NEWTON METHOD FAILURE");
      }
    }
    catch( std::runtime_error& e )  
    {
      cout << "\n" << e.what() << "\n";

      cout << "Switching to interval Krawczyk..\n";

      IMatrix approxDFInverse = approxInverse( DF );
      IMatrix Id( 2*pm_count, 2*pm_count );
      Id.setToIdentity();

      K = x0_vector- approxDFInverse*Fx0 + (Id - approxDFInverse*DF)*( X_vector - x0_vector ); 

      if( subsetInterior( K, X_vector ) )
      {
          cout << "Interval Krawczyk method succeeds! The periodic orbit is locally unique! \n";      
          cout << "norm(K) = " << vectorNorm(K) << " norm(X) = " << vectorNorm(X_vector) << " \n \n"; 
          wasKrawczykNeeded = 1;
      }
      else
      {
        cout << "Existence of periodic orbit using interval Krawczyk not proven\n"; 
        cout << "norm(K) = " << vectorNorm(K) << " norm(X) = " << vectorNorm(X_vector) << " \n"; 
        cout << "norm(Fx0) = " << vectorNorm( Fx0 ) << "\n";
        cout << "norm(Id-CdFX)" << vectorNorm(Id - approxDFInverse*DF) << "\n \n";

        throw "INTERVAL KRAWCZYK METHOD FAILURE \n";
      }
    }


  }

};




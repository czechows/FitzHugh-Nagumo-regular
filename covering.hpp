/* -------------------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing covering relations 
 * for finding periodic solutions of the FHN equations
 * ------------------------------------------------- */



class FhnCovering : public FhnFindPeriodicOrbit
{
  public:
  IMap *IVectorField;
  ITaylor ISolver;
  std::vector<IMatrix> IP_list;
  std::vector<IAffineSection> ISection;
  std::vector<IVector> X;
  DEuclNorm vectorNorm;
 
  FhnCovering( std::vector<DVector> initialGuess )
    : FhnFindPeriodicOrbit( initialGuess ),
      IVectorField( IFhn_vf ),
      ISolver( *IVectorField, rig_order),
      IP_list( pm_count ),
      ISection( pm_count, IAffineSection( IVector(dim), IVector(dim) ) ),
      X( pm_count )
  {
    for( int i = 0; i < pm_count; i++ )
    {
      (IP_list[ i % pm_count ]).resize( dim, dim );
      (IP_list[i % pm_count ]).setToIdentity();
      X[ i % pm_count ].resize( fast_dim );

      IVector IsectionCenter( dim );
      IVector IsectionNormal( dim );

      for( int j = 1; j <= dim; j++ )
      {
        IsectionCenter(j) = ( initialGuess[ i % pm_count ] )(j);
        IsectionNormal(j) = ( initialGuess[ (i+1) % pm_count ] - initialGuess[ i % pm_count ] )(j);
      }

      ISection[ i % pm_count ] = IAffineSection( IsectionCenter, IsectionNormal ); 
    }
  }


  DMatrix setCoordSystem( int coord_no, DVector _disp )            // similar to coordChange in fhn program, i < pm_count - 1 !       
  {
    DVector Gamma = section[ coord_no ].getOrigin() + P_list[ coord_no ]*DVector( _disp(1), _disp(2), 0. );

    int vdim( 3 );   // should be used only in dimension 3!
    DMatrix JacobianD( vdim, vdim );

    // a patch to set eps to 0 for the vector field for computing coordinates around slow manifold
    DMap vectorFieldZeroEps( *vectorField );          
    vectorFieldZeroEps.setParameter("eps", 0.);
    // WARNING: specific to the (type of) the vector field. Takes care of the problem that proof does not go for subintervals of parameter epsilon away from 0.
   
    for(int i=0; i<vdim; i++)              
    {
      for(int j=0; j<vdim; j++)
        JacobianD[i][j] = ( vectorFieldZeroEps[Gamma] )[i][j];
    }

    DVector tempvectRe( vdim );   
    DVector tempvectIm( vdim );

    DMatrix tempmatrixIm( vdim, vdim );
    
    DMatrix P( vdim, vdim );
    DMatrix P_result( vdim, vdim );

    computeEigenvaluesAndEigenvectors(JacobianD, tempvectRe, tempvectIm, P, tempmatrixIm);

    int i_min(42);  // i's initialized with any value
    int i_med(42);
    int i_max(42); 

    for( int i = 1; i <= vdim; i++ ) // sorting so we are sure we have stable coordinate first unstable second neutral third. ONLY FOR 3D vector fields!
    {
      if( tempvectRe(i) == min( tempvectRe(1), min(tempvectRe(2), tempvectRe(3)) ) )
        i_min = i;
      else if( tempvectRe(i) == max( tempvectRe(1), max(tempvectRe(2), tempvectRe(3)) ) )
        i_max = i;
      else 
        i_med = i;
    }

    for( int i = 1; i <=vdim; i++ )
    { 
      P_result(i,1) = P( i, i_min );
      P_result(i,2) = P( i, i_max );
      P_result(i,3) = P( i, i_med );
    }

    return P_result; 
  };



  void init( double _tolerance, double _radius )
  {
    std::vector<DVector> x0_double( newtonAlgorithm( _tolerance ) );

    for( int i = 0; i < pm_count; i++ )
      X[i] = IVector( { interval( -_radius, _radius ), interval( -_radius, _radius ) } );

    int breakPoint( pm_count + 1 );
    DMatrix P_breakPoint( dim, dim );

    for( int i = 0; i < pm_count; i++ )
    {
      if( ( (section[i].getOrigin() )(dim) > 0.05 ) && ( (section[i].getOrigin() )(dim)  < 0.06 ) ) 
      {
        P_breakPoint = setCoordSystem( i, x0_double[i] );
        breakPoint = i;
        break;
      }
    }

    if( breakPoint == pm_count + 1 )
      throw "COULD NOT INITIALIZE COORDINATE SYSTEMS AROUND APPROXIMATE ORBIT. NOT ENOUGH SUBSECTIONS";
  
    DVector unstableVect( P_breakPoint.column(1) );
    DVector stableVect( P_breakPoint.column(0) );

    double error(1.);

    while( error > _tolerance )
    {
      DVector unstableVectCp( unstableVect ),
              stableVectCp( stableVect );

      error = vectorNorm( unstableVectCp - unstableVect ) + vectorNorm( stableVectCp - stableVect );

      for( int i = 0; i < pm_count; i++ )
      {
        DMatrix uMonodromyMatrix( dim, dim ),
                sMonodromyMatrix( dim, dim );

        double uReturnTime(0.),
               sReturnTime(0.);

        int uSecCount( (breakPoint + i) % pm_count ),
            sSecCount( ( (breakPoint - i) + pm_count ) % pm_count );

        DVector uSetToIntegrate( correctedGuess[ uSecCount ] ),
                sSetToIntegrate( correctedGuess[ sSecCount ] );

        DPoincareMap uPM( solver, section[ (uSecCount + 1) % pm_count ] ), 
                     sPM( solverRev, section[ ((sSecCount - 1) + pm_count) % pm_count ] );

        DVector uResult( uPM( uSetToIntegrate, uMonodromyMatrix, uReturnTime ) ),
                sResult( sPM( sSetToIntegrate, sMonodromyMatrix, sReturnTime ) );

        DMatrix uDerivativeOfPm( uPM.computeDP( uResult, uMonodromyMatrix, uReturnTime ) ),
                sDerivativeOfPm( sPM.computeDP( sResult, sMonodromyMatrix, sReturnTime ) );

        unstableVect = uDerivativeOfPm * unstableVect;
        stableVect = sDerivativeOfPm * stableVect;

        unstableVect = unstableVect / vectorNorm( unstableVect );
        stableVect = stableVect / vectorNorm( stableVect );
        
        for( int j = 1; j <= dim; j++ )
        {
          (IP_list[ uSecCount ] )( j, 2 ) = unstableVect( j );
          (IP_list[ sSecCount ])( j, 1 ) = stableVect( j );
        }
      }

    }
 
    for( int i = 0; i < pm_count; i++)
    {
      DVector origin = DVector( section[i].getOrigin() + P_list[i] * DVector( (x0_double[i])(1), (x0_double[i])(2), 0. ) );
      (ISection[i]).setOrigin( IVector( origin ) );
    }

    for( int i = 0; i < pm_count; i++)
    {
      (ISection[i]).setNormalVector( ISection[ (i+1) % pm_count ].getOrigin() - (ISection[ i % pm_count ]).getOrigin() );

      for( int j = 1; j <= dim; j++ )
       ( IP_list[i] )( j, dim ) =( (ISection[i]).getNormalVector() )(j);
  
      IOrthogonalizeRelativeColumn( IP_list[ i % pm_count ], dim - 1 );
    }
  }

  IVector fi( int i, IVector _x0 )
  {
    IVector pm_result( dim );
    interval time(0.);

    IVector fi_eval( fast_dim ); 
   
    C0Rect2Set setToIntegrate( (ISection[ i % pm_count ]).getOrigin(), IP_list[ i % pm_count ], IVector({ _x0(1), _x0(2), 0. }) );

    IPoincareMap pm( ISolver, ISection[ (i+1) % pm_count ], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, (ISection[ (i+1) % pm_count ]).getOrigin(), inverseMatrix(IP_list[ (i+1) % pm_count ]), time ); 

    for( int j = 0; j < fast_dim; j++ )
      fi_eval[j] = pm_result[j];

    return fi_eval;
  }

  IVector leftU(const IVector &N)
  {
    IVector _leftU( N.dimension() ); 
    _leftU = N;
    _leftU[1] = N[1].leftBound();
    return _leftU;
  }

  IVector rightU(const IVector &N)
  {
    IVector _rightU( N.dimension() );
    _rightU = N;
    _rightU[1] =  N[1].rightBound();
    return _rightU;
  }
 
  IVector shrinkAndExpand(const IVector &N, interval factor )  // shrinks a rectangle in unstable direction and expands it in stable to get a covering (for example by original rectangle)
  {
      IVector result(N);
      result[0] = N[0]*factor;
      result[1] = N[1]/factor;
      return result;
  }

  IVector propagateToCoverByfi( int i, IVector setCovering, IVector unstableDirLimit ) // unstableDirLimit is the limit of how much the set to cover can grow in the unstable direction
  {
    const interval EPS = interval(1./1e15);
 
    IVector image = fi( i, setCovering );
    IVector unstablePart = IVector( { interval( fi( i, leftU( setCovering ) )[1].rightBound(), fi( i, rightU( setCovering ) )[1].leftBound() ) } );

    if( !( vectalg::containsZero( image ) && unstablePart[0].leftBound() < 0. && unstablePart[0].rightBound() > 0. ) )// these assumptions are not necessary but useful to check
      throw "COVERING ERROR! \n";

    unstablePart = intersection( unstableDirLimit, unstablePart );

    IVector setToCover = IVector( { image[0], unstablePart[0] } );

    return shrinkAndExpand( setToCover, interval( 1. + EPS ) );
  }


  bool isCovering( const IVector& setCovering, const IMatrix& setCoveringCoord, const IVector& setToCover ) 
                                                        // verifies covering between image of setCovering by a matrix setCoveringCoord over setToCover 
                                                        // first variable stable second unstable
  {
   bool leftCheck( (setCoveringCoord*( leftU(setCovering) ))[1] < setToCover[1].leftBound() );
   bool rightCheck( (setCoveringCoord*( rightU(setCovering) ))[1] > setToCover[1].rightBound() );

   if( leftCheck && rightCheck && subsetInterior( (setCoveringCoord*setCovering)[0], setToCover[0] ) ) 
    return 1;
   else
     return 0;
  };


  void proveExistenceOfOrbit( interval theta, interval eps, double _tolerance, double _radius )
  {
    double thetaD = ( theta.leftBound() + theta.rightBound() )/2. ;
    double epsD = ( eps.leftBound() + eps.rightBound() )/2. ;

    (*vectorField).setParameter( "theta", thetaD );
    (*vectorField).setParameter( "eps", epsD );
 
    (*vectorFieldRev).setParameter( "theta", thetaD );
    (*vectorFieldRev).setParameter( "eps", epsD );

    (*IVectorField).setParameter( "theta", theta );
    (*IVectorField).setParameter( "eps", eps );

    init( _tolerance, _radius );

    IVector vectorToCover( X[0] );
    IVector vectorCovering( X[0] );
    IVector unstableLimit( { 10*vectorCovering[1] } );

    for( int i = 0; i < pm_count; i++ )
    {
      //    cout << "test " << (( ISection[i].getOrigin() )(3)).leftBound() << "\n";
      vectorCovering = propagateToCoverByfi( i, vectorCovering, unstableLimit );
      //  cout << isCoveringByfi( i, X[i], X[i+1] ) << "\n"; 
      //  cout << fi( i, X[i] ) << "\n " << fi( i, leftU(X[i]) ) << " " << fi( i, rightU(X[i]) ) << "\n \n";
    }


    if( !isCovering( vectorCovering, IMatrix::Identity(2), vectorToCover ) )
      throw "COVERING ERROR! \n";
  }


};

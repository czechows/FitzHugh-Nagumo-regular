/* -----------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing Piotr Zgliczy{\'n}ski ideas
 * on how to compute the inverse derivative matrix
 * which occurs in problems of applying the interval Newton operator
 * to composition of Poincar{\'e} maps.
 * Methods used here are elements 
 * of the FhnIntervalNewton class.
 * ------------------------------------------ */


IMatrix setCoordSystem( int i )            // similar to coordChange in fhn program, i < pm_count - 1 !       
{
  IVector Gamma = ISection[i].getOrigin() + ( IMatrix(P_list[i]) )*IVector( (x0[i])(1), (x0[i])(2), 0. );

  int vdim( 3 );   // should be used only in dimension 3!
  DMatrix JacobianD( vdim, vdim );

  // a patch to set eps to 0 for the vector field for computing coordinates around slow manifold
  IMap vectorFieldZeroEps( IFhn_vf );          
  vectorFieldZeroEps.setParameter("eps", interval(0.) );
  // WARNING: specific to the (type of) the vector field. Takes care of the problem that proof does not go for subintervals of parameter epsilon away from 0.
 
  for(int i=0; i<vdim; i++)              
  {
    for(int j=0; j<vdim; j++)
      JacobianD[i][j] = ( ( vectorFieldZeroEps[Gamma] )[i][j] ).leftBound();
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

  IMatrix IP_result( P_result );

  for( int j = 1; j <= dim; j++ )
   IP_result( j, dim ) = ( IVector( IxPrecomputed[ (i+1) % pm_count]-IxPrecomputed[ i%pm_count ] ) )(j);

  IOrthogonalizeRelativeColumn( IP_result, dim - 1 );

  return IP_result; 
};



void Xx0init( double _tolerance, double _radius )
{
  std::vector<DVector> x0_double( newtonAlgorithm( _tolerance ) );

  for( int i = 0; i < pm_count; i++ )
  {
    for( int j = 1; j <= fast_dim; j++ )
      (x0[i])(j) = (x0_double[i])(j);
  }

  for( int i = 0; i < pm_count; i++ )
    X[i] = (x0[i])*interval( 1. - _radius, 1. + _radius );
 
  int breakPoint(7);

  for( int i = 0; i < pm_count; i++ )
  {

    if( ( ((ISection[i].getOrigin() )(dim)).leftBound() > 0.05 ) && ( ((ISection[i].getOrigin() )(dim)).rightBound()  < 0.06 ) ) 
    {
      IP_list[i] = setCoordSystem( i );
      breakPoint = i;
      break;
    }
  }
  
  for( int i = breakPoint; i < breakPoint + pm_count; i++ )
  {
    DMatrix monodromyMatrix( dim, dim );
    double returnTime(0.);
    DVector C1TempCenterSet( section[ i % pm_count ].getOrigin() + P_list[ i % pm_count ] * DVector( (x0_double[i % pm_count])(1), (x0_double[i % pm_count])(2), 0. ) );
    DPoincareMap tempPM( solver, section[ (i+1) % pm_count ] ); 

    DVector tempVect = tempPM( C1TempCenterSet, monodromyMatrix, returnTime );
    IEuclNorm matrixNorm;
    DMatrix derivativeOfPm = tempPM.computeDP( tempVect, monodromyMatrix, returnTime );
    
    IMatrix IDerivativeOfPm = IMatrix( derivativeOfPm );

    IP_list[ (i+1) % pm_count] = midVector( IDerivativeOfPm*( IP_list[ i % pm_count ] ) );
    IP_list[ (i+1) % pm_count] = IP_list[ (i+1) % pm_count] /( matrixNorm( IP_list[ (i+1) % pm_count] ) );
        
    for( int j = 1; j <= dim; j++ )
     (IP_list[ (i+1) % pm_count])( j, dim ) = ( IVector( IxPrecomputed[ (i+1) % pm_count]-IxPrecomputed[ i%pm_count ] ) )(j);

    IOrthogonalizeRelativeColumn( IP_list[ (i+1)%pm_count ], dim - 1 );

    cout << IP_list[ (i+1)%pm_count ] << "\n";
  }
}


IVector fi( int i, IVector _x0 )
{
  IVector pm_result( dim );
  interval time(0.);

  IVector fi_eval( fast_dim ); 
 
  if( i == pm_count - 1 )
  {
    C0Rect2Set setToIntegrate( (ISection[pm_count-1]).getOrigin(), IP_list[pm_count-1], IVector({ _x0(1), _x0(2), 0. }) );

    IPoincareMap pm( ISolver, ISection[0], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, (ISection[0]).getOrigin(), inverseMatrix(IP_list[0]), time );
  }
  else if( i < pm_count - 1 )
  {
    C0Rect2Set setToIntegrate( (ISection[i]).getOrigin(), IP_list[i], IVector({ _x0(1), _x0(2), 0. }) );

    IPoincareMap pm( ISolver, ISection[i+1], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, (ISection[i+1]).getOrigin(), inverseMatrix(IP_list[i+1]), time ); 
  }
  else
    throw "MAP INDEX OUT OF RANGE! \n";

  for( int i = 0; i < fast_dim; i++ )
    fi_eval[i] = pm_result[i];

  return fi_eval;
}





void proveExistenceOfOrbit( double _tolerance, double _radius )
{
  
  Xx0init( _tolerance, _radius );


}


std::vector<IVector> IReturnCorrectedOrbit( double _tolerance )
{
  std::vector<DVector> DxList( returnCorrectedOrbit( _tolerance ) );
  std::vector<IVector> xList( pm_count );

  for( int i = 0; i < pm_count; i++ )
    (xList[i]).resize( dim );

  for( int i = 0; i < pm_count; i++ )
  {
    for( int j = 1; j <= dim; j++ )
      (xList[i])(j) = (DxList[i])(j);
  }

  return xList;
}



void IintegrateOneTime( double _tolerance ) 
{
  std::vector<IVector> vect( IReturnCorrectedOrbit( _tolerance ) );
  for( int i = 0; i < pm_count-1; i++ )
  {

    IPoincareMap pm( ISolver, ISection[i+1], poincare::MinusPlus );
    interval temp_time = 0.;

    IMatrix monodromyMatrix(3,3);
    IMatrix poincareDer(3,3);

    C1Rect2Set setToIntegrate( vect[i] );

    IVector result = pm( setToIntegrate, monodromyMatrix, temp_time );
    poincareDer = pm.computeDP( result, monodromyMatrix, temp_time );
    

    if( ((ISection[i]).getOrigin())(3) > 0.05 && ((ISection[i]).getOrigin())(3) < 0.06 )
    {
      cout << (ISection[i]).getOrigin() << " " << maxDiam( result  - vect[i+1] )  << "\n";
      cout << poincareDer << "\n";
    }

  }
}



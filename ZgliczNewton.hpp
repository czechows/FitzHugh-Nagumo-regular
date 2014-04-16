/* -----------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing Piotr Zgliczy{\'n}ski ideas
 * on how to compute the inverse derivative matrix
 * which occurs in problems of applying the interval Newton operator
 * to composition of Poincar{\'e} maps.
 * Methods used here are elements 
 * of the FhnIntervalNewton class.
 * ------------------------------------------ */

void Xx0init( double _tolerance, double _radius )
{

  std::vector<DVector> x0_double( newtonAlgorithm( _tolerance ) );

  for( int i = 0; i < pm_count; i++ )
  {
    for( int j = 1; j <= fast_dim; j++ )
      (x0[i])(j) = (x0_double[i])(j);
  }

  for( int i = 0; i < pm_count; i++ )
    X[i] = (x0[i])*MpInterval( 1. - _radius, 1. + _radius );
}

MpIMatrix Ai( int i ) 
{
  MpIVector pm_result(3);
  MpInterval time(0.);

  MpIMatrix monodromyMatrix(3,3);
  MpIMatrix poincareDer(3,3);

  time = 0.;
  monodromyMatrix.setToIdentity(); // probably obsolete since it is done automatically (?) in the solver code

  if( i == 0 )
  {
    MpIVector tempvector(3);
    tempvector[0] = (X[pm_count-1])(1);
    tempvector[1] = (X[pm_count-1])(2);
    tempvector[2] = MpInterval( 0., 0. );

    MpC1Rect2Set setToIntegrate( (ISection[pm_count-1]).getOrigin(), IP_list[pm_count-1], tempvector );

    MpIPoincareMap pm( ISolver, ISection[0], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, monodromyMatrix, time );
    poincareDer = pm.computeDP( pm_result, monodromyMatrix, time );
  
    poincareDer = inverseMatrix( IP_list[0] )*poincareDer*( IP_list[pm_count-1] );
  }
  else if( i > 0 )
  {
    MpIVector tempvector(3);
    tempvector[0] = (X[i-1])(1);
    tempvector[1] = (X[i-1])(2);
    tempvector[2] = MpInterval( 0., 0. );

    MpC1Rect2Set setToIntegrate( (ISection[i-1]).getOrigin(), IP_list[i-1], tempvector );

    MpIPoincareMap pm( ISolver, ISection[i], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, monodromyMatrix, time );
    poincareDer = pm.computeDP( pm_result, monodromyMatrix, time );
  
    poincareDer = inverseMatrix( IP_list[i] )*poincareDer*( IP_list[i-1] );
  }
  else
    throw "DERIVATIVE MATRIX INDEX OUT OF RANGE! \n";

  MpIMatrix poincareDerRed(2,2);

  poincareDerRed[0][0] = poincareDer[0][0];
  poincareDerRed[0][1] = poincareDer[0][1];
  poincareDerRed[1][0] = poincareDer[1][0]; 
  poincareDerRed[1][1] = poincareDer[1][1];

  return poincareDerRed;
}

MpIMatrix Akl( int k, int l )
{
  MpIMatrix result(2,2);
  result.setToIdentity();

  if( l < k )
    throw "AKL INDEX ERROR! \n";

  for( int i = k; i < l ; i++ )
  {

    result = Ai( i % pm_count )*result;
    cout << i << "\n";
 //   if( vectalg::containsZero( MpIVector( {matrixAlgorithms::det( result )} )) )
  //    cout << i << " " << matrixAlgorithms::det( MpIMatrix::Identity( fast_dim ) - result ) << " " << matrixAlgorithms::det( result ) << " " << (ISection[i]).getOrigin() << "\n";
/*
    if( i % pm_count == 107 )
      cout << result << " " << matrixAlgorithms::det( result ) << "\n";


    if( i % pm_count == 107 )
      cout << Ai( i % pm_count ) << " " <<  matrixAlgorithms::det( Ai( i % pm_count ) )<< "\n" << result << " " << matrixAlgorithms::det( result ) << "\n ";// << matrixAlgorithms::det( MpIMatrix::Identity(fast_dim) - result )  << "\n";
  //  cout << result << "\n";*/
  }

  return result;
}


MpIVector fi( int i, MpIVector _x0 )
{
  MpIVector pm_result( dim );
  MpInterval time(0.);

  MpIVector fi_eval( fast_dim ); 
  
  MpIVector tempvector(3);
  tempvector[0] = _x0(1);
  tempvector[1] = _x0(2);
  tempvector[2] = MpInterval( 0., 0. );

  if( i == pm_count - 1 )
  {
    MpC0Rect2Set setToIntegrate( (ISection[pm_count-1]).getOrigin(), IP_list[pm_count-1], tempvector );

    MpIPoincareMap pm( ISolver, ISection[0], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, (ISection[0]).getOrigin(), inverseMatrix(IP_list[0]), time );
  }
  else if( i < pm_count - 1 )
  {
    MpC0Rect2Set setToIntegrate( (ISection[i]).getOrigin(), IP_list[i], tempvector );

    MpIPoincareMap pm( ISolver, ISection[i+1], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, (ISection[i+1]).getOrigin(), inverseMatrix(IP_list[i+1]), time ); 
  }
  else
    throw "MAP INDEX OUT OF RANGE! \n";

  for( int i = 0; i < fast_dim; i++ )
    fi_eval[i] = pm_result[i];

  return fi_eval;
}


MpIVector computeInvDFXtimesFx0( int i )   
{ 
  MpIVector sum( fast_dim );
  sum.clear();
  sum = MpIVector({ MpInterval(1., 1.), MpInterval(1., 1.) });

 // for( int j = 1; j <= pm_count; j++ )
  //  sum = sum + (Akl( i + j, i + pm_count ))*( x0[ (i + j)%pm_count ] - fi( (i + j - 1)%pm_count , x0[ (i + j - 1)%pm_count ] ) ); 

  MpIMatrix setToCout(  MpIMatrix::Identity( fast_dim ) - Akl( i, i + pm_count )  );
  /*
  for( int i = 1; i <= 2; i++ )
  {
    for( int j = 1; j <= 2; j++ )
    {
      setToCout( i, j ) = ( setToCout( i, j ) ).leftBound();
    }
  }

*/
  cout << setToCout << "  " << matrixAlgorithms::det(setToCout) << "\n";

  MpIVector result( matrixAlgorithms::gauss( MpIMatrix::Identity( fast_dim ) - Akl( i, i + pm_count ), sum ) );
  return result;
}


void proveExistenceOfOrbit( double _tolerance, double _radius )
{
  
  Xx0init( _tolerance, _radius );

  for( int i = 0; i < 1; i++ )
  {
    MpIVector N_i( x0[i] - computeInvDFXtimesFx0(i) ); 

    if( subsetInterior( x0[i], N_i ) )
     cout << "Warning! x0[" << i << "] is not contained in N[" << i << "]! \n";

    if( subsetInterior( N_i, X[i] ) )
      cout << "N[" << i << "] = " << N_i << " is a subset of X[" << i << "] = " << X[i] << " !\n";
    else
      cout << "EXISTENCE OF PERIODIC ORBIT NOT PROVEN" << " diam(N[" << i << "]) = " << maxDiam(N_i) << " diam(X[" << i << "]) = " << maxDiam(X[i]) << "! \n"; 
  }
}



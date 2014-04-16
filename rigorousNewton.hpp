/* -----------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing the rigorous Newton method
 * for finding periodic solutions
 * ------------------------------------------ */



class FhnIntervalNewton : public FhnFindPeriodicOrbit
{
  public:
  MpIMap MpIVectorField;
  MpITaylor ISolver;
  std::vector<MpIMatrix> IP_list;
  std::vector<MpIAffineSection> ISection;
  std::vector<MpIVector> x0;
  std::vector<MpIVector> X;
 
  FhnIntervalNewton( const int _pm_count )
    : FhnFindPeriodicOrbit( _pm_count ),
      MpIVectorField( IFhn_vf ),
      ISolver( MpIVectorField, order),
      IP_list( pm_count ),
      ISection( pm_count, MpIAffineSection(  IxPrecomputed[0],  IxPrecomputed[1]-IxPrecomputed[0] ) ),
      x0( pm_count ),
      X( pm_count )
  {
    for( int i = 0; i < pm_count - 1; i++ )
    {
      (IP_list[i]).resize( dim, dim );
      (IP_list[i]).setToIdentity();
      x0[i].resize( fast_dim );
      X[i].resize( fast_dim );

      for( int j = 1; j <= dim; j++ )
       (IP_list[i])( j, dim ) = ( MpIVector( IxPrecomputed[i+1]-IxPrecomputed[i] ) )(j);

      IOrthogonalizeRelativeColumn( IP_list[i], dim - 1 );

      ISection[i] = MpIAffineSection( IxPrecomputed[i], IxPrecomputed[i+1]-IxPrecomputed[i] ); // ? MpIVectors ?
    }
 
    (IP_list[ pm_count - 1 ]).resize( dim, dim );
    (IP_list[ pm_count - 1 ]).setToIdentity();
    x0[ pm_count - 1 ].resize( fast_dim );
    X[ pm_count - 1 ].resize( fast_dim );

    for( int j = 1; j <= dim; j++ )
      (IP_list[ pm_count - 1 ])( j, dim ) = ( IxPrecomputed[0]-IxPrecomputed[pm_count-1] )(j);
 
    ISection[ pm_count - 1 ] = MpIAffineSection( IxPrecomputed[pm_count - 1], IxPrecomputed[0]-IxPrecomputed[ pm_count - 1 ] );

    IOrthogonalizeRelativeColumn( IP_list[pm_count-1], dim - 1 );

  }
/*
  MpIMatrix computeDerMatrix( MpIVector x )   
  {
    MpIVector pm_result(3);
    MpInterval time(0.);

    MpIMatrix monodromyMatrix(3,3);
    MpIMatrix poincareDer(3,3);

    MpIMatrix derMatrix( x.dimension(), x.dimension() );
    derMatrix.setToIdentity();

    for( int i = 0; i < pm_count-1; i++ )
    {
      time = 0.;
      monodromyMatrix.setToIdentity(); // probably obsolete since it is done automatically (?) in the solver code
  
      MpC1Rect2Set setToIntegrate( (ISection[i]).getOrigin(), IP_list[i], MpIVector({ x[2*i], x[2*i + 1], 0. }) );

      MpIPoincareMap pm( ISolver, ISection[i+1], poincare::MinusPlus );
  
     // cout << "Integration to section " << i+1 << " in progress, setToIntegrate = " <<  MpIVector({ x[2*i], x[2*i + 1], 0. })  << " \n";

      pm_result = pm( setToIntegrate, monodromyMatrix, time );
      poincareDer = pm.computeDP( pm_result, monodromyMatrix, time );
      
      poincareDer = inverseMatrix( IP_list[i+1] )*poincareDer*( IP_list[i] );

      derMatrix[2*i+2][2*i] = -poincareDer[0][0];
      derMatrix[2*i+2][2*i+1] = -poincareDer[0][1];
      derMatrix[2*i+3][2*i] = -poincareDer[1][0]; 
      derMatrix[2*i+3][2*i+1] = -poincareDer[1][1];
    }

    time = 0.;
    monodromyMatrix.setToIdentity();
  
    MpC1Rect2Set setToIntegrate( (ISection[pm_count - 1]).getOrigin(), IP_list[pm_count - 1], MpIVector({ x[2*pm_count - 2], x[2*pm_count - 1], 0. }) );

    MpIPoincareMap pm( ISolver, ISection[0], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, monodromyMatrix, time );
    poincareDer = pm.computeDP( pm_result, monodromyMatrix, time);
 
    poincareDer = inverseMatrix( IP_list[0] )*poincareDer*( IP_list[pm_count - 1] );

    derMatrix[0][2*pm_count - 2] = -poincareDer[0][0];
    derMatrix[0][2*pm_count - 1] = -poincareDer[0][1];
    derMatrix[1][2*pm_count - 2] = -poincareDer[1][0];
    derMatrix[1][2*pm_count - 1] = -poincareDer[1][1];

    return derMatrix;
  }


  MpIVector computeF( MpIVector x )   
  { 
    MpIVector pm_result(3);
    MpInterval time(0.);

    MpIVector x_eval( x.dimension() ); // evaluates function xi - P_{i}(x_{i-1})

    for( int i = 0; i < pm_count-1; i++ )
    {
      time = 0.;
      MpC0Rect2Set setToIntegrate( (ISection[i]).getOrigin(), IP_list[i], MpIVector({ x[2*i], x[2*i + 1], 0 }) );
      
      MpIPoincareMap pm( ISolver, ISection[i+1], poincare::MinusPlus );
 
     // cout << "F Integration to section " << i+1 << " in progress \n";

      pm_result = pm(setToIntegrate, (ISection[i+1]).getOrigin(), inverseMatrix(IP_list[i+1]), time );
     
      x_eval[2*i + 2] = x[2*i + 2] - pm_result[0];
      x_eval[2*i + 3] = x[2*i + 3] - pm_result[1];
    }

    time = 0.;
    MpC0Rect2Set setToIntegrate( (ISection[pm_count - 1]).getOrigin(), IP_list[pm_count - 1], MpIVector({ x[2*pm_count - 2], x[2*pm_count - 1], 0 }) );
    
    MpIPoincareMap pm( ISolver, ISection[0], poincare::MinusPlus );

    pm_result = pm( setToIntegrate, (ISection[0]).getOrigin(), inverseMatrix(IP_list[0]), time );
 
    x_eval[0] = x[0] - pm_result[0];
    x_eval[1] = x[1] - pm_result[1];

    return x_eval;
  }
*/
  #include "ZgliczNewton.hpp"
/*
  void proveExistenceOfOrbit( double _tolerance, double _radius )
  {
    
    Xx0init( _tolerance, _radius );


 //   if( !subsetInterior(x0,X) )
  //    throw "x0 not in interior of X";

    MpIVector N( x0 - capd::matrixAlgorithms::gauss(  computeDerMatrix( X ) , computeF( x0 ) ) ); 

    if( subsetInterior( x0, N ) )
       cout << "ERROR";

    if( subsetInterior( N, X ) )
      cout << "Existence of periodic orbit proven \n";
    else
      cout << "Existence of periodic orbit not proven" << " diam(N) = " << N << " diam(X) = " << maxDiam(X) << " \n"; 
  }
*/
};




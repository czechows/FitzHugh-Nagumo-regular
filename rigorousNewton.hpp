/* -----------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing the rigorous Newton method
 * for finding periodic solutions
 * ------------------------------------------ */



class FhnIntervalNewton : public FhnFindPeriodicOrbit
{
  public:
  IMap IVectorField;
  ITaylor ISolver;
  std::vector<IMatrix> IP_list;
  std::vector<IAffineSection> ISection;
  std::vector<IVector> x0;
  std::vector<IVector> X;
 
  FhnIntervalNewton( const int _pm_count )
    : FhnFindPeriodicOrbit( _pm_count ),
      IVectorField( IFhn_vf ),
      ISolver( IVectorField, order),
      IP_list( pm_count ),
      ISection( pm_count, IAffineSection(  IxPrecomputed[0],  IxPrecomputed[1]-IxPrecomputed[0] ) ),
      x0( pm_count ),
      X( pm_count )
  {
    for( int i = 0; i < pm_count - 1; i++ )
    {
      (IP_list[i]).resize( dim, dim );
      (IP_list[i]).setToIdentity();
      x0[i].resize( fast_dim );
      X[i].resize( fast_dim );

      ISection[i] = IAffineSection( IxPrecomputed[i], IxPrecomputed[i+1]-IxPrecomputed[i] ); // ? IVectors ?
    }
 
    (IP_list[ pm_count - 1 ]).resize( dim, dim );
    (IP_list[ pm_count - 1 ]).setToIdentity();
    x0[ pm_count - 1 ].resize( fast_dim );
    X[ pm_count - 1 ].resize( fast_dim );

    ISection[ pm_count - 1 ] = IAffineSection( IxPrecomputed[pm_count - 1], IxPrecomputed[0]-IxPrecomputed[ pm_count - 1 ] );
  }

  #include "covering.hpp"

};




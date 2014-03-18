#include <iostream>
#include "capd/capdlib.h"

using std::cout;
using namespace capd;
//using namespace matrixAlgorithms;
//using namespace dynsys;


DMap Fhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
// FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)
DMap Fhn_fast("par:theta,v;var:u,w;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v);");
// fast subsystem of the FitzHugh-Nagumo system
  
const int pm_count(22); 
const int pm_count_up(8); // number of points on the upper branch

const double accuracy = 1e-12;           // accuracy for nonrigorous numerics (i.e. approximation of the slow manifold)
const int order = 6;

#include "numericsD.hpp"   // Warning! When changing the vector field, one needs to make manual changes in this header file (class FhnBifurcation)!


class FhnFindPeriodicOrbit
{
public:
  DMap vectorField;
  DTaylor solver;
  DVector GammaUL;
  DVector GammaUR;
  DVector GammaDL;
  DVector GammaDR;
  DVector [pm_count] x_start;

  FhnFindPeriodicOrbit( double eps, double theta, double _tolerance, const DVector& _GammaUL, const DVector& _GammaUR, const DVector& _GammaDL, const DVector& _GammaDR )
    : vectorField(Fhn_vf),
      solver(vectorField, order),
      GammaUL( _GammaUL ),
      GammaUR( _GammaUR ),
      GammaDL( _GammaDL ),
      GammaDR( _GammaDR )
  {
    vectorField.setParameter( "theta", _theta );
    vectorField.setParameter( "eps", _eps );

    GammaQuad_correct( _theta, GammaUL, GammaDL, GammaUR, GammaDR );     // we correct the initial guesses by nonrigorous Newtons methods (see numericsD.hpp)
  }

}


DVector Eq_correct(DVector& guess, double v)              // corrects initial guesses of u so they are closer to real equilibria, w is always 0
{
  vectorField.setParameter("v", v);
  double error = 1.;
  double result = guess[0];
  double oldresult = result;
  while(error > accuracy)
  {
    oldresult = result;
    result = result - ( vectorField(DVector({result, 0.}) )[1] / vectorField[ DVector({result, 0.}) ][1][0] );  // Newton algorithm to calculate zeroes of the vector field - w is always 0.,
                                                                                                                // derivative is of the second equation wr. to first variable - u
    error = abs(oldresult - result);
  }
  DVector new_Eq(2);
  new_Eq[0] = result;
  new_Eq[1] = 0.;
  return new_Eq;
}

DVector xPrecomputed [pm_count] =  
{
 0.0250442158334208;
 0.0988076360184288;
}; 
 


DVector oneNewtonStep( DVector x )   // one Newton step for function xi - P_{i}(x_{i-1}), i in mathbbZ
{
 
  DTaylor solver(Fhn_vf, order);

  DAffineSection *section [pm_count];
  DPoincareMap *pm [pm_count];

  for( int i = 0; i < pm_count; i++ )    // could have used coordinate sections but doesnt matter in the end
  {
    section[i] = new DAffineSection( xPrecomputed[i], DVector(0.,0.,1.) );
    if( i < pm_count_up )
      pm[i] = new DPoincareMap( solver, *section[i], poincare::MinusPlus );
    else
      pm[i] = new DPoincareMap( solver, *section[i], poincare::PlusMinus );

  }

  DVector pm_result(3);
  DVector tempvar(3);
  double time(0.);

  DMatrix monodromyMatrix(3,3);
  DMatrix poincareDer(3,3);

  DVector x_eval( x.dimension() ); // evaluates function xi - P_{i}(x_{i-1})
  DMatrix derMatrix( x.dimension(), x.dimension() );
  derMatrix.setToIdentity();

  for( int i = 0; i < pm_count-1; i++ )
  {
    monodromyMatrix.setToIdentity(); // probably obsolete since it is done automatically in the solver (?) code
    
    tempvar[0] = x[2*i];
    tempvar[1] = x[2*i + 1];
    tempvar[2] = xPrecomputed[i][2];

    pm_result = (*pm[i+1])(tempvar, monodromyMatrix, time);
    poincareDer = (*pm[i+1]).computeDP(pm_result, monodromyMatrix, time);

    x_eval[2*i + 2] = x[2*i + 2] - pm_result[0];
    x_eval[2*i + 3] = x[2*i + 3] - pm_result[1];

    derMatrix[2*i+2][2*i] = -poincareDer[0][0];
    derMatrix[2*i+2][2*i+1] = -poincareDer[0][1];
    derMatrix[2*i+3][2*i] = -poincareDer[1][0]; 
    derMatrix[2*i+3][2*i+1] = -poincareDer[1][1];
  }


  monodromyMatrix.setToIdentity();
 
  tempvar[0] = x[2*pm_count - 2];
  tempvar[1] = x[2*pm_count - 1];
  tempvar[2] = xPrecomputed[pm_count - 1][2];

  pm_result = (*pm[0])( tempvar, monodromyMatrix, time );
  poincareDer = (*pm[0]).computeDP(pm_result, monodromyMatrix, time);

  x_eval[0] = x[0] - pm_result[0];
  x_eval[1] = x[1] - pm_result[1];

  derMatrix[0][2*pm_count - 2] = -poincareDer[0][0];
  derMatrix[0][2*pm_count - 1] = -poincareDer[0][1];
  derMatrix[1][2*pm_count - 2] = -poincareDer[1][0];
  derMatrix[1][2*pm_count - 1] = -poincareDer[1][1];


  for( int i = 0; i < pm_count; i++ )
  {
    delete pm[i];
    delete section[i];
  }

  DVector result(x.dimension());
  result = x - capd::matrixAlgorithms::gauss( derMatrix, x_eval );

  return result;
};


DVector newtonAlgorithm( DVector x_start, double tolerance )
{
  double error = 1.;
  
  DVector x0( x_start );
  DVector x1( x_start );

  while( error > tolerance )
  {
    x0 = x1;
    x1 = oneNewtonStep( x1 );

    error = sqrt( scalarProduct(x1 - x0, x1 - x0) );
    cout << error << "\n";
  }

  return x1;
};





int main(){

  cout.precision(9);

  double theta = double(61.)/100.;
  double eps = double(0.001); 
 
  /*
  if( !(GammaUL[0] > GammaDL[0] && GammaUR[0] > GammaDR[0] && GammaUR[2] > GammaUL[2] && GammaDR[2] > GammaDL[2] ) )
    throw "NEWTON CORRECTION METHOD FOR CORNER POINTS ERROR! \n";
 */

  // some points on the orbit precomputed with MATCONT, available in file JJCperiodicorbit.dat
  // size of array x should be same as pm_count !!
  // there should be the same number of points on the upper branch as on the lower branch
  // ARTIFICIAL LINE - not from Matcont (not enough available) but just copypasted another line with different v value inserted to prevent blowups

  double time = 0.;

  DVector x_start( 2*pm_count );

  for(int i = 0; i < pm_count; i++ )
  {
    x_start[ 2*i ] = xPrecomputed[i][0];
    x_start[ 2*i + 1 ] = xPrecomputed[i][1];
  }


  DVector x( oneNewtonStep(oneNewtonStep( x_start )) );

  for(int i = 0; i < 2*pm_count - 1; i++)
    cout << x_start[i] << " ? " << x[i] << "\n";


 // newtonAlgorithm( x_start, 1e-5 );


  
  return 0;
} 


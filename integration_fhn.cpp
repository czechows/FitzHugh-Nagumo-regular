#include <iostream>
#include <vector>
#include "capd/capdlib.h"

using std::cout;
using namespace capd;
using matrixAlgorithms::inverseMatrix;

DMap Fhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
IMap IFhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
// FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)


  
double eps = double(1e-3); 
const int order = 6;
const int precomp_factor = 31; 

#include "matcontPrecomputedOrbit.hpp"

std::vector<DVector> xPrecomputed;
std::vector<IVector> IxPrecomputed;

#include "auxiliaries.hpp"
#include "nonrigorousNewton.hpp"
#include "rigorousNewton.hpp"

int main(){

  cout.precision(15);

  double theta = double(61.)/100.;
  double tolerance = 1e-14;
  double radius = double(1e-14);

  Fhn_vf.setParameter( "theta", theta );
  Fhn_vf.setParameter( "eps", eps );

  IFhn_vf.setParameter( "theta", theta );
  IFhn_vf.setParameter( "eps", eps );

  xPrecomputedFill();
  IxPrecomputedFill();

  const int pm_count = xPrecomputed.size();
  FhnFindPeriodicOrbit FindPeriodicOrbit( pm_count );

 
  try
  {
    FhnIntervalNewton IntervalNewton( pm_count );
    IntervalNewton.proveExistenceOfOrbit( tolerance, radius );
  }
  catch(const char* Message)
  {
    cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES THETA=" << theta << " AND EPS=" << eps << " NOT VERIFIED! \n";
  }

/*
  FindPeriodicOrbit.integrateOneTime( tolerance );

  std::vector<DVector> result = FindPeriodicOrbit.returnCorrectedOrbit( tolerance );
 // for( unsigned int i = 0; i < xPrecomputed.size(); i++ )
 //   cout << IxPrecomputed[i] << xPrecomputed[i] << "\n" << " " << result[i] << "\n";
*/
  return 0;
} 


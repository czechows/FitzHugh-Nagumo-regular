#include <iostream>
#include <vector>
#include "capd/capdlib.h"

using std::cout;
using namespace capd;
using matrixAlgorithms::inverseMatrix;

DMap Fhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
IMap IFhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 

DMap Fhn_vf_rev("par:theta,eps;var:u,w,v;fun:-w,(-2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(-eps/theta)*(u-v);"); 
IMap IFhn_vf_rev("par:theta,eps;var:u,w,v;fun:-w,(-2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(-eps/theta)*(u-v);"); 
// FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)


  
const int order = 6;
const int precomp_factor = 5; 

#include "matcontPrecomputedOrbit.hpp"

std::vector<DVector> xPrecomputed;
std::vector<IVector> IxPrecomputed;

#include "auxiliaries.hpp"
#include "nonrigorousNewton.hpp"
#include "covering.hpp"
#include "continuation.hpp"

int main(){

  cout.precision(15);

  interval theta = interval(61.)/100.;
  interval eps = interval(1e-4, 1e-3); 
  double tolerance = 1e-14;
  double radius = double(1e-8);

  xPrecomputedFill();
  IxPrecomputedFill();

 // const int pm_count = xPrecomputed.size();
  // FhnFindPeriodicOrbit FindPeriodicOrbit( pm_count );
 
 // FhnCovering cov( xPrecomputed );
 // cov.proveExistenceOfOrbit( theta, eps, tolerance, radius );

  bool isEpsIncreasing( 0 );

  try
  {
    continueOrbitWithEps( theta, eps, isEpsIncreasing, xPrecomputed, tolerance, radius, 1e-8 );
  }
  catch(const char* Message)
  {
    cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES THETA=" << theta << " AND EPS=" << eps << " NOT VERIFIED! \n";
  }

  return 0;
} 


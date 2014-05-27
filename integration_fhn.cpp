#include <iostream>
#include <vector>
#include <string.h>
#include "capd/capdlib.h"

using std::cout;
using namespace capd;
using matrixAlgorithms::inverseMatrix;

DMap *Fhn_vf;
IMap *IFhn_vf;

DMap *Fhn_vf_rev;
IMap *IFhn_vf_withEps;
  
int order = 6; // has to be low to avoid crossing sections in one step
int rig_order = 10; // higher better?

const int precomp_factor = 5; 

#include "matcontPrecomputedOrbit.hpp"

std::vector<DVector> xPrecomputed;
std::vector<IVector> IxPrecomputed;

#include "auxiliaries.hpp"
#include "nonrigorousNewton.hpp"
#include "covering.hpp"
#include "continuation.hpp"

int main(){

  Fhn_vf = new DMap("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
  IFhn_vf = new IMap("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 

  Fhn_vf_rev = new DMap("par:theta,eps;var:u,w,v;fun:-w,(-2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(-eps/theta)*(u-v);"); 
  IFhn_vf_withEps = new IMap("par:theta;var:u,w,v,eps;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v),0;");
  // FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)
 

  cout.precision(15);

  interval theta = interval(61.)/100.;
  interval eps = interval(1.5e-4, 1e-3); 
  double tolerance = 1e-13;
  double radius = double(1e-6);

  xPrecomputedFill();
  IxPrecomputedFill();


 // const int pm_count = xPrecomputed.size();
  // FhnFindPeriodicOrbit FindPeriodicOrbit( pm_count );
 
 // FhnCovering cov( xPrecomputed );
 // cov.proveExistenceOfOrbit( theta, eps, tolerance, radius );

  bool isEpsIncreasing( 0 );

  FhnValidatedContinuation cont( theta, eps, xPrecomputed, isEpsIncreasing, tolerance, radius );
  cont.continueOrbitWithEps();


 /* try
  {
    continueOrbitWithEps( theta, eps, isEpsIncreasing, xPrecomputed, tolerance, radius, 1e-8 );
  }
  catch(const char* Message)
  {
    cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES THETA=" << theta << " AND EPS=" << eps << " NOT VERIFIED! \n";
  }
*/

  delete Fhn_vf;
  delete IFhn_vf;

  delete Fhn_vf_rev;
  delete IFhn_vf_withEps;

  return 0;
} 


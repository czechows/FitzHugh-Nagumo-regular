#include <iostream>
#include <vector>
#include<fstream>
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
int rig_order = 10; // too low & too high are bad 

const int precomp_factor = 5; 
bool wasKrawczykNeeded = 0;

#include "matcontPrecomputedOrbit.hpp"
#include "savedOrbit.hpp"

std::vector<DVector> xPrecomputed;
std::vector<IVector> IxPrecomputed;

#include "auxiliaries.hpp"
#include "nonrigorousNewton.hpp"
#include "covering.hpp"
#include "krawczyk.hpp"
#include "continuation.hpp"

int main(){

  Fhn_vf = new DMap("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
  IFhn_vf = new IMap("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 

  Fhn_vf_rev = new DMap("par:theta,eps;var:u,w,v;fun:-w,(-2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(-eps/theta)*(u-v);"); 
  IFhn_vf_withEps = new IMap("par:theta;var:u,w,v,eps;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v),0;");
  // FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)
 

  cout.precision(15);

  interval theta = interval(61.)/100.;
  interval eps = interval("0.7e-4", "1e-3"); 
  double tolerance = 1e-13;
  double radius = double(5e-5);
  double startIncrementSize = 1e-6;
  interval integrationTimeBound = interval( 0.2, 4. );
  bool isEpsIncreasing( 0 );

/*  std::ifstream readSavedEpsRange;  // uncomment and switch isEpsIncreasing if necessary to start from where the proof for a given eps range ended
 *                                    // you need to run the program once and it has to terminate prematurely for this option to work
  readSavedEpsRange.open("savedEpsRange", std::ios::binary );
  interval savedEpsRange;
  hexRead( readSavedEpsRange, savedEpsRange );
  readSavedEpsRange.close();
  cout << "SAVED EPS RANGE READ: " << savedEpsRange << "\n";
  FhnValidatedContinuation cont( theta, savedEpsRange, savedOrbit, isEpsIncreasing, tolerance, radius, startIncrementSize, integrationTimeBound );
  cont.continueOrbitWithEps(); */

  xPrecomputedFill();
  IxPrecomputedFill();

 // FhnValidatedContinuation cont( theta, eps, xPrecomputed, isEpsIncreasing, tolerance, radius, startIncrementSize, integrationTimeBound );
 // cont.continueOrbitWithEps();

  isEpsIncreasing = 1;
 
  eps = interval("1e-3", "1.5e-3");

  FhnValidatedContinuation cont2( theta, eps, xPrecomputed, isEpsIncreasing, tolerance, radius, startIncrementSize, integrationTimeBound );
  cont2.continueOrbitWithEps();
 
 /*
  isEpsIncreasing = 0;  // speed test for interval already proven in singular range
  eps = interval("1.5e-4", "1.1e-4");
  FhnValidatedContinuation cont_test( theta, eps, savedOrbit, isEpsIncreasing, tolerance, radius, startIncrementSize, integrationTimeBound );
  cont_test.continueOrbitWithEps();
  */
  cout << wasKrawczykNeeded << "\n";
  delete Fhn_vf;
  delete IFhn_vf;

  delete Fhn_vf_rev;
  delete IFhn_vf_withEps;

  return 0;
} 


#include <iostream>
#include <vector>
#include<fstream>
#include <string.h>
#include "capd/capdlib.h"
#include "time.h"

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
#include "intervalNewton.hpp"
#include "continuation.hpp"

int main(){
  time_t start_cont,end_cont;
  time (&start_cont);

  Fhn_vf = new DMap("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
  IFhn_vf = new IMap("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 

  Fhn_vf_rev = new DMap("par:theta,eps;var:u,w,v;fun:-w,(-2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(-eps/theta)*(u-v);"); 
  IFhn_vf_withEps = new IMap("par:theta;var:u,w,v,eps;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v),0;");
  // FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)

  cout.precision(15);

  interval theta = interval(61.)/100.;
  interval eps = interval("1.5e-4", "1e-3"); 
  double tolerance = 1e-13;
  double radius = 1e-6;
  double startIncrementSize = 1e-6;
  interval integrationTimeBound = interval( 0.2, 4. );
  bool isEpsIncreasing( 0 );

/*  std::ifstream readSavedEpsRange;  // uncomment and switch isEpsIncreasing if necessary to start from where the proof for a given eps range ended
 *                                    // you need to run the program once for this option to work (it can terminate prematurely)
  readSavedEpsRange.open("savedEpsRange", std::ios::binary );
  interval savedEpsRange;
  hexRead( readSavedEpsRange, savedEpsRange );
  readSavedEpsRange.close();
  cout << "SAVED EPS RANGE READ: " << savedEpsRange << "\n";
  FhnValidatedContinuation cont( theta, savedEpsRange, savedOrbit, isEpsIncreasing, tolerance, radius, startIncrementSize, integrationTimeBound );
  cont.continueOrbitWithEps(); */

  xPrecomputedFill();
  IxPrecomputedFill();


  FhnValidatedContinuation cont_down( theta, eps, xPrecomputed, isEpsIncreasing, tolerance, radius, startIncrementSize, integrationTimeBound );
  cont_down.continueOrbitWithEps();
 
  
 /*  // speed test for interval already proved in the singular range
  time_t start_test,end_test;
  time( &start_test );

  std::vector<DVector> testStartApprox = cont_down.numericOrbitGuess;

  isEpsIncreasing = 0;    
  eps = interval("1.1e-4", "1.5e-4");
  FhnValidatedContinuation cont_test( theta, eps, testStartApprox, isEpsIncreasing, tolerance, radius, startIncrementSize, integrationTimeBound );
  cont_test.continueOrbitWithEps();
  time( &end_test );
  double dif_test = difftime( end_test, start_test );
  cout << "Elapsed time for the validated continuation in parameter range epsilon= " << eps << " is " << dif_test << " seconds. \n";
 */

  isEpsIncreasing = 1;
  eps = interval("1e-3", "1.5e-3");

  FhnValidatedContinuation cont_up( theta, eps, xPrecomputed, isEpsIncreasing, tolerance, radius, startIncrementSize, integrationTimeBound );
  cont_up.continueOrbitWithEps();

  time (&end_cont);
  double dif_cont = difftime( end_cont, start_cont );
  cout << "Elapsed time for the validated continuation proof is " << dif_cont << " seconds. \n \n";
 
 
  time_t start_newt,end_newt;
  time (&start_newt);

  eps = interval(15./10000.);
  FhnIntervalNewton newton( cont_up.numericOrbitGuess, 1 ); // no extra discretizations for integration
 
  newton.proveExistenceOfOrbitWithNewton( theta, eps, tolerance, radius );
  cout << "EXISTENCE OF THE PERIODIC ORBIT FOR EPSILON = " << eps << " PROVED!\n";
  cout << "\nWas Krawczyk operator necessary? (0 no, 1 yes): " << wasKrawczykNeeded << "\n";

  time(&end_newt);
  double dif_newt = difftime( end_newt, start_newt );
  cout << "Elapsed time for the interval Newton based proof for epsilon = " << eps << " is " << dif_newt << " seconds. \n";

  delete Fhn_vf;
  delete IFhn_vf;

  delete Fhn_vf_rev;
  delete IFhn_vf_withEps;

  return 0;
} 


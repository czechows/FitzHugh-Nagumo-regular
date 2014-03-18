#include <iostream>
#include "capd/capdlib.h"

using std::cout;
using namespace capd;
//using namespace matrixAlgorithms;
//using namespace dynsys;


DMap Fhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
// FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)

  
const int pm_count(22); 
const int pm_count_up(8); // number of points on the upper branch

int order(6);

DVector xPrecomputed [pm_count] =  
{
  // UPPER BRANCH
  { 0.556950455687874, 0.0841927072526861, 0.0261922195055318 }, // line 26
  { 0.646321880535920, 0.0741669547976277, 0.0272529736252821 }, // line 27
  { 0.852424959432130, 0.0308841164972789, 0.0320456651771511 }, // line 30
  { 0.938218000680795, -0.000512156800739995, 0.0463164555944456 }, // line 35
  { 0.909646535575510, -0.00234570819002980, 0.0678386186723165 }, // line 40
  { 0.815437941275343, -0.0171921631995435, 0.0945326937197397 }, // line 51
  { 0.793278903418992, -0.0233077299558358, 0.0958150707254695 }, // line 52
  { 0.729656580361698, -0.0401973081151965, 0.0980516028684074 }, // line 53
  // LOWER BRANCH
 // { -0.202908805592454, -0.0114221719391439, 0.100259409958824 }, // line 67   
  { -0.227432732980653, -0.00243159476850949, 0.0981212128291538 }, // line 69
  { -0.231543217009309, 0.000410400074557027, 0.0947645695930721 }, // line 70
  { -0.219312204383719, 0.000653723496723840, 0.0850446157980756 }, // line 73
  { -0.219312204383719, 0.000653723496723840, 0.078 }, // ARTIFICIAL LINE
  { -0.202697889812905, 0.000682372607947186, 0.0733638442208183 }, // line 74
  { -0.202697889812905, 0.000682372607947186, 0.066 }, // ARTIFICIAL LINE 
  { -0.186169211569390, 0.000659705066224590, 0.0627989263809506 }, // line 75
  { -0.186169211569390, 0.000659705066224590, 0.058 }, // ARTIFICIAL LINE
  { -0.170000649177500, 0.000655692431231735, 0.0532962953136077 }, // line 76
  { -0.147155956252173, 0.000632724076983211, 0.0413332845718417 }, // line 78
  { -0.126768568416608, 0.000610072527618383, 0.0320156651382415 }, // line 81
  { -0.121754579585000, 0.000604777717904711, 0.0299156971165403 }, // line 3
  { -0.112278484380553, 0.000760432317307061, 0.0263345944454681 },  // line 8
  { -0.105291191426533, 0.00193367348280320, 0.0250009165557845 }, // line 11
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

  Fhn_vf.setParameter( "theta", theta );
  Fhn_vf.setParameter( "eps", eps );
  


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


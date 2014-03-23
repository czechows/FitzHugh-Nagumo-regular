#include <iostream>
#include <vector>
#include "capd/capdlib.h"

using std::cout;
using namespace capd;
//using namespace matrixAlgorithms;
//using namespace dynsys;


DMap Fhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
// FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)
DMap Fhn_fast("par:theta,v;var:u,w;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v);");
// fast subsystem of the FitzHugh-Nagumo system
  
double eps = double(1e-3); 


const int pm_count(150); 
const int pm_count_up(75); // number of points on the upper branch

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
  std::vector<DVector> xStartlist;
  DVector xSectionVCoordVector;
  int dim;
  int fast_dim;
  double time;

  FhnFindPeriodicOrbit( double _theta, const DVector& _GammaUL, const DVector& _GammaUR, const DVector& _GammaDL, const DVector& _GammaDR )
    : vectorField(Fhn_vf),
      solver(vectorField, order),
      GammaUL( _GammaUL ),
      GammaUR( _GammaUR ),
      GammaDL( _GammaDL ),
      GammaDR( _GammaDR ),
      xStartlist( pm_count ),
      xSectionVCoordVector( pm_count ),
      dim(3),
      fast_dim(dim - 1),
      time(0.)
  {
    GammaQuad_correct( _theta, GammaUL, GammaDL, GammaUR, GammaDR );     // we correct the initial guesses by nonrigorous Newtons methods (see numericsD.hpp)
    cout << "GammaQuad correction done \n";

    for( unsigned int i = 0; i < xStartlist.size(); i++ )
      (xStartlist[i]).resize( fast_dim );
  }

  DVector convertToVector( std::vector<DVector> list )
  {
    DVector result( fast_dim*list.size() );

    for( unsigned int i = 0; i < list.size(); i++)
    {
      for( int j = 0; j < fast_dim; j++)
       result[ i*fast_dim+j ] = (list[i])[j];
    }

    return result;
  };

  std::vector<DVector> convertToList( DVector vect )
  {
    std::vector<DVector> result( vect.dimension()/fast_dim );

    for( unsigned int i = 0; i < result.size(); i++)
    {
      for( int j = 0; j < fast_dim; j++)
       result[i][j] = vect[ i*fast_dim+j ];
    }

    return result;
  };

  DVector Eq_correct(DVector guess)              // corrects initial guesses of (u,w,v) so they are closer to points on the slow manifold, w is always 0
  {
    double error(1.);
    double result(guess[0]);
    double oldresult(result);
    while(error > accuracy)
    {
      oldresult = result;
      result = result - (vectorField( DVector(result, 0., guess[2]) )[1])/(vectorField[ DVector(result, 0., guess[2]) ][1][0]);  
                                                                          // Newton algorithm to calculate zeroes of the vector field - w is always 0.,
                                                                          // derivative is of the second equation with respect to first variable - u
      error = abs(oldresult - result);
    }
    DVector new_Eq(3);
    new_Eq[0] = result;
    new_Eq[1] = double(0.);
    new_Eq[2] = guess[2];
    return new_Eq;
  }

  void setStartListAndVCoordVector() // sets start list by correcting the guesses
  {

    for( int i = 0; i < pm_count_up; i++ )
    {
      double ti( double(i+1)/double( pm_count_up + 1 ) );
      DVector Gamma_i =  GammaUL + (GammaUR-GammaUL)*ti;
      DVector slowManPoint( Eq_correct( Gamma_i ) );
      
      if( i == pm_count_up - 1 )
        slowManPoint = DVector( { 0.687025509271046, -0.0504824728367381, slowManPoint[2] } );  
        //  slowManPoint = slowManPoint + (coordChange( vectorField, slowManPoint ))*DVector( {0., eps, 0. } );
        //{ 0.687025509271046, -0.0504824728367381, 0.0989946073927174 }
        
      xStartlist[i] = DVector( { slowManPoint[0], slowManPoint[1] } );
      xSectionVCoordVector[i] = slowManPoint[2];     
    }


    for( unsigned int i = pm_count_up; i < xStartlist.size(); i++ )
    {
      double ti( double(i-pm_count_up+1)/double( xStartlist.size() - pm_count_up + 1 ) );
      DVector Gamma_i =  GammaDR + (GammaDL-GammaDR)*ti;
      DVector slowManPoint( Eq_correct( Gamma_i ) );
      
      if( i == xStartlist.size() - 1 )
        slowManPoint = DVector({ -0.108264186750857, 0.00129680306590373, slowManPoint[2] });
      //  slowManPoint = DVector({ -0.108264186750857, 0.00129680306590373, 0.0254084822724654 });
      //  slowManPoint =  slowManPoint + (coordChange( vectorField, slowManPoint ))*DVector( {0., 0., 0. } );

      xStartlist[i] = DVector( { slowManPoint[0], slowManPoint[1] } );
      xSectionVCoordVector[i] = slowManPoint[2]; 
    }
  }

        
  std::vector<DVector> oneNewtonStep( std::vector<DVector> xList )   // one Newton step for function xi - P_{i}(x_{i-1}), i in mathbbZ
  {
    DVector x( convertToVector( xList ) );
 
    DTaylor solver(Fhn_vf, order);

    DCoordinateSection *section [pm_count];
    DPoincareMap *pm [pm_count];

    for( int i = 0; i < pm_count; i++ )    // could have used coordinate sections but doesnt matter in the end
    {
      section[i] = new DCoordinateSection( dim, dim-1, xSectionVCoordVector[i] );
      if(  i < pm_count_up )
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

    cout << "test \n";

    for( int i = 0; i < pm_count-1; i++ )
    {
      monodromyMatrix.setToIdentity(); // probably obsolete since it is done automatically in the solver (?) code
    
      tempvar[0] = x[2*i];
      tempvar[1] = x[2*i + 1];
      tempvar[2] = xSectionVCoordVector[i];

      cout << "Integration to section " << i+1 << " in progress, tempvar = " << tempvar << " \n";
      
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
    tempvar[2] = xSectionVCoordVector[pm_count - 1];
 
    cout << tempvar << "\n";
      cout << xSectionVCoordVector[0] << "\n" ;


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

    std::vector<DVector> result( xList.size() );
 
    for( unsigned int i = 0; i < xList.size(); i++ )
      (result[i]).resize( fast_dim );

    cout << "AFTER INTEGRATION, BEFORE GAUSS \n";
    result = convertToList( x - capd::matrixAlgorithms::gauss( derMatrix, x_eval ) );
    cout << "AFTER GAUSS \n";
    
    for( unsigned int i = 0; i < result.size(); i++ )
      cout << result[i] << "\n";

    return result;
  };


  std::vector<DVector> newtonAlgorithm( double _tolerance )
  {
    setStartListAndVCoordVector();
    double error = 1.;
  
    DVector x0( convertToVector( xStartlist ) );
    DVector x1( convertToVector( xStartlist ) );
    int i = 0;

    while( error > _tolerance )
    {
      cout << i << " " << error << "\n";
      i++;

      x0 = x1;
      x1 = convertToVector( oneNewtonStep( convertToList(x1) ) );

      error = sqrt( scalarProduct(x1 - x0, x1 - x0) );
    }
    
    return convertToList( x1 );
  };


};



int main(){

  cout.precision(9);


  double theta = double(61.)/100.;
  double tolerance = 1e-4;

  Fhn_vf.setParameter( "theta", theta );
  Fhn_vf.setParameter( "eps", eps );

  cout << Fhn_vf(DVector({ -0.105291191426533, 0.00193367348280320, 0.0250009165557845 }) ) << "\n";
  cout << Fhn_vf(DVector({ 0.925215489440428, -0.00210506308805905, 0.0581098505368723 }) ) << "\n";
  cout << Fhn_vf(DVector({ 0.729656580361698, -0.0401973081151965, 0.0980516028684074 }) ) << "\n";


  DVector GammaUL(0.970345591417269, 0., 0.0250442158334208);                                   // some guesses for the corner points which are equilibria
  DVector GammaDL(-0.108412947498862, 0., 0.0250442158334208);                                  // of the fast subsystem for critical parameter v values (third variable)
                                                                                                 // where heteroclinics exist
  DVector GammaUR(0.841746280832201, 0., 0.0988076360184288);                                   // UR up right, DR down right, UL up left, DL down left
  DVector GammaDR(-0.237012258083933, 0., 0.0988076360184288);

 // FhnFindPeriodicOrbit FindPeriodicOrbit( theta, GammaUL, GammaUR, GammaDL, GammaDR );
  
  // FindPeriodicOrbit.newtonAlgorithm( tolerance );


 
  return 0;
} 


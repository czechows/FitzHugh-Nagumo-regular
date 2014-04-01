#include <iostream>
#include <vector>
#include "capd/capdlib.h"

using std::cout;
using namespace capd;
using matrixAlgorithms::inverseMatrix;

DMap Fhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
// FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)
DMap Fhn_fast("par:theta,v;var:u,w;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v);");
// fast subsystem of the FitzHugh-Nagumo system
  
double eps = double(1e-3); 
const int order = 6;
const int precomp_factor = 8; 

#include "matcontPrecomputedOrbit.hpp"

std::vector<DVector> xPrecomputed;

void xPrecomputedFill()
{
  std::vector<DVector> xPrecomputedTemp( precomp_factor * xPrecomputedBeforeCor.size() );

  for( unsigned int i = 0; i < xPrecomputedBeforeCor.size(); i++ )
  {
    double crit_eps = 0.01;
    double crit1 = 0.0250442158334208 + crit_eps;
    double crit2 = 0.0988076360184288 - crit_eps;

    if( (xPrecomputedBeforeCor[i])[2] > crit1 && (xPrecomputedBeforeCor[i])[2] < crit2 )
    {
      for( int j = 0; j < precomp_factor; j++ )
      {
        {
          DVector sum1 = ( double(precomp_factor - j)/double(precomp_factor) )*xPrecomputedBeforeCor[i];
          DVector sum2;

          if( i == xPrecomputedBeforeCor.size() - 1 )
            sum2 = ( double(j)/double(precomp_factor) )*xPrecomputedBeforeCor[0];
          else
            sum2 = ( double(j)/double(precomp_factor) )*xPrecomputedBeforeCor[i+1];

          xPrecomputedTemp[ i*precomp_factor + j ] = sum1 + sum2;
        }
      }
    }
    else
    {
      xPrecomputedTemp[ i*precomp_factor ] = xPrecomputedBeforeCor[i];

      for( int j = 1; j < precomp_factor; j++)
        xPrecomputedTemp[ i*precomp_factor + j ] = DVector({ 0., 0., 0. });
    }
  }

  for( unsigned int i = 0; i < xPrecomputedTemp.size(); i++ )
  {
    if( xPrecomputedTemp[i] != DVector({ 0., 0., 0. }) )
        xPrecomputed.push_back( xPrecomputedTemp[i] );
  }
}


void orthogonalizeRelativeColumn( DMatrix& matrixToOrthogonalize, unsigned int columnNo )
{
  for( unsigned int i = 0; i <= matrixToOrthogonalize.numberOfColumns() - 1; i++ ) 
  { 
    DVector vectorInvariant( matrixToOrthogonalize.column( columnNo ) );
    if( i != columnNo )
    {
      DVector vectorToOrthogonalize( matrixToOrthogonalize.column(i) );
      DVector projection = ( scalarProduct( vectorToOrthogonalize, vectorInvariant )/scalarProduct( vectorInvariant, vectorInvariant ) ) * vectorInvariant;

      for( unsigned int j = 1; j <= matrixToOrthogonalize.numberOfRows(); j++ )
      {
        matrixToOrthogonalize(j,i+1) = vectorToOrthogonalize(j) - projection(j);
      }
    }
  }
}



class FhnFindPeriodicOrbit
{
public:
  int pm_count;
  DMap vectorField;
  DTaylor solver;
  std::vector<DVector> xStartlist;
  std::vector<DMatrix> P_list;
  std::vector<DAffineSection> section;
  std::vector<DPoincareMap> pm;
  int dim;
  int fast_dim;
  double time;

  FhnFindPeriodicOrbit( const int _pm_count )
    : pm_count ( _pm_count ),
      vectorField( Fhn_vf ),
      solver(vectorField, order),
      xStartlist( pm_count ),
      P_list( pm_count ),
      section( pm_count, DAffineSection( xPrecomputed[0], xPrecomputed[1]-xPrecomputed[0]) ),
      dim(3),
      fast_dim(dim - 1),
      time(0.)
  {
    for( unsigned int i = 0; i < xStartlist.size(); i++ )
    {
      (xStartlist[i]).resize( fast_dim );
      
      for( int j = 1; j <= fast_dim; j++ )
        (xStartlist[i])(j) = 0.;
    }

    for( int i = 0; i < pm_count - 1; i++ )
    {
      (P_list[i]).resize( dim, dim );
      (P_list[i]).setToIdentity();

      for( int j = 1; j <= dim; j++ )
       (P_list[i])( j, dim ) = ( xPrecomputed[i+1]-xPrecomputed[i] )(j);

      orthogonalizeRelativeColumn( P_list[i], dim - 1 );

      section[i] = DAffineSection( xPrecomputed[i], xPrecomputed[i+1]-xPrecomputed[i] );
    }
 
    (P_list[ pm_count - 1 ]).resize( dim, dim );
    (P_list[ pm_count - 1 ]).setToIdentity();

    for( int j = 1; j <= dim; j++ )
      (P_list[ pm_count - 1 ])( j, dim ) = ( xPrecomputed[0]-xPrecomputed[pm_count-1] )(j);
 
    section[ pm_count - 1 ] = DAffineSection( xPrecomputed[pm_count - 1], xPrecomputed[0]-xPrecomputed[ pm_count - 1 ] );

    orthogonalizeRelativeColumn( P_list[pm_count-1], dim - 1 );

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


  std::vector<DVector> oneNewtonStep( std::vector<DVector> xList, double& error )   // one Newton step for function xi - P_{i}(x_{i-1}), i in mathbbZ
  {
    DVector x( convertToVector( xList ) );
 
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
  
      tempvar = (section[i]).getOrigin() + (P_list[i])*DVector({ x[2*i], x[2*i + 1], 0 });

 //     cout << "Integration to section " << i+1 << " in progress, tempvar = " << tempvar << " \n";
      
      DPoincareMap pm( solver, section[i+1], poincare::MinusPlus );

      pm_result = pm(tempvar, monodromyMatrix, time);
      poincareDer = pm.computeDP(pm_result, monodromyMatrix, time);
      
 //     cout << pm_result << "\n";

      pm_result = inverseMatrix( P_list[i+1] )*( pm_result - (section[i+1]).getOrigin() );
      poincareDer = inverseMatrix( P_list[i+1] )*poincareDer*(P_list[i]);

      x_eval[2*i + 2] = x[2*i + 2] - pm_result[0];
      x_eval[2*i + 3] = x[2*i + 3] - pm_result[1];

      derMatrix[2*i+2][2*i] = -poincareDer[0][0];
      derMatrix[2*i+2][2*i+1] = -poincareDer[0][1];
      derMatrix[2*i+3][2*i] = -poincareDer[1][0]; 
      derMatrix[2*i+3][2*i+1] = -poincareDer[1][1];
    }


    monodromyMatrix.setToIdentity();
  
    tempvar = (section[pm_count - 1]).getOrigin() + (P_list[pm_count - 1])*DVector({ x[2*pm_count - 2], x[2*pm_count - 1], 0 });

    DPoincareMap pm( solver, section[0], poincare::MinusPlus );

    pm_result = pm( tempvar, monodromyMatrix, time );
    poincareDer = pm.computeDP( pm_result, monodromyMatrix, time);
 
    pm_result = inverseMatrix( P_list[0] )*( pm_result - (section[0]).getOrigin() );
    poincareDer = inverseMatrix( P_list[0] )*poincareDer*(P_list[pm_count - 1]);

    x_eval[0] = x[0] - pm_result[0];
    x_eval[1] = x[1] - pm_result[1];

    derMatrix[0][2*pm_count - 2] = -poincareDer[0][0];
    derMatrix[0][2*pm_count - 1] = -poincareDer[0][1];
    derMatrix[1][2*pm_count - 2] = -poincareDer[1][0];
    derMatrix[1][2*pm_count - 1] = -poincareDer[1][1];

    error = scalarProduct( x_eval, x_eval );

    std::vector<DVector> result( xList.size() );
 
    for( unsigned int i = 0; i < xList.size(); i++ )
      (result[i]).resize( fast_dim );

 //   cout << "AFTER INTEGRATION, BEFORE GAUSS \n";
    result = convertToList( x - capd::matrixAlgorithms::gauss( derMatrix, x_eval ) );
 //   cout << "AFTER GAUSS \n";
    
    return result;
  };


  std::vector<DVector> newtonAlgorithm( double _tolerance )
  {
    double error = 1.;
  
    DVector x1( convertToVector( xStartlist ) );

    while( error > _tolerance )
    {
      x1 = convertToVector( oneNewtonStep( convertToList(x1), error ) );
      cout << error << "\n";
    }
    
    return convertToList( x1 );
  };


  std::vector<DVector> returnCorrectedOrbit( double _tolerance )
  {
    std::vector<DVector> xList = newtonAlgorithm( _tolerance );

    std::vector<DVector> result;

    for( int i = 0; i < pm_count; i++ )
      result.push_back( (section[i]).getOrigin() + (P_list[i]) * DVector({ (xList[i])[0], (xList[i])[1], 0. }) );

    return result;
  }

  void integrateOneTime( double _tolerance )
  {
    DVector vect = ( returnCorrectedOrbit( _tolerance ) )[0];
    cout << vect << "\n";
    
    DPoincareMap pm( solver, section[0], poincare::MinusPlus );
    double temp_time = 0;

    cout << pm( vect, temp_time ) << " " << temp_time << "\n";
  }

};






int main(){

  cout.precision(15);


  double theta = double(61.)/100.;
  double tolerance = 1e-28;

  Fhn_vf.setParameter( "theta", theta );
  Fhn_vf.setParameter( "eps", eps );

  xPrecomputedFill();

  const int pm_count = xPrecomputed.size();

  FhnFindPeriodicOrbit FindPeriodicOrbit( pm_count );
  
  std::vector<DVector> result = FindPeriodicOrbit.returnCorrectedOrbit( tolerance );
 

  for( unsigned int i = 0; i < xPrecomputed.size(); i++ )
    cout << xPrecomputed[i] << " " << result[i] << "\n";

  FindPeriodicOrbit.integrateOneTime( tolerance );

  return 0;
} 


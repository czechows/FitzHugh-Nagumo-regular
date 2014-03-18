
/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing numerical simulations necessary
 * for the computer assistance of proofs of periodic orbits in the FitzHugh-Nagumo system.
 * We find values of v for which there exist heteroclinic connections between equilibria
 * in the fast subsystem via shooting from stable/unstable manifolds of equilibria
 * to the Poincare section given in between.
 * The validity of the proofs does not depend on validity of numerics given below.
 * If the full FitzHugh-Nagumo vector field is alternated in any other way than changing
 * parameters theta or eps, one needs to manually readjust the fast vector field
 * and the fast reversed vector field given in the class FhnBifurcation
 * ----------------------------------------------------------------------------------------*/


/* ----------------------------------------------------------------------------------------- */
/* ---------------------------- FAST SUBSYSTEM NUMERICS ------------------------------------ */
/* ----------------------------------------------------------------------------------------- */

class FhnBifurcation
{
public:
  DMap vectorField;
  DMap vectorFieldRev;    // reversed (minus) vector field for backward integration
  DTaylor solver;
  DTaylor solverRev;
  DVector EqU;            // "upper" equilibrium (guess)
  DVector EqD;            // "lower" equilibrium (guess)
  double DISP;             // displacement from the equilibra in stable/unstable direction
  DAffineSection section;
  DPoincareMap pm;
  DPoincareMap pmRev;     // reversed Poincare map for backward integration
  bool dir;               // direction - do we go from EqU to EqD or other way round?
  FhnBifurcation (int order, double& _theta, const DVector& _EqU, const DVector& _EqD, double _DISP, bool _dir = 1) 
    : vectorField("par:theta,v;var:u,w;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v);"), // vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v is parameter
      vectorFieldRev("par:theta,v;var:u,w;fun:-w,(-2/10)*(theta*w+u*(u-1)*(u-(1/10))+v);"), //  minus vector field for reverse integration
      solver(vectorField,order),
      solverRev(vectorFieldRev,order),
      EqU(_EqU),
      EqD(_EqD),
      DISP(_DISP),
      section( DVector({0.2, 0. }), DVector({-1.,0.})), // an arbitrary choice of coordinate section here
      pm(solver,section),
      pmRev(solverRev,section),
      dir(_dir)
  {
    vectorField.setParameter("theta",_theta);
    vectorFieldRev.setParameter("theta",_theta);
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

  DMatrix J_correct(DVector Eq, double v)
  {
    vectorField.setParameter("v", v); 
    return( vectorField[Eq] );
  }

  double w_function(DVector guessEqU, DVector guessEqD, double v)  // returns distance in w variable on the v-poincare section between integrated displacement in unstable direction from EqU
                                                                   // and integrated displacement in stable direction from EqD if dir = 1 or vice-versa elsewise
  {
    double return_time = 1.;
    vectorField.setParameter("v", v); 
    vectorFieldRev.setParameter("v", v);

    DMatrix EigenvectU_real(2,2),
            EigenvectU_im(2,2);
    DVector EigenvalU_real(2),
            EigenvalU_im(2);
  
    DMatrix EigenvectD_real(2,2),
            EigenvectD_im(2,2);
    DVector EigenvalD_real(2),
            EigenvalD_im(2);   
    
    computeEigenvaluesAndEigenvectors( J_correct(guessEqU,v), EigenvalU_real, EigenvalU_im, EigenvectU_real, EigenvectU_im );
    computeEigenvaluesAndEigenvectors( J_correct(guessEqD,v), EigenvalD_real, EigenvalD_im, EigenvectD_real, EigenvectD_im );

    if(EigenvalU_real[0]*EigenvalU_real[1] >= 0.)
      throw "EIGENVALUES OF FAST SUBSYSTEM AT STATIONARY POINTS NOT OF OPPOSITE SIGNS! \n";
    if(EigenvalD_real[0]*EigenvalD_real[1] >= 0.)
      throw "EIGENVALUES OF FAST SUBSYSTEM AT STATIONARY POINTS NOT OF OPPOSITE SIGNS! \n";
    if(EigenvalU_im[0] != 0.)
      throw "EIGENVALUES OF FAST SUBSYSTEM AT STATIONARY POINTS NOT REAL! \n";
    if(EigenvalD_im[0] != 0.)
      throw "EIGENVALUES OF FAST SUBSYSTEM AT STATIONARY POINTS NOT REAL! \n";
    
    if(dir == 1) // here one has to play manually with pluses and minuses so we are on the right side of stable/unstable manifolds and we catch the right eigenvectors
        return pm(guessEqU + EigenvectU_real.column(1)*DISP, return_time)[1] - pmRev(guessEqD - EigenvectD_real.column(0)*DISP, return_time)[1]; 
    else
        return pmRev(guessEqU + EigenvectU_real.column(0)*DISP, return_time)[1] - pm(guessEqD - EigenvectD_real.column(1)*DISP, return_time)[1];   
  }

  double v_correct(double v) // secant method to correct v to the bifurcation point, as side effect corrects equilibria EqU and EqD to right positions 
  {
    double error = 1.;

    double v0 = v + 1e-4;
    double v1 = v;
    double v_temp = v;

    DVector EqU0 = Eq_correct(EqU, v0);
    DVector EqU1 = Eq_correct(EqU, v1);
    DVector EqD0 = Eq_correct(EqD, v0);
    DVector EqD1 = Eq_correct(EqD, v1);

   while(error > accuracy)
   {
    v_temp = v1;
    v1 = v1 - w_function(EqU1, EqD1, v1)*( (v1-v0) / ( w_function(EqU1, EqD1, v1) - w_function(EqU0, EqD0, v0) ) );
    EqU0 = Eq_correct( EqU0, v_temp );
    EqD0 = Eq_correct( EqD0, v_temp );
    EqU1 = Eq_correct( EqU1, v1 );
    EqD1 = Eq_correct( EqD1, v1 );
    v0 = v_temp;

    error = abs( w_function(EqU1, EqD1, v1) );
   }

   EqU = EqU1;
   EqD = EqD1;
    
   return v1;
  }
};


void GammaQuad_correct( const double& _theta, DVector& _GammaUL, DVector& _GammaDL, DVector& _GammaUR, DVector& _GammaDR ) // corrects original guesses of Gammas for given theta
{
  double theta( _theta.leftBound() );
  double DISP(1e-12);

  double vR( _GammaUR[2].leftBound() );
  double vL( _GammaUL[2].leftBound() );

  DVector EqUL(2),
          EqUR(2),
          EqDR(2),
          EqDL(2);      // equilibria are first two coordinates of each Gamma, they will be corrected by Newton inside the program

  EqUL[0] = _GammaUL[0].leftBound();
  EqUL[1] = _GammaUL[1].leftBound();
  EqUR[0] = _GammaUR[0].leftBound();
  EqUR[1] = _GammaUR[1].leftBound();

  EqDR[0] = _GammaDR[0].leftBound();
  EqDR[1] = _GammaDR[1].leftBound();
  EqDL[0] = _GammaDL[0].leftBound();
  EqDL[1] = _GammaDL[1].leftBound();


  FhnBifurcation BifR(order, theta, EqUR, EqDR, DISP);
  FhnBifurcation BifL(order, theta, EqUL, EqDL, DISP, 0); 

  double vR_c( BifR.v_correct(vR) );        // corrected v values & equlibria coordinates
  double vL_c( BifL.v_correct(vL) );

  DVector EqUL_c( BifL.EqU ),
          EqUR_c( BifR.EqU ),
          EqDR_c( BifR.EqD ), 
          EqDL_c( BifL.EqD );

  _GammaUL[0] = EqUL_c[0];                  // assignment back to Gammas
  _GammaUL[1] = EqUL_c[1];
 
  _GammaUR[0] = EqUR_c[0];
  _GammaUR[1] = EqUR_c[1];
 
  _GammaDL[0] = EqDL_c[0];
  _GammaDL[1] = EqDL_c[1];
 
  _GammaDR[0] = EqDR_c[0];
  _GammaDR[1] = EqDR_c[1];

  _GammaUL[2] = _GammaDL[2] = vL_c;
  _GammaUR[2] = _GammaDR[2] = vR_c;
}


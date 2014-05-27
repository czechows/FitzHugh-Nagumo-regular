/* -----------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing parameter continuation
 * of the computer assisted proof
 * of existence of periodic solutions
 * of the FitzHugh-Nagumo equations
 * ------------------------------------------ */

class FhnValidatedContinuation
{
  public:
  interval theta;
  interval epsRange;
  bool isFirstTry;
  bool epsIncreasing;
  double tolerance;
  double radius;
  std::vector<DVector> numericOrbitGuess;
  interval increment;
  double incrementFactor;
  interval currentEpsRange;
  interval integrationTimeBound;

  FhnValidatedContinuation( interval _theta, interval _eps, std::vector<DVector> _precomputedOrbit, bool _epsIncreasing, double _tolerance, double _radius )
    : theta( _theta ),
      epsRange( _eps ),
      isFirstTry( 1 ),
      epsIncreasing( _epsIncreasing ),
      tolerance( _tolerance ),
      radius( _radius ),
      numericOrbitGuess( _precomputedOrbit ),
      increment( 0., 1e-7 ),
      incrementFactor( 1.05 ),
      currentEpsRange( epsRange ),
      integrationTimeBound( 0.2, 4. )
  {
    if( epsIncreasing )
      currentEpsRange = left( epsRange ) + increment;
    else
      currentEpsRange = right( epsRange ) - increment;
  }

  void moveCurrentEpsRange()
  {
    if( epsIncreasing )
      currentEpsRange = right( currentEpsRange ) + increment;
    else
      currentEpsRange = left( currentEpsRange ) - increment;
  }

  void decreaseCurrentEpsRange()
  {
    increment = increment / incrementFactor;

    if( epsIncreasing )
      currentEpsRange = left( currentEpsRange ) + increment;
    else
      currentEpsRange = right( currentEpsRange ) - increment;
  }

  void tryToProveOrbit( FhnCovering *_cov )
  {
    while(true)
    {
      try
      {
        (*_cov).proveExistenceOfOrbit( theta, currentEpsRange, tolerance, radius ); 
        break;
      }
      catch(const char* Message)
      {
        cout << Message << "\n";
        isFirstTry = 0;
        decreaseCurrentEpsRange();
      }
      catch(std::domain_error)
      {
        cout << "DOMAIN ERROR! \n";
        isFirstTry = 0;
        decreaseCurrentEpsRange();
      }
    }
  }


  void continueOrbitWithEps()
  {
    while( ( !epsIncreasing && currentEpsRange.leftBound() >= epsRange.leftBound() ) ||  ( epsIncreasing && currentEpsRange.rightBound() <= epsRange.rightBound() ) )
    {
      moveCurrentEpsRange();

      FhnCovering *cov = new FhnCovering( numericOrbitGuess );
      tryToProveOrbit( cov );
      numericOrbitGuess = (*cov).getCorrectedGuess( integrationTimeBound );
      delete cov;

      if( isFirstTry )
        increment = increment * incrementFactor;


      cout << "Existence of periodic solution for parameter values eps = " << currentEpsRange << " and theta = " << theta << " proven. \n Increment size: " << increment.rightBound() << "\n";
      cout.flush();
    }

    cout << "\n EXISTENCE OF PERIODIC SOLUTION FOR PARAMETER VALUES EPS = " << epsRange << " AND THETA = " << theta << "PROVEN !. \n";
    
  }

};

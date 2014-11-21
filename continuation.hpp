/* -----------------------------------------
 * This is a header file for continuation_fhn.cpp
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
  int successCount;

  FhnValidatedContinuation( interval _theta, interval _eps, std::vector<DVector> _precomputedOrbit, bool _epsIncreasing, double _tolerance, double _radius, 
      double _startIncrementSize, interval _integrationTimeBound )
    : theta( _theta ),
      epsRange( _eps ),
      isFirstTry( 1 ),
      epsIncreasing( _epsIncreasing ),
      tolerance( _tolerance ),
      radius( _radius ),
      numericOrbitGuess( _precomputedOrbit ),
      increment( 0., _startIncrementSize ),
      incrementFactor( 1.05 ),
      currentEpsRange( epsRange ),
      integrationTimeBound( _integrationTimeBound ),
      successCount(0)
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

  void decreaseRadius()
  {
    radius = radius / incrementFactor;
  }

  void tryToProveOrbit( FhnCovering *_cov )
  {
    bool isFirstTryForCurrentEps(1);
    while(true)
    {
      try
      {
        (*_cov).proveExistenceOfOrbit( theta, currentEpsRange, tolerance, radius ); 
        break;
      }
      catch(const char* Message) // if the proof fails we try to decrease the epsilon range and h-sets size
      {
        cout << Message << "\n";
        isFirstTry = 0;
        isFirstTryForCurrentEps = 0;
        successCount = 0;
        decreaseCurrentEpsRange();
        decreaseCurrentEpsRange();
        decreaseRadius();
      }
      catch( capd::poincare::PoincareException< C0Rect2Set > )
      {
        cout << "POINCARE EXCEPTION! \n";
        isFirstTry = 0;
        isFirstTryForCurrentEps = 0;
        decreaseCurrentEpsRange();
        decreaseCurrentEpsRange();
        decreaseRadius();
      }
    }
    successCount = successCount + int( isFirstTryForCurrentEps );
  }


  void saveNewEpsRange() // we save the Eps range in binary (! exact value) so we can resume the program later for the saved range
  {
    interval savedEpsInterval;

    if( epsIncreasing )
      savedEpsInterval = interval( currentEpsRange.leftBound(), epsRange.rightBound() );
    else
      savedEpsInterval = interval( epsRange.leftBound(), currentEpsRange.rightBound() );

    std::ofstream savedEpsRange;
    savedEpsRange.open("savedEpsRange", std::ios::trunc | std::ios::binary );

    hexWrite( savedEpsRange, savedEpsInterval );
    savedEpsRange.close();

   /* savedEpsRange.precision(16);
    if( epsIncreasing )
      savedEpsRange << "interval savedEpsRange = interval( " << currentEpsRange.leftBound() << ", " << epsRange.rightBound() << ");";
    else
      savedEpsRange << "interval savedEpsRange = interval( " << epsRange.leftBound() << ", " << currentEpsRange.rightBound() << ");";*/
  }

  void continueOrbitWithEps()
  {
    interval oldEpsRange( currentEpsRange );

    while( ( !epsIncreasing && oldEpsRange.leftBound() >= epsRange.leftBound() ) ||  ( epsIncreasing && oldEpsRange.rightBound() <= epsRange.rightBound() ) )
    {
      FhnCovering *cov = new FhnCovering( numericOrbitGuess );
     // (*cov).initialize( theta, currentEpsRange, tolerance, radius );
      tryToProveOrbit( cov );

      numericOrbitGuess = (*cov).getCorrectedGuess( integrationTimeBound );

      delete cov;

      cout << "Existence of a periodic solution for parameter values eps = " << currentEpsRange << " and theta = " << theta << " proven. \nIncrement size: " << increment.rightBound() << "\n";
      cout.flush();

      oldEpsRange = currentEpsRange;
      moveCurrentEpsRange();
      saveNewEpsRange();

      if( isFirstTry || successCount > 100 )
      {
        increment = increment * incrementFactor;
        radius = radius * sqrt( incrementFactor );
        successCount = 0;
      }
    }

    cout << "\n EXISTENCE OF PERIODIC SOLUTION FOR PARAMETER VALUES EPS = " << epsRange << " AND THETA = " << theta << "PROVEN !. \n";
    
  }

};

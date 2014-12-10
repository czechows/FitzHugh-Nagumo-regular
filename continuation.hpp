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
  bool usingNewton;

  FhnValidatedContinuation( interval _theta, interval _eps, std::vector<DVector> _precomputedOrbit, bool _epsIncreasing, double _tolerance, double _radius, 
      double _startIncrementSize, interval _integrationTimeBound, bool _usingNewton=0 )
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
      successCount(0),
      usingNewton( _usingNewton )
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
 
  void tryToProveOrbit( FhnIntervalNewton *_newt )
  {
    bool isFirstTryForCurrentEps(1);
    while(true)
    {
      try
      {
        (*_newt).proveExistenceOfOrbitWithNewton( theta, currentEpsRange, tolerance, radius ); 
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
      catch( capd::poincare::PoincareException< C1Rect2Set > )
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
  }

  void continueOrbitWithEps()
  {
    if( usingNewton )
    {
      FhnFindPeriodicOrbit *num = new FhnFindPeriodicOrbit( numericOrbitGuess );   // one extra nonrigorous run before the continuation to adjust the number of sections
      (*num).setDParameters( theta, currentEpsRange );
      (*num).newtonAlgorithm( tolerance );
      numericOrbitGuess = (*num).getCorrectedGuess( integrationTimeBound );
      delete num; 
    }

    do
    {
      interval totalOrbitPeriod;

      if( usingNewton )
      {
        FhnIntervalNewton *method = new FhnIntervalNewton( numericOrbitGuess );
        tryToProveOrbit( method );
        numericOrbitGuess = (*method).getCorrectedGuess( integrationTimeBound );
        totalOrbitPeriod = (*method).totalPeriod;
        delete method;
      }
      else
      {
        FhnCovering *method = new FhnCovering( numericOrbitGuess );
        tryToProveOrbit( method );
        numericOrbitGuess = (*method).getCorrectedGuess( integrationTimeBound );
        totalOrbitPeriod = (*method).totalPeriod;
        delete method;
      }
    
      cout << "Existence of a periodic solution for parameter values eps = " << currentEpsRange << " and theta = " << theta << " proved. \nIncrement size: " << increment.rightBound() << "\n";
      cout << "Radius: " << radius << "\n";
      cout << "Bound for total period: " << totalOrbitPeriod << "\n";

      cout.flush();

      moveCurrentEpsRange();
      saveNewEpsRange();

     if( isFirstTry || successCount > 10 )
     {
      increment = increment * incrementFactor;
      radius = radius * sqrt( incrementFactor );
      successCount = 0;
     }
    }
    while( ( !epsIncreasing && currentEpsRange.rightBound() >= epsRange.leftBound() ) ||  ( epsIncreasing && currentEpsRange.leftBound() <= epsRange.rightBound() ) );

    cout << "\nEXISTENCE OF A PERIODIC SOLUTION FOR PARAMETER VALUES EPS = " << epsRange << " AND THETA = " << theta << " PROVED!. \n \n";
  }

};

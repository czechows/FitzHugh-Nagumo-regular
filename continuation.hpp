/* -----------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing parameter continuation
 * of the computer assisted proof
 * of existence of periodic solutions
 * of the FitzHugh-Nagumo equations
 * ------------------------------------------ */


void continueOrbitWithEps( interval theta, interval epsRange, bool epsIncreasing, std::vector<DVector> precomputedOrbit, 
    double tolerance, double radius, double startIntervalSize = 1e-5, double incrementFactor = 1.03  )
{

  std::vector<DVector> numericOrbitGuess( precomputedOrbit );
  double incrementSize( startIntervalSize );
  bool isFirstTry(1);

  if( epsIncreasing )
  {

    interval currentEpsRangeLeftBound( epsRange.leftBound() );
    interval currentEpsRangeRightBound( currentEpsRangeLeftBound );
 
    while( currentEpsRangeRightBound.leftBound() <= epsRange.rightBound() )
    {
      currentEpsRangeLeftBound = currentEpsRangeRightBound;
      currentEpsRangeRightBound = currentEpsRangeLeftBound + incrementSize;
      interval currentEpsRange( currentEpsRangeLeftBound.leftBound(), currentEpsRangeRightBound.rightBound() );
      
      FhnCovering cov( numericOrbitGuess );
      bool proofResult(0);
      
      try
      {
        bool covResult( cov.proveExistenceOfOrbit( theta, currentEpsRange, tolerance, radius ) ); 
        proofResult = covResult;
      }
      catch(const char* Message)
      {
        cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES THETA=" << theta << " AND EPS=" << currentEpsRange << " NOT PROVEN! \n";
        proofResult = 0;
      }
 
      while( !proofResult )
      {
        isFirstTry = 0;
        cout << "EXISTENCE OF PERIODIC SOLUTION FOR PARAMETER VALUES EPS = " << currentEpsRange << " AND THETA = " << theta << "NOT PROVEN. \n Increment size: " << incrementSize << "\n";
        incrementSize = incrementSize / incrementFactor;

        currentEpsRangeRightBound = currentEpsRangeLeftBound + incrementSize;
        currentEpsRange = interval( currentEpsRangeLeftBound.leftBound(), currentEpsRangeRightBound.rightBound() );

        try
        {
          bool covResult( cov.proveExistenceOfOrbit( theta, currentEpsRange, tolerance, radius ) ); 
          proofResult = covResult;
        }
        catch(const char* Message)
        {
          cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES THETA=" << theta << " AND EPS=" << currentEpsRange << " NOT PROVEN! \n";
          proofResult = 0;
        }
   
      }

      numericOrbitGuess = cov.getCorrectedGuess();
      
    // if( isFirstTry )
       incrementSize = incrementSize * incrementFactor;

      cout << "Existence of periodic solution for parameter values eps = " << currentEpsRange << " and theta = " << theta << " proven. \n Increment size: " << incrementSize << "\n";
    }

    cout << "\n EXISTENCE OF PERIODIC SOLUTION FOR PARAMETER VALUES EPS = " << epsRange << " AND THETA = " << theta << "PROVEN !. \n";
  }
  else
  {
    interval currentEpsRangeRightBound( epsRange.rightBound() );
    interval currentEpsRangeLeftBound( currentEpsRangeRightBound );

    while( currentEpsRangeLeftBound.rightBound() >= epsRange.leftBound() )
    {
      currentEpsRangeRightBound = currentEpsRangeLeftBound;
      currentEpsRangeLeftBound = currentEpsRangeRightBound - incrementSize;
      interval currentEpsRange( currentEpsRangeLeftBound.leftBound(), currentEpsRangeRightBound.rightBound() );
      
      FhnCovering cov( numericOrbitGuess );

      bool proofResult( 0 );

      try
      {
        bool covResult( cov.proveExistenceOfOrbit( theta, currentEpsRange, tolerance, radius ) ); 
        proofResult = covResult;
      }
      catch(const char* Message)
      {
        cout << Message << "\n";
        proofResult = 0;
      }
 
      while( !proofResult )
      {
        isFirstTry = 0;
        cout << "EXISTENCE OF PERIODIC SOLUTION FOR PARAMETER VALUES EPS = " << currentEpsRange << " AND THETA = " << theta << "NOT PROVEN. \n Increment size: " << incrementSize << "\n";
        incrementSize = incrementSize / incrementFactor;

        currentEpsRangeLeftBound = currentEpsRangeRightBound - incrementSize;
        currentEpsRange = interval( currentEpsRangeLeftBound.leftBound(), currentEpsRangeRightBound.rightBound() );

        try
        {
          bool covResult( cov.proveExistenceOfOrbit( theta, currentEpsRange, tolerance, radius ) ); 
          proofResult = covResult;
        }
        catch(const char* Message)
        {
          cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES THETA=" << theta << " AND EPS=" << currentEpsRange << " NOT PROVEN! \n";
          proofResult = 0;
        }
       }

      numericOrbitGuess = cov.getCorrectedGuess();
      
      if( isFirstTry && incrementSize < 1e-6 )
        incrementSize = incrementSize * incrementFactor;

      cout << "Existence of periodic solution for parameter values eps = " << currentEpsRange << " and theta = " << theta << " proven. \n Increment size: " << incrementSize << "\n";
    }

    cout << "\n EXISTENCE OF PERIODIC SOLUTION FOR PARAMETER VALUES EPS = " << epsRange << " AND THETA = " << theta << "PROVEN !. \n";
  }
}

/* -------------------------------------------------
 * This is a header file for nonrigorousNewton.hpp
 * implementing diagonalization of coordinates 
 * for optimization of the Newton algorithm 
 * * ------------------------------------------------- */


  DMatrix setCoordSystemD( int coord_no, DVector _disp )            // similar to coordChange in fhn program, i < pm_count - 1 !       
  {
    DVector Gamma = section[ coord_no ].getOrigin() + P_list[ coord_no ]*DVector( _disp(1), _disp(2), 0. );

    int vdim( 3 );   // should be used only in dimension 3!
    DMatrix JacobianD( vdim, vdim );

    // a patch to set eps to 0 for the vector field for computing coordinates around slow manifold
    DMap vectorFieldZeroEps( *vectorField );          
    vectorFieldZeroEps.setParameter("eps", 0.);
    // WARNING: specific to the (type of) the vector field 

    for(int i=0; i<vdim; i++)              
    {
      for(int j=0; j<vdim; j++)
        JacobianD[i][j] = ( vectorFieldZeroEps[Gamma] )[i][j];
    }

    DVector tempvectRe( vdim );   
    DVector tempvectIm( vdim );

    DMatrix tempmatrixIm( vdim, vdim );
    
    DMatrix P( vdim, vdim );
    DMatrix P_result( vdim, vdim );

    computeEigenvaluesAndEigenvectors(JacobianD, tempvectRe, tempvectIm, P, tempmatrixIm);

    int i_min(42);  // i's initialized with any value
    int i_med(42);
    int i_max(42); 

    for( int i = 1; i <= vdim; i++ ) // sorting so we are sure we have stable coordinate first unstable second neutral third. ONLY FOR 3D vector fields!
    {
      if( tempvectRe(i) == min( tempvectRe(1), min(tempvectRe(2), tempvectRe(3)) ) )
        i_min = i;
      else if( tempvectRe(i) == max( tempvectRe(1), max(tempvectRe(2), tempvectRe(3)) ) )
        i_max = i;
      else 
        i_med = i;
    }

    for( int i = 1; i <=vdim; i++ )
    { 
      P_result(i,1) = P( i, i_min );
      P_result(i,2) = P( i, i_max );
      P_result(i,3) = P( i, i_med );
    }

    return P_result; 
  };


      // SETTING GOOD COORDINATES FOR NONRIGOROUS NEWTON - DOES NOT HELP THAT MUCH
  void initD( double _tolerance )
  {
    int breakPoint( pm_count + 1 );
    DMatrix P_breakPoint( dim, dim );

    for( int i = 0; i < pm_count; i++ )
    {
      if( ( (section[i].getOrigin() )(dim) > 0.05 ) && ( (section[i].getOrigin() )(dim)  < 0.06 ) ) 
      {
        P_breakPoint = setCoordSystemD( i, section[i].getOrigin() );
        breakPoint = i;
        break;
      }
    }

    if( breakPoint == pm_count + 1 )
      throw "COULD NOT INITIALIZE COORDINATE SYSTEMS AROUND APPROXIMATE ORBIT. NOT ENOUGH SUBSECTIONS";
  
    DVector unstableVect( P_breakPoint.column(1) );
    DVector stableVect( P_breakPoint.column(0) );

    double error(1.);

    while( error > _tolerance )
    {
      DVector unstableVectCp( unstableVect ),
              stableVectCp( stableVect );

      error = vectorNorm( unstableVectCp - unstableVect ) + vectorNorm( stableVectCp - stableVect );

      for( int i = 0; i < pm_count; i++ )
      {
        DMatrix uMonodromyMatrix( dim, dim ),
                sMonodromyMatrix( dim, dim );

        double uReturnTime(0.),
               sReturnTime(0.);

        int uSecCount( (breakPoint + i) % pm_count ),
            sSecCount( ( (breakPoint - i) + pm_count ) % pm_count );

        DVector uSetToIntegrate( correctedGuess[ uSecCount ] ),
                sSetToIntegrate( correctedGuess[ sSecCount ] );

        DMatrix uDerivativeOfPm( dim, dim ),
                sDerivativeOfPm( dim, dim );

        int local_order( order );

      //  double maxStep(4.);

     //   while(true)
     //   {
     //     try
     //     {
        {
          DTaylor solver( *vectorField, local_order );
   //       solver.setMaxStep( maxStep );

          DTaylor solverRev( *vectorFieldRev, local_order );
    //      solverRev.setMaxStep( maxStep );

          DPoincareMap uPM( solver, section[ (uSecCount + 1) % pm_count ] ), 
                       sPM( solverRev, section[ ((sSecCount - 1) + pm_count) % pm_count ] );

          DVector uResult( uPM( uSetToIntegrate, uMonodromyMatrix, uReturnTime ) ),
                  sResult( sPM( sSetToIntegrate, sMonodromyMatrix, sReturnTime ) );

          uDerivativeOfPm =  uPM.computeDP( uResult, uMonodromyMatrix, uReturnTime ) ,
          sDerivativeOfPm = sPM.computeDP( sResult, sMonodromyMatrix, sReturnTime ) ;
        }

      //     break;
                                    // this catch would decrease the order or maxStep so nonrigorous integration would not cross a section in one step 
                                    // (then integrator throws domain_error); for the parameter range in the paper not necessary and did not work that well
                                    // so we disable it, however we leave it to uncomment for future purposes
     //   }
      //    catch(std::domain_error)
       //   {
       //     cout << "DOMAIN ERROR! \n";
         //   local_order = local_order - 1;
      //      maxStep = maxStep / 2.;
         //   cout << "ORDER DECREASED!";
       //     cout << "MAX STEP DECREASED;";
      //    }
      //  }
 

        unstableVect = uDerivativeOfPm * unstableVect;
        stableVect = sDerivativeOfPm * stableVect;

        unstableVect = unstableVect / vectorNorm( unstableVect );
        stableVect = stableVect / vectorNorm( stableVect );
        
        for( int j = 1; j <= dim; j++ )
        {
          (P_list[ uSecCount ] )( j, 2 ) = unstableVect( j );
          (P_list[ sSecCount ])( j, 1 ) = stableVect( j );
        }
      }

    }
  }; 
 


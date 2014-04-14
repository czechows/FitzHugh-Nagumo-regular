/* -----------------------------------------
 * This is a header file for integration_fhn.cpp
 * implementing auxiliaries for Newton algorithms
 * that is condensing the precomputed list of points
 * along the orbit and orthogonalizing matrices
 * to match the AffineSection implementation
 * (origin + normal)
 * ------------------------------------------ */




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

void IxPrecomputedFill()
{
  IxPrecomputed.resize( xPrecomputed.size(), IVector({ 0.,0.,0. }) );
  for( unsigned int i = 0; i < xPrecomputed.size(); i++ )
  {
    for( unsigned int j = 1; j <= ( xPrecomputed[i] ).dimension(); j++ )
      (IxPrecomputed[i])(j) = (xPrecomputed[i])(j);
  }
}


/* ------------------------------------------------------------------------------------------------------------------ */
/* ---------------------------- Matrix orthogonalization algorithms (should be done in a template) ------------------ */
/* ------------------------------------------------------------------------------------------------------------------ */


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

void IOrthogonalizeRelativeColumn( IMatrix& matrixToOrthogonalize, unsigned int columnNo )
{
  for( unsigned int i = 0; i <= matrixToOrthogonalize.numberOfColumns() - 1; i++ ) 
  { 
    IVector vectorInvariant( matrixToOrthogonalize.column( columnNo ) );
    if( i != columnNo )
    {
      IVector vectorToOrthogonalize( matrixToOrthogonalize.column(i) );
      IVector projection = ( scalarProduct( vectorToOrthogonalize, vectorInvariant )/scalarProduct( vectorInvariant, vectorInvariant ) ) * vectorInvariant;

      for( unsigned int j = 1; j <= matrixToOrthogonalize.numberOfRows(); j++ )
      {
        matrixToOrthogonalize(j,i+1) = vectorToOrthogonalize(j) - projection(j);
      }
    }
  }
}



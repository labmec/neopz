/**
 * @file pzvec_extra.h
 * @brief Extra utilities for TPZVec.
 */
// $Id: pzvec_extras.h,v 1.1.1.1 2003-02-04 16:45:27 cantao Exp $

#ifndef PZVEC_EXTRAS_H
#define PZVEC_EXTRAS_H

/**
 * Performs a saxpy operation: x <- x + s * y.
 *
 * @since Jan 9, 2002
 * @author Cantao!
 */
template< class T1, class T2, class Scalar >
void saxpy( TPZVec< T1 >& x, const TPZVec< T2 >& y, Scalar s )
{
   int size = x.NElements();

#ifndef NOTDEBUG
   if( size != y.NElements() ) {
      PZError << "SAXPY error!" << endl
	      << "Vectors with different sizes #x = " << size
	      << ", #y = " << y.NElements() << endl;

      PZError.flush();
   }
#endif

   for( int ii = 0; ii < size; ii++ )
   {
      x[ ii ] += s * y[ ii ];
   }
}

#endif //PZVEC_EXTRAS_H

//--| PZ |----------------------------------------------------------------------

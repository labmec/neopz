// -*- c++ -*-
/**
 * @file pzvec_extra.h
 * @brief Extra utilities for TPZVec.
 */
// $Id: pzvec_extras.h,v 1.4 2003-04-02 13:43:58 cantao Exp $

#ifndef PZVEC_EXTRAS_H
#define PZVEC_EXTRAS_H

#include <bits/stl_algo.h>

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

#ifdef DEBUG
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


/**
 * Performs a sscal operation: x <- x * s.
 *
 * @since Mar 20, 2003
 * @author Tiago Forti
 * @author Erick Santos
 */
template< class T1, class Scalar >
void sscal(TPZVec< T1 > & x, const Scalar s)
{
#ifdef BLAS
//The blas method's call shall be done here.
#endif  //BLAS

   int size = x.NElements();

   for( int ii = 0; ii < size; ii++)
   {
   	x[ii] *= s;
   }
}

/**
 * Performs a sdot operation: dot <- Transpose[x] * y
 *
 * @since Mar 20, 2003
 * @author Tiago Forti
 * @author Erick Santos
 */
template< class T1 >
double sdot(TPZVec< T1 > & x, TPZVec< T1 > & y)
{
#ifdef BLAS
//The blas method's call shall be done here.
#endif  //BLAS

   int size = x.NElements();
   double sum = 0.0;

#ifndef NOTDEBUG
   if( size != y.NElements() ) {
      PZError << "SDOT error!" << endl
	      << "Vectors with different sizes #x = " << size
	      << ", #y = " << y.NElements() << endl;

      PZError.flush();
   }
#endif

   for( int ii = 0; ii < size; ii++)
   {
   	sum += x[ii] * y[ii];
   }

   return sum;
}

//--| SORTING |-----------------------------------------------------------------

template< class T >
TPZVec< T >& Sort( TPZVec< T >& v )
{
   std::sort( static_cast< T* >( v ), v + v.NElements() );

   return v;
}

#endif //PZVEC_EXTRAS_H

//--| PZ |----------------------------------------------------------------------

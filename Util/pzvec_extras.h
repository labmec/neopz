// -*- c++ -*-
/**
 * @file pzvec_extra.h
 * @brief Extra utilities for TPZVec.
 */
// $Id: pzvec_extras.h,v 1.5 2003-05-05 15:09:37 phil Exp $

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


template< class T, int N >
void Intersect(const TPZVec< T > &one, const TPZVec< T > &two, TPZStack< T, N > &result)
{
	int firstc, secondc, nfirst, nsecond;
	nfirst = one.NElements();
	nsecond = two.NElements();
	firstc = 0;
	secondc = 0;
	while(firstc < nfirst && secondc < nsecond) {
		while(firstc < nfirst && one[firstc] < two[secondc]) 
		{
			firstc++;
		}
		if(firstc == nfirst) break;
		while(secondc < nsecond && two[secondc] < one[firstc])
		{
			secondc++;
		}
		if(firstc < nfirst && secondc < nsecond && one[firstc] == two[secondc])
		{
			result.Push(one[firstc]);
			firstc++;
			secondc++;
		}
	}
			
}

template< class T, int N >
void Intersect(const TPZVec< T > &one, const TPZVec< T > &two, const TPZVec< T > &three, TPZStack< T, N > &result)
{
	int firstc, secondc, thirdc, nfirst, nsecond, nthird;
	nfirst = one.NElements();
	nsecond = two.NElements();
	nthird = three.NElements();
	firstc = 0;
	secondc = 0;
	thirdc = 0;
	while(firstc < nfirst && secondc < nsecond && thirdc < nthird) {
		while(firstc < nfirst && (one[firstc] < two[secondc] || one[firstc] < three[thirdc])) 
		{
			firstc++;
		}
		if(firstc==nfirst)break;
		while(secondc < nsecond && (two[secondc] < one[firstc] || two[secondc] < three[thirdc]))
		{
			secondc++;
		}
		if(secondc==nsecond) break;
		while(thirdc < nthird && (three[thirdc] < one[firstc] || three[thirdc] < two[secondc]))
		{
			thirdc++;
		}
		if(firstc < nfirst && secondc < nsecond && thirdc < nthird && one[firstc] == two[secondc] && one[firstc] == three[thirdc])
		{
			result.Push(one[firstc]);
			firstc++;
			secondc++;
			thirdc++;
		}
	}
			
}

#endif //PZVEC_EXTRAS_H

//--| PZ |----------------------------------------------------------------------

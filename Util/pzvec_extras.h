// -*- c++ -*-
/**
 * @file
 * @brief Extra utilities for TPZVec. Implementations of the saxpy, sscal, sdot, intercept, max and min functions.
 */
// $Id: pzvec_extras.h,v 1.12 2011-03-24 19:58:12 phil Exp $

#ifndef PZVEC_EXTRAS_H
#define PZVEC_EXTRAS_H

#include <algorithm>
#include "pzstack.h"


/**
 * \addtogroup util
 * @{
 */

/**
 * @since Jan 9, 2002
 * @author Cantao!
 * @brief Performs a saxpy operation: x <- x + s * y.
 */
template< class T1, class T2, class Scalar >
void saxpy( TPZVec< T1 >& x, const TPZVec< T2 >& y, Scalar s )
{
	int size = x.NElements();
	
#ifdef DEBUG
	if( size != y.NElements() ) {
		PZError << "SAXPY error!" << std::endl
		<< "Vectors with different sizes #x = " << size
		<< ", #y = " << y.NElements() << std::endl;
		
		PZError.flush();
	}
#endif
	
	for( int ii = 0; ii < size; ii++ )
	{
		x[ ii ] += s * y[ ii ];
	}
}


/**
 * @brief Performs a sscal operation: x <- x * s.
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
 * @brief Performs a sdot operation: dot <- Transpose[x] * y
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
	
#ifndef NODEBUG
	if( size != y.NElements() ) {
		PZError << "SDOT error!" << std::endl
		<< "Vectors with different sizes #x = " << size
		<< ", #y = " << y.NElements() << std::endl;
		
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
/** @brief Sorting the elements into v */
template< class T >
TPZVec< T >& Sort( TPZVec< T >& v )
{
	std::sort( static_cast< T* >( v ), v + v.NElements() );
	
	return v;
}
/** @brief Finds if exists the element e into vector v */
template< class T >
int Find( TPZVec< T >& v, const T& e)
{
	T* found = std::find( static_cast< T* >( v ), v + v.NElements(), e );
	
	int dist = distance(static_cast< T* >( v ), found);
	
	return (dist == v.NElements()? -1 : dist);
}
/** @brief Returns the minimum element into v */
template< class T >
T Min( TPZVec< T >& v )
{
	int nel = v.NElements();
	
	T m = v[ 0 ];
	
	for( int ii = 1; ii < nel; ii++ )
	{
		m = min( m, v[ ii ] );
	}
	
	return m;
}
/** @brief Returns the maximum element into v */
template< class T >
T Max( TPZVec< T >& v )
{
	int nel = v.NElements();
	
	T m = v[ 0 ];
	
	for( int ii = 1; ii < nel; ii++ )
	{
		m = max( m, v[ ii ] );
	}
	
	return m;
}
/** @brief Gets commom elements into the one and two vectors */
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
/** @brief Gets commom elements into the one, two and three vectors */
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

/**
 * @}
 */

#endif //PZVEC_EXTRAS_H

//--| PZ |----------------------------------------------------------------------

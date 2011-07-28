// -*- c++ -*-
/******************************************************************************
 *
 * Data definition:     REAL
 *
 * Purpose:             Define what type of data the package is being compiled
 *                      for.  The choices are the following:
 *
 *                          double   real, double precision
 *                          float    real, single precision
 *
 *****************************************************************************/

#ifndef REALH
#define REALH

#include <math.h>
#include <iostream>

#ifdef WIN32
#define  __PRETTY_FUNCTION__ __FILE__
#endif

#ifndef ELLIPS

/** \addtogroup common 
 * @{
 */

/// Extern variable to control level of printting? ( or priority print? Jorge)
extern int gPrintLevel;


/// Operations to be counted: Sum, Product, Division, Square root, Power, \n
/// Sine, Cosine, Arc Sine, Arc Cosine, Arc Tangent, Exponencial and Logarithm.
/** 
 \brief Types of operations to be counted.
 */
enum MOptype {ESum,EProd,EDiv,ESqrt,EPow,ECos,ESin,EAcos,EAsin,EAtan,EExp,ELog};
const int gNumOp = 12;


/// @brief This class implements a counter
struct TPZCounter {
	/// Vector of counters by operation: sum, product, division, square root, power, \n
	/// Cosine, Sine, Arc cosine, arc Sine, arc Tangent, Exponencial and logarithm.
	int fCount[gNumOp];
	
	inline TPZCounter operator-(const TPZCounter &other)
	{
		TPZCounter result;
		int i;
		for(i=0; i<gNumOp; i++) result.fCount[i] = fCount[i]-other.fCount[i];
		return result;
	}
	
	void Print(std::ostream &out = std::cout) const;
	
};

std::ostream &operator<<(std::ostream &out,const TPZCounter &count);

// Note: How to sure the counter initializing from ZERO ??? I don't see where the counter vector is clean. (Jorge)

/**
 * By modifying the definition of the type of REAL, the number of floating point operations can be counted. \n
 * To modify you need to define "contar": \#define contar. (Jorge)
 *
 * @brief This class implements floating point number associated with a counter of the operations performed with it by type. \n
 * (arithmetic, trigonometric, exponencial and logarithmic operations performed with it)
 */
class TPZFlopCounter {
public:
	/// Floating point value
	double fVal;
	/// Containts the counter vector by operation performed
	static TPZCounter gCount;
	
	inline TPZFlopCounter()
	{
	}
	inline TPZFlopCounter(const double &val)
	{
		fVal = val;
	}
	/// Returns the product with the oth value and increments the counter of the products.
	inline TPZFlopCounter operator*(const TPZFlopCounter &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal*oth.fVal;
		gCount.fCount[EProd]++;
		return result;
	}
	/// Returns the division between oth value and increments the counter of the divisions.
	inline TPZFlopCounter operator/(const TPZFlopCounter &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal/oth.fVal;
		gCount.fCount[EDiv]++;
		return result;
	}
	
	/// Returns the sum with oth value and increments the counter of the sums.
	inline TPZFlopCounter operator+(const TPZFlopCounter &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal+oth.fVal;
		gCount.fCount[ESum]++;
		return result;
	}
	
	/// Returns the difference with oth value and increments the counter of the sums.
	inline TPZFlopCounter operator-(const TPZFlopCounter &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal-oth.fVal;
		gCount.fCount[ESum]++;
		return result;
	}
	/// Performs the sum with oth value on own value and increments the counter of the sums.
	inline TPZFlopCounter &operator+=(const TPZFlopCounter &oth) 
	{
		fVal += oth.fVal;
		gCount.fCount[ESum]++;
		return *this;
	}
	/// Performs the diference with oth value on own value and increments the counter of the sums.
	inline TPZFlopCounter &operator-=(const TPZFlopCounter &oth) 
	{
		fVal -= oth.fVal;
		gCount.fCount[ESum]++;
		return *this;
	}
	/// Performs the product with oth value on own value and increments the counter of the products.
	inline TPZFlopCounter &operator*=(const TPZFlopCounter &oth) 
	{
		fVal *= oth.fVal;
		gCount.fCount[EProd]++;
		return *this;
	}
	/// Performs the division between oth value on own value and increments the counter of the divisions.
	inline TPZFlopCounter &operator/=(const TPZFlopCounter &oth) 
	{
		fVal /= oth.fVal;
		gCount.fCount[EDiv]++;
		return *this;
	}
	/// Returns value with signal changed of its floating point value and increments the counter of the sums.
	inline TPZFlopCounter operator-() const
	{
		TPZFlopCounter result;
		result.fVal = -fVal;
		gCount.fCount[ESum]++;
		return result;
	}
	/// Returns the current floating point value and doesn't increments the counters.
	inline TPZFlopCounter operator+() const
	{
		return *this;
	}
	/// Compares the values and doesn't increments the counters.
	inline bool operator<(const TPZFlopCounter &sec) const
	{
		return (fVal < sec.fVal);
	}
	/// Compares the values and doesn't increments the counters.
	inline bool operator>(const TPZFlopCounter &sec) const
	{
		return (fVal > sec.fVal);
	}
	/// Compares the values and doesn't increments the counters.
	inline bool operator<=(const TPZFlopCounter &sec) const
	{
		return (fVal <= sec.fVal);
	}
	/// Compares the values and doesn't increments the counters.
	inline bool operator>=(const TPZFlopCounter &sec) const
	{
		return (fVal >= sec.fVal);
	}
	/// Compares the values and doesn't increments the counters.
	inline bool operator==(const TPZFlopCounter &sec) const
	{
		return (fVal == sec.fVal);
	}
	/// Compares the values and doesn't increments the counters.
	inline bool operator!=(const TPZFlopCounter &sec) const
	{
		return (fVal != sec.fVal);
	}
	
	/// Returns the square root of the value and increments the counter of the square root.
	friend TPZFlopCounter sqrt(const TPZFlopCounter &orig);
	/// Returns the power and increments the counter of the power.
	friend TPZFlopCounter pow(const TPZFlopCounter &orig,const TPZFlopCounter &xp);
	/// Returns the absolute value and doesn't increments the counters.
	friend TPZFlopCounter fabs(const TPZFlopCounter &orig);
	/// Returns the arc cosine in radians and increments the counter of the Arc Cosine.
	friend TPZFlopCounter acos(const TPZFlopCounter &orig);
	/// Returns the cosine in radians and increments the counter of the Cosine.
	friend TPZFlopCounter cos(const TPZFlopCounter &orig);
	/// Returns the arc sine in radians and increments the counter of the Arc Sine.
	friend TPZFlopCounter asin(const TPZFlopCounter &orig);
	/// Returns the sine in radians and increments the counter of the sine.
	friend TPZFlopCounter sin(const TPZFlopCounter &orig);
	/// Returns the arc tangent in radians and increments the counter of the Arc Tangent.
	friend TPZFlopCounter atan(const TPZFlopCounter &orig);
	/// Returns the arc tangent in radians and increments the counter of the Arc Tangent.
	/// ATAN2 returns the arc tangent of x(val1) and y(val2) coordinates as an angle expressed in radians.
	friend TPZFlopCounter atan2(const TPZFlopCounter &val1,const TPZFlopCounter &val2);
	/// Returns the exponencial and increments the counter of the Exponencial.
	friend TPZFlopCounter exp(const TPZFlopCounter &val);
	/// Returns the natural logarithm and increment the counter of the logarithm.
	friend TPZFlopCounter log(const TPZFlopCounter &val);
	/// Returns the decimal logarithm and increment the counter of the logarithm.
	friend TPZFlopCounter log10(const TPZFlopCounter &val);
};

inline TPZFlopCounter sqrt(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::sqrt(orig.fVal);
	TPZFlopCounter::gCount.fCount[ESqrt]++;
	return result;
}

inline TPZFlopCounter fabs(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::fabs(orig.fVal);
	return result;
}

inline TPZFlopCounter pow(const TPZFlopCounter &orig,const TPZFlopCounter &xp)
{
	TPZFlopCounter result;
	result.fVal = ::pow(orig.fVal,xp.fVal);
	TPZFlopCounter::gCount.fCount[EPow]++;
	return result;
}

inline TPZFlopCounter acos(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::acos(orig.fVal);
	TPZFlopCounter::gCount.fCount[EAcos]++;
	return result;
}

inline TPZFlopCounter asin(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::asin(orig.fVal);
	TPZFlopCounter::gCount.fCount[EAsin]++;
	return result;
}

inline TPZFlopCounter cos(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::cos(orig.fVal);
	TPZFlopCounter::gCount.fCount[ECos]++;
	return result;
}

inline TPZFlopCounter sin(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::sin(orig.fVal);
	TPZFlopCounter::gCount.fCount[ESin]++;
	return result;
}

inline TPZFlopCounter atan(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::atan(orig.fVal);
	TPZFlopCounter::gCount.fCount[EAtan]++;
	return result;
}

inline TPZFlopCounter atan2(const TPZFlopCounter &val1,const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result.fVal = ::atan2(val1.fVal,val2.fVal);
	TPZFlopCounter::gCount.fCount[EAtan]++;
	return result;
}

inline TPZFlopCounter exp(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::exp(orig.fVal);
	TPZFlopCounter::gCount.fCount[EExp]++;
	return result;
}

inline TPZFlopCounter log(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::log(orig.fVal);
	TPZFlopCounter::gCount.fCount[ELog]++;
	return result;
}

inline TPZFlopCounter log10(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::log10(orig.fVal);
	TPZFlopCounter::gCount.fCount[ELog]++;
	return result;
}
/// Performs val1 * val2. Doesn't increments counters.
inline TPZFlopCounter operator*(double val1, const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)*val2;
	return result;
	
}
/// Performs val1 / val2. Doesn't increments counters.	
inline TPZFlopCounter operator/(double val1, const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)/val2;
	return result;
	
}
/// Performs val1 + val2. Doesn't increments counters.	
inline TPZFlopCounter operator+(double val1, const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)+val2;
	return result;
	
}
/// Performs val1 + val2. Doesn't increments counters.	
inline TPZFlopCounter operator+(const TPZFlopCounter &val2, double val1 )
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)+val2;
	return result;
	
}
/// Performs val1 - val2. Doesn't increments counters.
inline TPZFlopCounter operator-(double val1, const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)-val2;
	return result;
	
}

/// Implements to write (output) only the floating point value
inline std::ostream &operator<<(std::ostream &out, const TPZFlopCounter &val)
{
	return out << val.fVal;
}
/// Implements to read (input) only the floating point value
inline std::istream &operator>>(std::istream &out, /*const*/ TPZFlopCounter &val)
{
	return out >> val.fVal;
}
#endif

#ifdef contar
/// PZ will use a double floating point number with a counter of the operations performed on (Jorge)
typedef std::TPZFlopCounter REAL;

#else
/// This is the type of floating point number PZ will use.
typedef double REAL;

#endif

/// Returns the tolerance to Zero value. Actually: 1e-10
inline REAL ZeroTolerance(){
	return 1e-10;
}
/** @brief Returns if the value a is close Zero as the allowable tolerance */
inline bool IsZero( REAL a ){
	return ( fabs( a ) < ZeroTolerance() );
}
/// Returns the maximum value between a and b
template <class T>
inline const T& Max( const T &a, const T &b ){
	return a > b ? a : b;
}
/// Returns the minimum value between a and b
template <class T>
inline const T& Min( const T & a, const T &b ){
	return a < b ? a : b;
}

/** @}
 */

#endif

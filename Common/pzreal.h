/**
 * @file
 * @brief Contains the declaration of TPZFlopCounter class and TPZCounter struct.
 */
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
#include <complex>
#include <config.h>

#ifdef WIN32
#define  __PRETTY_FUNCTION__ __FILE__
#endif

#ifndef ELLIPS

/** \addtogroup common 
 * @{
 */

/** @brief Extern variable to control level of printting (priority print?) */
extern int gPrintLevel;

/**
 * Operations to be counted: Sum, Product, Division, Square root, Power, \n
 * Sine, Cosine, Arc Sine, Arc Cosine, Arc Tangent, Exponencial and Logarithm.
 */
/** @brief Types of operations to be counted. */
enum MOptype {ESum,EProd,EDiv,ESqrt,EPow,ECos,ESin,EAcos,EAsin,EAtan,EExp,ELog};

/** @brief Number of type of the operations actually counted. */
const int gNumOp = 12;


/** @brief Implements a counter by operations. \ref common "Common" */
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

/** @brief Re-implements << operator to show the counter (count) data */
std::ostream &operator<<(std::ostream &out,const TPZCounter &count);

class TPZFlopCounter;

#ifdef contar
/** @brief PZ will use a double floating point number with a counter of the operations performed on (Jorge) */
typedef TPZFlopCounter REAL;
#endif

/** @brief This is the type of floating point number PZ will use. */
#ifdef REALfloat
typedef float REAL;
#endif
#ifdef REALdouble
typedef double REAL; //This is the default configuration
#endif
#ifdef REALlongdouble
typedef long double REAL;
#endif

/** @brief This is the type of State PZ will use. */
#ifdef STATEfloat
typedef float STATE;
#elif STATEdouble
typedef double STATE; //This is the default configuration
#elif STATElongdouble
typedef long double STATE;
#elif STATEcomplexf
typedef std::complex<float> STATE;
#elif STATEcomplexd
typedef std::complex<double> STATE;
#elif STATEcomplexld
typedef std::complex<long double> STATE;
#endif

// fabs function adapted to complex numbers.
inline float
fabs(std::complex <float> __x)
{
  std::complex <float> ret = abs(__x);
  return std::abs(ret.real());
}
inline double
fabs(std::complex <double> __x)
{
  std::complex <double> ret = abs(__x);
  return std::abs(ret.real());
}
inline long double
fabs(std::complex <long double> __x)
{
  std::complex <long double> ret = abs(__x);
  return std::abs(ret.real());
}


/**
 * By modifying the definition of the type of REAL, the number of floating point operations can be counted. \n
 * To modify you need to define "contar": \#define contar.
 */
/** @brief This class implements floating point number associated with a counter of the operations performed with it by type. \n
 * (arithmetic, trigonometric, exponencial and logarithmic operations performed with it). \ref common "Common"
 * @note How to sure the counter initializing from ZERO ??? I don't see where the counter vector is clean. (Jorge)
 */
class TPZFlopCounter {
public:
	/** @brief Floating point value */
#ifdef contar
	double fVal;
#else
	REAL fVal;
#endif
	/** @brief Containts the counter vector by operation performed */
	static TPZCounter gCount;
	
	inline TPZFlopCounter()
	{
	}
	inline TPZFlopCounter(const double &val)
	{
		fVal = val;
	}
	/** @brief Returns the product with the oth value and increments the counter of the products. */
	inline TPZFlopCounter operator*(const TPZFlopCounter &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal*oth.fVal;
		gCount.fCount[EProd]++;
		return result;
	}
	/** @brief Returns the division between oth value and increments the counter of the divisions. */
	inline TPZFlopCounter operator/(const TPZFlopCounter &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal/oth.fVal;
		gCount.fCount[EDiv]++;
		return result;
	}
	
	/** @brief Returns the sum with oth value and increments the counter of the sums. */
	inline TPZFlopCounter operator+(const TPZFlopCounter &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal+oth.fVal;
		gCount.fCount[ESum]++;
		return result;
	}
	
	/** @brief Returns the difference with oth value and increments the counter of the sums. */
	inline TPZFlopCounter operator-(const TPZFlopCounter &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal-oth.fVal;
		gCount.fCount[ESum]++;
		return result;
	}
	/** @brief Performs the sum with oth value on own value and increments the counter of the sums. */
	inline TPZFlopCounter &operator+=(const TPZFlopCounter &oth) 
	{
		fVal += oth.fVal;
		gCount.fCount[ESum]++;
		return *this;
	}
	/** @brief Performs the diference with oth value on own value and increments the counter of the sums. */
	inline TPZFlopCounter &operator-=(const TPZFlopCounter &oth) 
	{
		fVal -= oth.fVal;
		gCount.fCount[ESum]++;
		return *this;
	}
	/** @brief Performs the product with oth value on own value and increments the counter of the products. */
	inline TPZFlopCounter &operator*=(const TPZFlopCounter &oth) 
	{
		fVal *= oth.fVal;
		gCount.fCount[EProd]++;
		return *this;
	}
	/** @brief Performs the division between oth value on own value and increments the counter of the divisions. */
	inline TPZFlopCounter &operator/=(const TPZFlopCounter &oth) 
	{
		fVal /= oth.fVal;
		gCount.fCount[EDiv]++;
		return *this;
	}
	/** @brief Returns value with signal changed of its floating point value and increments the counter of the sums. */
	inline TPZFlopCounter operator-() const
	{
		TPZFlopCounter result;
		result.fVal = -fVal;
		gCount.fCount[ESum]++;
		return result;
	}
	/** @brief Returns the current floating point value and doesn't increments the counters. */
	inline TPZFlopCounter operator+() const
	{
		return *this;
	}
	/** @brief Compares the values and doesn't increments the counters. */
	inline bool operator<(const TPZFlopCounter &sec) const
	{
		return (fVal < sec.fVal);
	}
	/** @brief Compares the values and doesn't increments the counters. */
	inline bool operator>(const TPZFlopCounter &sec) const
	{
		return (fVal > sec.fVal);
	}
	/** @brief Compares the values and doesn't increments the counters. */
	inline bool operator<=(const TPZFlopCounter &sec) const
	{
		return (fVal <= sec.fVal);
	}
	/** @brief Compares the values and doesn't increments the counters. */
	inline bool operator>=(const TPZFlopCounter &sec) const
	{
		return (fVal >= sec.fVal);
	}
	/** @brief Compares the values and doesn't increments the counters. */
	inline bool operator==(const TPZFlopCounter &sec) const
	{
		return (fVal == sec.fVal);
	}
	/** @brief Compares the values and doesn't increments the counters. */
	inline bool operator!=(const TPZFlopCounter &sec) const
	{
		return (fVal != sec.fVal);
	}
	
	/** @brief Returns the square root of the value and increments the counter of the square root. */
	friend TPZFlopCounter sqrt(const TPZFlopCounter &orig);
	/** @brief Returns the power and increments the counter of the power. */
	friend TPZFlopCounter pow(const TPZFlopCounter &orig,const TPZFlopCounter &xp);
	/** @brief Returns the absolute value and doesn't increments the counters. */
	friend TPZFlopCounter fabs(const TPZFlopCounter &orig);
	/** @brief Returns the arc cosine in radians and increments the counter of the Arc Cosine. */
	friend TPZFlopCounter acos(const TPZFlopCounter &orig);
	/** @brief Returns the cosine in radians and increments the counter of the Cosine. */
	friend TPZFlopCounter cos(const TPZFlopCounter &orig);
	/** @brief Returns the arc sine in radians and increments the counter of the Arc Sine. */
	friend TPZFlopCounter asin(const TPZFlopCounter &orig);
	/** @brief Returns the sine in radians and increments the counter of the sine. */
	friend TPZFlopCounter sin(const TPZFlopCounter &orig);
	/** @brief Returns the arc tangent in radians and increments the counter of the Arc Tangent. */
	friend TPZFlopCounter atan(const TPZFlopCounter &orig);
	/** 
	 * @brief Returns the arc tangent in radians and increments the counter of the Arc Tangent. \n
	 * ATAN2 returns the arc tangent of x(val1) and y(val2) coordinates as an angle expressed in radians.
	 */
	friend TPZFlopCounter atan2(const TPZFlopCounter &val1,const TPZFlopCounter &val2);
	/** @brief Returns the exponencial and increments the counter of the Exponencial. */
	friend TPZFlopCounter exp(const TPZFlopCounter &val);
	/** @brief Returns the natural logarithm and increment the counter of the logarithm. */
	friend TPZFlopCounter log(const TPZFlopCounter &val);
	/** @brief Returns the decimal logarithm and increment the counter of the logarithm. */
	friend TPZFlopCounter log10(const TPZFlopCounter &val);
};

/** @brief Returns the square root of the value and increments the counter of the square root. */
inline TPZFlopCounter sqrt(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::sqrt(orig.fVal);
	TPZFlopCounter::gCount.fCount[ESqrt]++;
	return result;
}

/** @brief Returns the absolute value and doesn't increments the counters. */
inline TPZFlopCounter fabs(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::fabs(orig.fVal);
	return result;
}

/** @brief Returns the power and increments the counter of the power. */
inline TPZFlopCounter pow(const TPZFlopCounter &orig,const TPZFlopCounter &xp)
{
	TPZFlopCounter result;
	result.fVal = ::pow(orig.fVal,xp.fVal);
	TPZFlopCounter::gCount.fCount[EPow]++;
	return result;
}

/** @brief Returns the arc cosine in radians and increments the counter of the Arc Cosine. */
inline TPZFlopCounter acos(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::acos(orig.fVal);
	TPZFlopCounter::gCount.fCount[EAcos]++;
	return result;
}

/** @brief Returns the arc sine in radians and increments the counter of the Arc Sine. */
inline TPZFlopCounter asin(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::asin(orig.fVal);
	TPZFlopCounter::gCount.fCount[EAsin]++;
	return result;
}

/** @brief Returns the cosine in radians and increments the counter of the Cosine. */
inline TPZFlopCounter cos(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::cos(orig.fVal);
	TPZFlopCounter::gCount.fCount[ECos]++;
	return result;
}

/** @brief Returns the sine in radians and increments the counter of the sine. */
inline TPZFlopCounter sin(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::sin(orig.fVal);
	TPZFlopCounter::gCount.fCount[ESin]++;
	return result;
}

/** @brief Returns the arc tangent in radians and increments the counter of the Arc Tangent. */
inline TPZFlopCounter atan(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::atan(orig.fVal);
	TPZFlopCounter::gCount.fCount[EAtan]++;
	return result;
}

/** 
 * @brief Returns the arc tangent in radians and increments the counter of the Arc Tangent. \n
 * ATAN2 returns the arc tangent of x(val1) and y(val2) coordinates as an angle expressed in radians.
 */
inline TPZFlopCounter atan2(const TPZFlopCounter &val1,const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result.fVal = ::atan2(val1.fVal,val2.fVal);
	TPZFlopCounter::gCount.fCount[EAtan]++;
	return result;
}

/** @brief Returns the exponencial and increments the counter of the Exponencial. */
inline TPZFlopCounter exp(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::exp(orig.fVal);
	TPZFlopCounter::gCount.fCount[EExp]++;
	return result;
}

/** @brief Returns the natural logarithm and increment the counter of the logarithm. */
inline TPZFlopCounter log(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::log(orig.fVal);
	TPZFlopCounter::gCount.fCount[ELog]++;
	return result;
}

/** @brief Returns the decimal logarithm and increment the counter of the logarithm. */
inline TPZFlopCounter log10(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = ::log10(orig.fVal);
	TPZFlopCounter::gCount.fCount[ELog]++;
	return result;
}
/** @brief Performs \f$ val1 * val2\f$. Doesn't increments counters. */
inline TPZFlopCounter operator*(double val1, const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)*val2;
	return result;
	
}
/** @brief Performs \f$ val1 / val2\f$. Doesn't increments counters. */
inline TPZFlopCounter operator/(double val1, const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)/val2;
	return result;
	
}
/** @brief Performs \f$ val1 + val2\f$. Doesn't increments counters. */
inline TPZFlopCounter operator+(double val1, const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)+val2;
	return result;
	
}
/** @brief Performs \f$ val1 + val2\f$. Doesn't increments counters. */
inline TPZFlopCounter operator+(const TPZFlopCounter &val2, double val1 )
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)+val2;
	return result;
	
}
/** @brief Performs \f$ val1 - val2\f$. Doesn't increments counters. */
inline TPZFlopCounter operator-(double val1, const TPZFlopCounter &val2)
{
	TPZFlopCounter result;
	result = TPZFlopCounter(val1)-val2;
	return result;
	
}

/** @brief Implements to write (output) only the floating point value */
inline std::ostream &operator<<(std::ostream &out, const TPZFlopCounter &val)
{
	return out << val.fVal;
}
/** @brief Implements to read (input) only the floating point value */
inline std::istream &operator>>(std::istream &out, /*const*/ TPZFlopCounter &val)
{
	return out >> val.fVal;
}
#endif


/** @brief Returns the tolerance to Zero value. Actually: \f$ 1e-10 \f$ */
inline REAL ZeroTolerance(){
	return 1e-10;
}
#ifdef _AUTODIFF
/** @brief Returns if the value a is close Zero as the allowable tolerance */
template<class T>
inline bool IsZero( T a ){
	return ( fabs( a.val() ) < ZeroTolerance() );
}
#endif
/** @brief Returns if the value a is close Zero as the allowable tolerance */
//template<>
inline bool IsZero( long double a ){
	return ( fabs( a ) < ZeroTolerance() );
}
//template<>
inline bool IsZero( double a ){
	return ( fabs( a ) < ZeroTolerance() );
}
//template<>
inline bool IsZero( float a ){
	return ( fabs( a ) < ZeroTolerance() );
}
//template<>
inline bool IsZero( std::complex<long double> a ){
	return ( fabs( a ) < ZeroTolerance() );
}
//template<>
inline bool IsZero( std::complex<double> a ){
	return ( fabs( a ) < ZeroTolerance() );
}
//template<>
inline bool IsZero( std::complex<float> a ){
	return ( fabs( a ) < ZeroTolerance() );
}
//template<>
inline bool IsZero( int a ){
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

// In the math library (cmath.h) don't exist some overloading for some functions

inline float
pow(float __x, double __y)
{ return pow(__x, (float)(__y)); }

#endif

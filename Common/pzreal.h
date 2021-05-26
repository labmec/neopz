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

#include <pz_config.h>
#include <limits>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#endif // _USE_MATH_DEFINES
//additional math defines
#define M_SQRT3   1.732050807568877293527446341505872366
#define M_SQRT1_3 0.577350269189625764509148780501957455
#define M_SQRT6   2.449489742783178098197284074705891391
#define M_SQRT1_6 0.408248290463863016366214012450981898
#include <iostream>
#include <complex>
#include "fpo_exceptions.h"

template <int Num, class T> class TFad;
template <class T> class Fad;
template <typename Enumeration>
typename std::underlying_type<Enumeration>::type as_integer(const Enumeration value) {
    return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}

#ifdef __GNUC__ // GCC 4.8+, Clang, Intel and other compilers compatible with GCC (-std=c++0x or above)
[[noreturn]] inline __attribute__((always_inline)) void unreachable() {__builtin_unreachable();}
#elif defined(_MSC_VER) // MSVC
[[noreturn]] __forceinline void unreachable() {__assume(false);}
#else // ???
inline void unreachable() {}
#endif


class TPZFlopCounter;

//! Enables typename real_type<T>::type for the real type associated with T
template <typename, typename = void>
struct real_type;

//! Enabling it for arithmetic types
template <typename T>
struct real_type<T, std::enable_if_t<std::is_arithmetic_v<T>||std::is_integral_v<T>>>
 { using type = T; };

//! Enabling it for complex types
template <typename T>
struct real_type<std::complex<T>, void>
 { using type = T; };


//! Enabling it for TPZFlopCounter
template<>
struct real_type<TPZFlopCounter, void>
{ using type = TPZFlopCounter; };

//Convenience macros for real_type<T>::type
#define RTVar typename real_type<TVar>::type
#define RType(T) typename real_type<T>::type

//! Enabling it for Fad<T>
template <typename T>
struct real_type<Fad<T>, void>
{ using type =Fad<RType(T)>;};
//! Enabling it for TFad<6,T>
template <int N,typename T>
struct real_type<TFad<N,T>, void>
{ using type =TFad<N,RType(T)>;};

//! Enables typename complex_type<T>::type for the complex type associated with T
template <typename, typename = void>
struct complex_type;

//! Enabling it for arithmetic types
template <typename T>
struct complex_type<T, std::enable_if_t<std::is_arithmetic_v<T>||std::is_integral_v<T>>>
{ using type = std::complex<T>; };

//! Enabling it for complex types
template <typename T>
struct complex_type<std::complex<T>, void>
{ using type = std::complex<T>; };

//! Enabling it for TPZFlopCounter
template<>
struct complex_type<TPZFlopCounter, void>
{ using type = TPZFlopCounter; };

//Convenience macros for complex_type<T>::type
#define CTVar typename complex_type<TVar>::type
#define CType(T) typename complex_type<T>::type

//! Enabling it for Fad<T>
template <typename T>
struct complex_type<Fad<T>, void>
{ using type =Fad<CType(T)>;};
//! Enabling it for TFad<6,T>
template <int N,typename T>
struct complex_type<TFad<N,T>, void>
{ using type =TFad<N,CType(T)>;};

template<class T>
struct is_complex{ static constexpr bool value = false;};

template<class T>
struct is_complex<std::complex<T>> : 
    std::integral_constant<bool,
    std::is_integral<T>::value ||
    std::is_floating_point<T>::value>{};

template<class T>
struct is_arithmetic_pz :
    std::integral_constant<bool,
    std::is_integral<T>::value ||
    std::is_floating_point<T>::value ||
    is_complex<T>::value>{};

//! Enables typename complex_type<T>::type for the complex type associated with T
template <typename, typename = void>
struct is_fad{
	static constexpr bool value{false};
};

template<class T>
struct  is_fad<Fad<T>>{
	static constexpr bool value{true};
};
template<int N, class T>
struct  is_fad<TFad<N,T>>{
	static constexpr bool value{true};
};

/** @brief Gets maxime value between a and b */
#ifndef MAX
#define MAX( a, b )   ( (a) > (b) ? (a) : (b) )
#endif
/** @brief Gets minime value between a and b */
#ifndef MIN
#define MIN( a, b )   ( (a) < (b) ? (a) : (b) )
#endif

#ifdef WIN32
#define  __PRETTY_FUNCTION__ __FILE__
#endif


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
        uint64_t fCount[gNumOp];
	
        /// Counter constructor.
        TPZCounter() 
        {
	  clear();
	}

	inline TPZCounter operator-(const TPZCounter &other)
	{
	  TPZCounter result;
	  for(int i=0; i<gNumOp; i++) result.fCount[i] = fCount[i]-other.fCount[i];
	  return result;
	}

	inline TPZCounter& operator+=(const TPZCounter &other)
	{
	  for(int i=0; i<gNumOp; i++) fCount[i] += other.fCount[i];
	  return *this;
	}

	inline TPZCounter& operator-=(const TPZCounter &other)
	{
	  for(int i=0; i<gNumOp; i++) fCount[i] -= other.fCount[i];
	  return *this;
	}

	inline TPZCounter& copy(const TPZCounter &other)
	{
	  for (int i=0; i<gNumOp; i++) fCount[i] = other.fCount[i];
	  return *this;
	}

        inline TPZCounter& clear()
	{
	  for(int i=0; i<gNumOp; i++) fCount[i] = 0 ;
	  return *this;
	}
	
	void Print(std::ostream &out = std::cout) const;
	
};

/** @brief Re-implements << operator to show the counter (count) data */
std::ostream &operator<<(std::ostream &out,const TPZCounter &count);

/** @brief This is the type of floating point number PZ will use. */
#ifdef REALfloat
typedef float REAL;
#endif // REALfloat
#ifdef REALdouble
typedef double REAL; //This is the default configuration
#endif // REALdouble
#ifdef REALlongdouble
typedef long double REAL;
#endif // REALlongdouble
#ifdef REALpzfpcounter
typedef TPZFlopCounter REAL;
#endif // REALpzfpcounter

/** @brief This is the type of State PZ will use. */
#ifdef STATEfloat
typedef float STATE;
#endif // STATEfloat
#ifdef STATEdouble
//This is the default configuration
typedef double STATE;
#endif // STATEdouble
#ifdef STATElongdouble
typedef long double STATE;
#endif // STATElongdouble
typedef std::complex<STATE> CSTATE;
// set(STATE_COMPLEX "STATE_COMPLEX")

enum ESolType{ EReal=1,EComplex=2,EUndefined=3};
#ifdef VC
#include <io.h>
#ifndef NOMINMAX
#define NOMINMAX // Preventing the redefinition of min and max as macros
#endif // NOMINMAX
#include <Windows.h>
// sqrt function adapted to int numbers. required for VC
inline double
sqrt(int __x)
{
  return sqrt ((double)__x);
}
// fabs function adapted to int numbers. required for VC
inline int
fabs(int __x)
{
  return (int)fabs((double) __x);
}
inline double
sin(int __x)
{
  return sin((double) __x);
}
inline double
cos(int __x)
{
  return cos((double) __x);
}
inline double
atan(int __x)
{
  return atan((double) __x);
}
#endif //VC

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
 * To modify you need to define "REALpzfpcounter": \#define REALpzfpcounter
.
 */
/** @brief This class implements floating point number associated with a counter of the operations performed with it by type. \n
 * (arithmetic, trigonometric, exponencial and logarithmic operations performed with it). \ref common "Common"
 * @note How to sure the counter initializing from ZERO ??? I don't see where the counter vector is clean. (Jorge)
 */
class TPZFlopCounter {
public:
	/** @brief Floating point value */
#ifdef REALpzfpcounter
	double fVal;
#else
	REAL fVal;
#endif //REALpzfpcounter
	/** @brief Containts the counter vector by operation performed */
	static TPZCounter gCount;

	inline TPZFlopCounter()
	{
	}
	inline TPZFlopCounter(const double &val)
	{
		fVal = (REAL)val;
	}
    
    inline REAL val() const
    {
        return  fVal;
    }
    
    operator REAL() const{
        return fVal;
    }
    
    bool operator<=(const REAL &val) const
    {
        return fVal <= val;
    }
    
    bool operator<(const REAL &val) const
    {
        return fVal < val;
    }
    bool operator>(const REAL &val) const
    {
        return fVal > val;
    }
    bool operator>=(const REAL &val) const
    {
        return fVal >= val;
    }
    bool operator==(const REAL &val) const
    {
        return fVal == val;
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
	
	/** @brief Returns the product with the oth value and increments the counter of the products. */
	inline TPZFlopCounter operator*(const REAL &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal*oth;
		gCount.fCount[EProd]++;
		return result;
	}
	/** @brief Returns the division between oth value and increments the counter of the divisions. */
	inline TPZFlopCounter operator/(const double &oth) const
	{
		TPZFlopCounter result;
		result.fVal = fVal/((REAL)oth);
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
	friend TPZFlopCounter fabsFlop(const TPZFlopCounter &orig);
	/** @brief Returns the absolute value as REAL and doesn't increments the counters. */
	friend REAL fabs(const TPZFlopCounter &orig);
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
inline TPZFlopCounter fabsFlop(const TPZFlopCounter &orig)
{
	TPZFlopCounter result;
	result.fVal = std::abs(orig.fVal);
	return result;
}
/** @brief Returns the absolute value as REAL and doesn't increments the counters. */
inline REAL fabs(const TPZFlopCounter &orig)
{
    return std::abs(orig.fVal);
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



/** @brief Returns the tolerance to Zero value. Actually: \f$ 1e-10 \f$ */
inline REAL ZeroTolerance() {
    typedef std::numeric_limits< REAL > dbl;
	return (REAL)pow(10,(-1 * (dbl::max_digits10- 5)));
}

template<typename T,
        typename std::enable_if< is_arithmetic_pz<T>::value, int>::type* = nullptr>
inline void ZeroTolerance(T &tol) {
    typedef std::numeric_limits< T > dbl;
    tol = (T)pow(10,(-1 * (dbl::max_digits10- 5)));
}

template<typename T,
        typename std::enable_if<std::is_same<T,TPZFlopCounter>::value>::type* = nullptr>
inline void ZeroTolerance(T &tol) {
    typedef std::numeric_limits< REAL > dbl;
    tol = (T)pow(10,(-1 * (dbl::max_digits10- 5)));
}

template<int Num, typename T>
inline void ZeroTolerance(TFad<Num,T> &Tol) {
    ZeroTolerance(Tol.val());
}

/** @brief Returns if the value a is close Zero as the allowable tolerance */
template<class T>
inline bool IsZero( T a ) {
	return ( std::abs( a.val() ) < ZeroTolerance() );
}
/** @brief Returns if the value a is close Zero as the allowable tolerance */
//template<>
inline bool IsZero( long double a ) {
#ifdef WIN32
	return ( fabs( a ) < 1.e-12 );
#else
	return ( std::fabs( a ) < 1.e-16 );
#endif
}
//template<>
inline bool IsZero( double a ) {
#ifdef WIN32
	return ( fabs( a ) < 1.e-10 );
#else
	return ( fabs( a ) < 1.e-12 );
#endif
}
//template<>
inline bool IsZero( float a ) {
#ifdef WIN32
	return ( fabs( a ) < 1.e-6 );
#else
	return ( fabs( a ) < 1.e-6 );
#endif
}
//template<>
inline bool IsZero( std::complex<long double> a ) {
#ifdef WIN32
	return ( fabs( a ) < 1.e-12 );
#else
	return ( fabs( a ) < 1.e-16 );
#endif
}
//template<>
inline bool IsZero( std::complex<double> a ) {
#ifdef WIN32
	return ( fabs( a ) < 1.e-9 );
#else
	return ( fabs( a ) < 1.e-12 );
#endif
}
//template<>
inline bool IsZero( std::complex<float> a ) {
#ifdef WIN32
	return ( fabs( a ) < 1.e-5 );
#else
	return ( fabs( a ) < 1.e-7 );
#endif
}
//template<>
inline bool IsZero( int a ) {
	return ( a==0 );
}
inline bool IsZero( int64_t a ) {
	return ( a==0L );
}
/// Returns the maximum value between a and b
template <class T>
inline const T& Max( const T &a, const T &b ) {
	return a > b ? a : b;
}
/// Returns the minimum value between a and b
template <class T>
inline const T& Min( const T & a, const T &b ) {
	return a < b ? a : b;
}

/** @}
 */

// In the math library (cmath.h) don't exist some overloading for some functions

// SPECIAL FUNCTIONS NON STANDARD IN WINDOWS SYSTEM
#if (!defined(__cplusplus) || __cplusplus < 201103L) && (!defined(_MSC_VER) || _MSC_VER < 1900)// If we aren't using C++11.
/**
 * Function erf (Error function) implemented in 
 * http://www.johndcook.com/cpp_erf.html
 */
REAL erf(REAL arg);

#endif // not C++11

#if defined(_MSC_VER) && _MSC_VER < 1900 // Microsoft Visual Studio < 2015

#include <cfloat>
#define isnan(x) _isnan(x)

#endif

#endif

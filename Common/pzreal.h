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

extern int gPrintLevel;

enum MOptype {ESum,EProd,EDiv,ESqrt,EPow,ECos,ESin,EAcos,EAsin,EAtan,EExp,ELog};
const int gNumOp = 12;

/// This class implements a counter
/**
By modifying the definition of the type of REAL, the number of floating point operations can be counted
*/
struct TPZCounter {
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

namespace std {

/// This class implements floating point number associated with a counter
/**
By modifying the definition of the type of REAL, the number of floating point operations can be counted
*/
class TPZFlopCounter {
public:
    double fVal;
static TPZCounter gCount;
    
inline TPZFlopCounter()
{
}
inline TPZFlopCounter(const double &val)
{
  fVal = val;
}
inline TPZFlopCounter operator*(const TPZFlopCounter &oth) const
{
  TPZFlopCounter result;
  result.fVal = fVal*oth.fVal;
  gCount.fCount[EProd]++;
  return result;
}

inline TPZFlopCounter operator/(const TPZFlopCounter &oth) const
{
  TPZFlopCounter result;
  result.fVal = fVal/oth.fVal;
  gCount.fCount[EDiv]++;
  return result;
}

inline TPZFlopCounter operator+(const TPZFlopCounter &oth) const
{
  TPZFlopCounter result;
  result.fVal = fVal+oth.fVal;
  gCount.fCount[ESum]++;
  return result;
}

inline TPZFlopCounter operator-(const TPZFlopCounter &oth) const
{
  TPZFlopCounter result;
  result.fVal = fVal-oth.fVal;
  gCount.fCount[ESum]++;
  return result;
}

inline TPZFlopCounter &operator+=(const TPZFlopCounter &oth) 
{
  fVal += oth.fVal;
  gCount.fCount[ESum]++;
  return *this;
}

inline TPZFlopCounter &operator-=(const TPZFlopCounter &oth) 
{
  fVal -= oth.fVal;
  gCount.fCount[ESum]++;
  return *this;
}

inline TPZFlopCounter &operator*=(const TPZFlopCounter &oth) 
{
  fVal *= oth.fVal;
  gCount.fCount[EProd]++;
  return *this;
}

inline TPZFlopCounter &operator/=(const TPZFlopCounter &oth) 
{
  fVal /= oth.fVal;
  gCount.fCount[EDiv]++;
  return *this;
}

inline TPZFlopCounter operator-() const
{
  TPZFlopCounter result;
  result.fVal = -fVal;
  gCount.fCount[ESum]++;
  return result;
}

inline TPZFlopCounter operator+() const
{
  return *this;
}

inline bool operator<(const TPZFlopCounter &sec) const
{
  return (fVal < sec.fVal);
}
inline bool operator>(const TPZFlopCounter &sec) const
{
  return (fVal > sec.fVal);
}
inline bool operator<=(const TPZFlopCounter &sec) const
{
  return (fVal <= sec.fVal);
}
inline bool operator>=(const TPZFlopCounter &sec) const
{
  return (fVal >= sec.fVal);
}
inline bool operator==(const TPZFlopCounter &sec) const
{
  return (fVal == sec.fVal);
}
inline bool operator!=(const TPZFlopCounter &sec) const
{
  return (fVal != sec.fVal);
}

friend TPZFlopCounter sqrt(const TPZFlopCounter &orig);
friend TPZFlopCounter pow(const TPZFlopCounter &orig,const TPZFlopCounter &xp);
friend TPZFlopCounter fabs(const TPZFlopCounter &orig);
friend TPZFlopCounter acos(const TPZFlopCounter &orig);
friend TPZFlopCounter cos(const TPZFlopCounter &orig);
friend TPZFlopCounter asin(const TPZFlopCounter &orig);
friend TPZFlopCounter sin(const TPZFlopCounter &orig);
friend TPZFlopCounter atan(const TPZFlopCounter &orig);
friend TPZFlopCounter atan2(const TPZFlopCounter &val1,const TPZFlopCounter &val2);
friend TPZFlopCounter exp(const TPZFlopCounter &val);
friend TPZFlopCounter log(const TPZFlopCounter &val);
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
  //TPZFlopCounter::gCount.fCount[ESqrt]++;
  return result;
}

inline TPZFlopCounter pow(const TPZFlopCounter &orig,const TPZFlopCounter &xp)
{
  TPZFlopCounter result;
  result.fVal = ::pow(orig.fVal,xp.fVal);
  TPZFlopCounter::gCount.fCount[EPow]++;
  return result;
}

inline TPZFlopCounter operator*(double val1, const TPZFlopCounter &val2)
{
  TPZFlopCounter result;
  result = TPZFlopCounter(val1)*val2;
//  TPZFlopCounter::gCount.fCount[EProd]++;
  return result;
  
}

inline TPZFlopCounter operator/(double val1, const TPZFlopCounter &val2)
{
  TPZFlopCounter result;
  result = TPZFlopCounter(val1)/val2;
//  TPZFlopCounter::gCount.fCount[EProd]++;
  return result;
  
}

inline TPZFlopCounter operator+(double val1, const TPZFlopCounter &val2)
{
  TPZFlopCounter result;
  result = TPZFlopCounter(val1)+val2;
//  TPZFlopCounter::gCount.fCount[EProd]++;
  return result;
  
}

inline TPZFlopCounter operator+(const TPZFlopCounter &val2, double val1 )
{
  TPZFlopCounter result;
  result = TPZFlopCounter(val1)+val2;
//  TPZFlopCounter::gCount.fCount[EProd]++;
  return result;
  
}

inline TPZFlopCounter operator-(double val1, const TPZFlopCounter &val2)
{
  TPZFlopCounter result;
  result = TPZFlopCounter(val1)-val2;
//  TPZFlopCounter::gCount.fCount[EProd]++;
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

inline std::ostream &operator<<(std::ostream &out, const TPZFlopCounter &val)
{
  return out << val.fVal;
}
inline std::istream &operator>>(std::istream &out, const TPZFlopCounter &val)
{
  return out >> val.fVal;
}

}

#endif

#ifdef contar
typedef std::TPZFlopCounter REAL;
//typedef std::TPZFlopCounter* REALPtr;
#else
/// This is the type of floating point number PZ will use
typedef double REAL;
//typedef double *REALPtr;
#endif

#endif

//
//  ExtractVal.h
//  PZ
//
//  Created by Francisco Orlandini on 10/21/15.
//
//

#ifndef ExtractVal_h
#define ExtractVal_h

class TPZExtractVal
{
public:
  
  static REAL val( const int number)
  {
    return (REAL)number;
  }
  static REAL val( const int64_t number)
  {
    return (REAL)number;
  }
  static REAL val( const float number)
  {
    return (REAL)number;
  }
  static REAL val( const double number)
  {
    return (REAL)number;
  }
  static REAL val( const long double number)
  {
    return (REAL)number;
  }
  static REAL val( const std::complex<float> number)
  {
    return (REAL)number.real();
  }
  static REAL val( const std::complex<double> number)
  {
    return (REAL)number.real();
  }
  static REAL val( const std::complex<long double> number)
  {
    return (REAL)number.real();
  }
  
  template <class T>
  static REAL val(const T number)
  {
    return TPZExtractVal::val( number.val() ); // recursively downgrading until REAL type is reached
  }
  
  template<class T>
  static bool IsZero( const T & a ){
    return ::IsZero(TPZExtractVal::val(a));
  }
  
};


#endif /* ExtractVal_h */

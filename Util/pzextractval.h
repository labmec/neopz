//
//  ExtractVal.h
//  PZ
//
//  Created by Francisco Orlandini on 10/21/15.
//
//

#ifndef ExtractVal_h
#define ExtractVal_h

#include "pzreal.h"

#include "tfad.h"

class TPZExtractVal {
public:

    static REAL val(const int number) {
        return (REAL) number;
    }

    static REAL val(const int64_t number) {
        return (REAL) number;
    }

    static REAL val(const float number) {
        return (REAL) number;
    }

    static REAL val(const double number) {
        return (REAL) number;
    }

    static REAL val(const long double number) {
        return (REAL) number;
    }

    static REAL val(const std::complex<float> number) {
        return (REAL) number.real();
    }

    static REAL val(const std::complex<double> number) {
        return (REAL) number.real();
    }

    static REAL val(const std::complex<long double> number) {
        return (REAL) number.real();
    }

    template <class T>
    static REAL val(const T number) {
        return TPZExtractVal::val(number.val()); // recursively downgrading until REAL type is reached
    }

    static int& ref(int &number) {
        return number;
    }

    static int64_t& ref(int64_t &number) {
        return number;
    }

    static float& ref(float &number) {
        return number;
    }

    static double& ref(double &number) {
        return number;
    }

    static long double& ref(long double &number) {
        return number;
    }


    template <int Num, class T>
    static typename TFad<Num,T>::arithmetic_type& ref(TFad<Num,T> &number) {
        return TPZExtractVal::ref(number.val());
    }
    

    template<class T>
    static bool IsZero(const T & a) {
        return ::IsZero(TPZExtractVal::val(a));
    }

};


#endif /* ExtractVal_h */

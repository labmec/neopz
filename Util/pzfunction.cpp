/** 
 * @file 
 * @brief Contains the implementation of the methods to TPZFunction class. 
 */

#include "pzfunction.h"

template<>
int TPZFunction<float>::ClassId() {
    return Hash("TPZFunction") ^ Hash("float");
}

template<>
int TPZFunction<double>::ClassId() {
    return Hash("TPZFunction") ^ Hash("double");
}

template<>
int TPZFunction<long double>::ClassId() {
    return Hash("TPZFunction") ^ Hash("long double");
}

template<>
int TPZFunction<std::complex<float> >::ClassId() {
    return Hash("TPZFunction") ^ Hash("std::complex<float>");
}

template<>
int TPZFunction<std::complex<double> >::ClassId() {
    return Hash("TPZFunction") ^ Hash("std::complex<double>");
}

template<>
int TPZFunction<std::complex<long double> >::ClassId() {
    return Hash("TPZFunction") ^ Hash("std::complex<long double>");
}

template class TPZFunction<float>;
template class TPZFunction<double>;
template class TPZFunction<long double>;

template class TPZFunction<std::complex<float> >;
template class TPZFunction<std::complex<double> >;
template class TPZFunction<std::complex<long double> >;

template class TPZDummyFunction<float>;
template class TPZDummyFunction<double>;
template class TPZDummyFunction<long double>;

template class TPZDummyFunction<std::complex<float> >;
template class TPZDummyFunction<std::complex<double> >;
template class TPZDummyFunction<std::complex<long double> >;

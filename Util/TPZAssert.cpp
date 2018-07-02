/* 
 * File:   TPZAssert.cpp
 * Author: thiago
 * 
 * Created on 26 de Fevereiro de 2018, 13:49
 */

#include "TPZAssert.h"

template<> int& TPZAssert::NonNegative(int&);
template<> int64_t& TPZAssert::NonNegative(int64_t&);
template<> float& TPZAssert::NonNegative(float&);
template<> double& TPZAssert::NonNegative(double&);
template<> long double& TPZAssert::NonNegative(long double&);

template<> int TPZAssert::NonNegative(const int&);
template<> int64_t TPZAssert::NonNegative(const int64_t&);
template<> float TPZAssert::NonNegative(const float&);
template<> double TPZAssert::NonNegative(const double&);
template<> long double TPZAssert::NonNegative(const long double&);
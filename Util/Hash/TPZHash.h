/* 
 * File:   Hash.h
 * Author: quinelato
 *
 * Created on September 25, 2017, 12:23 PM
 */

#ifndef TPZHASH_H
#define TPZHASH_H

#include "pzreal.h"
#include "MurmurHash3.h"
#include <string>
#include <cstdint>

int32_t Hash(std::string str);

template <typename T>
int ClassIdOrHash(){
    return T().ClassId();
}

template <>
int ClassIdOrHash<TPZFlopCounter>();

template <>
int ClassIdOrHash<int>();

template <>
int ClassIdOrHash<long int>();

template <>
int ClassIdOrHash<long long>();

#ifndef __linux__
template <>
int ClassIdOrHash<int64_t>();
#endif

template <>
int ClassIdOrHash<u_int64_t>();

template <>
int ClassIdOrHash<float>();

template <>
int ClassIdOrHash<double>();

template <>
int ClassIdOrHash<long double>();

template <>
int ClassIdOrHash<std::complex<float>>();

template <>
int ClassIdOrHash<std::complex<double>>();

template <>
int ClassIdOrHash<std::complex<long double>>();


#endif /* TPZHASH_H */


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
#include <type_traits>

int32_t Hash(std::string str);

template <typename T>
typename std::enable_if<!std::is_pointer<T>::value && std::is_abstract<T>::value, int>::type ClassIdOrHash(){
    return T::StaticClassId();
}

template <typename T>
typename std::enable_if<!std::is_pointer<T>::value && !std::is_abstract<T>::value, int>::type ClassIdOrHash(){
    return T().ClassId();
}

template <typename T>
typename std::enable_if<std::is_pointer<T>::value, int>::type ClassIdOrHash(){
    return Hash("pointer") ^ ClassIdOrHash<typename std::remove_pointer<T>::type>() << 1;
}

//template <typename T>
//int ClassIdOrHash<typename std::enable_if<std::is_pointer<T>::value, T>::type>() {
//    return Hash("pointer") ^ ClassIdOrHash<T>() << 1;
//}

template <>
int ClassIdOrHash<TPZFlopCounter>();

template <>
int ClassIdOrHash<int>();

template <>
int ClassIdOrHash<long int>();

template <>
int ClassIdOrHash<long long>();

template <>
int ClassIdOrHash<uint64_t>();

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


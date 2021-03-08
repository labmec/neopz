/**
 * @file FadUnitTest.cpp
 * @brief Define a Unit Test using Boost for testing consistency of FAD
 * operations regardless of data types
 *
 */

#include "fad.h"
#include <math.h>

// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz fad_tests tests

#include "boost/test/tools/output_test_stream.hpp"
#include "boost/test/unit_test.hpp"

//#define NOISY //outputs operations' results

const long double tol = 1e-7;

template <typename T, typename U> void FadVsFadTest() {
    Fad<T> var_1(0.);
    Fad<T> var_2(0.);
    Fad<U> var_3(0.);
    Fad<U> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
#ifdef NOISY
    std::cout<< "op result (expected 13): "<< result.val() <<std::endl;
#endif
    BOOST_CHECK(checkResult);
    return;
}

template <typename T, typename U> void FadVsArithmeticTest() {

    Fad<T> var_1(0.);
    U var_2(0.);
    Fad<T> var_3(0.);
    Fad<T> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
#ifdef NOISY
    std::cout<< "op result (expected 13): "<< result.val() <<std::endl;
#endif
    BOOST_CHECK(checkResult);
    return;
}
template <typename T, typename U> void FadVsTempTest() {
    Fad<T> var_1(0.);
    Fad<T> var_3(0.);
    Fad<T> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result =
        var_1 * ((U)(2.0) * (var_1 + (U)(2.0) + (var_1 / (U)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
#ifdef NOISY
    std::cout<< "op result (expected 13): "<<result.val() <<std::endl;
#endif
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_SUITE(fad_fad_tests)

BOOST_AUTO_TEST_CASE(fad_float_fad_float_tests) {

    FadVsFadTest<float, float>();
}

BOOST_AUTO_TEST_CASE(fad_float_fad_double_tests) {
    
    FadVsFadTest<float, double>();
}

BOOST_AUTO_TEST_CASE(fad_float_fad_long_double_tests) {
    
    FadVsFadTest<float, long double>();
}

BOOST_AUTO_TEST_CASE(fad_double_fad_float_tests) {

    FadVsFadTest<double, float>();
}

BOOST_AUTO_TEST_CASE(fad_double_fad_double_tests) {
    
    FadVsFadTest<double, double>();
}

BOOST_AUTO_TEST_CASE(fad_double_fad_long_double_tests) {

    FadVsFadTest<double, long double>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_fad_float_tests) {
    
    FadVsFadTest<long double, float>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_fad_double_tests) {
    
    FadVsFadTest<long double, double>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_fad_long_double_tests) {
    
    FadVsFadTest<long double, long double>();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(fad_arithmetic_tests)

BOOST_AUTO_TEST_CASE(fad_float_float_tests) {

    FadVsArithmeticTest<float, float>();
}

BOOST_AUTO_TEST_CASE(fad_float_double_tests) {

    FadVsArithmeticTest<float, double>();
}

BOOST_AUTO_TEST_CASE(fad_float_long_double_tests) {

    FadVsArithmeticTest<float, long double>();
}

BOOST_AUTO_TEST_CASE(fad_double_float_tests) {

    FadVsArithmeticTest<double, float>();
}

BOOST_AUTO_TEST_CASE(fad_double_double_tests) {

    FadVsArithmeticTest<double, double>();
}

BOOST_AUTO_TEST_CASE(fad_double_long_double_tests) {

    FadVsArithmeticTest<double, long double>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_float_tests) {

    FadVsArithmeticTest<long double, float>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_double_tests) {

    FadVsArithmeticTest<long double, double>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_long_double_tests) {

    FadVsArithmeticTest<long double, long double>();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(fad_temp_tests)

BOOST_AUTO_TEST_CASE(fad_float_float_tests) { 
    
    FadVsTempTest<float, float>();
}

BOOST_AUTO_TEST_CASE(fad_float_double_tests) { 
    
    FadVsTempTest<float, double>(); 
}

BOOST_AUTO_TEST_CASE(fad_float_long_double_tests) {

    FadVsTempTest<float, long double>();
}

BOOST_AUTO_TEST_CASE(fad_double_float_tests) { 
    
    FadVsTempTest<double, float>(); 
}

BOOST_AUTO_TEST_CASE(fad_double_double_tests) {

    FadVsTempTest<double, double>();
}

BOOST_AUTO_TEST_CASE(fad_double_long_double_tests) {

    FadVsTempTest<double, long double>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_float_tests) {
    
    FadVsTempTest<long double, float>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_double_tests) {
    
    FadVsTempTest<long double, double>();
}

BOOST_AUTO_TEST_CASE(fad_long_double_long_double_tests) {
    
    FadVsTempTest<long double, long double>();
}

BOOST_AUTO_TEST_SUITE_END()

#endif

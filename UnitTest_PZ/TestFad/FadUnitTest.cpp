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

#include "boost/test/floating_point_comparison.hpp"
#include "boost/test/output_test_stream.hpp"
#include "boost/test/unit_test.hpp"

//#define NOISY //outputs operations' results

const long double tol = 1e-7;

BOOST_AUTO_TEST_SUITE(fad_fad_tests)

BOOST_AUTO_TEST_CASE(fad_float_fad_float_tests) {

    Fad<float> var_1(0.);
    Fad<float> var_2(0.);
    Fad<float> var_3(0.);
    Fad<float> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_fad_double_tests) {

    Fad<double> var_1(0.);
    Fad<double> var_2(0.);
    Fad<double> var_3(0.);
    Fad<double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_long_double_fad_long_double_tests) {

    Fad<long double> var_1(0.);
    Fad<long double> var_2(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_fad_float_tests) {

    Fad<float> var_1(0.);
    Fad<float> var_2(0.);
    Fad<double> var_3(0.);
    Fad<double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_fad_long_double_tests) {

    Fad<double> var_1(0.);
    Fad<double> var_2(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_float_fad_long_double_tests) {

    Fad<float> var_1(0.);
    Fad<float> var_2(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(fad_arithmetic_tests)

BOOST_AUTO_TEST_CASE(fad_float_float_tests) {

    Fad<float> var_1(0.);
    float var_2(0.);
    Fad<float> var_3(0.);
    Fad<float> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_float_double_tests) {

    Fad<float> var_1(0.);
    double var_2(0.);
    Fad<float> var_3(0.);
    Fad<float> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_float_long_double_tests) {

    Fad<float> var_1(0.);
    long double var_2(0.);
    Fad<float> var_3(0.);
    Fad<float> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_float_tests) {

    Fad<double> var_1(0.);
    float var_2(0.);
    Fad<double> var_3(0.);
    Fad<double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_double_tests) {

    Fad<double> var_1(0.);
    double var_2(0.);
    Fad<double> var_3(0.);
    Fad<double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_long_double_tests) {

    Fad<double> var_1(0.);
    long double var_2(0.);
    Fad<double> var_3(0.);
    Fad<double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_long_double_float_tests) {

    Fad<long double> var_1(0.);
    float var_2(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_long_double_double_tests) {

    Fad<long double> var_1(0.);
    double var_2(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_long_double_long_double_tests) {

    Fad<long double> var_1(0.);
    long double var_2(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_2 = 2.;
    var_3 = 3.;
    result = var_1 * (var_2 * (var_1 + var_2 + (var_1 / var_2 + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(fad_temp_tests)

BOOST_AUTO_TEST_CASE(fad_float_float_tests) {

    Fad<float> var_1(0.);
    Fad<float> var_3(0.);
    Fad<float> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result = var_1 * ((float)(2.0) *
                      (var_1 + (float)(2.0) + (var_1 / (float)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_float_double_tests) {

    Fad<float> var_1(0.);
    Fad<float> var_3(0.);
    Fad<float> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result = var_1 * ((double)(2.0) * (var_1 + (double)(2.0) +
                                       (var_1 / (double)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_float_long_double_tests) {

    Fad<float> var_1(0.);
    Fad<float> var_3(0.);
    Fad<float> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result =
        var_1 * ((long double)(2.0) * (var_1 + (long double)(2.0) +
                                       (var_1 / (long double)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_float_tests) {

    Fad<double> var_1(0.);
    Fad<double> var_3(0.);
    Fad<double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result = var_1 * ((float)(2.0) *
                      (var_1 + (float)(2.0) + (var_1 / (float)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_double_tests) {

    Fad<double> var_1(0.);
    Fad<double> var_3(0.);
    Fad<double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result = var_1 * ((double)(2.0) * (var_1 + (double)(2.0) +
                                       (var_1 / (double)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_double_long_double_tests) {

    Fad<double> var_1(0.);
    Fad<double> var_3(0.);
    Fad<double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result =
        var_1 * ((long double)(2.0) * (var_1 + (long double)(2.0) +
                                       (var_1 / (long double)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_long_double_float_tests) {

    Fad<long double> var_1(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result = var_1 * ((float)(2.0) *
                      (var_1 + (float)(2.0) + (var_1 / (float)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_long_double_double_tests) {

    Fad<long double> var_1(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result = var_1 * ((double)(2.0) * (var_1 + (double)(2.0) +
                                       (var_1 / (double)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_CASE(fad_long_double_long_double_tests) {

    Fad<long double> var_1(0.);
    Fad<long double> var_3(0.);
    Fad<long double> result(0.);
    bool checkResult = false;

    var_1 = 1.;
    var_3 = 3.;
    result =
        var_1 * ((long double)(2.0) * (var_1 + (long double)(2.0) +
                                       (var_1 / (long double)(2.0) + var_3)));
    checkResult = std::abs(result.val() - 13.) < tol;
    BOOST_CHECK(checkResult);
    return;
}

BOOST_AUTO_TEST_SUITE_END()

#endif

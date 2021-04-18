/**
 * @file FadUnitTest.cpp
 * @brief Define a Unit Test using Boost for testing consistency of FAD
 * operations regardless of data types
 *
 */

#include "fad.h"
#include <math.h>

#include <catch2/catch.hpp>
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
    REQUIRE(checkResult);
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
    REQUIRE(checkResult);
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
    REQUIRE(checkResult);
    return;
}


TEST_CASE("fad_float_fad_float_tests","[arithmetic_tests]") {

    FadVsFadTest<float, float>();
}

TEST_CASE("fad_float_fad_double_tests","[arithmetic_tests]") {
    
    FadVsFadTest<float, double>();
}

TEST_CASE("fad_float_fad_long_double_tests","[arithmetic_tests]") {
    
    FadVsFadTest<float, long double>();
}

TEST_CASE("fad_double_fad_float_tests","[arithmetic_tests]") {

    FadVsFadTest<double, float>();
}

TEST_CASE("fad_double_fad_double_tests","[arithmetic_tests]") {
    
    FadVsFadTest<double, double>();
}

TEST_CASE("fad_double_fad_long_double_tests","[arithmetic_tests]") {

    FadVsFadTest<double, long double>();
}

TEST_CASE("fad_long_double_fad_float_tests","[arithmetic_tests]") {
    
    FadVsFadTest<long double, float>();
}

TEST_CASE("fad_long_double_fad_double_tests","[arithmetic_tests]") {
    
    FadVsFadTest<long double, double>();
}

TEST_CASE("fad_long_double_fad_long_double_tests","[arithmetic_tests]") {
    
    FadVsFadTest<long double, long double>();
}


TEST_CASE("fad_float_float_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<float, float>();
}

TEST_CASE("fad_float_double_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<float, double>();
}

TEST_CASE("fad_float_long_double_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<float, long double>();
}

TEST_CASE("fad_double_float_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<double, float>();
}

TEST_CASE("fad_double_double_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<double, double>();
}

TEST_CASE("fad_double_long_double_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<double, long double>();
}

TEST_CASE("fad_long_double_float_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<long double, float>();
}

TEST_CASE("fad_long_double_double_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<long double, double>();
}

TEST_CASE("fad_long_double_long_double_tests","[fad_arithmetic_tests]") {

    FadVsArithmeticTest<long double, long double>();
}


TEST_CASE("tmp_fad_float_float_tests","[fad_temp_tests]") { 
    
    FadVsTempTest<float, float>();
}

TEST_CASE("tmp_fad_float_double_tests","[fad_temp_tests]") { 
    
    FadVsTempTest<float, double>(); 
}

TEST_CASE("tmp_fad_float_long_double_tests","[fad_temp_tests]") {

    FadVsTempTest<float, long double>();
}

TEST_CASE("tmp_fad_double_float_tests","[fad_temp_tests]") { 
    
    FadVsTempTest<double, float>(); 
}

TEST_CASE("tmp_fad_double_double_tests","[fad_temp_tests]") {

    FadVsTempTest<double, double>();
}

TEST_CASE("tmp_fad_double_long_double_tests","[fad_temp_tests]") {

    FadVsTempTest<double, long double>();
}

TEST_CASE("tmp_fad_long_double_float_tests","[fad_temp_tests]") {
    
    FadVsTempTest<long double, float>();
}

TEST_CASE("tmp_fad_long_double_double_tests","[fad_temp_tests]") {
    
    FadVsTempTest<long double, double>();
}

TEST_CASE("tmp_fad_long_double_long_double_tests","[fad_temp_tests]") {
    
    FadVsTempTest<long double, long double>();
}

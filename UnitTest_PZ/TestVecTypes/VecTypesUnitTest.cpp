/**
 * @file
 * @brief Contains Unit Tests for TPZVec<T> and TPZManVector<T>
 */

#include "pzmanvector.h"

#include <catch2/catch.hpp>


TEST_CASE("tpzvec_tests","[test_vectypes]") {
    //mv ctor
    TPZVec<REAL> vec1(std::move(TPZVec<REAL>(3,0.)));
    REQUIRE(true);
    //cp ctor
    TPZVec<REAL> vec2(vec1);
    REQUIRE(true);
    //mv assignment
    vec2= std::move(vec1);
    REQUIRE(true);
    //mv assignment with tpzmanvector
    TPZVec<REAL> vec3(std::move(TPZManVector<REAL,3>(3,0)));
    REQUIRE(true);
}


TEST_CASE("tpzmanvector_tests","[test_vectypes]") {

    //ctor and assignment tests
    //tpzvec
    {
        TPZVec<REAL> vec(2,0.);
        TPZManVector<REAL,3> manvec11(vec);
        TPZManVector<REAL,3> manvec12(TPZVec<REAL>(2,0.));
        REQUIRE(true);
        vec.Resize(3);
        TPZManVector<REAL,3> manvec21(vec);
        TPZManVector<REAL,3> manvec22(TPZVec<REAL>(3,0.));
        REQUIRE(true);
        vec.Resize(5);
        TPZManVector<REAL,3> manvec31(vec);
        TPZManVector<REAL,3> manvec32(TPZVec<REAL>(5,0.));
        REQUIRE(true);
        try{
            TPZManVector<REAL,3> manvec = std::move(vec);
        }catch(...){
            REQUIRE(false);
        }
        REQUIRE(true);
    }
    //N1==N2
    {
        TPZManVector<REAL,3> vec(2,0.);
        TPZManVector<REAL,3> manvec1(vec);
        REQUIRE(true);
        vec.Resize(3);
        TPZManVector<REAL,3> manvec2(vec);
        REQUIRE(true);
        vec.Resize(5);
        TPZManVector<REAL,3> manvec3(vec);
        TPZManVector<REAL,3> manvec4(std::move(TPZManVector<REAL,3>(2,0.)));
        TPZManVector<REAL,3> manvec5(std::move(TPZManVector<REAL,3>(3,0.)));
        TPZManVector<REAL,3> manvec6(std::move(TPZManVector<REAL,3>(4,0.)));
        REQUIRE(true);
        TPZManVector<REAL,3> manvec7(3,1);
        manvec7 = std::move(vec);
    }
    //N1>N2
    {
        TPZManVector<REAL,2> vec(2,0.);
        TPZManVector<REAL,3> manvec1(vec);
        REQUIRE(true);
        vec.Resize(3);
        TPZManVector<REAL,3> manvec2(vec);
        REQUIRE(true);
        vec.Resize(5);
        TPZManVector<REAL,3> manvec3(vec);
        TPZManVector<REAL,3> manvec4(TPZManVector<REAL,2>(1,0.));
        TPZManVector<REAL,3> manvec5(TPZManVector<REAL,2>(2,0.));
        TPZManVector<REAL,3> manvec6(TPZManVector<REAL,2>(3,0.));
        TPZManVector<REAL,3> manvec7(TPZManVector<REAL,2>(4,0.));
        REQUIRE(true);
    }
    //N1<N2
    {
        TPZManVector<REAL,5> vec(2,0.);
        TPZManVector<REAL,3> manvec1(vec);
        REQUIRE(true);
        vec.Resize(3);
        TPZManVector<REAL,3> manvec2(vec);
        REQUIRE(true);
        vec.Resize(5);
        TPZManVector<REAL,3> manvec3(vec);
        TPZManVector<REAL,3> manvec4(TPZManVector<REAL,5>(4,0.));
        TPZManVector<REAL,3> manvec5(TPZManVector<REAL,5>(5,0.));
        TPZManVector<REAL,3> manvec6(TPZManVector<REAL,5>(6,0.));
        REQUIRE(true);
    }
}
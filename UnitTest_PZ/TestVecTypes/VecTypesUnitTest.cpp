/**
 * @file
 * @brief Contains Unit Tests for TPZVec<T> and TPZManVector<T>
 */

#include "pzmanvector.h"
#include "fad.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>


TEMPLATE_TEST_CASE("tpzvec_tests","[test_vectypes]", REAL, Fad<REAL>) {
    const TestType val = (TestType) 0;
    //tpzvec  tpzvec tests
    {
        //mv ctor
        TPZVec<TestType> vec1(std::move(TPZVec<TestType>(3,val)));
        for(auto &v : vec1) {v+=1;}//let us test if the memory is accessible
        REQUIRE(true);
        //cp ctor
        TPZVec<TestType> vec2(vec1);
        for(auto &v : vec2) {v+=1;}//let us test if the memory is accessible
        REQUIRE(true);
        //mv assignment
        vec2= std::move(vec1);
        for(auto &v : vec2) {v+=1;}//let us test if the memory is accessible
        REQUIRE(true);
    }
    //tpzvec tpzmanvetor tests
    {
        //mv ctor with fExtAlloc = fStore
        TPZVec<TestType> vec1(std::move(TPZManVector<TestType,5>(3,val)));
        for(auto &v : vec1) {v+=1;}//let us test if the memory is accessible
        REQUIRE(true);
        //mv ctor with fExtAlloc != fStore
        TPZVec<TestType> vec2(std::move(TPZManVector<TestType,5>(10,val)));
        for(auto &v : vec2) {v+=1;}//let us test if the memory is accessible
        REQUIRE(true);
        //mv assignment with fExtAlloc = fStore
        TPZManVector<TestType,3> manvec1(3,val);
        TPZVec<TestType> vec3 = std::move(manvec1);
        for(auto &v : vec3) {v+=1;}//let us test if the memory is accessible
        REQUIRE(true);
        //mv assignment with fExtAlloc != fStore
        TPZManVector<TestType,1> manvec2(3,val);
        TPZVec<TestType> vec4 = std::move(manvec2);
        for(auto &v : vec4) {v+=1;}//let us test if the memory is accessible
        REQUIRE(true);
    }
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

TEST_CASE("mixedvecs_tests","[test_vectypes]") {

    auto myfunc = []() -> TPZVec<int> {
        //this vector has not allocated any memory
        TPZManVector<int,3> myvec = {1,2,3};
        return std::move(myvec);
    };


    //move ctor
    try{
        TPZVec<int> myextvec = myfunc();
        myextvec[0]++;//just to see if memory is still valid
    }catch(...){
        REQUIRE(false);
    }
    REQUIRE(true);

    
    //move assignment operator
    try{
        TPZManVector<int,3> manvec = {1,2,3};
        TPZVec<int> normalvec = {4,5,6};
        normalvec = std::move(manvec);
        normalvec[0]++;//just to see if memory is still valid
    }catch(...){
        REQUIRE(false);
    }
    REQUIRE(true);
}
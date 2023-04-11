/**
 * @file AutoPointerTest.cpp
 * @brief Define a Unit Test for TPZAutoPointer
 *
 */

#include "pzfmatrix.h"

#include <thread>
#include <chrono>
#include <catch2/catch.hpp>

namespace autoptrtest{
    //! stress test for autopointer
    void AutoPointerStressTest();
    //! conversion tests
    void AutoPointerConversionTest();
    //! comparison tests
    void AutoPointerComparisonTest();
};


TEST_CASE("stress_test","[autoptr_test]"){
    autoptrtest::AutoPointerStressTest();
}

TEST_CASE("conversion_test","[autoptr_test]"){
    autoptrtest::AutoPointerConversionTest();
}

TEST_CASE("comparison_test","[autoptr_test]"){
    autoptrtest::AutoPointerComparisonTest();
}

namespace autoptrtest{
    void ThreadTask(TPZAutoPointer<TPZFMatrix<REAL>> object) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
};

void autoptrtest::AutoPointerStressTest()
{
    const int nthreads = 10;
    std::thread allthreads[nthreads];
    {
        TPZAutoPointer<TPZFMatrix<REAL> > matrix = new TPZFMatrix<REAL>(10,10,1.);
        for(int i = 0; i<nthreads; i++)
        {
            allthreads[i] = std::thread(ThreadTask, matrix);
        }
    }
    for (int i=0; i<nthreads; i++) {
        allthreads[i].join();
    }
}

void autoptrtest::AutoPointerConversionTest(){

    SECTION("Testing TPZAutoPointerDynamicCast"){
        TPZAutoPointer<TPZMatrix<STATE>> mat_ptr(new TPZFMatrix<STATE>(1,1,0.));
        TPZAutoPointer<TPZFMatrix<STATE>> fmat_ptr =
            TPZAutoPointerDynamicCast<TPZFMatrix<STATE>>(mat_ptr);
        REQUIRE(fmat_ptr);
        TPZAutoPointer<TPZFNMatrix<9,STATE>> fnmat_ptr =
            TPZAutoPointerDynamicCast<TPZFNMatrix<9,STATE>>(mat_ptr);
        REQUIRE(!fnmat_ptr);
    }
    SECTION("Testing Conversion Constructor"){
        TPZAutoPointer<TPZFMatrix<STATE>> fmat_ptr(new TPZFMatrix<STATE>(1,1,0.));
        TPZAutoPointer<TPZMatrix<STATE>> mat_ptr(fmat_ptr);
        REQUIRE(mat_ptr);
    }
    SECTION("Testing Conversion Assignment"){
        TPZAutoPointer<TPZFMatrix<STATE>> fmat_ptr(new TPZFMatrix<STATE>(1,1,0.));
        TPZAutoPointer<TPZMatrix<STATE>> mat_ptr = nullptr;
        mat_ptr = fmat_ptr;
        REQUIRE(mat_ptr);
    }
}


void autoptrtest::AutoPointerComparisonTest(){
    /*
      In order for these tests to run as required, double
      parentheses are required, otherwise the
      TPZAutoPointer<T>::operator== will not be called
      See https://github.com/catchorg/Catch2/blob/devel/docs/assertions.md
     */
    
    SECTION("Testing nullptr nullptr comparison"){
        TPZAutoPointer<TPZFMatrix<STATE>> ap1, ap2;
        REQUIRE((ap1==ap2));
    }
    SECTION("Testing nullptr non-null comparison"){
        TPZAutoPointer<TPZFMatrix<STATE>> ap1, ap2;
        ap2 = new TPZFMatrix<STATE>(2,2,0.);
        REQUIRE((ap1!=ap2));
    }
    //same content, same pointer
    SECTION("Testing non-null non-null (equal ptr) comparison"){
        TPZAutoPointer<TPZFMatrix<STATE>> ap1, ap2;
        ap2 = new TPZFMatrix<STATE>(2,2,0.);
        ap1 = ap2;
        REQUIRE((ap1==ap2));
    }
    //same content, different pointers
    SECTION("Testing non-null non-null (diff ptr) comparison"){
        TPZAutoPointer<TPZFMatrix<STATE>> ap1, ap2;
        ap1 = new TPZFMatrix<STATE>(2,2,0.);
        ap2 = new TPZFMatrix<STATE>(2,2,0.);
        REQUIRE((ap1!=ap2));
    }
    //different content
    SECTION("Testing non-null non-null (diff) comparison"){
        TPZAutoPointer<TPZFMatrix<STATE>> ap1, ap2;
        ap1 = new TPZFMatrix<STATE>(1,1,0.);
        ap2 = new TPZFMatrix<STATE>(2,2,0.);
        REQUIRE((ap1!=ap2));
    }
}
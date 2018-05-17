/**
 * @file
 * @brief Contains Unit Tests for methods of the TPZChunkVector iterator.
 */

#include "pzadmchunk.h"

#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz iterator tests

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#endif

#ifdef USING_BOOST

#define EXP 1
#define N_CHUNKS 15
#define N_ELEM 30

void TestingFill() {
    TPZChunkVector<double, EXP> vector(N_CHUNKS);
    vector.Resize(N_ELEM);
    for (unsigned int i = 0; i < N_ELEM; ++i) {
        vector[i] = i;
    }
}

void ReadForeach(TPZChunkVector<double, EXP> &vector) {
    for (auto &elem : vector) {
        std::cout << elem << std::endl;
    }
}

void TestingFillReadForeach() {
    TPZChunkVector<double, EXP> vector(N_CHUNKS);
    vector.Resize(N_ELEM);
    int i = 0;
    for (auto &elem : vector) {
        elem = ++i;
    }

    ReadForeach(vector);
}

void ConstReadForeach(const TPZChunkVector<double, EXP> &vector) {
    for (auto &elem : vector) {
        std::cout << elem << std::endl;
    }
}

void TestingFillConstReadForeach() {
    TPZChunkVector<double, EXP> vector(N_CHUNKS);
    vector.Resize(N_ELEM);
    int i = 0;
    for (auto &elem : vector) {
        elem = ++i;
    }

    ConstReadForeach(vector);
}

BOOST_AUTO_TEST_SUITE(iterator_tests)


BOOST_AUTO_TEST_CASE(iterator_tests) {

    TestingFill();
    TestingFillReadForeach();
    TestingFillConstReadForeach();

}

BOOST_AUTO_TEST_SUITE_END()

#endif

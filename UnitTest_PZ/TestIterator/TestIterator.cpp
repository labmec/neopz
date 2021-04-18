/**
 * @file
 * @brief Contains Unit Tests for methods of the TPZChunkVector iterator.
 */

#include "pzadmchunk.h"
#include "catch2/catch.hpp"

#define EXP 2
#define N_CHUNKS 15
#define N_ELEM 30

void Fill() {
    TPZChunkVector<double, EXP> vector(N_CHUNKS);
    vector.Resize(N_ELEM);
    for (unsigned int i = 0; i < N_ELEM; ++i) {
        vector[i] = i;
    }
}

void TestingFill() {
    TPZChunkVector<double, EXP> vector(N_CHUNKS);
    vector.Resize(N_ELEM);
    unsigned int value = 0;
    for (auto &elem : vector) {
        elem = ++value;
    }
    for (unsigned int i = 0; i < N_ELEM; ++i) {
        REQUIRE(vector[i]== i+1);
    }
}

void ReadForeach(TPZChunkVector<double, EXP> &vector) {
    unsigned int i = 1;
    for (auto &elem : vector) {
        REQUIRE(elem== i++);
    }
}

void TestingFillReadForeach() {
    TPZChunkVector<double, EXP> vector(N_CHUNKS);
    vector.Resize(N_ELEM);
    unsigned int value = 0;
    for (auto &elem : vector) {
        elem = ++value;
    }

    ReadForeach(vector);
}

void ConstReadForeach(const TPZChunkVector<double, EXP> &vector) {
    unsigned int i = 1;
    for (auto &elem : vector) {
        REQUIRE(elem== i++);
    }
}

void TestingFillConstReadForeach() {
    TPZChunkVector<double, EXP> vector(N_CHUNKS);
    vector.Resize(N_ELEM);
    unsigned int value = 0;
    for (auto &elem : vector) {
        elem = ++value;
    }

    ConstReadForeach(vector);
}

TEST_CASE("iterator_tests_fill","[iterator_tests]") {
    TestingFill();
}
TEST_CASE("iterator_tests_fill_read","[iterator_tests]") {
    TestingFillReadForeach();
}
TEST_CASE("iterator_tests_const_read","[iterator_tests]") {
    TestingFillConstReadForeach();
}

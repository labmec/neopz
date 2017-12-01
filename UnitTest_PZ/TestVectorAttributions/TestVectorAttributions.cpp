//
//  TestVectorAttributions.cpp
//  PZ
//
//  Created by Thiago Quinelato on 1/dec/17.
//

#include "pzvec.h"
#include "pzmanvector.h"

#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz vector tests

#include "boost/test/output_test_stream.hpp"
#include "boost/test/unit_test.hpp"

TPZVec<int> BuildTPZVec(){
    TPZVec<int> vec(10);
    for (unsigned int i = 0; i < 10; ++i) {
        vec[i] = i;
    }
    return vec;
}

TPZManVector<int, 10> BuildTPZManVector(){
    TPZManVector<int, 10> vec(10);
    for (unsigned int i = 0; i < 10; ++i) {
        vec[i] = i;
    }
    return vec;
}

std::vector<int> BuildStdVector(){
    std::vector<int> vec(10);
    for (unsigned int i = 0; i < 10; ++i) {
        vec[i] = i;
    }
    return vec;
}

void TestingCopyTPZVecToTPZVec() {
    TPZVec<int> initialVec = BuildTPZVec();
    TPZVec<int> finalVec(initialVec);
    BOOST_ASSERT(finalVec.size() == 10);
    for (unsigned int i = 0; i < 10; ++i) {
        BOOST_ASSERT(finalVec[i] == i);
    }
}

void TestingMoveTPZVecToTPZVec() {
    TPZVec<int> finalVec(std::move(BuildTPZVec()));
    BOOST_ASSERT(finalVec.size() == 10);
    for (unsigned int i = 0; i < 10; ++i) {
        BOOST_ASSERT(finalVec[i] == i);
    }
}

void TestingCopyTPZManVectorToTPZManVector() {
    TPZManVector<int, 10> initialVec = BuildTPZManVector();
    TPZManVector<int, 10> finalVec(initialVec);
    BOOST_ASSERT(finalVec.size() == 10);
    for (unsigned int i = 0; i < 10; ++i) {
        BOOST_ASSERT(finalVec[i] == i);
    }
}

void TestingMoveTPZManVectorToTPZManVector() {
    TPZManVector<int, 10> finalVec(std::move(BuildTPZManVector()));
    BOOST_ASSERT(finalVec.size() == 10);
    for (unsigned int i = 0; i < 10; ++i) {
        BOOST_ASSERT(finalVec[i] == i);
    }
}

void TestingCopyTPZManVectorToTPZVec() {
    TPZManVector<int,10> initialVec = BuildTPZManVector();
    TPZVec<int> finalVec(initialVec);
    BOOST_ASSERT(finalVec.size() == 10);
    for (unsigned int i = 0; i < 10; ++i) {
        BOOST_ASSERT(finalVec[i] == i);
    }
}

void TestingMoveTPZManVectorToTPZVec() {
    TPZVec<int> finalVec(std::move(BuildTPZManVector()));
    BOOST_ASSERT(finalVec.size() == 10);
    for (unsigned int i = 0; i < 10; ++i) {
        BOOST_ASSERT(finalVec[i] == i);
    }
}
void TestingCopyStdVectorToTPZVec() {
    std::vector<int> initialVec = BuildStdVector();
    TPZVec<int> finalVec(initialVec);
    BOOST_ASSERT(finalVec.size() == 10);
    for (unsigned int i = 0; i < 10; ++i) {
        BOOST_ASSERT(finalVec[i] == i);
    }
}

BOOST_AUTO_TEST_SUITE(tpzvec_tests)

BOOST_AUTO_TEST_CASE(tpzvec_to_tpzvec_copy)
{
    TestingCopyTPZVecToTPZVec();
}

BOOST_AUTO_TEST_CASE(tpzvec_to_tpzvec_move)
{
    TestingMoveTPZVecToTPZVec();
}

BOOST_AUTO_TEST_CASE(tpzmanvector_to_tpzmanvector_copy)
{
    TestingCopyTPZManVectorToTPZManVector();
}

BOOST_AUTO_TEST_CASE(tpzmanvector_to_tpzmanvector_move)
{
    TestingMoveTPZManVectorToTPZManVector();
}

BOOST_AUTO_TEST_CASE(tpzmanvector_to_tpzvec_copy)
{
    TestingCopyTPZManVectorToTPZVec();
}

BOOST_AUTO_TEST_CASE(tpzmanvector_to_tpzvec_move)
{
    TestingMoveTPZManVectorToTPZVec();
}

BOOST_AUTO_TEST_CASE(stdvector_to_tpzvec_move)
{
    TestingCopyStdVectorToTPZVec();
}

BOOST_AUTO_TEST_SUITE_END()

#endif

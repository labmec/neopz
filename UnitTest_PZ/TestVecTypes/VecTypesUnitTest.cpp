/**
 * @file
 * @brief Contains Unit Tests for TPZVec<T> and TPZManVector<T>
 */

#include "pzmanvector.h"

#ifdef PZ_USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE pz test pzvec tests

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>



BOOST_AUTO_TEST_CASE(tpzvec_tests) {

    TPZVec<REAL> vec1(TPZVec<REAL>(3,0.));
    BOOST_CHECK(true);
    TPZVec<REAL> vec2(vec1);
    BOOST_CHECK(true);
    TPZVec<REAL> vec3 = std::move(vec2);
    BOOST_CHECK(true);
}


BOOST_AUTO_TEST_CASE(tpzmanvector_tests) {

    //ctor and assignment tests
    //tpzvec
    {
        TPZVec<REAL> vec(2,0.);
        TPZManVector<REAL,3> manvec11(vec);
        TPZManVector<REAL,3> manvec12(TPZVec<REAL>(2,0.));
        BOOST_CHECK(true);
        vec.Resize(3);
        TPZManVector<REAL,3> manvec21(vec);
        TPZManVector<REAL,3> manvec22(TPZVec<REAL>(3,0.));
        BOOST_CHECK(true);
        vec.Resize(5);
        TPZManVector<REAL,3> manvec31(vec);
        TPZManVector<REAL,3> manvec32(TPZVec<REAL>(5,0.));
        BOOST_CHECK(true);
        try{
            TPZManVector<REAL,3> manvec = std::move(vec);
        }catch(...){
            BOOST_CHECK(false);
        }
        BOOST_CHECK(true);
    }
    //N1==N2
    {
        TPZManVector<REAL,3> vec(2,0.);
        TPZManVector<REAL,3> manvec1(vec);
        BOOST_CHECK(true);
        vec.Resize(3);
        TPZManVector<REAL,3> manvec2(vec);
        BOOST_CHECK(true);
        vec.Resize(5);
        TPZManVector<REAL,3> manvec3(vec);
        TPZManVector<REAL,3> manvec4(TPZManVector<REAL,3>(2,0.));
        TPZManVector<REAL,3> manvec5(TPZManVector<REAL,3>(3,0.));
        TPZManVector<REAL,3> manvec6(TPZManVector<REAL,3>(4,0.));
        BOOST_CHECK(true);
    }
    //N1>N2
    {
        TPZManVector<REAL,2> vec(2,0.);
        TPZManVector<REAL,3> manvec1(vec);
        BOOST_CHECK(true);
        vec.Resize(3);
        TPZManVector<REAL,3> manvec2(vec);
        BOOST_CHECK(true);
        vec.Resize(5);
        TPZManVector<REAL,3> manvec3(vec);
        TPZManVector<REAL,3> manvec4(TPZManVector<REAL,2>(1,0.));
        TPZManVector<REAL,3> manvec5(TPZManVector<REAL,2>(2,0.));
        TPZManVector<REAL,3> manvec6(TPZManVector<REAL,2>(3,0.));
        TPZManVector<REAL,3> manvec7(TPZManVector<REAL,2>(4,0.));
        BOOST_CHECK(true);
    }
    //N1<N2
    {
        TPZManVector<REAL,5> vec(2,0.);
        TPZManVector<REAL,3> manvec1(vec);
        BOOST_CHECK(true);
        vec.Resize(3);
        TPZManVector<REAL,3> manvec2(vec);
        BOOST_CHECK(true);
        vec.Resize(5);
        TPZManVector<REAL,3> manvec3(vec);
        TPZManVector<REAL,3> manvec4(TPZManVector<REAL,5>(4,0.));
        TPZManVector<REAL,3> manvec5(TPZManVector<REAL,5>(5,0.));
        TPZManVector<REAL,3> manvec6(TPZManVector<REAL,5>(6,0.));
        BOOST_CHECK(true);
    }
}
#endif

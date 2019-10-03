/**
 * @file BlendUnitTest.cpp
 * @brief Define a Unit Test using Boost for blend elements over all kind of geometries available in the library
 *
 */
#include <iostream>

#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPattern3.h"

#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeoquad.h"

#include "pzgeotetrahedra.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"

//#include "pzgeoelrefless.h"
//
//#include "tpzquadraticline.h"
//
//#include "tpzquadratictrig.h"
//#include "tpzquadraticquad.h"
//
//#include "tpzquadratictetra.h"
//#include "tpzquadraticcube.h"
//#include "tpzquadraticprism.h"
//#include "tpzquadraticpyramid.h"
//
//#include "tpzgeoblend.h"
//
//#include "TPZVTKGeoMesh.h"
//
//#include "tpzgeoelrefpattern.h"
//
//#include "tpzarc3d.h"
//#include "TPZQuadSphere.h"
//#include "TPZTriangleSphere.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testgeom"));
#endif

#ifdef _AUTODIFF
#include "fad.h"
#endif

// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz blend_tests tests

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "boost/test/output_test_stream.hpp"

#endif

//#define NOISY_REFPATTERN //outputs x and grad comparisons
//#define NOISYVTK_REFPATTERN//prints all elements in .vtk format

std::string dirname = PZSOURCEDIR;

#ifdef USING_BOOST


BOOST_AUTO_TEST_SUITE(refpattern_tests)
    namespace refpatterntest{
        const REAL tol = 1e-8;
        template<class TGeo>
        void CompareWithOldRefPattern();
    }

    BOOST_AUTO_TEST_CASE(refpat_old_tests) {
        InitializePZLOG();
        refpatterntest::CompareWithOldRefPattern<pzgeom::TPZGeoLinear>();
        refpatterntest::CompareWithOldRefPattern<pzgeom::TPZGeoTriangle>();
        refpatterntest::CompareWithOldRefPattern<pzgeom::TPZGeoQuad>();
        refpatterntest::CompareWithOldRefPattern<pzgeom::TPZGeoTetrahedra>();
        refpatterntest::CompareWithOldRefPattern<pzgeom::TPZGeoCube>();
        refpatterntest::CompareWithOldRefPattern<pzgeom::TPZGeoPrism>();
        refpatterntest::CompareWithOldRefPattern<pzgeom::TPZGeoPyramid>();
    }

    template<class TGeo>
    void refpatterntest::CompareWithOldRefPattern(){
        MElementType elType = TGeo::Type();
        gRefDBase.InitializeUniformRefPattern(elType);
        auto oldRefPattern = gRefDBase.GetUniformRefPattern(elType);
        TPZAutoPointer<TPZRefPattern3> newRefPattern = new TPZRefPattern3(*(oldRefPattern.operator->()));
        oldRefPattern->PrintMore(std::cout);
        newRefPattern->PrintMore(std::cout);
}

BOOST_AUTO_TEST_SUITE_END()


#endif


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

#define NOISY_REFPATTERN //outputs x and grad comparisons
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
#ifdef NOISY_REFPATTERN
    std::cout<<"===================="<<std::endl;
    std::cout<<"TESTING ELEMENT TYPE "<<TGeo::TypeName()<<std::endl;
    std::cout<<"===================="<<std::endl;
#endif
        MElementType elType = TGeo::Type();
        gRefDBase.InitializeUniformRefPattern(elType);
        auto oldRefPattern = gRefDBase.GetUniformRefPattern(elType);
        TPZAutoPointer<TPZRefPattern3> newRefPattern = new TPZRefPattern3(oldRefPattern->RefPatternMesh());
//        TPZAutoPointer<TPZRefPattern3> newRefPattern = new TPZRefPattern3(*(oldRefPattern.operator->()));
//        oldRefPattern->PrintMore(std::cout);
//        newRefPattern->PrintMore(std::cout);

        const int nSubElOld = oldRefPattern->NSubElements();
        const int nSubElNew = newRefPattern->NSubElements();
        const bool isSameNSubEls = nSubElNew == nSubElOld;
        BOOST_CHECK(isSameNSubEls);
        for(int iSubEl = 0; (iSubEl < nSubElNew) && isSameNSubEls; iSubEl ++){
            TPZGeoEl *geoSubElNew = newRefPattern->Element(iSubEl + 1);
            TPZGeoEl *geoSubElOld = oldRefPattern->Element(iSubEl + 1);
            const int nSubElSidesOld = geoSubElOld->NSides();
            const int nSubElSidesNew = geoSubElNew->NSides();
            const bool isSameNSubElSides = nSubElSidesNew == nSubElSidesOld;
            BOOST_CHECK(isSameNSubEls);
            for(int iSubElSide = 0; isSameNSubElSides && (iSubElSide < nSubElSidesNew); iSubElSide++){
                const int fatherSideOld = oldRefPattern->FatherSide(iSubElSide,iSubEl);
                const int fatherSideNew = newRefPattern->FatherSide(iSubElSide,iSubEl);
                const bool isSameFatherSide = fatherSideNew == fatherSideOld;
                #ifdef NOISY_REFPATTERN
                std::cout<<"SUBEL: "<<iSubEl<<" SIDE: "<<iSubElSide<<std::endl;
                std::cout<<"\t----old----"<<std::endl;
                std::cout<<"\tfather side: "<<fatherSideOld<<std::endl;
                std::cout<<"\t----new----"<<std::endl;
                std::cout<<"\tfather side: "<<fatherSideNew<<std::endl;
                #endif
                BOOST_CHECK(isSameFatherSide);
                auto transformOld = oldRefPattern->Transform(iSubElSide,iSubEl);
                auto transformNew = newRefPattern->Transform(iSubElSide,iSubEl);
                bool checkTransform = !((bool)(transformOld.CompareTransform(transformNew)));
                if(!checkTransform){
                    int a = checkTransform;
                    a++;
                }
                BOOST_CHECK(checkTransform);
            }
        }

        for(int iFatherSide = 0; iFatherSide < TGeo::NSides; iFatherSide++){
            const int nSubElSideOld = oldRefPattern->NSideSubElements(iFatherSide);
            const int nSubElSideNew = newRefPattern->NSideSubGeoElSides(iFatherSide);
            const bool isSameNSubElSides = nSubElSideNew == nSubElSideOld;
            BOOST_CHECK(isSameNSubElSides);

            const int nSideNodesOld = oldRefPattern->NSideNodes(iFatherSide);
            const int nSideNodesNew = newRefPattern->NSideNodes(iFatherSide);
            const bool isSameNSideNodes = nSideNodesNew == nSideNodesOld;
            if(!isSameNSideNodes){
                int a = isSameNSideNodes;
                a++;
            }
            BOOST_CHECK(isSameNSideNodes);
            if(!isSameNSubElSides || iFatherSide < TGeo::NNodes) continue;
            for(int iSubEl = 0; iSubEl < nSubElSideNew; iSubEl ++){
                TPZGeoElSide subGeoElSide;
                newRefPattern->SideSubGeoElSide(iFatherSide, iSubEl, subGeoElSide);
                int subElIndexNew = subGeoElSide.Id();
                int subElSideNew = subGeoElSide.Side();
                int subElIndexOld = -1, subElSideOld = -1;
                oldRefPattern->SideSubElement(iFatherSide,iSubEl,subElIndexOld,subElSideOld);
                const bool isSameSubElIndex = subElIndexOld == subElIndexNew;
                const bool isSameSubElSide = subElSideOld == subElSideNew;
#ifdef NOISY_REFPATTERN
                std::cout<<"SIDE: "<<iFatherSide<<" SUBEL: "<<iSubEl<<std::endl;
                std::cout<<"\t----old----"<<std::endl;
                std::cout<<"\tindex: "<<subElIndexOld<<" side: "<<subElSideOld<<std::endl;
                std::cout<<"\t----new----"<<std::endl;
                std::cout<<"\tindex: "<<subElIndexNew<<" side: "<<subElSideNew<<std::endl;
#endif
                BOOST_CHECK(isSameSubElIndex);
                BOOST_CHECK(isSameSubElSide);
            }
        }

}

BOOST_AUTO_TEST_SUITE_END()

#ifdef NOISY_REFPATTERN
#undef NOISY_REFPATTERN
#endif
#ifdef NOISYVTK_REFPATTERN
#undef NOISYVTK_REFPATTERN
#endif
#endif


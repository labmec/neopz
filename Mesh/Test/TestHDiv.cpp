//
//  TestTopology.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/6/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzmanvector.h"
#include "pztrnsform.h"
#include "pzgengrid.h"
#include "tpzautopointer.h"
#include "pzpoisson3d.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "tpzpermutation.h"
#include "pzcompel.h"
#include "tpzintpoints.h"
#include "pztrnsform.h"
#include "pzintel.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testhdiv"));
#endif

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

TPZAutoPointer<TPZCompMesh> GenerateMesh(int type);
int CompareShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);
int CompareSideShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);

// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(mesh_tests)

BOOST_AUTO_TEST_CASE(sideshape_continuity)
{
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(0);
    TPZGeoMesh *gmesh = cmesh->Reference();
    TPZManVector<int,4> nodeids(4,0), nodesperm(4,0);
    TPZGeoEl *gel = gmesh->ElementVec()[4];
    for (int i = 0; i<4; i++) {
        nodeids[i] = gel->NodePtr(i)->Id();
    }
    TPZPermutation perm(4);
    perm++;
    while (!perm.IsFirst()) {
        perm.Permute(nodeids, nodesperm);
        for (int i = 0; i<4; i++) {
            gel->NodePtr(i)->SetNodeId(nodesperm[i]);
        }
        int nsides = gel->NSides();
        int is;
        int nwrong = 0;
        for (is=0; is<nsides; is++) {
            int dim = gel->SideDimension(is);
            if(dim != gel->Dimension()-1) continue;
            TPZStack<TPZCompElSide> connected;
            TPZGeoElSide gelside(gel,is);
            gelside.ConnectedCompElementList(connected, 1, 1);
            if (connected.NElements() != 1) {
                std::cout << "Number of elements connected along face side = " << connected.NElements() << std::endl;
                DebugStop();
            }
            TPZCompElSide celside = gelside.Reference();
            nwrong += CompareSideShapeFunctions(celside, connected[0]);
            nwrong += CompareShapeFunctions(celside, connected[0]);
            
        }
        if(nwrong)
        {
            std::cout << "Node ids " << nodesperm << " created incompatible shape functions\n";
        }
        BOOST_CHECK(nwrong == 0);
        // now compare the shape function value between neighbrouring elements
        // create a side integration rule
        // compute the transformation between the element and his neighbour
        // compute the vectors of shape function
        // compare
        perm++;
    }
}


BOOST_AUTO_TEST_SUITE_END()

TPZAutoPointer<TPZCompMesh> GenerateMesh(int type)
{
    TPZManVector<int,3> nx(2,3);
    TPZManVector<REAL,3> x0(3,0.),x1(3,3.);
    x1[2] = 0.;
    TPZGenGrid grid(nx,x0,x1);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    grid.Read(gmesh);
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    TPZAutoPointer<TPZMaterial> pois = new TPZMatPoisson3d(1, 2);
    cmesh->InsertMaterialObject(pois);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->SetDefaultOrder(3);
    cmesh->AutoBuild();
#ifdef LOG4CXX
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    return cmesh;
}

int CompareSideShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB)
{
    TPZGeoElSide gelsideA = celsideA.Reference();
    TPZGeoElSide gelsideB = celsideB.Reference();
    int sideA = gelsideA.Side();
    int sideB = gelsideB.Side();
    TPZCompEl *celA = celsideA.Element();
    TPZCompEl *celB = celsideB.Element();
    TPZInterpolatedElement *interA = dynamic_cast<TPZInterpolatedElement *>(celA);
    TPZInterpolatedElement *interB = dynamic_cast<TPZInterpolatedElement *>(celB);
    TPZGeoEl *gel = gelsideA.Element();
    int nshapeA = interA->NSideShapeF(sideA);
    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gelsideA.Side(), 4);
    int nwrong = 0;
    int npoints = intrule->NPoints();
    int ip;
    for (ip=0; ip<npoints; ip++) {
        TPZManVector<REAL,3> pointA(gelsideA.Dimension()),pointB(gelsideB.Dimension());
        REAL weight;
        intrule->Point(ip, pointA, weight);
        TPZTransform tr = gelsideA.NeighbourSideTransform(gelsideB);
        tr.Apply(pointA, pointB);
        TPZFNMatrix<200> phiA(nshapeA,1),dphiA(1,nshapeA),phiB(nshapeA,1),dphiB(1,nshapeA,1);
        interA->SideShapeFunction(sideA, pointA, phiA, dphiA);
        interB->SideShapeFunction(sideB, pointB, phiB, dphiB);
        int nshapeA = phiA.Rows();
        int nshapeB = phiB.Rows();
        BOOST_CHECK_EQUAL(nshapeA, nshapeB);
        int ish;
        for (ish=0; ish<nshapeA; ish++) {
            REAL Aval = phiA(ish,0);
            REAL Bval = phiB(ish,0);
            if(abs(Aval-Bval) > 1.e-6) 
            {
//                std::cout << "i " << ish << " " << Aval << " " << Bval << std::endl;
                nwrong++;   
            }
        }
//        if(nwrong)
//        {
//            std::cout << "\nNumber of different shape functions " << nwrong << std::endl;
//        }
//        BOOST_CHECK(nwrong == 0);
    }
    return nwrong;
}

int CompareShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB)
{
    TPZGeoElSide gelsideA = celsideA.Reference();
    TPZGeoElSide gelsideB = celsideB.Reference();
    int sideA = gelsideA.Side();
    int sideB = gelsideB.Side();
    TPZCompEl *celA = celsideA.Element();
    TPZCompEl *celB = celsideB.Element();
    TPZInterpolatedElement *interA = dynamic_cast<TPZInterpolatedElement *>(celA);
    TPZInterpolatedElement *interB = dynamic_cast<TPZInterpolatedElement *>(celB);
    TPZMaterialData dataA;
    TPZMaterialData dataB;
    interA->InitMaterialData(dataA);
    interB->InitMaterialData(dataB);
    TPZTransform tr = gelsideA.NeighbourSideTransform(gelsideB);
    TPZGeoEl *gelA = gelsideA.Element();
    TPZTransform trA = gelA->SideToSideTransform(gelsideA.Side(), gelA->NSides()-1);
    TPZGeoEl *gelB = gelsideB.Element();
    TPZTransform trB = gelB->SideToSideTransform(gelsideB.Side(), gelB->NSides()-1);
    
    int dimension = gelA->Dimension();
    
    int nshapeA = interA->NSideShapeF(sideA);
    int nshapeB = interB->NSideShapeF(sideB);
    int is;
    int firstShapeA = 0;
    int firstShapeB = 0;
    for (is=0; is<sideA; is++) {
        firstShapeA += interA->NSideShapeF(is);
    }
    for (is=0; is<sideB; is++) {
        firstShapeB += interB->NSideShapeF(is);
    }
    
    TPZIntPoints *intrule = gelA->CreateSideIntegrationRule(gelsideA.Side(), 4);
    int nwrong = 0;
    int npoints = intrule->NPoints();
    int ip;
    for (ip=0; ip<npoints; ip++) {
        TPZManVector<REAL,3> pointA(gelsideA.Dimension()),pointB(gelsideB.Dimension()), pointElA(gelA->Dimension()),pointElB(gelB->Dimension());
        REAL weight;
        intrule->Point(ip, pointA, weight);
        int sidedim = gelsideA.Dimension();
        TPZFNMatrix<9> jacobian(sidedim,sidedim),jacinv(sidedim,sidedim),axes(sidedim,3);
        REAL detjac;
        gelsideA.Jacobian(pointA, jacobian, jacinv, detjac, jacinv);
        TPZManVector<REAL,3> normal(3,0.);
        normal[0] = axes(0,1);
        normal[1] = -axes(0,0);
        tr.Apply(pointA, pointB);
        trA.Apply(pointA, pointElA);
        trB.Apply(pointB, pointElB);
        TPZFNMatrix<200> phiA(nshapeA,1),dphiA(dimension,nshapeA),phiB(nshapeB,1),dphiB(dimension,nshapeB);
        interA->Shape(pointElA, phiA, dphiA);
        interB->Shape(pointElB, phiB, dphiB);
        int nshapeA = phiA.Rows();
        int nshapeB = phiB.Rows();
        BOOST_CHECK_EQUAL(nshapeA, nshapeB);

        int i,j;
        for(i=firstShapeA,j=firstShapeB; i<firstShapeA+nshapeA; i++,j++)
        {
            int Avecind = dataA.fVecShapeIndex[i].first;
            int Ashapeind = dataA.fVecShapeIndex[i].second;
            int Bvecind = dataB.fVecShapeIndex[j].first;
            int Bshapeind = dataB.fVecShapeIndex[j].second;
            REAL vecnormalA = dataA.fNormalVec(0,Avecind)*normal[0]+dataA.fNormalVec(1,Avecind)*normal[1];
            REAL vecnormalB = dataB.fNormalVec(0,Bvecind)*normal[0]+dataB.fNormalVec(1,Bvecind)*normal[1];
            if(fabs(vecnormalA-vecnormalB) > 1.e-6)
            {
                nwrong++;
                LOGPZ_ERROR(logger, "normal vectors aren't equal")
            }
            if(fabs(phiA(Ashapeind,0)-phiB(Bshapeind,0)) > 1.-6)
            {
                nwrong ++;
                LOGPZ_ERROR(logger, "shape function values are different")
            }
        }
        //        if(nwrong)
        //        {
        //            std::cout << "\nNumber of different shape functions " << nwrong << std::endl;
        //        }
        //        BOOST_CHECK(nwrong == 0);
    }
    return nwrong;
}


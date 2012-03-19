//
//  TestCondensedElement.cpp
//  PZ
//
//  Created by Philippe Devloo on 1/2/12.
//  Copyright 2012 UNICAMP. All rights reserved.
//

//#include "TestCondensedElement.h"

#include "pzmanvector.h"
#include "pztrnsform.h"
#include "pzgengrid.h"
#include "tpzautopointer.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "tpzpermutation.h"
#include "pzcompel.h"
#include "tpzintpoints.h"
#include "pztrnsform.h"
#include "pzintel.h"
#include "pzstepsolver.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"

#include "pzcondensedcompel.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testmesh"));
#endif

#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN pz mesh tests

#include <boost/test/unit_test.hpp>

static TPZAutoPointer<TPZCompMesh> GenerateMesh(int type);


// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(meshcondensetests)

BOOST_AUTO_TEST_CASE(verifystiff)
{
    InitializePZLOG();
    std::cout << "Verifying creating and undoing condensed elements\n";
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(0);
    cmesh->ComputeNodElCon();
    std::ofstream out("TestCondensedMesh.txt");
    out << "Before condensing\n";
    cmesh->Print(out);
    TPZCreateApproximationSpace::CondenseLocalEquations(cmesh);
    cmesh->CleanUpUnconnectedNodes();
    out << "After condensing\n";
    cmesh->Print(out);
    TPZCreateApproximationSpace::UndoCondenseLocalEquations(cmesh);
    out << "After unwrapping the elements\n";
    cmesh->Print(out);
}

BOOST_AUTO_TEST_CASE(globalcompute)
{
    InitializePZLOG();
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(0);
    cmesh->ComputeNodElCon();
    TPZCreateApproximationSpace::CondenseLocalEquations(cmesh);
    cmesh->CleanUpUnconnectedNodes();
    TPZAnalysis an(cmesh);
    TPZSkylineStructMatrix skylstr(cmesh);
    an.SetStructuralMatrix(skylstr);
    TPZStepSolver<REAL> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
}



//void linpress(TPZVec<REAL> &x, TPZVec<REAL> &force)
//{
//    force[0] = x[0];
//}

BOOST_AUTO_TEST_SUITE_END()

static TPZAutoPointer<TPZCompMesh> GenerateMesh(int type)
{
    TPZManVector<int,3> nx(2,3);
    TPZManVector<REAL,3> x0(3,0.),x1(3,3.);
    x1[2] = 0.;
    TPZGenGrid grid(nx,x0,x1);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    grid.Read(gmesh);
    grid.SetBC(gmesh, 0, -1);
    grid.SetBC(gmesh, 1, -1);
    grid.SetBC(gmesh, 2, -1);
    grid.SetBC(gmesh, 3, -1);
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    TPZMatPoisson3d *matpois = new TPZMatPoisson3d(1, 2);
    TPZAutoPointer<TPZMaterial> pois(matpois);
    cmesh->InsertMaterialObject(pois);
    TPZFNMatrix<4> val1(1,1,0.),val2(1,1,0.);
    TPZBndCond *bnd = matpois->CreateBC(pois, -1, 0, val1, val2);
    TPZAutoPointer<TPZMaterial> matbnd(bnd);
    cmesh->InsertMaterialObject(matbnd);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(2);
    cmesh->SetDimModel(2);
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

// @TODO include tests with condensed elements and different type of solvers, substructed meshes etc. Verify the time gain of condensed elements


#endif

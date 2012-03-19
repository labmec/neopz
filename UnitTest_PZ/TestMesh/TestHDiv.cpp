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
#include "TPZParSkylineStructMatrix.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testhdiv"));
#endif

#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN pz meshHDiv tests

#include <boost/test/unit_test.hpp>

static TPZAutoPointer<TPZCompMesh> GenerateMesh(int type);
int CompareShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);
int CompareSideShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);

/// verify if the divergence of each vector function is included in the pressure space
static void CheckDRham(TPZInterpolatedElement *intel);


// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(mesh_tests)

BOOST_AUTO_TEST_CASE(sideshape_continuity)
{
    InitializePZLOG();
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

/// Check that the Div of the vector functions can be represented
BOOST_AUTO_TEST_CASE(drham_check)
{
    // generate a mesh
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(0);
    // for each computational element (not boundary) verify if the Div(vecspace) is included in the pressure space
    int nel = cmesh->NElements();
    int meshdim = cmesh->Dimension();
    int iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) continue;
        if(intel->Reference()->Dimension() != meshdim) continue;
        CheckDRham(intel);
    }
}

BOOST_AUTO_TEST_CASE(drham_permute_check)
{
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(0);
    TPZGeoMesh *gmesh = cmesh->Reference();
    int meshdim = cmesh->Dimension();

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
            TPZCompEl *cel = celside.Element();
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if(!intel) continue;
            if(intel->Reference()->Dimension() != meshdim) continue;
            CheckDRham(intel);            
        }
        if(nwrong)
        {
            std::cout << "Node ids " << nodesperm << " created drham incompatible shape functions\n";
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
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->SetDefaultOrder(3);
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
    
    int nSideshapeA = interA->NSideShapeF(sideA);
    int nSideshapeB = interB->NSideShapeF(sideB);
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
        int nshapeA = 0, nshapeB = 0;
        TPZFNMatrix<200> phiA(nshapeA,1),dphiA(dimension,nshapeA),phiB(nshapeB,1),dphiB(dimension,nshapeB);
        interA->Shape(pointElA, phiA, dphiA);
        interB->Shape(pointElB, phiB, dphiB);
        nshapeA = phiA.Rows();
        nshapeB = phiB.Rows();
        BOOST_CHECK_EQUAL(nshapeA, nshapeB);

        TPZManVector<REAL> shapesA(nSideshapeA), shapesB(nSideshapeB);
        int nwrongkeep(nwrong);
        int i,j;
        for(i=firstShapeA,j=firstShapeB; i<firstShapeA+nSideshapeA; i++,j++)
        {
            int Avecind = dataA.fVecShapeIndex[i].first;
            int Ashapeind = dataA.fVecShapeIndex[i].second;
            int Bvecind = dataB.fVecShapeIndex[j].first;
            int Bshapeind = dataB.fVecShapeIndex[j].second;
            REAL vecnormalA = dataA.fNormalVec(0,Avecind)*normal[0]+dataA.fNormalVec(1,Avecind)*normal[1];
            REAL vecnormalB = dataB.fNormalVec(0,Bvecind)*normal[0]+dataB.fNormalVec(1,Bvecind)*normal[1];
            shapesA[i-firstShapeA] = phiA(Ashapeind,0);
            shapesB[j-firstShapeB] = phiB(Bshapeind,0);
            if(fabs(vecnormalA-vecnormalB) > 1.e-6)
            {
                nwrong++;
                LOGPZ_ERROR(logger, "normal vectors aren't equal")
            }
            REAL diff = phiA(Ashapeind,0)-phiB(Bshapeind,0);
            REAL decision = fabs(diff)-1.e-6;
            if(decision > 0.)
            {
                nwrong ++;
                LOGPZ_ERROR(logger, "shape function values are different")
            }
        }
        if(nwrong != nwrongkeep)
        {
            std::cout << "shapeA " << shapesA << std::endl;
            std::cout << "shapeB " << shapesB << std::endl;
        }
        //        if(nwrong)
        //        {
        //            std::cout << "\nNumber of different shape functions " << nwrong << std::endl;
        //        }
        //        BOOST_CHECK(nwrong == 0);
    }
    delete intrule;
    return nwrong;
}

/// Generate the L2 matrix of the pressure space and the inner product of the divergence and the pressure shape functions
static void GenerateProjectionMatrix(TPZInterpolatedElement *intel, TPZAutoPointer<TPZMatrix<REAL> > L2, TPZFMatrix<REAL> &inner);

/// Given the multiplier coefficients of the pressure space, verify the correspondence of the divergence of the vector function and the L2 projection
static int VerifyProjection(TPZInterpolatedElement *intel, TPZFMatrix<REAL> &multiplier);

/// verify if the divergence of each vector function is included in the pressure space
static void CheckDRham(TPZInterpolatedElement *intel)
{
    TPZFMatrix<REAL> inner, multiplier;
    TPZAutoPointer<TPZMatrix<REAL> > L2 = new TPZFMatrix<REAL>;
    GenerateProjectionMatrix(intel, L2, inner);
    TPZStepSolver<REAL> step(L2);
    step.SetDirect(ELU);
    step.Solve(inner,multiplier);
    int nwrong = 0;
    nwrong = VerifyProjection(intel, multiplier);
    if(nwrong)
    {
        std::cout << "Number of points with wrong pressure projection " << nwrong << std::endl;
    }
    BOOST_CHECK(nwrong == 0);
    //return nwrong;

}

/// Generate the L2 matrix of the pressure space and the inner product of the divergence and the pressure shape functions
static void GenerateProjectionMatrix(TPZInterpolatedElement *intel, TPZAutoPointer<TPZMatrix<REAL> > L2, TPZFMatrix<REAL> &inner)
{
    TPZMaterialData dataA;
    intel->InitMaterialData(dataA);
    int dim = intel->Reference()->Dimension();
    const TPZIntPoints &intrule = intel->GetIntegrationRule();
    int np = intrule.NPoints();
    int npressure = dataA.numberdualfunctions;
    int nflux = dataA.fVecShapeIndex.NElements();
    L2->Redim(npressure,npressure);
    inner.Redim(npressure,nflux);
    int ip;
    for (ip=0; ip<np; ip++) {
        REAL weight;
        TPZManVector<REAL,3> pos(dim);
        intrule.Point(ip, pos, weight);
        intel->ComputeShape(pos, dataA.x, dataA.jacobian, dataA.axes, dataA.detjac, dataA.jacinv, dataA.phi, dataA.dphix);
        int firstpressure = dataA.phi.Rows()-npressure;
        int lastpressure = dataA.phi.Rows();
        int ish,jsh;
        for (ish=firstpressure; ish<lastpressure; ish++) {
            for (jsh=firstpressure; jsh<lastpressure; jsh++) {
                L2->s(ish-firstpressure,jsh-firstpressure) += dataA.phi(ish,0)*dataA.phi(jsh,0)*weight*fabs(dataA.detjac);
            }
            for (jsh=0; jsh<nflux; jsh++) {
                // compute the divergence of the shapefunction
                TPZManVector<REAL,3> vecinner(intel->Dimension(),0.);
                int vecindex = dataA.fVecShapeIndex[jsh].first;
                int phiindex = dataA.fVecShapeIndex[jsh].second;
                int j;
                int d;
                for (d=0; d<dim; d++) {
                    vecinner[d]=0;
                    for (j=0; j<3; j++) {
                        vecinner[d] += dataA.fNormalVec(j,vecindex)*dataA.axes(d,j);
                    }
                }
                REAL divphi = 0.;
                for (d=0; d<dim; d++) {
                    divphi += dataA.dphix(d,phiindex)*vecinner[d];                
                }
                inner(ish-firstpressure,jsh) += dataA.phi(ish,0)*divphi*weight*fabs(dataA.detjac);
            }
        }
    }
}

/// Given the multiplier coefficients of the pressure space, verify the correspondence of the divergence of the vector function and the L2 projection
static int VerifyProjection(TPZInterpolatedElement *intel, TPZFMatrix<REAL> &multiplier)
{
    TPZMaterialData dataA;
    intel->InitMaterialData(dataA);
    int dim = intel->Reference()->Dimension();
    const TPZIntPoints &intrule = intel->GetIntegrationRule();
    int np = intrule.NPoints();
    int npressure = dataA.numberdualfunctions;
    int nflux = dataA.fVecShapeIndex.NElements();
    TPZFNMatrix<30> pointpos(2,np);
    TPZFNMatrix<30> divergence(np,nflux);
    int ip;
    //std::cout << dataA.fVecShapeIndex << std::endl;
    int nwrong = 0;
    for (ip=0; ip<np; ip++) {
        REAL weight;
        TPZManVector<REAL,3> pos(dim);
        intrule.Point(ip, pos, weight);
        pointpos(0,ip) = pos[0];
        pointpos(1,ip) = pos[1];
        intel->ComputeShape(pos, dataA.x, dataA.jacobian, dataA.axes, dataA.detjac, dataA.jacinv, dataA.phi, dataA.dphix);
        int firstpressure = dataA.phi.Rows()-npressure;
        int lastpressure = dataA.phi.Rows();
        int ish,jsh;
        for (jsh=0; jsh<nflux; jsh++) {
            // compute the divergence of the shapefunction
            TPZManVector<REAL,3> vecinner(intel->Dimension(),0.);
            int vecindex = dataA.fVecShapeIndex[jsh].first;
            int phiindex = dataA.fVecShapeIndex[jsh].second;
            int j;
            int d;
            for (d=0; d<dim; d++) {
                vecinner[d]=0;
                for (j=0; j<3; j++) {
                    vecinner[d] += dataA.fNormalVec(j,vecindex)*dataA.axes(d,j);
                }
            }
            REAL divphi = 0.;
            for (d=0; d<dim; d++) {
                divphi += dataA.dphix(d,phiindex)*vecinner[d];                
            }
            divergence(ip,jsh) = divphi;
            REAL phival = 0;
            for (ish=firstpressure; ish<lastpressure; ish++) {
                phival += multiplier(ish-firstpressure,jsh)*dataA.phi(ish);
            }
            // the divergence of the vector function should be equal to the value of projected pressure space
            REAL diff = phival-divphi;
            if(fabs(diff) > 1.e-6) 
            {
                nwrong++;
                std::cout << "flux number " << jsh << " did not project\n";
            }
        }
    }
    /*
    int ifl;
    ofstream fluxes("fluxes.nb");
    for (ifl=0; ifl<nflux; ifl++) {
        fluxes << "flux" << ifl << " = {\n";
        for (ip=0; ip<np; ip++) {
            fluxes << "{ " << pointpos(0,ip) << " , " << pointpos(1,ip) << " , " << divergence(ip,ifl) << "} ";
            if(ip<np-1) fluxes << "," << std::endl;
        }
        fluxes << " };\n";
    }
     */
    return nwrong;
}

#endif
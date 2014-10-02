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
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "tpzpermutation.h"
#include "pzcompel.h"
#include "tpzintpoints.h"
#include "pztrnsform.h"
#include "pzintel.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZExtendGridDimension.h"

#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"

#include "pzlog.h"

#include "pzshapelinear.h"
#include "TPZRefPatternTools.h"
#include "pzshtmat.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"

#include "TPZVTKGeoMesh.h"

using namespace pzshape;


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testhdiv"));
#endif

#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz meshHDiv tests

#include <boost/test/unit_test.hpp>

static TPZAutoPointer<TPZCompMesh> GenerateMesh( TPZVec<TPZCompMesh *>  &meshvec,MElementType type, int nelem = 3, int fluxorder = 4, int ndiv = 0);
static int CompareShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);
static int CompareSideShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);

static TPZAutoPointer<TPZCompMesh> HDivMesh, PressureMesh;

template<class tshape>
static void CheckShapeOrder(int order);

template<class tshape>
void VectorDirections();
//static TPZCompMesh *HDivMesh, *PressureMesh;

static void ExactPressure(const TPZVec<REAL> &x, TPZVec<REAL> &force)
{
    force[0] = 5. ;//+ 3. * x[0] + 2. * x[1] + 4. * x[0] * x[1];
}
static void ExactNormalFluxTop(const TPZVec<REAL> &x, TPZVec<REAL> &force)
{
    force[0] = 0.;
}

static void ExactNormalFluxBottom(const TPZVec<REAL> &x, TPZVec<REAL> &force)
{
    force[0] = 0.;
}

/// verify if the divergence of each vector function is included in the pressure space
static void CheckDRham(TPZCompEl *cel);


// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(mesh_tests)


/// Check that the Div of the vector functions can be represented
BOOST_AUTO_TEST_CASE(bilinearsolution_check)
{
    InitializePZLOG();
    MElementType eltype = ETriangle; //ETriangle; //ECube;  //EQuadrilateral;
    int nelx = 1;
    int fluxorder = 1;
    int ndiv = 0; // para refinar a malha
     TPZVec<TPZCompMesh *>  meshvec(2);
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype,nelx,fluxorder, ndiv);

    
    {
        TPZMaterial *mat = cmesh->FindMaterial(-1);
        if(!mat) DebugStop();
        TPZDummyFunction<STATE> *dumforce = new TPZDummyFunction<STATE>(ExactPressure);
        TPZAutoPointer<TPZFunction<STATE> > autofunc (dumforce);
        mat->SetForcingFunction(autofunc);
    }
    {
        TPZMaterial *mat = cmesh->FindMaterial(-2);
        if(!mat) DebugStop();
        TPZDummyFunction<STATE> *dumforce = new TPZDummyFunction<STATE>(ExactNormalFluxTop);
        TPZAutoPointer<TPZFunction<STATE> > autofunc (dumforce);
        mat->SetForcingFunction(autofunc);
    }
    
    {
        TPZMaterial *mat = cmesh->FindMaterial(-3);
        if(!mat) DebugStop();
        TPZDummyFunction<STATE> *dumforce = new TPZDummyFunction<STATE>(ExactNormalFluxBottom);
        TPZAutoPointer<TPZFunction<STATE> > autofunc (dumforce);
        mat->SetForcingFunction(autofunc);
    }

    
//    {
//        //  Print Geometrical Base Mesh
//        std::ofstream Dummyfile("GeometricMesh.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(cmesh,Dummyfile, true);
//    }
    
   // cmesh->SaddlePermute();
    TPZAnalysis an(cmesh);
    // para resolver o sistema
    // escolhe entre isso 
//    TPZFStructMatrix str(cmesh);
//    an.SetStructuralMatrix(str);
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELU);
    /// ou isso
    TPZSkylineStructMatrix str(cmesh);
    an.SetStructuralMatrix(str);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    
    
    an.SetSolver(step);
    an.Run();

    
    an.Solution().Print("Solucao");
    
    std::string plotfile("GSaida.vtk");
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmesh.operator->());
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
	const int dim = cmesh->Dimension();
	int div =1;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);

    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Reference()->Print(sout);
        cmesh->Print(sout);
        an.Rhs().Print("Right Hand Side",sout);
        an.Solution().Print("Solution",sout);
        LOGPZ_DEBUG(logger, sout.str())
        std::ofstream test("onde_esta_este_arquivo.txt");
    }
#endif
    TPZCompEl *cel = cmesh->ElementVec()[0];
    TPZGeoEl *gel = cel->Reference();
    int ns = gel->NSides();
    TPZIntPoints *rule = gel->CreateSideIntegrationRule(ns-1, 4);//3
    int np = rule->NPoints();
    for(int ip=0; ip<np; ip++)
    {
        //TPZManVector<REAL,3> xi(2), xco(3), sol(1), exactsol(1);
        TPZManVector<REAL,3> xi(3), xco(3), sol(1), exactsol(1);
        REAL weight;
        rule->Point(ip, xi, weight);
        gel->X(xi, xco);
        cel->Solution(xi, 2, sol);
        ExactPressure(xco,exactsol);
//        autofunc->Execute(xco, exactsol);
        BOOST_CHECK(fabs(sol[0]-exactsol[0]) < 1.e-6);
    }
    
}


BOOST_AUTO_TEST_CASE(sideshape_continuity)
{
    InitializePZLOG();
    for (int i=0; i<8; i++) {
        TPZManVector<int,9> permgather(9);
        pztopology::TPZQuadrilateral::GetSideHDivPermutation(i,permgather);
        std::cout << "transform id " << i << " permgather " << permgather << std::endl;
        
    }
    MElementType eltype = ECube;
    TPZVec<TPZCompMesh *>  meshvec(2);
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype);
    TPZGeoMesh *gmesh = cmesh->Reference();
    TPZGeoEl *gel = gmesh->ElementVec()[0];
    const int ncorner = gel->NCornerNodes();
    TPZManVector<int,8> nodeids(ncorner,0), nodesperm(ncorner,0);
    for (int i = 0; i<ncorner; i++) {
        nodeids[i] = gel->NodePtr(i)->Id();
    }
    TPZPermutation perm(ncorner);
    perm++;
    long permcount = 1;
    while (!perm.IsFirst()) {
        perm.Permute(nodeids, nodesperm);
        for (int i = 0; i<ncorner; i++) {
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
            gelside.ConnectedCompElementList(connected, 0, 1);
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
        permcount++;
        if (!(permcount%1000)) {
            std::cout << "permcount = " << permcount << std::endl;
        }
    }
}
                
/// Check that the Div of the vector functions can be represented
BOOST_AUTO_TEST_CASE(drham_check)

{
		InitializePZLOG();
    // generate a mesh
    MElementType eltype = ECube;
    TPZVec<TPZCompMesh *>  meshvec(2);
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype);
    // for each computational element (not boundary) verify if the Div(vecspace) is included in the pressure space
    int nel = cmesh->NElements();
    int meshdim = cmesh->Dimension();
    int iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel)
        {
            DebugStop();
        }
        if(intel->Reference()->Dimension() != meshdim) continue;
        CheckDRham(intel);
    }
}

BOOST_AUTO_TEST_CASE(drham_permute_check)
{
    MElementType eltype = ECube;
    TPZVec<TPZCompMesh *>  meshvec(2);
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype);
    TPZGeoMesh *gmesh = cmesh->Reference();
    int meshdim = cmesh->Dimension();

    const long nel = gmesh->NElements();
    TPZGeoEl *gel = gmesh->ElementVec()[nel/2];
    const int ncorner = gel->NCornerNodes();
    TPZManVector<int,8> nodeids(ncorner,0), nodesperm(ncorner,0);
    for (int i = 0; i<ncorner; i++) {
        nodeids[i] = gel->NodePtr(i)->Id();
    }
    TPZPermutation perm(ncorner);
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
            gelside.ConnectedCompElementList(connected, 0, 1);
            if (connected.NElements() != 1) {
                std::cout << "Number of elements connected along face side = " << connected.NElements() << std::endl;
                DebugStop();
            }
            TPZCompElSide celside = gelside.Reference();
            TPZCompEl *cel = celside.Element();
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if(!intel)
            {
                DebugStop();
            }
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


BOOST_AUTO_TEST_CASE(shape_order)
{
    CheckShapeOrder<pzshape::TPZShapeQuad>(6);
    CheckShapeOrder<pzshape::TPZShapeTriang>(6);
    CheckShapeOrder<pzshape::TPZShapeCube>(6);
}

BOOST_AUTO_TEST_CASE(vector_direction)
{
    VectorDirections<pzshape::TPZShapeCube>();
    VectorDirections<pzshape::TPZShapeTriang>();
    VectorDirections<pzshape::TPZShapeQuad>();
}

BOOST_AUTO_TEST_SUITE_END()

static TPZAutoPointer<TPZCompMesh> GenerateMesh( TPZVec<TPZCompMesh *>  &meshvec, MElementType eltype, int nelem, int fluxorder, int ndiv)
{
    int dimmodel = 2;
    TPZManVector<int,3> nx(2,nelem);
    TPZManVector<REAL,3> x0(3,-1.),x1(3,1.);
    x1[2] = -1.;
    TPZGenGrid grid(nx,x0,x1);
    if (eltype == ETriangle || eltype == EPrisma) {
        grid.SetElementType(ETriangle);
    }
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    grid.Read(gmesh);
    grid.SetBC(gmesh, 4, -1);
    grid.SetBC(gmesh, 5, -1);
    grid.SetBC(gmesh, 6, -1);
    grid.SetBC(gmesh, 7, -1);
    
   
        for(int D = 0; D < ndiv; D++)
        {
            int nels = gmesh->NElements();
            for(int elem = 0; elem < nels; elem++)
            {
                TPZVec< TPZGeoEl * > filhos;
                TPZGeoEl * gel = gmesh->ElementVec()[elem];
                gel->Divide(filhos);
            }
        }

    


    
    if (eltype == EPrisma || eltype == ECube) {
        REAL thickness = 2.;
        TPZExtendGridDimension extend(gmesh,thickness);
        int numlayers = nelem;
        int bctop = -2;
        int bcbottom = -3 ;//normal negativa
        gmesh = extend.ExtendedMesh(numlayers,bcbottom,bctop);
        dimmodel = 3;
    }
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        gmesh->Print(sout);
        //mphysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZMatPoisson3d *matpoisP = new TPZMatPoisson3d(1, dimmodel);
    TPZMaterial *poisP(matpoisP);
    
    PressureMesh = new TPZCompMesh(gmesh);
    
    TPZFNMatrix<4,STATE> val1(1,1,0.),val2(1,1,0.);
    PressureMesh->InsertMaterialObject(poisP);
    if (eltype == ETriangle) {
        PressureMesh->SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        PressureMesh->SetAllCreateFunctionsContinuous();
        PressureMesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
//    TPZBndCond *bndP = matpoisP->CreateBC(poisP, -1, 0, val1, val2);
//    TPZMaterial *matbndP(bndP);
//    PressureMesh->InsertMaterialObject(matbndP);
//    bndP = matpoisP->CreateBC(poisP, -2, 1, val1, val2);
//    PressureMesh->InsertMaterialObject(bndP);
//    
//    TPZMaterial *matbndP2(bndP);
//    PressureMesh->InsertMaterialObject(matbndP2);
//    bndP = matpoisP->CreateBC(poisP, -3, 1, val1, val2);
//    PressureMesh->InsertMaterialObject(bndP);
    
    if (eltype != ETriangle) {
        PressureMesh->SetDefaultOrder(fluxorder);
    }
    else
    {
        PressureMesh->SetDefaultOrder(fluxorder-1);
    }
    PressureMesh->SetDimModel(dimmodel);
    std::set<int> matids;
    matids.insert(1);
    PressureMesh->AutoBuild(matids);
    PressureMesh->LoadReferences();
//    PressureMesh->ApproxSpace().CreateDisconnectedElements(false);
//    PressureMesh->AutoBuild();
    
    int ncon = PressureMesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = PressureMesh->ConnectVec()[i];
        //newnod.SetPressure(true);
        newnod.SetLagrangeMultiplier(1);
    }
    
    
    TPZMatPoisson3d *matpoisH = new TPZMatPoisson3d(1, dimmodel);
    TPZMaterial *poisH(matpoisH);
    HDivMesh = new TPZCompMesh(gmesh);
    HDivMesh->SetDimModel(dimmodel);
    TPZBndCond *bndh = matpoisH->CreateBC(poisH, -1, 0, val1, val2);
    HDivMesh->InsertMaterialObject(poisH);
    HDivMesh->SetAllCreateFunctionsHDiv();
    
    TPZMaterial *matbndh(bndh);
    HDivMesh->InsertMaterialObject(matbndh);
    
    bndh = matpoisH->CreateBC(poisH, -2, 1, val1, val2);
    HDivMesh->InsertMaterialObject(bndh);
    
    bndh = matpoisH->CreateBC(poisH, -3, 1, val1, val2);
    HDivMesh->InsertMaterialObject(bndh);
    
    HDivMesh->SetDefaultOrder(fluxorder);
    HDivMesh->SetDimModel(dimmodel);
    HDivMesh->AutoBuild();
    
    
    // malha multifisica
   // TPZVec<TPZCompMesh *>  meshvec(2);
    // static TPZAutoPointer<TPZCompMesh> HDivMesh, PressureMesh;
    meshvec[0] = HDivMesh.operator->();
    meshvec[1] = PressureMesh.operator->();
    
    gmesh->ResetReference();
    
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    TPZMixedPoisson *matpoisM = new TPZMixedPoisson(1, dimmodel);
    //TPZMatPoissonD3 *matpoisM = new TPZMatPoissonD3(1, dimmodel);
    TPZMaterial *poisM(matpoisM);
    mphysics->SetDimModel(dimmodel);
    mphysics->InsertMaterialObject(poisM);
    
    
    //Condicao de contorno Dirichlet
    TPZMaterial *bndm = matpoisM->CreateBC(poisM, -1, 0, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > BND;
//    BND = new TPZDummyFunction<STATE>(Force);//mudar force para SolExata
//    bndm->SetForcingFunction(BND);
    mphysics->InsertMaterialObject(bndm);
    
    
    
    //Condicao de contorno NEUMANN acima
    TPZMaterial *bnnm = matpoisM->CreateBC(poisM, -2, 1, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > BNN1;
//    BND = new TPZDummyFunction<STATE>(ExactNormalFlux);
    mphysics->InsertMaterialObject(bnnm);
    
    //Condicao de contorno NEUMANN abaixo
    TPZMaterial *bnnm2 = matpoisM->CreateBC(poisM, -3, 1, val1, val2);
//    TPZAutoPointer<TPZFunction<STATE> > BNN2;
//    BND = new TPZDummyFunction<STATE>(ExactNormalFlux);
    mphysics->InsertMaterialObject(bnnm2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, mphysics);
    
    {
        long nel = mphysics->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = mphysics->Element(el);
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mfcel) {
                DebugStop();
            }
            if (mfcel->NMeshes() != 2) {
                mfcel->AddElement(0, 1);
            }
        }
    }
    
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();

//    ofstream arq4("mphysics.txt");
//    mphysics->Print(arq4);
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
//        gmesh->Print(sout);
        mphysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    

    
    return mphysics;

    
    
    
    
    
  /**
    TPZManVector<int,3> nx(2,nelem);
    TPZManVector<REAL,3> x0(3,0.),x1(3,3.);
    x1[2] = 0.;
    TPZGenGrid grid(nx,x0,x1);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    grid.Read(gmesh);
    grid.SetBC(gmesh, 4, -1);
    grid.SetBC(gmesh, 5, -1);
    grid.SetBC(gmesh, 6, -1);
    grid.SetBC(gmesh, 7, -1);
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    TPZMatPoisson3d *matpois = new TPZMatPoisson3d(1, 2);
    TPZMaterial *pois(matpois);
    cmesh->InsertMaterialObject(pois);
    TPZFNMatrix<4,STATE> val1(1,1,0.),val2(1,1,0.);
    TPZBndCond *bnd = matpois->CreateBC(pois, -1, 0, val1, val2);
    TPZMaterial *matbnd(bnd);
    cmesh->InsertMaterialObject(matbnd);
    bnd = matpois->CreateBC(pois, -2, 1, val1, val2);
    cmesh->InsertMaterialObject(bnd);
    cmesh->SetAllCreateFunctionsHDivPressure();
    cmesh->SetDefaultOrder(fluxorder);
    cmesh->SetDimModel(2);
    cmesh->AutoBuild();
#ifdef LOG4CXX
    {
        std::stringstream sout;
        gmesh->Print(sout);
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    return cmesh;
   */
}

int CompareSideShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB)
{
    TPZGeoElSide gelsideA = celsideA.Reference();
    TPZGeoElSide gelsideB = celsideB.Reference();
    int sideA = gelsideA.Side();
    int sideB = gelsideB.Side();
    const int sidedimension = gelsideA.Dimension();
    TPZCompEl *celA = celsideA.Element();
    TPZCompEl *celB = celsideB.Element();
    TPZMultiphysicsElement *MFcelA = dynamic_cast<TPZMultiphysicsElement *>(celA);
    TPZMultiphysicsElement *MFcelB = dynamic_cast<TPZMultiphysicsElement *>(celB);
    TPZInterpolatedElement *interA = dynamic_cast<TPZInterpolatedElement *>(MFcelA->Element(0));
    TPZInterpolatedElement *interB = dynamic_cast<TPZInterpolatedElement *>(MFcelB->Element(0));
    if (!interA || !interB) {
        DebugStop();
    }
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
        TPZFNMatrix<200> phiA(nshapeA,1),dphiA(sidedimension,nshapeA),phiB(nshapeA,1),dphiB(sidedimension,nshapeA);
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
                std::cout << "i " << ish << " " << Aval << " " << Bval << std::endl;
                nwrong++;   
            }
        }
        if(nwrong)
        {
            std::cout << "\nNumber of different shape functions " << nwrong << std::endl;
        }
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
    TPZCompEl *celB = celsideB.Element();    TPZMultiphysicsElement *MFcelA = dynamic_cast<TPZMultiphysicsElement *>(celA);
    TPZMultiphysicsElement *MFcelB = dynamic_cast<TPZMultiphysicsElement *>(celB);
    TPZInterpolatedElement *interA = dynamic_cast<TPZInterpolatedElement *>(MFcelA->Element(0));
    TPZInterpolatedElement *interB = dynamic_cast<TPZInterpolatedElement *>(MFcelB->Element(0));

    TPZMaterialData dataA;
    TPZMaterialData dataB;
    interA->InitMaterialData(dataA);
    interB->InitMaterialData(dataB);
    TPZTransform tr = gelsideA.NeighbourSideTransform(gelsideB);
    TPZGeoEl *gelA = gelsideA.Element();
    TPZTransform trA = gelA->SideToSideTransform(gelsideA.Side(), gelA->NSides()-1);
    TPZGeoEl *gelB = gelsideB.Element();
    TPZTransform trB = gelB->SideToSideTransform(gelsideB.Side(), gelB->NSides()-1);
    
    int dimensionA = gelA->Dimension();
    int dimensionB = gelB->Dimension();
    
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
        TPZManVector<REAL,3> normal(3,0.), xA(3),xB(3);
        normal[0] = axes(0,1);
        normal[1] = -axes(0,0);
        tr.Apply(pointA, pointB);
        trA.Apply(pointA, pointElA);
        trB.Apply(pointB, pointElB);
        gelsideA.Element()->X(pointElA, xA);
        gelsideB.Element()->X(pointElB, xB);
        for (int i=0; i<3; i++) {
            BOOST_CHECK_CLOSE(xA[i], xB[i], 1.e-6);
        }
        int nshapeA = 0, nshapeB = 0;
        interA->ComputeRequiredData(dataA, pointElA);
        interB->ComputeRequiredData(dataB, pointElB);
        nshapeA = dataA.phi.Rows();
        nshapeB = dataB.phi.Rows();
        BOOST_CHECK_EQUAL(nSideshapeA, nSideshapeB);

        TPZManVector<REAL> shapesA(nSideshapeA), shapesB(nSideshapeB);
        int nwrongkeep(nwrong);
        int i,j;
        for(i=firstShapeA,j=firstShapeB; i<firstShapeA+nSideshapeA; i++,j++)
        {
            int Ashapeind = i;
            int Bshapeind = j;
            int Avecind = -1;
            int Bvecind = -1;
            // if A or B are boundary elements, their shapefunctions come in the right order
            if (dimensionA != sidedim) {
                Ashapeind = dataA.fVecShapeIndex[i].second;
                Avecind = dataA.fVecShapeIndex[i].first;
            }
            if (dimensionB != sidedim) {
                Bshapeind = dataB.fVecShapeIndex[j].second;
                Bvecind = dataB.fVecShapeIndex[j].first;
            }
            if (dimensionA != sidedim && dimensionB != sidedim) {
                // vefify that the normal component of the normal vector corresponds
                Avecind = dataA.fVecShapeIndex[i].first;
                Bvecind = dataB.fVecShapeIndex[j].first;
                REAL vecnormalA = dataA.fNormalVec(0,Avecind)*normal[0]+dataA.fNormalVec(1,Avecind)*normal[1];
                REAL vecnormalB = dataB.fNormalVec(0,Bvecind)*normal[0]+dataB.fNormalVec(1,Bvecind)*normal[1];
                if(fabs(vecnormalA-vecnormalB) > 1.e-6)
                {
                    nwrong++;
                    LOGPZ_ERROR(logger, "normal vectors aren't equal")
                }

            }
            shapesA[i-firstShapeA] = dataA.phi(Ashapeind,0);
            shapesB[j-firstShapeB] = dataB.phi(Bshapeind,0);
            REAL valA = dataA.phi(Ashapeind,0);
            REAL valB = dataB.phi(Bshapeind,0);
            REAL diff = valA-valB;
            REAL decision = fabs(diff)-1.e-6;
            if(decision > 0.)
            {
                nwrong ++;
                std::cout << "valA = " << valA << " valB = " << valB << " Avecind " << Avecind << " Bvecind " << Bvecind <<
                " Ashapeind " << Ashapeind << " Bshapeind " << Bshapeind <<
                " sideA " << sideA << " sideB " << sideB << std::endl;
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
static void GenerateProjectionMatrix(TPZCompEl *cel, TPZAutoPointer<TPZMatrix<STATE> > L2, TPZFMatrix<STATE> &inner);

/// Given the multiplier coefficients of the pressure space, verify the correspondence of the divergence of the vector function and the L2 projection
static int VerifyProjection(TPZCompEl *cel, TPZFMatrix<STATE> &multiplier);

/// verify if the divergence of each vector function is included in the pressure space
static void CheckDRham(TPZCompEl *cel)
{
    TPZFMatrix<STATE> inner, multiplier;
    TPZAutoPointer<TPZMatrix<STATE> > L2 = new TPZFMatrix<STATE>;
    GenerateProjectionMatrix(cel, L2, inner);
    TPZStepSolver<STATE> step(L2);
    step.SetDirect(ELU);
    step.Solve(inner,multiplier);
    int nwrong = 0;
    nwrong = VerifyProjection(cel, multiplier);
    if(nwrong)
    {
        std::cout << "Number of points with wrong pressure projection " << nwrong << std::endl;
    }
    BOOST_CHECK(nwrong == 0);
    //return nwrong;

}

/// Generate the L2 matrix of the pressure space and the inner product of the divergence and the pressure shape functions
static void GenerateProjectionMatrix(TPZCompEl *cel, TPZAutoPointer<TPZMatrix<STATE> > L2, TPZFMatrix<STATE> &inner)
{
    TPZMaterialData dataA,dataB;
    TPZMultiphysicsElement *celMF = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if (!celMF) {
        DebugStop();
    }
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(celMF->Element(0));
    TPZInterpolatedElement *intelP = dynamic_cast<TPZInterpolatedElement *>(celMF->Element(1));
    intel->InitMaterialData(dataA);
    intelP->InitMaterialData(dataB);
    int dim = intel->Reference()->Dimension();
    const TPZIntPoints &intrule = intel->GetIntegrationRule();
    int np = intrule.NPoints();
    int npressure = dataB.phi.Rows();
    int nflux = dataA.fVecShapeIndex.NElements();
    L2->Redim(npressure,npressure);
    inner.Redim(npressure,nflux);
    int ip;
    for (ip=0; ip<np; ip++) {
        REAL weight;
        TPZManVector<REAL,3> pos(dim);
        intrule.Point(ip, pos, weight);
//        intel->ComputeShape(pos, dataA.x, dataA.jacobian, dataA.axes, dataA.detjac, dataA.jacinv, dataA.phi, dataA.dphix);
        intel->ComputeShape(pos,dataA);
        intelP->ComputeShape(pos, dataB);
        int ish,jsh;
        for (ish=0; ish<npressure; ish++) {
            for (jsh=0; jsh<npressure; jsh++) {
                L2->s(ish,jsh) += dataB.phi(ish,0)*dataB.phi(jsh,0)*weight*fabs(dataB.detjac);
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
                inner(ish,jsh) += dataB.phi(ish,0)*divphi*weight*fabs(dataA.detjac);
            }
        }
    }
}

/// Given the multiplier coefficients of the pressure space, verify the correspondence of the divergence of the vector function and the L2 projection
static int VerifyProjection(TPZCompEl *cel, TPZFMatrix<STATE> &multiplier)
{
    TPZMaterialData dataA,dataB;
    TPZMultiphysicsElement *celMF = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if (!celMF) {
        DebugStop();
    }
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(celMF->Element(0));
    TPZInterpolatedElement *intelP = dynamic_cast<TPZInterpolatedElement *>(celMF->Element(1));
    intel->InitMaterialData(dataA);
    intelP->InitMaterialData(dataB);
    int dim = intel->Reference()->Dimension();
    const TPZIntPoints &intrule = intel->GetIntegrationRule();
    int np = intrule.NPoints();
    int npressure = dataB.phi.Rows();
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
//        intel->ComputeShape(pos, dataA.x, dataA.jacobian, dataA.axes, dataA.detjac, dataA.jacinv, dataA.phi, dataA.dphix);
        intel->ComputeShape(pos,dataA);
        intelP->ComputeShape(pos, dataB);
#ifdef LOG4CXX
				{
						std::stringstream sout;
						sout << "Phi's " << dataA.phi<< " dphix's "<< dataA.dphix<<std::endl;
						
						LOGPZ_DEBUG(logger,sout.str())
				}
#endif	
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
//#ifdef LOG4CXX
//						{
//								std::stringstream sout;
//								sout << "Div " << divphi<< std::endl;
//								
//								LOGPZ_DEBUG(logger,sout.str())
//						}
//#endif
            divergence(ip,jsh) = divphi;
            REAL phival = 0;

            for (ish=0; ish<npressure; ish++) {

                phival += multiplier(ish,jsh)*dataB.phi(ish);
            }
            // the divergence of the vector function should be equal to the value of projected pressure space
            REAL diff = phival-divphi;
#ifdef LOG4CXX
						{
						std::stringstream sout;
						sout << "phi: " << phival<<" dphi: "<< divphi <<"\n";
						sout << "flux number " << jsh << " diff: "<<diff<< "\n";
								LOGPZ_DEBUG(logger,sout.str())
						}
#endif
            if(fabs(diff) > 1.e-6) 
            {
                nwrong++;
                std::cout << "flux number " << jsh << " did not project: diff: "<<diff<<"\n";
            }
        }
    }
    
/*
    int ifl;
    std::ofstream fluxes("fluxes.nb");
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


template<class tshape>
void CheckShapeOrder(int order)
{
    TPZManVector<int,tshape::NSides> orders(tshape::NSides-tshape::NCornerNodes,order);
    TPZManVector<long,tshape::NSides> origids(tshape::NCornerNodes,0);
    const int nshape = tshape::NShapeF(orders);
    const int nsides = tshape::NSides;
    const int ncorner = tshape::NCornerNodes;
    const int numlegendre = (order+5) < 10 ? 10 : (order+5);
    const int dimension = tshape::Dimension;
    int numpermwrong = 0;
    
    // set up the original ids
    for (int ic = 0; ic<ncorner; ic++) {
        origids[ic] = ic;
    }
    
    // set up the permutations of the type
    TPZManVector<TPZVec<int> > permutation;
    MElementType thistype = tshape::Type();
    TPZRefPatternTools::GetElTypePermutations(thistype, permutation);
    const int nperm = permutation.size();
    
    // set up the integration rule
    TPZInt1d localrule(20);
    
    for (int iperm = 0; iperm<nperm; iperm++)
    {
        int numwrong = 0;
        TPZManVector<long,tshape::NCornerNodes> ids(tshape::NCornerNodes,0);
        for (int id = 0; id<ncorner; id++) {
            ids[permutation[iperm][id]] = id;
        }
        TPZGenMatrix<int> shapeorders(nshape,3);
        TPZGenMatrix<int> estimatedshapeorders(nshape,3);
        for (int ish=0; ish<nshape; ish++) {
            for (int d=0; d<3; d++) {
                estimatedshapeorders(ish,d) = 0;
            }
        }
        TPZVec<long> sides;
        tshape::ShapeOrder(ids, orders, shapeorders);

        int shapecounter = 0;
        for (int is=0; is<nsides; is++) {
            const int sidedim = tshape::SideDimension(is);
            if (sidedim < 1) {
                for (int dim=0; dim<1; dim++) {
                    estimatedshapeorders(shapecounter,dim) = 1;
                }
                shapecounter++;
                continue;
            }
            MElementType sidetype = tshape::Type(is);

            int nsideshape = tshape::NConnectShapeF(is, order);
            TPZStack<int> lowerdimensionsides;
            tshape::LowerDimensionSides(is, lowerdimensionsides);
            int firstshape = 0;
            for (int i=0; i<lowerdimensionsides.size(); i++) {
                firstshape += tshape::NConnectShapeF(lowerdimensionsides[i], order);
            }
            
            /// Compute the ids associated with the nodes of the side
            TPZManVector<long,8> locids(tshape::NSideNodes(is));
            for (int i=0; i<locids.size(); i++) {
                locids[i] = ids[tshape::SideNodeLocId(is, i)];
            }
            
            /// Estimate the shape orders in each direction
            for (int dim = 0; dim < sidedim; dim++) {
                TPZManVector<REAL,3> sidecenterel(dimension),sidecenter(sidedim),point(sidedim);
                TPZFNMatrix<9> jac(2,1);
                tshape::CenterPoint(is, sidecenterel);
                tshape::MapToSide(is, sidecenterel, sidecenter, jac);
                for (int i=0; i<sidedim; i++) {
                    sidecenter[i] += M_PI/20.;
                }
                TPZFNMatrix<500,REAL> integrationvals(numlegendre,nsideshape,0.);
                TPZManVector<REAL,1> pos(1);
                REAL weight;
                for (int ip=0; ip<localrule.NPoints(); ip++) {
                    localrule.Point(ip, pos, weight);
                    point = sidecenter;
                    point[dim]=pos[0];
//                    point[dim] = (pos[0]+1.)/2.;
                    if (sidetype == ETriangle) {
                        REAL a,b;
                        a = (point[0]-sidecenter[0])*cos(M_PI/20.)-(point[1]-sidecenter[1])*sin(M_PI/20.);
                        b = (point[0]-sidecenter[0])*sin(M_PI/20.)+(point[1]-sidecenter[1])*cos(M_PI/20.);
                        point[0] = a;
                        point[1] = b;
                    }
                    TPZFNMatrix<200,REAL> philegendre(numlegendre,1,0.);
                    TPZFNMatrix<200,REAL> dphi(sidedim,nsideshape+firstshape), dphilegendre(1,numlegendre,0.), phi(nsideshape+firstshape,1,0.);
                    tshape::SideShape(is, point, locids, orders, phi, dphi);
                    TPZShapeLinear::Legendre(pos[0],numlegendre,philegendre,dphilegendre);
                    for (int il=0; il<numlegendre; il++) {
                        for (int ishape = 0; ishape < nsideshape; ishape++) {
                            integrationvals(il,ishape) += philegendre(il)*phi(firstshape+ishape)*weight;
                        }
                    }
                }
                
//                integrationvals.Print("Integration values",std::cout,EMathematicaInput);
                
                
                for (int ishape = 0; ishape < nsideshape; ishape++) {
                    int ilegendre = numlegendre-1;
                    while (fabs(integrationvals(ilegendre,ishape)) < 1.e-6 && ilegendre > 0) {
                        ilegendre--;
                    }
                    estimatedshapeorders(shapecounter+ishape,dim) = ilegendre;
                }
                
            }
            shapecounter += nsideshape;
        }
        for (int ishape = 0; ishape < nshape; ishape++) {
            for( int idim = 0; idim<3; idim++)
            {
                if (shapeorders(ishape,idim) != estimatedshapeorders(ishape,idim)) {
                    numwrong++;
                }
            }
        }
        if (numwrong > 0) {
            numpermwrong++;
        }
        if (numwrong) {
            std::cout << "Permutation number " << iperm << " resulted in wrong shape order computation\n";
            shapeorders.Print("ShapeOrders as indicated by the method");
            estimatedshapeorders.Print("ShapeOrders estimated by fitting");
        }
    }
    BOOST_CHECK(numpermwrong == 0);

}

template<class tshape>
void VectorDirections()
{
    const int dimension = tshape::Dimension;
    const int numvectors = dimension * tshape::NumSides();
    int numnormalvectors = 0;
    int numfaces = tshape::NumSides(dimension-1);
    int numsides = tshape::NumSides();
    for (int side = numsides-numfaces-1; side < numsides-1; side++) {
        TPZStack<int> lower;
        tshape::LowerDimensionSides(side, lower);
        numnormalvectors += (lower.size()+1);
    }
    //    static void ComputeDirections(int is, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
    TPZFNMatrix<9> gradx(3,dimension,0.);
    for (int d=0; d<dimension; d++) gradx(d,d) = 1;
    TPZFNMatrix<3*3*27> directions(3,numvectors-numnormalvectors,0.);
    TPZManVector<int> sidevectors(numvectors-numnormalvectors);
    tshape::ComputeDirections(numsides-1, gradx, directions, sidevectors);

    TPZFNMatrix<27*3*3> NormalVectors(3,numvectors-numnormalvectors,0.);
    std::set<int> Mysides;
    int normalcount = 0;
    for (int is = 0; is < sidevectors.size(); is++) {
        int side = sidevectors[is];
        if (Mysides.find(side) != Mysides.end()) {
            continue;
        }
        Mysides.insert(side);
        int sidedim = tshape::SideDimension(side);
        TPZTransform tr = tshape::TransformSideToElement(side);
        for (int sd=0; sd < sidedim; sd++) {
            for (int d=0; d<dimension; d++) {
                NormalVectors(d,normalcount) = tr.Mult()(d,sd);
            }
            normalcount++;
        }
    }

//    directions.Print("Directions computed by topology",std::cout,EMathematicaInput);
//    NormalVectors.Print("Directions computed by actual test",std::cout,EMathematicaInput);
    
    // verify if the vectors are aligned
    for (int vec=0; vec < numvectors-numnormalvectors; vec++) {
        REAL inner = 0., norm = 0., normdir = 0.;
        for (int c=0; c<3; c++) {
            inner += directions(c,vec)*NormalVectors(c,vec);
            norm += NormalVectors(c,vec)*NormalVectors(c,vec);
            normdir += directions(c,vec)*directions(c,vec);
        }
        norm = sqrt(norm);
        normdir = sqrt(normdir);
        REAL diff = fabs(inner)/norm/normdir;
        diff -= 1.;
        BOOST_CHECK(fabs(diff) < 1.e-6);
    }
}

#endif
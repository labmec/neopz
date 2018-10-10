//
//  TestTopology.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/6/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "pzmanvector.h"
#include "pzvec_extras.h"
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
#include "TPZCompElDisc.h"
#include "tpzintpoints.h"
#include "pztrnsform.h"
#include "pzintel.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"

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
#include "pzshapeprism.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzshapepiramHdiv.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"


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

static int tetraedra_2[6][4]=
{
    {1,2,5,4},
    {4,7,3,2},
    {0,1,2,4},
    {0,2,3,4},
    {4,5,6,2},
    {4,6,7,2}
};

static bool MyDoubleComparer(REAL a, REAL b)
{
    if (IsZero(a-b)){
        return true;
    }
    else{
        return false;
    }
}

static void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem)
{
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (int64_t i=0; i<=nelem; i++) {
        for (int64_t j=0; j<=nelem; j++) {
            for (int64_t k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
}


static const int gfluxorder = 3;
static TPZAutoPointer<TPZCompMesh> GenerateMesh( TPZVec<TPZCompMesh *>  &meshvec,MElementType type, int nelem = 3, int fluxorder = gfluxorder, int ndiv = 0);

static TPZAutoPointer<TPZGeoMesh> /*TPZGeoMesh * */ CreateOneCuboWithTetraedrons(int nref);
static TPZAutoPointer<TPZGeoMesh> TetrahedralMeshCubo(int64_t nelem,int MaterialId);

static TPZAutoPointer<TPZGeoMesh> CreateGeoMeshHexaOfPir();
static TPZAutoPointer<TPZGeoMesh> CreateGeoMeshHexaOfPirTetra();


static int CompareShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);
static int CompareSideShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);

static TPZAutoPointer<TPZCompMesh> HDivMesh, PressureMesh;

static void CheckDRhamFacePermutations(MElementType Etype);
static void CheckDRhamPermutations(MElementType Etype);

template<class tshape>
static void CheckShapeOrder(int order);

static void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);

template<class tshape>
void VectorDirections();
//static TPZCompMesh *HDivMesh, *PressureMesh;

static void ExactPressure(const TPZVec<REAL> &x, TPZVec<STATE> &force)
{
    force[0] =  5. + 3. * x[0] + 2. * x[1] + 4. * x[0] * x[1];
}
static void ExactNormalFluxTop(const TPZVec<REAL> &x, TPZVec<STATE> &force)
{
    force[0] = 0.;
}

static void ExactNormalFluxBottom(const TPZVec<REAL> &x, TPZVec<STATE> &force)
{
    force[0] = 0.;
}

/// verify if the divergence of each vector function is included in the pressure space
static void CheckDRham(TPZCompEl *cel);

/// run a problem simulating a bilinear solution for the given element type
static void RunBilinear(MElementType eltype);

/// verify is the shape functions have continuity
static void VerifySideShapeContinuity(MElementType eltype);

/// verify if the pressure space is compatible with the flux space
static void VerifyDRhamCompatibility(MElementType eltype);

// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(mesh_tests)

BOOST_AUTO_TEST_CASE(vector_direction)
{
    std::cout << "Initializing vector_direction check\n";
    VectorDirections<pzshape::TPZShapePiram>();
    VectorDirections<pzshape::TPZShapeTetra>();
    VectorDirections<pzshape::TPZShapePrism>();
    VectorDirections<pzshape::TPZShapeCube>();
    VectorDirections<pzshape::TPZShapeTriang>();
    VectorDirections<pzshape::TPZShapeQuad>();
    std::cout << "Leaving vector_direction check\n";
}

BOOST_AUTO_TEST_CASE(sideshape_continuity)
{
    std::cout << "Initializing sideshape_continuity check\n";
    VerifySideShapeContinuity(EPiramide);
    VerifySideShapeContinuity(ETetraedro);
    VerifySideShapeContinuity(EPrisma);
    VerifySideShapeContinuity(ECube);
    VerifySideShapeContinuity(EQuadrilateral);
    VerifySideShapeContinuity(ETriangle);
    std::cout << "Leaving sideshape_continuity check\n";
}
    
    
BOOST_AUTO_TEST_CASE(shape_order)
{
    std::cout << "Initializing shape_order check\n";
    CheckShapeOrder<pzshape::TPZShapePiram>(6);
    CheckShapeOrder<pzshape::TPZShapeTetra>(6);
    CheckShapeOrder<pzshape::TPZShapeQuad>(6);
    CheckShapeOrder<pzshape::TPZShapeTriang>(6);
    CheckShapeOrder<pzshape::TPZShapeCube>(6);
    CheckShapeOrder<pzshape::TPZShapePrism>(6);
    std::cout << "Leaving shape_order check\n";
}
    

/// Check that the Div of the vector functions can be represented
BOOST_AUTO_TEST_CASE(drham_check)
{
    std::cout << "Initializing DRham consistency check\n";
    VerifyDRhamCompatibility(EPiramide);
    VerifyDRhamCompatibility(ETetraedro);
    VerifyDRhamCompatibility(EPrisma);
    VerifyDRhamCompatibility(ECube);
    VerifyDRhamCompatibility(EQuadrilateral);
    VerifyDRhamCompatibility(ETriangle);
    std::cout << "Leaving  DRham consistency check\n";
}

BOOST_AUTO_TEST_CASE(drham_permute_check)
{
    std::cout << "Initializing  DRham consistency under permutation check\n";
    CheckDRhamFacePermutations(EPiramide);
    CheckDRhamFacePermutations(ETetraedro);
    CheckDRhamFacePermutations(EPrisma);
    CheckDRhamFacePermutations(ECube);
    CheckDRhamPermutations(EQuadrilateral);
    CheckDRhamPermutations(ETriangle);
    std::cout << "Leaving  DRham consistency under permutation check\n";
}

/// Check that the Div of the vector functions can be represented
BOOST_AUTO_TEST_CASE(bilinearsolution_check)
{
    InitializePZLOG();
    std::cout << "Initializing solution check\n";
    RunBilinear(EPiramide);
    RunBilinear(ETetraedro);
    RunBilinear(EPrisma);
    RunBilinear(ETriangle);
    RunBilinear(EQuadrilateral);
    RunBilinear(ECube);
    std::cout << "Leaving solution check\n";
}

BOOST_AUTO_TEST_SUITE_END()

static TPZAutoPointer<TPZCompMesh> GenerateMesh( TPZVec<TPZCompMesh *>  &meshvec, MElementType eltype, int nelem, int fluxorder, int ndiv)
{
    int dimmodel = 2;
    TPZManVector<int,3> nx(2,nelem);
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = -1.;
    TPZGenGrid grid(nx,x0,x1);
    if (eltype == ETriangle|| eltype == EPrisma ) {
        grid.SetElementType(ETriangle);
    }
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    grid.Read(gmesh.operator->());
    grid.SetBC(gmesh, 4, -1);
    grid.SetBC(gmesh, 5, -1);
    grid.SetBC(gmesh, 6, -1);
    grid.SetBC(gmesh, 7, -1);
    
    if(eltype==ETriangle||eltype==EPrisma||eltype==ECube||eltype==EQuadrilateral )
    {
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
        
        
        {   // queria tanto ver a malha 2d
            std::ofstream Dummyfile("GeometricMesh2d.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(),Dummyfile, true);
        }
        
        
        
        if (eltype == EPrisma || eltype == ECube) {
            REAL thickness = 1.;//2.;
            TPZExtendGridDimension extend(gmesh,thickness);
            int numlayers = nelem;
            int bctop = -2;
            int bcbottom = -3 ;//normal negativa
            gmesh = extend.ExtendedMesh(numlayers,bcbottom,bctop);
            gmesh->SetDimension(3);
            dimmodel = 3;
        }
    }
    else if(eltype==ETetraedro)
    {
        // aqui
        dimmodel = 3;
        //gmesh = CreateOneCuboWithTetraedrons(ndiv); // AQUIDOUGLAS
        ndiv = 1;
        const int64_t NumberOfEl = ndiv;
        const int matid = 1;
        gmesh = TetrahedralMeshCubo(NumberOfEl, matid);
        gmesh->SetDimension(3);
        std::ofstream arg("gmesh.txt");
        gmesh->Print(arg);
        
    }
    else if(eltype==EPiramide)
    {
        // aqui
        dimmodel = 3;
        //gmesh = CreateOneCuboWithTetraedrons(ndiv); // AQUIDOUGLAS
        gmesh = CreateGeoMeshHexaOfPirTetra();
        std::ofstream arg("../gmesh.txt");
        gmesh->Print(arg);
        
    }
    else
    {
        // Elemento nao contemplado
        DebugStop();
    }
    
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream Dummyfile2("GeometricMesh3d.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile2, true);
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"Malha Geo FINAl \n\n";
        gmesh->Print(sout);
        //mphysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    int Axis;
    REAL theta, dump = 0.0;

    theta = 48.0;
    Axis = 1;
    RotateGeomesh(gmesh.operator->(), theta*dump, Axis);

    theta = -45.0;
    Axis = 2;
    RotateGeomesh(gmesh.operator->(), theta*dump, Axis);
    
    theta = 120.0;
    Axis = 3;
    RotateGeomesh(gmesh.operator->(), theta*dump, Axis);
    
    TPZMatPoisson3d *matpoisP = new TPZMatPoisson3d(1, dimmodel);
    TPZMaterial *poisP(matpoisP);
    
    PressureMesh = new TPZCompMesh(gmesh);
    
    TPZFNMatrix<4,STATE> val1(1,1,0.),val2(1,1,0.);
    PressureMesh->InsertMaterialObject(poisP);
    
    PressureMesh->SetAllCreateFunctionsContinuous();
    PressureMesh->ApproxSpace().CreateDisconnectedElements(true);
    
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
    
    PressureMesh->SetDefaultOrder(fluxorder);
    PressureMesh->SetDimModel(dimmodel);
    std::set<int> matids;
    matids.insert(1);
    const int64_t nel = gmesh->NElements();
    int64_t index;
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (gel->Type() != EPiramide){
            int64_t index;
            if (gel->Dimension() == gmesh->Dimension() && gel->MaterialId() == 1) {
                PressureMesh->ApproxSpace().CreateCompEl(gel, PressureMesh, index);
            }
        }
        else
        {
            new TPZIntelGen<TPZShapePiramHdiv>(PressureMesh,gel,index);
            //            cel->Print();
        }
        gel->ResetReference();
    }
    PressureMesh->ExpandSolution();
    int64_t ncon = PressureMesh->NConnects();
    for (int64_t ic=0; ic<ncon; ic++) {
        PressureMesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }

    if (eltype == ETriangle) { // Eprisma
        TPZCompElDisc::SetTotalOrderShape(PressureMesh.operator->());
    }
    
    PressureMesh->LoadReferences();
//    PressureMesh->ApproxSpace().CreateDisconnectedElements(false);
//    PressureMesh->AutoBuild();
    
    
    
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
    
    TPZCompEl::SetgOrder(fluxorder);
    HDivMesh->SetDefaultOrder(fluxorder);
    HDivMesh->SetDimModel(dimmodel);
    HDivMesh->AutoBuild();
    TPZCompMeshTools::AddHDivPyramidRestraints(HDivMesh.operator->());
    
//    {   int eqhdiv= HDivMesh->Solution().Rows();
//        TPZFMatrix<STATE> SolTriky(eqhdiv,1,1.0);
//        HDivMesh->LoadSolution(SolTriky);
//    }
    
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
        int64_t nel = mphysics->NElements();
        for (int64_t el=0; el<nel; el++) {
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
    
//    if (interface)
//    {
//        mphysics->Reference()->ResetReference();
//        mphysics->LoadReferences();
//        
//        int nel = mphysics->ElementVec().NElements();
//        for(int el = 0; el < nel; el++)
//        {
//            TPZCompEl * compEl = mphysics->ElementVec()[el];
//            if(!compEl) continue;
//            int index = compEl ->Index();
//            if(compEl->Dimension() == mphysics->Dimension())
//            {
//                TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(mphysics->ElementVec()[index]);
//                if(!InterpEl) continue;
//                InterpEl->CreateInterfaces();
//            }
//        }
//        
//    }
    
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();

//    ofstream arq4("mphysics.txt");
//    mphysics->Print(arq4);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
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
        TPZTransform<> tr = gelsideA.NeighbourSideTransform(gelsideB);
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
            if(fabs(Aval-Bval) > 1.e-6)
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
    TPZTransform<> tr = gelsideA.NeighbourSideTransform(gelsideB);
    TPZGeoEl *gelA = gelsideA.Element();
    TPZTransform<> trA = gelA->SideToSideTransform(gelsideA.Side(), gelA->NSides()-1);
    TPZGeoEl *gelB = gelsideB.Element();
    TPZTransform<> trB = gelB->SideToSideTransform(gelsideB.Side(), gelB->NSides()-1);
    
    int vecZeroA = -1;
    if(gelA->Type() == EPiramide)
    {
        TPZOneShapeRestraint restA = *interA->GetShapeRestraints().begin();
        int indexA = restA.fFaces[0].first;
        vecZeroA = 11 + (indexA - 1)*7;
    }
    
    int vecZeroB = -1;
    if(gelB->Type() == EPiramide)
    {
        TPZOneShapeRestraint restB = *interB->GetShapeRestraints().begin();
        int indexB = restB.fFaces[0].first;
        vecZeroB = 11 + (indexB - 1)*7;
    }
    
    
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
        TPZFNMatrix<9> jacobian(sidedim,sidedim),jacinv(sidedim,sidedim),axes(sidedim,3,0.);
        REAL detjac;
        gelsideA.Jacobian(pointA, jacobian, axes, detjac, jacinv);
        TPZManVector<REAL,3> normal(3,0.), xA(3),xB(3);
        if(axes.Rows() == 2)
        {
            
            for(int i=0; i<3; i++)
            {
                xA[i] = axes(0,i);
                xB[i] = axes(1,i);
            }
            Cross(xA,xB,normal);
        }
        else if(axes.Rows() == 1)
        {
            normal[1] = -axes(0,0);
            normal[0] = axes(0,1);
        }
        else
        {
            DebugStop();
        }
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
            REAL vecnormalA = 1.;
            REAL vecnormalB = 1.;
            // if A or B are boundary elements, their shapefunctions come in the right order
            if (dimensionA != sidedim) {
                Ashapeind = dataA.fVecShapeIndex[i].second;
                Avecind = dataA.fVecShapeIndex[i].first;
                vecnormalA = dataA.fNormalVec(0,Avecind)*normal[0]+dataA.fNormalVec(1,Avecind)*normal[1]+dataA.fNormalVec(2,Avecind)*normal[2];
            }
            if (dimensionB != sidedim) {
                Bshapeind = dataB.fVecShapeIndex[j].second;
                Bvecind = dataB.fVecShapeIndex[j].first;
                vecnormalB = dataB.fNormalVec(0,Bvecind)*normal[0]+dataB.fNormalVec(1,Bvecind)*normal[1]+dataB.fNormalVec(2,Bvecind)*normal[2];
                
            }
            if (dimensionA != sidedim && dimensionB != sidedim) {
                // vefify that the normal component of the normal vector corresponds
                Avecind = dataA.fVecShapeIndex[i].first;
                Bvecind = dataB.fVecShapeIndex[j].first;
                vecnormalA = dataA.fNormalVec(0,Avecind)*normal[0]+dataA.fNormalVec(1,Avecind)*normal[1]+dataA.fNormalVec(2,Avecind)*normal[2];
                vecnormalB = dataB.fNormalVec(0,Bvecind)*normal[0]+dataB.fNormalVec(1,Bvecind)*normal[1]+dataB.fNormalVec(2,Bvecind)*normal[2];
                if(fabs(vecnormalA-vecnormalB) > 1.e-6)
                {
                    nwrong++;
                    LOGPZ_ERROR(logger, "normal vectors aren't equal")
                }

            }
            shapesA[i-firstShapeA] = dataA.phi(Ashapeind,0);
            shapesB[j-firstShapeB] = dataB.phi(Bshapeind,0);
            REAL valA = dataA.phi(Ashapeind,0)*vecnormalA;
            REAL valB = dataB.phi(Bshapeind,0)*vecnormalB;
            REAL diff = valA-valB;
            REAL decision = fabs(diff)-1.e-6;
            if(decision > 0. && vecZeroA != Avecind && vecZeroB != Bvecind)
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
    int porder = cel->GetgOrder();
    std::string filename;
    {
        std::stringstream sout;
        sout << "../matrices" << porder << ".nb";
        filename = sout.str();
    }
    
    std::ofstream output(filename.c_str());
    {
        std::stringstream sout;
        sout << "L2" << porder << " = ";
        filename = sout.str();
    }
    L2->Print(filename.c_str(),output, EMathematicaInput);
    {
        std::stringstream sout;
        sout << "PressHDiv" << porder << " = ";
        filename = sout.str();
    }
    inner.Print(filename.c_str(),output,EMathematicaInput);
    TPZStepSolver<STATE> step(L2);
    step.SetDirect(ELU);
    step.Solve(inner,multiplier);
    {
        std::stringstream sout;
        sout << "multipl" << porder << " = ";
        filename = sout.str();
    }

    multiplier.Print(filename.c_str(),output,EMathematicaInput);
    output.close();
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
    TPZInterpolationSpace *intelP = dynamic_cast<TPZInterpolationSpace *>(celMF->Element(1));
    if (!intel || ! intelP) {
        DebugStop();
    }
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
        intel->ComputeRequiredData(dataA, pos);
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
    TPZInterpolationSpace *intelP = dynamic_cast<TPZInterpolationSpace *>(celMF->Element(1));
    
    if (!intelP || !intel) {
        DebugStop();
    }
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
        intel->ComputeRequiredData(dataA, pos);
        intelP->ComputeShape(pos, dataB);
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
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
            STATE phival = 0;

            for (ish=0; ish<npressure; ish++) {

                phival += multiplier(ish,jsh)*dataB.phi(ish);
            }
            // the divergence of the vector function should be equal to the value of projected pressure space
            STATE diff = phival-divphi;
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
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
                DebugStop();
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
    TPZManVector<int64_t,tshape::NSides> origids(tshape::NCornerNodes,0);
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
    TPZManVector<TPZManVector<int,8> > permutation;
    MElementType thistype = tshape::Type();
    TPZRefPatternTools::GetElTypePermutations(thistype, permutation);
    const int nperm = permutation.size();
    
    // set up the integration rule
    TPZInt1d localrule(20);
    
    for (int iperm = 0; iperm<nperm; iperm++)
    {
        int numwrong = 0;
        TPZManVector<int64_t,tshape::NCornerNodes> ids(tshape::NCornerNodes,0);
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
        TPZVec<int64_t> sides;
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
            TPZManVector<int64_t,8> locids(tshape::NSideNodes(is));
            for (int i=0; i<locids.size(); i++) {
                locids[i] = ids[tshape::SideNodeLocId(is, i)];
            }
            if(sidetype == EPiramide && sidedim == 3)
            {
                for(int shape = 0; shape < nsideshape; shape++)
                {
                    for(int d=0; d<3; d++) estimatedshapeorders(shapecounter+shape,d) = shapeorders(shapecounter+shape,d);
                }
                continue;
            }
            /// Estimate the shape orders in each direction
            for (int dim = 0; dim < sidedim; dim++) {
                TPZManVector<REAL,3> sidecenterel(dimension),sidecenter(sidedim),point(sidedim);
                TPZFNMatrix<9> jac(2,1);
                tshape::CenterPoint(is, sidecenterel);
//                if(sidetype==EPrisma && is>19){continue;} ///???????
//                else{tshape::MapToSide(is, sidecenterel, sidecenter, jac);}
                tshape::MapToSide(is, sidecenterel, sidecenter, jac);
                for (int i=0; i<sidedim; i++) {
                    sidecenter[i] += M_PI/230.;
                }
                TPZFNMatrix<500,REAL> integrationvals(numlegendre,nsideshape,0.);
                TPZManVector<REAL,1> pos(1);
                REAL weight;
                for (int ip=0; ip<localrule.NPoints(); ip++) {
                    localrule.Point(ip, pos, weight);
                    point = sidecenter;
                    point[dim]=pos[0];
//                    point[dim] = (pos[0]+1.)/2.;
                    if (sidetype == ETriangle|| (sidetype == EPrisma && dim<sidedim-1)  || sidetype == EPiramide ) {
                        REAL a,b;
                        a = (point[0]-sidecenter[0])*cos(M_PI/20.)-(point[1]-sidecenter[1])*sin(M_PI/20.);
                        b = (point[0]-sidecenter[0])*sin(M_PI/20.)+(point[1]-sidecenter[1])*cos(M_PI/20.);
                        point[0] = a;
                        point[1] = b;
                    }
                    if ( sidetype == ETetraedro ) {
                        REAL a,b,c,p0,p1,p2;
                        p0 = (point[0]-sidecenter[0])/3.;
                        p1 = (point[1]-sidecenter[1])/3.;
                        p2 = (point[2]-sidecenter[2])/3.;
                        // rotacao nos 3 eixos
                        REAL coss = cos(M_PI/20.), senn = sin(M_PI/20.);
                        a = p0*coss*coss  +  p1*(coss*senn - coss*senn*senn ) + p2*(coss*coss*senn + senn*senn);
                        b = -p0*(coss*senn) + p1*(coss*coss + senn*senn*senn) + p2*(coss*senn - coss*senn*senn);
                        c = -p0*senn - p1*(coss*senn) + p2*coss*coss;
                        // Rodar em z tambem
                        // aqui
                        point[0] = a+sidecenter[0];
                        point[1] = b+sidecenter[1];
                        point[2] = c+sidecenter[2];
                    }

                    TPZFNMatrix<200,REAL> philegendre(numlegendre,1,0.);
                    TPZFNMatrix<200,REAL> dphi(sidedim,nsideshape+firstshape), dphilegendre(1,numlegendre,0.), phi(nsideshape+firstshape,1,0.);
                    tshape::SideShape(is, point, locids, orders, phi, dphi);
                    
//                    for (int ishape = firstshape; ishape<firstshape+nsideshape; ishape++) {
//                        std::cout << phi(ishape,0) << " ";
//                    }
//                    std::cout << std::endl;
                    
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
    HDivPiola = 1;
    const int dimension = tshape::Dimension;
    int numvectors = dimension * tshape::NumSides();
    if(tshape::Type() == EPiramide) numvectors++;
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
    for (int d=0; d<dimension; d++) gradx(d,d) = 1.;
    TPZFNMatrix<3*3*27> directions(3,numvectors-numnormalvectors,0.);
    TPZFNMatrix<3*3*27> directionsAll(3,numvectors,0.);
    TPZManVector<int> sidevectors(numvectors-numnormalvectors);
    REAL detjac = 1.0;
//    tshape::ComputeDirections(numsides-1, gradx, directions, sidevectors);
//    // Nao tem que usar esse metodo?
//    tshape::ComputeDirections(gradx, detjac, directionsAll);
    tshape::ComputeDirections(numsides-1, gradx, directions, sidevectors);
    tshape::ComputeDirections(gradx, detjac, directionsAll);

    // copia dos vetroes internos por face
    int numintvec = numvectors-numnormalvectors;
    for (int j = 0; j< numintvec; j++) {
        directions(0,j) = directionsAll(0,j+numnormalvectors);
        directions(1,j) = directionsAll(1,j+numnormalvectors);
        directions(2,j) = directionsAll(2,j+numnormalvectors);
    }
    

    TPZFNMatrix<27*3*3> NormalVectors(3,numintvec,0.);
    std::set<int> Mysides;
    int normalcount = 0;
    for (int is = 0; is < sidevectors.size(); is++) {
        int side = sidevectors[is];
        // Repeated elements out of Mysides
        if (Mysides.find(side) != Mysides.end()) {
            continue;
        }
        Mysides.insert(side);
        int sidedim = tshape::SideDimension(side);
        TPZTransform<> tr = tshape::TransformSideToElement(side);
        for (int sd=0; sd < sidedim; sd++) {
            for (int d=0; d<dimension; d++) {
                NormalVectors(d,normalcount) = tr.Mult()(d,sd);
            }
            normalcount++;
        }
    }
    
//    cout << "sidevectors " << sidevectors << endl;
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

void CheckDRhamFacePermutations(MElementType eltype)
{
    TPZVec<TPZCompMesh *>  meshvec(2);
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype,3,4,0);
    TPZGeoMesh *gmesh = cmesh->Reference();
    int meshdim = cmesh->Dimension();

    const int64_t nel = gmesh->NElements();
    int64_t nel3d = 0;
    for (int64_t el=0; el<nel; el++) {
        if (gmesh->ElementVec()[el]->Dimension() == meshdim) {
            nel3d++;
        }
    }
    TPZGeoEl *gel = gmesh->ElementVec()[nel3d/2];
    TPZCompEl *cel = gel->Reference();
    TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if(!intel)
    {
        DebugStop();
    }
    const int gelcorner = gel->NCornerNodes();
    TPZManVector<int64_t,8> nodeids(gelcorner,0), nodesperm(gelcorner,0);
    for (int i = 0; i<gelcorner; i++) {
        nodeids[i] = gel->NodePtr(i)->Id();
    }
    int dimension = gel->Dimension();
    int nfaces = 0;
    for (int side = 0; side<dimension; side++) {
        if (gel->SideDimension(side) == dimension-1) {
            nfaces++;
        }
    }
    
    int nsides = gel->NSides();
    TPZManVector<std::set<int>, 27> verifiedperms(nsides);
    
    int64_t permcounter = 0;
    for (int side = 0; side < nsides; side++) {
        if (gel->SideDimension(side) != dimension-1) {
            continue;
        }
        for (int i = 0; i<gelcorner; i++) {
            gel->NodePtr(i)->SetNodeId(nodeids[i]);
        }
        int ncorner = gel->NSideNodes(side);
        TPZPermutation perm(ncorner);
        do {
            TPZManVector<int64_t> cornerids(ncorner);
            for (int ic = 0; ic<ncorner; ic++) {
                int locindex = gel->SideNodeLocIndex(side, ic);
                cornerids[ic] = nodeids[locindex];
            }
            perm.Permute(cornerids, nodesperm);
            for (int ic = 0; ic<ncorner; ic++) {
                int locindex = gel->SideNodeLocIndex(side, ic);
                int id = nodesperm[locindex];
                gel->NodePtr(locindex)->SetNodeId(id);
            }
            
            // verify the tranformation ids of all faces, starting with side
            int transid = 0;
            if (ncorner == 4) {
                transid = TPZShapeQuad::GetTransformId2dQ(nodesperm);
            }
            else if (ncorner == 3) {
                transid = TPZShapeTriang::GetTransformId2dT(nodesperm);
            }
            else if (ncorner == 2) {
                transid = TPZShapeLinear::GetTransformId1d(nodesperm);
            }
            else
            {
                DebugStop();
            }
            if (verifiedperms[side].find(transid) != verifiedperms[side].end()) {
                perm++;
                continue;
            }

            
            for (int is=0; is<nsides; is++) {
                if (gel->SideDimension(is) != dimension-1) {
                    continue;
                }
                int nc = gel->NSideNodes(is);
                TPZManVector<int64_t,4> cornids(nc);
                for (int ic=0; ic<nc; ic++) {
                    int locid = gel->SideNodeLocIndex(is, ic);
                    cornids[ic] = gel->NodePtr(locid)->Id();
                }
                if (nc == 4) {
                    int transid = TPZShapeQuad::GetTransformId2dQ(cornids);
                    verifiedperms[is].insert(transid);
                }
                else if(nc == 3)
                {
                    int transid = TPZShapeTriang::GetTransformId2dT(cornids);
                    verifiedperms[is].insert(transid);
                }
                else if(nc == 2)
                {
                    int transid = TPZShapeLinear::GetTransformId1d(cornids);
                    verifiedperms[is].insert(transid);
                }
                else
                {
                    DebugStop();
                }
            }
            
            CheckDRham(intel);
            // now compare the shape function value between neighbrouring elements
            // create a side integration rule
            // compute the transformation between the element and his neighbour
            // compute the vectors of shape function
            // compare
            perm++;
            permcounter++;
        } while (!perm.IsFirst() && permcounter < 150);
    }
}

void CheckDRhamPermutations(MElementType eltype)
{
    TPZVec<TPZCompMesh *>  meshvec(2);
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype);
    TPZGeoMesh *gmesh = cmesh->Reference();
    int meshdim = cmesh->Dimension();
    
    const int64_t nel = gmesh->NElements();
    int64_t nel3d = 0;
    for (int64_t el=0; el<nel; el++) {
        if (gmesh->ElementVec()[el]->Dimension() == meshdim) {
            nel3d++;
        }
    }
    TPZGeoEl *gel = gmesh->ElementVec()[nel3d/2];
    int dimension = gel->Dimension();
    if (dimension != meshdim) {
        DebugStop();
    }

    
    TPZCompEl *cel = gel->Reference();
    TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if(!intel)
    {
        DebugStop();
    }

    const int gelcorner = gel->NCornerNodes();
    TPZManVector<int,8> nodeids(gelcorner,0), nodesperm(gelcorner,0);
    for (int i = 0; i<gelcorner; i++) {
        nodeids[i] = gel->NodePtr(i)->Id();
    }
    int64_t permcounter = 0;
    TPZPermutation perm(gelcorner);
    do {
        TPZManVector<int> cornerids(gelcorner);
        perm.Permute(cornerids, nodesperm);
        for (int ic = 0; ic<gelcorner; ic++) {
            gel->NodePtr(ic)->SetNodeId(nodesperm[ic]);
        }
        CheckDRham(intel);
        // now compare the shape function value between neighbrouring elements
        // create a side integration rule
        // compute the transformation between the element and his neighbour
        // compute the vectors of shape function
        // compare
        perm++;
        permcounter++;
    } while (!perm.IsFirst() && permcounter < 150);
}

/// run a problem simulating a bilinear solution for the given element type
void RunBilinear(MElementType eltype)
{
    int nelx = 1;
    int fluxorder = 1;
    if (eltype == ETriangle) {
        fluxorder = 3;
    }
    if (eltype == EPrisma) {
        fluxorder = 2;
    }
    if (eltype == ETetraedro) {
        fluxorder = 2;
    }
    int ndiv = 0; // para refinar a malha
    TPZVec<TPZCompMesh *>  meshvec(2);
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype,nelx,fluxorder, ndiv);
    
    
    {
        TPZMaterial *mat = cmesh->FindMaterial(-1);
        if(!mat) DebugStop();
        TPZDummyFunction<STATE> *dumforce = new TPZDummyFunction<STATE>(ExactPressure, 5);
        TPZAutoPointer<TPZFunction<STATE> > autofunc (dumforce);
        mat->SetForcingFunction(autofunc);
    }
    {
        TPZMaterial *mat = cmesh->FindMaterial(-2);
        if(!mat) DebugStop();
        TPZDummyFunction<STATE> *dumforce = new TPZDummyFunction<STATE>(ExactNormalFluxTop, 5);
        TPZAutoPointer<TPZFunction<STATE> > autofunc (dumforce);
        mat->SetForcingFunction(autofunc);
    }
    
    {
        TPZMaterial *mat = cmesh->FindMaterial(-3);
        if(!mat) DebugStop();
        TPZDummyFunction<STATE> *dumforce = new TPZDummyFunction<STATE>(ExactNormalFluxBottom, 5);
        TPZAutoPointer<TPZFunction<STATE> > autofunc (dumforce);
        mat->SetForcingFunction(autofunc);
    }
    
    
//    {
//        //  Print Geometrical Base Mesh
//        std::ofstream Dummyfile("GeometricMesh.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(cmesh,Dummyfile, true);
//    }
    
    // cmesh->SaddlePermute();
    TPZAnalysis an(cmesh,false);
    // para resolver o sistema
    // escolhe entre isso
//    TPZFStructMatrix str(cmesh);
//    an.SetStructuralMatrix(str);
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELU);
    /// ou isso
    TPZSkylineStructMatrix str(cmesh);
    
//    TPZFMatrix<STATE> rhs, solteste;
//    TPZAutoPointer<TPZGuiInterface> guiInterface;
//    TPZAutoPointer<TPZMatrix<STATE> > matrix = str.CreateAssemble(rhs, guiInterface);
////    matrix->Print(std::cout,EMathematicaInput);
//    matrix->Print("EK = ", cout ,EMathematicaInput);
//    rhs.Print("EF = ", cout ,EMathematicaInput);
//    solteste = an.Solution();
//    solteste.Print("Sol = ", cout ,EMathematicaInput);
    
    an.SetStructuralMatrix(str);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    
    
    an.SetSolver(step);
    an.Run();
    
    
//    an.Solution().Print("Solucao");
    
    std::string plotfile("GSaida.vtk");
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmesh.operator->());
    TPZManVector<std::string,10> scalnames(1), vecnames(1);
    vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
    const int dim = cmesh->Dimension();
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Reference()->Print(sout);
        sout <<  "computacional" << std::endl;
        cmesh->Print(sout);
        //a matriz??
        an.Rhs().Print("Right Hand Side",sout);
        an.Solution().Print("Solution",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int nel = 1;//cmesh->NElements();
    for(int i = 0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension()<cmesh->Dimension()) continue;
        int ns = gel->NSides();
        TPZIntPoints *rule = gel->CreateSideIntegrationRule(ns-1, 4);//3
        int np = rule->NPoints();
        for(int ip=0; ip<np; ip++)
        {
            //TPZManVector<REAL,3> xi(2), xco(3), sol(1), exactsol(1);
            TPZManVector<REAL,3> xi(3), xco(3);
            TPZManVector<STATE,3> sol(1), exactsol(1);
            REAL weight;
            rule->Point(ip, xi, weight);
            gel->X(xi, xco);
            cel->Solution(xi, 2, sol);
            ExactPressure(xco,exactsol);
            //        autofunc->Execute(xco, exactsol);
            if (fabs(sol[0]-exactsol[0]) > 1.e-6) {
                std::cout << "xi = " << xi <<std::endl;
                std::cout << "xco = " << xco <<std::endl;
            }
            BOOST_CHECK(fabs(sol[0]-exactsol[0]) < 1.e-6);
        }
    }
}

/// verify is the shape functions have continuity
void VerifySideShapeContinuity(MElementType eltype)
{
    int64_t permcount = 0;
    TPZVec<TPZCompMesh *>  meshvec(2);
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype);
    TPZGeoMesh *gmesh = cmesh->Reference();
    TPZGeoEl *gel = gmesh->ElementVec()[0];
    if(gel->Type() != eltype)
    {
        DebugStop();
    }
    const int gelcorner = gel->NCornerNodes();
    TPZManVector<int,8> nodeids(gelcorner,0);
    for (int i = 0; i<gelcorner; i++) {
        nodeids[i] = gel->NodePtr(i)->Id();
    }
    int dimension = gel->Dimension();
    for (int side = 0; side < gel->NSides(); side++) {
        if (gel->SideDimension(side) != dimension-1) {
            continue;
        }
        for (int i = 0; i<gelcorner; i++) {
            gel->NodePtr(i)->SetNodeId(nodeids[i]);
        }
        int ncorner = gel->NSideNodes(side);
        TPZPermutation perm(ncorner);
        do {
            TPZManVector<int> cornerids(ncorner), nodesperm(ncorner,0);
            for (int ic = 0; ic<ncorner; ic++) {
                int locindex = gel->SideNodeLocIndex(side, ic);
                cornerids[ic] = nodeids[locindex];
            }
            perm.Permute(cornerids, nodesperm);
            for (int ic = 0; ic<ncorner; ic++) {
                int locindex = gel->SideNodeLocIndex(side, ic);
                int id = nodesperm[ic];
                gel->NodePtr(locindex)->SetNodeId(id);
            }
            TPZManVector<int64_t,8> gelids(gelcorner);
            for (int i=0; i< gelcorner; i++) {
                gelids[i] = gel->NodePtr(i)->Id();
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
        } while (!perm.IsFirst() && permcount < 150);
    }
}

/// verify if the pressure space is compatible with the flux space
void VerifyDRhamCompatibility(MElementType eltype)
{
    // generate a mesh
    TPZVec<TPZCompMesh *>  meshvec(2);
    int nelem =2;
    int fluxorder = gfluxorder;
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateMesh(meshvec,eltype,nelem,fluxorder);
    std::ofstream arg1("cmesh.txt");
    cmesh.operator->()->Print(arg1);
    // for each computational element (not boundary) verify if the Div(vecspace) is included in the pressure space
    int nel = cmesh->NElements();
    int meshdim = cmesh->Dimension();
    int iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!intel)
        {
            DebugStop();
        }
        if(intel->Reference()->Dimension() != meshdim) continue;
        CheckDRham(intel);
    }
}

TPZAutoPointer<TPZGeoMesh> TetrahedralMeshCubo(int64_t nelem,int MaterialId)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh,nelem);
    
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el=0; el<6; el++)
                {
                    TPZManVector<int64_t,4> elnodes(4);
                    int64_t index;
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra_2[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index,0);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    
    // Boundary Conditions
    const int numelements = gmesh->NElements();
    const int bczMinus = -3, bczplus = -2, bcids = -1;
//    const int bczMinus = -1, bczplus = -1, bcids = -1;
    
    for(int el=0; el<numelements; el++)
    {
        TPZManVector <TPZGeoNode,4> Nodefinder(4);
        TPZManVector <REAL,3> nodecoord(3);
        TPZGeoEl *tetra = gmesh->ElementVec()[el];
        // na face x = 0
        TPZVec<int64_t> ncoordVec(0); int64_t sizeOfVec = 0;
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bcids);	
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 1
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bcids);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 0
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bcids);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 1
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bcids);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face z = 0
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bczMinus);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face z = 1
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,bczplus);
        }
        
        
        
    }
    
    return gmesh;
}




TPZAutoPointer<TPZGeoMesh> /*TPZGeoMesh * */ CreateOneCuboWithTetraedrons(int nref)
{
    
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int nnodes = 8;
    //para que as condicoes de contorno fiquem como nos testes dos outros elementos, houve mudanca nos indices
    
    int idf0=-3;//-1;
    int idf1=-1;//-2;
    int idf2=-1;//-3;
    int idf3=-1;//-4;
    int idf4=-1;//-5;
    int idf5=-2;//-6;
    
    int matId = 1;
    
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    
    //cubo [0,1]ˆ3
    //c0
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c1
    coord[0] =  1.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c2
    coord[0] =  1.0;
    coord[1] =  1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c3
    coord[0] = 0.0;
    coord[1] =  1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    
    //c4
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c5
    coord[0] =  1.0;
    coord[1] = 0.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c6
    coord[0] =  1.0;
    coord[1] =  1.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c7
    coord[0] = 0.0;
    coord[1] =  1.0;
    coord[2] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    
    
    // cubo [-1,1]^3
    //    //c0
    //    coord[0] = -1.0;
    //    coord[1] = -1.0;
    //    coord[2] = -1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c1
    //    coord[0] =  1.0;
    //    coord[1] = -1.0;
    //    coord[2] = -1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c2
    //    coord[0] =  1.0;
    //    coord[1] =  1.0;
    //    coord[2] = -1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c3
    //    coord[0] = -1.0;
    //    coord[1] =  1.0;
    //    coord[2] = -1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c4
    //    coord[0] = -1.0;
    //    coord[1] = -1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c5
    //    coord[0] =  1.0;
    //    coord[1] = -1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c6
    //    coord[0] =  1.0;
    //    coord[1] =  1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c7
    //    coord[0] = -1.0;
    //    coord[1] =  1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    
    
    
    int index = 0;
    
    TPZManVector<int64_t,4> TopolTetra(4,0);
    TopolTetra[0] = 0;
    TopolTetra[1] = 1;
    TopolTetra[2] = 3;
    TopolTetra[3] = 4;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matId, *gmesh);
    index++;
    
    TopolTetra[0] = 1;
    TopolTetra[1] = 2;
    TopolTetra[2] = 3;
    TopolTetra[3] = 6;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matId, *gmesh);
    index++;
    
    TopolTetra[0] = 1;
    TopolTetra[1] = 5;
    TopolTetra[2] = 4;
    TopolTetra[3] = 6;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matId, *gmesh);
    index++;
    
    TopolTetra[0] = 3;
    TopolTetra[1] = 7;
    TopolTetra[2] = 6;
    TopolTetra[3] = 4;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matId, *gmesh);
    index++;
    
    TopolTetra[0] = 1;
    TopolTetra[1] = 3;
    TopolTetra[2] = 4;
    TopolTetra[3] = 6;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matId, *gmesh);
    index++;

    
    TPZVec<int64_t> TopolTriang(3);
    
    // bottom
    TopolTriang[0] = 0;
    TopolTriang[1] = 1;
    TopolTriang[2] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (index,TopolTriang,idf0,*gmesh);
    index++;
    
    TopolTriang[0] = 1;
    TopolTriang[1] = 2;
    TopolTriang[2] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (index,TopolTriang,idf0,*gmesh);
    index++;
    
    // Front
    TopolTriang[0] = 0;
    TopolTriang[1] = 1;
    TopolTriang[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (index,TopolTriang,idf1,*gmesh);
    index++;
    
    TopolTriang[0] = 1;
    TopolTriang[1] = 5;
    TopolTriang[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (index,TopolTriang,idf1,*gmesh);
    index++;
    
    // Rigth
    TopolTriang[0] = 1;
    TopolTriang[1] = 2;
    TopolTriang[2] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (index,TopolTriang,idf2,*gmesh);
    index++;
    
    TopolTriang[0] = 1;
    TopolTriang[1] = 5;
    TopolTriang[2] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (index,TopolTriang,idf2,*gmesh);
    index++;
    
    // Back
    TopolTriang[0] = 2;
    TopolTriang[1] = 3;
    TopolTriang[2] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (index,TopolTriang,idf3,*gmesh);
    index++;
    
    TopolTriang[0] = 3;
    TopolTriang[1] = 6;
    TopolTriang[2] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (index,TopolTriang,idf4,*gmesh);
    index++;
    
    // Left
    TopolTriang[0] = 0;
    TopolTriang[1] = 3;
    TopolTriang[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (index,TopolTriang,idf4,*gmesh);
    index++;
    
    TopolTriang[0] = 3;
    TopolTriang[1] = 4;
    TopolTriang[2] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (index,TopolTriang,idf4,*gmesh);
    index++;
    
    // Top
    TopolTriang[0] = 4;
    TopolTriang[1] = 5;
    TopolTriang[2] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (index,TopolTriang,idf5,*gmesh);
    index++;
    
    TopolTriang[0] = 4;
    TopolTriang[1] = 6;
    TopolTriang[2] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (index,TopolTriang,idf5,*gmesh);
    index++;
    
    
    
    gmesh->BuildConnectivity();
    
    /// gmesh para aqui
    
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("SingleCubeTetraWithBcs.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<REAL> RotationMatrix(3,3,0.0);

    switch (Axis) {
        case 1:
        {
            RotationMatrix(0,0) = 1.0;
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(1,2) =   -sin(theta);
            RotationMatrix(2,1) =   +sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 2:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,2) =   +sin(theta);
            RotationMatrix(1,1) = 1.0;
            RotationMatrix(2,0) =   -sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 3:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
        default:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
    }
    
    TPZVec<REAL> iCoords(3,0.0);
    TPZVec<REAL> iCoordsRotated(3,0.0);
    
    //RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

    TPZAutoPointer<TPZGeoMesh>  CreateGeoMeshHexaOfPir()
    {
        const int dim = 3;
        TPZGeoMesh *gmesh = new TPZGeoMesh;
        gmesh->SetDimension(dim);
        
        // Setando os nohs
        int nnodes = 9;
        gmesh->NodeVec().Resize(nnodes);
        int ino = 0;
        const int matid = 1;
        int64_t index = 0;
        
        // noh 0
        TPZManVector<REAL, 3> nodecoord(3,0.);
        nodecoord[0] = -1.;
        nodecoord[1] = -1.;
        nodecoord[2] = -1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 1
        nodecoord[0] = 1.;
        nodecoord[1] = -1.;
        nodecoord[2] = -1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 2
        nodecoord[0] = 1.;
        nodecoord[1] = 1.;
        nodecoord[2] = -1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 3
        nodecoord[0] = -1.;
        nodecoord[1] = 1.;
        nodecoord[2] = -1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 4
        nodecoord[0] = -1.;
        nodecoord[1] = -1.;
        nodecoord[2] = 1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 5
        nodecoord[0] = 1.;
        nodecoord[1] = -1.;
        nodecoord[2] = 1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 6
        nodecoord[0] = 1.;
        nodecoord[1] = 1.;
        nodecoord[2] = 1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 7
        nodecoord[0] = -1.;
        nodecoord[1] = 1.;
        nodecoord[2] = 1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 8
        nodecoord[0] = 0.;
        nodecoord[1] = 0.;
        nodecoord[2] = 0.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        
        // Criando elemento
        TPZManVector<int64_t,5> topolPyr(5);
        int myels[6][5] = {{0,1,5,4,8},{1,2,6,5,8},{2,3,7,6,8},{0,3,7,4,8},{0,1,2,3,8},{4,5,6,7,8}};
        //int myels[6][5] = {{0,1,5,4,8},{6,5,1,2,8},{2,3,7,6,8},{7,4,0,3,8},{0,1,2,3,8},{4,5,6,7,8}}; //Sequencia trocada soh para funcionar o AddHDivPyramidRestraints
        for (int iel = 0; iel < 6; iel++) {
            for (int i = 0; i < 5; i++) {
                topolPyr[i] = myels[iel][i];
            }
            gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index,0);
        }
        
        const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
        
        const int64_t nel = gmesh->NElements();
        for (int64_t iel = 0; iel < nel; iel++) {
            gmesh->Element(iel)->CreateBCGeoEl(13, bc0);
        }
        
        gmesh->BuildConnectivity();
        
        std::ofstream out("HexaPyrGmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        
        return gmesh;
    }

    TPZAutoPointer<TPZGeoMesh> CreateGeoMeshHexaOfPirTetra()
    {
        const int dim = 3;
        TPZGeoMesh *gmesh = new TPZGeoMesh;
        gmesh->SetDimension(dim);
        
        // Setando os nohs
        int nnodes = 9;
        gmesh->NodeVec().Resize(nnodes);
        int ino = 0;
        const int matid = 1;
        int64_t index = 0;
        
        // noh 0
        TPZManVector<REAL, 3> nodecoord(3,0.);
        nodecoord[0] = -1.;
        nodecoord[1] = -1.;
        nodecoord[2] = -1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 1
        nodecoord[0] = 1.;
        nodecoord[1] = -1.;
        nodecoord[2] = -1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 2
        nodecoord[0] = 1.;
        nodecoord[1] = 1.;
        nodecoord[2] = -1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 3
        nodecoord[0] = -1.;
        nodecoord[1] = 1.;
        nodecoord[2] = -1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 4
        nodecoord[0] = -1.;
        nodecoord[1] = -1.;
        nodecoord[2] = 1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 5
        nodecoord[0] = 1.;
        nodecoord[1] = -1.;
        nodecoord[2] = 1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 6
        nodecoord[0] = 1.;
        nodecoord[1] = 1.;
        nodecoord[2] = 1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 7
        nodecoord[0] = -1.;
        nodecoord[1] = 1.;
        nodecoord[2] = 1.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        ino++;
        
        // noh 8
        nodecoord[0] = 0.;
        nodecoord[1] = 0.;
        nodecoord[2] = 0.;
        gmesh->NodeVec()[ino].SetCoord(nodecoord);
        gmesh->NodeVec()[ino].SetNodeId(ino);
        //    ino++;
        gmesh->SetNodeIdUsed(ino);
        
        // Criando elemento
        TPZManVector<int64_t,5> topolPyr(5), topolTet(4), topolTri(3);
        int myelsp[2][5] = {{4,0,2,6,1},{4,0,2,6,7}};
        int myelst[2][4] = {{4,6,5,1},{0,2,3,7}};
        //                          front           right          top             back            left            bottom
        int triangles[12][3] = {{0,1,4},{1,5,4},{1,2,6},{1,6,5},{4,5,6},{4,6,7},{2,6,7},{2,7,3},{0,3,7},{0,7,4},{0,1,2},{0,2,3} };
        //int myels[6][5] = {{0,1,5,4,8},{6,5,1,2,8},{2,3,7,6,8},{7,4,0,3,8},{0,1,2,3,8},{4,5,6,7,8}}; //Sequencia trocada soh para funcionar o AddHDivPyramidRestraints
        for (int iel = 0; iel < 2; iel++) {
            for (int i = 0; i < 5; i++) {
                topolPyr[i] = myelsp[iel][i];
            }
            gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index,0);
        }
        for (int iel = 0; iel < 2; iel++) {
            for (int i = 0; i < 4; i++) {
                topolTet[i] = myelst[iel][i];
            }
            gmesh->CreateGeoElement(ETetraedro, topolTet, matid, index,0);
        }
        
        const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
        
        for (int64_t iel = 0; iel < 12; iel++) {
            for (int i = 0; i < 3; i++) {
                topolTri[i] = triangles[iel][i];
            }
            gmesh->CreateGeoElement(ETriangle, topolTri, bc0, index,0);
        }
        
        gmesh->BuildConnectivity();
        
        std::ofstream out("../HexaPyrTetGmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        
        return gmesh;
    }


#endif

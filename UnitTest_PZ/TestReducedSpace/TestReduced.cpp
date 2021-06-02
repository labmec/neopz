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
#include "TPZGenGrid2D.h"
#include "tpzautopointer.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZBndCond.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "tpzpermutation.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "tpzintpoints.h"
#include "pztrnsform.h"
#include "pzintel.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"

#include "TPZExtendGridDimension.h"

#include "TPZLinearAnalysis.h"

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

#include "TPZAnalyticSolution.h"

using namespace pzshape;


#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testhdiv");
#endif

#include<catch2/catch.hpp>

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

static TPZAutoPointer<TPZGeoMesh> /*TPZGeoMesh * */ CreateOneCuboWithTetraedrons(int nref);
static TPZAutoPointer<TPZGeoMesh> TetrahedralMeshCubo(int64_t nelem,int MaterialId);

static TPZAutoPointer<TPZGeoMesh> GenerateMesh(MElementType eltype, int nelem, int ndiv);
static TPZCompMesh *GenerateH1Mesh(TPZAutoPointer<TPZGeoMesh> gmesh);

static TPZCompMesh *GenerateReducedMesh(TPZCompMesh *cmesh_orig);

void LoadSolutions(TPZCompMesh *cmesh);

static void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);

static TPZFMatrix<STATE> RunConfig(TPZCompMesh *cmesh, int sol);

std::vector<MElementType> allmeshes = {EQuadrilateral,ETriangle,ETetraedro,EPrisma,ECube};

// Tests for the 'voidflux' class.
TEST_CASE("reduced_space","[reduced_space]")
{
    int nelem = 2;
    int ndiv = 0;
    int ntypes = allmeshes.size();
    for (int meshtype = 0; meshtype<ntypes; meshtype++)
    {
        auto gmesh = GenerateMesh(allmeshes[meshtype], nelem, ndiv);
        auto cmeshH1 = GenerateH1Mesh(gmesh);
        LoadSolutions(cmeshH1);
        auto cmeshreduced = GenerateReducedMesh(cmeshH1);
        int nsol = cmeshreduced->Solution().Rows();

        auto oldPrecision = Catch::StringMaker<REAL>::precision;
        for (int sol = 0; sol < nsol; sol++) {
            auto result = RunConfig(cmeshreduced, sol);
            TPZVec<REAL> correct(result.Rows(),0.);
            correct[sol] = 1.;
            for (int i = 0; i< result.Rows(); i++) {
                CAPTURE(i);
                REQUIRE_THAT(result(i), Catch::Matchers::WithinAbs(correct[i],1.e-8));
            }
        }
        Catch::StringMaker<REAL>::precision = oldPrecision;
    }

}

static TPZAutoPointer<TPZGeoMesh> GenerateMesh(MElementType eltype, int nelem, int ndiv)
{
    int dimmodel = 2;
    TPZManVector<int,3> nx(2,nelem);
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = -1.;
    TPZGenGrid2D grid(nx,x0,x1);
    if (eltype == ETriangle|| eltype == EPrisma ) {
        grid.SetElementType(MMeshType::ETriangular);
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
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
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
    return gmesh;
}



template<class TMat = TPZMatPoisson<STATE>>
static void InsertMaterials(TPZCompMesh *cmesh)
{
    int dimmodel = cmesh->Dimension();
    TMat *matpoisP = new TMat(1, dimmodel);
    {
        TPZMaterial *poisP(matpoisP);
        cmesh->InsertMaterialObject(matpoisP);
        cmesh->SetAllCreateFunctionsContinuous();
        TPZFNMatrix<4,STATE> val1(1,1,0.);
        TPZManVector<STATE,1> val2 = {0};
        TPZBndCond *bndP = matpoisP->CreateBC(poisP, -1, 0, val1, val2);
        cmesh->InsertMaterialObject(bndP);
        bndP = matpoisP->CreateBC(poisP, -2, 1, val1, val2);
        cmesh->InsertMaterialObject(bndP);
        bndP = matpoisP->CreateBC(poisP, -3, 1, val1, val2);
        cmesh->InsertMaterialObject(bndP);
    }

}

static TPZCompMesh *GenerateH1Mesh(TPZAutoPointer<TPZGeoMesh> gmesh)
{
    int dimmodel = gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    InsertMaterials<TPZMatPoisson<STATE>>(cmesh);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    // reorder the equations
    TPZLinearAnalysis an(cmesh);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
#endif
    return cmesh;
}

static TPZCompMesh *GenerateReducedMesh(TPZCompMesh *cmesh_orig)
{
    TPZConnect &c0 = cmesh_orig->ConnectVec()[0];
    int nstate = c0.NState();
    int dimmodel = cmesh_orig->Dimension();
    int64_t nshape = cmesh_orig->Solution().Cols();
    TPZCompMesh *cmesh = new TPZCompMesh(cmesh_orig->Reference());
    InsertMaterials<TPZMatPoisson<STATE>>(cmesh);
    TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmesh);
    cmesh->AllocateNewConnect(nshape, nstate, 1);
    auto nc = cmesh->NConnects();
    if(nc != 1) DebugStop();
    TPZConnect &c = cmesh->ConnectVec()[0];
    c.SetNShape(nshape);
    c.SetNState(nstate);
    cmesh->AutoBuild();
    cmesh->Reference()->ResetReference();
    cmesh_orig->LoadReferences();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZReducedSpace *reduce = dynamic_cast<TPZReducedSpace *>(cel);
        if(!reduce) DebugStop();
        TPZGeoEl *gel = reduce->Reference();
        reduce->SetReferredElement(gel->Reference());
    }
    return cmesh;
}

void SetExactSolution(TLaplaceExample1 &config, TPZCompMesh *cmesh);

TPZFMatrix<STATE> ComputeSolution(TLaplaceExample1 &config, TPZCompMesh *cmesh)
{
    SetExactSolution(config, cmesh);
    // avoid renumbering the equations
    TPZLinearAnalysis an(cmesh,false);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
    return an.Solution();
}

typedef TLaplaceExample1 TL;

std::vector<TL::EExactSol> allexact = {TL::EConst, TL::EX, TL::ESinSin,TL::ESinDist, TL::E10SinSin,TL::ECosCos};

void LoadSolutions(TPZCompMesh *cmesh)
{
    int64_t neq = cmesh->Solution().Rows();
    TPZFMatrix<STATE> allsol(neq,allexact.size());
    for (int iex=0; iex<allexact.size(); iex++) {
        TL config;
        config.fExact = allexact[iex];
        ComputeSolution(config, cmesh);
        TPZFMatrix<STATE> &sol = cmesh->Solution();
        for (int ieq=0; ieq<neq; ieq++) {
            allsol(ieq,iex) = sol(ieq,0);
        }
    }
    TPZFMatrix<STATE> &sol = cmesh->Solution();
    sol = allsol;
}

void SetExactSolution(TLaplaceExample1 &config, TPZCompMesh *cmesh)
{
    int matids[] = {1,-1,-2,-3};
    {
        auto *mat =
            dynamic_cast<TPZMaterialT<STATE> *>(cmesh->FindMaterial(1));
        if(!mat) DebugStop();
        
        auto forcingFunction = [&config](const TPZVec<REAL>&x, TPZVec<STATE>&u){
            config.ForcingFunction()->Execute(x, u);
        };
        const auto pOrderForcingFunction = config.ForcingFunction()->PolynomialOrder();
        mat->SetForcingFunction(forcingFunction,pOrderForcingFunction);
    }
    for(int i=1; i<4; i++)
    {
        auto *mat =
            dynamic_cast<TPZBndCondT<STATE> *>(cmesh->FindMaterial(matids[i]));
        if(!mat) DebugStop();
        auto exact = [&config](const TPZVec<REAL>&x, TPZVec<STATE>&u,
                               TPZFMatrix<STATE>&du){
            config.Exact()->Execute(x, u, du);
        };
        mat->SetForcingFunctionBC(exact);
    }
}

static TPZFMatrix<STATE> RunConfig(TPZCompMesh *cmesh, int sol)
{
    TLaplaceExample1 config;
    config.fExact = allexact[sol];
    SetExactSolution(config, cmesh);
    TPZLinearAnalysis an(cmesh);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
    return an.Solution();
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
#ifdef PZ_LOG
                if(logger.isDebugEnabled())
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

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
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZBndCondT.h"
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
#include "pzstepsolver.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
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


using namespace pzshape;


#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testhdivconstant");
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
static TPZAutoPointer<TPZGeoMesh> GenerateMesh(MElementType type, int nelem = 3, int ndiv = 0);
static TPZAutoPointer<TPZGeoMesh> TetrahedralMeshCubo(int64_t nelem,int MaterialId);


static void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);

template<class TSHAPE>
static void IntegralNormal();

// Tests for the 'voidflux' class.
TEST_CASE("integral_normal","[hdivconstant_tests]")
{
    std::cout << "Initializing vector_direction check\n";
//    IntegralNormal<pzshape::TPZShapePrism>();
    IntegralNormal<pzshape::TPZShapeCube>();
    IntegralNormal<pzshape::TPZShapeTetra>();
    IntegralNormal<pzshape::TPZShapeTriang>();
    IntegralNormal<pzshape::TPZShapeQuad>();
//    IntegralNormal<pzshape::TPZShapePrism>();
    std::cout << "Leaving integral_normal check\n";
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
    

    
    return gmesh;


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

#include "TPZShapeHDivConstant.h"

template<class TSHAPE>
void IntegralNormal()
{
    auto gmesh = GenerateMesh(TSHAPE::Type(),1,0);
    int64_t nel = gmesh->NElements();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Type() != TSHAPE::Type()) continue;
        TPZManVector<int64_t> ids(gel->NCornerNodes());
        int nnodes = gel->NCornerNodes();
        for(int no=0; no<nnodes; no++) ids[no] = gel->Node(no).Id();
        TPZShapeHDivConstant<TSHAPE> shape;
        TPZShapeData data;
        TPZManVector<int> orders(TSHAPE::NFacets+1,1);
        TPZManVector<int> sideorient(TSHAPE::NFacets,1);
        shape.Initialize(ids, orders, sideorient, data);
        data.fSideTransformationId.Resize(TSHAPE::NFacets, 0);
        data.fSideOrient.Resize(TSHAPE::NFacets, 1);
        int nshape = shape.NHDivShapeF(data);
        TPZManVector<REAL> pt(3,0.);
        TPZFNMatrix<8,REAL> phi(TSHAPE::Dimension,nshape),divphi(nshape,1);
        shape.Shape(pt,data,phi,divphi);
        // integrate the divergence. It should be one
        divphi *= TSHAPE::RefElVolume();
        phi.Print("phi  = ",std::cout);
        divphi.Print("divphi  = ",std::cout);
    }
}

template void IntegralNormal<pzshape::TPZShapeQuad>();

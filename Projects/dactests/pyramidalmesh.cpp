//
//  piramide.cpp
//  PZ
//
//  Created by Douglas Castro on 1/27/15.
//
//


#include "pyramidalmesh.h"
#include "pzreal.h"
#include "pzvec.h"

const int matId = 1;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;
const int bc4 = -5;
const int bc5 = -6;


static int piramide_2[6][5]=
{
    {0,1,2,3,8},
    {0,1,5,4,8},
    {1,2,6,5,8},
    {3,2,6,7,8},
    {0,3,7,4,8},
    {4,5,6,7,8}
};

PyramidalMesh::PyramidalMesh(TPZGeoMesh *gmesh, int ndiv)
{
    gmesh = new TPZGeoMesh;
    REAL dndiv = ndiv;
    int nref = (int) pow(2., dndiv);
    gmesh = CreateGMeshCubeWithPyramids( nref, matId);
    
}

PyramidalMesh::~PyramidalMesh()
{
    
}

bool PyramidalMesh::DoubleComparer(REAL a, REAL b)
{
    if (IsZero(a-b)){
        return true;
    }
    else{
        return false;
    }
}

void PyramidalMesh::GenerateNodesforPyramidalMesh(TPZGeoMesh *gmesh, int64_t nelem)
{
    int64_t sizenodevec = (nelem+1)*(nelem+1)*(nelem+1)+(nelem*nelem*nelem);
    gmesh->NodeVec().Resize(sizenodevec);
    int posicao = 0;
    for (int64_t i=0; i<=nelem; i++) {
        for (int64_t j=0; j<=nelem; j++) {
            for (int64_t k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                posicao = i*(nelem+1)*(nelem+1)+j*(nelem+1)+k;
                gmesh->NodeVec()[posicao].Initialize(x, *gmesh);
                
            }
        }
    }
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<REAL,3> xc(3);
                xc[0] = k*1./nelem + (0.5)*1./nelem;
                xc[1] = j*1./nelem + (0.5)*1./nelem;
                xc[2] = i*1./nelem + (0.5)*1./nelem;
                posicao++;
                gmesh->NodeVec()[posicao].Initialize(xc, *gmesh);
            }
        }
    }
}

TPZGeoMesh *PyramidalMesh::CreateGMeshCubeWithPyramids(int64_t nelem, int MaterialId)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodesforPyramidalMesh(gmesh,nelem);
    
    int posicao = 0;
    for (int64_t i=0; i<=nelem; i++) {
        for (int64_t j=0; j<=nelem; j++) {
            for (int64_t k=0; k<=nelem; k++) {
                posicao = i*(nelem+1)*(nelem+1)+j*(nelem+1)+k;
            }
        }
    }
    
    
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,9> nodes(9,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[8] = posicao + (k)*(nelem)*(nelem)+(j)*(nelem)+i+1;

                for (int el=0; el<6; el++)
                {
                    TPZManVector<int64_t,5> elnodes(5);
                    int64_t index;
                    for (int il=0; il<5; il++) {
                        elnodes[il] = nodes[piramide_2[el][il]];
                    }
                    gmesh->CreateGeoElement(EPiramide, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    
    // Boundary Conditions
    const int numelements = gmesh->NElements();
    //    const int bczMinus = -3, bczplus = -2, bcids = -1;
    //    const int bczMinus = -1, bczplus = -1, bcids = -1;
    
    for(int el=0; el<numelements; el++)
    {
        TPZManVector <TPZGeoNode,5> Nodefinder(5);
        TPZManVector <REAL,3> nodecoord(3);
        TPZGeoEl *piramide = gmesh->ElementVec()[el];
        TPZVec<int64_t> ncoordVec(0); int64_t sizeOfVec = 0;
        
        // na face z = 0
        for (int i = 0; i < 5; i++)
        {
            int64_t pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (DoubleComparer(nodecoord[2],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc0);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 0
        for (int i = 0; i < 5; i++)
        {
            int64_t pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (DoubleComparer(nodecoord[1],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc1);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 1
        for (int i = 0; i < 5; i++)
        {
            int64_t pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (DoubleComparer(nodecoord[0],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc2);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 1
        for (int i = 0; i < 5; i++)
        {
            int64_t pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (DoubleComparer(nodecoord[1],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc3);
        }
        
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 0
        for (int i = 0; i < 5; i++)
        {
            int64_t pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (DoubleComparer(nodecoord[0],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc4);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face z = 1
        for (int i = 0; i < 5; i++)
        {
            int64_t pos = piramide->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (DoubleComparer(nodecoord[2],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 4)
        {
            int lado = piramide->WhichSide(ncoordVec);
            TPZGeoElSide piramideSide(piramide, lado);
            TPZGeoElBC(piramideSide,bc5);
        }
        
        
        
    }
    
    
    std::ofstream out("CubeWithPyramidswithBcs.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}
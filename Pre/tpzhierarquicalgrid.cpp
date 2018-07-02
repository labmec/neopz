/*
 <one line to give the program's name and a brief idea of what it does.>
 Copyright (C) 2014  <copyright holder> <email>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "tpzhierarquicalgrid.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeotetrahedra.h"

TPZHierarquicalGrid::TPZHierarquicalGrid()
{
    
    fFileName   = "untitled";
    fNonAffineQ = false;
    fIsQuad = true;
    fIsPrism = false;
    fIsTetrahedron = false;
    fComputedGeomesh = NULL;
    ffrontMatID = -999;
    fbackMatID = -1000;
    fSubBases.Resize(0);
    fBase = NULL;
}

TPZHierarquicalGrid::TPZHierarquicalGrid(TPZGeoMesh *Geomesh)
{
    if(Geomesh->NElements() == 0)
    {
        if(fBase) std::cout << "Number of elements" << fBase->NElements() << std::endl;
        DebugStop();
    }
    
    if(Geomesh->NNodes() == 0)
    {
        std::cout << "Number of nodes" << fBase->NElements() << std::endl;
        DebugStop();
    }
    
    fFileName   = "untitled";
    fNonAffineQ = false;
    fIsQuad = true;
    fIsPrism = false;
    fIsTetrahedron = false;
    fComputedGeomesh = NULL;
    ffrontMatID = -999;
    fbackMatID = -1000;
    fBase       = Geomesh;
    fSubBases.Resize(0);
}



TPZHierarquicalGrid::TPZHierarquicalGrid(const TPZHierarquicalGrid& other)
{
    
    fFileName   = other.fFileName;
    fNonAffineQ = other.fNonAffineQ;
    fIsQuad = other.fIsQuad;
    fIsPrism = other.fIsPrism;
    fIsTetrahedron = other.fIsTetrahedron;
    fComputedGeomesh = other.fComputedGeomesh;
    ffrontMatID = other.ffrontMatID;
    fbackMatID = other.fbackMatID;
    fSubBases = other.fSubBases;
    fBase = other.fBase;
    
}

TPZHierarquicalGrid::~TPZHierarquicalGrid()
{
    
}

TPZHierarquicalGrid& TPZHierarquicalGrid::operator=(const TPZHierarquicalGrid& other)
{
    fFileName   = other.fFileName;
    fNonAffineQ = other.fNonAffineQ;
    fIsQuad = other.fIsQuad;
    fIsPrism = other.fIsPrism;
    fIsTetrahedron = other.fIsTetrahedron;
    fComputedGeomesh = other.fComputedGeomesh;
    ffrontMatID = other.ffrontMatID;
    fbackMatID = other.fbackMatID;
    fSubBases = other.fSubBases;
    fBase = other.fBase;
    return *this;
}

bool TPZHierarquicalGrid::operator==(const TPZHierarquicalGrid& other) const
{
    DebugStop();
    ///TODO: return ...;
    return 0; //to fix WIN compiler error
}

TPZGeoMesh * TPZHierarquicalGrid::ComputeExtrusion(REAL t, REAL dt, int n)
{
    fSubBases.Resize(n+1);
    fComputedGeomesh = new TPZGeoMesh;
    int NNodesBase = fBase->NNodes();
    fComputedGeomesh->NodeVec().Resize( (n+1) * NNodesBase);
    
    // Creating new elements
    int nodeId = 0;
    REAL sing = 1.0;
    for(int il = 0; il < (n+1); il++ )
    {
        // copying l extrusions
        // For a while all of them are equal
        //        fSubBases[il] = new TPZGeoMesh(fBase);
        
        for(int inode = 0; inode < NNodesBase; inode++)
        {
            TPZVec<REAL> tpara(1),NewCoordinates_r(3,0.0);
            TPZVec<REAL> NewCoordinates(3,0.0);
            tpara[0] = dt*il+t;
            
            TPZVec<REAL> Coordinates(3,0.0);
            fComputedGeomesh->NodeVec()[inode + il * NNodesBase] =  fBase->NodeVec()[inode];
            //            NewGeomesh->NodeVec()[inode + il*fBase->NNodes()] =  fSubBases[il]->NodeVec()[inode];
            fComputedGeomesh->NodeVec()[inode + il * NNodesBase].GetCoordinates(Coordinates);
            
            fParametricFunction->Execute(tpara,NewCoordinates);
            NewCoordinates_r[0] = REAL(NewCoordinates[0]);
            NewCoordinates_r[1] = REAL(NewCoordinates[1]);
            NewCoordinates_r[2] = REAL(NewCoordinates[2]);
            
            if(((il+1)%2==0 && fNonAffineQ) && fBase->Dimension() == 2){
                Coordinates[0]+=NewCoordinates_r[0];
                Coordinates[1]+=NewCoordinates_r[1];
                Coordinates[2]+=NewCoordinates_r[2];
                Coordinates[2]+= sing*dt/2.0;
                sing *= -1.0;
            }
            else{
                Coordinates[0]+=NewCoordinates_r[0];
                Coordinates[1]+=NewCoordinates_r[1];
                Coordinates[2]+=NewCoordinates_r[2];
            }
            
            fComputedGeomesh->NodeVec()[inode + il * NNodesBase].SetCoord(Coordinates);
            fComputedGeomesh->NodeVec()[inode + il * NNodesBase].SetNodeId(nodeId);
            nodeId++;
        }
    }
    
    fComputedGeomesh->SetMaxNodeId(nodeId);
    
//    int fbasedim = fBase->Dimension();
    
    int elid=0;
    for(int iel = 0; iel < fBase->NElements(); iel++)
    {
        int ielDim      = fBase->ElementVec()[iel]->Dimension();
        int ielMatId    = fBase->ElementVec()[iel]->MaterialId();
        CreateGeometricElement(n,iel,ielDim,ielMatId,elid);
    }
    
    fComputedGeomesh->SetMaxElementId(elid);
    fComputedGeomesh->BuildConnectivity();
    return fComputedGeomesh;
    
}

void TPZHierarquicalGrid::CreateGeometricElement(int n, int iel,int eldim, int elmatid, int &elid)
{
    int jump =  fBase->NNodes();
    int dim = fBase->Dimension();
    
    TPZGeoEl *gel =  fBase->ElementVec()[iel];
    int gelNodes = gel->NNodes();

    // Computing  current topology
    TPZManVector<int64_t,10> CTopology(gelNodes);
    for(int inode = 0; inode < CTopology.size(); inode++)
    {
        CTopology[inode] = gel->Node(inode).Id();
    }
    
    bool Is1D = false;
    bool Is2D = false;
    bool Is3D = false;
    
    for(int il = 1; il < (n+1); il++ )
    {
        
        switch (eldim) {
            case 0:
            {
                
                //Defining boundaries
                if (dim==0) {
                    if (il==1){
                        TPZVec<int64_t> Topology(gelNodes);
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elid++, Topology, ffrontMatID,*fComputedGeomesh);
                    }
                    if (il==n)
                    {
                        TPZVec<int64_t> Topology(gelNodes);
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elid++, Topology, fbackMatID,*fComputedGeomesh);
                    }
                }
                
                //                    if (dim==1) {
                //                        if (il==1){
                //                            TPZVec<int64_t> Topology(gelNodes+1);
                //                            Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                //                            Topology[1]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                //                            new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elid++, Topology, elmatid,*fComputedGeomesh);
                //                        }
                //                        if (il==n)
                //                        {
                //                            TPZVec<int64_t> Topology(gelNodes+1);
                //                            Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                //                            Topology[1]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                //                            new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elid++, Topology, elmatid,*fComputedGeomesh);
                //
                //                        }
                //                    }
                
                TPZVec<int64_t> Topology(gelNodes+1);
                Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                Topology[1]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elid++, Topology, elmatid,*fComputedGeomesh);
                Is1D = true;
            }
                break;
            case 1:
            {
                if (dim==1) {
                    if (il==1){
                        TPZVec<int64_t> Topology(gelNodes+1);
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elid++, Topology, ffrontMatID,*fComputedGeomesh);
                    }
                    if (il==n)
                    {
                        TPZVec<int64_t> Topology(gelNodes+1);
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elid++, Topology, fbackMatID,*fComputedGeomesh);
                    }
                }
                
                if (fIsQuad) {
                    // quadrilateras
                    TPZVec<int64_t> Topology(gelNodes+2);
                    Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                    Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                    Topology[2]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                    Topology[3]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                    new TPZGeoElRefPattern < pzgeom::TPZGeoQuad > (elid++, Topology, elmatid,*fComputedGeomesh);
                }
                else
                {
                    // triangles
                    TPZVec<int64_t> Topology(gelNodes+1);
                    Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                    Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                    Topology[2]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle > (elid++, Topology, elmatid,*fComputedGeomesh);
                    
                    Topology[0]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                    Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                    Topology[2]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle > (elid++, Topology, elmatid,*fComputedGeomesh);
                    
                }
                Is2D = true;
                
            }
                break;
            case 2:
            {
                
                if (dim==2) {
                    if (il==1){
                        
                        if (gel->Type() == EQuadrilateral) {
                            // quadrilateras
                            TPZVec<int64_t> Topology(gelNodes+2);
                            Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                            Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                            Topology[2]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                            Topology[3]=fComputedGeomesh->NodeVec()[CTopology[3] + (il - 1) * jump].Id();
                            new TPZGeoElRefPattern < pzgeom::TPZGeoQuad > (elid++, Topology, ffrontMatID,*fComputedGeomesh);
                        }
                        else{
                            // triangles
                            TPZVec<int64_t> Topology(gelNodes+1);
                            Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                            Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                            Topology[2]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                            new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle > (elid++, Topology, ffrontMatID,*fComputedGeomesh);
                            
                        }
                    }
                    if (il==n)
                    {
                        if (gel->Type() == EQuadrilateral) {
                            // quadrilateras
                            TPZVec<int64_t> Topology(gelNodes+2);
                            Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                            Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                            Topology[2]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 0) * jump].Id();
                            Topology[3]=fComputedGeomesh->NodeVec()[CTopology[3] + (il - 0) * jump].Id();
                            new TPZGeoElRefPattern < pzgeom::TPZGeoQuad > (elid++, Topology, fbackMatID,*fComputedGeomesh);
                        }
                        else{
                            // triangles
                            TPZVec<int64_t> Topology(gelNodes+1);
                            Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                            Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                            Topology[2]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 0) * jump].Id();
                            new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle > (elid++, Topology, fbackMatID,*fComputedGeomesh);
                            
                        }
                    }
                }
                
                if (!fIsTetrahedron) {

                    if (fIsPrism) {
                        // Prisms
                        TPZVec<int64_t> Topology(2*gelNodes);
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[2]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                        Topology[3]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        Topology[4]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[5]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoPrism > (elid++, Topology, elmatid,*fComputedGeomesh);
                    }
                    else
                    {
                        // Cubes
                        TPZVec<int64_t> Topology(gelNodes+4);
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[2]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                        Topology[3]=fComputedGeomesh->NodeVec()[CTopology[3] + (il - 1) * jump].Id();
                        Topology[4]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        Topology[5]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[6]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 0) * jump].Id();
                        Topology[7]=fComputedGeomesh->NodeVec()[CTopology[3] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoCube > (elid++, Topology, elmatid,*fComputedGeomesh);
                    }
                }
                else
                {
                    if (fIsPrism) {
                        // Prisms
                        TPZVec<int64_t> Topology(2*gelNodes);
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[2]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                        Topology[3]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        Topology[4]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[5]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoPrism > (elid++, Topology, elmatid,*fComputedGeomesh);
                    }
                    else{
                        // Tetrahedron
                        TPZVec<int64_t> Topology(gelNodes+1);
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[2]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        Topology[3]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoTetrahedra > (elid++, Topology, elmatid,*fComputedGeomesh);
                        
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[2]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        Topology[3]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoTetrahedra > (elid++, Topology, elmatid,*fComputedGeomesh);
                        
                        Topology[0]=fComputedGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        Topology[1]=fComputedGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[2]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 0) * jump].Id();
                        Topology[3]=fComputedGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoTetrahedra > (elid++, Topology, elmatid,*fComputedGeomesh);
                        
                    }
                    
                }
                Is3D = true;
                
            }
                break;
            default:
            {
                std::cout << "Connection not implemented " << std::endl;
                DebugStop();
            }
                break;
        }
    }
    
    if (Is1D) {
        fComputedGeomesh->SetDimension(1);
    }
    
    if (Is2D) {
        fComputedGeomesh->SetDimension(2);
    }
    
    if (Is3D) {
        fComputedGeomesh->SetDimension(3);
    }
    
}

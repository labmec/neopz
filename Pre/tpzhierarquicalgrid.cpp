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
#include "pzgeotetrahedra.h"

TPZHierarquicalGrid::TPZHierarquicalGrid()
{
    std::string Name = "untitled";
    fFileName   = Name;
    fSubBases.Resize(0);
}

TPZHierarquicalGrid::TPZHierarquicalGrid(TPZGeoMesh *Geomesh)
{
    if(Geomesh->NElements() == 0)
    {
        std::cout << "Number of elements" << fBase->NElements() << std::endl;
        DebugStop();
    }
    
    if(Geomesh->NNodes() == 0)
    {
        std::cout << "Number of nodes" << fBase->NElements() << std::endl;
        DebugStop();
    }    
    
    std::string Name = "untitled";
    fFileName   = Name;
    fBase       = Geomesh;
    fSubBases.Resize(0);    
}



TPZHierarquicalGrid::TPZHierarquicalGrid(const TPZHierarquicalGrid& other)
{

}

TPZHierarquicalGrid::~TPZHierarquicalGrid()
{

}

TPZHierarquicalGrid& TPZHierarquicalGrid::operator=(const TPZHierarquicalGrid& other)
{
return *this;
}

bool TPZHierarquicalGrid::operator==(const TPZHierarquicalGrid& other) const
{
    DebugStop();
///TODO: return ...;
}

TPZGeoMesh * TPZHierarquicalGrid::ComputeExtrusion(STATE t, STATE dt, int n)
{
    fSubBases.Resize(n+1);
    TPZGeoMesh * NewGeomesh = new TPZGeoMesh;
    NewGeomesh->NodeVec().Resize( (n+1) * fBase->NNodes());    
    
    // Creating new elements
    int nodeId = 0;
    
    for(int il = 0; il < (n+1); il++ )
    {
        // copying l extrusions
        // For a while all of them are equal
        fSubBases[il] = new TPZGeoMesh(fBase);
        
        for(int inode = 0; inode < fBase->NNodes(); inode++)
        {
            TPZVec<REAL> tpara(1),NewCoordinates(3,0.0);
            tpara[0] = dt*il+t;
            
            TPZVec<REAL> Coordinates(3,0.0);
            NewGeomesh->NodeVec()[inode + il*fBase->NNodes()] =  fSubBases[il]->NodeVec()[inode];
            NewGeomesh->NodeVec()[inode + il*fBase->NNodes()].GetCoordinates(Coordinates);

            fParametricFunction->Execute(tpara,NewCoordinates);
            Coordinates[0]+=NewCoordinates[0];
            Coordinates[1]+=NewCoordinates[1];
            Coordinates[2]+=NewCoordinates[2];
            
            NewGeomesh->NodeVec()[inode + il*fBase->NNodes()].SetCoord(Coordinates);
            NewGeomesh->NodeVec()[inode + il*fBase->NNodes()].SetNodeId(nodeId++);
        }        
    }
    
    int elid=0;
    for(int iel = 0; iel < fBase->NElements(); iel++)
    {
        int ielDim      = fBase->ElementVec()[iel]->Dimension();
        int ielMatId    = fBase->ElementVec()[iel]->MaterialId();
        CreateGeometricElement(n,iel,ielDim,ielMatId,elid,NewGeomesh);
    }
    
//        TPZVec<long> Topology(1);
//        
//        Topology.Resize(1);
//        Topology[0]=NewGeomesh->NodeVec()[0].Id();
//        new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elInd++, Topology, -1,*NewGeomesh);         
//        Topology[0]=NewGeomesh->NodeVec()[NewGeomesh->NNodes()-1].Id();
//        new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elInd++, Topology, -2,*NewGeomesh);
    
        NewGeomesh->BuildConnectivity();
        return NewGeomesh;
        
}

void TPZHierarquicalGrid::CreateGeometricElement(int n, int iel,int eldim, int elmatid, int &elid, TPZGeoMesh * NewGeomesh)
{
    int jump =  fBase->NNodes();
    
    TPZGeoEl *gel =  fBase->ElementVec()[iel];
    int gelNodes = gel->NNodes();
    
    // Computing  current topology
    TPZVec<long> CTopology(gelNodes);
    for(int inode = 0; inode < CTopology.size(); inode++)
    {
        TPZGeoNode GelNode = gel->Node(inode);
        CTopology[inode] = GelNode.Id();
    }
    
    for(int il = 1; il < (n+1); il++ )
    {

            switch (eldim) {
                case 0:
                {
                    TPZVec<long> Topology(gelNodes+1);
                    Topology[0]=NewGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                    Topology[1]=NewGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elid++, Topology, elmatid,*NewGeomesh);
                }
                    break;
                case 1:
                {
                    if (true) {
                        // quadrilateras
                        TPZVec<long> Topology(gelNodes+2);
                        Topology[0]=NewGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=NewGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[2]=NewGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[3]=NewGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoQuad > (elid++, Topology, elmatid,*NewGeomesh);
                    }
                    else
                    {
                        // triangles
                        TPZVec<long> Topology(gelNodes+1);
                        Topology[0]=NewGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=NewGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[2]=NewGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle > (elid++, Topology, elmatid,*NewGeomesh);
                        
                        Topology[0]=NewGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[1]=NewGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[2]=NewGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle > (elid++, Topology, elmatid,*NewGeomesh);
                        
                    }

                }
                    break;
                case 2:
                {
                    if (true) {
                        // Cubes
                        TPZVec<long> Topology(gelNodes+4);
                        Topology[0]=NewGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=NewGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[2]=NewGeomesh->NodeVec()[CTopology[2] + (il - 1) * jump].Id();
                        Topology[3]=NewGeomesh->NodeVec()[CTopology[3] + (il - 1) * jump].Id();
                        Topology[4]=NewGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        Topology[5]=NewGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[6]=NewGeomesh->NodeVec()[CTopology[2] + (il - 0) * jump].Id();
                        Topology[7]=NewGeomesh->NodeVec()[CTopology[3] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoCube > (elid++, Topology, elmatid,*NewGeomesh);
                    }
                    else
                    {
                        // Tetrahedros
                        TPZVec<long> Topology(gelNodes+1);
                        Topology[0]=NewGeomesh->NodeVec()[CTopology[0] + (il - 1) * jump].Id();
                        Topology[1]=NewGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[2]=NewGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle > (elid++, Topology, elmatid,*NewGeomesh);
                        
                        Topology[0]=NewGeomesh->NodeVec()[CTopology[1] + (il - 1) * jump].Id();
                        Topology[1]=NewGeomesh->NodeVec()[CTopology[1] + (il - 0) * jump].Id();
                        Topology[2]=NewGeomesh->NodeVec()[CTopology[0] + (il - 0) * jump].Id();
                        new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle > (elid++, Topology, elmatid,*NewGeomesh);
                        
                    }
                    
                }
                    break;
                default:
                {
                    std::cout << "connection not implemented " << std::endl;
                    DebugStop();
                }
                    break;
            }
    }
    

}

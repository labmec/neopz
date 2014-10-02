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
    int elInd = 0;
    int matId = 1;    
    
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
    
    // chose how to make the connections
    
        TPZVec<long> Topology(2);
        
        for(int inode = 1; inode < NewGeomesh->NNodes(); inode++)
        {
            Topology[0]=NewGeomesh->NodeVec()[inode-1].Id();
            Topology[1]=NewGeomesh->NodeVec()[inode].Id();
            new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elInd++, Topology, matId,*NewGeomesh); 
        }        
        
        Topology.Resize(1);
        Topology[0]=NewGeomesh->NodeVec()[0].Id();
        new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elInd++, Topology, -1,*NewGeomesh);         
        Topology[0]=NewGeomesh->NodeVec()[NewGeomesh->NNodes()-1].Id();
        new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elInd++, Topology, -2,*NewGeomesh);
        
        NewGeomesh->BuildConnectivity();
        return NewGeomesh;
        
}

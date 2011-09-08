/**
 * @file
 * @brief Contains the implementation of the TPZChangeEl methods. 
 */
#include "tpzchangeel.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>

#include "tpzmathtools.h"

#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticpyramid.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticcube.h"

#include "TPZGeoElement.h.h"
#include "pzgeoelside.h"
#include "pzstack.h"
#include "tpzgeoelrefpattern.h"

#include <sstream>
using namespace std;
using namespace pzgeom;

TPZChangeEl::TPZChangeEl()
{
}

TPZChangeEl::~TPZChangeEl()
{
}

TPZGeoEl * TPZChangeEl::ChangeToQuadratic(TPZGeoMesh *Mesh, int ElemIndex)
{
	TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
    
    MElementType elType = OldElem->Type();
    int oldId = OldElem->Id();
    int oldIndex = OldElem->Index();
    int oldMatId = OldElem->MaterialId();
    int nsides = OldElem->NSides();
    
    TPZVec<TPZGeoElSide> oldNeigh(nsides);
    for(int s = 0; s < nsides; s++)
    {   
        TPZGeoElSide mySide(OldElem, s);
        oldNeigh[s] = mySide.Neighbour();
    }
    
	TPZGeoEl * NewElem = NULL;
    
    switch(elType) 
    {
        case(EOned) :
        {
            /** Creating Midnodes */
            TPZGeoNode Node2;
            TPZVec<REAL> Coord2(3);
            
            /** Setting Midnodes Coordinates */
            TPZGeoNode Node0 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(0)];
            TPZGeoNode Node1 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(1)];
            for(int i = 0; i < 3; i++)
            {
                Coord2[i] = (Node0.Coord(i) + Node1.Coord(i)) / 2.;
            }
            Node2.SetCoord(Coord2);
            
            /** Setting Midnodes Id's */
            int NewNodeId = Mesh->NNodes();
            Mesh->SetNodeIdUsed(NewNodeId);
            Node2.SetNodeId(NewNodeId);
            
            /** Allocating Memory for MidNodes and Pushing Them */
            TPZVec <int> NodesSequence(3) ;
            NodesSequence[0]  = OldElem->NodeIndex(0);
            NodesSequence[1]  = OldElem->NodeIndex(1);
            NodesSequence[2]  = Mesh->NodeVec().AllocateNewElement();
            Mesh->NodeVec()[NodesSequence[2]] = Node2;               
            
            /** Deleting OldElem */
            Mesh->DeleteElement(OldElem);
            
            /** Inserting New Element in Mesh */            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticLine>(oldId,NodesSequence,oldMatId,*Mesh);
            for(int s = 0; s < nsides; s++)
            {   
                TPZGeoElSide neigh = oldNeigh[s];
                NewElem->SetNeighbour(s, neigh);
            }
            
            break;
        }  
        case(ETriangle) :
        {
            int oldNNodes = 3;
            int newNNodes = 6;
            
            TPZVec<int> NodesSequence(newNNodes);
            
            /** Creating Midnodes */
            TPZVec<TPZGeoNode> node(newNNodes);
            for(int n = 0; n < oldNNodes; n++)
            {
                node[n] = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(n)];
                NodesSequence[n]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(n);
            }
            TPZVec<REAL> Coord(3);
            
            /** Setting Midnodes Coordinates */
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[1].Coord(i)) / 2.;
            node[3].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[2].Coord(i)) / 2.;
            node[4].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[0].Coord(i)) / 2.;
            node[5].SetCoord(Coord);
            
            /** Setting Midnodes Id's */
            int NewNodeId = Mesh->NNodes();
            Mesh->SetNodeIdUsed(NewNodeId + (newNNodes - oldNNodes - 1));
            int nCount = 0;
            for(int n = oldNNodes; n < newNNodes; n++)
            {
                node[n].SetNodeId(NewNodeId + nCount);
                
                /** Allocating Memory for MidNodes and Pushing Them */
                NodesSequence[n]  = Mesh->NodeVec().AllocateNewElement();
                Mesh->NodeVec()[NodesSequence[n]] = node[n];
                
                nCount++;
            }
            
            /** Deleting OldElem */
            Mesh->DeleteElement(OldElem);
            
            /** Inserting New Element in Mesh */            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticTrig>(oldId,NodesSequence,oldMatId,*Mesh);
            for(int s = 0; s < nsides; s++)
            {   
                TPZGeoElSide neigh = oldNeigh[s];
                NewElem->SetNeighbour(s, neigh);
            }
            
            break;
        }
        case(EQuadrilateral) :
        {
            int oldNNodes = 4;
            int newNNodes = 8;
            
            TPZVec<int> NodesSequence(newNNodes);
            
            /** Creating Midnodes */
            TPZVec<TPZGeoNode> node(newNNodes);
            for(int n = 0; n < oldNNodes; n++)
            {
                node[n] = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(n)];
                NodesSequence[n]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(n);
            }
            TPZVec<REAL> Coord(3);
            
            /** Setting Midnodes Coordinates */
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[1].Coord(i)) / 2.;
            node[4].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[2].Coord(i)) / 2.;
            node[5].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[3].Coord(i)) / 2.;
            node[6].SetCoord(Coord);

            for(int i = 0; i < 3; i++) Coord[i] = (node[3].Coord(i) + node[0].Coord(i)) / 2.;
            node[7].SetCoord(Coord);
            
            /** Setting Midnodes Id's */
            int NewNodeId = Mesh->NNodes();
            Mesh->SetNodeIdUsed(NewNodeId + (newNNodes - oldNNodes - 1));
            int nCount = 0;
            for(int n = oldNNodes; n < newNNodes; n++)
            {
                node[n].SetNodeId(NewNodeId + nCount);
                
                /** Allocating Memory for MidNodes and Pushing Them */
                NodesSequence[n]  = Mesh->NodeVec().AllocateNewElement();
                Mesh->NodeVec()[NodesSequence[n]] = node[n];
                
                nCount++;
            }
            
            /** Deleting OldElem */
            Mesh->DeleteElement(OldElem);
            
            /** Inserting New Element in Mesh */            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticQuad>(oldId,NodesSequence,oldMatId,*Mesh);
            for(int s = 0; s < nsides; s++)
            {   
                TPZGeoElSide neigh = oldNeigh[s];
                NewElem->SetNeighbour(s, neigh);
            }
            
            break;
        }
        case(ETetraedro) :
        {
            int oldNNodes = 4;
            int newNNodes = 10;
            
            TPZVec<int> NodesSequence(newNNodes);
            
            /** Creating Midnodes */
            TPZVec<TPZGeoNode> node(newNNodes);
            for(int n = 0; n < oldNNodes; n++)
            {
                node[n] = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(n)];
                NodesSequence[n]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(n);
            }
            TPZVec<REAL> Coord(3);
            
            /** Setting Midnodes Coordinates */
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[1].Coord(i)) / 2.;
            node[4].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[2].Coord(i)) / 2.;
            node[5].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[0].Coord(i)) / 2.;
            node[6].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[3].Coord(i)) / 2.;
            node[7].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[3].Coord(i)) / 2.;
            node[8].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[3].Coord(i)) / 2.;
            node[9].SetCoord(Coord);
            
            /** Setting Midnodes Id's */
            int NewNodeId = Mesh->NNodes();
            Mesh->SetNodeIdUsed(NewNodeId + (newNNodes - oldNNodes - 1));
            int nCount = 0;
            for(int n = oldNNodes; n < newNNodes; n++)
            {
                node[n].SetNodeId(NewNodeId + nCount);
                
                /** Allocating Memory for MidNodes and Pushing Them */
                NodesSequence[n]  = Mesh->NodeVec().AllocateNewElement();
                Mesh->NodeVec()[NodesSequence[n]] = node[n];
                
                nCount++;
            }
            
            /** Deleting OldElem */
            Mesh->DeleteElement(OldElem);
            
            /** Inserting New Element in Mesh */            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticTetra>(oldId,NodesSequence,oldMatId,*Mesh);
            for(int s = 0; s < nsides; s++)
            {   
                TPZGeoElSide neigh = oldNeigh[s];
                NewElem->SetNeighbour(s, neigh);
            }
            
            break;
        }
        case(EPiramide) :
        {
            int oldNNodes = 5;
            int newNNodes = 13;
            
            TPZVec<int> NodesSequence(newNNodes);
            
            /** Creating Midnodes */
            TPZVec<TPZGeoNode> node(newNNodes);
            for(int n = 0; n < oldNNodes; n++)
            {
                node[n] = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(n)];
                NodesSequence[n]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(n);
            }
            TPZVec<REAL> Coord(3);
            
            /** Setting Midnodes Coordinates */
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[1].Coord(i)) / 2.;
            node[5].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[2].Coord(i)) / 2.;
            node[6].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[3].Coord(i)) / 2.;
            node[7].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[3].Coord(i) + node[0].Coord(i)) / 2.;
            node[8].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[4].Coord(i)) / 2.;
            node[9].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[4].Coord(i)) / 2.;
            node[10].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[4].Coord(i)) / 2.;
            node[11].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[3].Coord(i) + node[4].Coord(i)) / 2.;
            node[12].SetCoord(Coord);
            
            /** Setting Midnodes Id's */
            int NewNodeId = Mesh->NNodes();
            Mesh->SetNodeIdUsed(NewNodeId + (newNNodes - oldNNodes - 1));
            int nCount = 0;
            for(int n = oldNNodes; n < newNNodes; n++)
            {
                node[n].SetNodeId(NewNodeId + nCount);
                
                /** Allocating Memory for MidNodes and Pushing Them */
                NodesSequence[n]  = Mesh->NodeVec().AllocateNewElement();
                Mesh->NodeVec()[NodesSequence[n]] = node[n];
                
                nCount++;
            }
            
            /** Deleting OldElem */
            Mesh->DeleteElement(OldElem);
            
            /** Inserting New Element in Mesh */            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticPyramid>(oldId,NodesSequence,oldMatId,*Mesh);
            for(int s = 0; s < nsides; s++)
            {   
                TPZGeoElSide neigh = oldNeigh[s];
                NewElem->SetNeighbour(s, neigh);
            }
            
            break;
        }    
        case(EPrisma) :
        {
            int oldNNodes = 6;
            int newNNodes = 15;
            
            TPZVec<int> NodesSequence(newNNodes);
            
            /** Creating Midnodes */
            TPZVec<TPZGeoNode> node(newNNodes);
            for(int n = 0; n < oldNNodes; n++)
            {
                node[n] = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(n)];
                NodesSequence[n]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(n);
            }
            TPZVec<REAL> Coord(3);
            
            /** Setting Midnodes Coordinates */
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[1].Coord(i)) / 2.;
            node[6].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[2].Coord(i)) / 2.;
            node[7].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[0].Coord(i)) / 2.;
            node[8].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[3].Coord(i)) / 2.;
            node[9].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[4].Coord(i)) / 2.;
            node[10].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[5].Coord(i)) / 2.;
            node[11].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[3].Coord(i) + node[4].Coord(i)) / 2.;
            node[12].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[4].Coord(i) + node[5].Coord(i)) / 2.;
            node[13].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[5].Coord(i) + node[3].Coord(i)) / 2.;
            node[14].SetCoord(Coord);
            
            /** Setting Midnodes Id's */
            int NewNodeId = Mesh->NNodes();
            Mesh->SetNodeIdUsed(NewNodeId + (newNNodes - oldNNodes - 1));
            int nCount = 0;
            for(int n = oldNNodes; n < newNNodes; n++)
            {
                node[n].SetNodeId(NewNodeId + nCount);
                
                /** Allocating Memory for MidNodes and Pushing Them */
                NodesSequence[n]  = Mesh->NodeVec().AllocateNewElement();
                Mesh->NodeVec()[NodesSequence[n]] = node[n];
                
                nCount++;
            }
            
            /** Deleting OldElem */
            Mesh->DeleteElement(OldElem);
            
            /** Inserting New Element in Mesh */            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticPrism>(oldId,NodesSequence,oldMatId,*Mesh);
            for(int s = 0; s < nsides; s++)
            {   
                TPZGeoElSide neigh = oldNeigh[s];
                NewElem->SetNeighbour(s, neigh);
            }
            
            break;
        }
        case(ECube) :
        {
            int oldNNodes = 8;
            int newNNodes = 20;
            
            TPZVec<int> NodesSequence(newNNodes);
            
            /** Creating Midnodes */
            TPZVec<TPZGeoNode> node(newNNodes);
            for(int n = 0; n < oldNNodes; n++)
            {
                node[n] = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(n)];
                NodesSequence[n]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(n);
            }
            TPZVec<REAL> Coord(3);
            
            /** Setting Midnodes Coordinates */
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[1].Coord(i)) / 2.;
            node[8].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[2].Coord(i)) / 2.;
            node[9].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[3].Coord(i)) / 2.;
            node[10].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[3].Coord(i) + node[0].Coord(i)) / 2.;
            node[11].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[0].Coord(i) + node[4].Coord(i)) / 2.;
            node[12].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[1].Coord(i) + node[5].Coord(i)) / 2.;
            node[13].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[2].Coord(i) + node[6].Coord(i)) / 2.;
            node[14].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[3].Coord(i) + node[7].Coord(i)) / 2.;
            node[15].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[4].Coord(i) + node[5].Coord(i)) / 2.;
            node[16].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[5].Coord(i) + node[6].Coord(i)) / 2.;
            node[17].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[6].Coord(i) + node[7].Coord(i)) / 2.;
            node[18].SetCoord(Coord);
            
            for(int i = 0; i < 3; i++) Coord[i] = (node[7].Coord(i) + node[4].Coord(i)) / 2.;
            node[19].SetCoord(Coord);
            
            /** Setting Midnodes Id's */
            int NewNodeId = Mesh->NNodes();
            Mesh->SetNodeIdUsed(NewNodeId + (newNNodes - oldNNodes - 1));
            int nCount = 0;
            for(int n = oldNNodes; n < newNNodes; n++)
            {
                node[n].SetNodeId(NewNodeId + nCount);
                
                /** Allocating Memory for MidNodes and Pushing Them */
                NodesSequence[n]  = Mesh->NodeVec().AllocateNewElement();
                Mesh->NodeVec()[NodesSequence[n]] = node[n];
                
                nCount++;
            }
            
            /** Deleting OldElem */
            Mesh->DeleteElement(OldElem);
            
            /** Inserting New Element in Mesh */            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticCube>(oldId,NodesSequence,oldMatId,*Mesh);
            for(int s = 0; s < nsides; s++)
            {   
                TPZGeoElSide neigh = oldNeigh[s];
                NewElem->SetNeighbour(s, neigh);
            }
            
            break;
        }
        default :
        {
            DebugStop();
            break;
        }
    }
    
	return NewElem;
}

TPZGeoEl* TPZChangeEl::QuarterPoints(TPZGeoMesh *Mesh, int ElemIndex, int side)
{
	TPZGeoEl *OldElem = Mesh->ElementVec()[ElemIndex];
    if(!OldElem)
    {
        DebugStop();
    }
	if(OldElem->TypeName() == "Triangle" && OldElem->NNodes() == 3)
	{
		if(side < 0 || side > 6) { cout << "Invalid Side to Compute Quarter Points!\nSee QuarterPoints Method!\n"; exit(-1);}
		Mesh->SetNodeIdUsed(2);
		TPZGeoEl * newEl = ChangeToQuadratic(Mesh,ElemIndex);
		switch(side)
		{
			case 0:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(3)].
					SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
					SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
				}
				break;
			case 1:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(3)].
					SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
				}
				break;
			case 2:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
				}
				break;
			case 3:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
				}
				break;
			case 4:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(3)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
				}
				break;
			case 5:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(3)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
				}
			case 6:
				// Do nothing!
				break;
		}
		return newEl;
	}
	else if(OldElem->TypeName() == "Quad" && OldElem->NNodes() == 4)
	{
		if(side < 0 || side > 8) { cout << "Invalid Side to Compute Quarter Points!\nSee QuarterPoints Method!\n"; exit(-1);}
		Mesh->SetNodeIdUsed(3);
		TPZGeoEl* newEl = ChangeToQuadratic(Mesh,ElemIndex);
		switch(side)
		{
			case 0:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(7)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
				}
				break;
			case 1:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
				}
				break;
			case 2:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(6)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
				}
				break;
			case 3:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(6)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(7)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
				}
				break;
			case 4:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(7)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
				}
				break;
			case 5:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(6)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
				}
				break;
			case 6:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(7)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
				}
				break;
			case 7:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
					
                    Mesh->NodeVec()[newEl->NodeIndex(6)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
				}
			case 8:
				// Do nothing!
				break;
		}
        return newEl;
	}
	else if(OldElem->NNodes() == 2)
	{
		if(side < 0 || side > 2) { cout << "Invalid Side to Compute Quarter Points!\nSee QuarterPoints Method!\n"; exit(-1);}
		TPZGeoEl* newEl = ChangeToQuadratic(Mesh,ElemIndex);
		switch(side)
		{
			case 0:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(2)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
				}
				break;
			case 1:
				for(int i = 0; i < 3; i++)
				{
                    Mesh->NodeVec()[newEl->NodeIndex(2)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
				}
				break;
		}
        return newEl;
	}
	else { cout << "Element type don't recognized!\nSee QuarterPoints Method!\n"; exit(-1); }
}

TPZGeoEl * TPZChangeEl::ChangeToGeoBlend(TPZGeoMesh *Mesh, int ElemIndex)
{
	TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
	if(!OldElem){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " - NULL geometric element.\n";
		return NULL;
	}
	
	const int nnodes = OldElem->NNodes();
	TPZManVector<int> nodeindexes(nnodes);
	for(int i = 0; i < nnodes; i++) nodeindexes[i] = OldElem->NodeIndex(i);
	int newindex;
	TPZGeoEl * NewElem = Mesh->CreateGeoBlendElement(OldElem->Type(),nodeindexes,OldElem->MaterialId(),newindex);
	TPZChangeEl::AdjustNeighbourhood(OldElem,NewElem);
	NewElem->BuildBlendConnectivity();
	
	delete OldElem;
	
	return NewElem;
}//method

void TPZChangeEl::AdjustNeighbourhood(TPZGeoEl* OldElem, TPZGeoEl*NewElem){
	for(int j = 0; j < OldElem->NSides(); j++)
	{ 
		TPZGeoElSide currSide(OldElem,j);
		TPZGeoElSide newElSide(NewElem,j);
		TPZGeoElSide neighbour(currSide.Neighbour());
		
		if (neighbour == currSide)
		{
			newElSide.SetNeighbour(newElSide);
			continue;
		}
		newElSide.SetNeighbour(neighbour);
		
		while (neighbour != currSide)
		{
			TPZGeoElSide next = neighbour.Neighbour();
			if (next == currSide)
			{
				neighbour.SetNeighbour(newElSide);
				break;
			}
			neighbour = neighbour.Neighbour();
		}
		
	}
}//method



// #include "pznoderep.h.h"

//          template class pzgeom::TPZNodeRep<8,TPZGeoBlend<TPZGeoCube> >;
//          template class pzgeom::TPZNodeRep<6,TPZGeoBlend<TPZGeoPrism> >;


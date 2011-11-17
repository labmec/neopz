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
using namespace pztopology;

TPZChangeEl::TPZChangeEl()
{
}

TPZChangeEl::~TPZChangeEl()
{
}

TPZGeoEl * TPZChangeEl::ChangeToQuadratic(TPZGeoMesh *Mesh, int ElemIndex)
{
	TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
    
    #ifdef DEBUG
    if(!OldElem)
    {
        std::cout << "Null geoel on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    #endif
    
    MElementType elType = OldElem->Type();
    int oldId = OldElem->Id();
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

TPZGeoEl* TPZChangeEl::DragQuarterPoints(TPZGeoMesh *Mesh, int ElemIndex, int targetSide)
{
    TPZGeoEl * gel = Mesh->ElementVec()[ElemIndex];
    
    #ifdef DEBUG
    if(!gel)
    {
        std::cout << "Null geoel on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    #endif
    
    if(gel->IsLinearMapping())
    {
        gel = TPZChangeEl::ChangeToQuadratic(Mesh, ElemIndex);
    }
    
    MElementType gelType = gel->Type();
    
    TPZStack<int> targetSideNodes, targetSideEdges;
    
    switch(gelType)
    {
        case (EPoint):
        {
            TPZPoint::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZPoint::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (EOned):
        {
            TPZLine::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZLine::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (ETriangle):
        {
            TPZTriangle::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZTriangle::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (EQuadrilateral):
        {
            TPZQuadrilateral::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZQuadrilateral::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (ETetraedro):
        {
            TPZTetrahedron::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZTetrahedron::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (EPiramide):
        {
            TPZPyramid::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZPyramid::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (EPrisma):
        {
            TPZPrism::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZPrism::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
		case (ECube):
        {
            TPZCube::LowerDimensionSides(targetSide,targetSideNodes,0);
            TPZCube::LowerDimensionSides(targetSide,targetSideEdges,1);
            
            break;
        }
        default:
        {
            std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
    }
    
    std::set<int> targetSideNodes_set, targetSideEdges_set, NOTtargetSideEdges_set;
    TPZGeoElSide targetSideEl(gel,targetSide);
    for(int nd = 0; nd < targetSideNodes.NElements(); nd++)
    {
        int tsNode = targetSideNodes[nd];
        targetSideNodes_set.insert(tsNode);
    }
    if(targetSideEl.Dimension() == 0)
    {
        targetSideNodes_set.insert(targetSide);
    }
    for(int ed = 0; ed < targetSideEdges.NElements(); ed++)
    {
        int tsEdge = targetSideEdges[ed];
        targetSideEdges_set.insert(tsEdge);
    }
    if(targetSideEl.Dimension() == 1)
    {
        targetSideNodes_set.insert(targetSide);
    }
    for(int sd = gel->NCornerNodes(); sd < gel->NSides(); sd++)
    {
        TPZGeoElSide gelSide(gel,sd);
        if(gelSide.Dimension() == 1 && targetSideEdges_set.find(sd) == targetSideEdges_set.end())
        {
            NOTtargetSideEdges_set.insert(sd);
        }
        else if(gelSide.Dimension() > 1)
        {
            break;
        }
    }
    
    std::set<int>::iterator it;
    for(it = NOTtargetSideEdges_set.begin(); it != NOTtargetSideEdges_set.end(); it++)
    {
        int edg = *it;
        TPZGeoElSide gelSide(gel,edg);
        
        /**
         * Embora o elemento quadratico possua mais nohs (NNodes), a topologia segue igual ao
         * elemento linear no qual o elemento quadratico foi baseado.
         *
         * Isso quer dizer que Linear(Geom::NNodes) < Quadratic(Geom::NNodes), mas Linear(Geom::NSides) = Quadratic(Geom::NSides)
         *
         * Com isso a numeracao dos nohs nos meios das arestas coincide com a numeracao dos lados unidimensionais do elemento.
         */
        int middleNodeId = edg; 
        //********************
        
std::cout << "Coords antes:\n";
Mesh->NodeVec()[middleNodeId].Print();
        
        int initNode = gelSide.SideNodeIndex(0);
        int finalNode = gelSide.SideNodeIndex(1); 
        int meshInitNode = gel->NodeIndex(initNode);
        int meshFinalNode = gel->NodeIndex(finalNode);
        double coordNear, coordFar;
        
        if(targetSideNodes_set.find(initNode) != targetSideNodes_set.end())//drag middle node to node 0 of this edge
        {
            for(int c = 0; c < 3; c++)
            {
                coordNear = Mesh->NodeVec()[meshInitNode].Coord(c);
                coordFar = Mesh->NodeVec()[meshFinalNode].Coord(c);
                Mesh->NodeVec()[middleNodeId].SetCoord(c, 0.75*coordNear + 0.25*coordFar);
            }
        }
        else if(targetSideNodes_set.find(finalNode) != targetSideNodes_set.end())//drag middle node to node 1 of this edge
        {
            for(int c = 0; c < 3; c++)
            {
                coordNear = Mesh->NodeVec()[meshFinalNode].Coord(c);
                coordFar = Mesh->NodeVec()[meshInitNode].Coord(c);
                Mesh->NodeVec()[middleNodeId].SetCoord(c, 0.75*coordNear + 0.25*coordFar);
            }
        }

std::cout << "Coords depois:\n";
Mesh->NodeVec()[middleNodeId].Print();

    }				
    
    return gel;
}

TPZGeoEl * TPZChangeEl::ChangeToGeoBlend(TPZGeoMesh *Mesh, int ElemIndex)
{
	TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
	if(!OldElem)
    {
		PZError << "Error at " << __PRETTY_FUNCTION__ << " - NULL geometric element.\n";
		return NULL;
	}
    MElementType oldType = OldElem->Type();
    int oldId = OldElem->Id();
    int oldMatId = OldElem->MaterialId();
    int nsides = OldElem->NSides();
    
    TPZVec<TPZGeoElSide> oldNeigh(nsides);
    for(int s = 0; s < nsides; s++)
    {   
        TPZGeoElSide mySide(OldElem, s);
        oldNeigh[s] = mySide.Neighbour();
    }
	
	const int nnodes = OldElem->NCornerNodes();
	TPZManVector<int> nodeindexes(nnodes);
	for(int i = 0; i < nnodes; i++)
    {
        nodeindexes[i] = OldElem->NodeIndex(i);
    }
    
    Mesh->DeleteElement(OldElem);
    
	TPZGeoEl * NewElem = Mesh->CreateGeoBlendElement(oldType, nodeindexes, oldMatId, oldId);

    for(int s = 0; s < nsides; s++)
    {   
        TPZGeoElSide neigh = oldNeigh[s];
        NewElem->SetNeighbour(s, neigh);
    }
    
	NewElem->BuildBlendConnectivity();
	
	return NewElem;
}



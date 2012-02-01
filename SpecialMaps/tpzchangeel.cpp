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
//------------------------------------------------------------------------------------------------------------

TPZChangeEl::~TPZChangeEl()
{
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl * TPZChangeEl::ChangeToQuadratic(TPZGeoMesh *Mesh, int ElemIndex)
{
	TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
    
    #ifdef DEBUG
    if(!OldElem)
    {
        std::cout << "Null geoel on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    if(!OldElem->IsLinearMapping())
    {
        return OldElem;
    }
    #endif
    
    TPZGeoEl * father = OldElem->Father();
    
    int midN, nsides = OldElem->NSides();
    
    //backingup oldElem neighbourhood
    TPZVec<REAL> Coord(3);
    TPZVec< std::vector<TPZGeoElSide> > neighbourhood(nsides);
    TPZVec<int> NodesSequence(0);
    for(int s = 0; s < nsides; s++)
    {
        neighbourhood[s].resize(0);
        TPZGeoElSide mySide(OldElem,s);
        TPZGeoElSide neighS = mySide.Neighbour();
        if(mySide.Dimension() == 0)
        {
            int oldSz = NodesSequence.NElements();
            NodesSequence.resize(oldSz+1);
            NodesSequence[oldSz] = OldElem->NodeIndex(s);
        }
        if(CreateMiddleNodeAtEdge(Mesh, ElemIndex, s, midN))
        {
            int oldSz = NodesSequence.NElements();
            NodesSequence.resize(oldSz+1);
            NodesSequence[oldSz] = midN;
        }
        while(mySide != neighS)
        {
            neighbourhood[s].push_back(neighS);
            neighS = neighS.Neighbour();
        }
    }
    
    MElementType elType = OldElem->Type();
    int oldId = OldElem->Id();
    int oldMatId = OldElem->MaterialId();
    
	TPZGeoEl * NewElem = NULL;
    
    /** Deleting OldElem */
    Mesh->DeleteElement(OldElem);

    switch(elType) /** Inserting New Element in Mesh */
    {
        case(EOned) :
        {             
            NewElem = new TPZGeoElRefPattern<TPZQuadraticLine>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }  
        case(ETriangle) :
        {            
            NewElem = new TPZGeoElRefPattern<TPZQuadraticTrig>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(EQuadrilateral) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticQuad>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(ETetraedro) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticTetra>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(EPiramide) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticPyramid>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }    
        case(EPrisma) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticPrism>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(ECube) :
        {
            NewElem = new TPZGeoElRefPattern<TPZQuadraticCube>(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        default :
        {
            DebugStop();
            break;
        }
    }
    
    if(father)
    {
        NewElem->SetFather(father);
    }
    
    for(int s = 0; s < nsides; s++)
    {
        TPZGeoEl * tempEl = NewElem;
        TPZGeoElSide tempSide(NewElem,s);
        int byside = s;
        for(unsigned int n = 0; n < neighbourhood[s].size(); n++)
        {
            TPZGeoElSide neighS = neighbourhood[s][n];
            tempEl->SetNeighbour(byside, neighS);
            tempEl = neighS.Element();
            byside = neighS.Side();
        }
        tempEl->SetNeighbour(byside, tempSide);
    }
        
    if(NewElem->HasSubElement())
    {
        //Mudar subelementos para TPZGeoElMapped
    }
    
	return NewElem;
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl* TPZChangeEl::ChangeToQuarterPoint(TPZGeoMesh *Mesh, int ElemIndex, int targetSide)
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
        targetSideEdges_set.insert(targetSide);
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
        
        int initNode = gelSide.SideNodeLocIndex(0);
        int finalNode = gelSide.SideNodeLocIndex(1);; 
        int meshInitNode = gelSide.SideNodeIndex(0);
        int meshFinalNode = gelSide.SideNodeIndex(1);
        /**
         * Embora o elemento quadratico possua mais nohs (NNodes), a topologia segue igual ao
         * elemento linear no qual o elemento quadratico foi baseado.
         *
         * Isso quer dizer que Linear(Geom::NNodes) < Quadratic(Geom::NNodes), mas Linear(Geom::NSides) = Quadratic(Geom::NSides)
         *
         * Com isso a numeracao dos nohs nos meios das arestas coincide com a numeracao dos lados unidimensionais do elemento.
         */
        int meshMiddleNode = gel->NodeIndex(edg);
        //********************
                                            
        double coordNear, coordFar;
        
        if(targetSideNodes_set.find(initNode) != targetSideNodes_set.end())//drag middle node to node 0 of this edge
        {
            for(int c = 0; c < 3; c++)
            {
                coordNear = Mesh->NodeVec()[meshInitNode].Coord(c);
                coordFar = Mesh->NodeVec()[meshFinalNode].Coord(c);
                Mesh->NodeVec()[meshMiddleNode].SetCoord(c, 0.75*coordNear + 0.25*coordFar);
            }
        }
        else if(targetSideNodes_set.find(finalNode) != targetSideNodes_set.end())//drag middle node to node 1 of this edge
        {
            for(int c = 0; c < 3; c++)
            {
                coordNear = Mesh->NodeVec()[meshFinalNode].Coord(c);
                coordFar = Mesh->NodeVec()[meshInitNode].Coord(c);
                Mesh->NodeVec()[meshMiddleNode].SetCoord(c, 0.75*coordNear + 0.25*coordFar);
            }
        }
    }				
    
    return gel;
}
//------------------------------------------------------------------------------------------------------------

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
//------------------------------------------------------------------------------------------------------------

bool TPZChangeEl::NearestNode(TPZGeoEl * gel, TPZVec<REAL> &x, int &meshNode, double tol)
{    
	meshNode = -1;
	bool IsNearSomeNode = false;
	
	TPZVec<REAL> nodeCoord(3);
	int nnodes = gel->NNodes();
	
	for(int n = 0; n < nnodes; n++)
	{
		double dist = 0.;
		gel->NodePtr(n)->GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			dist += (x[c] - nodeCoord[c])*(x[c] - nodeCoord[c]);
		}
		dist = sqrt(dist);
		
		if(dist <= tol)
		{
			meshNode = gel->NodeIndex(n);
			IsNearSomeNode = true;
			break;
		}
	}
	
	return IsNearSomeNode;
}
//------------------------------------------------------------------------------------------------------------


int TPZChangeEl::NearestNode(TPZGeoMesh * gmesh, TPZVec<REAL> &x, double tol)
{
	int meshNode = -1;
	
	TPZVec<REAL> nodeCoord(3);
	int nnodes = gmesh->NNodes();
	
	for(int n = 0; n < nnodes; n++)
	{
		double dist = 0.;
		gmesh->NodeVec()[n].GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			dist += (x[c] - nodeCoord[c])*(x[c] - nodeCoord[c]);
		}
		dist = sqrt(dist);
		
		if(dist <= tol)
		{
			meshNode = n;
			break;
		}
	}
	
	if(meshNode == -1)
	{
		std::cout << "Node not found for coordinates ( " << x[0] << " , " << x[1] << " , " << x[2] << " )" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		
		DebugStop();
	}
	
	return meshNode;
}
//------------------------------------------------------------------------------------------------------------

bool TPZChangeEl::CreateMiddleNodeAtEdge(TPZGeoMesh *Mesh, int ElemIndex, int edge, int &middleNodeId)
{
    TPZGeoEl * gel = Mesh->ElementVec()[ElemIndex];
    
    TPZGeoElSide myEdge(gel,edge);
    TPZGeoElSide neighEdge(myEdge.Neighbour());
    if(myEdge.Dimension() != 1)
    {
        return false;
    }
    
    TPZVec<REAL> n0Coord(3), n1Coord(3), middleCoordLocal(1), middleCoord(3);
    middleCoordLocal[0] = 0.;//middle node in edge (1D master element) is on qsi=0
    
    int nearestN;
    myEdge.X(middleCoordLocal,middleCoord);
    if(NearestNode(gel, middleCoord, nearestN, 1.E-8))
    {
        middleNodeId = nearestN;//a malha jah contem um noh na coordenada desejada (middleCoord)
        return true;
    }
    //else...
    while(neighEdge != myEdge)
    {
        neighEdge.X(middleCoordLocal,middleCoord);
        if(NearestNode(neighEdge.Element(), middleCoord, nearestN, 1.E-8))
        {
            middleNodeId = nearestN;//a malha jah contem um noh na coordenada desejada (middleCoord)
            return true;
        }
        neighEdge = neighEdge.Neighbour();
    }
    
    //if not returned true...
    TPZGeoNode midNode;
    midNode.SetCoord(middleCoord);
    
    /** Setting Midnodes Id's */
    int NewNodeId = Mesh->NNodes();
    Mesh->SetNodeIdUsed(NewNodeId);
    midNode.SetNodeId(NewNodeId);
    
    /** Allocating Memory for MidNodes and Pushing Them */
    middleNodeId = Mesh->NodeVec().AllocateNewElement();
    Mesh->NodeVec()[middleNodeId] = midNode;

    return true;
}



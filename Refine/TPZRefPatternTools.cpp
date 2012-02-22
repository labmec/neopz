/**
 * @file
 * @brief Contains the implementation of the TPZRefPatternTools methods. 
 */
/*
 *  Created by Cesar Lucci on 10/03/10.
 *  Copyright 2010 LabMeC. All rights reserved.
 *
 */

#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"

using namespace std;

TPZRefPatternTools::TPZRefPatternTools()
{
	
}

TPZRefPatternTools::~TPZRefPatternTools()
{
	
}

void TPZRefPatternTools::GetCompatibleRefPatterns(TPZGeoEl *gel, std::list<TPZAutoPointer<TPZRefPattern> > &refs)
{
	if(!gel) return;
	
	refs.clear();
	if(gel->HasSubElement())
	{
		refs.push_back(gel->GetRefPattern());
		return;
	}
	
	// first we build the refinement patterns associated with the neighbours of the current element
	int side, nsides, nnodes;
	nsides = gel->NSides();
	nnodes = gel->NCornerNodes();
	TPZManVector<TPZAutoPointer<TPZRefPattern>, 27> NeighSideRefPatternVec(nsides,0);
	
	for(side = nnodes; side < nsides; side++)
	{
		TPZGeoElSide gelside(gel, side);
		TPZGeoElSide neighside = gelside.Neighbour();
		if(!neighside.Element()){
			DebugStop();
		}
		while(neighside != gelside)
		{
			if(neighside.NSubElements() > 1)
			{
				TPZAutoPointer<TPZRefPattern> neighRefp = neighside.Element()->GetRefPattern();
				if(neighRefp)
				{
					TPZTransform trans = neighside.NeighbourSideTransform(gelside);
					NeighSideRefPatternVec[side] = neighRefp->SideRefPattern(neighside.Side(),trans);
					break;
				}
			}
			neighside = neighside.Neighbour();
		}
	}
	
	// having the refinement patterns associated with the sides, look for compatible refinement patterns
	std::list< TPZAutoPointer<TPZRefPattern> > gelReflist = gRefDBase.RefPatternList(gel->Type());
	std::list< TPZAutoPointer<TPZRefPattern> >::iterator gelReflistIt;
	
	for(gelReflistIt = gelReflist.begin(); gelReflistIt != gelReflist.end(); gelReflistIt++)
	{
		// compare the side refinement patterns
		for(side = nnodes; side < nsides; side++)
		{
			TPZAutoPointer<TPZRefPattern> GelSideRefPattern = (*gelReflistIt)->SideRefPattern(side);
			TPZAutoPointer<TPZRefPattern> NeighSideRefPattern = NeighSideRefPatternVec[side];
			
			if(GelSideRefPattern && NeighSideRefPattern)
			{				
				if(  !( (TPZRefPattern &)GelSideRefPattern == NeighSideRefPattern)  )
				{			
					break;
				}
			}
		}
		
		// if all refinement patterns are equal
		if(side == nsides)
		{
			refs.push_back(*gelReflistIt);
		}
	}
}

TPZAutoPointer<TPZRefPattern> TPZRefPatternTools::ModelRefPattern(TPZGeoEl *gel, std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > > &neighCorresp)
{
	neighCorresp.clear();
	TPZAutoPointer<TPZRefPattern> modelPat;//NULL
	
	int nnodes = gel->NNodes();
	int nsides = gel->NSides();
	
	bool neigh_isrefined, topolAreCompatibles;
	
	std::list< TPZAutoPointer<TPZRefPattern> > candidates = gRefDBase.RefPatternList(gel->Type());
	std::list< TPZAutoPointer<TPZRefPattern> >::iterator it;
	
	int side;
	for(it = candidates.begin(); it != candidates.end(); it++)
	{
		modelPat = *it;		
		for(side = nnodes; side < nsides; side++)
		{
			TPZGeoElSide gelside(gel, side);
			TPZGeoElSide neighside = gelside.Neighbour();
			
			neigh_isrefined = false;
			topolAreCompatibles = false;
			
			//Se meu lado eh aresta e nao tem vizinho, significa que nao precisa ser refinado!
			//Mas se o candidato a modelRefPattern quer refinar por este lado, jah nao serve... vamos para outro candidato!
			if(neighside == gelside && gelside.Dimension() == 1 && modelPat->SideRefPattern(side) && side != nsides-1)
			{
				neigh_isrefined = true;
				topolAreCompatibles = false;
				break;
			}
			//Se eu nao tenho vizinho pelo lado (nsides-1), o candidato a modelRefPattern pode refina-lo!!!
			else if(neighside == gelside && side == nsides-1)
			{
				continue;
			}
			
			while(neighside != gelside)
			{
				if(neighside.NSubElements() > 1)
				{
					neigh_isrefined = true;
					TPZAutoPointer<TPZRefPattern> neighRefp = neighside.Element()->GetRefPattern()->SideRefPattern(neighside.Side());
					TPZAutoPointer<TPZRefPattern> modelsideRefp = modelPat->SideRefPattern(side);
					
					if(neighRefp && modelsideRefp)
					{
						TPZTransform transBetweenNeigh = neighside.NeighbourSideTransform(gelside);
						TPZTransform transBetweenSides = gel->SideToSideTransform(side, nsides-1);
						
						std::map<int, int> pairedNodes;
						
						//correspondencias entre nohs entre LADO DO PADRAO MODELO e LADO DO PADRAO VIZINHO
						topolAreCompatibles = CompareTopologies(neighRefp, modelsideRefp, transBetweenNeigh, pairedNodes);
						
						if(topolAreCompatibles)
						{
							neighCorresp[side] = make_pair(neighside.Element(), pairedNodes);
							break;
						}
					}
				}
				neighside = neighside.Neighbour();
			}
			if( (neigh_isrefined && !topolAreCompatibles) || (gelside.Dimension() == 1 && !neigh_isrefined && modelPat->SideRefPattern(side)) )
			{
				neighCorresp.clear();
				break;
			}
		}		
		if(side == nsides)
		{
			return modelPat;	
		}
	}
	
	return NULL;
}

TPZAutoPointer<TPZRefPattern> TPZRefPatternTools::PerfectMatchRefPattern(TPZGeoEl *gel, TPZVec<int> &sidestorefine)
{
	if(!gel)
	{
		return NULL;
	}
	
	std::list<TPZAutoPointer<TPZRefPattern> > patlist;
	TPZRefPatternTools::GetCompatibleRefPatterns(gel, patlist);
	
	std::list<TPZAutoPointer<TPZRefPattern> >::iterator it;
	for(it = patlist.begin(); it != patlist.end(); it++)
	{
		if( !(*it) )
		{
			continue;	
		}
		TPZGeoEl * gelTemp = (*it)->Element(0);
		int ncorners = gelTemp->NCornerNodes();
		int nsides = gelTemp->NSides();
		if(gelTemp->Dimension() == 3)
		{
			nsides--;
		}
		
		int is;
		for(is = ncorners; is < nsides; is++)
		{
			if( sidestorefine[is] != (*it)->NSideNodes(is) )
			{
				break;	
			}
		}
		if(is == nsides)
		{
			return (*it);	
		}
	}
	
	return NULL;
}

TPZAutoPointer<TPZRefPattern> TPZRefPatternTools::PerfectMatchRefPattern(TPZGeoEl *gel)
{
	if(!gel)
	{
		return NULL;
	}
	TPZVec<int> sidestorefine;
	
	bool thereIsAnyNeighbourRefined = TPZRefPatternTools::SidesToRefine(gel, sidestorefine);
	if(!thereIsAnyNeighbourRefined)
	{
		return 0;
	}
	TPZAutoPointer<TPZRefPattern> pat = PerfectMatchRefPattern(gel, sidestorefine);
	
	if(!pat)
	{
		//if no refpattern matches perfectly...
		std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > > neighCorresp;
		TPZAutoPointer<TPZRefPattern> modelPat = TPZRefPatternTools::ModelRefPattern(gel, neighCorresp);	
		if(modelPat)
		{
			//modelPat->Print();
			pat = TPZRefPatternTools::DragModelPatNodes(gel, modelPat, neighCorresp);
		}
		else
		{
			std::cout << "\nModelRefPattern NOT FOUND in " << __PRETTY_FUNCTION__ << std::endl;
			std::cout << "You should create it and add in Refinement Patterns Folder!" << std::endl;
			std::cout << "Open file ModelRefPatternNOTFOUND.vtk in Paraview to see the neighbourhood\n";
            
            std::ofstream outNotFound("ModelRefPatternNOTFOUND.vtk");
            TPZVTKGeoMesh::PrintGMeshVTKneighbourhood(gel->Mesh(), gel->Id(), outNotFound);
			
            DebugStop();
		}
	}
	
	return pat;
}

TPZAutoPointer<TPZRefPattern> TPZRefPatternTools::DragModelPatNodes(TPZGeoEl * gel, TPZAutoPointer<TPZRefPattern> modelPat, std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > > &neighCorresp)
{
	TPZAutoPointer<TPZRefPattern> modelPat_copy = new TPZRefPattern(modelPat);
	int nnodes = gel->NNodes();
	int nsides = gel->NSides();
	
	for(int side = nnodes; side < nsides; side++)
	{
		std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > >::iterator correspIt = neighCorresp.find(side);
		if(correspIt == neighCorresp.end())
		{
			continue;
		}
		TPZGeoEl * neighStored = correspIt->second.first;
		TPZGeoElSide gelside(gel, side);
		TPZGeoElSide neighside = gelside.Neighbour();
		
		while(neighside != gelside)
		{
			if(neighside.Element() == neighStored)
			{
				//int SideOfNeighSide = neighside.Side();
				TPZAutoPointer<TPZRefPattern> neighRefp = neighside.Element()->GetRefPattern()->SideRefPattern(neighside.Side());
				TPZAutoPointer<TPZRefPattern> modelsideRefp = modelPat_copy->SideRefPattern(side);
				
				TPZTransform transBetweenNeigh = neighside.NeighbourSideTransform(gelside);
				TPZTransform transBetweenSides = gel->SideToSideTransform(side, nsides-1);
				
				//correspondencias entre nohs entre LADO DO PADRAO MODELO e LADO DO PADRAO VIZINHO
				std::map<int, int> pairedNodes = correspIt->second.second;
				
				//correspondencias entre nohs entre LADO DO PADRAO MODELO e o proprio PADRAO MODELO
				std::map<int, int> pairSideNodes;
				PairMeshesNodesMatchingCoordinates(modelsideRefp->RefPatternMesh(), modelPat_copy->RefPatternMesh(), transBetweenSides, pairSideNodes);
				std::map<int, int>::iterator it, jt;
				for(it = pairedNodes.begin(); it != pairedNodes.end(); it++)
				{
					TPZVec<REAL> NodeCoord(3,0.), sideNodeCoord(3,0.), neighNodeCoord(3,0.);
					
					neighRefp->RefPatternMesh().NodeVec()[it->first].GetCoordinates(neighNodeCoord);
					transBetweenNeigh.Apply(neighNodeCoord, sideNodeCoord);
					transBetweenSides.Apply(sideNodeCoord, NodeCoord);
					
					jt = pairSideNodes.find(it->second);
					if(jt != pairSideNodes.end())
					{
						int gNode = jt->second;
						for(int n = 0; n < 3; n++)
						{
							modelPat_copy->RefPatternMesh().NodeVec()[gNode].SetCoord(n, NodeCoord[n]);
						}
					}
				}
				
				break;
			}
			neighside = neighside.Neighbour();
		}
	}
	modelPat_copy->ComputeTransforms();
	modelPat_copy->ComputePartition();
	modelPat_copy->GenerateSideRefPatterns(); 
	if(!gRefDBase.FindRefPattern(modelPat_copy))
	{
		gRefDBase.InsertRefPattern(modelPat_copy);
		modelPat_copy->InsertPermuted();
	}
	
	return modelPat_copy;
}

bool TPZRefPatternTools::CompareTopologies(TPZAutoPointer<TPZRefPattern> refA, TPZAutoPointer<TPZRefPattern> refB, TPZTransform &fromAtoB, std::map<int, int> &pairedNodes)
{	
	pairedNodes.clear();
	
	if(refA->Type() != refB->Type())
	{
		//Their are NOT topologicaly compatibles
		return false;
	}	
	
	TPZGeoMesh meshA = refA->RefPatternMesh();
	TPZGeoMesh meshB = refB->RefPatternMesh();
	
	if(meshA.NNodes() != meshB.NNodes() || meshA.NElements() != meshB.NElements())
	{
		//Their are NOT topologicaly compatibles
		return false;
	}
	
	//Subgrouping sub-elements by types
	std::map<MElementType, std::list<TPZGeoEl*> > eltypeA, eltypeB;
	for(int el = 1; el < meshA.NElements(); el++)//Father is NOT included!
	{
		TPZGeoEl *gelA = meshA.ElementVec()[el];
		TPZGeoEl *gelB = meshB.ElementVec()[el];
		
		eltypeA[gelA->Type()].push_back(gelA);
		eltypeB[gelB->Type()].push_back(gelB);
	}
	if(eltypeA.size() != eltypeB.size())
	{
		//Their are NOT topologicaly compatibles
		return false;
	}
	std::map<MElementType, std::list<TPZGeoEl*> >::iterator eltypeAit, eltypeBit;
	for(eltypeAit = eltypeA.begin(); eltypeAit != eltypeA.end(); eltypeAit++)
	{
		eltypeBit = eltypeB.find(eltypeAit->first);
		if( (eltypeBit == eltypeB.end()) || (eltypeAit->second.size() != eltypeBit->second.size()) )
		{
			//Their are NOT topologicaly compatibles
			return false;
		}
	}
	
	//PARTY BEGINS!!!!!!!!!	
	PairMeshesCornerNodesMatchingCoordinates(meshA, meshB, fromAtoB, pairedNodes);
	
	TPZGeoEl* fatherA = meshA.ElementVec()[0];
	int unpairedNNodes = meshA.NNodes() - pairedNodes.size();
	
#ifdef DEBUG
	if((int) pairedNodes.size() < fatherA->NNodes())
	{
		std::cout << "\nThere is something going wrong with meshA and/or meshB in " << __PRETTY_FUNCTION__ << std::endl;
		std::cout << "father->ConnerNodes should be paired at least!!!\n\n";
		DebugStop();
	}
#endif
	
	std::map<int, int>::iterator pairedNodesIT;
	
	if(unpairedNNodes > 0)
	{
		for(int s = fatherA->NNodes(); s < fatherA->NSides(); s++)//pairing fatherA~fatherB (Edges)Nodes
		{
			TPZGeoElSide Aside(fatherA, s);
			if(Aside.Dimension() == 1)
			{
				pairedNodesIT = pairedNodes.find(Aside.SideNodeIndex(0));
				int nodeIniA = pairedNodesIT->first;
				int nodeIniB = pairedNodesIT->second;
				
				pairedNodesIT = pairedNodes.find(Aside.SideNodeIndex(1));
				int nodeFinA = pairedNodesIT->first;
				int nodeFinB = pairedNodesIT->second;			
				
				TPZManVector<int> NodesHuntedA;
				NodesHunter(meshA, NodesHuntedA, nodeIniA, nodeFinA);
				
				TPZManVector<int> NodesHuntedB;
				NodesHunter(meshB, NodesHuntedB, nodeIniB, nodeFinB);
				
				for(int n = 1; n < NodesHuntedA.size()-1; n++)
				{
					int nA = NodesHuntedA[n];
					int nB = NodesHuntedB[n];
					
					if(pairedNodes.find(nA) == pairedNodes.end())
					{
						pairedNodes[nA] = nB;
						unpairedNNodes--;					
					}
				}
			}
		}
	}
	
	//pairing father->(Faces&Volume)Nodes AND SUB-ELEMENTS!!!
	int nRemainingSubels = refA->NSubElements();
	TPZVec<bool> AsubelementIsAlreadyPaired(meshA.NElements(), false);
	TPZVec<bool> BsubelementIsAlreadyPaired(meshB.NElements(), false);
	
#define IwantPairedElementsToo
#ifdef IwantPairedElementsToo
	std::map<int, int> pairedElements;
	pairedElements[0] = 0;//fathers are already paired at this line of code
#endif
	
	std::list<TPZGeoEl*>::iterator elAit, elBit;
	std::list< std::list<TPZGeoEl*>::iterator > pointedB;
	TPZVec<int> Anodes(0);
	TPZVec< TPZVec<int> > Bnodes(0);
	std::list< TPZVec<int> > BthatMatch;
	int pairedNodesCount, unsuccessfulLoop = 0;
	while(unsuccessfulLoop < 2 && nRemainingSubels > 0)
	{
		for(eltypeAit = eltypeA.begin(); eltypeAit != eltypeA.end(); eltypeAit++)
		{
			eltypeBit = eltypeB.find(eltypeAit->first);
			for(elAit = eltypeAit->second.begin(); elAit != eltypeAit->second.end(); elAit++)
			{
				TPZGeoEl * A = *elAit;
				Anodes.Resize(A->NNodes());
				pairedNodesCount = 0;
				for(int n = 0; n < A->NNodes(); n++)
				{
					Anodes[n] = A->NodeIndex(n);
					pairedNodesIT = pairedNodes.find(Anodes[n]);
					if(pairedNodesIT != pairedNodes.end())
					{
						pairedNodesCount++;
					}
				}
				if(pairedNodesCount > 1)
				{
					BthatMatch.clear();
					pointedB.clear();
					for(elBit = eltypeBit->second.begin(); elBit != eltypeBit->second.end(); elBit++)
					{
						TPZGeoEl * B = *elBit;
						TPZRefPatternTools::GetGelPermutations(B, Bnodes);
						for(int p = 0; p < Bnodes.NElements(); p++)
						{
							bool goodCandidate = true;
							for(int n = 0; n < Anodes.NElements(); n++)
							{
								int nA = Anodes[n];
								int nB = Bnodes[p][n];
								
								pairedNodesIT = pairedNodes.find(nA);
								if(pairedNodesIT != pairedNodes.end() && pairedNodesIT->second != nB)
								{
									goodCandidate = false;
									break;
								}
							}
							if(goodCandidate)
							{
								BthatMatch.push_back(Bnodes[p]);
								pointedB.push_back(elBit);  
							}
						}
					}//for elBit
				}//if(pairedNodesCount > 1)
				if(pairedNodesCount > 1 && BthatMatch.size() == 1)//se temos apenas 1 candidato da lista de B
				{
					for(int n = 0; n < Anodes.NElements(); n++)//pareando os nohs ainda nao pareados e deletando da lista os elementos pareados
					{
						pairedNodesIT = pairedNodes.find(Anodes[n]);
						if(pairedNodesIT == pairedNodes.end() && unpairedNNodes > 0)
						{
							pairedNodes[Anodes[n]] = (*BthatMatch.begin())[n];
							unpairedNNodes--;
						}
					}
					nRemainingSubels--;
					
#ifdef IwantPairedElementsToo
					int Aid = (*(elAit))->Id();
					int Bid = (*(*(pointedB.begin())))->Id();
					pairedElements[Aid] = Bid;
#endif
					
					eltypeAit->second.erase(elAit);
					eltypeBit->second.erase(*pointedB.begin());
					
					unsuccessfulLoop = 0;
					break;
				}
			}//for elAit
		}//for eltypeAit
		unsuccessfulLoop++;
	}//while
	
	if(nRemainingSubels == 0 && unpairedNNodes == 0)
	{
		//Their ARE topologicaly compatibles
		return true;
	}
	
	return false;
}

void TPZRefPatternTools::PairMeshesCornerNodesMatchingCoordinates(TPZGeoMesh meshA, TPZGeoMesh meshB, TPZTransform fromAtoB, std::map<int, int> &pairedNodes)
{
	TPZVec<REAL> ANodeCoord(3,0.), BNodeCoord(3,0.);
	
	int Annodes = meshA.ElementVec()[0]->NNodes();
	int Bnnodes = meshB.ElementVec()[0]->NNodes();
	for(int nA = 0; nA < Annodes; nA++)
	{
		meshA.NodeVec()[meshA.ElementVec()[0]->NodeIndex(nA)].GetCoordinates(ANodeCoord);
		fromAtoB.Apply(ANodeCoord, BNodeCoord);
		
		TPZVec<REAL> BNodeCoord_compare(3,0.);
		for(int nB = 0; nB < Bnnodes; nB++)
		{
			double dif = 0.;
			meshB.NodeVec()[meshB.ElementVec()[0]->NodeIndex(nB)].GetCoordinates(BNodeCoord_compare);
			for(int c = 0; c < 3; c++)
			{
				dif += fabs(BNodeCoord[c] - BNodeCoord_compare[c]);
			}
			if(dif < 1.E-10)
			{
				int nA_Id = meshA.NodeVec()[meshA.ElementVec()[0]->NodeIndex(nA)].Id();
				int nB_Id = meshB.NodeVec()[meshB.ElementVec()[0]->NodeIndex(nB)].Id();
				pairedNodes[nA_Id] = nB_Id;
				
				break;
			}
		}
	}
}

void TPZRefPatternTools::PairMeshesNodesMatchingCoordinates(TPZGeoMesh meshA, TPZGeoMesh meshB, TPZTransform fromAtoB, std::map<int, int> &pairedNodes)
{
	TPZVec<REAL> ANodeCoord(3,0.), BNodeCoord(3,0.);
	
	int Annodes = meshA.NNodes();
	int Bnnodes = meshB.NNodes();
	for(int nA = 0; nA < Annodes; nA++)
	{
		meshA.NodeVec()[nA].GetCoordinates(ANodeCoord);
		fromAtoB.Apply(ANodeCoord, BNodeCoord);
		
		TPZVec<REAL> BNodeCoord_compare(3,0.);
		for(int nB = 0; nB < Bnnodes; nB++)
		{
			double dif = 0.;
			meshB.NodeVec()[nB].GetCoordinates(BNodeCoord_compare);
			for(int c = 0; c < 3; c++)
			{
				dif += fabs(BNodeCoord[c] - BNodeCoord_compare[c]);
			}
			if(dif < 1.E-10)
			{
				int nA_Id = meshA.NodeVec()[nA].Id();
				int nB_Id = meshB.NodeVec()[nB].Id();
				pairedNodes[nA_Id] = nB_Id;
				
				break;
			}
		}
	}
}

std::string TPZRefPatternTools::BuildRefPatternModelName(TPZRefPattern &refp)
{
	std::string refpTypeName;
	std::stringstream rpName;
	
	if(&refp == NULL)
	{
		std::cout << "Null refpattern parameter on " << __PRETTY_FUNCTION__ << std::endl;
		return refpTypeName;
	}
	
	TPZGeoEl *gel = refp.Element(0);
	int ncorners = gel->NCornerNodes();
	int nsides = gel->NSides();
	
	std::string perfix = gel->TypeName();
	for(int i = 0; i < 3; i++)
	{
		rpName <<  perfix[i];
	}
	for(int s = 0; s < nsides; s++)
	{
		if(s < ncorners)
		{
			rpName << "0";
		}
		else 
		{
			rpName << refp.NSideNodes(s);	
		}
	}
	rpName >> refpTypeName;
	
	if(refpTypeName.length() == 0)
	{
		refpTypeName = nonInitializedName;
	}
	
	return refpTypeName;
}

std::string TPZRefPatternTools::BuildRefPatternModelName(TPZAutoPointer<TPZRefPattern> refp)
{
	std::string refpTypeName;
	std::stringstream rpName;
	
	if(!refp)
	{
		std::cout << "Null refpattern parameter on " << __PRETTY_FUNCTION__ << std::endl;
		return refpTypeName;
	}
	
	TPZGeoEl *gel = refp->Element(0);
	int ncorners = gel->NCornerNodes();
	int nsides = gel->NSides();
	
	std::string perfix = gel->TypeName();
	for(int i = 0; i < 3; i++)
	{
		rpName <<  perfix[i];
	}
	for(int s = 0; s < nsides; s++)
	{
		if(s < ncorners)
		{
			rpName << "0";
		}
		else 
		{
			rpName << refp->NSideNodes(s);	
		}
	}
	rpName >> refpTypeName;
	
	if(refpTypeName.length() == 0)
	{
		refpTypeName = nonInitializedName;
	}
	
	return refpTypeName;
}

std::string TPZRefPatternTools::BuildRefPatternModelName(TPZGeoEl *gel)
{
	std::string refpTypeName = "";
	std::stringstream rpName;
	
	if(!gel)
	{
		std::cout << "Null geoelement parameter on " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	
	int nsides = gel->NSides();
	TPZVec<int> sidesToRefine;
	bool thereIsAnyNeighbourRefined = TPZRefPatternTools::SidesToRefine(gel, sidesToRefine);
	if(!thereIsAnyNeighbourRefined)
	{
		return refpTypeName;
	}
	
	std::string perfix = gel->TypeName();
	for(int i = 0; i < 3; i++)
	{
		rpName <<  perfix[i];
	}
	
	for(int s = 0; s < nsides; s++)
	{
		rpName << sidesToRefine[s];
	}
	rpName >> refpTypeName;
	
	if(refpTypeName.length() == 0)
	{
		refpTypeName = nonInitializedName;
	}
	
	return refpTypeName;
}

bool TPZRefPatternTools::SidesToRefine(TPZGeoEl *gel, TPZVec<int> &sidestoref)
{
	bool thereIsAnyNeighbourRefined = false;
	
	int ncorners = gel->NCornerNodes();
	int nsides = gel->NSides();
	sidestoref.Resize(nsides, 0);
	
	for(int s = ncorners; s < nsides; s++)
	{
		TPZGeoElSide gelside(gel, s);
		TPZGeoElSide neighside = gelside.Neighbour();
		if(!neighside.Exists())
		{
			break;
		}
		while(neighside != gelside)
		{
			if(neighside.Element()->HasSubElement() && neighside.Element()->NSideSubElements2(neighside.Side()) > 1)
			{
				thereIsAnyNeighbourRefined = true;
				
				int ns = neighside.Side();
				TPZVec<int> MidNodesIndexes;
				
				TPZAutoPointer<TPZRefPattern> elrefpattern = neighside.Element()->GetRefPattern();				
				TPZAutoPointer<TPZRefPattern> refSide = elrefpattern->SideRefPattern(ns);
				if (!refSide)
				{
					std::cout << "An element is refined but the pattern has no side refpattern along side " << ns << "\n";
					neighside.Element()->GetRefPattern()->Print(cout);
					
					DebugStop();
				}
				elrefpattern->SideNodes(ns, MidNodesIndexes);
				
				sidestoref[s] = MidNodesIndexes.NElements();
				
				break;
			} 
			neighside = neighside.Neighbour();
		}
	}
	
	return thereIsAnyNeighbourRefined;
}

void TPZRefPatternTools::RefineDirectional(TPZGeoEl *gel, std::set<int> &matids)
{
	if(gel->HasSubElement())
	{
		return;
	}
	int matid = gel->MaterialId();
	if(matids.count(matid)) return;
	TPZManVector<int,27> sidestorefine(gel->NSides(), 0);
	TPZManVector<int,27> cornerstorefine(gel->NSides(), 0);
	
	//Look for corners which are on the boundary
	int numrefribs = 0;
	for(int in=0; in<gel->NCornerNodes(); in++)
	{
		TPZGeoElSide gels(gel,in);
		TPZGeoElSide neigh(gels.Neighbour());
		while(gels != neigh)
		{
			if(matids.count(neigh.Element()->MaterialId()))
			{
				cornerstorefine[in] = 1;
				break;
			}
			neigh = neigh.Neighbour();
		}
	}
	
	// look for ribs which touch the boundary but which do no lay on the boundary
	for(int is=gel->NCornerNodes(); is<gel->NSides(); is++)
	{
		// we are only interested in ribs
		if(gel->SideDimension(is) != 1)
		{
			continue;
		}
		
		// the side is a candidate if it contains a corner which is neighbour of the boundary condition
		if(cornerstorefine[gel->SideNodeLocIndex(is,0)] || cornerstorefine[gel->SideNodeLocIndex(is,1)])
		{
			sidestorefine[is] = 1;
			numrefribs++;
			TPZGeoElSide gels(gel,is);
			TPZGeoElSide neigh(gels.Neighbour());
			while(neigh != gels)
			{
				// if this condition is true the rib lies on the boundary
				if(matids.count(neigh.Element()->MaterialId()))
				{
					sidestorefine[is] = 0;
					numrefribs--;
					break;
				}
				neigh = neigh.Neighbour();
			}
		}
	}
	if(!numrefribs)
	{
		return;
	}
	
	TPZAutoPointer<TPZRefPattern> patt = TPZRefPatternTools::PerfectMatchRefPattern(gel, sidestorefine);
	if(patt)
	{
		gel->SetRefPattern(patt);
		TPZManVector<TPZGeoEl *> subel;
		gel->Divide(subel);
	}
	else
	{		
#ifdef WIN32
		DebugStop();
#endif
		std::cout << "|"; std::cout.flush();
		std::ofstream arquivo ("NotListedPatterns.txt",std::ios::app);
		std::list<TPZAutoPointer<TPZRefPattern> >::iterator it;
		arquivo << "Compatible refinement patterns\n";
		
		arquivo << std::endl;
		arquivo << "Element Type : " << gel->Type() << ' ' << gel->TypeName() << std::endl;
		arquivo << "Sides selected for refinement :" << std::endl;
		int i;
		for (i=0 ; i<gel->NSides() ; i++)
		{
			/*		
            if(cornerstorefine[i] == 1)
			{
             arquivo << " " << i << " ";
            }     
            */
			if (sidestorefine[i] == 1) {
				arquivo << " " << i << " " ;
			}
		}
		gel->Print(arquivo);
		int in;
		arquivo << std::endl;
		arquivo << "Neighbouring information \n";
		for(in=0; in<gel->NSides(); in++)
		{
			arquivo << "Side : " << in << " ";
			TPZGeoElSide gels(gel,in);
			arquivo << "Dim " << gels.Dimension() << " ";
			TPZGeoElSide neigh(gels.Neighbour());
			while(gels != neigh)
			{
				if(matids.count(neigh.Element()->MaterialId()))
				{
					arquivo << neigh.Element()->Id() << "-l-" << neigh.Side() << " ";
					if (neigh.Side() == 9 && gel->Type() == ETetraedro)
					{
						arquivo << "Teje pego meliante..." << std::endl;
						neigh.Element()->Print(arquivo);
					}
				}
				neigh = neigh.Neighbour();
			}
			arquivo << std::endl;
		}
		arquivo << std::endl;
		
		arquivo << "Element information : " << gel->Index() << std::endl;
		arquivo << "Vizinhos dos lados maraados para refinamento:" << std::endl;
		for (i=0 ; i<gel->NSides() ; i++)
		{
			if(cornerstorefine[i] == 1 || sidestorefine[i] == 1)
			{
				TPZGeoElSide gelside (gel,i);
				TPZGeoElSide neigh = gelside.Neighbour();
				while (neigh != gelside)
				{
					arquivo << "*********** my side = " << i << " neighside " << neigh.Side() << std::endl;
					neigh.Element()->Print(arquivo);
					neigh = neigh.Neighbour();
				}
			}
		}
		arquivo << std::endl << std::endl << std::endl << std::endl;
	}
	
	return;
}

void TPZRefPatternTools::RefineDirectional(TPZGeoEl *gel, std::set<int> &matids, int gelMat)
{
	if(gel->HasSubElement())
	{
		return;	
	}
	int matid = gel->MaterialId();
	if(matids.count(matid)) return;
	TPZManVector<int,27> sidestorefine(gel->NSides(), 0);
	TPZManVector<int,27> cornerstorefine(gel->NSides(), 0);
	
	//Look for corners which are on the boundary
	int numrefribs = 0;
	for(int in=0; in<gel->NCornerNodes(); in++)
	{
		TPZGeoElSide gels(gel,in);
		TPZGeoElSide neigh(gels.Neighbour());
		while(gels != neigh)
		{
			if(matids.count(neigh.Element()->MaterialId()))
			{
				cornerstorefine[in] = 1;
				break;
			}
			neigh = neigh.Neighbour();
		}
	}
	
	// look for ribs which touch the boundary but which do no lay on the boundary
	for(int is = gel->NCornerNodes(); is < gel->NSides(); is++)
	{
		// we are only interested in ribs
		if(gel->SideDimension(is) != 1)
		{
			continue;
		}
		
		// the side is a candidate if it contains a corner which is neighbour of the boundary condition
		if(cornerstorefine[gel->SideNodeLocIndex(is,0)] || cornerstorefine[gel->SideNodeLocIndex(is,1)])
		{
			sidestorefine[is] = 1;
			numrefribs++;
			TPZGeoElSide gels(gel,is);
			TPZGeoElSide neigh(gels.Neighbour());
			while(neigh != gels)
			{
				// if this condition is true the rib lies on the boundary
				if(matids.count(neigh.Element()->MaterialId()))
				{
					sidestorefine[is] = 0;
					numrefribs--;
					break;
				}
				neigh = neigh.Neighbour();
			}
		}
	}
	if(!numrefribs)
	{
		return;
	}
	
	TPZAutoPointer<TPZRefPattern> patt = TPZRefPatternTools::PerfectMatchRefPattern(gel, sidestorefine);
	if(patt)
	{
		gel->SetMaterialId(gelMat+1);
		gel->SetRefPattern(patt);
		TPZManVector<TPZGeoEl *> subel;
		gel->Divide(subel);
		gel->SetMaterialId(gelMat);
	}
	else
	{		
		std::cout << "|"; std::cout.flush();
		std::ofstream arquivo ("NotListedPatterns.txt",std::ios::app);
		std::list<TPZAutoPointer<TPZRefPattern> >::iterator it;
		arquivo << "Compatible refinement patterns\n";
		
		arquivo << std::endl;
		arquivo << "Element Type : " << gel->Type() << ' ' << gel->TypeName() << std::endl;
		arquivo << "Sides selected for refinement :" << std::endl;
		int i;
		for (i=0 ; i<gel->NSides() ; i++)
		{
			if(cornerstorefine[i] == 1)
			{
				arquivo << " " << i << " ";
			}
			if (sidestorefine[i] == 1) {
				arquivo << " " << i << " " ;
			}
		}
		gel->Print(arquivo);
		int in;
		arquivo << std::endl;
		arquivo << "Neighbouring information \n";
		for(in=0; in<gel->NSides(); in++)
		{
			arquivo << "Side : " << in << " ";
			TPZGeoElSide gels(gel,in);
			arquivo << "Dim " << gels.Dimension() << " ";
			TPZGeoElSide neigh(gels.Neighbour());
			while(gels != neigh)
			{
				if(matids.count(neigh.Element()->MaterialId()))
				{
					arquivo << neigh.Element()->Id() << "-l-" << neigh.Side() << " ";
					if (neigh.Side() == 9 && gel->Type() == ETetraedro)
					{
						arquivo << "Teje pego meliante..." << std::endl;
						neigh.Element()->Print(arquivo);
					}
				}
				neigh = neigh.Neighbour();
			}
			arquivo << std::endl;
		}
		arquivo << std::endl;
		
		arquivo << "Element information : " << gel->Index() << std::endl;
		arquivo << "Vizinhos dos lados maraados para refinamento:" << std::endl;
		for (i=0 ; i<gel->NSides() ; i++)
		{
			if(cornerstorefine[i] == 1 || sidestorefine[i] == 1)
			{
				TPZGeoElSide gelside (gel,i);
				TPZGeoElSide neigh = gelside.Neighbour();
				while (neigh != gelside)
				{
					arquivo << "*********** my side = " << i << " neighside " << neigh.Side() << std::endl;
					neigh.Element()->Print(arquivo);
					neigh = neigh.Neighbour();
				}
			}
		}
		arquivo << std::endl << std::endl << std::endl << std::endl;
	}	
	
	return;
}

void TPZRefPatternTools::RefineUniformIfNeighMat(TPZGeoEl *gel, std::set<int> &matids)
{
	int nsides = gel->NSides();
	for(int s = 0; s < nsides; s++)
	{
		TPZGeoElSide gelside(gel,s);
		TPZGeoElSide neighside(gelside.Neighbour());
		while(gelside != neighside)
		{
			if(matids.count(neighside.Element()->MaterialId()))
			{
				TPZAutoPointer<TPZRefPattern> refP = gRefDBase.GetUniformRefPattern(gel->Type());
				if(refP)
				{
					TPZVec<TPZGeoEl*> sons;
					gel->SetRefPattern(refP);
					gel->Divide(sons);
					
					return;
				}
				else
				{
					std::cout << "Uniform refpattern was not found!" << std::endl;
					std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
					return;
				}
				
			}
			neighside = neighside.Neighbour();
		}
	}
}

bool TPZRefPatternTools::ConstJacobian(TPZGeoElSide gelside, REAL tol)
{
	TPZIntPoints * rule = gelside.Element()->CreateSideIntegrationRule(gelside.Side(),5);
	int dim = gelside.Dimension();
	TPZManVector<REAL,4> pt(dim,0.);
	REAL weight;
	int np = rule->NPoints();
	rule->Point(0, pt, weight);
	TPZFNMatrix<9> jacobian(dim,dim,0.),jacinv(dim,dim,0.);
	TPZFNMatrix<9> axes(3,3,0.);
	REAL detjac;
	gelside.Jacobian(pt, jacobian, axes , detjac, jacinv);
	int ip;
	for(ip=1; ip<np; ip++)
	{
		rule->Point(ip, pt, weight);
		TPZFNMatrix<9> jacobian2(dim,dim,0.),jacinv2(dim,dim,0.);
		REAL detjac2;
		TPZFNMatrix<9> axes2(3,3,0.);
		gelside.Jacobian(pt, jacobian2, axes2 , detjac2, jacinv2);
		jacobian2 -= jacobian;
		axes2 -= axes;
		detjac2 -= detjac;
		jacinv2 -= jacinv;
		REAL diff = Norm(jacobian2)+Norm(axes2)+Norm(jacinv2)+fabs(detjac2);
		if(diff > 1.e-8) 
		{
			delete rule;
			return false;
		}
	}
	delete rule;
	
	return true;
}

void TPZRefPatternTools::TransformationTest(TPZRefPattern * refp)
{
	int isub,sd,ip;
	
	//x1 no filho deformado, x2 no pai deformado
	TPZManVector<REAL,3> x1(3),pf(3),x2(3),xpf(3,0.);
	REAL weight;
	TPZGeoEl *father = refp->Element(0);//pai
	int dimfatside,fatside,nsides;
	int dimfat = father->Dimension();
	
	TPZGeoEl *subel;
	int nsubs = refp->NSubElements();
	for(isub=0;isub<nsubs;isub++)
	{
		subel  = refp->Element(isub+1);
		nsides = subel->NSides();
		for(sd=0; sd<nsides; sd++)
		{
			TPZGeoElSide subside(subel,sd);
			if(!ConstJacobian(subside,1.e-6))
			{
				continue;
			}
			if( sd<subel->NNodes() && refp->IsFatherNeighbour(TPZGeoElSide(subel,sd),father) ) continue;
			
			int dims = subel->SideDimension(sd);
			int dimsub = subel->Dimension();
			
			//regra de integracao para o espaco paramcorico do lado do sub-elemento
			TPZIntPoints *rule = subel->CreateSideIntegrationRule(sd,5);
			TPZVec<int> order(dims,2);
			TPZManVector<REAL,3> point(dims,0.),point2(dimsub),pointparamfather(dimfat,0.);
			rule->SetOrder(order);
			
			for(ip=0;ip<rule->NPoints();ip++)
			{
				//ponto no espaco paramcorico do lado do filho
				rule->Point(ip,point,weight);
				TPZTransform sidet(dimsub);//transformacao unitcoia
				if(dims < dimsub)
				{
					sidet = subel->SideToSideTransform(sd,nsides-1);//transf. no elemento metre
				}
				
				//transformacao para o interior do mestre do sub-elemento
				sidet.Apply(point,point2);
				
				//ponto no lado do filho deformado
				subel->X(point2,x1);
				
				//transformacao: espaco paramcorico do filho/lado  -> espaco paramcorico do pai/lado
				father->ComputeXInverse(x1, pointparamfather, 1.e-10);
				
				//no elemento mestre do pai, father coo deformado
				fatside = father->WhichSide(pointparamfather);
				dimfatside = father->SideDimension(fatside);
				
				//transformacoo calculada pelo TPZRefPattern
				TPZTransform trans = refp->Transform(sd,isub);
				
				TPZVec<REAL> pf(dimfatside);
				
				//ponto pf no espaco parametrico do lado do pai
				trans.Apply(point,pf);
				
				//unitcoia do lado do pai
				TPZTransform sidetf(dimfatside);
				if(dimfatside < dimfat)
				{
					//para o interior do sub-elemento
					sidetf = father->SideToSideTransform(fatside,father->NSides()-1);
					
					//sidetf =  father->SideToElemShapeT(fatside);
				}
				
				//do espaco parametrico do lado do pai para o interior do pai
				sidetf.Apply(pf,xpf);
				
				//ponto no lado do pai deformado
				father->X(xpf,x2);
				if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 )
				{
					PZError << "\nTransformacao errada\n";
					PZError << "son    = " << (subel->Id()) << std::endl;
					PZError << "father = " << (father->Id()) << std::endl;
					PZError << "side   = " << sd << std::endl << std::endl;
					
					DebugStop();
				}
				else
				{
					//					std::cout << "Transformacao OK!\n";
					//					std::cout << "Filho/lado : " << subel->Id() << "/" << sd << std::endl;
					//					std::cout << "Pai : " << father->Id() << std::endl << std::endl;
				}//fim if sqrt..
			}//fim rule
			
			delete rule;
		}//fim for sd
	}//fim for isub
}

void TPZRefPatternTools::NodesHunter(TPZGeoMesh &gMesh, TPZVec<int>& NodesHunted, int IdIni, int IdFin, double Tol)
{
    /** Although this method considers coordinates in R3, Coord-Z must be constant for all nodes */
    /** i.e., nodes coordinates belong to XY parallel plane */
    int dim = 3;
	
    int VecSize = gMesh.NodeVec().NElements();
    int posIni = -1;
    int posFin = -1;
	
    /** Changing NodesCoords from Cartesian Notation to Vectorial Notation */
    /** with respect to IniNode coordinate */
    TPZFMatrix VectorialNotation(VecSize,dim);
    for(int pos = 0; pos < VecSize; pos++)
    {
        if(gMesh.NodeVec()[pos].Id() == IdIni) posIni = pos;
        if(gMesh.NodeVec()[pos].Id() == IdFin) posFin = pos;
		
        for(int coord = 0; coord < dim; coord++)
        {
            double val = gMesh.NodeVec()[pos].Coord(coord) - gMesh.NodeVec()[IdIni].Coord(coord);
            VectorialNotation.PutVal(pos,coord,val);
        }
    }
	
#ifdef DEBUG
    if(posIni == -1 || posFin == -1)
    {
        std::cout << "Initial Node index or Final Node index doesn't belong to the given TPZGeoNode TPZVec!\n";
        std::cout << "See NodesHunter method!\n";
		DebugStop();
    }
#endif
	
	// Computing BasisChange Matrix
	// Where NewBase X_axis is defined by IniNode->FinNode orientation
	// and   NewBase Y_axis is perpendicular to X_axis in XY plane counter-clockwise
    TPZFMatrix IfromCntoBase(dim,dim,0.);
    IfromCntoBase(dim-1, dim-1) = 1.;
    double norm = 0.;
	
	//computing X_axis
    for(int coord = 0; coord < dim; coord++)
    {
        double val = VectorialNotation.GetVal(posFin,coord) - VectorialNotation.GetVal(posIni,coord);
        IfromCntoBase.PutVal(0,coord,val);
        norm += val*val;
    }
	//normalizing X_axis
    for(int coord = 0; coord < dim; coord++)
	{
        IfromCntoBase.PutVal(0,coord,IfromCntoBase.GetVal(0,coord)/sqrt(norm));
    }
	//computing Y_axis, i.e., if X_axis=(x,y) then Y_axis=(-y,x)
    IfromCntoBase.PutVal(1,0,-IfromCntoBase.GetVal(0,1));
    IfromCntoBase.PutVal(1,1, IfromCntoBase.GetVal(0,0));
	
    /** Changing VectorialNotation from Canonic Base to NewBase */
    for(int i = 0; i < VecSize; i++)
    {
		TPZManVector <double> temp; temp.Resize(3);
		for(int j = 0; j < dim; j++)
		{
			double val = 0.;
			for(int k = 0; k < dim; k++)
			{
				val += IfromCntoBase.GetVal(j,k)*VectorialNotation.GetVal(i,k);
			}
			temp[j] = val;
		}
		for(int p = 0; p < dim; p++) VectorialNotation.PutVal(i,p,temp[p]);
    }
	
    /** Hunting Nodes */
    std::map <double,int> mymap;
    for(int h = 0; h < VecSize; h++)
    {
        if( VectorialNotation.GetVal(h,0) >= VectorialNotation.GetVal(posIni,0) &&
		   VectorialNotation.GetVal(h,0) <= VectorialNotation.GetVal(posFin,0) &&
		   (fabs(VectorialNotation.GetVal(h,1)) + fabs(VectorialNotation.GetVal(h,2))) < Tol )
		{
			std::pair< double , int> Item(VectorialNotation.GetVal(h,0), gMesh.NodeVec()[h].Id());
			mymap.insert(Item);
		}
    }
    NodesHunted.Resize(mymap.size());
    std::map<double, int>::iterator it = mymap.begin();
    int i = 0;
    for(it = mymap.begin(); it!= mymap.end(); it++, i++)
    {
        NodesHunted[i] = it->second;
    }
}

void TPZRefPatternTools::GetGelPermutations(TPZGeoEl * gel, TPZVec< TPZVec<int> > &permutation)
{
	int id;
	GetElTypePermutations(gel->Type(), permutation);
	for(int p = 0; p < permutation.NElements(); p++)
	{
		for(int n = 0; n < permutation[p].NElements(); n++)
		{
			id = gel->NodeIndex(permutation[p][n]);
			permutation[p][n] = id;
		}
	}
}

void TPZRefPatternTools::GetElTypePermutations(MElementType elType, TPZVec< TPZVec<int> > &permutation)
{
	int nperm, nnodes, p;
	switch (elType)
	{
		case (EOned) :
		{
			nperm = 2;
			nnodes = 2;
			permutation.Resize(nperm);
			for(p = 0; p < nperm; p++)
			{
				permutation[p].Resize(nnodes);
			}
			p = 0;
			permutation[p][0] = 0; permutation[p][1] = 1; p++;
			permutation[p][0] = 1; permutation[p][1] = 0;
			
			break;
		}
		case (ETriangle) :
		{
			nperm = 6;
			nnodes = 3;
			permutation.Resize(nperm);
			for(p = 0; p < nperm; p++)
			{
				permutation[p].Resize(nnodes);
			}
			p = 0;
			permutation[p][0] = 0; permutation[p][1] = 1; permutation[p][2] = 2; p++;
			permutation[p][0] = 1; permutation[p][1] = 2; permutation[p][2] = 0; p++;
			permutation[p][0] = 2; permutation[p][1] = 0; permutation[p][2] = 1; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 2; permutation[p][2] = 1; p++;
			permutation[p][0] = 2; permutation[p][1] = 1; permutation[p][2] = 0; p++;
			permutation[p][0] = 1; permutation[p][1] = 0; permutation[p][2] = 2;
			
			break;
		}
		case (EQuadrilateral) :
		{
			nperm = 8;
			nnodes = 4;
			permutation.Resize(nperm);
			for(p = 0; p < nperm; p++)
			{
				permutation[p].Resize(nnodes);
			}
			p = 0;
			permutation[p][0] = 0; permutation[p][1] = 1; permutation[p][2] = 2; permutation[p][3] = 3; p++;
			permutation[p][0] = 1; permutation[p][1] = 2; permutation[p][2] = 3; permutation[p][3] = 0; p++;
			permutation[p][0] = 2; permutation[p][1] = 3; permutation[p][2] = 0; permutation[p][3] = 1; p++;
			permutation[p][0] = 3; permutation[p][1] = 0; permutation[p][2] = 1; permutation[p][3] = 2; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 3; permutation[p][2] = 2; permutation[p][3] = 1; p++;
			permutation[p][0] = 3; permutation[p][1] = 2; permutation[p][2] = 1; permutation[p][3] = 0; p++;
			permutation[p][0] = 2; permutation[p][1] = 1; permutation[p][2] = 0; permutation[p][3] = 3; p++;
			permutation[p][0] = 1; permutation[p][1] = 0; permutation[p][2] = 3; permutation[p][3] = 2;
			
			break;
		}
		case (ETetraedro) :
		{
			nperm = 24;
			nnodes = 4;
			permutation.Resize(nperm);
			for(p = 0; p < nperm; p++)
			{
				permutation[p].Resize(nnodes);
			}
			p = 0;
			permutation[p][0] = 0; permutation[p][1] = 1; permutation[p][2] = 2;  permutation[p][3] = 3; p++;
			permutation[p][0] = 1; permutation[p][1] = 2; permutation[p][2] = 0;  permutation[p][3] = 3; p++;
			permutation[p][0] = 2; permutation[p][1] = 0; permutation[p][2] = 1;  permutation[p][3] = 3; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 2; permutation[p][2] = 1;  permutation[p][3] = 3; p++;
			permutation[p][0] = 2; permutation[p][1] = 1; permutation[p][2] = 0;  permutation[p][3] = 3; p++;
			permutation[p][0] = 1; permutation[p][1] = 0; permutation[p][2] = 2;  permutation[p][3] = 3; p++;
			//
			permutation[p][0] = 0; permutation[p][1] = 1; permutation[p][2] = 3;  permutation[p][3] = 2; p++;
			permutation[p][0] = 1; permutation[p][1] = 3; permutation[p][2] = 0;  permutation[p][3] = 2; p++;
			permutation[p][0] = 3; permutation[p][1] = 0; permutation[p][2] = 1;  permutation[p][3] = 2; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 3; permutation[p][2] = 1;  permutation[p][3] = 2; p++;
			permutation[p][0] = 3; permutation[p][1] = 1; permutation[p][2] = 0;  permutation[p][3] = 2; p++;
			permutation[p][0] = 1; permutation[p][1] = 0; permutation[p][2] = 3;  permutation[p][3] = 2; p++;
			//
			permutation[p][0] = 0; permutation[p][1] = 3; permutation[p][2] = 2;  permutation[p][3] = 1; p++;
			permutation[p][0] = 3; permutation[p][1] = 2; permutation[p][2] = 0;  permutation[p][3] = 1; p++;
			permutation[p][0] = 2; permutation[p][1] = 0; permutation[p][2] = 3;  permutation[p][3] = 1; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 2; permutation[p][2] = 3;  permutation[p][3] = 1; p++;
			permutation[p][0] = 2; permutation[p][1] = 3; permutation[p][2] = 0;  permutation[p][3] = 1; p++;
			permutation[p][0] = 3; permutation[p][1] = 0; permutation[p][2] = 2;  permutation[p][3] = 1; p++;
			//
			permutation[p][0] = 3; permutation[p][1] = 1; permutation[p][2] = 2;  permutation[p][3] = 0; p++;
			permutation[p][0] = 1; permutation[p][1] = 2; permutation[p][2] = 3;  permutation[p][3] = 0; p++;
			permutation[p][0] = 2; permutation[p][1] = 3; permutation[p][2] = 1;  permutation[p][3] = 0; p++;
			
			permutation[p][0] = 3; permutation[p][1] = 2; permutation[p][2] = 1;  permutation[p][3] = 0; p++;
			permutation[p][0] = 2; permutation[p][1] = 1; permutation[p][2] = 3;  permutation[p][3] = 0; p++;
			permutation[p][0] = 1; permutation[p][1] = 3; permutation[p][2] = 2;  permutation[p][3] = 0;
			
			break;
		}
		case (EPiramide) :
		{
			nperm = 8;
			nnodes = 5;
			permutation.Resize(nperm);
			for(p = 0; p < nperm; p++)
			{
				permutation[p].Resize(nnodes);
			}
			p = 0;
			permutation[p][0] = 0; permutation[p][1] = 1; permutation[p][2] = 2; permutation[p][3] = 3; permutation[p][4] = 4; p++;
			permutation[p][0] = 1; permutation[p][1] = 2; permutation[p][2] = 3; permutation[p][3] = 0; permutation[p][4] = 4; p++;
			permutation[p][0] = 2; permutation[p][1] = 3; permutation[p][2] = 0; permutation[p][3] = 1; permutation[p][4] = 4; p++;
			permutation[p][0] = 3; permutation[p][1] = 0; permutation[p][2] = 1; permutation[p][3] = 2; permutation[p][4] = 4; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 3; permutation[p][2] = 2; permutation[p][3] = 1; permutation[p][4] = 4; p++;
			permutation[p][0] = 3; permutation[p][1] = 2; permutation[p][2] = 1; permutation[p][3] = 0; permutation[p][4] = 4; p++;
			permutation[p][0] = 2; permutation[p][1] = 1; permutation[p][2] = 0; permutation[p][3] = 3; permutation[p][4] = 4; p++;
			permutation[p][0] = 1; permutation[p][1] = 0; permutation[p][2] = 3; permutation[p][3] = 2;
			
			break;
		}
		case (EPrisma) :
		{
			nperm = 12;
			nnodes = 6;
			permutation.Resize(nperm);
			for(p = 0; p < nperm; p++)
			{
				permutation[p].Resize(nnodes);
			}
			p = 0;
			permutation[p][0] = 0; permutation[p][1] = 1; permutation[p][2] = 2; permutation[p][3] = 3; permutation[p][4] = 4; permutation[p][5] = 5; p++;
			permutation[p][0] = 1; permutation[p][1] = 2; permutation[p][2] = 0; permutation[p][3] = 4; permutation[p][4] = 5; permutation[p][5] = 3; p++;
			permutation[p][0] = 2; permutation[p][1] = 0; permutation[p][2] = 1; permutation[p][3] = 5; permutation[p][4] = 3; permutation[p][5] = 4; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 2; permutation[p][2] = 1; permutation[p][3] = 3; permutation[p][4] = 5; permutation[p][5] = 4; p++;
			permutation[p][0] = 2; permutation[p][1] = 1; permutation[p][2] = 0; permutation[p][3] = 5; permutation[p][4] = 4; permutation[p][5] = 3; p++;
			permutation[p][0] = 1; permutation[p][1] = 0; permutation[p][2] = 2; permutation[p][3] = 4; permutation[p][4] = 3; permutation[p][5] = 5; p++;
			//
			permutation[p][0] = 3; permutation[p][1] = 4; permutation[p][2] = 5; permutation[p][3] = 0; permutation[p][4] = 1; permutation[p][5] = 2; p++;
			permutation[p][0] = 4; permutation[p][1] = 5; permutation[p][2] = 3; permutation[p][3] = 1; permutation[p][4] = 2; permutation[p][5] = 0; p++;
			permutation[p][0] = 5; permutation[p][1] = 3; permutation[p][2] = 4; permutation[p][3] = 2; permutation[p][4] = 0; permutation[p][5] = 1; p++;
			
			permutation[p][0] = 3; permutation[p][1] = 5; permutation[p][2] = 4; permutation[p][3] = 0; permutation[p][4] = 2; permutation[p][5] = 1; p++;
			permutation[p][0] = 5; permutation[p][1] = 4; permutation[p][2] = 3; permutation[p][3] = 2; permutation[p][4] = 1; permutation[p][5] = 0; p++;
			permutation[p][0] = 4; permutation[p][1] = 3; permutation[p][2] = 5; permutation[p][3] = 1; permutation[p][4] = 0; permutation[p][5] = 2;
			
			break;
		}
		case (ECube) :
		{
			nperm = 48;
			nnodes = 8;
			permutation.Resize(nperm);
			for(p = 0; p < nperm; p++)
			{
				permutation[p].Resize(nnodes);
			}
			p = 0;
			permutation[p][0] = 0; permutation[p][1] = 1; permutation[p][2] = 2; permutation[p][3] = 3; permutation[p][4] = 4; permutation[p][5] = 5; permutation[p][6] = 6; permutation[p][7] = 7; p++;
			permutation[p][0] = 1; permutation[p][1] = 2; permutation[p][2] = 3; permutation[p][3] = 0; permutation[p][4] = 5; permutation[p][5] = 6; permutation[p][6] = 7; permutation[p][7] = 4; p++;
			permutation[p][0] = 2; permutation[p][1] = 3; permutation[p][2] = 0; permutation[p][3] = 1; permutation[p][4] = 6; permutation[p][5] = 7; permutation[p][6] = 4; permutation[p][7] = 5; p++;
			permutation[p][0] = 3; permutation[p][1] = 0; permutation[p][2] = 1; permutation[p][3] = 2; permutation[p][4] = 7; permutation[p][5] = 4; permutation[p][6] = 5; permutation[p][7] = 6; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 3; permutation[p][2] = 2; permutation[p][3] = 1; permutation[p][4] = 4; permutation[p][5] = 7; permutation[p][6] = 6; permutation[p][7] = 5; p++;
			permutation[p][0] = 3; permutation[p][1] = 2; permutation[p][2] = 1; permutation[p][3] = 0; permutation[p][4] = 7; permutation[p][5] = 6; permutation[p][6] = 5; permutation[p][7] = 4; p++;
			permutation[p][0] = 2; permutation[p][1] = 1; permutation[p][2] = 0; permutation[p][3] = 3; permutation[p][4] = 6; permutation[p][5] = 5; permutation[p][6] = 4; permutation[p][7] = 7; p++;
			permutation[p][0] = 1; permutation[p][1] = 0; permutation[p][2] = 3; permutation[p][3] = 2; permutation[p][4] = 5; permutation[p][5] = 4; permutation[p][6] = 7; permutation[p][7] = 6; p++;
			//
			permutation[p][0] = 0; permutation[p][1] = 1; permutation[p][2] = 5; permutation[p][3] = 4; permutation[p][4] = 3; permutation[p][5] = 2; permutation[p][6] = 6; permutation[p][7] = 7; p++;
			permutation[p][0] = 1; permutation[p][1] = 5; permutation[p][2] = 4; permutation[p][3] = 0; permutation[p][4] = 2; permutation[p][5] = 6; permutation[p][6] = 7; permutation[p][7] = 3; p++;
			permutation[p][0] = 5; permutation[p][1] = 4; permutation[p][2] = 0; permutation[p][3] = 1; permutation[p][4] = 6; permutation[p][5] = 7; permutation[p][6] = 3; permutation[p][7] = 2; p++;
			permutation[p][0] = 4; permutation[p][1] = 0; permutation[p][2] = 1; permutation[p][3] = 5; permutation[p][4] = 7; permutation[p][5] = 3; permutation[p][6] = 2; permutation[p][7] = 6; p++;
			
			permutation[p][0] = 0; permutation[p][1] = 4; permutation[p][2] = 5; permutation[p][3] = 1; permutation[p][4] = 3; permutation[p][5] = 7; permutation[p][6] = 6; permutation[p][7] = 2; p++;
			permutation[p][0] = 4; permutation[p][1] = 5; permutation[p][2] = 1; permutation[p][3] = 0; permutation[p][4] = 7; permutation[p][5] = 6; permutation[p][6] = 2; permutation[p][7] = 3; p++;
			permutation[p][0] = 5; permutation[p][1] = 1; permutation[p][2] = 0; permutation[p][3] = 4; permutation[p][4] = 6; permutation[p][5] = 2; permutation[p][6] = 3; permutation[p][7] = 7; p++;
			permutation[p][0] = 1; permutation[p][1] = 0; permutation[p][2] = 4; permutation[p][3] = 5; permutation[p][4] = 2; permutation[p][5] = 3; permutation[p][6] = 7; permutation[p][7] = 6; p++;
			//
			permutation[p][0] = 1; permutation[p][1] = 2; permutation[p][2] = 6; permutation[p][3] = 5; permutation[p][4] = 0; permutation[p][5] = 3; permutation[p][6] = 7; permutation[p][7] = 4; p++;
			permutation[p][0] = 2; permutation[p][1] = 6; permutation[p][2] = 5; permutation[p][3] = 1; permutation[p][4] = 3; permutation[p][5] = 7; permutation[p][6] = 4; permutation[p][7] = 0; p++;
			permutation[p][0] = 6; permutation[p][1] = 5; permutation[p][2] = 1; permutation[p][3] = 2; permutation[p][4] = 7; permutation[p][5] = 4; permutation[p][6] = 0; permutation[p][7] = 3; p++;
			permutation[p][0] = 5; permutation[p][1] = 1; permutation[p][2] = 2; permutation[p][3] = 6; permutation[p][4] = 4; permutation[p][5] = 0; permutation[p][6] = 3; permutation[p][7] = 7; p++;
			
			permutation[p][0] = 1; permutation[p][1] = 5; permutation[p][2] = 6; permutation[p][3] = 2; permutation[p][4] = 0; permutation[p][5] = 4; permutation[p][6] = 7; permutation[p][7] = 3; p++;
			permutation[p][0] = 5; permutation[p][1] = 6; permutation[p][2] = 2; permutation[p][3] = 1; permutation[p][4] = 4; permutation[p][5] = 7; permutation[p][6] = 3; permutation[p][7] = 0; p++;
			permutation[p][0] = 6; permutation[p][1] = 2; permutation[p][2] = 1; permutation[p][3] = 5; permutation[p][4] = 7; permutation[p][5] = 3; permutation[p][6] = 0; permutation[p][7] = 4; p++;
			permutation[p][0] = 2; permutation[p][1] = 1; permutation[p][2] = 5; permutation[p][3] = 6; permutation[p][4] = 3; permutation[p][5] = 0; permutation[p][6] = 4; permutation[p][7] = 7; p++;
			//
			permutation[p][0] = 2; permutation[p][1] = 3; permutation[p][2] = 7; permutation[p][3] = 6; permutation[p][4] = 1; permutation[p][5] = 0; permutation[p][6] = 4; permutation[p][7] = 5; p++;
			permutation[p][0] = 3; permutation[p][1] = 7; permutation[p][2] = 6; permutation[p][3] = 2; permutation[p][4] = 0; permutation[p][5] = 4; permutation[p][6] = 5; permutation[p][7] = 1; p++;
			permutation[p][0] = 7; permutation[p][1] = 6; permutation[p][2] = 2; permutation[p][3] = 3; permutation[p][4] = 4; permutation[p][5] = 5; permutation[p][6] = 1; permutation[p][7] = 0; p++;
			permutation[p][0] = 6; permutation[p][1] = 2; permutation[p][2] = 3; permutation[p][3] = 7; permutation[p][4] = 5; permutation[p][5] = 1; permutation[p][6] = 0; permutation[p][7] = 4; p++;
			
			permutation[p][0] = 2; permutation[p][1] = 6; permutation[p][2] = 7; permutation[p][3] = 3; permutation[p][4] = 1; permutation[p][5] = 5; permutation[p][6] = 4; permutation[p][7] = 0; p++;
			permutation[p][0] = 6; permutation[p][1] = 7; permutation[p][2] = 3; permutation[p][3] = 2; permutation[p][4] = 5; permutation[p][5] = 4; permutation[p][6] = 0; permutation[p][7] = 1; p++;
			permutation[p][0] = 7; permutation[p][1] = 3; permutation[p][2] = 2; permutation[p][3] = 6; permutation[p][4] = 4; permutation[p][5] = 0; permutation[p][6] = 1; permutation[p][7] = 5; p++;
			permutation[p][0] = 3; permutation[p][1] = 2; permutation[p][2] = 6; permutation[p][3] = 7; permutation[p][4] = 0; permutation[p][5] = 1; permutation[p][6] = 5; permutation[p][7] = 4; p++;			
			//
			permutation[p][0] = 3; permutation[p][1] = 0; permutation[p][2] = 4; permutation[p][3] = 7; permutation[p][4] = 1; permutation[p][5] = 5; permutation[p][6] = 6; permutation[p][7] = 4; p++;
			permutation[p][0] = 0; permutation[p][1] = 4; permutation[p][2] = 7; permutation[p][3] = 3; permutation[p][4] = 5; permutation[p][5] = 6; permutation[p][6] = 4; permutation[p][7] = 1; p++;
			permutation[p][0] = 4; permutation[p][1] = 7; permutation[p][2] = 3; permutation[p][3] = 0; permutation[p][4] = 6; permutation[p][5] = 4; permutation[p][6] = 1; permutation[p][7] = 5; p++;
			permutation[p][0] = 7; permutation[p][1] = 3; permutation[p][2] = 0; permutation[p][3] = 4; permutation[p][4] = 4; permutation[p][5] = 1; permutation[p][6] = 5; permutation[p][7] = 6; p++;
			
			permutation[p][0] = 3; permutation[p][1] = 7; permutation[p][2] = 4; permutation[p][3] = 0; permutation[p][4] = 2; permutation[p][5] = 6; permutation[p][6] = 5; permutation[p][7] = 1; p++;
			permutation[p][0] = 7; permutation[p][1] = 4; permutation[p][2] = 0; permutation[p][3] = 3; permutation[p][4] = 6; permutation[p][5] = 5; permutation[p][6] = 1; permutation[p][7] = 2; p++;
			permutation[p][0] = 4; permutation[p][1] = 0; permutation[p][2] = 3; permutation[p][3] = 7; permutation[p][4] = 5; permutation[p][5] = 1; permutation[p][6] = 2; permutation[p][7] = 6; p++;
			permutation[p][0] = 0; permutation[p][1] = 3; permutation[p][2] = 7; permutation[p][3] = 4; permutation[p][4] = 1; permutation[p][5] = 2; permutation[p][6] = 6; permutation[p][7] = 5; p++;
			//
			permutation[p][0] = 4; permutation[p][1] = 5; permutation[p][2] = 6; permutation[p][3] = 7; permutation[p][4] = 0; permutation[p][5] = 1; permutation[p][6] = 2; permutation[p][7] = 3; p++;
			permutation[p][0] = 5; permutation[p][1] = 6; permutation[p][2] = 7; permutation[p][3] = 4; permutation[p][4] = 1; permutation[p][5] = 2; permutation[p][6] = 3; permutation[p][7] = 0; p++;
			permutation[p][0] = 6; permutation[p][1] = 7; permutation[p][2] = 4; permutation[p][3] = 5; permutation[p][4] = 2; permutation[p][5] = 3; permutation[p][6] = 0; permutation[p][7] = 1; p++;
			permutation[p][0] = 7; permutation[p][1] = 4; permutation[p][2] = 5; permutation[p][3] = 6; permutation[p][4] = 3; permutation[p][5] = 0; permutation[p][6] = 1; permutation[p][7] = 2; p++;
			
			permutation[p][0] = 4; permutation[p][1] = 7; permutation[p][2] = 6; permutation[p][3] = 5; permutation[p][4] = 0; permutation[p][5] = 3; permutation[p][6] = 2; permutation[p][7] = 1; p++;
			permutation[p][0] = 7; permutation[p][1] = 6; permutation[p][2] = 5; permutation[p][3] = 4; permutation[p][4] = 3; permutation[p][5] = 2; permutation[p][6] = 1; permutation[p][7] = 0; p++;
			permutation[p][0] = 6; permutation[p][1] = 5; permutation[p][2] = 4; permutation[p][3] = 7; permutation[p][4] = 2; permutation[p][5] = 1; permutation[p][6] = 0; permutation[p][7] = 3; p++;
			permutation[p][0] = 5; permutation[p][1] = 4; permutation[p][2] = 7; permutation[p][3] = 6; permutation[p][4] = 1; permutation[p][5] = 0; permutation[p][6] = 3; permutation[p][7] = 2; p++;
			
			break;
		}
		default:
		{
			cout << "Cant return permutation because MElementType was not found on " << __PRETTY_FUNCTION__ << endl;
			DebugStop();
		}
	}
}




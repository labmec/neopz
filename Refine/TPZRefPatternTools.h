#ifndef TPZREFPATTERNTOOLSH
#define TPZREFPATTERNTOOLSH 

/*
 *  TPZRefPatternTools.h
 *  Crack
 *
 *  Created by Cesar Lucci on 10/03/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

using namespace std;

#include <set>

#include "pzgeoel.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternTools.h"
#include "tpzintpoints.h"

class TPZRefPatternTools
{
	
public:
	
	TPZRefPatternTools();
  ~TPZRefPatternTools();
	
	/**
	 * Search for refpatterns that could be used by a given element with respect to their neighbours.
	 * @param gel - input data: geometric element for which the list of compatible refpatterns will be filled
	 * @param refs - output data: list of compatible refpatterns with respect to their neighbours
	 */
	static void GetCompatibleRefPatterns(TPZGeoEl *gel, std::list<TPZAutoPointer<TPZRefPattern> > &refs);
	
	/**
	 * Returns the refpattern that matches the sides refinement by neighbours
	 * @param gel - input data: geometric element for which the model refpattern will be returned
	 * @param neighCorresp - output data: map that group (neighbour geoelement) and nodes correspondences between (neighbour->SideRefpattern) and (gel->SideRefpattern), indexed by (gel->Side)
	 * IF THERE IS NO NEIGHBOUR ALREADY REFINED, IT RETURNS NULL
	 */
	static TPZAutoPointer<TPZRefPattern> ModelRefPattern(TPZGeoEl *gel, std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > > &neighCorresp);
	
	/**
	 * This methos is used by RefineDirectional method!!!
	 * Returns the refpattern that matches the sides refinement intensity and midnodes coordinates with respect to sidestorefine vector
	 * @param gel - input data: geometric element for which the perfect match refpattern will be returned
	 * @param sidestorefine - input data: vector filled with sides refinement intensity
	 */
	static TPZAutoPointer<TPZRefPattern> PerfectMatchRefPattern(TPZGeoEl *gel, TPZVec<int> &sidestorefine);
	
	/**
	 * Returns the refpattern that matches the sides refinement intensity and midnodes coordinates with respect to sideNeighbours
	 * @param gel - input data: geometric element for which the perfect match refpatterns will be returned
	 * IF THERE IS NO NEIGHBOUR ALREADY REFINED, IT RETURNS NULL
	 */
	static TPZAutoPointer<TPZRefPattern> PerfectMatchRefPattern(TPZGeoEl *gel);
	
	/**
	 * For a given refpattern model, its midnodes will be dragged to match with a geoel neighbourhood refinement
	 * @param gel - input data: geometric element for which the model refpattern nodes will be dragged
	 * @param modelPat - input data: Model RefPattern that is topologicaly compatible with neighbourhood
	 * @param neighCorresp - input data: map that group (neighbour geoelement) and nodes correspondences between (neighbour->SideRefpattern) and (gel->SideRefpattern), indexed by (gel->Side)
	 */
	static void DragModelPatNodes(TPZGeoEl * gel, TPZAutoPointer<TPZRefPattern> modelPat, std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > > &neighCorresp);
	
	/**
	 * Returns if the given refPatterns (refA and refB) are topologicaly compatibles.
	 * If they are, pairNodes represents the correspondence between nodesIds from refAmesh to refBmesh.
	 * @param refA - input data: first refpattern to be compared
	 * @param refB - input data: second refpattern to be compared
	 * @param fromAtoB - input data: Linear Transformation from refA->RefPatternMesh to refB->RefPatternMesh() (its needed to pair Nodes correctly in permuted cases)
	 * @param pairNodes - output data: correspondence between refA and refB nodes, in case they are topologicaly compatibles
	 */
	static bool CompareTopologies(TPZAutoPointer<TPZRefPattern> refA, TPZAutoPointer<TPZRefPattern> refB, TPZTransform &fromAtoB, std::map<int, int> &pairNodes);
	
	/**
	 *	This method pair nodes from meshA to meshB, using the givem transformation from meshA to meshB to match coordinates.
	 * The output is the map pairNodes, thar represents the A_nodeId paired with B_nodeId.
	 * Obs.: Be careful with the output interpretation! It contains the nodeIds, NOT the nodes positions in mesh.NodeVec()!!!
	 * @param meshA - input data: first mesh that will be considered its nodes coordinates
	 * @param meshB - input data: second mesh that will be considered its nodes coordinates
	 * @param fromAtoB - input data: Linear Transformation from meshA to meshB (its needed to pair Nodes correctly in permuted cases)
	 * @param pairedNodes - output data: correspondence between meshA and meshB nodes that matches its coordinates
	 */
	static void PairMeshesNodesMatchingCoordinates(TPZGeoMesh &meshA, TPZGeoMesh &meshB, TPZTransform fromAtoB, std::map<int, int> &pairedNodes);
	
	/**
	 * Returns the the name of refpattern model.
	 * To do this, it starts with the 3 initial characters of element nametype,
	 * followed by the quantity of midnodes for each side of element.
	 */
	static std::string BuildRefPatternModelName(TPZRefPattern &refp);
	static std::string BuildRefPatternModelName(TPZAutoPointer<TPZRefPattern> refp);
	static std::string BuildRefPatternModelName(TPZGeoEl *gel);
	
	/**
	 * Returns if there is any neigbour already refined
	 * @param gel - input data: geometric element whose refinements of the neighbors will define the refinement of its sides
	 * @sidestoref - output data: vector whose positions mention the sides of gel, and its contents mention the intensity of refinement of the respective side
	 */
	static bool SidesToRefine(TPZGeoEl *gel, TPZVec<int> &sidestoref);
	
	/**
     * Refine the element if it touches an element with a material id included in matids
     */
	static void RefineDirectional(TPZGeoEl *gel, std::set<int> &matids);
	static void RefineDirectional(TPZGeoEl *gel, std::set<int> &matids, int gelMat);
	
	/**
	 * Method to test if the jacobian of a TPZGeoElSide element is constant
	 */
	static bool ConstJacobian(TPZGeoElSide gelside, REAL tol = 1.e-6);
	
	/**
	 * Algorithm that evaluates the veracity of the hashings between sides
     * of the elements children and corresponding sides of the father. A
     * point p in the parametric space of the side of the sub-element is
     * overcome and is calculated it mentioned hashing getting pf point in
     * the element father. One calculates for p and pf the corresponding
     * deformed point. Itself the hashing is consistent the deformed point
     * must the same be.
     */
	static void TransformationTest(TPZRefPattern * refp);
	
	/**
	 * NodesHunted vector is the sequential nodesIds that belongs (i.e.: "Tol" far) to InitialNode(IdIni)~FinalNode(IdFin) alignment of gMesh.NodeVec()
	 * Obs.: InitialNode and FinalNode are also included!!!
	 */
	static void NodesHunter(TPZGeoMesh &gMesh, std::vector<int>& NodesHunted, int IdIni, int IdFin, double Tol = 1.E-1);

	/**
	 * Fill the TPZVec "permutation" with the valid permutations of "gel"
	 * Note: The permutations is with respect to Master Element nodes, NOT Gel nodes in a geomesh context (i.e.: NOT geomesh nodes ids)
	 */
	static void GetGelPermutations(TPZGeoEl * gel, TPZVec< TPZVec<int> > &permutation);
	
	/**
	 * Fill the TPZVec "permutation" with the valid permutations of a given element type
	 */
	static void GetElTypePermutations(MElementType elType, TPZVec< TPZVec<int> > &permutation);
	
	/**
	 * Generate an output of all geomesh to VTK
	 */
	static void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, bool matColor = false)
	{
		file.clear();
		int nelements = gmesh->NElements();
		
		stringstream node, connectivity, type, material;
		
		//Header
		file << "# vtk DataFile Version 3.0" << std::endl;
		file << "TPZGeoMesh VTK Visualization" << std::endl;
		file << "ASCII" << std::endl << std::endl;
		
		file << "DATASET UNSTRUCTURED_GRID" << std::endl;
		file << "POINTS ";
		
		int actualNode = -1, size = 0, nVALIDelements = 0;
		
		for(int el = 0; el < nelements; el++)
		{				
			if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines, Arc3D and Ellipse3D
			{
				continue;
			}
			if(gmesh->ElementVec()[el]->HasSubElement())
			{
				continue;
			}
			
			int elNnodes = gmesh->ElementVec()[el]->NNodes();
			size += (1+elNnodes);
			connectivity << elNnodes;
			
			for(int t = 0; t < elNnodes; t++)
			{
				for(int c = 0; c < 3; c++)
				{
					double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
					node << coord << " ";
				}			
				node << std::endl;
				
				actualNode++;
				connectivity << " " << actualNode;
			}
			connectivity << std::endl;
			
			int elType;
			switch (gmesh->ElementVec()[el]->Type())
			{
				case (ETriangle):
				{
					elType = 5;
					break;				
				}
				case (EQuadrilateral):
				{
					elType = 9;
					break;				
				}
				case (ETetraedro):
				{
					elType = 10;
					break;				
				}
				case (EPiramide):
				{
					elType = 14;
					break;				
				}
				case (EPrisma):
				{
					elType = 13;
					break;				
				}
				case (ECube):
				{
					elType = 12;
					break;				
				}
				default:
				{
					elType = -1;//ElementType NOT Found!!!
					DebugStop();
					break;	
				}
			}
			type << elType << std::endl;
			if(matColor == true)
			{
				material << gmesh->ElementVec()[el]->MaterialId() << std::endl;
			}
			else
			{
				material << elType << std::endl;
			}

			nVALIDelements++;
		}
		node << std::endl;
		actualNode++;
		file << actualNode << " float" << std::endl << node.str();
		
		file << "CELLS " << nVALIDelements << " ";
		
		file << size << std::endl;
		file << connectivity.str() << std::endl;
		
		file << "CELL_TYPES " << nVALIDelements << std::endl;
		file << type.str() << std::endl;
		
		file << "CELL_DATA" << " " << nVALIDelements << std::endl;
		file << "FIELD FieldData 1" << std::endl;
		if(matColor == true)
		{
			file << "material 1 " << nVALIDelements << " int" << std::endl;
		}
		else
		{
			file << "ElementType 1 " << nVALIDelements << " int" << std::endl;
		}
		file << material.str();
		
		file.close();
	}
	
	/**
	 * Generate an output of all geomesh to VTK, associating to each one the given data
	 */
	static void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, TPZVec<int> &elData)
	{
		if(gmesh->NElements() != elData.NElements())
		{
			std::cout << "Wrong vector size of elements data!" << std::endl;
			std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		}
		file.clear();
		int nelements = gmesh->NElements();
		
		stringstream node, connectivity, type, material;
		
		//Header
		file << "# vtk DataFile Version 3.0" << std::endl;
		file << "TPZGeoMesh VTK Visualization" << std::endl;
		file << "ASCII" << std::endl << std::endl;
		
		file << "DATASET UNSTRUCTURED_GRID" << std::endl;
		file << "POINTS ";
		
		int actualNode = -1, size = 0, nVALIDelements = 0;
		
		for(int el = 0; el < nelements; el++)
		{				
			if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines, Arc3D and Ellipse3D
			{
				continue;
			}
			if(gmesh->ElementVec()[el]->HasSubElement())
			{
				continue;
			}
			
			int elNnodes = gmesh->ElementVec()[el]->NNodes();
			size += (1+elNnodes);
			connectivity << elNnodes;
			
			for(int t = 0; t < elNnodes; t++)
			{
				for(int c = 0; c < 3; c++)
				{
					double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
					node << coord << " ";
				}			
				node << std::endl;
				
				actualNode++;
				connectivity << " " << actualNode;
			}
			connectivity << std::endl;
			
			int elType;
			switch (gmesh->ElementVec()[el]->Type())
			{
				case (ETriangle):
				{
					elType = 5;
					break;				
				}
				case (EQuadrilateral):
				{
					elType = 9;
					break;				
				}
				case (ETetraedro):
				{
					elType = 10;
					break;				
				}
				case (EPiramide):
				{
					elType = 14;
					break;				
				}
				case (EPrisma):
				{
					elType = 13;
					break;				
				}
				case (ECube):
				{
					elType = 12;
					break;				
				}
				default:
				{
					elType = -1;//ElementType NOT Found!!!
					DebugStop();
					break;	
				}
			}
			type << elType << std::endl;
			
			material << elData[el] << std::endl;
			
			nVALIDelements++;
		}
		node << std::endl;
		actualNode++;
		file << actualNode << " float" << std::endl << node.str();
		
		file << "CELLS " << nVALIDelements << " ";
		
		file << size << std::endl;
		file << connectivity.str() << std::endl;
		
		file << "CELL_TYPES " << nVALIDelements << std::endl;
		file << type.str() << std::endl;
		
		file << "CELL_DATA" << " " << nVALIDelements << std::endl;
		file << "FIELD FieldData 1" << std::endl;

		file << "Substructure 1 " << nVALIDelements << " int" << std::endl;

		file << material.str();
		
		file.close();
	}
	
	/**
	 * Based on a given geomesh, just the elements that have an neighbour with a given material id will be exported to an VTK file
	 */
	static void PrintGMeshVTKneighbour_material(TPZGeoMesh * gmesh, std::ofstream &file, int neighMaterial, bool matColor = false)
	{
		file.clear();
		int nelements = gmesh->NElements();
		
		stringstream node, connectivity, type, material;
		
		//Header
		file << "# vtk DataFile Version 3.0" << std::endl;
		file << "TPZGeoMesh VTK Visualization" << std::endl;
		file << "ASCII" << std::endl << std::endl;
		
		file << "DATASET UNSTRUCTURED_GRID" << std::endl;
		file << "POINTS ";
		
		int actualNode = -1, size = 0, nVALIDelements = 0;
		
		for(int el = 0; el < nelements; el++)
		{				
			if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines, Arc3D and Ellipse3D
			{
				continue;
			}
			if(gmesh->ElementVec()[el]->HasSubElement())
			{
				continue;
			}
			
			bool matFound = false;
			for(int s = 0; s < gmesh->ElementVec()[el]->NSides(); s++)
			{
				TPZGeoElSide thisSide(gmesh->ElementVec()[el], s);
				TPZGeoElSide neighSide = thisSide.Neighbour();
				
				while(thisSide != neighSide)
				{
					if(neighSide.Element()->MaterialId() == neighMaterial)
					{
						matFound = true;
						break;
					}
					neighSide = neighSide.Neighbour();
				}
				if(matFound)
				{
					break;
				}
			}
			if(!matFound)
			{
				continue;
			}
			
			int elNnodes = gmesh->ElementVec()[el]->NNodes();
			size += (1+elNnodes);
			connectivity << elNnodes;
			
			for(int t = 0; t < elNnodes; t++)
			{
				for(int c = 0; c < 3; c++)
				{
					double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
					node << coord << " ";
				}			
				node << std::endl;
				
				actualNode++;
				connectivity << " " << actualNode;
			}
			connectivity << std::endl;
			
			int elType;
			switch (gmesh->ElementVec()[el]->Type())
			{
				case (ETriangle):
				{
					elType = 5;
					break;				
				}
				case (EQuadrilateral):
				{
					elType = 9;
					break;				
				}
				case (ETetraedro):
				{
					elType = 10;
					break;				
				}
				case (EPiramide):
				{
					elType = 14;
					break;				
				}
				case (EPrisma):
				{
					elType = 13;
					break;				
				}
				case (ECube):
				{
					elType = 12;
					break;				
				}
				default:
				{
					elType = -1;//ElementType NOT Found!!!
					DebugStop();
					break;	
				}
			}
			type << elType << std::endl;
			if(matColor == true)
			{
				material << gmesh->ElementVec()[el]->MaterialId() << std::endl;
			}
			else
			{
				material << elType << std::endl;
			}
			
			nVALIDelements++;
		}
		node << std::endl;
		actualNode++;
		file << actualNode << " float" << std::endl << node.str();
		
		file << "CELLS " << nVALIDelements << " ";
		
		file << size << std::endl;
		file << connectivity.str() << std::endl;
		
		file << "CELL_TYPES " << nVALIDelements << std::endl;
		file << type.str();
		
		file << "CELL_DATA" << " " << nVALIDelements << std::endl;
		file << "FIELD FieldData 1" << std::endl;
		if(matColor == true)
		{
			file << "material 1 " << nVALIDelements << " int" << std::endl;
		}
		else
		{
			file << "ElementType 1 " << nVALIDelements << " int" << std::endl;
		}
		file << material.str();
		
		file.close();
	}
	
	/**
	 * Based on a given geomesh, just the elements that have the given material id will be exported to an VTK file
	 */
	static void PrintGMeshVTKmy_material(TPZGeoMesh * gmesh, std::ofstream &file, std::set<int> myMaterial, bool matColor = false)
	{
		file.clear();
		int nelements = gmesh->NElements();
		
		stringstream node, connectivity, type, material;
		
		//Header
		file << "# vtk DataFile Version 3.0" << std::endl;
		file << "TPZGeoMesh VTK Visualization" << std::endl;
		file << "ASCII" << std::endl << std::endl;
		
		file << "DATASET UNSTRUCTURED_GRID" << std::endl;
		file << "POINTS ";
		
		int actualNode = -1, size = 0, nVALIDelements = 0;
		
		for(int el = 0; el < nelements; el++)
		{				
			if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines, Arc3D and Ellipse3D
			{
				continue;
			}
			int mat = gmesh->ElementVec()[el]->MaterialId();
			bool found = !(myMaterial.find(mat) == myMaterial.end() );
			if(gmesh->ElementVec()[el]->HasSubElement() || !found)
			{
				continue;
			}
			
			int elNnodes = gmesh->ElementVec()[el]->NNodes();
			size += (1+elNnodes);
			connectivity << elNnodes;
			
			for(int t = 0; t < elNnodes; t++)
			{
				for(int c = 0; c < 3; c++)
				{
					double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
					node << coord << " ";
				}			
				node << std::endl;
				
				actualNode++;
				connectivity << " " << actualNode;
			}
			connectivity << std::endl;
			
			int elType;
			switch (gmesh->ElementVec()[el]->Type())
			{
				case (ETriangle):
				{
					elType = 5;
					break;				
				}
				case (EQuadrilateral):
				{
					elType = 9;
					break;				
				}
				case (ETetraedro):
				{
					elType = 10;
					break;				
				}
				case (EPiramide):
				{
					elType = 14;
					break;				
				}
				case (EPrisma):
				{
					elType = 13;
					break;				
				}
				case (ECube):
				{
					elType = 12;
					break;				
				}
				default:
				{
					elType = -1;//ElementType NOT Found!!!
					DebugStop();
					break;	
				}
			}
			type << elType << std::endl;
			if(matColor == true)
			{
				material << gmesh->ElementVec()[el]->MaterialId() << std::endl;
			}
			else
			{
				material << elType << std::endl;
			}
			nVALIDelements++;
		}
		node << std::endl;
		actualNode++;
		file << actualNode << " float" << std::endl << node.str();
		
		file << "CELLS " << nVALIDelements << " ";
		
		file << size << std::endl;
		file << connectivity.str() << std::endl;
		
		file << "CELL_TYPES " << nVALIDelements << std::endl;
		file << type.str();
		
		file << "CELL_DATA" << " " << nVALIDelements << std::endl;
		file << "FIELD FieldData 1" << std::endl;
		if(matColor == true)
		{
			file << "material 1 " << nVALIDelements << " int" << std::endl;
		}
		else
		{
			file << "ElementType 1 " << nVALIDelements << " int" << std::endl;
		}
		file << material.str();
		
		file.close();
	}
};

#endif

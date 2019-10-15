/**
 * @file
 * @brief Contains the TPZRefPatternTools class which defines tools of pattern.
 */

#ifndef TPZREFPATTERNTOOLSH
#define TPZREFPATTERNTOOLSH 

/*
 *  TPZRefPatternTools.h
 *  Crack
 *
 *  Created by Cesar Lucci on 10/03/10.
 *  Copyright 2010 LabMeC. All rights reserved.
 */


#include <set>

#include "pzgeoel.h"
#include "tpzintpoints.h"

class TPZRefPattern;
/**
 * @ingroup refine
 * @brief Defines tools of pattern. \ref refine "Refine"
 */
class TPZRefPatternTools
{
	
public:
	
	TPZRefPatternTools();
	~TPZRefPatternTools();
	
	/**
	 * @brief Search for refpatterns that could be used by a given element with respect to their neighbours.
	 * @param gel - input data: geometric element for which the list of compatible refpatterns will be filled
	 * @param refs - output data: list of compatible refpatterns with respect to their neighbours
	 */
	static void GetCompatibleRefPatterns(TPZGeoEl *gel,
                                         std::list<TPZAutoPointer<TPZRefPattern> > &refs);
	
	/**
	 * @brief Returns the refpattern that matches the sides refinement by neighbours
	 * @param gel - input data: geometric element for which the model refpattern will be returned
	 * @param neighCorresp - output data:	map that group (neighbour geoelement) \n
	 *									and nodes correspondences between (neighbour->SideRefpattern) and (gel->SideRefpattern),
	 *									indexed by (gel->Side)
	 * @note IF THERE IS NO NEIGHBOUR ALREADY REFINED, IT RETURNS NULL
	 */
	static TPZAutoPointer<TPZRefPattern> ModelRefPattern(TPZGeoEl *gel,
                                                         std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > > &neighCorresp);
	
	/**
	 * @note This methos is used by RefineDirectional method!!!
	 * @brief Returns the refpattern that matches the sides refinement intensity and midnodes coordinates with respect to sidestorefine vector
	 * @param gel - input data: geometric element for which the perfect match refpattern will be returned
	 * @param sidestorefine - input data: vector filled with sides refinement intensity
	 */
	static TPZAutoPointer<TPZRefPattern> PerfectMatchRefPattern(TPZGeoEl *gel,
                                                                TPZVec<int> &sidestorefine);
	
	/**
	 * @brief Returns the refpattern that matches the sides refinement intensity and midnodes coordinates with respect to sideNeighbours
	 * @param gel - input data: geometric element for which the perfect match refpatterns will be returned
	 * @note IF THERE IS NO NEIGHBOUR ALREADY REFINED, IT RETURNS NULL
	 */
	static TPZAutoPointer<TPZRefPattern> PerfectMatchRefPattern(TPZGeoEl *gel);
    
    /**
     * This method receives an vector of real mesh elements (index0=FATHER , index1...=SONS) and returns the corresponding refinement pattern.
     * ONCE THE REFPATTERN COULD PRESENT PERMUTATIONS OF SONS, THIS METHOD FIX THE TOPOLOGICAL SEQUENCE OF REAL MESH ELEMENTS.
     * This implies the call BuildConnectivity() method after this whole operation (which should be done by the programmer).
     */
    static TPZAutoPointer<TPZRefPattern> GetRefPatternBasedOnRealMeshElements(TPZVec<TPZGeoEl *> & realMeshElementVec);

    
    static void GenerateGMeshFromElementVec(const TPZVec<TPZGeoEl *> & elementVec, TPZGeoMesh & refGMesh);    
	
    /**
     * A partir de um padrao encontrado, corrige a sequencia topologica dos elementos originais
     */
    static void ModifyElementsBasedOnRefpFound(TPZAutoPointer<TPZRefPattern> & refpFound,
                                               TPZAutoPointer<TPZRefPattern> & refp,
                                               TPZVec<TPZGeoEl *> &elementVec);    
    
	/**
	 * @brief Return an refpattern based on a gived one (modelPat), whose midnodes was dragged to match with a geoel neighbourhood refinement
	 * @param gel - input data: geometric element for which the model refpattern nodes will be dragged
	 * @param modelPat - input data: Model RefPattern that is topologicaly compatible with neighbourhood
	 * @param neighCorresp - input data: map that group (neighbour geoelement) and nodes correspondences between (neighbour->SideRefpattern) and (gel->SideRefpattern), indexed by (gel->Side)
	 */
	static TPZAutoPointer<TPZRefPattern> DragModelPatNodes(TPZGeoEl * gel,
                                                           TPZAutoPointer<TPZRefPattern> modelPat,
                                                           std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > > &neighCorresp);
	
	/**
	 * @brief Returns if the given refPatterns (refA and refB) are topologicaly compatibles.
	 * @note If they are, pairNodes represents the correspondence between nodesIds from refAmesh to refBmesh.
	 * @param refA - input data: first refpattern to be compared
	 * @param refB - input data: second refpattern to be compared
	 * @param fromAtoB - input data: Linear Transformation from refA->RefPatternMesh to refB->RefPatternMesh() (its needed to pair Nodes correctly in permuted cases)
	 * @param pairNodes - output data: correspondence between refA and refB nodes, in case they are topologicaly compatibles
	 */
	static bool CompareTopologies(TPZAutoPointer<TPZRefPattern> refA,
                                  TPZAutoPointer<TPZRefPattern> refB,
                                  TPZTransform<> &fromAtoB,
                                  std::map<int, int> &pairNodes);
	
	/**
	 * @brief This method pair CORNER nodes from refA->mesh.father to refB->mesh.father, using the givem transformation from refA->mesh to refB->mesh to match coordinates.
	 * @param meshA - input data: mesh of first refpattern that will be considered its nodes coordinates
	 * @param meshB - input data: mesh of second refpattern that will be considered its nodes coordinates
	 * @param fromAtoB - input data: Linear Transformation from refA->mesh to refB->mesh (its needed to pair Nodes correctly in permuted cases)
	 * @param pairedNodes - output data: correspondence between refA->mesh and refB->mesh nodes that matches its coordinates in fathers elements
	 */
	/**
	 * The output is the map pairNodes, thar represents the A_nodeId paired with B_nodeId. \n
	 * Obs.: Be careful with the output interpretation! It contains the nodeIds, NOT the nodes positions in mesh.NodeVec()!!!
	 */
	static void PairMeshesCornerNodesMatchingCoordinates(TPZGeoMesh &meshA,
                                                         TPZGeoMesh &meshB,
                                                         TPZTransform<> &fromAtoB,
                                                         std::map<int, int> &pairedNodes);
	
	/**
	 * @brief This method pair nodes from refA->mesh to refB->mesh, using the givem transformation from refA->mesh to refB->mesh to match coordinates.
	 * @note The output is the map pairNodes, thar represents the A_nodeId paired with B_nodeId. \n
	 * Obs: Be careful with the output interpretation! It contains the nodeIds, NOT the nodes positions in mesh.NodeVec()!!!
	 * @param meshA - input data: mesh of first refpattern that will be considered its nodes coordinates
	 * @param meshB - input data: mesh of second refpattern that will be considered its nodes coordinates
	 * @param fromAtoB - input data: Linear Transformation from refA->mesh to refB->mesh (its needed to pair Nodes correctly in permuted cases)
	 * @param pairedNodes - output data: correspondence between refA->mesh and refB->mesh nodes that matches its coordinates
	 */
	static void PairMeshesNodesMatchingCoordinates(TPZGeoMesh &meshA,
                                                   TPZGeoMesh &meshB,
                                                   TPZTransform<> &fromAtoB,
                                                   std::map<int, int> &pairedNodes);
	
	/** @brief Returns the the name of refpattern model. */
	/** 
	 * To do this, it starts with the 3 initial characters of element nametype,
	 * followed by the quantity of midnodes for each side of element.
	 */
	static std::string BuildRefPatternModelName(TPZRefPattern &refp);
	static std::string BuildRefPatternModelName(TPZAutoPointer<TPZRefPattern> refp);
	static std::string BuildRefPatternModelName(TPZGeoEl *gel);
	
	/**
	 * @brief Returns if there is any neigbour already refined
	 * @param gel - input data: geometric element whose refinements of the neighbors will define the refinement of its sides
	 * @param sidestoref - output data: vector whose positions mention the sides of gel, and its contents mention the intensity of refinement of the respective side
	 */
	static bool SidesToRefine(TPZGeoEl *gel, TPZVec<int> &sidestoref);
	
	/** @brief Refines the element if it touches an element with a material id included in matids */
	static void RefineDirectional(TPZGeoEl *gel, std::set<int> &matids);
	static void RefineDirectional(TPZGeoEl *gel, std::set<int> &matids, int gelMat);
    
    static void RefineDirectional(TPZGeoMesh *gmesh, std::set<int> &matids);
    static void RefineDirectional(TPZGeoMesh *gmesh, std::set<int> &matids, int gelmat);
    
	
	static void RefineUniformIfNeighMat(TPZGeoEl *gel, std::set<int> &matids);
	
	/** @brief Method to test if the jacobian of a TPZGeoElSide element is constant */
	static bool ConstJacobian(TPZGeoElSide gelside, REAL tol = 1.e-6);
	
	/**
	 * @brief Algorithm that evaluates the veracity of the hashings between sides
     * of the elements children and corresponding sides of the father.
	 */
	/** 
	 * A point p in the parametric space of the side of the sub-element is overcome and is calculated \n
	 * it mentioned hashing getting pf point in the element father. \n
	 * One calculates for p and pf the corresponding deformed point. \n
	 * Itself the hashing is consistent the deformed point must the same be.
     */
	static void TransformationTest(TPZRefPattern * refp);
	
	/**
	 * @brief NodesHunted vector is the sequential nodesIds that belongs (i.e.: "Tol" far) to InitialNode(IdIni)~FinalNode(IdFin) alignment of gMesh.NodeVec()
	 * @note Obs.: InitialNode and FinalNode are also included!!!
	 */
	static void NodesHunter(TPZGeoMesh &gMesh,
                            TPZVec<int>& NodesHunted,
                            int IdIni,
                            int IdFin,
                            double Tol = 1.E-1);
	
	/**
	 * @brief Fill the TPZVec "permutation" with the valid permutations of "gel"
	 * @note Note: The permutations is with respect to Master Element nodes, NOT Gel nodes in a geomesh context (i.e.: NOT geomesh nodes ids)
	 */
	static void GetGelPermutations(TPZGeoEl * gel,
                                   TPZVec< TPZManVector<int,8> > &permutation);
	
	/** @brief Fill the TPZVec "permutation" with the valid permutations of a given element type */
	static void GetElTypePermutations(MElementType elType,
                                      TPZVec< TPZManVector<int, 8> > &permutation);
};

#endif

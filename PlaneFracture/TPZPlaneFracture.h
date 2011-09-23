/**
 * @file
 * @brief Contains the TPZPlaneFracture class which defines a plane of the fracture.
 */
#ifndef TPZPLANEFRACTUREH
#define TPZPLANEFRACTUREH

/*
 *  TPZPlaneFracture.h
 *  Crack
 *
 *  Created by Cesar Lucci on 09/08/10.
 *  Copyright 2010 LabMeC. All rights reserved.
 *
 */

using namespace std;

#include <set>
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzvec.h"

/** @brief Initializing for refpattern modulation */
const int __minTrimQTD = 4;
/** @brief Initializing for soft deviation of node (will be use to multiply the __minTrimQTD) */
const int __TrimQTDmultiplier = 5;

/** 
 * @brief Plane of the fracture. 
 * @author Cesar Lucci
 * @since 09/08/2010
 */
class TPZPlaneFracture
{
	public:
	
	/**
	 * @brief Constructor
	 * @param planeMesh [in] bidimensional mesh in R3 space
	 * @param nodeOrigin [in] Id of node of mesh that will be considered origin of local R2 coordinate system
	 * @param nodeX [in] Id of node of mesh that will define the X-axis direction of local R2 coordinate system
	 * @param nodeY [in] Id of node of mesh that will define the Y-axis direction of local R2 coordinate system
	 * @param TrimQTD [in] quantity of edges stretches that defines possibles trimCoords
	 * @note Obs.1 : The local R2 coordinate system, defined by 3 nodes described above,
	 *             will be orthogonalized and normalized (Gram-Schmidt) in this method, so dont worry too much!
	 *
	 * Obs.2 : The given planeMesh will be preserved as class atribute, so many poligonalChain
	 *             geometry could be created without change the original planeMesh
	 */
	TPZPlaneFracture(TPZGeoMesh * planeMesh, int nodeOrigin, int nodeX, int nodeY, int TrimQTD = __minTrimQTD);
	~TPZPlaneFracture();
	
	/**
	 * @brief Returns an GeoMesh based on original planeMesh, contemplating the poligonalChains geometry by refined elements
	 * @param poligonalChain [in] vector of boundary points coordinates
	 *
	 * @note Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]\n
	 *		y coordinate of first point of crack boundary: poligonalChain[1]\n
	 *		z coordinate of first point of crack boundary: poligonalChain[2]\n
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]\n
	 *		y coordinate of second point of crack boundary: poligonalChain[4]\n
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 */
	TPZGeoMesh * GetFractureMesh(TPZVec<REAL> &poligonalChain);
		
//-----------------------------------------------------------------------------------------------------------------------------------------------------
	
	private:
	
	/*
	 * @brief Computes the edges of elements of fractMesh that are intercepted by the crack tip defined by poligonalChain points (defined by a vector coordinates)
	 * @param poligonalChain [in] vector of boundary points coordinates
	 *
	 * @note Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]\n
	 *		y coordinate of first point of crack boundary: poligonalChain[1]\n
	 *		z coordinate of first point of crack boundary: poligonalChain[2]\n
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]\n
	 *		y coordinate of second point of crack boundary: poligonalChain[4]\n
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 *
	 * @param fractMesh [in] geomesh of considered elements
	 * @param elId_TrimCoords [out] map that contains 1D element Id and a set of it trim 1D coordinates
	 * @param elIdSequence [out] the same of elId_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 */
	void DetectEdgesCrossed(TPZVec<REAL> &poligonalChain, TPZGeoMesh * fractMesh,
							std::map< int, std::set<double> > &elId_TrimCoords, std::list< std::pair<int,double> > &elIdSequence);
	
	/**
	 * @brief Returns (gel->Neighbour) relative to the side intersepted by the x+alpha.dx line
	 * @param gel [in] gel crossed by the line
	 * @param x [input and output data]  
	 *				x [as input] start point of line \n
	 *				x [as output] end point of line in gel and gel->Neighbour interface
	 * @param dx [in] direction of line from point x (input)
	 * @param alphaMin [in] if an start point (x) is in already in one edge of gel, it might be included or not in the intersections \n
	 *				        so, using alphaMin=0, in this case the first intersection (the x itself) is included...
	 *							   using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 * @param elId_TrimCoords [out] map that contains the trim coordinates of 1D element, indexed by its Id (once 1D element was inserted in gel->Mesh)\n
	 *							    obs.: elId_TrimCoords was idealized to work in accumulative conception, i.e.,
	 *								    each time this method is called, this map grows!
	 * @param elIdSequence [out] the same of elId_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 * @param pushback [in] set if dots on element edges will be inserted at the end, or not (i.e.: at beggining), of fCrackBoundary list
	 */
	TPZGeoEl * CrossToNextNeighbour(TPZGeoEl * gel, TPZFMatrix &x, TPZFMatrix dx, double alphaMin, std::map< int,
									std::set<double> > &elId_TrimCoords, std::list< std::pair<int,double> > &elIdSequence, bool pushback);
	
	/**
	 * @brief For a given element and internal point and an direction, computes the intersection
	 * coordinates with respect to it edges, and the respective intersected edge.
	 *
	 * @param gel [in] 2D geometric element whose edge will be intersected by (x+alphaX.dx) line
	 * @param x [in] element internal point coordinates
	 * @param dx [in] direction from point p
	 * @param edge [out] side Id of element edge that will be intersected by (x+alphaX.dx) line
	 * @param ExactIntersect [out] exact intersection coordinates with respect to edges parameter
	 * @param ModulatedIntersect [out] exact intersection coordinates, dragged to the nearest module (defined by fTrimQTD atribute)
	 * @param alphaMin [in] if an start point (x) is in already in one edge of gel, it might be included or not in the intersections\n
	 *					    so, using alphaMin=0, in this case the first intersection (the x itself) is included...\n
	 *						using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 */
	bool EdgeIntersection(TPZGeoEl * gel, TPZFMatrix &x, TPZFMatrix &dx, TPZVec<int> &edge,
						  TPZVec< TPZFMatrix > &ExactIntersect, TPZVec< TPZFMatrix > &ModulatedIntersect, double alphaMin);
	
	/*
	 * @brief Returns an pointer to element of gMesh that contains the given point p (i.e.: point "p" belongs to element domain)
	 * @param p [in] point whose elements is going to be localized
	 * @param fractMesh [in] geomesh of elements candidates
	 * @param poligonalChain [in] vector of boundary points coordinates
	 *
	 * @note Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]\n
	 *		y coordinate of first point of crack boundary: poligonalChain[1]\n
	 *		z coordinate of first point of crack boundary: poligonalChain[2]\n
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]\n
	 *		y coordinate of second point of crack boundary: poligonalChain[4]\n
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 */
	static TPZGeoEl * PointElement(int p, TPZGeoMesh * fractMesh, TPZVec<REAL> &poligonalChain);
	
	// alphaNode eh uma das solucoes do sistema: {x + alphaX.dx == node + alphaNode.dnode}, ou seja,
	// a norma que multiplica o vetor dnode e cruza a reta (x+alphaX.dx)
	/**
	 * @brief Given two vectorial lines \f$ x + alphaX.dx\f$ and \f$ node + alphaNode.dnode\f$, \n
	 * this method returns the alphaNode (norm that multiplies the unit vector dnode to intersect the line (x + alphax.dx) )
	 * @param x given point
	 * @param dx direction of line from point x
	 * @param node connect
	 * @param dnode
	 * @param norm [in] norm of edge that belongs to (node + alphaNode.dnode) line
	 * @param modulate [in] set if alphaNode will be modulated by stretches
	 * @param smooth [in] if alphaNode will be modulated, set if the stretches will be (norm/fTrimQTD) or smallest stretches \n
	 *							defined by (norm/(fTrimQTD*__TrimQTDmultiplier)).
	 * @note Obs.: alphaNode modulation is useful to reduce the possibilities of non regular refpatterns.
 	 * @note OBS.: dx and dnode MUST BE UNIT VECTORS!!!
	 */
	double ComputeAlphaNode(TPZFMatrix &x, TPZFMatrix &dx, TPZFMatrix &node, TPZFMatrix &dnode, double norm, bool modulate, bool smooth);
	
	// alphaX eh uma das solucoes do sistema: {x + alphaX.dx == node + alphaNode.dnode}, ou seja,
	// a norma que multiplica o vetor dx e cruza a reta (node+alphaNode.dnode)
	/**
	 * @brief Given two vectorial lines \f$ x + alphaX.dx\f$ and \f$ node + alphaNode.dnode\f$, \n
	 * this method returns the alphaX (norm that multiplies the unit vector dx to intersect the line (none + alphaNode.dnode) )
	 * @note dx and dnode MUST BE UNIT VECTORS!!!
	 */
	double ComputeAlphaX(TPZFMatrix &x, TPZFMatrix &dx, TPZFMatrix &node, TPZFMatrix &dnode);
	
	/**
	 * @brief Return if a given point x is near to some node of a given geo element
	 * @param gel [in] given geo element
	 * @param x [in] given point
	 * @param node [out] id of node that is in the x range
	 * @param tol [in] x range radius
	 */
	static bool NearestNode(TPZGeoEl * gel, TPZFMatrix &x, int &node, double tol);
	
	static int NearestNode(TPZGeoMesh * gmesh, TPZVec<REAL> &x, double tol);
	
	/**
	 * @brief Given 2 nodes (n0 and n1) and one point (x) in \f$ n0->n1 \f$ line, returns the point x in the line parametric space \f$ [-1,+1]\f$
	 */
	static double LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1);
	
	/**
	 * @brief This method return a reffpattern of an unidimentional element that matches with the trim coordinates.
	 * @param TrimCoord Set of 1D element trimmed coordinates ( \f$[ -1 , +1 ]\f$ domain )
	 */
	static TPZAutoPointer<TPZRefPattern> Generate1DRefPatt(std::set<double> &TrimCoord);
	
	/**
	 * @brief Updates poligonal chain.
	 * @note The original Poligonal Chain (input data on GetFractureMesh method) is dots coordinates in the 2D mesh. This points normally are inside elements domain.\n
	 * The edges intersections of the original Poligonal Chain originate a new Poligonal Chain named poligonalChainUpdated 
	 */
	static void UpdatePoligonalChain(TPZGeoMesh * gmesh, std::list< std::pair<int,double> > &elIdSequence,
							  TPZVec<REAL> &poligonalChainUpdated);
	
	/**
	 * @param gmesh geometric mesh
	 * @param elIdSequence - output data: list that contains 1D element Id and it trim 1D coordinates in generation sequence order
	 */
	void GenerateCrackBoundary(TPZGeoMesh * gmesh, std::list< std::pair<int,double> > &elIdSequence);
	
	
//--------------------------------------------------------------------------------------------------------------------------------------------------
	
	protected:
	
	/** @brief Original mesh (keeped intact for any poligonalChain configuration) */
	const TPZGeoMesh * fplaneMesh;
	/** @brief It limits the amount of possible points in the edge of the elements */
	int fTrimQTD;
	/** @brief Transformation from crack plane in R3 to crack plane in R2 */
	TPZFMatrix fFromR3toR2;
};

#endif
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

const int __minTrimQTD = 4;//for refpattern modulation
const int __TrimQTDmultiplier = 5;//for soft deviation of node (will be use to multiply the __minTrimQTD)

class TPZPlaneFracture
{
	public:
	
	/**
	 * @param planeMesh - input data: bidimensional mesh in R3 space
	 * @param nodeOrigin - input data: Id of node of mesh that will be considered origin of local R2 coordinate system
	 * @param nodeX - input data: Id of node of mesh that will define the X-axis direction of local R2 coordinate system
	 * @param nodeY - input data: Id of node of mesh that will define the Y-axis direction of local R2 coordinate system
	 * @param TrimQTD - inout data: quantity of edges stretches that defines possibles trimCoords
	 * Obs.1 : The local R2 coordinate system, defined by 3 nodes described above,
	 *             will be orthogonalized and normalized (Gram-Schmidt) in this method, so dont worry too much!
	 *
	 * Obs.2 : The given planeMesh will be preserved as class atribute, so many poligonalChain
	 *             geometry could be created without change the original planeMesh
	 */
	TPZPlaneFracture(TPZGeoMesh * planeMesh, int nodeOrigin, int nodeX, int nodeY, int TrimQTD = __minTrimQTD);
	~TPZPlaneFracture();
	
	/**
	 * Return an GeoMesh based on original planeMesh, contemplating the poligonalChains geometry by refined elements
	 * @param poligonalChain - input data: vector of boundary points coordinates
	 *
	 * Note: Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]
	 *		y coordinate of first point of crack boundary: poligonalChain[1]
	 *		z coordinate of first point of crack boundary: poligonalChain[2]
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]
	 *		y coordinate of second point of crack boundary: poligonalChain[4]
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 */
	TPZGeoMesh * GetFractureMesh(TPZVec<REAL> &poligonalChain);
		
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private:
	
	/*
	 * Computes the edges of elements of fractMesh that are intercepted by the crack tip defined by poligonalChain points (defined by a vector coordinates)
	 * @param poligonalChain - input data: vector of boundary points coordinates
	 *
	 * Note: Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]
	 *		y coordinate of first point of crack boundary: poligonalChain[1]
	 *		z coordinate of first point of crack boundary: poligonalChain[2]
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]
	 *		y coordinate of second point of crack boundary: poligonalChain[4]
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 *
	 * @param fractMesh - input data: geomesh of considered elements
	 * @param elId_TrimCoords - output data: map that contains 1D element Id and a set of it trim 1D coordinates
	 * @param elIdSequence - output data: the same of elId_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 */
	void DetectEdgesCrossed(TPZVec<REAL> &poligonalChain, TPZGeoMesh * fractMesh,
							std::map< int, std::set<double> > &elId_TrimCoords, std::list< std::pair<int,double> > &elIdSequence);
	
	/**
	 * Returns (gel->Neighbour) relative to the side intersepted by the x+alpha.dx line
	 * @param gel - input data: gel crossed by the line
	 * @param x - input and output data: 
	 *				x - as input: start point of line
	 *				x - as output: end point of line in gel and gel->Neighbour interface
	 * @param dx - input data: direction of line from point x (input)
	 * @param alphaMin - input data: if an start point (x) is in already in one edge of gel, it might be included or not in the intersections
	 *							   so, using alphaMin=0, in this case the first intersection (the x itself) is included...
	 *							   using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 * @param elId_TrimCoords - output data: map that contains the trim coordinates of 1D element, indexed by its Id (once 1D element was inserted in gel->Mesh)
	 *							   obs.: elId_TrimCoords was idealized to work in accumulative conception, i.e.,
	 *								    each time this method is called, this map grows!
	 * @param elIdSequence - output data: the same of elId_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 * @param pushback - input data: set if dots on element edges will be inserted at the end, or not (i.e.: at beggining), of fCrackBoundary list
	 */
	TPZGeoEl * CrossToNextNeighbour(TPZGeoEl * gel, TPZFMatrix &x, TPZFMatrix dx, double alphaMin, std::map< int,
									std::set<double> > &elId_TrimCoords, std::list< std::pair<int,double> > &elIdSequence, bool pushback);
	
	/**
	 * For a given element and internal point and an direction, computes the intersection
	 * coordinates with respect to it edges, and the respective intersected edge.
	 *
	 * @param gel - input data: 2D geometric element whose edge will be intersected by (x+alphaX.dx) line
	 * @param x - input data: element internal point coordinates
	 * @param dx - input data: direction from point p
	 * @param edge - output data: side Id of element edge that will be intersected by (x+alphaX.dx) line
	 * @param ExactIntersect - output data: exact intersection coordinates with respect to edges parameter
	 * @param ModulatedIntersect - output data: exact intersection coordinates, dragged to the nearest module (defined by fTrimQTD atribute)
	 * @param alphaMin - input data: if an start point (x) is in already in one edge of gel, it might be included or not in the intersections
	 *							   so, using alphaMin=0, in this case the first intersection (the x itself) is included...
	 *							   using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 */
	bool EdgeIntersection(TPZGeoEl * gel, TPZFMatrix &x, TPZFMatrix &dx, TPZVec<int> &edge,
						  TPZVec< TPZFMatrix > &ExactIntersect, TPZVec< TPZFMatrix > &ModulatedIntersect, double alphaMin);
	
	/*
	 * Return an pointer to element of gMesh that contains the given point p (i.e.: point "p" belongs to element domain)
	 * @param p - input data: point whose elements is going to be localized
	 * @param fractMesh - input data: geomesh of elements candidates
	 * @param poligonalChain - input data: vector of boundary points coordinates
	 *
	 * Note: Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]
	 *		y coordinate of first point of crack boundary: poligonalChain[1]
	 *		z coordinate of first point of crack boundary: poligonalChain[2]
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]
	 *		y coordinate of second point of crack boundary: poligonalChain[4]
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 */
	static TPZGeoEl * PointElement(int p, TPZGeoMesh * fractMesh, TPZVec<REAL> &poligonalChain);
	
	// alphaNode eh uma das solucoes do sistema: {x + alphaX.dx == node + alphaNode.dnode}, ou seja,
	// a norma que multiplica o vetor dnode e cruza a reta (x+alphaX.dx)
	/**
	 * Given two vectorial lines x + alphaX.dx and node + alphaNode.dnode,
	 * this method returns the alphaNode (norm that multiplies the unit vector dnode to intersect the line (x + alphax.dx) )
	 * @param norm - input data: norm of edge that belongs to (node + alphaNode.dnode) line
	 * @param modulate - input data: set if alphaNode will be modulated by stretches
	 * @param smooth - input data: if alphaNode will be modulated, set if the stretches will be (norm/fTrimQTD) or smallest stretches
	 *							defined by (norm/(fTrimQTD*__TrimQTDmultiplier)).
	 * Obs.: alphaNode modulation is useful to reduce the possibilities of non regular refpatterns.
 	 * OBS.: dx and dnode MUST BE UNIT VECTORS!!!
	 */
	double ComputeAlphaNode(TPZFMatrix &x, TPZFMatrix &dx, TPZFMatrix &node, TPZFMatrix &dnode, double norm, bool modulate, bool smooth);
	
	// alphaX eh uma das solucoes do sistema: {x + alphaX.dx == node + alphaNode.dnode}, ou seja,
	// a norma que multiplica o vetor dx e cruza a reta (node+alphaNode.dnode)
	/**
	 * Given two vectorial lines x + alphaX.dx and node + alphaNode.dnode,
	 * this method returns the alphaX (norm that multiplies the unit vector dx to intersect the line (none + alphaNode.dnode) )
	 * OBS.: dx and dnode MUST BE UNIT VECTORS!!!
	 */
	double ComputeAlphaX(TPZFMatrix &x, TPZFMatrix &dx, TPZFMatrix &node, TPZFMatrix &dnode);
	
	/**
	 * Return if a given point x is near to some node of a given geo element
	 * @param gel - input data: given geo element
	 * @param x - input data: given point
	 * @param node - output data: id of node that is in the x range
	 * @param tol - input data: x range radius
	 */
	static bool NearestNode(TPZGeoEl * gel, TPZFMatrix &x, int &node, double tol);
	
	static int NearestNode(TPZGeoMesh * gmesh, TPZVec<REAL> &x, double tol);
	
	/**
	 * Given 2 nodes (n0 and n1) and one point (x) in n0->n1 line, returns the point x in the line parametric space [-1,+1]
	 */
	static double LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1);
	
	/**
	 * This method return a reffpattern of an unidimentional element that matches with the trim coordinates.
	 * @param TrimCoord - set of 1D element trimmed coordinates ( [ -1 , +1 ] domain )
	 */
	static TPZAutoPointer<TPZRefPattern> Generate1DRefPatt(std::set<double> &TrimCoord);
	
	/**
	 * The original Poligonal Chain (input data on GetFractureMesh method) is dots coordinates in the 2D mesh. This points normally are inside elements domain.
	 * The edges intersections of the original Poligonal Chain originate a new Poligonal Chain named poligonalChainUpdated 
	 */
	static void UpdatePoligonalChain(TPZGeoMesh * gmesh, std::list< std::pair<int,double> > &elIdSequence,
							  TPZVec<REAL> &poligonalChainUpdated);
	
	/**
	 * @param elIdSequence - output data: list that contains 1D element Id and it trim 1D coordinates in generation sequence order
	 */
	void GenerateCrackBoundary(TPZGeoMesh * gmesh, std::list< std::pair<int,double> > &elIdSequence);
	
	
//---------------------------------------------------------------------------------------------------------------------------------------------------------------	
	
	protected:
	
	TPZGeoMesh * fplaneMesh;//original mesh (keeped intact for any poligonalChain configuration)
	int fTrimQTD;//It limits the amount of possible points in the edge of the elements
	TPZFMatrix fFromR3toR2;//Linear transformation from R3 to R2
};

#endif
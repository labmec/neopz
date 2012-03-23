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

class TPZPlaneFracture
{
	public:
	
	/**
	 * @param planeMesh - input data: bidimensional mesh in R3 space
	 * @param nodeOrigin - input data: Id of node of mesh that will be considered origin of local R2 coordinate system
	 * @param nodeX - input data: Id of node of mesh that will define the X-axis direction of local R2 coordinate system
	 * @param nodeY - input data: Id of node of mesh that will define the Y-axis direction of local R2 coordinate system
	 * Obs.1 : The local R2 coordinate system, defined by 3 nodes described above,
	 *             will be orthogonalized and normalized (Gram-Schmidt) in this method, so dont worry too much!
	 *
	 * Obs.2 : The given planeMesh will be preserved as class atribute, so many crackTip
	 *             geometry could be created without change the original planeMesh
	 */
	TPZPlaneFracture(TPZGeoMesh * planeMesh, int nodeOrigin, int nodeX, int nodeY);
	~TPZPlaneFracture();
	
	/**
	 * This method return a reffpattern of an unidimentional element, folowing the "ContourIntersections" instructions.
	 * @param TrimCoord - set of 1D element trimmed coordinates ( [ -1 , +1 ] domain )
	 * @param gmesh - F.E.M. geometric mesh where the element refined will be inserted
	 */
	TPZAutoPointer<TPZRefPattern> Generate1DRefPatt(set<double> &TrimCoord);
	
	double LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1);
	
	/**
	 * For a given element and internal point and an direction, returns the intersection
	 * coordinates with respect to it edges, and the respective intersected edge.
	 *
	 * @param gel - input data: 2D geometric element whose edge will be intersected by (p+alphaP.dp) line
	 * @param p - input data: element internal point coordinates
	 * @param dp - input data: direction from point p
	 * @param edge - output data: side Id of element edge that will be intersected by (p+alphaP.dp) line
	 */
	TPZVec<REAL> EdgeIntersection(TPZGeoEl * gel, const TPZVec<REAL> &p, const TPZVec<REAL> &dp, int &edge);
	
	/**
	 * Return an GeoMesh based on original planeMesh, contemplating the crackTips geometry by refined elements
	 * @param crackTip - input data: vector of crack boundaries with vector of boundary points coordinates
	 *
	 * Note: Each vector position is another vector that store x, y and z coordinates IN SEQUENCE of an crackTip geometry.
	 *
	 * Example:
	 *
	 *	First crack boundary: crackTip[0];
	 *		x coordinate of first point of first crack boundary: crackTip[0][0]
	 *		y coordinate of first point of first crack boundary: crackTip[0][1]
	 *		z coordinate of first point of first crack boundary: crackTip[0][2]
	 *		//
	 *		x coordinate of second point of first crack boundary: crackTip[0][3]
	 *		y coordinate of second point of first crack boundary: crackTip[0][4]
	 *		z coordinate of second point of first crack boundary: crackTip[0][5]
	 *
	 *	Second crack boundary (islands): crackTip[1];
	 *		x coordinate of first point of second crack boundary: crackTip[1][0]
	 *		y coordinate of first point of second crack boundary: crackTip[1][1]
	 *		z coordinate of first point of second crack boundary: crackTip[1][2]
	 *		//
	 *		x coordinate of second point of second crack boundary: crackTip[1][3]
	 *		y coordinate of second point of second crack boundary: crackTip[1][4]
	 *		z coordinate of second point of second crack boundary: crackTip[1][5]
	 */
	TPZGeoMesh * GetFractureMesh(TPZVec< TPZVec<REAL> > &crackTip);
	
	protected:
	TPZGeoMesh * fplaneMesh;
	
	TPZFMatrix fFromR3toR2;
};

#endif
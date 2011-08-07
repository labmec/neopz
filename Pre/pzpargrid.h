/**
 * @file
 * @brief Contains the TPZGenPartialGrid class which implements the generation of a geometric grid.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO/Philippe Devloo
//
// File:   tpargrid.h
//
// Class:  tparrid
//
// Obs.:   Gera uma malha retangular:
//
// Versao: 06 / 1996.
//


#ifndef _TPZPARGRIDHH_
#define _TPZPARGRIDHH_

class TPZCompMesh;
class TPZGeoMesh;

#include <stdio.h>
#include <iostream>
#include "pzreal.h"
#include "pzvec.h"

/**
 * @ingroup pre
 * @brief Implements the generation of a geometric grid. \ref pre "Getting Data"
 */
/** Implements the generation of part of the grid
 * This class uses DEPRECATED objects, but can be easily updated
 */
class TPZGenPartialGrid{
public:
	/**
	 @param x0 lower left coordinate
	 @param x1 upper right coordinate
	 @param nx number of nodes in x and y
	 @param rangex range of nodes which need to be created
	 @param rangey range of nodes which need to be created
	 */
	TPZGenPartialGrid(TPZVec<int> &nx, TPZVec<int> &rangex, TPZVec<int> &rangey, TPZVec<REAL> &x0, TPZVec<REAL> &x1);
	
	~TPZGenPartialGrid();
	
	short Read (TPZGeoMesh & malha);
	
	void SetBC(TPZGeoMesh &gr, int side, int bc);
	
	void Print( char *name = NULL, std::ostream &out = std::cout );
	
	void SetElementType(int type) {
		fElementType = type;
	}
	
protected:
	
	void Coord(int i, TPZVec<REAL> &coord);
	
	int NodeIndex(int i, int j);
	
	int ElementIndex(int i, int j);
	
	void ElementConnectivity(int iel, TPZVec<int> &nodes);
	
	TPZVec<int> fNx;
	TPZVec<int> fRangex,fRangey;
	TPZVec<REAL> fX0,fX1,fDelx;
	int fNumNodes;
	int fElementType;
	
};

#endif // _TPZGENGRIDHH_

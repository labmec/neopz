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


/*************/
class TPZGenPartialGrid{
public:
	TPZGenPartialGrid(TPZVec<int> &nx, TPZVec<int> &rangex, TPZVec<int> &rangey, TPZVec<REAL> &x0, TPZVec<REAL> &x1);

	~TPZGenPartialGrid();

	short Read (TPZGeoMesh & malha);

void SetBC(TPZGeoMesh &gr, int side, int bc);

	void Print( char *name = NULL, ostream &out = cout );

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





//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   TGenGrid.h
//
// Class:  TGenGrid
//
// Obs.:   Gera uma malha retangular:
//
// Versao: 06 / 1996.
//


#ifndef _TPZGENGRIDHH_
#define _TPZGENGRIDHH_

class TPZCompMesh;
class TPZGeoMesh;
#include "pzvec.h"

#include <stdio.h>
#include <iostream>
#include "pzreal.h"

#include <fstream>
/*************/
class TPZGenGrid{

public:

TPZGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numl = 1, REAL rot = 0.5);

virtual ~TPZGenGrid();

virtual short Read (TPZGeoMesh & malha);

virtual void SetBC(TPZGeoMesh *gr, int side, int bc);

virtual void SetBC(TPZGeoMesh *g, TPZVec<REAL> &start, TPZVec<REAL> &end, int bc);

virtual void Print( char *name = NULL, ostream &out = cout );

virtual void SetElementType(int type);

   int ElemId(int iel,int jel, int layer);

   int SuperElemId(int iblock,int jblock);

   void CreateSuperElements(TPZCompMesh &malha,int length);

	REAL Distance(TPZVec<REAL> &x1,TPZVec<REAL> &x2);




protected:

virtual void Coord(int i, TPZVec<REAL> &coord);

virtual int GlobalI(int ix, int iy, int layer);

	void ElementConnectivity(int iel, TPZVec<int> &nodes);

virtual void GenerateNodes(TPZGeoMesh &grid);

virtual void GenerateElements(TPZGeoMesh &grid);


	TPZVec<int> fNx;
	TPZVec<REAL> fX0,fX1,fDelx;
	int fNumNodes;
	int fElementType;

   int fSuperElemlength;
   int 	fNumBlocks[2];

   int fNumLayers;
   REAL fRotAngle;



};

#endif // _TGENGRIDHH_

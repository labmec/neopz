/**
 * @file
 * @brief Contains the implementation of the TPZGraphElT2dMapped methods. 
 */
//
// C++ Implementation: tpzgraphelt2dmapped
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "tpzgraphelt2dmapped.h"
#include "pzshapequad.h"

static REAL cornerco[4][2] = 
{
	{0.,0.},
	{1.,0.},
	{0.,1.},
	{0.,1.}
};

TPZGraphElT2dMapped::~TPZGraphElT2dMapped()
{
}

void TPZGraphElT2dMapped::QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta)
{
	TPZGraphElQ2dd::QsiEta(i,imax,qsieta);
	TPZFNMatrix<8> phi(4,1,0.),dphi(2,4,0.);
	pzshape::TPZShapeQuad::ShapeCorner(qsieta,phi,dphi);
	REAL temp[2] = {0.,0.};
	int is;
	for(is=0; is<4; is++)
	{
		temp[0] += cornerco[is][0]*phi(is,0);
		temp[1] += cornerco[is][1]*phi(is,0);
	}
	qsieta[0] = temp[0];
	qsieta[1] = temp[1];
	
}


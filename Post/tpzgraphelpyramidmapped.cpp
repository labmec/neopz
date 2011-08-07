/**
 * @file
 * @brief Contains the implementation of the TPZGraphElPyramidMapped methods. 
 */
//
// C++ Implementation: tpzgraphelpyramidmapped
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "tpzgraphelpyramidmapped.h"
#include "pzshapecube.h"

static REAL cornerco[8][3] = 
{
	{-1.,-1.,0.},
	{1.,-1.,0.},
	{1.,1.,0.},
	{-1.,1.,0.},
	{0.,0.,1.},
	{0.,0.,1.},
	{0.,0.,1.},
	{0.,0.,1.}
};

TPZGraphElPyramidMapped::TPZGraphElPyramidMapped(TPZCompEl* cel, TPZGraphMesh* gmesh): TPZGraphElQ3dd(cel, gmesh)
{
}


TPZGraphElPyramidMapped::~TPZGraphElPyramidMapped()
{
}

void TPZGraphElPyramidMapped::QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta)
{
	TPZGraphElQ3dd::QsiEta(i,imax,qsieta);
	TPZFNMatrix<24> phi(8,1,0.),dphi(3,8,0.);
	pzshape::TPZShapeCube::ShapeCorner(qsieta,phi,dphi);
	REAL temp[3] = {0.,0.,0.};
	int is;
	for(is=0; is<8; is++)
	{
		temp[0] += cornerco[is][0]*phi(is,0);
		temp[1] += cornerco[is][1]*phi(is,0);
		temp[2] += cornerco[is][2]*phi(is,0);
	}
	qsieta[0] = temp[0];
	qsieta[1] = temp[1];
	qsieta[2] = temp[2];
	
}




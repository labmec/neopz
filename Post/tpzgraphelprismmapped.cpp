/**
 * @file
 * @brief Contains the implementation of the TPZGraphElPrismMapped methods. 
 */

#include "tpzgraphelprismmapped.h"
#include "pzshapecube.h"

static REAL cornerco[8][3] = 
{
	{0.,0.,-1.},
	{1.,0.,-1.},
	{0.,1.,-1.},
	{0.,1.,-1.},
	{0.,0.,1.},
	{1.,0.,1.},
	{0.,1.,1.},
	{0.,1.,1.}
};

TPZGraphElPrismMapped::TPZGraphElPrismMapped(TPZCompEl* cel, TPZGraphMesh* gmesh): TPZGraphElQ3dd(cel, gmesh)
{
}


TPZGraphElPrismMapped::~TPZGraphElPrismMapped()
{
}

void TPZGraphElPrismMapped::QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta)
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


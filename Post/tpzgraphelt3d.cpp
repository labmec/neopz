/**
 * @file
 * @brief Contains the implementation of the TPZGraphElT3d methods. 
 */

#include "tpzgraphelt3d.h"


TPZGraphElT3d::~TPZGraphElT3d()
{
}

/**
 * This method maps the index of a point to parameter space as a function
 * of the number of divisions
 */
void TPZGraphElT3d::QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta)
{
	TPZGraphElQ3dd::QsiEta(i,imax,qsieta);
	REAL temp[3];
	temp[0] = (-1.+qsieta[1])*(1.+qsieta[0])*(-1.+qsieta[2])/8.;
	temp[1] = -(1.+qsieta[1])*(-1.+qsieta[2])/4.;
	temp[2] = (1.+qsieta[2])/2.;
	qsieta[0] = temp[0];
	qsieta[1] = temp[1];
	qsieta[2] = temp[2];
}



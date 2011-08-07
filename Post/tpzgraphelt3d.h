/**
 * @file
 * @brief Contains the TPZGraphElT3d class which implements the graphical representation of a tetrahedra element.
 */
//
// C++ Interface: tpzgraphelt3d
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef TPZGRAPHELT3D_H
#define TPZGRAPHELT3D_H

#include "pzgraphelq3dd.h"

/**
 * @ingroup post
 * @brief Implements the graphical representation of a tetrahedra element. \ref post "Post processing"
 * @author Philippe R. B. Devloo
 */
class TPZGraphElT3d : public TPZGraphElQ3dd
{
public:
	TPZGraphElT3d(TPZCompEl *cel, TPZGraphMesh *gmesh) : TPZGraphElQ3dd(cel,gmesh){
	}
	
    ~TPZGraphElT3d();
	
	/**
	 * @brief This method maps the index of a point to parameter space as a function
	 * of the number of divisions
	 */
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);
	
};

#endif

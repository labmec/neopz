/**
 * @file
 * @brief Contains the TPZGraphElT3d class which implements the graphical representation of a tetrahedra element.
 */

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
	/** @brief Constructor for graphical tetrahedra element */
	TPZGraphElT3d(TPZCompEl *cel, TPZGraphMesh *gmesh) : TPZGraphElQ3dd(cel,gmesh){
	}
	/** @brief Default destructor */
    ~TPZGraphElT3d();
	
	/** @brief This method maps the index of a point to parameter space as a function of the number of divisions */
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta) override;
	
};

#endif

/**
 * @file
 * @brief Contains the TPZGraphElT2dMapped class which implements a graphical element for a triangle mapped into de quadrilateral element.
 */

#ifndef TPZGRAPHELT2DMAPPED_H
#define TPZGRAPHELT2DMAPPED_H

#include "pzgraphelq2dd.h"

/**
 * @ingroup post
 * @brief Implements a graphical element for a triangle mapped into de quadrilateral element. \ref post "Post processing"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZGraphElT2dMapped : public TPZGraphElQ2dd
{
public:
	/** @brief Constructor for graphical element */
	TPZGraphElT2dMapped(TPZCompEl *cel, TPZGraphMesh *gmesh) : TPZGraphElQ2dd(cel,gmesh){
	}
	/** @brief Default destructor */
    ~TPZGraphElT2dMapped();
	
	/** @brief This method maps the index of a point to parameter space as a function of the number of divisions */
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta) override;
	
};

#endif

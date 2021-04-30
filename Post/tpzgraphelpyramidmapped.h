/**
 * @file
 * @brief Contains the TPZGraphElPyramidMapped class which implements the graphical element for a pyramid using a map to the cube element.
 */

#ifndef TPZGRAPHELPYRAMIDMAPPED_H
#define TPZGRAPHELPYRAMIDMAPPED_H

#include "pzgraphelq3dd.h"

/**
 * @ingroup post
 * @brief Implements the graphical element for a pyramid using a map to the cube element. \ref post "Post processing"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZGraphElPyramidMapped : public TPZGraphElQ3dd
{
public:
	/** @brief Constructor for graphical pyramid element */
    TPZGraphElPyramidMapped(TPZCompEl* cel, TPZGraphMesh* gmesh);
	/** @brief Simple destructor */
    ~TPZGraphElPyramidMapped();
	
	/** @brief This method maps the index of a point to parameter space as a function of the number of divisions */
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta) override;
	
};

#endif

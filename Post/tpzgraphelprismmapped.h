/**
 * @file
 * @brief Contains the TPZGraphElPrismMapped class which implements the graphical element for a prism using a degenerated cube element.
 */

#ifndef TPZGRAPHELPRISMMAPPED_H
#define TPZGRAPHELPRISMMAPPED_H

#include "pzgraphelq3dd.h"

/**
 * @ingroup post
 * @brief Implements the graphical element for a prism using a degenerated cube element. \ref post "Post processing"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZGraphElPrismMapped : public TPZGraphElQ3dd
{
public:
	/** @brief Constructor for graphical prism element */
    TPZGraphElPrismMapped(TPZCompEl* cel, TPZGraphMesh* gmesh);
	/** @brief Simple destructor */
    ~TPZGraphElPrismMapped();
	
	/** @brief This method maps the index of a point to parameter space as a function of the number of divisions */
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta) override;
	
};

#endif

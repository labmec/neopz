/**
 * @file
 * @brief Contains the TPZGraphElPrismMapped class which implements the graphical element for a prism using a degenerated cube element.
 */
//
// C++ Interface: tpzgraphelprismmapped
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//

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
    TPZGraphElPrismMapped(TPZCompEl* cel, TPZGraphMesh* gmesh);
	
    ~TPZGraphElPrismMapped();
	
	/**
	 * @brief This method maps the index of a point to parameter space as a function
	 * of the number of divisions
	 */
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);
	
};

#endif

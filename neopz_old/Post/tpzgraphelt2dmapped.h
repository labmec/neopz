//
// C++ Interface: tpzgraphelt2dmapped
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZGRAPHELT2DMAPPED_H
#define TPZGRAPHELT2DMAPPED_H

#include "pzgraphelq2dd.h"

/**
This class implements a graphical element for a triangle mapped into de quadrilateral element

	@author Philippe R. B. Devloo <phil@fec.unicamp.br>
*/
class TPZGraphElT2dMapped : public TPZGraphElQ2dd
{
public:
  TPZGraphElT2dMapped(TPZCompEl *cel, TPZGraphMesh *gmesh) : TPZGraphElQ2dd(cel,gmesh){
  }

    ~TPZGraphElT2dMapped();

 /**
     * This method maps the index of a point to parameter space as a function
     * of the number of divisions
  */
  virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);

};

#endif

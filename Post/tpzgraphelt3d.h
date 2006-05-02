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
//
#ifndef TPZGRAPHELT3D_H
#define TPZGRAPHELT3D_H

#include <pzgraphelq3dd.h>

/**
This class implements the graphical representation of a tetrahedra elemennt

@author Philippe R. B. Devloo
*/
class TPZGraphElT3d : public TPZGraphElQ3dd
{
public:
  TPZGraphElT3d::TPZGraphElT3d(TPZCompEl *cel, TPZGraphMesh *gmesh) : TPZGraphElQ3dd(cel,gmesh){
  }


    ~TPZGraphElT3d();

 /**
     * This method maps the index of a point to parameter space as a function
     * of the number of divisions
  */
  virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);


};

#endif

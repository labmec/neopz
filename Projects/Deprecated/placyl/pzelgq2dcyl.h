/***************************************************************************
                          pzelgq2dcyl.h  -  description
                             -------------------
    begin                : Thu Jun 15 2000
    copyright            : (C) 2000 by Edimar Cesar Rylo
    email                : ecrylo@fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef ELGQ2DCYL
#define ELGQ2DCYL

#include "pzreal.h"
#include "pzelgq2d.h"
#include "pzcosys.h"

class TPZGeoElQ2dCyl : public TPZGeoElQ2d {

public:
  /**Constructors. Parameters: id - element id, nodeindexes - vector containing node indexes,
     matind - material index, refind - index of the node which indicates the reference direction*/
  TPZGeoElQ2dCyl(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh, int cosysindex);
  TPZGeoElQ2dCyl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh, int cosysindex);

  //Destructor
  ~TPZGeoElQ2dCyl();
  TPZGeoEl* CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh);


private :
  REAL fRadius[4];
  REAL fTheta[4];
  REAL fZ[4];
  int fCosysIndex;
public:
  TPZCosys *fCosys;
  void Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);
  void X(TPZVec<REAL> & loc,TPZVec<REAL> &result);
  void NormalVector(int side,TPZVec<REAL> &param,TPZVec<REAL> &normal,TPZFMatrix &axes,TPZFMatrix &jac1d);
  void VerifyTheta();
  void Print(ostream & out);
};

#endif

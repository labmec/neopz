/**************************************************************************
                          tpzmatplacaiso.h  -  description
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

#ifndef TPZMATPLACAISO_H
#define TPZMATPLACAISO_H

#include "pzmatplaca2.h"

/**
  *@author Edimar Cesar Rylo
  */

class TPZMatPlacaIso : public TPZMatPlaca2  {
public: 
	TPZMatPlacaIso(int num, REAL h, REAL f, REAL E1 , REAL E2 ,
           REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
           REAL G23 , TPZFMatrix &naxes, TPZVec<REAL> &xf);

	~TPZMatPlacaIso();
  /** Solver */
  void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
  /**  Contribute*/
  void Contribute(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &, double weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);


};

#endif

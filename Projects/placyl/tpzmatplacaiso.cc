/***************************************************************************
                          tpzmatplacaiso.cpp  -  description
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

#include "tpzmatplacaiso.h"

TPZMatPlacaIso::TPZMatPlacaIso(int num, REAL h, REAL f, REAL E1 , REAL E2 ,
           REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
           REAL G23 , TPZFMatrix &naxes, TPZVec<REAL> &xf):TPZMatPlaca2(num,h,f,E1,E2,ni1,ni2,G12,G13,
           G23 , naxes ,xf){
}


TPZMatPlacaIso::~TPZMatPlacaIso(){
}


/**Contribute */
void TPZMatPlacaIso::Contribute(TPZVec<REAL> &x , TPZVec<REAL> &sol , TPZFMatrix &mat , double weight , TPZFMatrix &axes , TPZFMatrix &phi , TPZFMatrix &dphi , TPZFMatrix &ek , TPZFMatrix &ef){
	
  SetNAxes(axes);
  TPZMatPlaca2::Contribute(x,sol,mat,weight,axes,phi,dphi,ek,ef);
}


/** Solver */
void TPZMatPlacaIso::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){

  SetNAxes(axes);
  TPZMatPlaca2::Solution(Sol,DSol,axes,var,Solout);
}

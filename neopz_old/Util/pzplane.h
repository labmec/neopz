/***************************************************************************
                          TPZPlane.h  -  description
                             -------------------
    begin                : Wed Mar 27 2002
    copyright            : (C) 2002 by Renato Gomes Damas
    email                : rgdamas@fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TPZPLANE_H
#define TPZPLANE_H

#include <cmath>
#include <iostream>
#include "pzvec.h"
#include "pznumeric.h"

/**Esta classe cria e guarda a equação de um plano.
Verifica também se determinado elemento pertence ao plano
  *@author Renato Gomes Damas
  */

class TPZPlane {

 public:
  
  TPZPlane();
  ~TPZPlane();
  /** Dado três pontos calcula a equação do plano que os contém. */
  int SetPlane(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2, const TPZVec<REAL> &p3);
  /** Verifica se o ponto[3] pertence ao plano. Se pertencer retorna 1, caso contrário 0.*/
  bool Belongs(const TPZVec<REAL> &ponto);
  /** Verifica se o plano coincide com plano formado pelos três pontos passados. Se pertencer retorna 1, caso contrário 0. */
  bool Belongs(const TPZVec<REAL> &ponto1, const TPZVec<REAL> &ponto2, const TPZVec<REAL> &ponto3);
  
 private:
  
  /** Coeficientes da equação do plano: fCoef[0]*x + fCoef[1]*y + fCoef[2]*z + fCoef[3] = 0. */
  TPZVec<REAL> fCoef;

};

#endif

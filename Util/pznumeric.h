/***************************************************************************
                          pznumeric.h  -  description
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

#ifndef TPZNUMERIC_H
#define TPZNUMERIC_H


/**Tem implementados vários metodos de cálculo.
  *@author Renato Gomes Damas
  */
#include "pzvec.h"
using namespace std;

class TPZNumeric {

public:
	
	TPZNumeric();
	~TPZNumeric();
  	/** Retorna o determinante da matriz em &det. */
  	static void MatrixDet(REAL matrix[3][3], REAL &det);
  	/** Retorna o determinante da matriz. */
  	static REAL MatrixDet(REAL matrix[3][3]);
  	/** Dada a array[3] armazena sua ordem decrescente, em valor absoluto, em ordem[3]. */
  	static void SortArray3(const TPZVec<REAL> &array, int ordem[3]);
  	/** Dada a array[3], retorna-a em ordem decrescente, em valor absoluto. */
  	static void SortArray3(TPZVec<REAL> &array);
  	/** dados dois vetores calcula o produto vetorial. */
  	static void ProdVetorial(TPZVec<REAL> &u, TPZVec<REAL> &v, TPZVec<REAL> &result);

	static void NormalizeVetor3(TPZVec<REAL> &vetor);

    static REAL Norma(const TPZVec<REAL> &vetor);
};

#endif

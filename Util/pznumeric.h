/**
 * @file
 * @brief Contains declaration of the TPZNumeric class which implements several methods to calculation.
 */
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

#include "pzvec.h"

/** 
 * @ingroup util
 * @brief Implements several methods to calculation. \ref util "Utility"
 * @author Renato Gomes Damas
 */
class TPZNumeric {
	
public:
	
	TPZNumeric();
	~TPZNumeric();
  	/** Retorna o determinante da matriz em &det. */
	/** @brief Compute the 3x3-matrix determinant */
  	static void MatrixDet(REAL matrix[3][3], REAL &det);
  	/** Retorna o determinante da matriz. */
	/** @brief Returns the 3x3-matrix determinant */
  	static REAL MatrixDet(REAL matrix[3][3]);
  	/** Dada a array[3] armazena sua ordem decrescente, em valor absoluto, em ordem[3]. */
	/** @brief Sorts in descending order, in absolute value and stores the indexes in ordem. */
  	static void SortArray3(const TPZVec<REAL> &array, int ordem[3]);
  	/** Dada a array[3], retorna-a em ordem decrescente, em valor absoluto. */
	/** @brief Sorts in descending order, in absolute value on self storage vector */
  	static void SortArray3(TPZVec<REAL> &array);
  	/** dados dois vetores calcula o produto vetorial. */
	/** @brief Computes the vectorial product u x v */
  	static void ProdVetorial(TPZVec<REAL> &u, TPZVec<REAL> &v, TPZVec<REAL> &result);
	/** @brief Normalizes the vector */
	static void NormalizeVetor3(TPZVec<REAL> &vetor);
	/** @brief Returns the L2-norm of the vector */
    static REAL Norma(const TPZVec<REAL> &vetor);
};

#endif

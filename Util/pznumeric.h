/**
 * @file
 * @brief Contains declaration of the TPZNumeric class which implements several methods to calculation.
 */

#ifndef TPZNUMERIC_H
#define TPZNUMERIC_H

#include "pzvec.h"
#include "fad.h"

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
    template <class Tvar>
  	static void SortArray3(const TPZVec<Tvar> &array, int ordem[3]);
  	/** Dada a array[3], retorna-a em ordem decrescente, em valor absoluto. */
	/** @brief Sorts in descending order, in absolute value on self storage vector */
    template <class Tvar>
  	static void SortArray3(TPZVec<Tvar> &array);
  	/** dados dois vetores calcula o produto vetorial. */
	/** @brief Computes the vectorial product u x v */
    template <class Tvar>
  	static void ProdVetorial(TPZVec<Tvar> &u, TPZVec<Tvar> &v, TPZVec<Tvar> &result);
	/** @brief Normalizes the vector with 3 elements */
    template <class Tvar>
	static void NormalizeVetor3(TPZVec<Tvar> &vetor);
    /** @brief Normalizes the vector */
    template <class Tvar>
    static void NormalizeVetor(TPZVec<Tvar> &vetor);
	/** @brief Returns the L2-norm of the vector */
    template <class Tvar>
    static Tvar Norm(const TPZVec<Tvar> &vetor);
};

#endif

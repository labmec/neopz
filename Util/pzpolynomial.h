/**
 * @file
 * @brief Contains declaration of the TPZPolynomial class which implements a polynomial.
 */
#ifndef TPZPOLYNOMIAL_H
#define TPZPOLYNOMIAL_H
#include <iostream>
#include <cmath>
#include "pzvec.h"
#include "pznumeric.h"

/**
 * @ingroup util
 * @brief Implements a polynomial. \ref util "Utility"
 */
class TPZPolynomial {
public:
	/** @brief Default constructor */
    TPZPolynomial();
	
	/** @brief Constructor based on coef as coefficients and tol as tolerance value */
    TPZPolynomial(const TPZVec<REAL> &coef, const REAL &tol);
	
	/** @brief Constructor based on coef as coefficients */
    TPZPolynomial(const TPZVec<REAL> &coef);
	
	/** @brief Sets the tolerance value to computes */
    void SetTolerance(const REAL &tol);
	
	/** @brief Gets the tolerance value */
    void GetTolerance(REAL & tol);
	
	/** @brief Sets up four coefficients to polynomial */
    void SetCoef(const REAL &c0, const REAL &c1, const REAL &c2, const REAL &c3);
	
	/** @brief Sets the coefficients of the polynomial */
    void SetCoef(const TPZVec<REAL> &coef);
	
	/** @brief On given coefficients computes the roots of the polynomial in r */ 
    int GetRoots(const TPZVec<REAL> &coef, TPZVec<REAL> &r);
	
	/** @brief Computes the roots of the polynomial in r */
    int GetRoots(TPZVec<REAL> &r);
	
	/** @brief Computes the roots of the cubic polynomial using Tartaglia method (until 3 degree) */
    int Tartaglia(const TPZVec<REAL> &coef, TPZVec<REAL> &real, REAL &imagem);
	
private:
	/** @brief Computes the roots of the polynomial and stores into the fReal */
    int SetRoots();
	
	/** @brief Tolerance value to computes */
    REAL fTolerance;
	
	/** @brief Polynomial coefficients */
    TPZVec <REAL> fCo;
	
	/** @brief Roots of the polynomial */
    TPZVec <REAL> fReal;
    REAL fImagem;
};

#endif

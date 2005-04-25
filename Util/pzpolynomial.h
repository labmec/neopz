#ifndef TPZPOLYNOMIAL_H
#define TPZPOLYNOMIAL_H
#include <iostream>
#include <cmath>
#include "pzvec.h"
#include "pznumeric.h"
class TPZPolynomial {
public:
    /** Construtor da classe. */
    TPZPolynomial();

    /** Construtor da classe. */
    TPZPolynomial(const TPZVec<REAL> &coef, const REAL &tol);

    /** Construtor da classe.* */
    TPZPolynomial(const TPZVec<REAL> &coef);

    /** Define a tolerância dos cálculos. */
    void SetTolerance(const REAL &tol);

    /** Fornece a tolerância do cálculo, armazenando em tol. */
    void GetTolerance(REAL & tol);

    /** Armazena os coeficientes do polinomio em fCo. */
    void SetCoef(const REAL &c0, const REAL &c1, const REAL &c2, const REAL &c3);

    /** Armazena os coeficientes do polinomio em fCo. */
    void SetCoef(const TPZVec<REAL> &coef);

    /** Dados os coeficientes do polinomio, retorna as raízes do polinomio em r. */
    int GetRoots(const TPZVec<REAL> &coef, TPZVec<REAL> &r);

    /** Retorna as raízes do polinomio em r. */
    int GetRoots(TPZVec<REAL> &r);

    int Tartaglia(const TPZVec<REAL> &coef, TPZVec<REAL> &real, REAL &imagem);

private:
    /** Calcula as raízes do polinomio e armazena em r. */
    int SetRoots();

    /** Tolerância para o cálculo. */
    REAL fTolerance;

    /** Coeficientes do polinômio. */
    TPZVec <REAL> fCo;

    /** Raízes do polinômio. */

    TPZVec <REAL> fReal;
    REAL fImagem;
};
#endif

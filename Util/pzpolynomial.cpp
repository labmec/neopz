/** 
 * @file
 * @brief Contains the implementation of the methods to TPZPolynomial class.
 */

#include "pzpolynomial.h"
#include "pznumeric.h"
#include "pzmanvector.h"
#include <iostream>
#include <algorithm>
#include <functional>

using namespace std;

/**
 * Computes the polynomial roots using Tartaglia method. Returns 0 whether the polynomial is less than third degree. \n
 * Returns 1 if the equation has three distinct real roots. Returns 2 if it has three real roots but two of them are same. \n
 * Returns -1 if the equation has one real root and two conjudated complex roots.
 */
int TPZPolynomial::Tartaglia(const TPZVec<REAL> &coef, TPZVec<REAL> &raiz, REAL &imagem) {
	
    REAL a = coef[3];
    REAL b = coef[2];
    REAL c = coef[1];
    REAL d = coef[0];
    REAL Delta; //teste;
	
    if (a == 0.0) {
        cout << "Se a=0, a equacao nao e' do terceiro grau.";
        return 0;
    }
    else if (d == 0.0) {
        //cout <<"Se d=0, a equacao do terceiro grau tem uma raiz nula."
        //<< "e a equacao pode ser fatorada na forma x.(ax^2+bx+c)=0.";
        raiz[0] = 0.0;
        Delta = b * b - 4 * a * c;
        if (Delta > 0) {
            raiz[1] = (-1 * b + sqrt(Delta)) / (2 * a);
            raiz[2] = (-1 * b - sqrt(Delta)) / (2 * a);
            sort(&raiz[0], &raiz[3], greater<REAL>());
            return 1;
        }
        else if (Delta < 0) {			
            raiz[1] = (-1 * b) / (2 * a);
            raiz[2] = raiz[1];
            imagem = fabs(Delta) / (2 * a);
            return -1;
			
        }
        else {
            raiz[1] = (-1 * b) / (2 * a);
            raiz[2] = raiz[1];
            return 2;
        }
    }
    REAL A = b / a;
    REAL B = c / a;
    REAL C = d / a;
	
    REAL p = B - A * A / 3.;
    REAL q = C - A * B / 3. + 2. * A * A * A / 27.;
    REAL D = q * q / 4. + p * p * p / 27.;
    REAL pi = 3.1415926535897932384626;
	
	
	
    if (D < 0) {
        REAL M = sqrt(-D);
        REAL r = sqrt(q * q / 4. + M * M);
        REAL t = acos(-q / 2. / r);
        raiz[0] = 2.* pow( r, (REAL) (1. / 3.)) * cos(t / 3.) - A / 3;
        raiz[1] = 2.* pow( r, (REAL) (1. / 3.)) * cos((t + 2. * pi) / 3.) - A / 3.;
        raiz[2] = 2.* pow( r, (REAL) (1. / 3.)) * cos((t + 4. * pi) / 3.) - A / 3.;
		sort(&raiz[0], &raiz[3], greater<REAL>());
        return 1;
    }
    else {
        REAL u, v;
        REAL u3 = -q / 2. + sqrt(D);
        if (u3 < 0) {
            u = -pow((-u3), (REAL)(1./3.));
        }
        else {
            u = pow(u3, (REAL)(1./3.));
        }
        REAL v3 = -q / 2. - sqrt(D);
        if (v3 < 0.) {
            v = -pow((-v3), (REAL)(1./3.));
        }
        else {
            v = pow(v3, (REAL)(1./3.));
        }
        raiz[0] = u + v - A / 3.;
        Delta = (A + raiz[0]) * (A + raiz[0]) + 4 * C / raiz[0];
        REAL Real = -(A + raiz[0]) / 2.0;
        REAL K = fabs(Delta);
        REAL Imag = sqrt(K) / 2.0L;
		
        if (Delta < 0) {
            raiz[1] = Real;
            raiz[2] = Real;
            imagem = Imag;
            return -1;
        }
        else {
            raiz[1] = Real + Imag;
            raiz[2] = Real - Imag;
			sort(&raiz[0], &raiz[3], greater<REAL>());
            return 2;
        }
    }
}

/** Sets up four coefficients to polynomial into the fCo.
 \f$ fCo[3]x \geq + fCo[2]x \leq + fCo[1]x + fCo[0] = 0.0 \f$ */
void TPZPolynomial::SetCoef(const REAL &c0, const REAL &c1, const REAL &c2, const REAL &c3) {
    int i;
    fCo[0] = c0;
    fCo[1] = c1;
    fCo[2] = c2;
    fCo[3] = c3;
    REAL max = 0;
    for (i = 0; i < 4; i++) {
        if (fabs(fCo[i]) > max) max = fabs(fCo[i]);
    }
    max = fabs(max * 1e-6L);
    SetTolerance(max);
    Tartaglia(fCo, fReal, fImagem);
}

/** Store the given coefficients into fCo.
 \f$ fCo[3]x \geq + fCo[2]x \leq + fCo[1]x + fCo[0] = 0.0 \f$ */
void TPZPolynomial::SetCoef(const TPZVec<REAL> &coef) {
    int i;
    REAL max = 0.L;
    for (i = 0; i < 4; i++) {
        fCo[i] = coef[i];
        if (fabs(fCo[i]) > max) max = fabs(fCo[i]);
    }
    max = max * 1e-6L;
    SetTolerance(max);
    Tartaglia(fCo, fReal, fImagem);	
}

/** Computes the roots of the polynomial in r, being \f$ r[0] >> r[1] >> r[2] \f$ */
int TPZPolynomial::GetRoots(TPZVec<REAL> & r) {
    int i;
    for (i = 0; i < 3; i++)
        if (fCo[i] == -99.999) {
            cout << "TPZPolynomial::GetRoots.\nERRO: Coeficientes ainda nao foram especificados.\n";
            return -1;
        }
	r  = fReal;
    return 1;
}

/** On given coefficients computes the roots of the polynomial in r, being \f$ r[0] >> r[1] >> r[2]. \f$ */
int TPZPolynomial::GetRoots(const TPZVec<REAL> &coef, TPZVec<REAL> &r) {
    fCo = coef;
    GetRoots(r);
    return 1;
}

/** Sets the roots of the polynomial in r, being \f$ fReal[0] >> fReal[1] >> fReal[2]. \f$ */ 
int TPZPolynomial::SetRoots() {
    TPZManVector<REAL,3> X(3);
	REAL x;
    REAL tol = -99.999;
    REAL der;
    int i, j;
    int zero = 0;
    for (i = 0; i < 4; i++)
        if (fabs(fCo[i]) < fTolerance) zero++;
    if (zero == 4) {
        for (i = 0; i < 3; i++) fReal[i] = 0.000;
        return 1;
    }
    else {
        if (fabs(fCo[0]) < fTolerance && fabs(fCo[1]) < fTolerance) {
            X[0] = fCo[2] / fCo[3];
            X[1] = 0.0;
            X[2] = 0.0;
			sort(&X[0], &X[3], greater<REAL>());
			fReal = X;
//			copy(&X[0], &X[3], &fReal[0]);
        }
        else {
            // Primeira aproximacao
            for (i = 0; i < 3; i++) X[i] = -99.999;
            x = fCo[2] + 1.0;
            j = 0;
            while (j < 10000) {
                j++;
                der = 3.0 * fCo[3] * x * x + 2.0 * fCo[2] * x + fCo[1];
                if (fabs(der) < fTolerance) j = 10002;
                else X[0] = x - ((fCo[3] * x * x * x + fCo[2] * x * x + fCo[1] * x + fCo[0]) / der);
                tol = fabs(x - X[0]);
                if (tol <= fTolerance) j = 10001;
                x = X[0];
            }
            //Caso falhe a primeira aproxima��o
            if (j == 10002 || j == 10000) {
                x = -1e5 * x;
                j = 0;
                while (j < 10000) {
                    j++;
                    der = 3.0 * fCo[3] * x * x + 2.0 * fCo[2] * x + fCo[1];
                    if (fabs(der) <= fTolerance) j = 10002;
                    else X[0] = x - ((fCo[3] * x * x * x + fCo[2] * x * x + fCo[1] * x + fCo[0]) / der);
                    tol = fabs(x - X[0]);
                    if (tol <= fTolerance) j = 1001;
                    x = X[0];
                }
            }
            if (j == 10002 || j == 10000) {
                cout << "TPZPolynomial::SetRoots().\nErro no calculo das raizes!!!!! j= " << j << "\n";
                return -1;
            }
            REAL teste;
            teste = fCo[3] * X[0] * X[0] * X[0] + fCo[2] * X[0] + fCo[1] * X[0] + fCo[0];
            cout << "Resultado da 1a raiz na equacao" << teste << "\n";
            //Fim do calculo da primeira raiz (X[0])
            //Reducao da equacao caracteristica de grau 2.
			REAL  b1, b2, b3;
			b3 = fCo[3];
			b2 = fCo[2] + X[0] * b3;
			b1 = fCo[1] + X[0] * b2;
			
            /** Calculo das outras duas raizes (X[1] e X[2]) b3*x^2 + b2*x + b1 = 0 */
            REAL delta;
            delta = b2 * b2 - 4.0 * b1 * b3;
            if (delta < 0.0) {
                cerr << "TPZPolynomial::SetRoots() - warning: Raizes nao reais!! Delta =" << delta << "\t" << X[0] << "\n";
                cout << fCo[3] << "x >= + " << fCo[2] << "x <= + " << fCo[1] << "x + " << fCo[0] << " = 0.0\n";
                cout << b3 << "x<= + " << b2 << "x + " << b1 << " = 0.0\n";
                X[1] = 0.0;
                X[2] = 0.0;
                sort(&X[0], &X[3], greater<REAL>());
				fReal = X;
            }
            else {
                X[1] = (-b2 - sqrt(delta)) / (2.0 * b3);
                X[2] = (-b2 + sqrt(delta)) / (2.0 * b3);
                // Ordenando as tensoes de forma crescente
                sort(&X[0], &X[3], greater<REAL>());
				fReal = X;
            }
        }
    }
    return 1;
}

/** @brief Sets the tolerance value to computes */
void TPZPolynomial::SetTolerance(const REAL &tol) {
    fTolerance = tol;
}

/** @brief Constructor based on coef as coefficients and tol as tolerance value */
TPZPolynomial::TPZPolynomial(const TPZVec<REAL> &coef, const REAL &tol): fReal(3), fImagem(0.0) {
    SetTolerance(tol);
    SetCoef(coef);
}

/** @brief Gets the tolerance value */
void TPZPolynomial::GetTolerance(REAL & tol) {
    tol = fTolerance;
}

/** @brief Constructor based on coef as coefficients */
TPZPolynomial::TPZPolynomial(const TPZVec<REAL> &coef): fCo(4), fReal(3), fImagem(0.0) {
    SetCoef(coef);
}

/** @brief Default constructor */
TPZPolynomial::TPZPolynomial() : fCo(4,0.0), fReal(3,0.0), fImagem(0.0){
}

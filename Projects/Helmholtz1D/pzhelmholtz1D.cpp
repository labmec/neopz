/*
 * @file Helmhotz1D.cpp
 * @brief Contains the implementation of the TPZHelmholtz1D methods.
 */

#include "pzhelmholtz1D.h"

TPZHelmholtz1D::TPZHelmholtz1D(int id,int dim) : TPZMat1dLin(id) {
}

TPZHelmholtz1D::TPZHelmholtz1D(const TPZHelmholtz1D &helm) : TPZMat1dLin(helm) {
}

TPZHelmholtz1D::~TPZHelmholtz1D() {
}

void TPZHelmholtz1D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {

    TPZManVector<REAL,2> alphaval(2), betaval(2), phiaval(2);
    fAlpha->Execute(data.x, alphaval);
    fBeta->Execute(data.x, betaval);
    fPhi->Execute(data.x, phiaval);
    TPZFNMatrix<4,REAL> xk(2,2),xb(2,2),xc(2,2,0.),xf(2,1);
    xk(0,0) = alphaval[0];
    xk(0,1) = -alphaval[1];
    xk(1,0) = alphaval[1];
    xk(1,1) = alphaval[0];
    xb(0,0) = betaval[0];
    xb(0,1) = -betaval[1];
    xb(1,0) = betaval[1];
    xb(1,1) = betaval[0];
    xf(0,0) = phiaval[0];
    xf(1,0) = phiaval[1];
    
    
  
    // atualizar xf, adicionar uma funcao para tal.
    
    SetMaterial(xk, xc, xb, xf);
    TPZMat1dLin::Contribute(data, weight, ek, ef);
}

/*
void TPZHelmholtz1D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc) {

    
}
*/
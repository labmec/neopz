/*
 * @file Helmhotz1D.cpp
 * @brief Contains the implementation of the TPZHelmholtz1D methods.
 */

#include "pzhelmholtz1D.h"

TPZHelmholtz1D::TPZHelmholtz1D(int id,int dim) : TPZMaterial(1) {
}

TPZHelmholtz1D::TPZHelmholtz1D(const TPZHelmholtz1D &helm) {
}

TPZHelmholtz1D::~TPZHelmholtz1D() {
}

void TPZHelmholtz1D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {

}

void TPZHelmholtz1D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc) {

}

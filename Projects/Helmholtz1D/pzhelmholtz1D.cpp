/*
 * @file Helmhotz1D.cpp
 * @brief Contains the implementation of the TPZHelmholtz1D methods.
 */

#include "Helmhotz1D.h"

TPZHelmholtz1D::TPZHelmholtz1D(int id,int dim) : TPZMaterial(1) {
}

TPZHelmholtz1D::TPZHelmholtz1D(TPZHelmholtz1D &helm) {
}

TPZHelmholtz1D::~TPZHelmholtz1D() {
}

virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {

}

virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc) {

}

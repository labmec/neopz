/**
 * \file
 * @brief Contains implementations of the TPZDiscontinuousGalerkin methods.
 */
// $Id: pzdiscgal.cpp,v 1.12 2010-10-08 19:18:15 fortiago Exp $

#include "pzdiscgal.h"
#include "pzmaterialdata.h"
#include "pzmaterialid.h"

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin() : TPZMaterial(){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(int nummat) : TPZMaterial(nummat){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy) : TPZMaterial(copy) {}

TPZDiscontinuousGalerkin::~TPZDiscontinuousGalerkin(){}

std::string TPZDiscontinuousGalerkin::Name() { return "TPZDiscontinuousGalerkin"; }

void TPZDiscontinuousGalerkin::FillDataRequirementsInterface(TPZMaterialData &data){
	data.SetAllRequirements(true);
	data.fNeedsSol = false;
	if(fLinearContext == false){
		data.fNeedsNeighborSol = true;
	}
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                                   REAL weight, TPZFMatrix &ef){
	TPZFMatrix fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeInterface(data, dataleft, dataright, weight, fakeek, ef);
}

void TPZDiscontinuousGalerkin::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ef,TPZBndCond &bc){
	TPZFMatrix fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeBCInterface(data, dataleft, weight, fakeek, ef, bc);
}

int TPZDiscontinuousGalerkin::IsInterfaceConservative(){
	return 0;
}

void TPZDiscontinuousGalerkin::InterfaceJump(TPZVec<REAL> &x, 
											 TPZSolVec &leftu,
											 TPZSolVec &rightu,
											 TPZSolVec &jump){
    int numbersol = leftu.size();
    for (int is=0; is<numbersol; is++) {
        const int n = leftu[is].NElements();
        jump[is].Resize(n);
        for(int i = 0; i < n; i++){
            jump[is][i] = leftu[is][i] - rightu[is][i];
        }
    }
}

void TPZDiscontinuousGalerkin::BCInterfaceJump(TPZVec<REAL> &x, 
                                               TPZSolVec &leftu,
                                               TPZBndCond &bc,
                                               TPZSolVec & jump){
	PZError << __PRETTY_FUNCTION__ << " - method not implemented in derived class" << std::endl;
	DebugStop();
}

int TPZDiscontinuousGalerkin::ClassId() const{
	return TPZDISCONTINUOUSGALERKIN;
}

void TPZDiscontinuousGalerkin::Write(TPZStream &buf, int withclassid){
	TPZMaterial::Write(buf, withclassid);
}

void TPZDiscontinuousGalerkin::Read(TPZStream &buf, void *context){
	TPZMaterial::Read(buf, context);
}

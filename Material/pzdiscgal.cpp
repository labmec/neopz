/**
 * @file
 * @brief Contains implementations of the TPZDiscontinuousGalerkin methods.
 */

#include "pzdiscgal.h"
#include "pzmaterialdata.h"
#include <algorithm>



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
												   REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    std::cout << __PRETTY_FUNCTION__ << " please implement me\n";
    DebugStop();
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                                   REAL weight, TPZFMatrix<STATE> &ef){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeInterface(data, dataleft, dataright, weight, fakeek, ef);
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, 
                                                   REAL weight, TPZFMatrix<STATE> &ef){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeInterface(data, dataleft, dataright, weight, fakeek, ef);
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeInterface(data, dataleft, dataright, weight, fakeek, ef);

}
void TPZDiscontinuousGalerkin::ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                                 REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    std::cout << __PRETTY_FUNCTION__ << " please implement me\n";
    DebugStop();    
}

void TPZDiscontinuousGalerkin::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
	this->ContributeBCInterface(data, dataleft, weight, fakeek, ef, bc);
}



void TPZDiscontinuousGalerkin::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
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

void TPZDiscontinuousGalerkin::Errors(TPZMaterialData &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    TPZMaterial::Errors(data,u_exact,du_exact,errors);
}

void TPZDiscontinuousGalerkin::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    int nspace = data.size();
    int ispace = 0;
    for(; ispace <nspace; ispace++)
    {
        if(data[ispace].fShapeType != TPZMaterialData::EEmpty) break;
    }
    TPZMaterial::Errors(data[ispace],u_exact,du_exact,errors);
}


void TPZDiscontinuousGalerkin::BCInterfaceJump(TPZVec<REAL> &x, 
                                               TPZSolVec &leftu,
                                               TPZBndCond &bc,
                                               TPZSolVec & jump){
	PZError << __PRETTY_FUNCTION__ << " - method not implemented in derived class" << std::endl;
	DebugStop();
}

int TPZDiscontinuousGalerkin::ClassId() const{
    return Hash("TPZDiscontinuousGalerkin") ^ TPZMaterial::ClassId() << 1;
}

void TPZDiscontinuousGalerkin::Write(TPZStream &buf, int withclassid) const{
	TPZMaterial::Write(buf, withclassid);
}

void TPZDiscontinuousGalerkin::Read(TPZStream &buf, void *context){
	TPZMaterial::Read(buf, context);
}

/// return the integration order as a function of interpolation orders of the left and right elements
int TPZDiscontinuousGalerkin::GetIntegrationOrder(TPZVec<int> &porder_left, TPZVec<int> &porder_right) const
{
    int maxl = 0, maxr = 0;
    for (auto porder: porder_left) {
        maxl = std::max(maxl,porder);
    }
    for (auto porder: porder_right) {
        maxr = std::max(maxr,porder);
    }
    return maxl+maxr;
}

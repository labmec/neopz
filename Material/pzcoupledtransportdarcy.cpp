/**
 * \file
 * @brief Contains implementations of the TPZCoupledTransportDarcy methods.
 */

//$Id: pzcoupledtransportdarcy.cpp,v 1.8 2009-09-01 19:44:46 phil Exp $

#include "pzcoupledtransportdarcy.h"
#include "pzcoupledtransportdarcyBC.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmanvector.h"
#include <math.h>

using namespace std;
int TPZCoupledTransportDarcy::gCurrentEq = 0;

void TPZCoupledTransportDarcy::SetCurrentMaterial(const int i){
	if (i == 0 || i == 1) TPZCoupledTransportDarcy::gCurrentEq = i;
	else PZError << "Error! - " << __PRETTY_FUNCTION__ << endl;
}

int TPZCoupledTransportDarcy::CurrentEquation(){ return TPZCoupledTransportDarcy::gCurrentEq; }

TPZCoupledTransportDarcy::TPZCoupledTransportDarcy(int nummat, int nummat0, int nummat1, int dim) : 
TPZDiscontinuousGalerkin(nummat), fAlpha(1.) {
	this->fMaterials[0] = new TPZMatPoisson3d(nummat0, dim);
	fMaterialRefs[0] = fMaterials[0];
	this->fMaterials[1] = new TPZMatPoisson3d(nummat1, dim);
	fMaterialRefs[1] = fMaterials[1];
}

void TPZCoupledTransportDarcy::SetAlpha(REAL alpha){
	this->fAlpha = alpha;
}

TPZCoupledTransportDarcy::~TPZCoupledTransportDarcy() {
}

int TPZCoupledTransportDarcy::NStateVariables() {
	return this->GetCurrentMaterial()->NStateVariables();
}

void TPZCoupledTransportDarcy::Print(ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Base Class properties : \n";
	TPZMaterial::Print(out);
}

void TPZCoupledTransportDarcy::ContributeInterface(TPZMaterialData &data,
												   REAL weight,
												   TPZFMatrix &ek,
												   TPZFMatrix &ef){
	this->UpdateConvectionDirInterface(data.dsoll, data.dsolr);
	this->GetCurrentMaterial()->ContributeInterface(data, weight, ek, ef);
}

void TPZCoupledTransportDarcy::ContributeBCInterface(TPZMaterialData &data,
													 REAL weight, 
													 TPZFMatrix &ek,
													 TPZFMatrix &ef,
													 TPZBndCond &bc){
	this->UpdateConvectionDir(data.dsoll); 
	this->GetCurrentMaterial()->ContributeBCInterface( data, weight, ek, ef, bc );
}

void TPZCoupledTransportDarcy::Contribute(TPZMaterialData &data,
                                          REAL weight,
                                          TPZFMatrix &ek,TPZFMatrix &ef) {
	this->UpdateConvectionDir(data.dsol);
	this->GetCurrentMaterial()->Contribute(data, weight, ek, ef);
}


void TPZCoupledTransportDarcy::ContributeBC(TPZMaterialData &data,
                                            REAL weight,
                                            TPZFMatrix &ek,TPZFMatrix &ef,
                                            TPZBndCond &bc) {
	this->GetCurrentMaterial()->ContributeBC(data, weight, ek, ef, bc);
}

/** returns the variable index associated with the name*/
int TPZCoupledTransportDarcy::VariableIndex(const std::string &name){
	return this->GetCurrentMaterial()->VariableIndex(name);
}

int TPZCoupledTransportDarcy::NSolutionVariables(int var){
	return this->GetCurrentMaterial()->NSolutionVariables(var);
}

void TPZCoupledTransportDarcy::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,
										int var,TPZVec<REAL> &Solout){
	TPZMaterialData data;
	data.sol = Sol;
	data.dsol = DSol;
	data.axes = axes;
	return this->GetCurrentMaterial()->Solution(data, var, Solout);
}//method

void TPZCoupledTransportDarcy::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
	//Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

void TPZCoupledTransportDarcy::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
									  TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
									  TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {
	this->GetCurrentMaterial()->Errors(x, u, dudx, axes, flux, u_exact, du_exact, values);
}


TPZCoupledTransportDarcyBC * TPZCoupledTransportDarcy::CreateBC(int id){
	return new TPZCoupledTransportDarcyBC(this, id);
}

void TPZCoupledTransportDarcy::UpdateConvectionDir(TPZFMatrix &dsol){
	if (TPZCoupledTransportDarcy::CurrentEquation() == 1){
		//It is necessary to set Beta1 = alpha * (-K0 Grad[p] )
		REAL K, C;
		const int dim = this->Dimension();
		TPZManVector<REAL, 3> dir(dim);
		this->GetMaterial(0)->GetParameters(K, C, dir);
		const REAL K0 = K;
		this->GetMaterial(1)->GetParameters(K, C, dir);    
		
		int i;
		for(i = 0; i < dim; i++){
			dir[i] = dsol(i,0);
			dir[i] *= -1. * K0 * this->fAlpha;
		}
		
		this->GetMaterial(1)->SetParameters(K, 1., dir);
	}
}

void TPZCoupledTransportDarcy::UpdateConvectionDirInterface(TPZFMatrix &dsolL, TPZFMatrix &dsolR){
	if (TPZCoupledTransportDarcy::CurrentEquation() == 1){
		int i, j;
		int nrows = dsolL.Rows();
		int ncols = dsolL.Cols();
		TPZFNMatrix<100> dsol(nrows, ncols);
		for(i = 0; i < nrows; i++) for(j = 0; j < ncols; j++) dsol(i,j) = 0.5 * ( dsolL(i,j) + dsolR(i,j) );
		this->UpdateConvectionDir(dsol);
	}
}

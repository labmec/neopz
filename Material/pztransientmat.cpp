/**
 * \file
 * @brief Contains implementations of the TPZTransientMaterial methods.
 */

//$Id: pztransientmat.cpp,v 1.8 2009-05-06 20:22:18 fortiago Exp $

#include "pztransientmat.h"

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::SetExplicit(){
	this->fTemporalIntegrator = EExplicit;
}
template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::SetImplicit(){
	this->fTemporalIntegrator = EImplicit;
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::TPZTransientMaterial(int nummat, int dim, REAL TimeStep):TBASEMAT(nummat, dim){
	this->SetTimeStep(TimeStep);
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::TPZTransientMaterial(const TPZTransientMaterial &cp):TBASEMAT(cp){
	this->fTemporalIntegrator = cp.fTemporalIntegrator;
	this->fStep = cp.fStep;
	this->fTimeStep = cp.fTimeStep;
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::~TPZTransientMaterial(){
	//NOTHING TO BE DONE
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::Contribute(TPZMaterialData &data,
                                                  REAL weight,
                                                  TPZFMatrix<REAL> &ek,
                                                  TPZFMatrix<REAL> &ef){
	
	// Mostly for implicit
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	if (this->fStep == ECurrent){
		TBASEMAT::Contribute(data,weight,ek,ef);
		this->ContributeSolutionRhs(data.sol[0], data.phi, weight, ef);
		this->ContributeTangent(data.sol[0], data.phi, weight, ek);
		return;
	}
	
	if (this->fStep == ELast){
		this->ContributeSolutionRhs(data.sol[0], data.phi, weight, ef);
		return;
	}
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		this->ContributeTangent(data.sol[0], data.phi, weight, ek);
		return;
	}
	
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TBASEMAT::Contribute(data,weight,ek,ef);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeBC(TPZMaterialData &data,
                                                    REAL weight,
                                                    TPZFMatrix<REAL> &ek,
                                                    TPZFMatrix<REAL> &ef,
                                                    TPZBndCond &bc){
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeBC(data,weight,ek,ef,bc);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		TPZFNMatrix<1000> fakeef(ek.Rows(),1,0.);
		TBASEMAT::ContributeBC(data,weight,ek,fakeef,bc);
		return;
	}
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TPZFNMatrix<1000> fakeef(ef.Rows(),ef.Rows(),0.);
		TBASEMAT::ContributeBC(data,weight,ek,ef,bc);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                                           REAL weight,
                                                           TPZFMatrix<REAL> &ek,
                                                           TPZFMatrix<REAL> &ef){
	
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeInterface(data,dataleft,dataright, weight, ek, ef);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		return;
	}
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TBASEMAT::ContributeInterface(data,dataleft,dataright, weight, ek, ef);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                                             REAL weight, 
                                                             TPZFMatrix<REAL> &ek,
                                                             TPZFMatrix<REAL> &ef,
                                                             TPZBndCond &bc){
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeBCInterface(data,dataleft, weight,ek, ef, bc);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		return;
	}
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TBASEMAT::ContributeBCInterface(data, dataleft, weight,  ek, ef, bc);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeSolutionRhs(TPZVec<REAL> &sol, TPZFMatrix<REAL> &phi, REAL weight, TPZFMatrix<REAL> &ef){
	REAL Mult = +1.; 
	//Last solution is added to residual
	if (this->fStep == ECurrent) Mult = -1.; 
	//Current solution is subtracted from residual
	const int phr = phi.Rows();
	const int nstate = this->NStateVariables();
	const REAL DeltaT = this->TimeStep();
	int i, k;
	for(i = 0; i < phr; i++) {
		for(k = 0; k < nstate; k++){
			ef(i*nstate+k, 0) += Mult * weight * sol[k] * phi(i,0) / DeltaT;
		}//k
	}//i
}//method

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeTangent(TPZVec<REAL> &sol, TPZFMatrix<REAL> &phi, REAL weight, TPZFMatrix<REAL> &ek){
	const int phr = phi.Rows();
	const int nstate = this->NStateVariables();
	const REAL DeltaT = this->TimeStep();
	int i, j, k;
	for(i = 0; i < phr; i++) {
		for(j = 0; j < phr; j++){
			for(k = 0; k < nstate; k++){
				ek(i*nstate+k, j*nstate+k) += weight * phi(i,0) * phi(j,0) / DeltaT;
			}//k
		}//j
	}//i
}//method

#include "pzpoisson3d.h"
template class TPZTransientMaterial< TPZMatPoisson3d >;

#include "pznonlinearpoisson3d.h"
template class TPZTransientMaterial< TPZNonLinearPoisson3d >;

#include "pzburger.h"
template class TPZTransientMaterial< TPZBurger >;

/** @brief Instantiations to TPZMatPoisson3d, TPZNonLinearPoisson3d and TPZBurger. */
void TestInstantiations(){
	TPZTransientMaterial< TPZMatPoisson3d > A(1,1,1.);
	TPZTransientMaterial< TPZNonLinearPoisson3d > B(1,1,1.);
	TPZTransientMaterial< TPZBurger > C(1,1,1.);
}



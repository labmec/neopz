/**
 * \file
 * @brief Contains implementations of the TPZL2Projection methods.
 */
//$Id: pzl2projection.cpp,v 1.13 2011-05-26 03:28:57 phil Exp $

#include "pzl2projection.h"
#include "pzbndcond.h"
#include "tpzintpoints.h"

TPZL2Projection::TPZL2Projection(int id, int dim, int nstate, TPZVec<REAL> &sol,
                                 int IntegrationOrder) :TPZDiscontinuousGalerkin(id) , fScale(1.)
{
    this->fDim = dim;
    this->fNStateVars = nstate;
    this->fSol = sol;
    this->fIntegrationOrder = IntegrationOrder;
    this->SetIsReferred(false);
}

TPZL2Projection::TPZL2Projection(const TPZL2Projection &cp):TPZDiscontinuousGalerkin(cp), fScale(cp.fScale)
{
    this->fDim = cp.fDim;
    this->fNStateVars = cp.fNStateVars;
    this->fSol = cp.fSol;
    this->fIntegrationOrder = cp.fIntegrationOrder;
    this->SetIsReferred(cp.fIsReferred);
}

TPZL2Projection::~TPZL2Projection()
{
}

void TPZL2Projection::SetIsReferred(bool val){
	this->fIsReferred = val;
}

TPZAutoPointer<TPZMaterial> TPZL2Projection::NewMaterial(){
	return new TPZL2Projection(*this);
}

void TPZL2Projection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
	
	if (this->HasForcingFunction()){
		this->fForcingFunction->Execute(data.x, this->fSol);
	}
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	const int nvars = this->fNStateVars;
	if (this->fIsReferred){
		this->fSol.Resize(nvars);
		if (data.sol.NElements() < 2*nvars){//it means the referred element does not exist or it is an interface element. In that case, I ASSUME the referred solution is ZERO
			this->fSol.Fill(0.);
		}//if
		else{
			for(int i = 0; i < nvars; i++){
				this->fSol[i] = data.sol[0][i+nvars];
			}//for
		}//else
	}//if
	
	const int nshape = data.phi.Rows();
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
			for(int ivi = 0; ivi < nvars; ivi++){
				for(int ivj = 0; ivj < nvars; ivj++){
					const int posI = nvars*i+ivi;
					const int posJ = nvars*j+ivj;
					ek(posI, posJ) += weight*fScale*data.phi(i,0)*data.phi(j,0);
				}//ivj
			}//ivi
		}//for j
		for(int ivi = 0; ivi < nvars; ivi++){
			const int posI = nvars*i+ivi;
			ef(posI,0) += weight*fScale*data.phi(i,0)*this->fSol[ivi];
		}//ivi
	}//for i
}

void TPZL2Projection::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
	
	const int nvars = this->fNStateVars;
	TPZFMatrix<REAL> &phi = data.phi;
	const int phr = phi.Rows();
	int in, jn, iv;
	
	switch (bc.Type()){
			
			// Dirichlet condition
		case 0 : {      
			for(iv = 0; iv < nvars; iv++){
				for(in = 0 ; in < phr; in++) {
					ef(nvars*in+iv,0) += TPZMaterial::gBigNumber * bc.Val2()(iv,0) * phi(in,0) * weight;      
					for (jn = 0 ; jn < phr; jn++) {
						ek(nvars*in+iv,nvars*jn+iv) += TPZMaterial::gBigNumber * phi(in,0) * phi(jn,0) * weight;
					}//jn
				}//in
			}//iv
			break;
		}
			
			// Neumann condition
		case 1 : {
			for(iv = 0; iv < nvars; iv++){
				for(in = 0 ; in < phr; in++) {
					ef(nvars*in+iv,0) += bc.Val2()(iv,0) * fScale * phi(in,0) * weight;
				}//in
			}//iv
			break;
		}
			
		default:{
			std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
		}
	}//switch
	
}

int TPZL2Projection::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str())) return ESolution;
	return TPZMaterial::VariableIndex(name);
}

int TPZL2Projection::NSolutionVariables(int var){
	const int nvars = this->NStateVariables();
	if(var == ESolution) return nvars;
	
	return 0;
}

void TPZL2Projection::Solution(TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol,
                               TPZFMatrix<REAL> &axes, int var, TPZVec<REAL> &Solout){
	if (var == ESolution){
		Solout = Sol;
		return;
	}
	
	Solout.Resize(0);
}

/** 
 * Get the order of the integration rule necessary to integrate an
 * element with polinomial order p
 */
int TPZL2Projection::IntegrationRuleOrder(int elPMaxOrder) const
{
    if (this->fIntegrationOrder == -1) {
        return TPZDiscontinuousGalerkin::IntegrationRuleOrder(elPMaxOrder);
    }
    else
    {
        int order = (fIntegrationOrder > (2*elPMaxOrder) ) ? fIntegrationOrder : 2*elPMaxOrder;
        return order;
    }
}

/*
 void TPZL2Projection::SetIntegrationRule(TPZAutoPointer<TPZIntPoints> rule,
 int elPMaxOrder,
 int elDimension){
 if(this->fIntegrationOrder == -1){
 TPZDiscontinuousGalerkin::SetIntegrationRule(rule,elPMaxOrder,elDimension);
 }
 else{
 const int order = (fIntegrationOrder > (2*elPMaxOrder) ) ? fIntegrationOrder : 2*elPMaxOrder;
 TPZManVector<int,3> p2(elDimension,order);
 rule->SetOrder(p2);
 }
 }
 */

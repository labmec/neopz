/**
 * \file
 * @brief Contains implementations of the TPZL2Projection methods.
 */

#include "pzl2projection.h"
#include "pzbndcond.h"
#include "tpzintpoints.h"

TPZL2Projection::TPZL2Projection(int id, int dim, int nstate, TPZVec<STATE> &sol,
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

TPZMaterial * TPZL2Projection::NewMaterial(){
	return new TPZL2Projection(*this);
}

void TPZL2Projection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype==data.EVecShape){
        ContributeVecShape(data,weight,ek, ef);
        return;
    }
    
    TPZManVector<STATE> solloc(fSol);
	if (this->HasForcingFunction()){
        solloc.Resize(fForcingFunction->NFunctions());
		this->fForcingFunction->Execute(data.x, solloc);
	}
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	const int nvars = this->fNStateVars;
	if (this->fIsReferred){
		solloc.Resize(nvars);
		if (data.sol[0].NElements() < 2*nvars){//it means the referred element does not exist or it is an interface element. In that case, I ASSUME the referred solution is ZERO
			solloc.Fill(0.);
		}//if
		else{
			for(int i = 0; i < nvars; i++){
				solloc[i] = data.sol[0][i+nvars];
			}//for
		}//else
	}//if
	
	const int nshape = data.phi.Rows();
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
			for(int ivi = 0; ivi < nvars; ivi++){
                const int posI = nvars*i+ivi;
                const int posJ = nvars*j+ivi;
                ek(posI, posJ) += weight*fScale*data.phi(i,0)*data.phi(j,0);
			}//ivi
		}//for j
		for(int ivi = 0; ivi < nvars; ivi++){
			const int posI = nvars*i+ivi;
			ef(posI,0) += (STATE)weight*(STATE)fScale*(STATE)data.phi(i,0)*solloc[ivi];
		}//ivi
	}//for i
}

void TPZL2Projection::ContributeVecShape(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZManVector<STATE> solloc(fSol);
	if (this->HasForcingFunction()){
        solloc.Resize(fForcingFunction->NFunctions());
		this->fForcingFunction->Execute(data.x, solloc);
	}
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	const int nvariables = this->fNStateVars;
    const int vecFuncSize = data.phi.Rows();
    const int nshape = data.phi.Cols();
    if(nvariables != vecFuncSize)
    {
        //Nao implementado neste material ainda
        DebugStop();
    }
    
	if (this->fIsReferred){
        //Nao implementado neste material ainda
        DebugStop();
	}//if
	
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
			for(int ivi = 0; ivi < nvariables; ivi++){
                ek(i,j) += weight*fScale*data.phi(ivi,i)*data.phi(ivi,j);
			}//ivi
		}//for j
		for(int ivi = 0; ivi < nvariables; ivi++){
			ef(i,0) += (STATE)weight*(STATE)fScale*(STATE)data.phi(ivi,i)*solloc[ivi];
		}//ivi
	}//for i
}

void TPZL2Projection::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	
	const int nvars = this->fNStateVars;
	TPZFMatrix<REAL> &phi = data.phi;
	const int phr = phi.Rows();
	int in, jn, iv;
	
	switch (bc.Type()){
			
			// Dirichlet condition
		case 0 : {      
			for(iv = 0; iv < nvars; iv++){
				for(in = 0 ; in < phr; in++) {
					ef(nvars*in+iv,0) += (STATE)TPZMaterial::gBigNumber * bc.Val2()(iv,0) * (STATE)phi(in,0) * (STATE)weight;      
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
					ef(nvars*in+iv,0) += bc.Val2()(iv,0) * (STATE)fScale * (STATE)phi(in,0) * (STATE)weight;
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

void TPZL2Projection::Solution(TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol,
                               TPZFMatrix<REAL> &axes, int var, TPZVec<STATE> &Solout){
	if (var == ESolution){
#ifndef STATE_COMPLEX
        Solout.Resize(Sol.size());
        for (int i=0; i<Sol.size(); i++) {
            Solout[i] = Sol[i];
        }
#else
        
#endif
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

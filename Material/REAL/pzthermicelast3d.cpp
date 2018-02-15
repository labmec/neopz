/**
 * @file
 * @brief Contains implementations of the TPZThermicElast3D methods.
 */

#include "pzthermicelast3d.h"

TPZThermicElast3D::TPZThermicElast3D(int nummat, STATE ThermalCoeff, STATE RefTemp, STATE E, STATE poisson, TPZVec<STATE> &force):
TPZRegisterClassId(&TPZThermicElast3D::ClassId),
TPZElasticity3D(nummat,E,poisson,force){  
	this->SetReferredTemperatureField();  
	this->fThermalCoeff = ThermalCoeff;
	this->fRefTemperature = RefTemp;  
}

TPZThermicElast3D::~TPZThermicElast3D(){}

void TPZThermicElast3D::Contribute(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,
								   TPZFMatrix<STATE> &ef){
	
	TPZElasticity3D::Contribute(data, weight, ek, ef);
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	this->ContributeThermalStress(data.sol[0], data.phi, data.dphix, weight, ef);
}

void TPZThermicElast3D::ContributeThermalStress(TPZVec<STATE> &sol, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, REAL weight, TPZFMatrix<STATE> &ef){
	
	const int nshape = phi.Rows();  
	const STATE E  = this->fE;
	const STATE nu = this->fPoisson;
	const STATE RefTemp = this->fRefTemperature;
	
	STATE FinalTemp;
	if (this->IsReferredTemperatureField()) FinalTemp = sol[0];
	else FinalTemp = this->fFinalTemperature;
	
	const STATE DeltaT = FinalTemp - RefTemp; 
	const STATE Gamma = this->fThermalCoeff;
	
	const STATE ThermalStress = DeltaT * E * Gamma / (1.-2.*nu);
	
	int i, k;
	for(i = 0; i < nshape; i++) {
		for(k = 0; k < 3; k++){
			ef(i*3+k, 0) += weight * dphi(k, i) * ThermalStress;
		}//kd
	}//in
	
}//method

void TPZThermicElast3D::Solution(TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol,
                                 TPZFMatrix<REAL> &axes, int var, TPZVec<STATE> &Solout){
	
	if (var == EPrincipalStress || var == EVonMisesStress || var == EStress || var == EStress1){
		//When computing a stress solution, the thermal stress must be subtracted from the elastic stress
		//For that purpose the thermal strain will be subtracted before stress computations
		
		const STATE RefTemp = this->fRefTemperature;  
		STATE FinalTemp;
		if (this->IsReferredTemperatureField()){
			FinalTemp = Sol[3];
		}
		else FinalTemp = this->fFinalTemperature;  
		
		const STATE DeltaT = FinalTemp - RefTemp; 
		const STATE Gamma = this->fThermalCoeff;
		
		const STATE ThermalStrain = DeltaT * Gamma;
		for(int it = 0; it < 3; it++) DSol(it,it) += -1. * ThermalStrain;
	}//if var
    
	TPZElasticity3D::Solution(Sol, DSol, axes, var, Solout);
	
}

int TPZThermicElast3D::ClassId() const{
    return Hash("TPZThermicElast3D") ^ TPZElasticity3D::ClassId() << 1;
}

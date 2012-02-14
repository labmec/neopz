/**
 * \file
 * @brief Contains implementations of the TPZThermicElast3D methods.
 */
//$Id: pzthermicelast3d.cpp,v 1.2 2007-05-11 19:15:18 joao Exp $

#include "pzthermicelast3d.h"

TPZThermicElast3D::TPZThermicElast3D(int nummat, REAL ThermalCoeff, REAL RefTemp, REAL E, REAL poisson, TPZVec<REAL> &force)
:TPZElasticity3D(nummat,E,poisson,force){  
	this->SetReferredTemperatureField();  
	this->fThermalCoeff = ThermalCoeff;
	this->fRefTemperature = RefTemp;  
}

TPZThermicElast3D::~TPZThermicElast3D(){}

void TPZThermicElast3D::Contribute(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix &ek,
								   TPZFMatrix &ef){
	
	TPZElasticity3D::Contribute(data, weight, ek, ef);
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	this->ContributeThermalStress(data.sol[0], data.phi, data.dphix, weight, ef);
}

void TPZThermicElast3D::ContributeThermalStress(TPZVec<REAL> &sol, TPZFMatrix &phi, TPZFMatrix &dphi, REAL weight, TPZFMatrix &ef){
	
	const int nshape = phi.Rows();  
	const REAL E  = this->fE;
	const REAL nu = this->fPoisson;
	const REAL RefTemp = this->fRefTemperature;
	
	REAL FinalTemp;
	if (this->IsReferredTemperatureField()) FinalTemp = sol[0];
	else FinalTemp = this->fFinalTemperature;
	
	const REAL DeltaT = FinalTemp - RefTemp; 
	const REAL Gamma = this->fThermalCoeff;
	
	const REAL ThermalStress = DeltaT * E * Gamma / (1.-2.*nu);
	
	int i, k;
	for(i = 0; i < nshape; i++) {
		for(k = 0; k < 3; k++){
			ef(i*3+k, 0) += weight * dphi(k, i) * ThermalStress;
		}//kd
	}//in
	
}//method

void TPZThermicElast3D::Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
                                 TPZFMatrix &axes, int var, TPZVec<REAL> &Solout){
	
	if (var == EPrincipalStress || var == EVonMisesStress || var == EStress || var == EStress1){
		//When computing a stress solution, the thermal stress must be subtracted from the elastic stress
		//For that purpose the thermal strain will be subtracted before stress computations
		
		const REAL RefTemp = this->fRefTemperature;  
		REAL FinalTemp;
		if (this->IsReferredTemperatureField()){
			FinalTemp = Sol[3];
		}
		else FinalTemp = this->fFinalTemperature;  
		
		const REAL DeltaT = FinalTemp - RefTemp; 
		const REAL Gamma = this->fThermalCoeff;
		
		const REAL ThermalStrain = DeltaT * Gamma;
		for(int it = 0; it < 3; it++) DSol(it,it) += -1. * ThermalStrain;
	}//if var
    
	TPZElasticity3D::Solution(Sol, DSol, axes, var, Solout);
	
}


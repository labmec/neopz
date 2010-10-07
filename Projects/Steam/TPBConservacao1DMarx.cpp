/*
 *  TPBConservacao1D.cpp
 *  IP3D_v4
 *
 *  Created by Philippe Devloo on 27/03/10.
 *  Copyright 2010 UNICAMP. All rights reserved.
 *
 */

#include "TPBConservacao1DMarx.h"
#include "ThermalMethodsTables.h"


int main()
{
	TPBrCellMarx first;
	REAL MaterialPermeability = 1.e-12;
	TPZManVector<REAL> Viscosity(3,1.e-4);
	TPZManVector<REAL> Density(3,1000.);
	TPZManVector<REAL,3> SpecificHeat(3,1.);
	REAL TimeStep = 1.;
	
	REAL CellVolume = 1.;
	REAL LeftArea = 1.;
	REAL RightArea = 1.;
	REAL CellSize = 0.1;
	REAL Porosity = 0.8;
	REAL DensityRock = 3000.;
	
	TPZManVector<REAL> Pressure(3,1.e6);
	TPZManVector<REAL> Massflux(3,100.);
	TPZManVector<REAL> DarcyVelocity(3,1.);
	TPZManVector<REAL> Saturation(3,0.);
	Saturation[0] = 1.;
	TPZManVector<REAL,3> Energy(3,5.);
	REAL TotalEnergy(15.);
	REAL Temperature(100.);
	
	
	TPZManVector<REAL> initial(TPBrCellMarx::NUMVARS,0.), residual(TPBrCellMarx::NUMVARS,0.);
	TPZManVector<REAL> leftstate(TPBrCellMarx::NUMVARS,0.),rightstate(TPBrCellMarx::NUMVARS,0.);
	
	first.SetMaterialProperty(MaterialPermeability, Viscosity,SpecificHeat,TimeStep);
	first.SetGeometry(CellVolume,LeftArea,RightArea,CellSize,Porosity, DensityRock);
	first.SetState(Pressure,Massflux,DarcyVelocity,Density,Saturation,Energy,Temperature);
	initial[TPBrCellMarx::EPressureOil] = Pressure[TPBrCellMarx::EOil];
	initial[TPBrCellMarx::EPressureWater] = Pressure[TPBrCellMarx::EWater];
	initial[TPBrCellMarx::EPressureSteam] = Pressure[TPBrCellMarx::ESteam];
	initial[TPBrCellMarx::EDarcyVelocityOil] = DarcyVelocity[TPBrCellMarx::EOil];
	initial[TPBrCellMarx::EDarcyVelocityWater] = DarcyVelocity[TPBrCellMarx::EWater];
	initial[TPBrCellMarx::EDarcyVelocitySteam] = DarcyVelocity[TPBrCellMarx::ESteam];
	initial[TPBrCellMarx::EMassFluxOil] = Massflux[TPBrCellMarx::EOil];
	initial[TPBrCellMarx::EMassFluxWater] = Massflux[TPBrCellMarx::EWater];
	initial[TPBrCellMarx::EMassFluxSteam] = Massflux[TPBrCellMarx::ESteam];
	initial[TPBrCellMarx::ESaturationOil] = Saturation[TPBrCellMarx::EOil];
	initial[TPBrCellMarx::ESaturationWater] = Saturation[TPBrCellMarx::EWater];
	initial[TPBrCellMarx::ESaturationSteam] = Saturation[TPBrCellMarx::ESteam];
	initial[TPBrCellMarx::EEnthalpyOil] = Energy[TPBrCellMarx::EOil];
	initial[TPBrCellMarx::EEnthalpyWater] = Energy[TPBrCellMarx::EWater];
	initial[TPBrCellMarx::EEnthalpySteam] = Energy[TPBrCellMarx::ESteam];
	initial[TPBrCellMarx::ETotalEnergy] = TotalEnergy;
	initial[TPBrCellMarx::ETemperature] = Temperature;
	initial[TPBrCellMarx::EPhaseChange] = 0.;
	
	leftstate = initial;
	rightstate = initial;
	
	first.TotalResidual(leftstate,rightstate,residual);
	
	TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > tangent(TPBrCellMarx::NUMVARS,0.);
	first.TotalResidual(leftstate, rightstate, tangent);
	
	
	TPZFMatrix tangentmatrix,residualmatrix,statematrix;
	first.ExtractMatrix(tangent,tangentmatrix);
	first.ExtractMatrix(residual,residualmatrix);

	std::cout << "Residual " << residual << std::endl;
	std::cout << "Tangent " << tangent << std::endl;
	
	tangentmatrix.Print("tangent = ",cout, EMathematicaInput);
	
	return 0;
	
	while (Norm(residualmatrix) > 1.e-6) {
		tangentmatrix.SolveDirect(residualmatrix, ELU);
		statematrix -= residualmatrix;
		first.ConvertState(statematrix,rightstate);
		first.TotalResidual(leftstate, rightstate, residual);
		first.TotalResidual(leftstate, rightstate, tangent);
	}
	
	std::cout << "Residual " << residual << std::endl;
	std::cout << "Tangent " << tangent << std::endl;
	
	tangentmatrix.Print("tangent = ",cout, EMathematicaInput);
	return 0;
	
};

// metodos para recuperar os dados tabulados em funcao
template<class T>
//[kg/m3]
T TPBrCellMarx::DensityOil(T temp) {
	OilData oil;
	return oil.getDensityToOil(temp);
}
template<class T>
//[kg/m3]
T TPBrCellMarx::DensityWater(T temp) {
	WaterDataInStateOfSaturation t;
	return t.getSaturationStateDensityToLiquidWater(temp);
}
template<class T>
//[kg/m3]
T TPBrCellMarx::DensitySteam(T temp) {
	WaterDataInStateOfSaturation t;
	return t.getSaturationStateDensityToSteam(temp);
}

//[Pa*sec]
REAL TPBrCellMarx::ViscosityOil(REAL temp) {
	OilData oil;
	return oil.getDynamicViscosityToOil(temp);
}

//[Pa*sec]
REAL TPBrCellMarx::ViscosityWater(REAL temp) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateViscosityToLiquidWater(temp);
}

//[Pa*sec]
REAL TPBrCellMarx::ViscositySteam(REAL temp) {
	WaterDataInStateOfSaturation steam;
	return steam.getSaturationStateViscosityToSteam(temp);
}

template<class T>
//dimensionless
void TPBrCellMarx::ComputeRelativePermeability(TPZManVector<T> &saturation,TPZManVector<T> &relativepermeability) {
	relativepermeability = saturation;
}
template<class T>
T TPBrCellMarx::TemperatureSaturation(T p) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateTemperature(p);
}

template<class T>
//[kJoule/kg]
T TPBrCellMarx::EnthalpyWater(T temperature) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
}
template<class T>
//[kJoule/kg]
T TPBrCellMarx::EnthalpySteam(T temperature) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateSpecificEnthalpyToSteam(temperature);
}
template<class T>
//[kJoule/kg]
T TPBrCellMarx::EnthalpyOil(T temperature) {
	OilData oil;
#warning Verificar se specific heat e enthalply??
	return  oil.getSpecificHeatToOil(temperature);
}


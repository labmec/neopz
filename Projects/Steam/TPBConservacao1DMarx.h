#ifndef TPBCONSERVACAO1DMARX
#define TPBCONSERVACAO1DMARX

/*
 *  TPBConservacao1D.h
 *  IP3D_v4
 *
 *  Created by Philippe Devloo on 27/03/10.
 *  Copyright 2010 UNICAMP. All rights reserved.
 *
 */

#include "PropertiesTable.h"
#include "ThermalMethodsTables.h"

#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"

#ifdef _AUTODIFF
using namespace std;
#include "fadType.h"


//@TODO COLOCAR ENERGIA

/// implements a first order method for one dimensional fluid flow
class TPBrCellMarx
{
	// Passo de tempo
	
	// Célula
	/**
	 *Volume [m3]; Area [m2] e Size [m]
	 */
	REAL fCellVolume;
	REAL fLeftArea;
	REAL fRightArea;
	REAL fCellSize;
	
	// propriedades da rocha
	/**
	 *Porisity [%];
	 *SpecificHeat [kJ/(kg*Kelvin)]
	 *Density [kg/m3]
	 *ThermalConductivity [ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
	 */
	REAL fPorosityRock;
	REAL fSpecificHeatRock;
	REAL fDensityRock;
	
	REAL fResidualOil;
	
	REAL fMaterialPermeability;//[m2]
	
	
	TPZManVector<REAL> fInitialState;
public:
	enum VARNAME {EMassFluxOil,EMassFluxWater,EMassFluxSteam,EDarcyVelocityOil,EDarcyVelocityWater,EDarcyVelocitySteam,
		EPressureOil,EPressureWater,EPressureSteam,ESaturationOil,ESaturationWater,ESaturationSteam,
		EEnthalpyOil,EEnthalpyWater,EEnthalpySteam,ETotalEnergy,ETemperature,EPhaseChange, EDeltaT};
	
	//Original
	enum NUMBEREQUATIONS {EqMassConservationOil, EqMassConservationWater,EqMassConservationSteam, 
		EqRelationMassFlowVDarcyOil, EqRelationMassFlowVDarcyWater, EqRelationMassFlowVDarcySteam, 
		EqRelationDpDxDarcyOil, EqRelationDpDxDarcyWater, EqRelationDpDxDarcySteam, 
		EqCapilaryPressureOilWater, EqCapilaryPressureOilSteam, EqSumSaturation,
		EqEnergyOil, EqEnergyWater, EqEnergySteam,
		EqTotalEnergy, EqEnergyConservation, EqZeroSaturation, EqTemperSaturation};
	
	
	//Only Swap
	//enum NUMBEREQUATIONS {EqRelationMassFlowVDarcyOil, EqRelationMassFlowVDarcyWater, EqRelationMassFlowVDarcySteam,
	//		EqEnergyConservation,EqMassConservationOil,EqRelationDpDxDarcySteam, 
	//		EqCapilaryPressureOilWater, EqCapilaryPressureOilSteam, EqRelationDpDxDarcyWater,EqTotalEnergy,
	//		EqRelationDpDxDarcyOil,EqMassConservationWater,EqEnergyOil,EqEnergyWater, EqEnergySteam,EqSumSaturation,
	//		EqZeroSaturation,EqMassConservationSteam
	//	};
	
	
	
	//Swap and factor scale
	//enum NUMBEREQUATIONS {EqMassConservationOil, EqMassConservationWater,EqMassConservationSteam,
	//		EqRelationDpDxDarcyOil,EqRelationMassFlowVDarcyWater,EqRelationDpDxDarcySteam, 
	//		EqCapilaryPressureOilWater, EqCapilaryPressureOilSteam,EqRelationDpDxDarcyWater,
	//		EqRelationMassFlowVDarcyOil,EqSumSaturation,EqRelationMassFlowVDarcySteam, 
	//		EqEnergyOil, EqEnergyWater, EqEnergySteam,
	//		EqTotalEnergy, EqEnergyConservation, EqZeroSaturation
	//	};
	
	
	enum VARINDEX {EOil, EWater, ESteam};
	
	static const int NUMVARS = 19;
	
	TPBrCellMarx();
	
	virtual ~TPBrCellMarx();
	
	
	void SetMaterialProperty(REAL materialpermeability, PhysicalProperties &rock);
		
	void SetGeometry(REAL cellvolume, REAL leftarea, REAL rightarea, REAL cellsize);
		
	void SetInjectionState(REAL pressurewater, TPZVec<REAL> &massflux, TPZManVector<REAL> &leftstate);
		
	void SetCellState(REAL pressurewater, TPZVec<REAL> &saturation, REAL temperature, REAL deltIni );
		
	template<class T>
	//[kJ] = [1000 kg (m2/s2)]
	T EnergySolid(T temperature);
	
	template<class T>
	//[kJ] = [1000 kg (m2/s2)]
	T EnergySaturatedSteam(T pressuresteam);
	
	template<class T>
	//pressuresteam->Pa
	//[Celsius]
	T TemperatureSaturatedSteam(T pressuresteam);
	
	template<class T>
	void InitializeState(TPZManVector<T> &state);
	
	template<class T>
	void Inflow(TPZVec<REAL> &leftval, TPZVec<T> &flux, T timestep);
	
	template<class T>
	void Outflow(TPZVec<T> &rightval, TPZVec<T> &flux);
	
	void ReferenceResidualValues(TPZManVector<REAL> &state, TPZManVector<REAL> &scalevalues);
	void ReferenceStateValues(TPZManVector<REAL> &state, TPZManVector<REAL> &statescalevalues);
	
	//Retorna a razao entre o volume de agua expulsa e o volume da celula
	template<class T>
	void RazaoVexpVcel(TPZVec<T> &state, REAL &razao);		
	
	template<class T>
	void InternalEquations(TPZVec<T> &state, TPZVec<T> &residual);
	
	template<class T>
	void TotalResidual(TPZVec<REAL> &leftval, TPZVec<T> &state, TPZVec<T> &residual);
	
	static void ExtractMatrix(TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > &input, TPZFMatrix<REAL> &output);
	
	static void ExtractMatrix(TPZManVector<REAL> &input, TPZFMatrix<REAL> &output);
	
	template<class T>
	void ConvertState(TPZFMatrix<REAL> &input, TPZManVector<T> &output);
	
	// metodos para recuperar os dados tabulados em funcao
	template<class T>
	T DensityOil(T temp);//[kg/m3]
	template<class T>
	T DensityWater(T temp);//[kg/m3]
	template<class T>
	T DensitySteam(T temp);//[kg/m3]
	
	REAL ViscosityOil(REAL temp);//[Pa*sec]
	REAL ViscosityWater(REAL temp);//[Pa*sec]
	REAL ViscositySteam(REAL temp);//[Pa*sec]
	
	template<class T>
	void ComputeRelativePermeability(TPZManVector<T> &saturation,TPZManVector<T> &relativepermeability);//dimensionless
	
	template<class T>
	T TemperatureSaturation(T p);//[Celsius]
	template<class T>
	T EnthalpyWater(T temperature);//[kJ/kg]
	template<class T>
	T EnthalpySteam(T temperature);//[kJ/kg]
	template<class T>
	T EnthalpyOil(T temperature);//[kJ/kg]
	
};


template<class T>
inline void TPBrCellMarx::Inflow(TPZVec<REAL> &leftval, TPZVec<T> &flux, T timestep)
{		
	flux[EqMassConservationOil] = leftval[EMassFluxOil]*timestep;
	flux[EqMassConservationWater] = leftval[EMassFluxWater]*timestep;
	flux[EqMassConservationSteam] = leftval[EMassFluxSteam]*timestep;
	//( [kJ/kg]*[kg/s] )*[s] -> [kJ]
	flux[EqEnergyConservation] = (leftval[EEnthalpyOil]*leftval[EMassFluxOil]
								  +leftval[EEnthalpyWater]*leftval[EMassFluxWater]
								  +leftval[EEnthalpySteam]*leftval[EMassFluxSteam]
								  )*timestep;
}

template<class T>
inline void TPBrCellMarx::Outflow(TPZVec<T> &state, TPZVec<T> &flux)
{		
	//T QWaterSaida = state[EMassFluxWater]*state[EDeltaT];
	//cout << " \nQWaterSaida = "<< QWaterSaida<< endl;
	//cout << " \nVolumeCelula = "<< fCellVolume <<endl;
	
	flux[EqMassConservationOil] =state[EMassFluxOil]*state[EDeltaT];
	flux[EqMassConservationWater] = state[EMassFluxWater]*state[EDeltaT];
	flux[EqMassConservationSteam] = state[EMassFluxSteam]*state[EDeltaT];
	//( [kJ/kg]*[m/s]*[kg/m3] )*[m2]*[s] -> [kJ]
	flux[EqEnergyConservation] = (fInitialState[EEnthalpyOil]*state[EDarcyVelocityOil]*DensityOil(fInitialState[ETemperature])
								  +fInitialState[EEnthalpyWater]*state[EDarcyVelocityWater]*DensityWater(fInitialState[ETemperature])
								  +fInitialState[EEnthalpySteam]*state[EDarcyVelocitySteam]*DensitySteam(fInitialState[ETemperature])
								  )*fRightArea*state[EDeltaT];
}

inline void TPBrCellMarx::ReferenceResidualValues(TPZManVector<REAL> &state, TPZManVector<REAL> &scalevalues)
{
	REAL maxmassflux = max(state[EMassFluxOil],state[EMassFluxWater]);
	maxmassflux = max(state[EMassFluxSteam],maxmassflux);
	scalevalues[EqMassConservationOil] = maxmassflux *1000.;
	scalevalues[EqMassConservationWater] = maxmassflux*1000.;
	scalevalues[EqMassConservationSteam] = maxmassflux *1000.;
	
	scalevalues[EqRelationMassFlowVDarcyOil] = maxmassflux;
	scalevalues[EqRelationMassFlowVDarcyWater] = maxmassflux;
	scalevalues[EqRelationMassFlowVDarcySteam] = maxmassflux;
	
	REAL mindesity = max( DensityWater(state[ETemperature]), DensitySteam(state[ETemperature]) );
	REAL vDarcy = maxmassflux/(mindesity*fRightArea);
	scalevalues[EqRelationDpDxDarcyOil] = vDarcy*fCellSize;
	scalevalues[EqRelationDpDxDarcyWater] = vDarcy*fCellSize;
	scalevalues[EqRelationDpDxDarcySteam] = vDarcy*fCellSize;
	
	scalevalues[EqCapilaryPressureOilWater] = 1.;
	scalevalues[EqCapilaryPressureOilSteam] = 1.;
	scalevalues[EqSumSaturation] = 1.;
	
	EnthalpyOil(state[ETemperature]);
	EnthalpyWater(state[ETemperature]);
	EnthalpySteam(state[ETemperature]);
	scalevalues[EqEnergyOil] = EnthalpyOil(state[ETemperature]);
	scalevalues[EqEnergyWater] = EnthalpyWater(state[ETemperature]);
	scalevalues[EqEnergySteam] = EnthalpySteam(state[ETemperature]);
	
	scalevalues[EqTotalEnergy] = EnthalpyWater(state[ETemperature])*DensityWater(state[ETemperature])*fCellVolume
	+ (1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*fCellVolume*state[ETemperature];
	scalevalues[EqEnergyConservation] = EnthalpyWater(state[ETemperature])*DensityWater(state[ETemperature])*fCellVolume
	+ (1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*fCellVolume*state[ETemperature];
	
	scalevalues[EqZeroSaturation] = 1.;
	
	scalevalues[EqTemperSaturation] = 200.;
}

inline void TPBrCellMarx::ReferenceStateValues(TPZManVector<REAL> &state, TPZManVector<REAL> &statescalevalues)
{
	REAL maxmassflux = max(state[EMassFluxOil],state[EMassFluxWater]);
	maxmassflux = max(state[EMassFluxSteam],maxmassflux);
	
	statescalevalues[EMassFluxOil] = maxmassflux;
	statescalevalues[EMassFluxWater] = maxmassflux;
	statescalevalues[EMassFluxSteam] = maxmassflux;
	
	REAL pressure = max(state[EPressureOil],state[EPressureWater]);
	pressure = max(pressure,state[EPressureSteam]);
	statescalevalues[EPressureOil] = 100.;
	statescalevalues[EPressureSteam] = 100.;
	statescalevalues[EPressureWater] = 100.;
	
	REAL mindensity = max(DensityWater(state[ETemperature]), DensitySteam(state[ETemperature]));
	REAL vDarcy = maxmassflux/(mindensity*fRightArea);
	statescalevalues[EDarcyVelocityOil] = vDarcy;
	statescalevalues[EDarcyVelocitySteam] = vDarcy;
	statescalevalues[EDarcyVelocityWater] = vDarcy;
	
	statescalevalues[ESaturationOil] = 1.;	
	statescalevalues[ESaturationWater] = 1.;
	statescalevalues[ESaturationSteam] = 1.;
	
	statescalevalues[EEnthalpyOil] = state[EEnthalpyOil];
	statescalevalues[EEnthalpyWater] = state[EEnthalpyWater];
	statescalevalues[EEnthalpySteam] = state[EEnthalpySteam];
	
	statescalevalues[ETotalEnergy] = state[ETotalEnergy];
	
	statescalevalues[ETemperature] = state[ETemperature];
	
	REAL phasechange = maxmassflux;
	statescalevalues[EPhaseChange] = phasechange;
	
	statescalevalues[EDeltaT] = 1000.;
	
}

//--------------
// metodos para recuperar os dados tabulados em funcao
template<class T>
//[kg/m3]
inline T TPBrCellMarx::DensityOil(T temp) {
	OilData oil;
	return oil.getDensityToOil(temp);
}
template<class T>
//[kg/m3]
inline T TPBrCellMarx::DensityWater(T temp) {
	WaterDataInStateOfSaturation t;
	return t.getSaturationStateDensityToLiquidWater(temp);
}
template<class T>
//[kg/m3]
inline T TPBrCellMarx::DensitySteam(T temp) {
	WaterDataInStateOfSaturation t;
	return t.getSaturationStateDensityToSteam(temp);
}

//[Pa*sec]
inline REAL TPBrCellMarx::ViscosityOil(REAL temp) {
	OilData oil;
	return oil.getDynamicViscosityToOil(temp);
}

//[Pa*sec]
inline REAL TPBrCellMarx::ViscosityWater(REAL temp) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateViscosityToLiquidWater(temp);
}

//[Pa*sec]
inline REAL TPBrCellMarx::ViscositySteam(REAL temp) {
	WaterDataInStateOfSaturation steam;
	return steam.getSaturationStateViscosityToSteam(temp);
}

template<class T>
////dimensionless
inline void TPBrCellMarx::ComputeRelativePermeability(TPZManVector<T> &saturation,TPZManVector<T> &relativepermeability) {
	relativepermeability = saturation;
}

template<class T>
inline T TPBrCellMarx::TemperatureSaturation(T p) {
	WaterDataInStateOfSaturation water;
	T p1000 = p/1000.;
	return water.getSaturationStateTemperature(p1000);
}

template<class T>
//[kJoule/kg]
inline T TPBrCellMarx::EnthalpyWater(T temperature) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
}
template<class T>
//[kJoule/kg]
inline T TPBrCellMarx::EnthalpySteam(T temperature) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateSpecificEnthalpyToSteam(temperature);
}
template<class T>
//[kJoule/kg]
inline T TPBrCellMarx::EnthalpyOil(T temperature) {
	OilData oil;
// #warning Verificar se specific heat e enthalply??
	return  oil.getSpecificHeatToOil(temperature);
}
//----------------

template<class T>
inline void TPBrCellMarx::InternalEquations(TPZVec<T> &state, TPZVec<T> &residual)
{
	
	// conservation of mass of oil in the cell -> ESaturationOil
	// [m3]*[Kg/m3] -> [kg]
	residual[EqMassConservationOil] = fCellVolume*fPorosityRock*(DensityOil(state[ETemperature])*state[ESaturationOil]-DensityOil(fInitialState[ETemperature])*fInitialState[ESaturationOil]);
	
	// conservation of mass of water in the cell ->ESaturationWater
	// [m3]*[Kg/m3] -> [kg]
	residual[EqMassConservationWater] = fCellVolume*fPorosityRock*(DensityWater(state[ETemperature])*state[ESaturationWater]-DensityWater(fInitialState[ETemperature])*fInitialState[ESaturationWater])+state[EPhaseChange];
	
	// conservation of mass of steam in the cell ->ESaturationSteam
	// [m3]*[Kg/m3] -> [kg]
	//T  dS = DensitySteam(fInitialState[ETemperature]);
	residual[EqMassConservationSteam] = fCellVolume*fPorosityRock*(DensitySteam(state[ETemperature])*state[ESaturationSteam]-DensitySteam(fInitialState[ETemperature])*fInitialState[ESaturationSteam])-state[EPhaseChange];
	
	// relation between massflow rate and darcy velocity ->EMassFluxOil
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [kg/s]
	residual[EqRelationMassFlowVDarcyOil] = state[EMassFluxOil]-state[EDarcyVelocityOil]*DensityOil(fInitialState[ETemperature])*fRightArea;
	
	// relation between massflow rate and darcy velocity ->EMassFluxWater
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [gk/s]
	residual[EqRelationMassFlowVDarcyWater] = state[EMassFluxWater]-state[EDarcyVelocityWater]*DensityWater(fInitialState[ETemperature])*fRightArea;
	
	// relation between massflow rate and darcy velocity ->EMassFluxSteam
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [gk/s]
	residual[EqRelationMassFlowVDarcySteam] = state[EMassFluxSteam]-state[EDarcyVelocitySteam]*DensitySteam(fInitialState[ETemperature])*fRightArea;
	
	TPZManVector<T> saturation(3), relatpermeability(3);
	saturation[EOil]=state[ESaturationOil];
	saturation[EWater]=state[ESaturationWater];
	saturation[ESteam]=state[ESaturationSteam];
	ComputeRelativePermeability(saturation, relatpermeability);
	// relation between dp/dx and the darcy velocity ->EDarcyVelocityOil
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s] 
	residual[EqRelationDpDxDarcyOil] = (fInitialState[EPressureOil]-state[EPressureOil])*fMaterialPermeability*relatpermeability[EOil]/ViscosityOil(fInitialState[ETemperature])-state[EDarcyVelocityOil]*fCellSize;
	
	// relation between dp/dx and the darcy velocity -> EDarcyVelocityWater
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s]
	residual[EqRelationDpDxDarcyWater] = (fInitialState[EPressureWater]-state[EPressureWater])*fMaterialPermeability*relatpermeability[EWater]/ViscosityWater(fInitialState[ETemperature])-state[EDarcyVelocityWater]*fCellSize;
	
	// relation between dp/dx and the darcy velocity -> EDarcyVelocitySteam
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s]
	residual[EqRelationDpDxDarcySteam] = (fInitialState[EPressureSteam]-state[EPressureSteam])*fMaterialPermeability*relatpermeability[ESteam]/ViscositySteam(fInitialState[ETemperature])-state[EDarcyVelocitySteam]*fCellSize;
	
	// relation between the pressure of oil and water (capillary pressure) -> EPressureWater
	//[Pa] - [Pa] -> [Pa]
	if(state[ESaturationOil]<=fResidualOil) {
		residual[EqCapilaryPressureOilWater] = state[ESaturationOil]-fResidualOil;
	}
	else {
		residual[EqCapilaryPressureOilWater] = state[EPressureOil]-state[EPressureWater]; 
	}
	
	// relation between the pressure of oil and steam (capillary pressure) -> EPressureSteam
	//[Pa] - [Pa] -> [Pa]
	residual[EqCapilaryPressureOilSteam] = state[EPressureWater]-state[EPressureSteam];
	
	// the sum of saturations is equal 1 -> EPressureOil
	//dimensionless
	residual[EqSumSaturation] = state[ESaturationOil]+state[ESaturationWater]+state[ESaturationSteam] - 1.;
	
	// ENERGIA **************
	T TempSaturation = TemperatureSaturation(T(fInitialState[EPressureWater]));
	//([kJ/kg][kg/m3] + [kJ/kg][kg/m3])*[m3] - ([kJ/(kg*Kelvin)][Celsius][kg/m3])*(m3) ->...->[kJ]
	T EnergiaT = fCellVolume*fPorosityRock*(EnthalpyOil(TempSaturation)*DensityOil(TempSaturation)*state[ESaturationOil]
											+EnthalpyWater(TempSaturation)*DensityWater(TempSaturation)*(state[ESaturationWater]/*+state[ESaturationSteam]*/)
											+EnthalpySteam(TempSaturation)*DensitySteam(TempSaturation)*state[ESaturationSteam])
	+(1. - fPorosityRock)*fCellVolume*fSpecificHeatRock*fDensityRock*state[ETemperature];
	//cout << "\nEnergiaT --> " <<EnergiaT<<endl;
	
	// energy of the different fases : oil -> EEnthalpyOil
	//EEnthalpyOil ->[kJ/(kg*Kelvin)][Celsius] -> [kJ/kg]
	residual[EqEnergyOil] = state[EEnthalpyOil]-EnthalpyOil(state[ETemperature]); 
	
	// energy of the different fases : water -> EEnthalpyWater
	//EEnthalpyOil ->[kJ/(kg*Kelvin)][Celsius] -> [kJ/kg]
	residual[EqEnergyWater] = state[EEnthalpyWater]-EnthalpyWater(state[ETemperature]);
	
	// energy of the different fases : steam -> EEnthalpySteam
	//EEnthalpyOil ->[kJ/(kg*Kelvin)][Celsius] -> [kJ/kg]
	residual[EqEnergySteam] = state[EEnthalpySteam]-EnthalpySteam(state[ETemperature]);
	
	// total energy -> ETotalEnergy
	//ETotalEnergy -> ([kJ/kg]*[kg/m3])*[m3] - ([kJ/(kg*Kelvin)]*[kg/m3]*[m3]*[Celsius]) -> [kJ] 
	residual[EqTotalEnergy] = state[ETotalEnergy]-fCellVolume*fPorosityRock*(EnthalpyOil(state[ETemperature])*DensityOil(state[ETemperature])*state[ESaturationOil]
																			 +EnthalpyWater(state[ETemperature])*DensityWater(state[ETemperature])*state[ESaturationWater]
																			 +EnthalpySteam(state[ETemperature])*DensitySteam(state[ETemperature])*state[ESaturationSteam])
	-(1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*fCellVolume*state[ETemperature];
	
	// conservation of energy internal contribution -> ETemperature -> EPhaseChange
	//[kJ]
	residual[EqEnergyConservation] = (state[ETotalEnergy]-fInitialState[ETotalEnergy]);
	
	//T TotalEnergy = state[ETotalEnergy]; 
	//if(state[ETotalEnergy] <= EnergiaT/*EnergyCondensedSystem(state)*/){
	// Energy < EMin
	residual[EqZeroSaturation] = state[ESaturationSteam];
	//} else {
	//		residual[EqZeroSaturation] = state[ETemperature] - TempSaturation/*CondensationTemperature(state)*/;
	//		cout <<" Entrei ======="<<endl;
	//		cout << " Energia_min = "<< EnergiaT<<endl;
	//		
	//	}
	
	//equacao para obter o valor de DeltaT, para o qual a teperatura eh igual a temperatura de saturacao de vapor
	residual[EqTemperSaturation] = state[ETemperature] - TempSaturation;
	
}

template<class T>
inline void TPBrCellMarx::TotalResidual(TPZVec<REAL> &leftval, TPZVec<T> &state, TPZVec<T> &residual)
{
	TPZManVector<T> fluxl(NUMVARS,0.);
	TPZManVector<T> fluxr(NUMVARS,0.);
	TPZManVector<T> internal(NUMVARS,0.);
	Inflow(leftval, fluxl, state[EDeltaT]);
	Outflow(state, fluxr);
	InternalEquations(state, internal);
	int i;
	for (i=0; i<NUMVARS; i++) {
		residual[i] = internal[i]-(fluxl[i]-fluxr[i]);
	}
}

template<>
inline void TPBrCellMarx::InitializeState(TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > &state)
{		
	state.Resize(NUMVARS);
	state[EPressureOil] = fInitialState[EPressureOil];
	state[EPressureOil].fastAccessDx(EPressureOil) = 1.;
	state[EPressureWater] = fInitialState[EPressureWater];
	state[EPressureWater].fastAccessDx(EPressureWater) = 1.;
	state[EPressureSteam] = fInitialState[EPressureSteam];
	state[EPressureSteam].fastAccessDx(EPressureSteam) = 1.;
	state[EMassFluxOil] = fInitialState[EMassFluxOil];
	state[EMassFluxOil].fastAccessDx(EMassFluxOil) = 1.;
	state[EMassFluxWater] = fInitialState[EMassFluxWater];
	state[EMassFluxWater].fastAccessDx(EMassFluxWater) = 1.;
	state[EMassFluxSteam] = fInitialState[EMassFluxSteam];
	state[EMassFluxSteam].fastAccessDx(EMassFluxSteam) = 1.;
	state[EDarcyVelocityOil] = fInitialState[EDarcyVelocityOil];
	state[EDarcyVelocityOil].fastAccessDx(EDarcyVelocityOil) = 1.;
	state[EDarcyVelocityWater] = fInitialState[EDarcyVelocityWater];
	state[EDarcyVelocityWater].fastAccessDx(EDarcyVelocityWater) = 1.;
	state[EDarcyVelocitySteam] = fInitialState[EDarcyVelocitySteam];
	state[EDarcyVelocitySteam].fastAccessDx(EDarcyVelocitySteam) = 1.;
	state[ESaturationOil] = fInitialState[ESaturationOil];
	state[ESaturationOil].fastAccessDx(ESaturationOil) = 1.;
	state[ESaturationWater] = fInitialState[ESaturationWater];
	state[ESaturationWater].fastAccessDx(ESaturationWater) = 1.;
	state[ESaturationSteam] = fInitialState[ESaturationSteam];
	state[ESaturationSteam].fastAccessDx(ESaturationSteam) = 1.;
	state[EEnthalpyOil] = fInitialState[EEnthalpyOil];
	state[EEnthalpyOil].fastAccessDx(EEnthalpyOil) = 1.;
	state[EEnthalpyWater] = fInitialState[EEnthalpyWater];
	state[EEnthalpyWater].fastAccessDx(EEnthalpyWater) = 1.;
	state[EEnthalpySteam] = fInitialState[EEnthalpySteam];
	state[EEnthalpySteam].fastAccessDx(EEnthalpySteam) = 1.;
	state[ETotalEnergy] = fInitialState[ETotalEnergy];
	state[ETotalEnergy].fastAccessDx(ETotalEnergy) = 1.;
	state[ETemperature] = fInitialState[ETemperature];
	state[ETemperature].fastAccessDx(ETemperature) = 1.;
	state[EPhaseChange] = fInitialState[EPhaseChange];
	state[EPhaseChange].fastAccessDx(EPhaseChange) = 1.;
	state[EDeltaT] = fInitialState[EDeltaT];
	state[EDeltaT].fastAccessDx(EDeltaT) = 1.;
	
}

template<>
inline void TPBrCellMarx::InitializeState(TPZManVector<REAL> &state)
{
	state=fInitialState;
	
	//state.Resize(NUMVARS);
	//	state[EPressureOil] = 1.9755862983093257e6;
	//	state[EPressureWater] = 1.9755862983093257e6;
	//	state[EPressureSteam] = 1.9755862983093257e6;
	//	state[EMassFluxOil] = 0.18397689836721653;
	//	state[EMassFluxWater] = 0.45115809177932803;
	//	state[EMassFluxSteam] = 0.;
	//	state[EDarcyVelocityOil] = 0.00003337274897335301;
	//	state[EDarcyVelocityWater] = 0.00006651848185443374;
	//	state[EDarcyVelocitySteam] = 0.;
	//	state[ESaturationOil] = 0.017084997912007262;
	//	state[ESaturationWater] = 0.9829150020879928;
	//	state[ESaturationSteam] = 0.;
	//	state[EEnthalpyOil] = fInitialState[EEnthalpyOil];
	//	state[EEnthalpyWater] = fInitialState[EEnthalpyWater];
	//	state[EEnthalpySteam] = fInitialState[EEnthalpySteam];
	//	state[ETotalEnergy] = fInitialState[ETotalEnergy];
	//	state[ETemperature] = fInitialState[ETemperature];
	//	state[EPhaseChange] = -0.00555556;
}

template<>
inline void TPBrCellMarx::ConvertState(TPZFMatrix<REAL> &input, TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > &output)
{
	int nrows = input.Rows();
	int i;
	for (i=0; i<nrows; i++) {
		output[i] = input(i,0);
		output[i].fastAccessDx(i) = 1.;
	}
}
template<class T>
void TPBrCellMarx::ConvertState(TPZFMatrix<REAL> &input, TPZManVector<T> &output)
{
	InitializeState(output);
	int nrows = input.Rows();
	int i;
	for (i=0; i<nrows; i++) {
		output[i] = input(i,0);
	}
}

template<class T>
void  TPBrCellMarx::RazaoVexpVcel(TPZVec<T> &state, REAL &razao){
	T densidAgua = DensityWater(state[ETemperature]);
	T quantAguaExp = state[EMassFluxWater]*state[EDeltaT];
	REAL VexpAgua = quantAguaExp/densidAgua;
	REAL Vcel = fCellVolume*fPorosityRock;
	
	razao = VexpAgua/Vcel;
}

// nothing is compiled if _AUTODIFF isnt defined
#endif

#endif
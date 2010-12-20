/*
 *  TPBConservacao1D.h
 *  IP3D_v4
 *
 *  Created by Philippe Devloo on 27/03/10.
 *  Copyright 2010 UNICAMP. All rights reserved.
 *
 */
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"

//#ifdef _AUTODIFF
#include "fadType.h"
//#endif
#include "PropertiesTable.h"

#include <iostream>


//@TODO COLOCAR ENERGIA

/// implements a first order method for one dimensional fluid flow
class TPBrCellMarx
{
	// Passo de tempo
	REAL fTimeStep;
	
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
		EEnthalpyOil,EEnthalpyWater,EEnthalpySteam,ETotalEnergy,ETemperature,EPhaseChange};
	
	//Original
	enum NUMBEREQUATIONS {EqMassConservationOil, EqMassConservationWater,EqMassConservationSteam, 
		EqRelationMassFlowVDarcyOil, EqRelationMassFlowVDarcyWater, EqRelationMassFlowVDarcySteam, 
		EqRelationDpDxDarcyOil, EqRelationDpDxDarcyWater, EqRelationDpDxDarcySteam, 
		EqCapilaryPressureOilWater, EqCapilaryPressureOilSteam, EqSumSaturation,
		EqEnergyOil, EqEnergyWater, EqEnergySteam,
		EqTotalEnergy, EqEnergyConservation, EqZeroSaturation};
	
	
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
	
	static const int NUMVARS = 18;
	
	TPBrCellMarx();
	
	virtual ~TPBrCellMarx();
	
	void SetTimeStep(REAL &timestep);
	
	
	void SetMaterialProperty(REAL materialpermeability, REAL delt, PhysicalProperties &rock)
	{
		fMaterialPermeability = materialpermeability;
		fTimeStep = delt;
		fPorosityRock = rock.porosity;
		fDensityRock = rock.density;
		fSpecificHeatRock = rock.specificheat;
	}
	
	void SetGeometry(REAL cellvolume, REAL leftarea, REAL rightarea, REAL cellsize)
	{
		fCellVolume = cellvolume;
		fLeftArea = leftarea;
		fRightArea = rightarea;
		fCellSize = cellsize;
	}
	
	void SetInjectionState(REAL pressurewater, TPZVec<REAL> &massflux, TPZManVector<REAL> &leftstate)
	{
		REAL temperature = TemperatureSaturation (pressurewater);
		leftstate[ETemperature] = temperature;
		leftstate[EEnthalpyOil] = EnthalpyOil(temperature);
		leftstate[EEnthalpyWater] = EnthalpyWater(temperature);
		leftstate[EEnthalpySteam] = EnthalpySteam(temperature);
		leftstate[EMassFluxOil] = massflux[EOil];
		leftstate[EMassFluxWater] = massflux[EWater];
		leftstate[EMassFluxSteam] = massflux[ESteam];
		
	}
	
void SetCellState(REAL pressurewater, TPZVec<REAL> &saturation, REAL temperature ){
								
		fInitialState[EMassFluxOil] = 0.;
		fInitialState[EMassFluxWater] = 0.;
		fInitialState[EMassFluxSteam] = 0.;
		fInitialState[EDarcyVelocityOil] = 0.;
		fInitialState[EDarcyVelocityWater] = 0.;
		fInitialState[EDarcyVelocitySteam] = 0.;
		fInitialState[EPressureOil] = pressurewater;
		fInitialState[EPressureWater] = pressurewater;
		fInitialState[EPressureSteam] = pressurewater;
		fInitialState[ESaturationOil] = saturation[EOil];
		fInitialState[ESaturationWater] = saturation[EWater];
		fInitialState[ESaturationSteam] = saturation[ESteam];
		fInitialState[EEnthalpyOil] = EnthalpyOil(temperature);
		fInitialState[EEnthalpyWater] = EnthalpyWater(temperature);
		fInitialState[EEnthalpySteam] = EnthalpySteam(temperature);
		//([kJ/kg]*[kg/m3])*[m3] + [kJ/(kg*Kelvin)][kg/m3]*[m3]*[Celsius] -> [kJ]
		fInitialState[ETotalEnergy] = fPorosityRock*(EnthalpyOil(temperature)*DensityOil(temperature)*saturation[EOil] +
													 EnthalpyWater(temperature)*DensityWater(temperature)*saturation[EWater] +
													 EnthalpySteam(temperature)*DensitySteam(temperature)*saturation[ESteam])*fCellVolume 
		+ (1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*fCellVolume*temperature;
		fInitialState[ETemperature] = temperature;
	}
	
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
	void Inflow(TPZVec<REAL> &leftval, TPZVec<T> &flux);
	
	template<class T>
	void Outflow(TPZVec<T> &rightval, TPZVec<T> &flux);
	
	void ReferenceResidualValues(TPZManVector<REAL> &state, TPZManVector<REAL> &scalevalues);
	void ReferenceStateValues(TPZManVector<REAL> &state, TPZManVector<REAL> &statescalevalues);
	
	template<class T>
	void InternalEquations(TPZVec<T> &state, TPZVec<T> &residual);
	
	template<class T>
	void TotalResidual(TPZVec<REAL> &leftval, TPZVec<T> &state, TPZVec<T> &residual);
	
	static void ExtractMatrix(TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > &input, TPZFMatrix &output);
	
	static void ExtractMatrix(TPZManVector<REAL> &input, TPZFMatrix &output);
	
	template<class T>
	void ConvertState(TPZFMatrix &input, TPZManVector<T> &output);
	
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



TPBrCellMarx::TPBrCellMarx() : fInitialState(NUMVARS,0.)
{
	fInitialState[ESaturationOil] = 1.;
	fResidualOil = 0.01;
}

TPBrCellMarx::~TPBrCellMarx() {
}

void TPBrCellMarx::SetTimeStep(REAL &timestep)
{
	fTimeStep = timestep;
}

template<class T>
inline void TPBrCellMarx::Inflow(TPZVec<REAL> &leftval, TPZVec<T> &flux)
{		
	flux[EqMassConservationOil] = leftval[EMassFluxOil]*fTimeStep;
	flux[EqMassConservationWater] = leftval[EMassFluxWater]*fTimeStep;
	flux[EqMassConservationSteam] = leftval[EMassFluxSteam]*fTimeStep;
	//( [kJ/kg]*[kg/s] )*[s] -> [kJ]
	flux[EqEnergyConservation] = (leftval[EEnthalpyOil]*leftval[EMassFluxOil]
			   +leftval[EEnthalpyWater]*leftval[EMassFluxWater]
			   +leftval[EEnthalpySteam]*leftval[EMassFluxSteam]
			   )*fTimeStep;
	//cout << "\nFluxo de entrada --> "<< flux[EqEnergyConservation] <<endl; 
}

template<class T>
inline void TPBrCellMarx::Outflow(TPZVec<T> &state, TPZVec<T> &flux)
{		
	flux[EqMassConservationOil] =state[EMassFluxOil]*fTimeStep;
	flux[EqMassConservationWater] = state[EMassFluxWater]*fTimeStep;
	flux[EqMassConservationSteam] = state[EMassFluxSteam]*fTimeStep;
	//( [kJ/kg]*[m/s]*[kg/m3] )*[m2]*[s] -> [kJ]
	flux[EqEnergyConservation] = (fInitialState[EEnthalpyOil]*state[EDarcyVelocityOil]*DensityOil(fInitialState[ETemperature])
				+fInitialState[EEnthalpyWater]*state[EDarcyVelocityWater]*DensityWater(fInitialState[ETemperature])
				+fInitialState[EEnthalpySteam]*state[EDarcyVelocitySteam]*DensitySteam(fInitialState[ETemperature])
				)*fRightArea*fTimeStep;
	
	//cout << "\nFluxo de saida --> "<< flux[EqEnergyConservation] <<endl; 
	
}

inline void TPBrCellMarx::ReferenceResidualValues(TPZManVector<REAL> &state, TPZManVector<REAL> &scalevalues)
{
	REAL maxmassflux = max(state[EMassFluxOil],state[EMassFluxWater]);
	maxmassflux = max(state[EMassFluxSteam],maxmassflux);
	scalevalues[EqMassConservationOil] = maxmassflux *fTimeStep;
	scalevalues[EqMassConservationWater] = maxmassflux*fTimeStep;
	scalevalues[EqMassConservationSteam] = maxmassflux * fTimeStep;
	
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

	REAL mindensity = max( DensityWater(state[ETemperature]), DensitySteam(state[ETemperature]) );
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
	
}
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
	T Temp = TemperatureSaturation(fInitialState[EPressureWater]);
	//([kJ/kg][kg/m3] + [kJ/kg][kg/m3])*[m3] - ([kJ/(kg*Kelvin)][Celsius][kg/m3])*(m3) ->...->[kJ]
	T EnergiaT = fCellVolume*fPorosityRock*(EnthalpyOil(Temp)*DensityOil(Temp)*state[ESaturationOil]
											+EnthalpyWater(Temp)*DensityWater(Temp)*(state[ESaturationWater]/*+state[ESaturationSteam]*/)
											+EnthalpySteam(Temp)*DensitySteam(Temp)*state[ESaturationSteam])
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
	//cout << "\nETotalInicial --> " << fInitialState[ETotalEnergy] <<endl;
	//cout << "\nETotalFinal --> " << state[ETotalEnergy] <<endl;
	
	T TotalEnergy = state[ETotalEnergy]; 
	if(state[ETotalEnergy] < EnergiaT/*EnergyCondensedSystem(state)*/){
		// Energy < EMin
		residual[EqZeroSaturation] = state[ESaturationSteam];
	} else {
		residual[EqZeroSaturation] = state[ETemperature] - Temp/*CondensationTemperature(state)*/;
	}

}

template<class T>
inline void TPBrCellMarx::TotalResidual(TPZVec<REAL> &leftval, TPZVec<T> &state, TPZVec<T> &residual)
{
	TPZManVector<T> fluxl(NUMVARS,0.);
	TPZManVector<T> fluxr(NUMVARS,0.);
	TPZManVector<T> internal(NUMVARS,0.);
	Inflow(leftval, fluxl);
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

void TPBrCellMarx::ExtractMatrix(TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > &input, TPZFMatrix &output)
{
	output.Resize(NUMVARS, NUMVARS);
	int i,j;
	for (i=0; i<NUMVARS; i++) {
		for (j=0; j<NUMVARS; j++) {
			output(i,j) = input[i].d(j);
		}
	}
}

void TPBrCellMarx::ExtractMatrix(TPZManVector<REAL> &input, TPZFMatrix &output)
{
	output.Resize(input.NElements(), 1);
	int i;
	for (i=0; i<input.NElements(); i++) {
		output(i,0) = input[i];
	}
}

template<>
inline void TPBrCellMarx::ConvertState(TPZFMatrix &input, TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > &output)
{
	int nrows = input.Rows();
	int i;
	for (i=0; i<nrows; i++) {
		output[i] = input(i,0);
		output[i].fastAccessDx(i) = 1.;
	}
}
template<class T>
void TPBrCellMarx::ConvertState(TPZFMatrix &input, TPZManVector<T> &output)
{
	InitializeState(output);
	int nrows = input.Rows();
	int i;
	for (i=0; i<nrows; i++) {
		output[i] = input(i,0);
	}
}

/**
 template<class T>
 inline void TPBrCellMarx::InternalEquations(TPZVec<REAL> &leftval, TPZVec<T> &residual)
 {
 TPZManVector<T> state;
 InitializeState(state);
 // conservation of mass of oil in the cell -> ESaturationOil
 residual[0] = fCellVolume*fPorosity*fDensity[EOil]*(state[ESaturationOil]-fInitialState[ESaturationOil]);
 // conservation of mass of water in the cell ->ESaturationWater
 residual[1] = fCellVolume*fPorosity*fDensity[EWater]*(state[ESaturationWater]-fInitialState[ESaturationWater])+state[EPhaseChange];
 // conservation of mass of steam in the cell ->ESaturationSteam
 residual[2] = fCellVolume*fPorosity*fDensity[ESteam]*(state[ESaturationSteam]-fInitialState[ESaturationSteam])-state[EPhaseChange];
 
 // relation between massflow rate and darcy velocity ->EMassFluxOil
 residual[3] = state[EMassFluxOil]-state[EDarcyVelocityOil]*fDensity[EOil]*fRightArea;
 // relation between massflow rate and darcy velocity ->EMassFluxWater
 residual[4] = state[EMassFluxWater]-state[EDarcyVelocityWater]*fDensity[EWater]*fRightArea;
 // relation between massflow rate and darcy velocity ->EMassFluxSteam
 residual[5] = state[EMassFluxSteam]-state[EDarcyVelocitySteam]*fDensity[ESteam]*fRightArea;
 
 // relation between dp/dx and the darcy velocity ->EDarcyVelocityOil
 residual[6] = (state[EPressureOil]-leftval[EPressureOil])*fPermeability[EOil]/fViscosity[EOil]-state[EDarcyVelocityOil]*fCellSize;
 // relation between dp/dx and the darcy velocity -> EDarcyVelocityWater
 residual[7] = (state[EPressureWater]-leftval[EPressureWater])*fPermeability[EWater]/fViscosity[EWater]-state[EDarcyVelocityWater]*fCellSize;
 // relation between dp/dx and the darcy velocity -> EDarcyVelocitySteam
 residual[8] = (state[EPressureSteam]-leftval[EPressureSteam])*fPermeability[ESteam]/fViscosity[ESteam]-state[EDarcyVelocitySteam]*fCellSize;
 
 // relation between the pressure of oil and water (capillary pressure) -> EPressureWater
 residual[9] = state[EPressureOil]-state[EPressureWater]; 
 // relation between the pressure of oit and steam (capillary pressure) -> EPressureSteam
 residual[10] = state[EPressureOil]-state[EPressureSteam];
 // the sum of saturations is equal 1 -> EPressureOil
 residual[11] = state[ESaturationOil]+state[ESaturationWater]+state[ESaturationSteam]-1.;
 
 // ENERGIA **************
 
 // energy of the different fases : oil -> EEnthalpyOil
 residual[12] = state[EEnthalpyOil]-fSpecificHeat[EOil]*state[ETemperature]; 
 // energy of the different fases : water -> EEnthalpyWater
 residual[13] = state[EEnthalpyWater]-fSpecificHeat[EWater]*state[ETemperature];
 // energy of the different fases : steam -> EEnthalpySteam
 residual[14] = state[EEnthalpySteam]-fSpecificHeat[ESteam]*state[ETemperature];
 // total energy -> ETotalEnergy
 residual[15] = state[ETotalEnergy]-fPorosity*(state[EEnthalpyOil]*fDensity[EOil]*state[ESaturationOil]
 +state[EEnthalpyWater]*fDensity[EWater]*state[ESaturationWater]+state[EEnthalpySteam]*fDensity[ESteam]*state[ESaturationSteam])
 -(1.-fPorosity)*fSpecificHeatRock*state[ETemperature];
 ;
 // conservation of energy internal contribution -> ETemperature -> EPhaseChange
 residual[16] = (state[ETotalEnergy]-fInitialState[ETotalEnergy])*fCellVolume;
 if(state[ETotalEnergy] < 0./*EnergyCondensedSystem(state)* /)
{
	// Energy < EMin
	residual[17] = state[ESaturationSteam];
} else {
	residual[17] = state[ETemperature] - 100./*CondensationTemperature(state)* /;
}

}
*/
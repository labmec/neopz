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
	REAL fThermalConductivityRock;
	
	REAL fResidualOil;
	
	REAL fMaterialPermeability;//[m2]
	
	// propriedades dos fluidos
	/**
	 *fViscosity [Pa*sec]
	 *fDensity [kg/m3]
	 *fSpecificHeat [kJ/(kg*Kelvin)]
	 *fPressure [Pa]
	 *fDarcyVelocity [m/s]
	 *fMassflux [kg/s]
	 */
	TPZManVector<REAL,3> fViscosity;
	TPZManVector<REAL,3> fDensity;
	TPZManVector<REAL,3> fSpecificHeat;
	TPZManVector<REAL,3> fPressure;
	TPZManVector<REAL,3> fDarcyVelocity;
	TPZManVector<REAL,3> fMassflux;
	TPZManVector<REAL,3> fSaturation;
	
	/**
	 *fEnergy [kJoule] = [1000 kg*(m2/s2)]
	 *fTemperatura [Celsius]
	*/
	TPZManVector<REAL,3> fEnergy;
	REAL fTotalEnergy;
	REAL fTemperature;
	
	REAL fPhaseChange;
	
	TPZManVector<REAL> fInitialState;
public:
	enum VARNAME {EMassFluxOil,EMassFluxWater,EMassFluxSteam,EDarcyVelocityOil,EDarcyVelocityWater,EDarcyVelocitySteam,
		EPressureOil,EPressureWater,EPressureSteam,ESaturationOil,ESaturationWater,ESaturationSteam,
		EEnergyOil,EEnergyWater,EEnergySteam,ETotalEnergy,ETemperature,EPhaseChange};
	enum VARINDEX {EOil, EWater, ESteam};
	
	static const int NUMVARS = 18;
	
	TPBrCellMarx();
	
	virtual ~TPBrCellMarx();
	
	void SetMaterialProperty(REAL &materialpermeability, TPZVec<REAL> &viscosity, TPZVec<REAL> &specificheat, REAL delt)
	{
		fMaterialPermeability = materialpermeability;
		fViscosity = viscosity;
		fSpecificHeat = specificheat;
		fTimeStep = delt;
	}
	
	void SetGeometry(REAL cellvolume, REAL leftarea, REAL rightarea, REAL cellsize, REAL porosityrock, REAL densityrock)
	{
		fCellVolume = cellvolume;
		fLeftArea = leftarea;
		fRightArea = rightarea;
		fCellSize = cellsize;
		fPorosityRock = porosityrock;
		fDensityRock = densityrock;
	}
	
	void SetState(TPZVec<REAL> &pressure, TPZVec<REAL> massflux, TPZVec<REAL> darcyvelocity, TPZVec<REAL> &density,
				  TPZVec<REAL> &saturation, TPZVec<REAL> &enthalpy, REAL temperature )
	{
		fPressure = pressure;
		fMassflux = massflux;
		fDarcyVelocity = darcyvelocity;
		fDensity = density;
		fSaturation = saturation;
		fEnthalpy = enthalpy;
		fTemperature = temperature;
		fTotalEnergy = fPorosityRock*(enthalpy[0]*density[0]*saturation[0] + enthalpy[1]*density[1]*saturation[1] + enthalpy[2]*density[2]*saturation[2])
			+(1-fPorosityRock)*fSpecificHeatRock*fDensityRock*fTemperature;
	}
	
	template<class T>
	//[kJ] = [1000 kg (m2/s2)]
	T EnergySolid(T temperature);
	
	template<class T>
	//[kJ] = [1000 kg (m2/s2)]
	T EnergySaturatedSteam(T pressuresteam);
	
	template<class T>
	//[Celsius]
	T TemperatureSaturatedSteam(T pressuresteam);

	template<class T>
	void InitializeState(TPZManVector<T> &state);
	
	template<class T>
	void Inflow(TPZVec<REAL> &leftval, TPZVec<T> &flux);
	
	template<class T>
	void Outflow(TPZVec<REAL> &rightval, TPZVec<T> &flux);
	
	template<class T>
	void InternalEquations(TPZVec<REAL> &leftval, TPZVec<T> &residual);
	
	template<class T>
	void TotalResidual(TPZVec<REAL> &leftval, TPZVec<REAL> &rightval, TPZVec<T> &residual);
	
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
}
TPBrCellMarx::~TPBrCellMarx() {
}

template<class T>
inline void TPBrCellMarx::Inflow(TPZVec<REAL> &leftval, TPZVec<T> &flux)
{
	flux[0] = leftval[EMassFluxOil]*fTimeStep;
	flux[1] = leftval[EMassFluxWater]*fTimeStep;
	flux[2] = leftval[EMassFluxSteam]*fTimeStep;
	//
	flux[16] = (leftval[EEnergyOil]*leftval[EDarcyVelocityOil]*fDensity[EOil]
			   +leftval[EEnergyWater]*leftval[EDarcyVelocityWater]*fDensity[EWater]
			   +leftval[EEnergySteam]*leftval[EDarcyVelocitySteam]*fDensity[ESteam]
			   )*fLeftArea*fTimeStep;
}

template<class T>
inline void TPBrCellMarx::Outflow(TPZVec<REAL> &rightval, TPZVec<T> &flux)
{
	TPZManVector<T> state;
	InitializeState(state);
	flux[0] = state[EMassFluxOil]*fTimeStep;
	flux[1] = state[EMassFluxWater]*fTimeStep;
	flux[2] = state[EMassFluxSteam]*fTimeStep;
	flux[16] = (state[EEnergyOil]*state[EDarcyVelocityOil]*fDensity[EOil]
				+state[EEnergyWater]*state[EDarcyVelocityWater]*fDensity[EWater]
				+state[EEnergySteam]*state[EDarcyVelocitySteam]*fDensity[ESteam]
				)*fRightArea*fTimeStep;
}

template<class T>
inline void TPBrCellMarx::InternalEquations(TPZVec<REAL> &leftval, TPZVec<T> &residual)
{
	TPZManVector<T> state;
	InitializeState(state);
	
	// conservation of mass of oil in the cell -> ESaturationOil
	// [m3]*[Kg/m3] -> [kg]
	residual[0] = fCellVolume*fPorosityRock*fDensity[EOil]*(state[ESaturationOil]-fInitialState[ESaturationOil]);
	
	// conservation of mass of water in the cell ->ESaturationWater
	// [m3]*[Kg/m3] -> [kg]
	residual[1] = fCellVolume*fPorosityRock*fDensity[EWater]*(state[ESaturationWater]-fInitialState[ESaturationWater])+state[EPhaseChange];
	
	// conservation of mass of steam in the cell ->ESaturationSteam
	// [m3]*[Kg/m3] -> [kg]
	residual[2] = fCellVolume*fPorosityRock*fDensity[ESteam]*(state[ESaturationSteam]-fInitialState[ESaturationSteam])-state[EPhaseChange];
	
	// relation between massflow rate and darcy velocity ->EMassFluxOil
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [gk/s]
	residual[3] = state[EMassFluxOil]-state[EDarcyVelocityOil]*DensityOil(fInitialState[ETemperature])*fRightArea;
	
	// relation between massflow rate and darcy velocity ->EMassFluxWater
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [gk/s]
	residual[4] = state[EMassFluxWater]-state[EDarcyVelocityWater]*DensityWater(fInitialState[ETemperature])*fRightArea;
	
	// relation between massflow rate and darcy velocity ->EMassFluxSteam
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [gk/s]
	residual[5] = state[EMassFluxSteam]-state[EDarcyVelocitySteam]*DensitySteam(fInitialState[ETemperature])*fRightArea;
	
	TPZManVector<REAL> saturation(3), relatpermeability(3);
	saturation[EOil]=fInitialState[ESaturationOil];
	saturation[EWater]=fInitialState[ESaturationWater ];
	saturation[ESteam]=fInitialState[ESaturationSteam ];
	ComputeRelativePermeability(saturation, relatpermeability);
	// relation between dp/dx and the darcy velocity ->EDarcyVelocityOil
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s]  
	residual[6] = (state[EPressureOil]-leftval[EPressureOil])*fMaterialPermeability*relatpermeability[EOil]/ViscosityOil(fInitialState[ETemperature])-state[EDarcyVelocityOil]*fCellSize;
	
	// relation between dp/dx and the darcy velocity -> EDarcyVelocityWater
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s]
	residual[7] = (state[EPressureWater]-leftval[EPressureWater])*fMaterialPermeability*relatpermeability[EWater]/ViscosityWater(fInitialState[ETemperature])-state[EDarcyVelocityWater]*fCellSize;
	
	// relation between dp/dx and the darcy velocity -> EDarcyVelocitySteam
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s]
	residual[8] = (state[EPressureSteam]-leftval[EPressureSteam])*fMaterialPermeability*relatpermeability[ESteam]/ViscositySteam(fInitialState[ETemperature])-state[EDarcyVelocitySteam]*fCellSize;

	// relation between the pressure of oil and water (capillary pressure) -> EPressureWater
	//[Pa] - [Pa] -> [Pa]
	if(state[ESaturationOil]<fResidualOil) {
		residual[9] = state[ESaturationOil]-fResidualOil;
	}
	else {
		residual[9] = state[EPressureOil]-state[EPressureWater]; 
	}
	
	// relation between the pressure of oit and steam (capillary pressure) -> EPressureSteam
	//[Pa] - [Pa] -> [Pa]
	residual[10] = state[EPressureWater]-state[EPressureSteam];
	
	// the sum of saturations is equal 1 -> EPressureOil
	//dimensionless
	residual[11] = state[ESaturationOil]+state[ESaturationWater]+state[ESaturationSteam] - 1.;

	// ENERGIA **************
	T Temp = TemperatureSaturation(fPressure[EWater]);
	//([kJ/kg][kg/m3] + [kJ/kg][kg/m3]) - ([kJ/(kg*Kelvin)][Celsius][kg/m3]) ->...->[kJ/m3]
	T EnergiaT = fPorosityRock*(EnthalpyOil(Temp)*DensityOil(Temp)*state[ESaturationOil]
					+EnthalpyWater(Temp)*DensityWater(Temp)*(state[ESaturationWater]+state[ESaturationSteam]))
											-(1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*state[ETemperature];
	
	// energy of the different fases : oil -> EEnergyOil
	residual[12] = state[EEnergyOil]-fSpecificHeat[EOil]*state[ETemperature]; 
	// energy of the different fases : water -> EEnergyWater
	residual[13] = state[EEnergyWater]-fSpecificHeat[EWater]*state[ETemperature];
	// energy of the different fases : steam -> EEnergySteam
	residual[14] = state[EEnergySteam]-fSpecificHeat[ESteam]*state[ETemperature];
	
	// total energy -> ETotalEnergy
	//ETotalEnergy -> ([kJ/kg]*[kg/m3])*[m3] - ([kJ/(kg*Kelvin)]*[kg/m3]*[m3]*[Celsius]) -> [kJ] 
	residual[15] = state[ETotalEnergy]-fPorosityRock*(EnthalpyOil(state[ETemperature])*DensityOil(state[ETemperature])*state[ESaturationOil]
		+EnthalpyWater(state[ETemperature])*DensityWater(state[ETemperature])*state[ESaturationWater]
		+EnthalpySteam(state[ETemperature])*DensitySteam(state[ETemperature])*state[ESaturationSteam])*fCellVolume
	-(1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*fCellVolume*state[ETemperature];
	
	// conservation of energy internal contribution -> ETemperature -> EPhaseChange
	//[kJ]
	residual[16] = (state[ETotalEnergy]-fInitialState[ETotalEnergy]);
	
	if(state[ETotalEnergy] < EnergiaT/*EnergyCondensedSystem(state)*/)
	{
		// Energy < EMin
		residual[17] = state[ESaturationSteam];
	} else {
		residual[17] = state[ETemperature] - Temp/*CondensationTemperature(state)*/;
	}

}

template<class T>
inline void TPBrCellMarx::TotalResidual(TPZVec<REAL> &leftval, TPZVec<REAL> &rightval, TPZVec<T> &residual)
{
	TPZManVector<T> fluxl(NUMVARS,0.);
	TPZManVector<T> fluxr(NUMVARS,0.);
	TPZManVector<T> internal(NUMVARS,0.);
	Inflow(leftval, fluxl);
	Outflow(rightval, fluxr);
	InternalEquations(leftval, internal);
	int i;
	for (i=0; i<NUMVARS; i++) {
		residual[i] = fluxl[i]-fluxr[i]+internal[i];
	}
}

template<>
inline void TPBrCellMarx::InitializeState(TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > &state)
{
	state.Resize(NUMVARS);
	state[EPressureOil] = fPressure[EOil];
	state[EPressureOil].fastAccessDx(EPressureOil) = 1.;
	state[EPressureWater] = fPressure[EWater];
	state[EPressureWater].fastAccessDx(EPressureWater) = 1.;
	state[EPressureSteam] = fPressure[ESteam];
	state[EPressureSteam].fastAccessDx(EPressureSteam) = 1.;
	state[EMassFluxOil] = fMassflux[EOil];
	state[EMassFluxOil].fastAccessDx(EMassFluxOil) = 1.;
	state[EMassFluxWater] = fMassflux[EWater];
	state[EMassFluxWater].fastAccessDx(EMassFluxWater) = 1.;
	state[EMassFluxSteam] = fMassflux[ESteam];
	state[EMassFluxSteam].fastAccessDx(EMassFluxSteam) = 1.;
	state[EDarcyVelocityOil] = fDarcyVelocity[EOil];
	state[EDarcyVelocityOil].fastAccessDx(EDarcyVelocityOil) = 1.;
	state[EDarcyVelocityWater] = fDarcyVelocity[EWater];
	state[EDarcyVelocityWater].fastAccessDx(EDarcyVelocityWater) = 1.;
	state[EDarcyVelocitySteam] = fDarcyVelocity[ESteam];
	state[EDarcyVelocitySteam].fastAccessDx(EDarcyVelocitySteam) = 1.;
	state[ESaturationOil] = fSaturation[EOil];
	state[ESaturationOil].fastAccessDx(ESaturationOil) = 1.;
	state[ESaturationWater] = fSaturation[EWater];
	state[ESaturationWater].fastAccessDx(ESaturationWater) = 1.;
	state[ESaturationSteam] = fSaturation[ESteam];
	state[ESaturationSteam].fastAccessDx(ESaturationSteam) = 1.;
	state[EEnergyOil] = fEnergy[EOil];
	state[EEnergyOil].fastAccessDx(EEnergyOil) = 1.;
	state[EEnergyWater] = fEnergy[EWater];
	state[EEnergyWater].fastAccessDx(EEnergyWater) = 1.;
	state[EEnergySteam] = fEnergy[ESteam];
	state[EEnergySteam].fastAccessDx(EEnergySteam) = 1.;
	state[ETotalEnergy] = fTotalEnergy;
	state[ETotalEnergy].fastAccessDx(ETotalEnergy) = 1.;
	state[ETemperature] = fTemperature;
	state[ETemperature].fastAccessDx(ETemperature) = 1.;
	state[EPhaseChange] = fPhaseChange;
	state[EPhaseChange].fastAccessDx(EPhaseChange) = 1.;

}

template<>
inline void TPBrCellMarx::InitializeState(TPZManVector<REAL> &state)
{
	state.Resize(NUMVARS);
	state[EPressureOil] = fPressure[EOil];
	state[EPressureWater] = fPressure[EWater];
	state[EPressureSteam] = fPressure[ESteam];
	state[EMassFluxOil] = fMassflux[EOil];
	state[EMassFluxWater] = fMassflux[EWater];
	state[EMassFluxSteam] = fMassflux[ESteam];
	state[EDarcyVelocityOil] = fDarcyVelocity[EOil];
	state[EDarcyVelocityWater] = fDarcyVelocity[EWater];
	state[EDarcyVelocitySteam] = fDarcyVelocity[ESteam];
	state[ESaturationOil] = fSaturation[EOil];
	state[ESaturationWater] = fSaturation[EWater];
	state[ESaturationSteam] = fSaturation[ESteam];
	state[EEnergyOil] = fEnergy[EOil];
	state[EEnergyWater] = fEnergy[EWater];
	state[EEnergySteam] = fEnergy[ESteam];
	state[ETotalEnergy] = fTotalEnergy;
	state[ETemperature] = fTemperature;
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
 
 // energy of the different fases : oil -> EEnergyOil
 residual[12] = state[EEnergyOil]-fSpecificHeat[EOil]*state[ETemperature]; 
 // energy of the different fases : water -> EEnergyWater
 residual[13] = state[EEnergyWater]-fSpecificHeat[EWater]*state[ETemperature];
 // energy of the different fases : steam -> EEnergySteam
 residual[14] = state[EEnergySteam]-fSpecificHeat[ESteam]*state[ETemperature];
 // total energy -> ETotalEnergy
 residual[15] = state[ETotalEnergy]-fPorosity*(state[EEnergyOil]*fDensity[EOil]*state[ESaturationOil]
 +state[EEnergyWater]*fDensity[EWater]*state[ESaturationWater]+state[EEnergySteam]*fDensity[ESteam]*state[ESaturationSteam])
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
#ifndef STEAMFLUXHPP
#define STEAMFLUXHPP

#include "pzvec.h"

#ifdef _AUTODIFF
using namespace std;
#include "fadType.h"

class TPBrSteamSimulation;
class TPBrSteamMesh;
#include "ThermalMethodsTables.h"

/*
 *  tpbrsteamflux.h
 *  PZ
 *
 *  Created by Philippe Devloo on 3/8/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */

//extern WaterDataInStateOfSaturation waterdata;
//extern OilData oildata;

class TPBrSteamFlux
{
    friend class TPBrSteamSimulation;
    friend class TPBrSteamMesh;
    
	REAL fMaterialPermeability;
	
	static REAL fFarfieldPressureOil;
	static REAL fFarfieldPressureWater;
	static REAL fFarfieldPressureSteam;
	static REAL fFarfieldTemperature;
	static REAL fFarfieldSaturationOil;
	static REAL fFarfieldSaturationWater;
	static REAL fFarfieldSaturationSteam;
	
	static REAL fInletEnergyFlux; //[KJ/s]
	static REAL fInletMassFlux; //[Kg/s];
	
	
public:
	
	/// empty constructor initializing all variables
	TPBrSteamFlux();
	
	enum ESteamFluxEq {
		EMassFluxWaterEq, EMassFluxSteamEq, EMassFluxOilEq, EDarcyVelocityWaterEq, EDarcyVelocitySteamEq, EDarcyVelocityOilEq, EEnergyFluxEq
	};
	enum ESteamFluxVars {
		EMassFluxWater, EMassFluxSteam, EMassFluxOil, EDarcyVelocityWater, EDarcyVelocitySteam, EDarcyVelocityOil, EEnergyFlux
	};

	enum EInletEq { EInletEnergyFlux, EInletStateEqs, EInletMassFlux };
	
	enum EInletVars { EInletPressure, EInletSteamSaturation, EInletTemperature };
	
	
	enum VARINDEX {EOil, EWater, ESteam};

	
	static const int NumFluxEq = 7;
	
	/// this variable allows to incorporate inlet variable specifications such as energy flux, and mass flux
	static const int NumInletVars = 3;

	/// Compute the relative permeabilities between fases
	template<class T>
	void ComputeRelativePermeability(TPZManVector<T> &saturation,TPZManVector<T> &relativepermeability);//dimensionless

	template<class T>
	static T EnthalpyWater(T temperature);//[kJ/kg]
	template<class T>
	static T EnthalpySteam(T temperature);//[kJ/kg]
	template<class T>
	static T EnthalpyOil(T temperature);//[kJ/kg]
	
	template<class T>
	static T ViscosityOil(T temp);//[Pa*sec]
	template<class T>
	static T ViscosityWater(T temp);//[Pa*sec]
	template<class T>
	static T ViscositySteam(T temp);//[Pa*sec]
	
	// metodos para recuperar os dados tabulados em funcao
	template<class T>
	static T DensityOil(T temp);//[kg/m3]
	template<class T>
	static T DensityWater(T temp);//[kg/m3]
	template<class T>
	static T DensitySteam(T temp);//[kg/m3]
	
    template<class T>
	static T TemperatureSaturation(T p);//[Celsius]

	
	/// calcula o fluxo entre duas celulas
	template<class T>
	void FluxResidual(TPZVec<T> &leftstate, TPZVec<T> &rightstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual);
	
	/// calcula a contribuicao para a matriz de rigidez
	void CalcStiff(TPZVec<REAL> &leftstate, TPZVec<REAL> &rightstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
				   TPZFMatrix &ek, TPZFMatrix &ef);
	
	/// calcula a contribuicao para a matriz de rigidez das equacoes de entrada
	void InletCalcStiff(TPZVec<REAL> &rightstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
				   TPZFMatrix &ek, TPZFMatrix &ef);
	
	/// calcula a contribuicao para a matriz de rigidez das equacoes de entrada
	void OutletCalcStiff(TPZVec<REAL> &leftstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
						TPZFMatrix &ek, TPZFMatrix &ef);
	
    /// Compute a limit for correcting the solution
    static REAL LimitRangeInlet(REAL scale,TPZVec<REAL> &inletstate,TPZVec<REAL> &cellstate, TPZVec<REAL> &inletcorrection, TPZVec<REAL> &cellcorrection);
	
    /// Compute a limit for correcting the solution
    static REAL LimitRange(REAL scale,TPZVec<REAL> &interfacestate,TPZVec<REAL> &interfacecorrection);
	
	/// Incorporate the partial derivatives in the state variables
	template<int N>
	static void Initialize(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset);
	
	/// Incorporate the partial derivatives in the state variables
	template<int N>
	static void InitializeInlet(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset);
	
	/// complete residual vector as a function of the inletstate and rightstate
	template<class T>
	void InletFluxResidual(TPZVec<T> &rightstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual);

	/// complete residual vector as a function of the inletstate and rightstate
	template<class T>
	void OutletFluxResidual(TPZVec<T> &leftstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual);
	
	/// compute left state as a function of the inletstate
	template<class T>
	void ComputeLeftState(TPZVec<T> &inletstate, TPZVec<T> &leftstate);
    
    void Print(std::ostream &out = std::cout);
    
    /// associate scale factors with the equations and state variables
    static void Scales(TPZVec<REAL> &eqscales, TPZVec<REAL> &statescales);

    /// associate scale factors with the equations and state variables for the inlet equations
    static void InletScales(TPZVec<REAL> &eqscales, TPZVec<REAL> &statescales);

};


/// Incorporate the partial derivatives in the state variables
template<int N>
inline void TPBrSteamFlux::Initialize(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset)
{
	// EMassFluxWater, EMassFluxSteam, EMassFluxOil, EDarcyVelocityWater, EDarcyVelocitySteam, EDarcyVelocityOil, EEnergyFlux
	
	fadstate[EMassFluxWater].val() = state[EMassFluxWater];
	fadstate[EMassFluxWater].fastAccessDx(EMassFluxWater+offset) = 1.;
	fadstate[EMassFluxSteam] = state[EMassFluxSteam];
	fadstate[EMassFluxSteam].fastAccessDx(EMassFluxSteam+offset) = 1.;
	fadstate[EMassFluxOil].val() = state[EMassFluxOil];
	fadstate[EMassFluxOil].fastAccessDx(EMassFluxOil+offset) = 1.;
	fadstate[EDarcyVelocityWater].val() = state[EDarcyVelocityWater];
	fadstate[EDarcyVelocityWater].fastAccessDx(EDarcyVelocityWater+offset) = 1.;
	fadstate[EDarcyVelocitySteam].val() = state[EDarcyVelocitySteam];
	fadstate[EDarcyVelocitySteam].fastAccessDx(EDarcyVelocitySteam+offset) = 1.;
	fadstate[EDarcyVelocityOil].val() = state[EDarcyVelocityOil];
	fadstate[EDarcyVelocityOil].fastAccessDx(EDarcyVelocityOil+offset) = 1.;
	fadstate[EEnergyFlux].val() = state[EEnergyFlux];
	fadstate[EEnergyFlux].fastAccessDx(EEnergyFlux+offset) = 1.;
	
}


// nothing is compiled if _AUTODIFF isnt defined
#endif

#endif
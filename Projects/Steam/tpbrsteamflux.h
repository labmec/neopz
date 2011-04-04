#ifndef STEAMFLUXHPP
#define STEAMFLUXHPP

#include "pzvec.h"

#ifdef _AUTODIFF
using namespace std;
#include "fadType.h"
/*
 *  tpbrsteamflux.h
 *  PZ
 *
 *  Created by Philippe Devloo on 3/8/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */

class TPBrSteamFlux
{
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

	enum EInletEq { EInletMassFlux, EInletEnergyFlux };
	
	enum EInletVars { EInletPressure, EInletSteamSaturation };
	
	
	enum VARINDEX {EOil, EWater, ESteam};


	/// Compute the relative permeabilities between fases
	template<class T>
	void ComputeRelativePermeability(TPZManVector<T> &saturation,TPZManVector<T> &relativepermeability);//dimensionless

	template<class T>
	T EnthalpyWater(T temperature);//[kJ/kg]
	template<class T>
	T EnthalpySteam(T temperature);//[kJ/kg]
	template<class T>
	T EnthalpyOil(T temperature);//[kJ/kg]
	
	template<class T>
	T ViscosityOil(T temp);//[Pa*sec]
	template<class T>
	T ViscosityWater(T temp);//[Pa*sec]
	template<class T>
	T ViscositySteam(T temp);//[Pa*sec]
	
	// metodos para recuperar os dados tabulados em funcao
	template<class T>
	T DensityOil(T temp);//[kg/m3]
	template<class T>
	T DensityWater(T temp);//[kg/m3]
	template<class T>
	T DensitySteam(T temp);//[kg/m3]
	
	
	/// calcula o fluxo entre duas celulas
	template<class T>
	void FluxResidual(TPZVec<T> &leftstate, TPZVec<T> &rightstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual);
	
	/// calcula a contribuicao para a matriz de rigidez
	void CalcStiff(TPZVec<REAL> &leftstate, TPZVec<REAL> &rightstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
				   TPZFMatrix &ek, TPZFMatrix &ef);
	
	/// calcula a contribuicao para a matriz de rigidez das equacoes de entrada
	void InletCalcStiff(TPZVec<REAL> &inletstate, TPZVec<REAL> &rightstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
				   TPZFMatrix &ek, TPZFMatrix &ef);
	
	/// calcula a contribuicao para a matriz de rigidez das equacoes de entrada
	void OutletCalcStiff(TPZVec<REAL> &leftstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
						TPZFMatrix &ek, TPZFMatrix &ef);
	
	
	static const int NumFluxEq = 7;
	
	/// this variable allows to incorporate inlet variable specifications such as energy flux, and mass flux
	static const int NumInletVars = 2;
	
	/// Incorporate the partial derivatives in the state variables
	template<int N>
	static void Initialize(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset);
	
	/// Incorporate the partial derivatives in the state variables
	template<int N>
	static void InitializeInlet(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset);
	
	/// complete residual vector as a function of the inletstate and rightstate
	template<class T>
	void InletFluxResidual(TPZVec<T> &inletstate, TPZVec<T> &rightstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual);

	/// complete residual vector as a function of the inletstate and rightstate
	template<class T>
	void OutletFluxResidual(TPZVec<T> &leftstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual);
	
	/// compute left state as a function of the inletstate
	template<class T>
	void ComputeLeftState(TPZVec<T> &inletstate, TPZVec<T> &leftstate);
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
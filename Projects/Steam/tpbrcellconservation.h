#ifndef TPBRCELLCONSERVATION
#define TPBRCELLCONSERVATION
/*
 *  tpbrcellconservation.h
 *  PZ
 *
 *  Created by Philippe Devloo on 3/8/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */

#include "pzvec.h"
#include "pzfmatrix.h"
#include "ThermalMethodsTables.h"

#ifdef _AUTODIFF

using namespace std;
#include "fadType.h"

extern WaterDataInStateOfSaturation waterdata;
extern OilData oildata;

class TPBrCellConservation
{
	static REAL fPorosityRock;
	static REAL fDensityRock;
	static REAL fSpecificHeatRock;
	static REAL fResidualOil;
	
public:
	
	// metodos para recuperar os dados tabulados em funcao
	template<class T>
	T DensityOil(T temp);//[kg/m3]
	template<class T>
	T DensityWater(T temp);//[kg/m3]
	template<class T>
	T DensitySteam(T temp);//[kg/m3]

	template<class T>
	static T TemperatureSaturation(T p);//[Celsius]
	template<class T>
	T EnthalpyWater(T temperature);//[kJ/kg]
	template<class T>
	T EnthalpySteam(T temperature);//[kJ/kg]
	template<class T>
	T EnthalpyOil(T temperature);//[kJ/kg]
	
	
	
public:
	/// Equation associated with a cell
	enum ECellEq { EMassWater, EMassOil, EMassSteam, ECapillaryPressureWO, ECapillaryPressureOS, ESumSaturation, EEnergyCons,
		EZeroSaturation
	};
	
	/// State variables associated with a cell
	enum ECellState {
		ESaturationWater, ESaturationOil, ESaturationSteam, ETemperature, EPressureWater, EPressureSteam, EPressureOil,EPhaseChange
	};
	
	static const int NumCellEq = 8;
	
	/// calcula o residuo das equacoes associadas a celula
	template<class T>
	void CellResidual(TPZVec<T> &leftflux, TPZVec<T> &cellstate, TPZVec<T> &rightflux, TPZVec<REAL> &initialstate, 
					  REAL volume, REAL delt, TPZVec<T> &cellresidual);
	
	/// calcula a contribuicao para o residuo e matriz tangente
	void CalcStiff(TPZVec<REAL> &leftflux, TPZVec<REAL> &cellstate, TPZVec<REAL> &rightflux, TPZVec<REAL> &initialstate, 
				   REAL volume, REAL delt, TPZFMatrix &ek, TPZFMatrix &ef);
	
	/// Incorporate the partial derivatives in the state variables
	template<int N>
	static void Initialize(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset);
};


// metodos para recuperar os dados tabulados em funcao
template<class T>
inline T TPBrCellConservation::DensityOil(T temp)//[kg/m3]
{
	return oildata.getDensityToOil(temp);
}
template<class T>
inline T TPBrCellConservation::DensityWater(T temp)//[kg/m3]
{
	return waterdata.getSaturationStateDensityToLiquidWater(temp);
}
template<class T>
inline T TPBrCellConservation::DensitySteam(T temp)//[kg/m3]
{
	return waterdata.getSaturationStateDensityToSteam(temp);
}

template<class T>
inline T TPBrCellConservation::TemperatureSaturation(T pressuresteam)//[Celsius]
{
	T  val_log, temp;
	T temp_de_saturac;
	val_log = log(pressuresteam*0.0001450377438972831);
	temp=561.435 + 33.8866*val_log + 2.18893*(val_log*val_log) + 0.0808998*(val_log*val_log*val_log) +
	0.0342030*(val_log*val_log*val_log*val_log);
	temp_de_saturac = (temp - 32. - 459.67)/1.8;
	
	return temp_de_saturac;
}
template<class T>
inline T TPBrCellConservation::EnthalpyWater(T temperature)//[kJ/kg]
{
	return waterdata.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
}
template<class T>
inline T TPBrCellConservation::EnthalpySteam(T temperature)//[kJ/kg]
{
	return waterdata.getSaturationStateSpecificEnthalpyToSteam(temperature);
}
template<class T>
inline T TPBrCellConservation::EnthalpyOil(T temperature)//[kJ/kg]
{
	return  oildata.getSpecificHeatToOil(temperature);	
}

/// Incorporate the partial derivatives in the state variables
template<int N>
inline void TPBrCellConservation::Initialize(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset)
{
	//ESaturationWater, ESaturationOil, ESaturationSteam, ETemperature, EPressureWater, EPressureSteam, EPressureOil,EPhaseChange
	fadstate[ESaturationWater] = state[ESaturationWater];
	fadstate[ESaturationWater].fastAccessDx(ESaturationWater+offset) = 1.;
	fadstate[ESaturationOil] = state[ESaturationOil];
	fadstate[ESaturationOil].fastAccessDx(ESaturationOil+offset) = 1.;
	fadstate[ESaturationSteam] = state[ESaturationSteam];
	fadstate[ESaturationSteam].fastAccessDx(ESaturationSteam+offset) = 1.;
	fadstate[ETemperature] = state[ETemperature];
	fadstate[ETemperature].fastAccessDx(ETemperature+offset) = 1.;
	fadstate[EPressureOil] = state[EPressureOil];
	fadstate[EPressureOil].fastAccessDx(EPressureOil+offset) = 1.;
	fadstate[EPressureWater] = state[EPressureWater];
	fadstate[EPressureWater].fastAccessDx(EPressureWater+offset) = 1.;
	fadstate[EPressureSteam] = state[EPressureSteam];
	fadstate[EPressureSteam].fastAccessDx(EPressureSteam+offset) = 1.;
	fadstate[EPhaseChange] = state[EPhaseChange];
	fadstate[EPhaseChange].fastAccessDx(EPhaseChange+offset) = 1.;
	
}

// Nothing is compiled if _AUTODIFF isnt defined
#endif

#endif
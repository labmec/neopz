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


class TPBrCellConservation
{
	static REAL fPorosityRock;
	static REAL fDensityRock;
	static REAL fSpecificHeatRock;
	static REAL fResidualOil;
	
public:
	
	// metodos para recuperar os dados tabulados em funcao
	template<class T>
	static T DensityOil(T temp);//[kg/m3]
	template<class T>
	static T DensityWater(T temp);//[kg/m3]
	template<class T>
	static T DensitySteam(T temp);//[kg/m3]

	template<class T>
	static T TemperatureSaturation(T pressure);//[Celsius]
	template<class T>
	static T EnthalpyWater(T temperature);//[kJ/kg]
	template<class T>
	static T EnthalpySteam(T temperature);//[kJ/kg]
	template<class T>
	static T EnthalpyOil(T temperature);//[kJ/kg]
	
	
	
public:
	/// Equation associated with a cell
	enum ECellEq { EMassWater, EMassOil, EEnergyCons, ECapillaryPressureOS, ECapillaryPressureWO, EMassSteam,
		EZeroSaturation, ESumSaturation
	};
	
	/// State variables associated with a cell
	enum ECellState {
		ESaturationWater, ESaturationOil, ETemperature, EPressureSteam, EPressureOil,EPhaseChange, ESaturationSteam, EPressureWater
	};
	
	static const int NumCellEq = 8;
	
	/// calcula o residuo das equacoes associadas a celula
	template<class T>
	void CellResidual(TPZVec<T> &leftflux, TPZVec<T> &cellstate, TPZVec<T> &rightflux, TPZVec<REAL> &initialstate, 
					  REAL volume, REAL delt, TPZVec<T> &cellresidual);
	
	/// calcula a contribuicao para o residuo e matriz tangente
	void CalcStiff(TPZVec<REAL> &leftflux, TPZVec<REAL> &cellstate, TPZVec<REAL> &rightflux, TPZVec<REAL> &initialstate, 
				   REAL volume, REAL delt, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
	/// Incorporate the partial derivatives in the state variables
	template<int N>
	static void Initialize(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset);
    
    /// limit the correction of the Newton Iteration
    static REAL LimitRange(REAL scale,TPZVec<REAL> &cellstate,TPZVec<REAL> &cellcorrection);
    
    /// Print the data of the cell
    void Print(std::ostream &out = std::cout);
    
    /// Compute the total energy of the cell
    static REAL Energy(TPZVec<REAL> &cellstate, REAL volume);
    
    /// associate scale factors with the equations and state variables
    static void Scales(TPZVec<REAL> &eqscales, TPZVec<REAL> &statescales);
};


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

extern WaterDataInStateOfSaturation waterdata;

template<class T>
inline T TPBrCellConservation::TemperatureSaturation(T pressuresteam)//[Celsius]
{
    T press1000 = pressuresteam/T(1000.0);
    T temp_de_saturac = waterdata.getSaturationStateTemperature(press1000);
	return temp_de_saturac;
}


// Nothing is compiled if _AUTODIFF isnt defined
#endif

#endif
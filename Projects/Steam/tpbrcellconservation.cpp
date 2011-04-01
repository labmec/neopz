/*
 *  tpbrcellconservation.cpp
 *  PZ
 *
 *  Created by Philippe Devloo on 3/8/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */

#include "tpbrcellconservation.h"
#include "tpbrsteamflux.h"
#include "fadType.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("br.steam.cell"));
#endif


REAL TPBrCellConservation::fPorosityRock = 0.196;
REAL TPBrCellConservation::fDensityRock = 1.055;
REAL TPBrCellConservation::fSpecificHeatRock = 0.766;
REAL TPBrCellConservation::fResidualOil = 0.05;


/// calcula a contribuicao para o residuo e matriz tangente

void TPBrCellConservation::CalcStiff(TPZVec<REAL> &leftflux, TPZVec<REAL> &cellstate, TPZVec<REAL> &rightflux, TPZVec<REAL> &initialstate,
									 REAL volume, REAL delt, TPZFMatrix &ek, TPZFMatrix &ef)
{
	const int totaleq = 2*TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq;
	TPZManVector<TFad<totaleq,REAL> , TPBrCellConservation::NumCellEq> cellfad(TPBrCellConservation::NumCellEq); 
	TPZManVector<TFad<totaleq,REAL> , TPBrSteamFlux::NumFluxEq> leftfad(TPBrSteamFlux::NumFluxEq), rightfad(TPBrSteamFlux::NumFluxEq);
	
	TPBrSteamFlux::Initialize<totaleq>(leftflux,leftfad,0);
	TPBrCellConservation::Initialize<totaleq>(cellstate,cellfad,TPBrSteamFlux::NumFluxEq);
	TPBrSteamFlux::Initialize<totaleq>(rightflux,rightfad,TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq);

	TPZManVector<TFad<totaleq,REAL> , totaleq> cellresidualfad(NumCellEq); 

	CellResidual(leftfad, cellfad, rightfad, initialstate, volume, delt, cellresidualfad );
	
	ek.Redim(NumCellEq, totaleq);
	ef.Redim(NumCellEq, 1);
	int i,j;
	for (i=0; i<NumCellEq; i++) 
	{
		ef(i,0) = cellresidualfad[i].val();
		for (j=0; j<totaleq; j++) {
			ek(i,j) = cellresidualfad[i].d(j);
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		ek.Print("Cell stiffness",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

/*
 /// Equation associated with a cell
 enum ECellEq { EMassWater, EMassOil, EMassSteam, EEnergy, ECapillaryPressureWO, ECapillaryPressureOS, ESumSaturation, EEnergyCons,
 EZeroSaturation
 };
 
 /// State variables associated with a cell
 enum ECellState {
 ESaturationWater, ESaturationOil, ESaturationSteam, ETemperature, EPressureWater, EPressureSteam, EPressureOil,EPhaseChange
 };
*/ 
/*
 enum ESteamFluxVars {
 EMassFluxWater, EMassFluxSteam, EMassFluxOil, EDarcyVelocityWater, EDarcyVelocitySteam, EDarcyVelocityOil, EEnergyFlux
 };
*/ 

/// calcula o residuo das equacoes associadas a celula
template<class T>
void TPBrCellConservation::CellResidual(TPZVec<T> &leftflux, TPZVec<T> &cellstate, TPZVec<T> &rightflux, TPZVec<REAL> &initialstate, 
				  REAL volume, REAL delt, TPZVec<T> &cellresidual)
{
	// conservation of mass of oil in the cell -> ESaturationOil
	// [m3]*[Kg/m3] -> [kg]
	cellresidual[EMassOil] = 
		volume*fPorosityRock*(DensityOil(cellstate[ETemperature])*cellstate[ESaturationOil]-DensityOil(initialstate[ETemperature])*initialstate[ESaturationOil]);
	cellresidual[EMassOil] += leftflux[TPBrSteamFlux::EMassFluxOil]-rightflux[TPBrSteamFlux::EMassFluxOil];
	// conservation of mass of water in the cell ->ESaturationWater
	// [m3]*[Kg/m3] -> [kg]
	cellresidual[EMassWater] = volume*fPorosityRock*(DensityWater(cellstate[ETemperature])*cellstate[ESaturationWater]-DensityWater(initialstate[ETemperature])*initialstate[ESaturationWater])+cellstate[EPhaseChange];
	cellresidual[EMassWater] += leftflux[TPBrSteamFlux::EMassFluxWater]-rightflux[TPBrSteamFlux::EMassFluxWater];

	// conservation of mass of steam in the cell ->ESaturationSteam
	// [m3]*[Kg/m3] -> [kg]
	cellresidual[EMassSteam] = volume*fPorosityRock*(DensitySteam(cellstate[ETemperature])*cellstate[ESaturationSteam]-DensitySteam(initialstate[ETemperature])*initialstate[ESaturationSteam])-cellstate[EPhaseChange];
	cellresidual[EMassSteam] += leftflux[TPBrSteamFlux::EMassFluxSteam]-rightflux[TPBrSteamFlux::EMassFluxSteam];
	
	
	// relation between the pressure of oil and water (capillary pressure) -> EPressureWater
	//[Pa] - [Pa] -> [Pa]
	if(cellstate[ESaturationOil]<fResidualOil) {
		cellresidual[ECapillaryPressureWO] = cellstate[ESaturationOil]-fResidualOil;
	}
	else {
		cellresidual[ECapillaryPressureWO] = cellstate[EPressureOil]-cellstate[EPressureWater]; 
	}
	
	// relation between the pressure of oil and steam (capillary pressure) -> EPressureSteam
	//[Pa] - [Pa] -> [Pa]
	cellresidual[ECapillaryPressureOS] = cellstate[EPressureWater]-cellstate[EPressureSteam];
	
	// the sum of saturations is equal 1 -> EPressureOil
	//dimensionless
	cellresidual[ESumSaturation] = cellstate[ESaturationOil]+cellstate[ESaturationWater]+cellstate[ESaturationSteam] - 1.;
	
	// ENERGIA **************
	T Temp = TemperatureSaturation(cellstate[EPressureWater]);
	//([kJ/kg][kg/m3] + [kJ/kg][kg/m3])*[m3] - ([kJ/(kg*Kelvin)][Celsius][kg/m3])*(m3) ->...->[kJ]
	T EnergiaMax = volume*
	(
	 fPorosityRock*(EnthalpyOil(Temp)*DensityOil(Temp)*cellstate[ESaturationOil]
					+EnthalpyWater(Temp)*DensityWater(Temp)*(cellstate[ESaturationWater]+cellstate[ESaturationSteam]))
	 +(1. - fPorosityRock)*fSpecificHeatRock*fDensityRock*cellstate[ETemperature]
	 );
		
	// total energy -> ETotalEnergy
	//ETotalEnergy -> ([kJ/kg]*[kg/m3])*[m3] - ([kJ/(kg*Kelvin)]*[kg/m3]*[m3]*[Celsius]) -> [kJ] 
	REAL InitialTotalEnergy = fPorosityRock*(EnthalpyOil(initialstate[ETemperature])*DensityOil(initialstate[ETemperature])*initialstate[ESaturationOil]
											 +EnthalpyWater(initialstate[ETemperature])*DensityWater(initialstate[ETemperature])*initialstate[ESaturationWater]
											 +EnthalpySteam(initialstate[ETemperature])*DensitySteam(initialstate[ETemperature])*initialstate[ESaturationSteam])*volume
	+(1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*volume*initialstate[ETemperature];
	
	T TotalEnergy = fPorosityRock*(EnthalpyOil(cellstate[ETemperature])*DensityOil(cellstate[ETemperature])*cellstate[ESaturationOil]
																 +EnthalpyWater(cellstate[ETemperature])*DensityWater(cellstate[ETemperature])*cellstate[ESaturationWater]
																 +EnthalpySteam(cellstate[ETemperature])*DensitySteam(cellstate[ETemperature])*cellstate[ESaturationSteam])*volume
					+(1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*volume*cellstate[ETemperature];
	
	// conservation of energy internal contribution -> ETemperature -> EPhaseChange
	//[kJ]
	cellresidual[EEnergyCons] = TotalEnergy-InitialTotalEnergy;
	cellresidual[EEnergyCons] += leftflux[TPBrSteamFlux::EEnergyFlux]-rightflux[TPBrSteamFlux::EEnergyFlux];
	
	if(cellstate[ESaturationSteam] == T(0.) && TotalEnergy < EnergiaMax)
	{
		// Energy < EMin
		cellresidual[EZeroSaturation] = cellstate[ESaturationSteam];
	} 
	// the energy is larger or equal to the maximum energy
	else
	{
		cellresidual[EZeroSaturation] = cellstate[ETemperature] - Temp/*CondensationTemperature(state)*/;
	}
	
	
	// more steam condensed than steam present in the cell
	if (cellstate[ESaturationSteam] < T(0.)) {
		// Energy < EMin
		cellresidual[EZeroSaturation] = cellstate[ESaturationSteam];
	}
}


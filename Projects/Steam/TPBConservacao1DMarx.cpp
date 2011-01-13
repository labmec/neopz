/*
 *  TPBConservacao1D.cpp
 *  IP3D_v4
 *
 *  Created by Philippe Devloo on 27/03/10.
 *  Copyright 2010 UNICAMP. All rights reserved.
 *
 */

#include "TPBConservacao1DMarx.h"
//#include "ThermalMethodsTables.h"
//#include "PropertiesTable.h"
//#include "pzseqsolver.h"
#include <math.h>

TPBrCellMarx::TPBrCellMarx() : fInitialState(NUMVARS,0.)
{
	fInitialState[ESaturationOil] = 1.;
	fResidualOil = 0.01;
}

TPBrCellMarx::~TPBrCellMarx() {
}

void TPBrCellMarx::SetMaterialProperty(REAL materialpermeability, PhysicalProperties &rock)
{
	fMaterialPermeability = materialpermeability;
	fPorosityRock = rock.porosity;
	fDensityRock = rock.density;
	fSpecificHeatRock = rock.specificheat;
}

void TPBrCellMarx::SetGeometry(REAL cellvolume, REAL leftarea, REAL rightarea, REAL cellsize)
{
	fCellVolume = cellvolume;
	fLeftArea = leftarea;
	fRightArea = rightarea;
	fCellSize = cellsize;
}

void TPBrCellMarx::SetInjectionState(REAL pressurewater, TPZVec<REAL> &massflux, TPZManVector<REAL> &leftstate)
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

void TPBrCellMarx::SetCellState(REAL pressurewater, TPZVec<REAL> &saturation, REAL temperature, REAL deltIni ){
	
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
	fInitialState[EDeltaT] = deltIni;
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

//template<class T>
//void  TPBrCellMarx::RazaoVexpVcel(TPZVec<T> &state, REAL &razao){
//	T densidAgua = DensityWater(state[ETemperature]);
//	T quantAguaExp = state[EMassFluxWater]*state[EDeltaT];
//	REAL VexpAgua = quantAguaExp/densidAgua;
//	REAL Vcel = fCellVolume;
//	razao = VexpAgua/Vcel;
//}

template<class T>
// by Agnaldo
//uso a formula dada na tese SteamFlood, pagina 34
//pressuresteam-> Pa
//[Celsius]
T TPBrCellMarx::TemperatureSaturatedSteam(T pressuresteam){
	
	T  val_log, temp;
	T temp_de_saturac;
	val_log = log(pressuresteam*0.0001450377438972831);
	temp=561.435 + 33.8866*val_log + 2.18893*(val_log*val_log) + 0.0808998*(val_log*val_log*val_log) +
	0.0342030*(val_log*val_log*val_log*val_log);
	temp_de_saturac = (temp-32. - 459.67)/1.8;
	
	return temp_de_saturac;
}

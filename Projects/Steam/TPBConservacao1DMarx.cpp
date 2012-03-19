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
#include "PropertiesTable.h"
#include "tpbrsteammesh.h"
#include "pzlog.h"
#include "pzseqsolver.h"
#include "tpbrthermaldisc.h"
#include "tpbrsolutionlist.h"
#include <math.h>

#ifdef _AUTODIFF

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

//WaterDataInStateOfSaturation waterdata;
OilData oildata;

void ScaleFactor(TPZFMatrix<REAL> &tangentmatrix, TPZFMatrix<REAL> &residualmatrix, TPZManVector<REAL> &scalevalues, TPZManVector<REAL> &statescalevalues );
void ScaleFactorSol(TPZFMatrix<REAL> &residualmatrix, TPZManVector<REAL> &scalevalues );

/// compute the flux and energy timestepping with step delt
void FluxEvolution(REAL tinlet, REAL delt, REAL Tfinal, const std::string &fluxfilename, const std::string &energyfilename);

/// compute the flux and energy as the domain expands at a linear rate
void ExpandingDomain(REAL tinlet, REAL DADt, REAL deltatime, REAL TimeFinal, const std::string &fluxfilename, const std::string &energyfilename);

int main()
{
    
#ifdef LOG4CXX
	InitializePZLOG();
#endif
    /*
	int numcells = 1;
	REAL temperature = 98.;
	REAL pressure = 2.e6;
	REAL WellRadius = 0.13;
	REAL ReservoirRadius = 100.;
	REAL oilsaturation = 0.7;
	
	TPBrSteamMesh mesh(numcells,temperature,pressure,WellRadius,ReservoirRadius,oilsaturation);
	mesh.TimeStep(10.);
	return 0;
	*/

	REAL tinlet(100.);
    REAL deltatime = 1.;
    REAL timefinal(20.);
//    FluxEvolution(tinlet, deltatime, timefinal, "fluxdelt1.txt","energy1.txt");
//    FluxEvolution(tinlet, deltatime/10., timefinal, "fluxdelt01.txt","energy01.txt");
    REAL DADt(1.);
    ExpandingDomain(tinlet, DADt, deltatime, timefinal, "fluxdelt1.txt","energy1.txt");
    ExpandingDomain(tinlet, DADt, deltatime/10., timefinal, "fluxdelt01.txt","energy01.txt");
	REAL domainsize = 100.;
	int nelements = 50;
	REAL cp = 1.;
	REAL K = 1.;
    REAL density = 1.;
	REAL initialtemp = 100.;
	TPBRThermalDiscretization discrete(domainsize,nelements,cp,K,density,initialtemp);
	TPZFNMatrix<101> sol(nelements+1,1,initialtemp), nextsol(nelements+1,1,0.);
	discrete.SetTimeStep(1.);
	discrete.ComputeStiffness();
	REAL flux1,flux2;
    REAL energy1, energy2;
    REAL dQdT = discrete.DQDT();
	discrete.NextSolution(100., sol,nextsol,flux1);
    energy1 = discrete.Energy(nextsol);
	discrete.NextSolution(101., sol,nextsol,flux2);
    energy2 = discrete.Energy(nextsol);
    REAL residual = flux2-flux1-dQdT;
    std::cout << "flux2 " << flux2 << " flux1 " << flux1 << " dQdT " << dQdT << " residual " << residual << std::endl;
    std::cout << "energy1 " << energy1 << " energy2 " << energy2 << std::endl;
	nextsol.Print("Next Solution", std::cout);
	discrete.NextSolution(1., nextsol,nextsol,flux1);
	nextsol.Print("Next Solution", std::cout);
	return 0;
// #warning "this should be resolved"
}

void FluxEvolution(REAL tinlet, REAL delt, REAL Tfinal, const std::string &fluxfilename, const std::string &energyfilename)
{
	REAL domainsize = 100.;
	int nelements = 100;
	REAL cp = 1.;
	REAL K = 1.;
    REAL density = 1.;
	REAL initialtemp = 0.;
	TPBRThermalDiscretization discrete(domainsize,nelements,cp,K,density,initialtemp);
	TPZFNMatrix<11> sol(nelements+1,1,0.);
    //sol(0,0) = tinlet;
	discrete.SetTimeStep(delt);
	discrete.ComputeStiffness();
    std::ofstream outflux(fluxfilename.c_str());
    std::ofstream outenergy(energyfilename.c_str());
    REAL flux;
    REAL t;
    for (t=0; t<= Tfinal; t+=delt) {
        discrete.NextSolution(tinlet, sol, sol, flux);
        outflux << t+delt << " " << flux <<  std::endl;
        outenergy << t+delt << " " << discrete.Energy(sol) <<  std::endl;
    }
    
}

/// compute the flux and energy as the domain expands at a linear rate
void ExpandingDomain(REAL tinlet, REAL DADt, REAL deltatime, REAL TimeFinal, const std::string &fluxfilename, const std::string &energyfilename)
{
    TPBRSolutionList locallist;
    REAL domainsize = 100.;
	int nelements = 100;
	REAL cp = 1.;
	REAL K = 1.;
    REAL density = 1.;
	REAL initialtemp = 0.;
	TPBRThermalDiscretization discrete(domainsize,nelements,cp,K,density,initialtemp);
    locallist.SetDiscretization(discrete);
    std::ofstream outflux(fluxfilename.c_str());
    std::ofstream outenergy(energyfilename.c_str());
    REAL flux, DQDT;
    REAL t;
    for (t=0; t<= TimeFinal; t+=deltatime) 
    {
        TPBRThermalSolution sol(deltatime*DADt);
        locallist.AddSolution(sol);
        locallist.AdvanceAllSolutions(deltatime, tinlet, flux, DQDT, true);
        REAL energy = locallist.Energy();
        outflux << t+deltatime << " " << flux <<  std::endl;
        outenergy << t+deltatime << " " << energy  <<  std::endl;
    }


}

void TPBrCellMarx::SetInjectionState(REAL pressurewater, TPZVec<REAL> &massflux, TPZManVector<REAL> &leftstate)
{
//	REAL PI = 4*atan(1.);	
	TPBrCellMarx first;
		
	//----------------- dados de entrada ------------------------------
	//dados numerico
//	REAL TimeStep =10.;//(1260. s = tempo para atingir a Energia máxima)
	//REAL TimeStep_target =1400.;
	REAL temperature = TemperatureSaturation(pressurewater);
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


void TPBrCellMarx::ExtractMatrix(TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > &input, TPZFMatrix<REAL> &output)
{
	output.Resize(NUMVARS, NUMVARS);
	int i,j;
	for (i=0; i<NUMVARS; i++) {
		for (j=0; j<NUMVARS; j++) {
			output(i,j) = input[i].d(j);
		}
	}
}

void TPBrCellMarx::ExtractMatrix(TPZManVector<REAL> &input, TPZFMatrix<REAL> &output)
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

#else
int main()
{
    return 0;
}
// Nothing is compiled if _AUTODIFF isnt defined
#endif
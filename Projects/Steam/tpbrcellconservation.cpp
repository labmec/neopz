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
#include "PropertiesTable.h"

//WaterDataInStateOfSaturation waterdata;
static OilData oildata;

#ifdef _AUTODIFF

#include "fadType.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("br.steam.cell"));
#endif


REAL TPBrCellConservation::fPorosityRock = 0.196;
REAL TPBrCellConservation::fDensityRock = 1.055;
REAL TPBrCellConservation::fSpecificHeatRock = 0.766;
REAL TPBrCellConservation::fResidualOil = 0.0;


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
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "left flux " << leftfad << std::endl;
        sout << "cell state " << cellfad << std::endl;
        sout << "right flux " << rightfad << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
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
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Cell residual Mass Oil first part " << cellresidual[EMassOil];
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	cellresidual[EMassOil] -= delt*(leftflux[TPBrSteamFlux::EMassFluxOil]-rightflux[TPBrSteamFlux::EMassFluxOil]);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Cell residual Mass Oil accumulated " << cellresidual[EMassOil];
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	// conservation of mass of water in the cell ->ESaturationWater
	// [m3]*[Kg/m3] -> [kg]
	cellresidual[EMassWater] = volume*fPorosityRock*(DensityWater(cellstate[ETemperature])*cellstate[ESaturationWater]-DensityWater(initialstate[ETemperature])*initialstate[ESaturationWater])+cellstate[EPhaseChange];
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Cell residual Mass Water first part " << cellresidual[EMassWater] << std::endl;
        sout << "left flux " << leftflux[TPBrSteamFlux::EMassFluxWater] << std::endl;
        sout << "right flux " << rightflux[TPBrSteamFlux::EMassFluxWater] << std::endl;
        sout << "phase change " << cellstate[EPhaseChange] << std::endl;
        T cellvolprev = volume*fPorosityRock*DensityWater(initialstate[ETemperature])*initialstate[ESaturationWater];
        T cellvolnext = volume*fPorosityRock*DensityWater(cellstate[ETemperature])*cellstate[ESaturationWater];
        sout << "water volume before " << cellvolprev << std::endl;
        sout << "water volume next " << cellvolnext << std::endl;
        sout << "saturation before " << initialstate[ESaturationWater] << std::endl;
        sout << "saturation next " << cellstate[ESaturationWater] << std::endl;
        sout << "Density before " << DensityWater(initialstate[ETemperature]) << std::endl;
        sout << "Density after " << DensityWater(cellstate[ETemperature]) << std::endl;
        sout << "temperature before " << initialstate[ETemperature] << std::endl;
        sout << "temperature next " << cellstate[ETemperature] << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	cellresidual[EMassWater] -= delt*(leftflux[TPBrSteamFlux::EMassFluxWater]-rightflux[TPBrSteamFlux::EMassFluxWater]);

#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Cell residual Mass Water accumulated " << cellresidual[EMassWater];
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	// conservation of mass of steam in the cell ->ESaturationSteam
	// [m3]*[Kg/m3] -> [kg]
	cellresidual[EMassSteam] = volume*fPorosityRock*(DensitySteam(cellstate[ETemperature])*cellstate[ESaturationSteam]-DensitySteam(initialstate[ETemperature])*initialstate[ESaturationSteam])-cellstate[EPhaseChange];
	cellresidual[EMassSteam] -= delt*(leftflux[TPBrSteamFlux::EMassFluxSteam]-rightflux[TPBrSteamFlux::EMassFluxSteam]);
	
	
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
    T press = cellstate[EPressureWater]+T(TPBrScales::fReferencePressure);
	T Temp = TemperatureSaturation(press);
	//([kJ/kg][kg/m3] + [kJ/kg][kg/m3])*[m3] - ([kJ/(kg*Kelvin)][Celsius][kg/m3])*(m3) ->...->[kJ]
	T EnergiaMax = volume*
	(
	 fPorosityRock*(EnthalpyOil(Temp)*DensityOil(Temp)*cellstate[ESaturationOil]
					+EnthalpyWater(Temp)*DensityWater(Temp)*(cellstate[ESaturationWater]+cellstate[ESaturationSteam]))
	 +(((REAL)1.0) - fPorosityRock)*fSpecificHeatRock*fDensityRock*cellstate[ETemperature]
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
					+(((REAL)1.0)-fPorosityRock)*fSpecificHeatRock*fDensityRock*volume*cellstate[ETemperature];
	
	// conservation of energy internal contribution -> ETemperature -> EPhaseChange
	//[kJ]
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Initial total energy " << InitialTotalEnergy << std::endl;
        sout << "Next total energy " << TotalEnergy << std::endl;
        sout << "delt " << delt << std::endl;
        T deltInlet = delt*leftflux[TPBrSteamFlux::EEnergyFlux];
        sout << "Inlet energy flux " << deltInlet << std::endl;
        T deltright = delt*rightflux[TPBrSteamFlux::EEnergyFlux];
        sout << "Right energy flux " << deltright << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	cellresidual[EEnergyCons] = TotalEnergy-InitialTotalEnergy;
	cellresidual[EEnergyCons] -= delt*(leftflux[TPBrSteamFlux::EEnergyFlux]-rightflux[TPBrSteamFlux::EEnergyFlux]);
	
	if(cellstate[ESaturationSteam] == T(0.) && TotalEnergy < EnergiaMax)
	{
		// Energy < EMin
		cellresidual[EZeroSaturation] = cellstate[ESaturationSteam];
	} 
	// the energy is larger or equal to the maximum energy
	else
	{
		cellresidual[EZeroSaturation] = cellstate[ETemperature] - Temp;//CondensationTemperature(state)
	}
	
	
	// more steam condensed than steam present in the cell
	if (cellstate[ESaturationSteam] < T(0.)) {
		// Energy < EMin
		cellresidual[EZeroSaturation] = cellstate[ESaturationSteam];
	}
 
}

/*
static REAL fPorosityRock;
static REAL fDensityRock;
static REAL fSpecificHeatRock;
static REAL fResidualOil;
*/

/// Print the data of the cell
void TPBrCellConservation::Print(std::ostream &out)
{
    out << "Cell conservation\n";
    out << "Rock porosity " << fPorosityRock << std::endl;
    out << "Rock density " << fDensityRock << std::endl;
    out << "Specific heat " << fSpecificHeatRock << std::endl;
    out << "Residual saturation oil "  << fResidualOil << std::endl;
}
/*
enum ECellEq { EMassWater, EMassOil, EMassSteam, ECapillaryPressureWO, ECapillaryPressureOS, ESumSaturation, EEnergyCons,
    EZeroSaturation
};

/// State variables associated with a cell
enum ECellState {
    ESaturationWater, ESaturationOil, ESaturationSteam, ETemperature, EPressureWater, EPressureSteam, EPressureOil,EPhaseChange
};
*/
/// associate scale factors with the equations and state variables
void TPBrCellConservation::Scales(TPZVec<REAL> &eqscales, TPZVec<REAL> &statescales)
{
    eqscales.Resize(NumCellEq);
    statescales.Resize(NumCellEq);
    eqscales[EMassWater] = 1.;
    eqscales[EMassOil] = 1.;
    eqscales[EMassSteam] = 1.;
    eqscales[ECapillaryPressureWO] = TPBrScales::fPressureScale;
    eqscales[ECapillaryPressureOS] = TPBrScales::fPressureScale;
    eqscales[ESumSaturation] = 1.;
    eqscales[EEnergyCons] = TPBrScales::fEnergyScale;
    eqscales[EZeroSaturation] = 1.;
    statescales[ESaturationWater] = 1.;
    statescales[ESaturationOil] = 1.;
    statescales[ESaturationSteam] = 1.;
    statescales[ETemperature] = 1.;
    statescales[EPressureWater] = TPBrScales::fPressureScale;
    statescales[EPressureSteam] = TPBrScales::fPressureScale;
    statescales[EPressureOil] = TPBrScales::fPressureScale;
    statescales[EPhaseChange] = 1.;
}

template<class T>
T TPBrCellConservation::DensityWater(T temp)//[kg/m3]
{
	return waterdata.getSaturationStateDensityToLiquidWater(temp);
}
template<class T>
T TPBrCellConservation::DensitySteam(T temp)//[kg/m3]
{
	return waterdata.getSaturationStateDensityToSteam(temp);
}

//template<>
//REAL TPBrCellConservation::TemperatureSaturation<REAL>(REAL pressuresteam); //[Celsius]

//template<>
//REAL TPBrCellConservation::TemperatureSaturation(REAL pressuresteam) //[Celsius]
//{
 //   double press1000 = pressuresteam/1000.;
   // double temp_de_saturac = waterdata.getSaturationStateTemperature(press1000);
//	return temp_de_saturac;
//}

template<class T>
T TPBrCellConservation::EnthalpyWater(T temperature)//[kJ/kg]
{
	return waterdata.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
}
template<class T>
T TPBrCellConservation::EnthalpySteam(T temperature)//[kJ/kg]
{
	return waterdata.getSaturationStateSpecificEnthalpyToSteam(temperature);
}

// metodos para recuperar os dados tabulados em funcao
template<class T>
T TPBrCellConservation::DensityOil(T temp)//[kg/m3]
{
	return oildata.getDensityToOil(temp);
}
template<class T>
T TPBrCellConservation::EnthalpyOil(T temperature)//[kJ/kg]
{
	return  oildata.getSpecificHeatToOil(temperature);	
}

/// limit the correction of the Newton Iteration
REAL TPBrCellConservation::LimitRange(REAL scale,TPZVec<REAL> &cellstate,TPZVec<REAL> &cellcorrection)
{
    return scale;
}

/// Compute the total energy of the cell
REAL TPBrCellConservation::Energy(TPZVec<REAL> &cellstate, REAL volume)
{
    REAL TotalEnergy = fPorosityRock*(EnthalpyOil(cellstate[ETemperature])*DensityOil(cellstate[ETemperature])*cellstate[ESaturationOil]
                                   +EnthalpyWater(cellstate[ETemperature])*DensityWater(cellstate[ETemperature])*cellstate[ESaturationWater]
                                   +EnthalpySteam(cellstate[ETemperature])*DensitySteam(cellstate[ETemperature])*cellstate[ESaturationSteam])*volume
                                +(1.-fPorosityRock)*fSpecificHeatRock*fDensityRock*volume*cellstate[ETemperature];
    return TotalEnergy;
    
}


#endif
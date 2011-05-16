/*
 *  tpbrsteamflux.cpp
 *  PZ
 *
 *  Created by Philippe Devloo on 3/8/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */

#include "tpbrsteamflux.h"
#include "tpbrcellconservation.h"

#include "ThermalMethodsTables.h"

#ifdef _AUTODIFF

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("br.steam.steamflux"));
#endif

extern WaterDataInStateOfSaturation waterdata;
extern OilData oildata;


REAL TPBrSteamFlux::fFarfieldPressureOil = 2.e6;
REAL TPBrSteamFlux::fFarfieldPressureWater = 2.e6;
REAL TPBrSteamFlux::fFarfieldPressureSteam = 2.e6;
REAL TPBrSteamFlux::fFarfieldTemperature = 80;
REAL TPBrSteamFlux::fFarfieldSaturationOil = 0.7;
REAL TPBrSteamFlux::fFarfieldSaturationWater = 0.3;
REAL TPBrSteamFlux::fFarfieldSaturationSteam = 0.;

REAL TPBrSteamFlux::fInletEnergyFlux = 500.; //[KJ/s]
REAL TPBrSteamFlux::fInletMassFlux = 0.7; //[Kg/s];


/// empty constructor initializing all variables
TPBrSteamFlux::TPBrSteamFlux()
{
	fMaterialPermeability = 0.8e-12;
}

/// calcula a contribuicao para a matriz de rigidez
void TPBrSteamFlux::CalcStiff(TPZVec<REAL> &leftstate, TPZVec<REAL> &rightstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
			   TPZFMatrix &ek, TPZFMatrix &ef)
{
	const int totaleq = 2*TPBrCellConservation::NumCellEq+TPBrSteamFlux::NumFluxEq;
	TPZManVector<TFad<totaleq,REAL> , TPBrCellConservation::NumCellEq> leftcellfad(TPBrCellConservation::NumCellEq),rightcellfad(TPBrCellConservation::NumCellEq);
	TPZManVector<TFad<totaleq,REAL> , TPBrSteamFlux::NumFluxEq> interfacefad(TPBrSteamFlux::NumFluxEq);
	
	TPBrCellConservation::Initialize<totaleq>(leftstate,leftcellfad,0);
	TPBrSteamFlux::Initialize<totaleq>(interfacestate,interfacefad,TPBrCellConservation::NumCellEq);
	TPBrCellConservation::Initialize<totaleq>(rightstate,rightcellfad,TPBrCellConservation::NumCellEq+TPBrSteamFlux::NumFluxEq);
	
	TPZManVector<TFad<totaleq,REAL> , NumFluxEq> cellresidualfad(NumFluxEq); 
	
	FluxResidual(leftcellfad, interfacefad, rightcellfad, delx, area, delt, cellresidualfad );
	
	ek.Redim(NumFluxEq, totaleq);
	ef.Redim(NumFluxEq, 1);
	int i,j;
	for (i=0; i<NumFluxEq; i++) 
	{
		ef(i,0) = cellresidualfad[i].val();
		for (j=0; j<totaleq; j++) {
			ek(i,j) = cellresidualfad[i].d(j);
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		ek.Print("Flux stiffness",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

/// calcula a contribuicao para a matriz de rigidez das equacoes de entrada
void TPBrSteamFlux::InletCalcStiff(TPZVec<REAL> &inletstate, TPZVec<REAL> &rightstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
					TPZFMatrix &ek, TPZFMatrix &ef)
{
	const int totaleq = NumInletVars+NumFluxEq+TPBrCellConservation::NumCellEq;
	TPZManVector<TFad<totaleq,REAL> , NumInletVars> inletfad(NumInletVars);
	TPZManVector<TFad<totaleq,REAL> , TPBrSteamFlux::NumFluxEq> interfacefad(TPBrSteamFlux::NumFluxEq);
	TPZManVector<TFad<totaleq,REAL> , TPBrCellConservation::NumCellEq> rightcellfad(TPBrCellConservation::NumCellEq);
	
	TPBrSteamFlux::InitializeInlet<totaleq>(inletstate,inletfad,0);
	TPBrSteamFlux::Initialize<totaleq>(interfacestate,interfacefad,NumInletVars);
	TPBrCellConservation::Initialize<totaleq>(rightstate,rightcellfad,NumInletVars+TPBrSteamFlux::NumFluxEq);
#ifdef LOG4CXX
	{
		std::stringstream sout;
        sout << "before calling inlet flux residual\n";
		sout << "inletfad " << inletfad << std::endl;
		sout << "interfacefad " << interfacefad << std::endl;
        sout << "rightcellfad " << rightcellfad << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	TPZManVector<TFad<totaleq,REAL> , NumFluxEq+NumInletVars> cellresidualfad(NumFluxEq+NumInletVars); 
	
	InletFluxResidual(inletfad, interfacefad, rightcellfad, delx, area, delt, cellresidualfad );
	
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "cellresidual " << cellresidualfad;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	ek.Redim(NumFluxEq+NumInletVars, totaleq);
	ef.Redim(NumFluxEq+NumInletVars, 1);
	int i,j;
	for (i=0; i<NumFluxEq+NumInletVars; i++) 
	{
		ef(i,0) = cellresidualfad[i].val();
		for (j=0; j<totaleq; j++) {
			ek(i,j) = cellresidualfad[i].d(j);
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		ek.Print("Inlet stiffness",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

/// calcula a contribuicao para a matriz de rigidez das equacoes de entrada
void TPBrSteamFlux::OutletCalcStiff(TPZVec<REAL> &leftstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
					 TPZFMatrix &ek, TPZFMatrix &ef)
{
	const int totaleq = NumFluxEq+TPBrCellConservation::NumCellEq;
	TPZManVector<TFad<totaleq,REAL> , TPBrCellConservation::NumCellEq> leftcellfad(TPBrCellConservation::NumCellEq);
	TPZManVector<TFad<totaleq,REAL> , TPBrSteamFlux::NumFluxEq> interfacefad(TPBrSteamFlux::NumFluxEq);
	
	TPBrCellConservation::Initialize<totaleq>(leftstate,leftcellfad,0);
	TPBrSteamFlux::Initialize<totaleq>(interfacestate,interfacefad,TPBrCellConservation::NumCellEq);
	
	TPZManVector<TFad<totaleq,REAL> , NumFluxEq> cellresidualfad(NumFluxEq); 
	
	OutletFluxResidual(leftcellfad, interfacefad, delx, area, delt, cellresidualfad );
	
	ek.Redim(NumFluxEq, totaleq);
	ef.Redim(NumFluxEq, 1);
	int i,j;
	for (i=0; i<NumFluxEq; i++) 
	{
		ef(i,0) = cellresidualfad[i].val();
		for (j=0; j<totaleq; j++) {
			ek(i,j) = cellresidualfad[i].d(j);
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		ek.Print("Outlet flux stiffness",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}


/// Incorporate the partial derivatives in the state variables
template<int N>
void TPBrSteamFlux::InitializeInlet(TPZVec<REAL> &state, TPZVec<TFad<N,REAL> > &fadstate, int offset)
{
/*
	enum EInletVars { EInletPressure, EInletSteamSaturation };
*/
	
	fadstate[EInletPressure] = state[EInletPressure];
	fadstate[EInletPressure].fastAccessDx(EInletPressure+offset) = 1.;
	fadstate[EInletSteamSaturation] = state[EInletSteamSaturation];
	fadstate[EInletSteamSaturation].fastAccessDx(EInletSteamSaturation+offset) = 1.;
	fadstate[EInletTemperature] = state[EInletTemperature];
	fadstate[EInletTemperature].fastAccessDx(EInletTemperature+offset) = 1.;
}


/*
enum ESteamFluxEq {
	EMassFluxWaterEq, EMassFluxSteamEq, EMassFluxOilEq, EDarcyVelocityWaterEq, EDarcyVelocitySteamEq, EDarcyVelocityOilEq, EEnergyFluxEq
};
enum ESteamFluxVars {
	EMassFluxWater, EMassFluxSteam, EMassFluxOil, EDarcyVelocityWater, EDarcyVelocitySteam, EDarcyVelocityOil, EEnergyFlux
};
*/


/// calcula o fluxo entre duas celulas
template<class T>
void TPBrSteamFlux::FluxResidual(TPZVec<T> &leftstate, TPZVec<T> &rightstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual)
{
	TPZManVector<T> saturation(3), relatpermeability(3);
	T temperature;
	if (leftstate[TPBrCellConservation::EPressureWater] > rightstate[TPBrCellConservation::EPressureWater]) 
	{
		saturation[EOil]=leftstate[TPBrCellConservation::ESaturationOil];
		saturation[EWater]=leftstate[TPBrCellConservation::ESaturationWater];
		saturation[ESteam]=leftstate[TPBrCellConservation::ESaturationSteam];
		temperature = leftstate[TPBrCellConservation::ETemperature];
	}
	else 
	{
		saturation[EOil]=rightstate[TPBrCellConservation::ESaturationOil];
		saturation[EWater]=rightstate[TPBrCellConservation::ESaturationWater];
		saturation[ESteam]=rightstate[TPBrCellConservation::ESaturationSteam];
		temperature = rightstate[TPBrCellConservation::ETemperature];
	}

	ComputeRelativePermeability(saturation, relatpermeability);
	// relation between dp/dx and the darcy velocity ->EDarcyVelocityOil
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s]  
	fluxresidual[EDarcyVelocityOilEq] = (leftstate[TPBrCellConservation::EPressureOil]-rightstate[TPBrCellConservation::EPressureOil])*fMaterialPermeability*relatpermeability[EOil]/ViscosityOil(temperature)
							-interfacestate[EDarcyVelocityOil]*delx;
	T DarcyVelocityOil = (
						  (leftstate[TPBrCellConservation::EPressureOil]-rightstate[TPBrCellConservation::EPressureOil])*fMaterialPermeability*relatpermeability[EOil]/ViscosityOil(temperature)
						  )/delx;
	//----deletar depois
	T viscoil = ViscosityOil(temperature);
	T KrOil = relatpermeability[EOil];
	T pressaoIniOil =leftstate[TPBrCellConservation::EPressureOil];
	T presscal6 = fMaterialPermeability*relatpermeability[EOil]/viscoil;
	//---------------
	
	// relation between dp/dx and the darcy velocity -> EDarcyVelocityWater
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s]
	fluxresidual[EDarcyVelocityWaterEq] = (leftstate[TPBrCellConservation::EPressureWater]-rightstate[TPBrCellConservation::EPressureWater])*fMaterialPermeability*relatpermeability[EWater]/ViscosityWater(temperature)
									-interfacestate[EDarcyVelocityWater]*delx;
	T DarcyVelocityWater = (
							(leftstate[TPBrCellConservation::EPressureWater]-rightstate[TPBrCellConservation::EPressureWater])*fMaterialPermeability*relatpermeability[EWater]/ViscosityWater(temperature)
							)/delx;
	//----deletar depois
	T viscWater = ViscosityWater(temperature);
	T KrWater = relatpermeability[EWater];
	T pressaoIniWater =leftstate[TPBrCellConservation::EPressureWater];
	//--------------
	
	// relation between dp/dx and the darcy velocity -> EDarcyVelocitySteam
	//[Pa]*[m2]/[Pa*sec] - [m/sec]*[m] -> [m2/s] - [m2/s] -> ... ->[m/s]
	fluxresidual[EDarcyVelocitySteamEq] = (leftstate[TPBrCellConservation::EPressureSteam]-rightstate[TPBrCellConservation::EPressureSteam])*fMaterialPermeability*relatpermeability[ESteam]/ViscositySteam(temperature)
									-interfacestate[EDarcyVelocitySteam]*delx;
	T DarcyVelocitySteam = (
							(leftstate[TPBrCellConservation::EPressureSteam]-rightstate[TPBrCellConservation::EPressureSteam])*fMaterialPermeability*relatpermeability[ESteam]/ViscositySteam(temperature)
							)/delx;
	//----deletar depois
	T viscSteam = ViscositySteam(temperature);
	T KrSteam = relatpermeability[ESteam];
	T pressaoIniSteam =leftstate[TPBrCellConservation::EPressureSteam];
	//--------------
	
	// relation between massflow rate and darcy velocity ->EMassFluxOil
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [gk/s]
	fluxresidual[EMassFluxOilEq] = interfacestate[EMassFluxOil]-DarcyVelocityOil*DensityOil(temperature)*area*delt;
	
	// relation between massflow rate and darcy velocity ->EMassFluxWater
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [gk/s]
	fluxresidual[EMassFluxWaterEq] = interfacestate[EMassFluxWater]-DarcyVelocityWater*DensityWater(temperature)*area*delt;
	
	// relation between massflow rate and darcy velocity ->EMassFluxSteam
	//[gk/s] - [m/s]*[gk/m3]*[m2] -> [gk/s]
	fluxresidual[EMassFluxSteamEq] = interfacestate[EMassFluxSteam]-DarcyVelocitySteam*DensitySteam(temperature)*area*delt;
	
	T vEnthalpyOil = EnthalpyOil(temperature);
	T vEnthalpyWater = EnthalpyWater(temperature);
	T vEnthalpySteam = EnthalpySteam(temperature);
	//( [kJ/kg]*[m/s]*[kg/m3] )*[m2]*[s] -> [kJ]
	fluxresidual[EEnergyFluxEq] = (vEnthalpyOil*DarcyVelocityOil*DensityOil(temperature)
								  +vEnthalpyWater*DarcyVelocityWater*DensityWater(temperature)
								  +vEnthalpySteam*DarcyVelocitySteam*DensitySteam(temperature)
								  )*area*delt;
	
	
}

/*
 /// State variables associated with a cell
 enum ECellState {
 ESaturationWater, ESaturationOil, ESaturationSteam, ETemperature, EPressureWater, EPressureSteam, EPressureOil,EPhaseChange
 };
*/

/// compute left state as a function of the inletstate
template<class T>
void TPBrSteamFlux::ComputeLeftState(TPZVec<T> &inletstate, TPZVec<T> &leftstate)
{
	T pressure = inletstate[EInletPressure];
	T saturationSteam = inletstate[EInletSteamSaturation];
    T temperature = inletstate[EInletTemperature];

    leftstate[TPBrCellConservation::ESaturationWater] = 1.-saturationSteam;
	leftstate[TPBrCellConservation::ESaturationOil] = 0.;
    leftstate[TPBrCellConservation::ESaturationSteam] = saturationSteam;
	leftstate[TPBrCellConservation::ETemperature] = temperature;
	leftstate[TPBrCellConservation::EPressureWater] = pressure;
	leftstate[TPBrCellConservation::EPressureSteam] = pressure;
	leftstate[TPBrCellConservation::EPressureOil] = pressure;
	leftstate[TPBrCellConservation::EPhaseChange] = 0.;
	
}

/// complete residual vector as a function of the inletstate and rightstate
template<class T>
void TPBrSteamFlux::InletFluxResidual(TPZVec<T> &inletstate, TPZVec<T> &rightstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual)
{
	TPZManVector<T> leftstate(TPBrCellConservation::NumCellEq);
    ComputeLeftState(inletstate,leftstate);
#ifdef LOG4CXX
	{
		std::stringstream sout;
        sout << "Transferring the inletstate to a leftstate\n";
		sout << "\ninletstate " << inletstate;
		sout << "\nleftstate " << leftstate;
        sout << "\nrightstate " << rightstate;
        sout << "\ninterfacestate " << interfacestate;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	FluxResidual(leftstate,rightstate,interfacestate,delx,area,delt,fluxresidual);
	int i;
	for (i=NumInletVars; i<NumFluxEq+NumInletVars; i++) {
		fluxresidual[i] = fluxresidual[i-NumInletVars];
	}

	fluxresidual[EInletMassFlux] = fInletMassFlux*delt - interfacestate[EMassFluxSteam] - interfacestate[EMassFluxWater];
	fluxresidual[EInletEnergyFlux] = fInletEnergyFlux*delt - interfacestate[EEnergyFlux];
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "residual of the inlet mass flux " << fluxresidual[EInletMassFlux] << std::endl;
        sout << "residual of the inlet energy " << fluxresidual[EInletEnergyFlux] << std::endl;
        sout << "fInletMassFlux " << fInletMassFlux << std::endl;
        sout << "interfacestate[EMassFluxSteam] " << interfacestate[EMassFluxSteam] << std::endl;
        sout << "interfacestate[EMassFluxWater] " << interfacestate[EMassFluxWater] << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    T saturationtemperature = TemperatureSaturation(inletstate[EInletPressure]);
    std::stringstream sout;
    if (inletstate[EInletTemperature] > saturationtemperature) 
    {
        sout << "The inlet temperature is larger than the saturation temperature at line " << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletTemperature]-saturationtemperature;
        inletstate[EInletTemperature] = saturationtemperature;
    }
    if (inletstate[EInletSteamSaturation] < T(0.))
    {
        sout << "The steam saturation is less than zero at line " << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletSteamSaturation];
        inletstate[EInletSteamSaturation] = 0.;
    }
    if (inletstate[EInletSteamSaturation] > T(0.))
    {
        sout << "The inlet steam saturation is larger than zero " << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletTemperature]-saturationtemperature;
    }
    else if(inletstate[EInletTemperature] < saturationtemperature)
    {
        sout << "The temperature is lower than the saturation temperature " << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletSteamSaturation];
    }
    else if(fluxresidual[EInletEnergyFlux] > T(0.))
    {
        sout << "Very unusual border case steam saturation zero AND inlet temperature is steam temperature AND positive Flux balance" << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletTemperature]-saturationtemperature;
    }
    else
    {
        sout << "Very unusual border case steam saturation zero AND inlet temperature is steam temperature AND non positive Flux balance" << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletSteamSaturation];        
    }
#ifdef LOG4CXX
    sout << std::endl;
    LOGPZ_DEBUG(logger, sout.str())
#endif
}

/// complete residual vector as a function of the inletstate and rightstate
template<class T>
void TPBrSteamFlux::OutletFluxResidual(TPZVec<T> &leftstate, TPZVec<T> &interfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual)
{
/*
	static REAL fFarfieldPressureOil;
	static REAL fFarfieldPressureWater;
	static REAL fFarfieldPressureSteam;
	static REAL fFarfieldTemperature;
	static REAL fFarfieldSaturationOil;
	static REAL fFarfieldSaturationWater;
	static REAL fFarfieldSaturationSteam;
*/	
	TPZManVector<T> rightstate(TPBrCellConservation::NumCellEq);
	rightstate[TPBrCellConservation::ESaturationWater] = T(fFarfieldSaturationWater);
	rightstate[TPBrCellConservation::ESaturationOil] = T(fFarfieldSaturationOil);
	rightstate[TPBrCellConservation::ETemperature] = T(fFarfieldTemperature);
	rightstate[TPBrCellConservation::EPressureWater] = T(fFarfieldPressureWater);
	rightstate[TPBrCellConservation::EPressureSteam] = T(fFarfieldPressureSteam);
	rightstate[TPBrCellConservation::EPressureOil] = T(fFarfieldPressureOil);
	rightstate[TPBrCellConservation::EPhaseChange] = T(0.);
	FluxResidual(leftstate,rightstate,interfacestate,delx,area,delt,fluxresidual);
}

template<class T>
T TPBrSteamFlux::EnthalpyWater(T temperature)//[kJ/kg]
{
	return waterdata.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
}
template<class T>
T TPBrSteamFlux::EnthalpySteam(T temperature)//[kJ/kg]
{
	return waterdata.getSaturationStateSpecificEnthalpyToSteam(temperature);
}
template<class T>
T TPBrSteamFlux::EnthalpyOil(T temperature)//[kJ/kg]
{
	return  oildata.getSpecificHeatToOil(temperature);	
}
template<class T>
T TPBrSteamFlux::ViscosityOil(T temp)//[Pa*sec]
{
	return oildata.getDynamicViscosityToOil(temp);
}
template<class T>
T TPBrSteamFlux::ViscosityWater(T temp)//[Pa*sec]
{
	return waterdata.getSaturationStateViscosityToLiquidWater(temp);
}
template<class T>
T TPBrSteamFlux::ViscositySteam(T temp)//[Pa*sec]
{
	return waterdata.getSaturationStateViscosityToSteam(temp);
}

// metodos para recuperar os dados tabulados em funcao
template<class T>
T TPBrSteamFlux::DensityOil(T temp)//[kg/m3]
{
	return oildata.getDensityToOil(temp);
}
template<class T>
T TPBrSteamFlux::DensityWater(T temp)//[kg/m3]
{
	return waterdata.getSaturationStateDensityToLiquidWater(temp);
}
template<class T>
T TPBrSteamFlux::DensitySteam(T temp)//[kg/m3]
{
	return waterdata.getSaturationStateDensityToSteam(temp);
}

template<class T>
//dimensionless
void TPBrSteamFlux::ComputeRelativePermeability(TPZManVector<T> &saturation,TPZManVector<T> &relativepermeability) {
	relativepermeability = saturation;
}

/*	REAL fMaterialPermeability;
 
 static REAL fFarfieldPressureOil;
 static REAL fFarfieldPressureWater;
 static REAL fFarfieldPressureSteam;
 static REAL fFarfieldTemperature;
 static REAL fFarfieldSaturationOil;
 static REAL fFarfieldSaturationWater;
 static REAL fFarfieldSaturationSteam;
 
 static REAL fInletEnergyFlux; //[KJ/s]
 static REAL fInletMassFlux; //[Kg/s];

 */

void TPBrSteamFlux::Print(std::ostream &out)
{
    out << "Steam flux data\n";
    out << "Material permeability " << fMaterialPermeability << std::endl;
    out << "Farfield pressure oil " << fFarfieldPressureOil << std::endl;
    out << "Farfield pressure water " << fFarfieldPressureWater << std::endl;
    out << "Farfield pressure steam " << fFarfieldPressureSteam << std::endl;
    out << "Farfield temperature " << fFarfieldTemperature << std::endl;
    out << "Farfield saturation oil " << fFarfieldSaturationOil << std::endl;
    out << "Farfield saturation water " << fFarfieldSaturationWater << std::endl;
    out << "Farfield saturation steam " << fFarfieldSaturationSteam << std::endl;
    out << "Energy inlet flux " << fInletEnergyFlux << std::endl;
    out << "Inlet mass flux " << fInletMassFlux << std::endl;
}

#endif
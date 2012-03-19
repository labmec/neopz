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
#include "PropertiesTable.h"

#ifdef _AUTODIFF

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("br.steam.steamflux"));
static LoggerPtr logger2(Logger::getLogger("br.steam.steamflux"));
#endif

//WaterDataInStateOfSaturation waterdata;
//extern WaterDataInStateOfSaturation waterdata;
static OilData oildata;


REAL TPBrSteamFlux::fFarfieldPressureOil = 2.e6;
REAL TPBrSteamFlux::fFarfieldPressureWater = 2.e6;
REAL TPBrSteamFlux::fFarfieldPressureSteam = 2.e6;
REAL TPBrSteamFlux::fFarfieldTemperature = 80;
REAL TPBrSteamFlux::fFarfieldSaturationOil = 0.;
REAL TPBrSteamFlux::fFarfieldSaturationWater = 1.;
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
			   TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
	const int totaleq = 2*TPBrCellConservation::NumCellEq+TPBrSteamFlux::NumFluxEq;
	TPZManVector<TFad<totaleq,REAL> , TPBrCellConservation::NumCellEq> leftcellfad(TPBrCellConservation::NumCellEq),rightcellfad(TPBrCellConservation::NumCellEq);
	TPZManVector<TFad<totaleq,REAL> , TPBrSteamFlux::NumFluxEq> interfacefad(TPBrSteamFlux::NumFluxEq);
	
	TPBrCellConservation::Initialize<totaleq>(leftstate,leftcellfad,0);
	TPBrSteamFlux::Initialize<totaleq>(interfacestate,interfacefad,TPBrCellConservation::NumCellEq);
	TPBrCellConservation::Initialize<totaleq>(rightstate,rightcellfad,TPBrCellConservation::NumCellEq+TPBrSteamFlux::NumFluxEq);
	
	TPZManVector<TFad<totaleq,REAL> , NumFluxEq> cellresidualfad(NumFluxEq); 
	
	FluxResidual(leftcellfad, rightcellfad, interfacefad, delx, area, delt, cellresidualfad );
	
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
void TPBrSteamFlux::InletCalcStiff(TPZVec<REAL> &rightstate, TPZVec<REAL> &interfacestate, REAL delx, REAL area, REAL delt, 
					TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
	const int totaleq = NumInletVars+NumFluxEq+TPBrCellConservation::NumCellEq;
//	TPZManVector<TFad<totaleq,REAL> , NumInletVars> inletfad(NumInletVars);
	TPZManVector<TFad<totaleq,REAL> , TPBrSteamFlux::NumFluxEq+TPBrSteamFlux::NumInletVars> interfacefad(TPBrSteamFlux::NumFluxEq+TPBrSteamFlux::NumInletVars);
	TPZManVector<TFad<totaleq,REAL> , TPBrCellConservation::NumCellEq> rightcellfad(TPBrCellConservation::NumCellEq);
	
	//TPBrSteamFlux::InitializeInlet<totaleq>(inletstate,inletfad,0);
	TPBrSteamFlux::InitializeInlet<totaleq>(interfacestate,interfacefad,0);
	TPBrCellConservation::Initialize<totaleq>(rightstate,rightcellfad,NumInletVars+TPBrSteamFlux::NumFluxEq);
#ifdef LOG4CXX
	{
		std::stringstream sout;
        sout << "before calling inlet flux residual\n";
//		sout << "inletfad " << inletfad << std::endl;
		sout << "interfacefad " << interfacefad << std::endl;
        sout << "rightcellfad " << rightcellfad << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	TPZManVector<TFad<totaleq,REAL> , NumFluxEq+NumInletVars> cellresidualfad(NumFluxEq+NumInletVars); 
	
	InletFluxResidual(rightcellfad, interfacefad, delx, area, delt, cellresidualfad );
	
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
					 TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
	const int totaleq = NumFluxEq+TPBrCellConservation::NumCellEq;
	TPZManVector<TFad<totaleq,REAL> , TPBrCellConservation::NumCellEq> leftcellfad(TPBrCellConservation::NumCellEq);
	TPZManVector<TFad<totaleq,REAL> , TPBrSteamFlux::NumFluxEq> interfacefad(TPBrSteamFlux::NumFluxEq);
	
	TPBrCellConservation::Initialize<totaleq>(leftstate,leftcellfad,0);
	TPBrSteamFlux::Initialize<totaleq>(interfacestate,interfacefad,TPBrCellConservation::NumCellEq);
	
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "leftstate " << leftstate << std::endl;
        sout << "leftcellfad " << leftcellfad << std::endl;
        sout << "interfacestate " << interfacestate << std::endl;
        sout << "interfacefad " << interfacefad << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
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
    
    TPZManVector<REAL,NumFluxEq> fluxstate(NumFluxEq);
    TPZManVector<TFad<N,REAL>, NumFluxEq> fadfluxstate(NumFluxEq);
    int is;
    for (is=0; is<NumFluxEq; is++) {
        fluxstate[is] = state[is+NumInletVars];
    }
    Initialize(fluxstate, fadfluxstate, NumInletVars);
    for (is=0; is<NumFluxEq; is++) {
        fadstate[is+NumInletVars] = fadfluxstate[is];
    }
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
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "leftstate " << leftstate << std::endl;
        sout << "rightstate " << rightstate << std::endl;
        sout << "interfacestate " << interfacestate << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	TPZManVector<T> saturation(3), relatpermeability(3);
	T temperature;
	if (leftstate[TPBrCellConservation::EPressureWater] >= rightstate[TPBrCellConservation::EPressureWater]) 
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
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Darcy velocity oil " << interfacestate[EDarcyVelocityOil] << std::endl;
        sout << "flux residual # " << EDarcyVelocityOilEq << " residual " << fluxresidual[EDarcyVelocityOilEq];
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
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
#ifdef LOG4CXX
    {
        using namespace std;
        std::stringstream sout;
        sout << "Data to compute the mass flux\n";
        sout << "material permeability " << fMaterialPermeability << std::endl;
        sout << "relative permeability water " << relatpermeability[EWater] << endl;
        sout << "viscosity water " << ViscosityWater(temperature) << endl;
        sout << "delx " << delx << endl;
        sout << "water density " << DensityWater(temperature) << endl;
        sout << "area " << area << endl;
        T mult = DensityWater(temperature)*area*fMaterialPermeability*relatpermeability[EWater]/ViscosityWater(temperature);
        sout << "delp multiplier coefficient " << mult << endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
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
	//[kg/s] - [m/s]*[kg/m3]*[m2] -> [kg/s]
	fluxresidual[EMassFluxOilEq] = interfacestate[EMassFluxOil]-DarcyVelocityOil*DensityOil(temperature)*area;
	
	// relation between massflow rate and darcy velocity ->EMassFluxWater
	//[kg/s] - [m/s]*[kg/m3]*[m2] -> [kg/s]
	fluxresidual[EMassFluxWaterEq] = interfacestate[EMassFluxWater]-DarcyVelocityWater*DensityWater(temperature)*area;
	
	// relation between massflow rate and darcy velocity ->EMassFluxSteam
	//[kg/s] - [m/s]*[kg/m3]*[m2] -> [kg/s]
	fluxresidual[EMassFluxSteamEq] = interfacestate[EMassFluxSteam]-DarcyVelocitySteam*DensitySteam(temperature)*area;
	
	T vEnthalpyOil = EnthalpyOil(temperature);
	T vEnthalpyWater = EnthalpyWater(temperature);
    T enthalpy2 = waterdata.getSpecificEnthalpyToLiquidWater(temperature);
	T vEnthalpySteam = EnthalpySteam(temperature);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Temperature " << temperature << std::endl;
        sout << "Darcy velocity oil " << DarcyVelocityOil << std::endl;
        sout << "Darcy velocity water " << DarcyVelocityWater << std::endl;
        sout << "Darcy velocity steam " << DarcyVelocitySteam << std::endl;
        sout << "enthalpy oil " << vEnthalpyOil << std::endl;
        sout << "enthalpy water " << vEnthalpyWater << std::endl;
        sout << "enthalpy computed differently " << enthalpy2 << std::endl;
        sout << "enthalpy steam " << vEnthalpySteam << std::endl;
        sout << "Density Oil " << DensityOil(temperature) << std::endl;
        sout << "Density Water " << DensityWater(temperature) << std::endl;
        sout << "Density Steam " << DensitySteam(temperature) << std::endl;
        sout << "area " << area << std::endl;
        T massfluxwater = DarcyVelocityWater*DensityWater(temperature)*area;
        sout << "mass flux water " << massfluxwater << std::endl;
        sout << "mass flux water by interface " << interfacestate[EMassFluxWater] << std::endl;
        T massfluxsteam = DarcyVelocitySteam*DensitySteam(temperature)*area;
        sout << "mass flux steam " << massfluxsteam << std::endl;
        sout << "mass flux steam by interface " << interfacestate[EMassFluxSteam] << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	//( [kJ/kg]*[m/s]*[kg/m3] )*[m2]*[s] -> [kJ]
	fluxresidual[EEnergyFluxEq] = interfacestate[EEnergyFlux]-(vEnthalpyOil*DarcyVelocityOil*DensityOil(temperature)
								  +vEnthalpyWater*DarcyVelocityWater*DensityWater(temperature)
								  +vEnthalpySteam*DarcyVelocitySteam*DensitySteam(temperature)
								  )*area;
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "flux residual energy " << fluxresidual[EEnergyFluxEq];
        LOGPZ_DEBUG(logger, sout.str())
    }	
#endif
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

    leftstate[TPBrCellConservation::ESaturationWater] = T(1.0)-saturationSteam;
	leftstate[TPBrCellConservation::ESaturationOil] = T(0.0);
    leftstate[TPBrCellConservation::ESaturationSteam] = saturationSteam;
	leftstate[TPBrCellConservation::ETemperature] = temperature;
	leftstate[TPBrCellConservation::EPressureWater] = pressure;
	leftstate[TPBrCellConservation::EPressureSteam] = pressure;
	leftstate[TPBrCellConservation::EPressureOil] = pressure;
	leftstate[TPBrCellConservation::EPhaseChange] = T(0.0);
	
}


/// complete residual vector as a function of the inletstate and rightstate
template<class T>
void TPBrSteamFlux::InletFluxResidual(TPZVec<T> &rightstate, TPZVec<T> &suminterfacestate, REAL delx, REAL area, REAL delt, TPZVec<T> &fluxresidual)
{
    TPZManVector<T,NumFluxEq> interfacestate(NumFluxEq);
    TPZManVector<T,NumInletVars> inletstate(NumInletVars);
    int is;
    for (is=0; is<NumInletVars; is++) {
        inletstate[is] = suminterfacestate[is];
    }
    for (is=0; is<NumFluxEq; is++) {
        interfacestate[is] = suminterfacestate[is+NumInletVars];
    }
    T difference = (inletstate[EInletPressure] - rightstate[TPBrCellConservation::EPressureWater]);
    while(difference < T(10.))
    {
        std::cout << "NEEDED TO CORRECT THE INLET PRESSURE\n";
        inletstate[EInletPressure] += T(1000.);
        difference = (inletstate[EInletPressure] - rightstate[TPBrCellConservation::EPressureWater]);
    }

	fluxresidual[EInletMassFlux] = fInletMassFlux - interfacestate[EMassFluxSteam] - interfacestate[EMassFluxWater];
	fluxresidual[EInletEnergyFlux] = fInletEnergyFlux - interfacestate[EEnergyFlux];
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "inlet mass flux " << fInletMassFlux << std::endl;
        sout << "residual of the inlet mass flux " << fluxresidual[EInletMassFlux] << std::endl;
        sout << "residual of the inlet energy " << fluxresidual[EInletEnergyFlux] << std::endl;
        sout << "fInletMassFlux " << fInletMassFlux << std::endl;
        sout << "interfacestate[EMassFluxSteam] " << interfacestate[EMassFluxSteam] << std::endl;
        sout << "interfacestate[EMassFluxWater] " << interfacestate[EMassFluxWater] << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    T press = inletstate[EInletPressure]+T(TPBrScales::fReferencePressure);

    if(press < T(0.))
    {
        LOGPZ_DEBUG(logger,"negative pressure")
    }
//    T press = inletstate[EInletPressure];
    T saturationtemperature = TemperatureSaturation(press);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "pressure " << press << std::endl;
        sout << "temperature " << inletstate[EInletTemperature] << std::endl;
        sout << "saturation temperature " << saturationtemperature << std::endl;
        sout << "inlet saturation " << inletstate[EInletSteamSaturation] << std::endl;
        LOGPZ_DEBUG(logger2, sout.str())
    }
#endif

    std::stringstream sout;
    if (inletstate[EInletTemperature] > saturationtemperature) 
    {
        sout << "The inlet temperature is larger than the saturation temperature at line " << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletTemperature]-saturationtemperature;
//        inletstate[EInletSteamSaturation] += T(0.01);
        //        inletstate[EInletTemperature] = saturationtemperature;
    }
    if (inletstate[EInletSteamSaturation] < T(0.))
    {
        sout << "The steam saturation is less than zero at line " << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletSteamSaturation];
//        inletstate[EInletSteamSaturation] = 0.;
    }
    if (inletstate[EInletSteamSaturation] > T(0.))
    {
        sout << "The inlet steam saturation is larger than zero at line " << __LINE__ << std::endl;
        sout << "The inlet temperature " << inletstate[EInletTemperature] << std::endl;
        sout << "saturation temperature " << saturationtemperature << std::endl;
        fluxresidual[EInletStateEqs] = inletstate[EInletTemperature]-saturationtemperature;
    }
    else if(inletstate[EInletTemperature] < saturationtemperature)
    {
        sout << "The temperature is lower than the saturation temperature at line " << __LINE__;
        fluxresidual[EInletStateEqs] = inletstate[EInletSteamSaturation];
    }

#ifdef LOG4CXX
    sout << std::endl;
    sout << "flux residual inlet equation " << EInletStateEqs << " " << fluxresidual[EInletStateEqs] << std::endl;
    LOGPZ_DEBUG(logger2, sout.str())
#endif
	TPZManVector<T> leftstate(TPBrCellConservation::NumCellEq);
    ComputeLeftState(inletstate,leftstate);
#ifdef LOG4CXX
	{
		std::stringstream sout;
        sout << "Transferring the inletstate to a leftstate\n";
        sout << "\nsuminterfacestate " << suminterfacestate;
        sout << "\ninletstate " << inletstate;
		sout << "\nleftstate " << leftstate;
        sout << "\nrightstate " << rightstate;
        sout << "\ninterfacestate " << interfacestate;
        sout << "\nMass flux water " << interfacestate[EMassFluxWater];
        sout << "\nMass flux steam " << interfacestate[EMassFluxSteam];
		LOGPZ_DEBUG(logger2,sout.str())
	}
#endif
    TPZManVector<T,NumFluxEq> localfluxresidual(NumFluxEq);
	FluxResidual(leftstate,rightstate,interfacestate,delx,area,delt,localfluxresidual);
	int i;
	for (i=NumInletVars+NumFluxEq-1; i>= NumInletVars; i--) {
		fluxresidual[i] = localfluxresidual[i-NumInletVars];
	}
    fluxresidual[EInletMassFlux] += fluxresidual[NumInletVars+EMassFluxWaterEq];
    fluxresidual[EInletEnergyFlux] += fluxresidual[NumInletVars+EEnergyFluxEq];


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
	rightstate[TPBrCellConservation::EPressureWater] = T(fFarfieldPressureWater-TPBrScales::fReferencePressure);
	rightstate[TPBrCellConservation::EPressureSteam] = T(fFarfieldPressureSteam-TPBrScales::fReferencePressure);
	rightstate[TPBrCellConservation::EPressureOil] = T(fFarfieldPressureOil-TPBrScales::fReferencePressure);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "farfield (relative) pressure water" << fFarfieldPressureWater-TPBrScales::fReferencePressure << std::endl;
        sout << "farfield temperature " << fFarfieldTemperature;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
//	rightstate[TPBrCellConservation::EPressureWater] = T(fFarfieldPressureWater);
//	rightstate[TPBrCellConservation::EPressureSteam] = T(fFarfieldPressureSteam);
//	rightstate[TPBrCellConservation::EPressureOil] = T(fFarfieldPressureOil);
	rightstate[TPBrCellConservation::EPhaseChange] = T(0.);
	FluxResidual(leftstate,rightstate,interfacestate,delx,area,delt,fluxresidual);
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

/*
 enum ESteamFluxEq {
 EMassFluxWaterEq, EMassFluxSteamEq, EMassFluxOilEq, EDarcyVelocityWaterEq, EDarcyVelocitySteamEq, EDarcyVelocityOilEq, EEnergyFluxEq
 };
 enum ESteamFluxVars {
 EMassFluxWater, EMassFluxSteam, EMassFluxOil, EDarcyVelocityWater, EDarcyVelocitySteam, EDarcyVelocityOil, EEnergyFlux
 };
 
 enum EInletEq { EInletMassFlux, EInletEnergyFlux, EInletStateEqs };
 
 enum EInletVars { EInletPressure, EInletSteamSaturation, EInletTemperature };
 

*/
/// associate scale factors with the equations and state variables
void TPBrSteamFlux::Scales(TPZVec<REAL> &eqscales, TPZVec<REAL> &statescales)
{
    eqscales.Resize(NumFluxEq);
    statescales.Resize(NumFluxEq);
    eqscales[EMassFluxWaterEq] = 1.;
    eqscales[EMassFluxSteamEq] = 1.;
    eqscales[EMassFluxOilEq] = 1.;
    eqscales[EDarcyVelocityWaterEq] = 1.;
    eqscales[EDarcyVelocitySteamEq] = 1.;
    eqscales[EDarcyVelocityOilEq] = 1.;
    eqscales[EEnergyFluxEq] = 1.;
    statescales[EMassFluxWater] = 1.;
    statescales[EMassFluxSteam] = 1.;
    statescales[EMassFluxOil] = 1.;
    statescales[EDarcyVelocityWater] = 1.;
    statescales[EDarcyVelocitySteam] = 1.;
    statescales[EDarcyVelocityOil] = 1.;
    statescales[EEnergyFlux] = 1.;
}

/// associate scale factors with the equations and state variables
void TPBrSteamFlux::InletScales(TPZVec<REAL> &eqscales, TPZVec<REAL> &statescales)
{
    eqscales.Resize(NumInletVars+NumFluxEq);
    statescales.Resize(NumInletVars+NumFluxEq);
    eqscales[EInletMassFlux] = 1.;
    eqscales[EInletEnergyFlux] = 1.;
    eqscales[EInletStateEqs] = 1.;
    statescales[EInletPressure] = TPBrScales::fReferencePressure;
    statescales[EInletSteamSaturation] = 1.;
    statescales[EInletTemperature] = 1.;
    TPZManVector<REAL> loceqscale,locstatescale;
    Scales(loceqscale, locstatescale);
    int ist;
    int nst = loceqscale.NElements();
    for (ist=0; ist<nst; ist++) {
        eqscales[NumInletVars+ist] = loceqscale[ist];
        statescales[NumInletVars+ist] = locstatescale[ist];
    }
}

template<class T>
T TPBrSteamFlux::TemperatureSaturation(T pressuresteam)//[K]
{
    T pressadapt = pressuresteam/(T(1000.0));
    T temp_de_saturac = waterdata.getSaturationStateTemperature(pressadapt);
    /*
     T  val_log, temp;
     T temp_de_saturac;
     val_log = log(pressuresteam*0.0001450377438972831);
     temp=561.435 + 33.8866*val_log + 2.18893*(val_log*val_log) + 0.0808998*(val_log*val_log*val_log) +
     0.0342030*(val_log*val_log*val_log*val_log);
     temp_de_saturac = (temp - 32. - 459.67)/1.8;
     */
	return temp_de_saturac;
}



template<class T>
T TPBrSteamFlux::EnthalpyWater(T temperature)//[kJ/kg]
{
	return waterdata.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
}

template<>
REAL TPBrSteamFlux::EnthalpyWater(REAL temperature)//[kJ/kg]
{
	return waterdata.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
}

/// Compute a limit for correcting the solution
REAL TPBrSteamFlux::LimitRangeInlet(REAL scale,TPZVec<REAL> &inletstate,TPZVec<REAL> &cellstate, TPZVec<REAL> &inletcorrection, TPZVec<REAL> &cellcorrection)
{
    REAL delp = inletstate[EInletPressure]-cellstate[TPBrCellConservation::EPressureWater];
    REAL corrdelp = inletcorrection[EInletPressure]-cellcorrection[TPBrCellConservation::EPressureWater];
    if(delp < 0.)
    {
        DebugStop();
    }
    if(delp+scale*corrdelp < 0.)
    {
        scale = -delp/corrdelp;
    }
    REAL press = inletstate[EInletPressure]+TPBrScales::fReferencePressure;
    REAL tempsat = TemperatureSaturation(press);
    if(inletstate[EInletSteamSaturation] > 0. && inletstate[EInletTemperature] > tempsat*1.1)
    {
        std::stringstream sout;
        REAL sat = inletstate[EInletSteamSaturation];
        REAL temp = inletstate[EInletTemperature];
        sout << "Saturation " << sat << " temperature " << temp << " sat temperature " << tempsat;
        LOGPZ_ERROR(logger, sout.str())
        std::cout << sout.str() << std::endl;
    }
    if(IsZero(inletstate[EInletSteamSaturation]) && inletstate[EInletTemperature] < tempsat 
        && inletstate[EInletTemperature]+scale*inletcorrection[EInletTemperature] > 1.1*tempsat)
    {
        scale = (1.1*tempsat-inletstate[EInletTemperature])/inletcorrection[EInletTemperature];
    }
    return scale;
}

/// Compute a limit for correcting the solution
REAL TPBrSteamFlux::LimitRange(REAL scale,TPZVec<REAL> &interfacestate,TPZVec<REAL> &interfacecorrection)
{
    return scale;
}



#endif
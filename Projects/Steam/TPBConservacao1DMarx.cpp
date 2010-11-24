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

#include <math.h>

void ScaleFactor(TPZFMatrix &tangentmatrix, TPZFMatrix &residualmatrix, TPZManVector<REAL> &scalevalues, TPZManVector<REAL> &statescalevalues );
void ScaleFactorSol(TPZFMatrix &residualmatrix, TPZManVector<REAL> &scalevalues );
int main()
{
	REAL PI = 4*atan(1.);	
	TPBrCellMarx first;
	
	//WaterDataInStateOfSaturation aux;
//	double temp =0.;
//	double pres = 2000.;
//	temp= aux.getSaturationStateTemperature(pres);
//	pres = aux.getSaturationStatePressure(temp);
	
	//----------------- dados de entrada ------------------------------
	//dados numerico
	REAL TimeStep = 0.01;
	
	//dados da celula
	REAL rint = 0.15;
	REAL rext = 1.5;
	REAL CellSize = 1.;
	REAL LeftArea = PI*rint*rint;
	REAL RightArea = PI*rext*rext;
	REAL CellVolume = PI*(rext-rint)*(rext-rint)*CellSize;// 5.725552611167399 ;
	
	//dados da rocha
	REAL MaterialPermeability = 0.8e-12;
	PhysicalProperties minharocha(0,1);
	//dados da injecao
	REAL PressureWater(2.e6);
	REAL tempReservtorio = 98.0;
	REAL Temperature(tempReservtorio);
	TPZManVector<REAL> InitialSaturation(3,0.);
	InitialSaturation[TPBrCellMarx::EOil] = 0.0170871;
	InitialSaturation[TPBrCellMarx::EWater] = 1.-0.0170871;
	InitialSaturation[TPBrCellMarx::ESteam] = 0. ;
	TPZManVector<REAL> Massflux(3,0.);
	Massflux[TPBrCellMarx::EOil] = 0.0;
	Massflux[TPBrCellMarx::EWater] = 0.121951;
	Massflux[TPBrCellMarx::ESteam] = 0.555556 ;
	//-----------------------------------------------------------------------------------------
	
	TPZManVector<REAL> initial(TPBrCellMarx::NUMVARS,0.), residual(TPBrCellMarx::NUMVARS,0.);
	TPZManVector<REAL> leftstate(TPBrCellMarx::NUMVARS,0.),rightstate(TPBrCellMarx::NUMVARS,0.);
	
	first.SetMaterialProperty(MaterialPermeability, TimeStep,minharocha);
	first.SetGeometry(CellVolume,LeftArea,RightArea,CellSize);
	first.SetCellState(PressureWater,InitialSaturation,Temperature);
	first.SetInjectionState(PressureWater, Massflux, leftstate);
	
	first.InitializeState(initial);
		
	first.TotalResidual(leftstate,initial,residual);
			
	TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > tangent(TPBrCellMarx::NUMVARS,0.), state(TPBrCellMarx::NUMVARS,0.);
	first.InitializeState(state);
	first.TotalResidual(leftstate, state, tangent);
	
	
	//-------------------------- FATORES DE ESCALA  ---------------------------------------------------
	TPZManVector <REAL> vecState(initial), scalevalues(18,0.), statescalevalues(18,0.);
	TPZFMatrix scalevaluesmatrix,statescalevaluesmatrix;
	vecState[TPBrCellMarx::EMassFluxOil] = leftstate[TPBrCellMarx::EMassFluxOil];
	vecState[TPBrCellMarx::EMassFluxWater] = leftstate[TPBrCellMarx::EMassFluxWater];
	vecState[TPBrCellMarx::EMassFluxSteam] = leftstate[TPBrCellMarx::EMassFluxSteam];
	first.ReferenceResidualValues(leftstate, scalevalues);
	first.ReferenceStateValues(vecState,statescalevalues);
	
	first.ExtractMatrix(scalevalues,scalevaluesmatrix);
	first.ExtractMatrix(statescalevalues,statescalevaluesmatrix);
	scalevaluesmatrix.Print("scalevalues = ",cout, EMathematicaInput);
	statescalevaluesmatrix.Print("statescalevalues = ",cout, EMathematicaInput);
	
	//----------------------------- corrigir residuo e tangente -----------------------------
	TPZFMatrix tangentmatrix,residualmatrix,statematrix;
	first.ExtractMatrix(tangent,tangentmatrix);
	first.ExtractMatrix(residual,residualmatrix);
	first.ExtractMatrix(initial,statematrix);
	
	ScaleFactor(tangentmatrix, residualmatrix, scalevalues, statescalevalues );
	//---------------------------------------------------------------------------------------------	
			
	tangentmatrix.Print("tangentmatrix = ",cout, EMathematicaInput);
	residualmatrix.Print("Residualmatrix = ",cout, EMathematicaInput);
	
	cout << "\n ======= Metodo de Newton =======\n\n";
	while (Norm(residualmatrix) > 1.e-6) {
		tangentmatrix.SolveDirect(residualmatrix, ELU);
		residualmatrix.Print("ResidualMatrix",cout);
		// multiplicar o valor pelo scalestate
		ScaleFactorSol(residualmatrix, statescalevalues);
		
		statematrix += residualmatrix;
		statematrix.Print("statematrix", cout);
		first.ConvertState(statematrix,rightstate);
		first.ConvertState(statematrix,state);
		first.TotalResidual(leftstate, rightstate, residual);
		first.TotalResidual(leftstate, state, tangent);
		tangentmatrix.Zero();
		first.ExtractMatrix(tangent,tangentmatrix);
		
		// aplicar os fatores de escala novamente no tangentmatrix e residual
		ScaleFactor(tangentmatrix, residualmatrix, scalevalues, statescalevalues );
	}
	
	
	tangentmatrix.Print("tangentmatrix = ",cout, EMathematicaInput);
	residualmatrix.Print("Residualmatrix = ",cout, EMathematicaInput);
	
	return 0;
};

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
	temp_de_saturac = (temp-32 - 459.67)/1.8;
	
	return temp_de_saturac;
}


// metodos para recuperar os dados tabulados em funcao
template<class T>
//[kg/m3]
T TPBrCellMarx::DensityOil(T temp) {
	OilData oil;
	return oil.getDensityToOil(temp);
}
template<class T>
//[kg/m3]
T TPBrCellMarx::DensityWater(T temp) {
	WaterDataInStateOfSaturation t;
	return t.getSaturationStateDensityToLiquidWater(temp);
}
template<class T>
//[kg/m3]
T TPBrCellMarx::DensitySteam(T temp) {
	WaterDataInStateOfSaturation t;
	return t.getSaturationStateDensityToSteam(temp);
}

//[Pa*sec]
REAL TPBrCellMarx::ViscosityOil(REAL temp) {
	OilData oil;
	return oil.getDynamicViscosityToOil(temp);
}

//[Pa*sec]
REAL TPBrCellMarx::ViscosityWater(REAL temp) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateViscosityToLiquidWater(temp);
}

//[Pa*sec]
REAL TPBrCellMarx::ViscositySteam(REAL temp) {
	WaterDataInStateOfSaturation steam;
	return steam.getSaturationStateViscosityToSteam(temp);
}

template<class T>
//dimensionless
void TPBrCellMarx::ComputeRelativePermeability(TPZManVector<T> &saturation,TPZManVector<T> &relativepermeability) {
	relativepermeability = saturation;
}
template<class T>
T TPBrCellMarx::TemperatureSaturation(T p) {
	WaterDataInStateOfSaturation water;
	T p1000 = p/1000.;
	return water.getSaturationStateTemperature(p1000);
}

template<class T>
//[kJoule/kg]
T TPBrCellMarx::EnthalpyWater(T temperature) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateSpecificEnthalpyToLiquidWater(temperature);
}
template<class T>
//[kJoule/kg]
T TPBrCellMarx::EnthalpySteam(T temperature) {
	WaterDataInStateOfSaturation water;
	return water.getSaturationStateSpecificEnthalpyToSteam(temperature);
}
template<class T>
//[kJoule/kg]
T TPBrCellMarx::EnthalpyOil(T temperature) {
	OilData oil;
#warning Verificar se specific heat e enthalply??
	return  oil.getSpecificHeatToOil(temperature);
}


void ScaleFactor(TPZFMatrix &tangentmatrix, TPZFMatrix &residualmatrix, TPZManVector<REAL> &scalevalues, TPZManVector<REAL> &statescalevalues ){
	
	int numvar = TPBrCellMarx ::NUMVARS;
	for (int ir =0; ir<numvar; ir++) {
		residualmatrix(ir,0) = residualmatrix(ir,0)/scalevalues[ir] ;
	}
	
	TPZFMatrix  tangentmatrixtr;
	for (int it = 0; it<numvar; it++) {
		for (int jt = 0; jt<numvar; jt++) {
			tangentmatrix(it,jt) = tangentmatrix(it,jt)/scalevalues[it];  
		}
	}
	tangentmatrix.Transpose(&tangentmatrixtr);
	for (int itr = 0; itr<numvar; itr++) {
		for (int jtr = 0; jtr<numvar; jtr++) {
			tangentmatrixtr(itr,jtr) = tangentmatrixtr(itr,jtr)*statescalevalues[itr];  
		}
	}
	tangentmatrixtr.Transpose(&tangentmatrix);
}

void ScaleFactorSol(TPZFMatrix &residualmatrix, TPZManVector<REAL> &scalevalues ){
	
	int numvar = TPBrCellMarx ::NUMVARS;
	for (int ir =0; ir<numvar; ir++) {
		residualmatrix(ir,0) = residualmatrix(ir,0)*scalevalues[ir] ;
	}	
}
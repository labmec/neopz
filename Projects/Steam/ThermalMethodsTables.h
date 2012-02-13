/*
 *  ThermalMethodsTables.h
 *  FrenteVapor
 *
 *  Created by Jorge Calle on 4/2/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ThermalMethodsTablesToSteamHHH
#define ThermalMethodsTablesToSteamHHH

#include "pzerror.h"
#include "pzreal.h"

#ifdef _AUTODIFF
using namespace std;
#include "tfad.h"
#endif

#include <vector>


// Funções auxiliares para interpolar o valor para var (variável) desde os valores tabulados em Table.
template<class T>
T getInterpolatedValue(std::vector<std::vector<double> > &Table,T var,int index);

// Devolve o valor da variável interpolada na coluna de variáveis independentes a partir do valor
// de uma variável dependente dada na coluna indexproperty
//double getIndexerFromValuesInterpolating(std::vector<std::vector<double> > &fTable,int indexproperty,double value);

// Classe que armazena a tabela de valores das propriedades físicas para algum material
// Pode ser indicado o indice da variavel independente, assumindo então todas as demais dependentes dela
// sugere-se utilizar: O índice 1 (coluna 1) para armazenar a variável da qual as propriedades físicas das outras colunas são dependentes
// O número de linhas é o total de entradas tabuladas
// O número de colunas-1 é o total de propriedades físicas tabuladas e dependentes do valor armazenado na coluna 1
class DataInTable {
	
protected:
	std::vector<std::vector<double> > fTable;
	
	// Indice da variavel independente
	int fIndependentVarIndex;
	
	// Método para interpolar o valor Temp, entre os valores tabulados no indice (i-1), Value(i-1) e no indice (i), Value(i)
	template<class T>
	T getInterpolatedValue(T temp,int indexproperty);
	
	// Função que a partir dos pares (initialvar, initialvalue) e (finalvar, finalvalue) interpola o valor da variável dependente sobre o valor da independente - var.
	template<class T>
	T getInterpolatedValue(double initialvar,double intialvalue,double finalvar,double finalvalue,T var);
public:
	// Construtor com: ndata -> numero de dados associados (linhas) na tabela, nitens -> numero de valores relacionados para cada linha (colunas)
	// Primeira e segunda coluna correspondem a Temperatura e Pressão, os demais dependem do tipo de dado armazenado. Assim  nitens > 2.
	DataInTable(int ndata,int nitens);
	
	// Método que para uma dada temperatura obtemos o indice da linha na tabela (a temperatura é armazenada em intervalos iguais)
	template<class T>
	int getIndex(T temp);
	
	// Métodos que dada uma temperatura obtem o indice (linha) e pega o correspondente valor do item
	//double getValue(int ndata,int nitem);
	//double getValue(double indexer,int nitem);
	
	// Método que devolve o valor da variável independente para o índice fornecido
	//double getIndependentValue(int index);
	//double getDependentValue(int index,int indexproperty);
};
//-------------------------------------------------------------------------------------------------------------------

// Classe que construi uma tabela com as propriedades físicas da agua liquida e gasosa (vapor) em estado de saturação:
// Na tabela, o índice 1 (coluna 1) contem as temperaturas do estado de saturação [ C - graus Celsius ]
// e o índice 2 (coluna 2) contem a pressão associada para o estado de saturação  [ kPa - Kilo Pascal ]
// O índice 3 (coluna 3) armazena o calor latente para o correspondente estado de saturação. [ kJ/kg ]
// O calor latente para temperaturas superiores a 374 C é nulo, assim a água passa do estado líquido para o vapor diretamente.
// Os índices do 4 ao 8 armazenam propriedades físicas da água no estado aquoso (líquido)
// Os índices do 9 ao 13 armazenam propriedades físicas da água no estado gasoso (vapor)
// O índice 4(líquido) e 9(vapor) armazenam as densidades correspondentes [ kg/(m3) ]
// O índice 5(líquido) e 10(vapor) armazenam a viscosidade dinâmica [ Pa * sec ]
// O índice 6(líquido) e 11(vapor) armazenam o calor específico [ kJ/(kg * K) ]
// O índice 7(líquido) e 12(vapor) armazenam a condutividade térmica [ kJ/(sec*m*K) ]
// O índice 8(líquido) e 13(vapor) armazenam a entalphia específica [ kJ/kg ]
class  WaterDataInStateOfSaturation : public DataInTable {
	
public:
	WaterDataInStateOfSaturation();
	
	// Método que devolve a temperatura associada a pressão fornecida para um correspondente estado de saturação
	template<class T>
	T getSaturationStateTemperature(T pressure);//[ C - graus Celsius pressure kPa]
	
	// Método que devolve a pressão e o calor latente associados a temperatura (temp) fornecida para um correspondente estado de saturação
	template<class T>
	T getSaturationStatePressure(T temp);//[ kPa - Kilo Pascal ]
	
	template<class T>
	T getSaturationStateLatentHeat(T temp);//[ kJ/kg ]
	
	// Métodos que devolvem as propriedades físicas: densidade, viscosidade, calor específico, condutividade térmica e entalphia específica 
	// a partir dos dados tabulados - para a água líquida em estado de saturação
	template<class T>
	T getSaturationStateDensityToLiquidWater(T temp);//[ kg/(m3) ]
	template<class T>
	T getSaturationStateViscosityToLiquidWater(T temp);//[ Pa * sec ]
	template<class T>
	T getSaturationStateSpecificHeatToLiquidWater(T temp);//[ kJ/(kg*Kelvin) ]
	template<class T>
	T getSaturationStateThermalConductivityToLiquidWater(T temp);//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
	template<class T>
	T getSaturationStateSpecificEnthalpyToLiquidWater(T temp);//[ kJ/kg ]
	
	// Métodos que devolvem as propriedades físicas: densidade, viscosidade, calor específico, condutividade térmica e entalphia específica 
	// a partir dos dados tabulados - para a água gasosa (vapor) em estado de saturação
	template<class T>
	T getSaturationStateDensityToSteam(T temp);//[ kg/(m3) ]
	template<class T>
	T getSaturationStateViscosityToSteam(T temp);//[ Pa * sec ]
	template<class T>
	T getSaturationStateSpecificHeatToSteam(T temp);//[ kJ/(kg*Kelvin) ]
	template<class T>
	T getSaturationStateThermalConductivityToSteam(T temp);//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
	template<class T>
	T getSaturationStateSpecificEnthalpyToSteam(T temp);//[ kJ/kg ]
	
	// Métodos que devolvem as propriedades físicas: densidade, viscosidade, calor específico, condutividade térmica e entalphia específica 
	// a partir de uma fórmula - para a água líquida (aquosa)
	template<class T>
	T getDensityToLiquidWater(T temp);//[ kg/(m3) ]
	template<class T>
	T getViscosityToLiquidWater(T temp);//[ Pa * sec ]
	template<class T>
	T getSpecificHeatToLiquidWater(T temp);//[ kJ/(kg*Kelvin) ]
	template<class T>
	T getThermalConductivityToLiquidWater(T temp);//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
	template<class T>
	T getSpecificEnthalpyToLiquidWater(T temp);//[ kJ/kg ]
	
	// Métodos que devolvem as propriedades físicas: densidade, viscosidade, calor específico, condutividade térmica e entalphia específica 
	// a partir de uma fórmula - para o vapor (gasoso)
	template<class T>
	T getDensityToSteam(T temp);//[ kg/(m3) ]
	template<class T>
	T getViscosityToSteam(T temp);//[ Pa * sec ]
	template<class T>
	T getSpecificHeatToSteam(T temp);//[ kJ/(kg*Kelvin) ]
	template<class T>
	T getThermalConductivityToSteam(T temp);//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
	template<class T>
	T getSpecificEnthalpyToSteam(T temp);//[ kJ/kg ]
	
	// Obtendo valores das saturações da água ou vapor dependendo da qualidade (titularidade) do vapor
	template<class T>
	T getSaturationWater(int i,T quality,T So);
};
//-------------------------------------------------------------------------------------------------------------------------

// Physical Properties tabuled for some crude oils
class OilData : public DataInTable {
public:
	OilData();
	
	// Métodos que retornam as propriedades físicas do óleo
	template<class T>
	T getDensityToOil(T temp);//[ kg/(m3) ]
	template<class T>
	T getSpecificGravityToOil(T temp);//dimensionless (kg/m3)/(kg/m3): (density of the material)/(density of water) at a specified temperature
	
	template<class T>
	T getThermalConductivityToOil(T temp);//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
	template<class T>
	T getSpecificHeatToOil(T temp);//[ kJ/(kg*Kelvin) ]
	template<class T>
	T getDynamicViscosityToOil(T temp);//[ Pa * sec ]
};
//------------------------------------------------------------------------------------------------------------

// Métodos que retornam as propriedades físicas do óleo
template<class T>
inline T OilData::getDensityToOil(T temp) {
	// Estado inicial conhecido
	REAL rho_ini = fTable[4][2];
//	REAL T_ini = fTable[4][0], p_ini = fTable[4][1], ThermalExpansion_ini = 0.001, Compressibility_ini = 0.1;
	//	double rho = rho_ini*(1.-(ThermalExpansion_ini*(temp-T_ini))+(Compressibility_ini*(getPression(temp)-p_ini)));
	//	return rho;
	return rho_ini;
}

template<class T>
//dimensionless
inline T OilData::getSpecificGravityToOil(T temp) {
	return 0.0;
}

template<class T>
//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
inline T OilData::getThermalConductivityToOil(T temp) {
	return 0.0;
}

template<class T>
//[ kJ/(kg*Kelvin) ]
inline T OilData::getSpecificHeatToOil(T temp) {
	return fTable[4][4];
}

template<class T>
//[ Pa * sec ]
inline T OilData::getDynamicViscosityToOil(T temp) {
	return 1.e-5;
}

// Apenas pelo fato de a temperatura e a pressão serem relacionadas para um estado de saturação, assume-se
// que a temperatura e a pressão foram tabulados em forma crescente
// O fato de ser armazenado em forma crescente ou decrescente pode ser definido olhando as primeiras entradas da tabela
template<class T>
//[ C - graus Celsius ]
//[ pressure kpa]
inline T WaterDataInStateOfSaturation::getSaturationStateTemperature(T pressure) {
	fIndependentVarIndex = 1;
	T temp = getInterpolatedValue(pressure,0);
	fIndependentVarIndex = 0;
	return temp;
}

template<class T>
//[ kPa - Kilo Pascal ]
inline T WaterDataInStateOfSaturation::getSaturationStatePressure(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,1);
}

template<class T>
//[ kJ/kg ]
inline T WaterDataInStateOfSaturation::getSaturationStateLatentHeat(T temp) {
	fIndependentVarIndex = 0;
	if(temp > 374.15) return 0.0;
	return getInterpolatedValue(temp,2);
}

// Physical properties to liquid water (fase aquosa) from tabulated data
template<class T>
//[kg/m3]
inline T WaterDataInStateOfSaturation::getSaturationStateDensityToLiquidWater(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,3);
}

template<class T>
//[Pa*sec]
inline T WaterDataInStateOfSaturation::getSaturationStateViscosityToLiquidWater(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,4);
}

template<class T>
//[ kJ/(kg*Kelvin) ]
inline T WaterDataInStateOfSaturation::getSaturationStateSpecificHeatToLiquidWater(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,5);
}

template<class T>
//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
inline T WaterDataInStateOfSaturation::getSaturationStateThermalConductivityToLiquidWater(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,6);
}

template<class T>
//[kJoule/kg]
inline T WaterDataInStateOfSaturation::getSaturationStateSpecificEnthalpyToLiquidWater(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,7);
}

// Physical properties to steam (fase gasosa) from tabulated data
template<class T>
//[kg/m3]
inline T WaterDataInStateOfSaturation::getSaturationStateDensityToSteam(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,8);
}

template<class T>
//[Pa*sec]
inline T WaterDataInStateOfSaturation::getSaturationStateViscosityToSteam(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,9);
}

template<class T>
//[ kJ/(kg*Kelvin) ]
inline T WaterDataInStateOfSaturation::getSaturationStateSpecificHeatToSteam(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,10);
}

template<class T>
//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
inline T WaterDataInStateOfSaturation::getSaturationStateThermalConductivityToSteam(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,11);
}

template<class T>
//[kJoule/kg]
inline T WaterDataInStateOfSaturation::getSaturationStateSpecificEnthalpyToSteam(T temp) {
	fIndependentVarIndex = 0;
	return getInterpolatedValue(temp,12);
}

// Physical properties to steam (fase aquosa) from formula
template<class T>
//[kg/m3]
inline T WaterDataInStateOfSaturation::getDensityToLiquidWater(T temp) {
	if(temp > 0.0 && temp < 330)
		return (997 - (.046*temp) - (0.00306*temp*temp));   // Formula com aprox. 1% Burger - pag 47
	// Utilizar dados a partir de um conhecido relacionamento de densidade, pressao temperatura e os coeficientes
	// PENDENTE
	return getSaturationStateDensityToLiquidWater(temp);
}

template<class T>
//[Pa*sec]
inline T WaterDataInStateOfSaturation::getViscosityToLiquidWater(T temp) {
	if(temp > 0.0 && temp < 310) {
		T visc = temp - 8.435;
		T value = sqrt(visc*visc + 8078.4);
		value += visc;
		value *= 0.021482;
		return (1./(value - 1.2));
	}
	return getSaturationStateViscosityToLiquidWater(temp);
}

template<class T>
//[ kJ/(kg*Kelvin) ]
inline T WaterDataInStateOfSaturation::getSpecificHeatToLiquidWater(T temp) {
	return getSaturationStateSpecificHeatToLiquidWater(temp);
}

template<class T>
//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
inline T WaterDataInStateOfSaturation::getThermalConductivityToLiquidWater(T temp) {
	if(temp > 0.0 && temp < 330)
		return (0.56 + 0.00185*temp - (0.00000634*temp*temp));   // Burger 51 - Para agua liquida o maximo valor é para temp = 130 C
	return getSaturationStateThermalConductivityToLiquidWater(temp);
}

template<class T>
//[kJoule/kg]
inline T WaterDataInStateOfSaturation::getSpecificEnthalpyToLiquidWater(T temp) {
	if((temp - T(0.01 - ZeroTolerance())) < 0. && (-temp + T(0.01) - T(ZeroTolerance())) < 0.)
		return T(0.0);
	else if(temp < T(130.))
		return (T(4.2)*temp);
	else if(temp < T(310))
		return (T(4.2)*temp + (T(0.003)*(temp-T(130.))*(temp-T(130))));
	else if(temp < T(373.))
		return (T(4.2)*temp + (T(0.003)*(temp-T(130.))*(temp-T(130.))) + (0.0008*(temp-T(310.))*(temp-T(310.))*(temp-T(310.))));
	return getSaturationStateSpecificEnthalpyToLiquidWater(temp);
}

// Physical properties to steam (fase gasosa) from formula
template<class T>
//[kg/m3]
inline T WaterDataInStateOfSaturation::getDensityToSteam(T temp) {
	double Z = 1.012 - 0.0004461*temp + (0.00000298*temp*temp) - (0.00000001663*temp*temp*temp);
	if(!IsZero(temp) && !IsZero(Z)) {
		return ((18.015/(8.3143*Z))*(getSaturationStatePression(temp)/temp));
	}
	return getSaturationStateDensityToSteam(temp);
}

template<class T>
//[Pa*sec]
inline T WaterDataInStateOfSaturation::getViscosityToSteam(T temp) {
	if(temp > 0.0 && temp < 330)
		return (88. + 0.36*temp);
	else if(temp > 0.0 && temp < 374.)
		return getInterpolatedValue(3,temp);
	else if(temp > 0.0 && temp < 700)
		return (88. + 0.38*temp);
	return getSaturationStateViscosityToSteam(temp);
}

template<class T>
//[ kJ/(kg*Kelvin) ]
inline T WaterDataInStateOfSaturation::getSpecificHeatToSteam(T temp) {
	return getSaturationStateSpecificHeatToSteam(temp);
}

template<class T>
//[ (kJ/sec)/(m*Kelvin) ] = [Watt/m*Kelvin]
inline T WaterDataInStateOfSaturation::getThermalConductivityToSteam(T temp) {
	return getSaturationStateThermalConductivityToSteam(temp);
}

template<class T>
//[kJoule/kg]
inline T WaterDataInStateOfSaturation::getSpecificEnthalpyToSteam(T temp) {
	if(temp > 0.0 && temp < 320)
		return (2500.+(1.88*temp)-(0.0000037*pow(temp,3.2)));
	return getSaturationStateSpecificEnthalpyToSteam(temp);
}

template<class T>
//dimensionless
inline T WaterDataInStateOfSaturation::getSaturationWater(int i,T quality,T So) {
	T Sw = 1. - (quality+So);
	if(Sw < 0.0)
		return 0.0;
	if(Sw > 1.0)
		return 1.0;
	return Sw;
	
}

//-------------------------------------------------------------------------------------------------------------------------------------
template<class T>
inline T getInterpolatedValue(std::vector<std::vector<double> > &fTable,T var,int index) {
	// Determinando primeiro em qual intervalo de dados da tabela está o valor no qual interpolaremos - value
	// os valores referencia para a tabela (dado inicial - ex. temperatura) devem estar dados em intervalo constante
	int last = fTable.size();
	if(last < 2 || index < 0) DebugStop();
	// Se index é zero será retornado o mesmo valor, pois está se refirindo a própria variável independente
	if(!index) return var;
	
	// Verificando se a tabela tem o número suficiente de colunas para devolver o valor interpolado na propriedade index
	int indexes = fTable[0].size(), i;
	if(indexes < 2 || index >= indexes)
		DebugStop();
	// Verificando se var está no intervalo fTable[0][0] até fTable[0][last-1]
	if(var < fTable[0][0])
		return fTable[0][index];
	else if(var > fTable[last-1][0])
		return fTable[last-1][index];
	
	// Procuramos o intervalo entre os quais está o valor var na coluna da propriedade independente
	for(i=1;i<last;i++) {
		if(IsZero(var - fTable[i-1][0]))
			return fTable[i-1][index];
		if(var < fTable[i][0])
			break;
	}
	if(i==last)
		return fTable[last-1][index];
	
	// i e (i-1) indicam as linhas correspondente entre os quais deveremos interpolar a variável var
	double dx = (fTable[i][0] - fTable[i-1][0]);
	double dy = (fTable[i][index] - fTable[i-1][index]);
	
	return (fTable[i-1][index] + (var - fTable[i-1][0])*(dy/dx));
}

template<class T>
inline T DataInTable::getInterpolatedValue(double initialvar,double initialvalue,double finalvar,double finalvalue,T var) {
    T zero = T(0.0);
	if(T(initialvar) > var || T(finalvar) < var)
		return zero;
	return (T(initialvalue) + ((var-T(initialvar))*T((finalvalue-initialvalue)/(finalvar-initialvar))));
}

// Devolve o valor interpolado da propriedade i entre as linhas correspondentes da variável independente que contem o valor temp no meio
template<class T>
inline T DataInTable::getInterpolatedValue(T temp,int indexproperty) {
	int index = getIndex(temp);
	if(index < 0 || !fTable.size() || indexproperty < 0 || (indexproperty+1) > fTable[0].size()) return T(0.0);
	// aqui ja temos valores validos de index e de indexproperty
	if(!index) return T(fTable[0][indexproperty]);
	T value = getInterpolatedValue(fTable[index-1][fIndependentVarIndex],fTable[index-1][indexproperty],fTable[index][fIndependentVarIndex],fTable[index][indexproperty],temp);
	return value;
}

// Retorna o indice da linha na tabela para uma temperatura temp, considerando que foram armazenadas as temperaturas a intervalos de 10 em 10 graus Celsius (Centigrados)
template<class T>
inline int DataInTable::getIndex(T temp) {
	int indice = -1;
	int last = fTable.size();
	if(last < 2 || temp < T(fTable[0][fIndependentVarIndex]) || temp > T(fTable[last-1][fIndependentVarIndex]))
	{
		DebugStop();
	}
	// Procurando a linha i até a qual superamos o valor de temp
	for(indice=1;indice<last;indice++)
		if(temp < T(fTable[indice][fIndependentVarIndex]))
			return indice;
	return -1;
}




#endif 
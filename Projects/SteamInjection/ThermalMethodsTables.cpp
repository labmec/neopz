/*
 *  ThermalMethodsTables.cpp
 *  FrenteVapor
 *
 *  Created by Jorge Calle on 4/2/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "ThermalMethodsTables.h"
#include "ConvertionMethods.h"
#include <math.h>

double DataInTable::getInterpolatedValue(double initialvar,double initialvalue,double finalvar,double finalvalue,double var)
{
	if(initialvar > var || finalvar < var)
	{
		return 0.0;
	}

	double val = (initialvalue + ((var-initialvar)*((finalvalue-initialvalue)/(finalvar-initialvar))));

	return val;
}

DataInTable::DataInTable(int ndata,int nitens)
{
	if(nitens < 2)
	{
		nitens = 2;
	}
	if(ndata < 1)
	{
		ndata = 1;
	}

	fTable.resize(ndata);

	for(int i = 0; i < ndata; i++)
	{
		fTable[i].resize(nitens);
	}
	// Primeira coluna será por default a coluna para a variavel independente
	fIndependentVarIndex = 0;
}

double DataInTable::getInterpolatedValue(double temp,int indexproperty)
{
	int index = getIndex(temp);
	if(index < 0 || !fTable.size() || indexproperty < 2 || (indexproperty+1) > fTable[0].size())
	{
		return 0.0;
	}

	if(!index)
	{
		return fTable[0][indexproperty];
	}

	double value = getInterpolatedValue(fTable[index-1][fIndependentVarIndex],fTable[index-1][indexproperty],fTable[index][fIndependentVarIndex],fTable[index][indexproperty],temp);

	return value;
}

// Retorna o indice da linha na tabela para uma temperatura temp, considerando que foram armazenadas as temperaturas a intervalos de 10 em 10 graus Celsius (Centigrados)
int DataInTable::getIndex(double temp)
{
	int indice = -1;
	int last = fTable.size();
	if(last < 2 || temp < fTable[0][fIndependentVarIndex] || temp > fTable[last-1][fIndependentVarIndex])
	{
		return indice;
	}

	// Procurando a linha i até a qual superamos o valor de temp
	for(indice=1;indice<last;indice++)
	{
		if(temp < fTable[indice][fIndependentVarIndex])
		{
			return indice;
		}
	}

	return -1;
}


// PHYSICAL PROPERTIES TO WATER

// Data of water in state of saturation. Na tabela existem os valores das propriedades físicas da água, para as fases: --> gasosa (vapor), aquosa (líquida)
WaterDataInStateOfSaturation::WaterDataInStateOfSaturation() : DataInTable(39,13)
{
	// Valores relacionados da temperatura e da pressão para o estado de saturação do fluido
	// Temperatura -> indice 0; Pressao -> indice 1; Calor Latente -> indice 7
	double Temperature[39] = { 0., 10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,300.,310., 320., 330., 340., 350., 360., 370., 374.15 };
	double Pression[39] = {0.6108, 1.227, 2.337, 4.241, 7.375, 12.335, 19.92, 31.16, 47.36, 70.11, 101.33, 143.27, 198.54, 270.13, 361.4, 476., 618.1, 792., 1002.7, 1255.1, 1554.9, 1907.7, 2319.8, 2797.6, 3347.8, 3977.6, 4694.3, 5505.8, 6420.2, 7446.1, 8592.7, 9870., 11289., 12863., 14605., 16535., 18675., 21054., 22120.};
	double LatentHeat[39] = {2501.6, 2477.9, 2454.3, 2430.7, 2406.9, 2382.9, 2358.6, 2334., 2308.8, 2283.2, 2256.9, 2230., 2202.2, 2173.6, 2144., 2113.2, 2081.3, 2047.9, 2013.1, 1976.7, 1938.6, 1898.5, 1856.2, 1811.7, 1764.6, 1714.6, 1661.5, 1604.6, 1543.6, 1477.6, 1406., 1327.6, 1241.1, 1143.6, 1030.7, 895.7, 721.3, 452.6, 0.0};
	// Caso a temperatura seja T = 372 o calor latente eh 351.4
	
	// Valores tabulados para estado de saturação - propriedades na fase gasosa (vapor)
	// Densidade -> indice 2; Viscosidade -> indice 3; Calor especifico -> indice 4; Condutividade termica -> indice 5; Entalpia especifica -> indice 6
	double Density[39] = {0.004847, .009396, .01729, .03037, .05116, .08302, .1302, .1982, .2933, .4235, .5977, .8265, 1.122, 1.497, 1.967, 2.548, 3.26, 4.123, 5.16, 6.397, 7.864, 9.593, 11.62, 14., 16.76, 19.99, 23.73, 28.1, 33.19, 39.16, 46.19, 54.54, 64.6, 76.99, 92.76, 113.6, 144., 201.1, 315.5};
	double DynamicViscosity[39] = {8.84, 9.17, 9.52, 9.86, 10.18, 10.51, 10.88, 11.28, 11.65, 12.02, 12.37, 12.71, 13.04, 13.37, 13.7, 14.02, 14.36, 14.72, 15.06, 15.4, 15.74, 16.09, 16.44, 16.81, 17.18, 17.57, 17.9, 18.28, 18.65, 19.07, 19.53, 20.1, 20.8, 21.6, 22.5, 23.7, 25.7, 29.7, 40.6};
	double SpecificHeat[39] = {1.1854, 1.86, 1.866, 1.875, 1.885, 1.899, 1.916, 1.936, 1.962, 1.992, 2.028, 2.07, 2.12, 2.176, 2.241, 2.314, 2.398, 2.491, 2.596, 2.713, 2.843, 2.988, 3.15, 3.331, 3.536, 3.772, 4.047, 4.373, 4.767, 5.253, 5.863, 6.65, 7.722, 9.361, 12.21, 17.15, 25.12, 76.92,1000.};
	double ThermalConductivity[39] = {1.67, 1.74, 1.81, 1.9, 1.97, 2.04, 2.12, 2.22, 2.31, 2.4, 2.5, 2.57, 2.68, 2.8, 2.94, 3.1, 3.25, 3.4, 3.55, 3.71, 3.88, 4.08, 4.3, 4.53, 4.8, 5.1, 5.4, 5.75, 6.15, 6.65, 7.32, 8.2, 9.2, 10.4, 11.9, 13.8, 17.4, 29.3, 91.4};
	double SpecificEnthalpy[39] = {2501.6, 2519.9, 2538.2, 2556.4, 2574.4, 2592.2, 2609.7, 2626.9, 2643.8, 2660.1, 2676., 2691.3, 2706., 2719.9, 2733.1, 2745.4, 2756.7, 2767.1, 2776.3, 2784.3, 2790.9, 2796.2, 2799.9, 2802., 2802.2, 2800.4, 2796.4, 2789.9, 2780.4, 2767.6, 2751., 2730., 2703.7, 2670.2, 2626.2, 2567.7, 2485.4, 2342.8, 2107.4};
	// Caso a temperatura seja T = 372 a entalpia especifica eh 2286.9  e se   T = 0.01 a entalpia especifica eh 2501.6
	
	// Valores tabulados para estado de saturação - propriedades na fase líquida (água)
	// Densidade -> indice 2; Viscosidade -> indice 3; Calor especifico -> indice 4; Condutividade termica -> indice 5; Entalpia especifica -> indice 6
	double LWDensity[39] = {999.8, 999.7, 998.3, 995.7, 992.3, 988., 983.2, 977.7, 971.6, 965.2, 958.1, 950.7, 942.9, 934.6, 925.8, 916.8, 907.3, 897.3, 886.9, 876., 864.7, 852.8, 840.3, 827.3, 813.6, 799.2, 783.9, 767.8, 750.5, 732.1, 712.2, 690.6, 666.9, 640.4, 610.2, 574.4, 527.5, 451.8, 315.5};
	double LWDynamicViscosity[39] = {1792., 1305., 1003., 798., 654.,547., 466., 403., 354., 315., 282., 255., 233., 214., 197., 183., 170., 160., 150., 142., 134., 128., 122., 116., 111., 107., 102., 98., 94., 90., 85.6, 81.6, 77.6, 73.6, 69.5, 65.2, 59.7, 52.2, 40.6};
	double LWSpecificHeat[39] = {4.217, 4.193, 4.182, 4.179, 4.179, 4.181, 4.185, 4.19, 4.197, 4.205, 4.216, 4.229, 4.245, 4.263, 4.285, 4.31, 4.339, 4.371, 4.408, 4.449, 4.497, 4.551, 4.613, 4.685, 4.769, 4.867, 4.983, 5.122, 5.29, 5.499, 5.762, 6.104, 6.565, 7.219, 8.233, 10.11, 14.58, 43.17, 1000.0};
	double LWThermalConductivity[39] = {56.5, 58.4, 60.2, 61.7, 63.1, 64.2, 65.2, 66., 66.9, 67.5, 67.9, 68.1, 68.5, 68.6, 68.6, 68.6, 68.3, 68., 67.6, 67., 66.4, 65.6, 64.7, 63.7, 62.6, 61.5, 60.2, 59.1, 57.8, 56.3, 54.7, 53., 51.2, 49.1, 47., 44.7, 42.5, 41.8, 91.4};
	double LWSpecificEnthalpy[39] = {-0.04, 41.99, 83.86, 125.66, 167.45, 209.26, 251.09, 292.97, 334.92, 376.94, 419.06, 461.32, 503.72, 546.31, 589.1, 632.15, 675.47, 719.12, 763.12, 807.52, 852.37, 897.74, 943.67, 990.26, 1037.6, 1085.8, 1134.9, 1185.2, 1236.8, 1290., 1345., 1402.4, 1462.6, 1526.5, 1595.5, 1671.9, 1764.2, 1890.2, 2107.4};
	// Caso a temperatura seja T = 0.01 a entalpia especifica é ZERO para a aguar em fase liquida, e se T = 372 a entalpia especifica eh 1935.6
	
	vector<vector<double> >::iterator it;
	int i = 0;
	for(it = fTable.begin(); it < fTable.end(); it++, i++)
	{
		vector<double> vec;
		vec.push_back(Temperature[i]);
		vec.push_back(Pression[i]);
		vec.push_back(LatentHeat[i]);
		vec.push_back(LWDensity[i]);
		vec.push_back(0.000001*LWDynamicViscosity[i]);
		vec.push_back(LWSpecificHeat[i]);
		vec.push_back(ConvertWTokJperSecond(0.01*LWThermalConductivity[i]));
		vec.push_back(LWSpecificEnthalpy[i]);
		vec.push_back(Density[i]);
		vec.push_back(0.000001*DynamicViscosity[i]);
		vec.push_back(SpecificHeat[i]);
		vec.push_back(ConvertWTokJperSecond(0.01*ThermalConductivity[i]));
		vec.push_back(SpecificEnthalpy[i]);
		// Copiando o vetor preenchido na linha correspondente
		*it = vec;
	}
}

// Physical properties to liquid water (fase aquosa) from tabulated data
double WaterDataInStateOfSaturation::getSaturationStateDensityToLiquidWater(double temp)
{
	fIndependentVarIndex = 0;
	double val = getInterpolatedValue(temp,3);
	return val;
}

double WaterDataInStateOfSaturation::getSaturationStateSpecificHeatToLiquidWater(double temp)
{
	fIndependentVarIndex = 0;
	double val = getInterpolatedValue(temp,5);

	return val;
}

// Physical properties to steam (fase gasosa) from tabulated data
double WaterDataInStateOfSaturation::getSaturationStateDensityToSteam(double temp)
{
	fIndependentVarIndex = 0;
	double val = getInterpolatedValue(temp,8);

	return val;
}

double WaterDataInStateOfSaturation::getSaturationStateViscosityToSteam(double temp)
{
	fIndependentVarIndex = 0;
	double val = getInterpolatedValue(temp,9);

	return val;
}

double WaterDataInStateOfSaturation::getSaturationStateSpecificHeatToSteam(double temp)
{
	fIndependentVarIndex = 0;
	double val = getInterpolatedValue(temp,10);

	return val;
}

double WaterDataInStateOfSaturation::getSaturationWater(int i, double quality, double So)
{
	double Sw = 1. - (quality + So);

	if(Sw < 0.0)
	{
		Sw = 0.0;
	}
	if(Sw > 1.0)
	{
		Sw = 1.0;
	}

	return Sw;
}





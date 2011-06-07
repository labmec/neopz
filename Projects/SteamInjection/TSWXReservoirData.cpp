//---------------------------------------------------------------------------


#pragma hdrstop

#include "TSWXReservoirData.h"

//---------------------------------------------------------------------------

/**
 * Metodos relacionados ao objeto reservatorio
 */
TSwxReservoirData::TSwxReservoirData()
{
		frw = -1.;
		fre = -1.;
		fHr = -1.;
		fPe = -1.;
		fKh = -1.;
		fPoros = -1.;
		fSo = -1.;
		fTr = -1.;
		fRockDensity = -1.;
		fRockSpecificHeat = -1.;
		fOilDensity = -1.;
		fOilSpecificHeat = -1.;

		// Preenchendo a tabela de valores para a permeabilidade relativa
		// Valores relacionados da temperatura e da pressão para o estado de saturação do fluido
		// Temperatura -> indice 0; Pressao -> indice 1; Calor Latente -> indice 7

		// Observar: A saturação de água não pode ser ZERO, pois existe a Saturação
		// irredutível de água, e também não pode ser UM pois existe a saturação mínima de óleo
		const int nData = 17;
		fTable.resize(nData);
		double WaterSaturation[nData] =
		{
			0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
			0.75, 0.8, 0.85, 0.9
		};
		double RelativePermeability_Vapour[nData] =
		{
			0.85, 0.71, 0.6, 0.51, 0.43, 0.36, 0.3, 0.25, 0.215, 0.187, 0.166,
			0.137, 0.11, 0.1, 0.065, 0.03, 0.001
		};
		double RelativePermeability_Water[nData] =
		{
			0.008, 0.011, 0.015, 0.023, 0.035, 0.048, 0.065, 0.083, 0.105, 0.127,
			0.152, 0.185, 0.23, 0.28, 0.39, 0.56, 0.91
		};
		for(int d = 0; d < nData; d++)
		{
				std::vector<double> vec(3);
				vec[0] = WaterSaturation[d];
				vec[1] = RelativePermeability_Vapour[d];
				vec[2] = RelativePermeability_Water[d];
				fTable[d] = vec;
		}
}

	void TSwxReservoirData::SetData(double rw, double re, double Hr, double Pe,
																	double Kh, double Poros, double So, double Tr,
																	double RockDensity, double RockSpecificHeat,
																	double OilDensity, double OilSpecificHeat, double OilViscosity)
	{
		frw = rw;
		fre = re;
		fHr = Hr;
		fPe = Pe;
		fKh = Kh;
		fPoros = Poros;
		fSo = So;
		fTr = Tr;
		fRockDensity = RockDensity;
		fRockSpecificHeat = RockSpecificHeat;
		fOilDensity = OilDensity;
		fOilSpecificHeat = OilSpecificHeat;
		fOilVisc = OilViscosity;
	}

double TSwxReservoirData::Rho_C_Estrela(WaterDataInStateOfSaturation &watertable, double temp, double Sw, double So)
{
	double oilDensity = fOilDensity;
	double oilSpecificHeat = fOilSpecificHeat;

	double Rho_C_Fluidos = So * oilDensity * oilSpecificHeat;

	Rho_C_Fluidos += (Sw * watertable.getSaturationStateDensityToLiquidWater(temp) * watertable.getSaturationStateSpecificHeatToLiquidWater(temp));
	Rho_C_Fluidos += ((1 - (Sw + So)) * watertable.getSaturationStateDensityToSteam(temp) * watertable.getSaturationStateSpecificHeatToSteam(temp));

	double val = (((1 - fPoros) * fRockDensity * fRockSpecificHeat)	+ fPoros * Rho_C_Fluidos);

	return val;
}

#pragma package(smart_init)

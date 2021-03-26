//---------------------------------------------------------------------------

#ifndef TSWXInputDataH
#define TSWXInputDataH
//---------------------------------------------------------------------------

#include "ThermalMethodsTables.h"
#include "ConvertionMethods.h"

#include "PropertiesTable.h"

#include "TSWXReservoirData.h"

#include "TSWXConfinementData.h"

#include "TSWXInjectionData.h"

#include "TSWXMaxSigmaTheta.h"

/** Dados de entrada */
class TSWXInputData
{

	public:

	TSWXInputData() : fReservoirData(), fInjectionData(), fConfinementData(), fWaterInSaturationState()
	{

		fMaxSigmaThetaData = new TSWXMaxSigmaTheta;

	}

	double getLiquidWaterViscosity()
	{
		double temp = getReservoirData().Temperature();
		if(temp > 100.)
		{
			temp = 100;
		}
		double waterVisc = getWaterInSaturationState().getViscosityOfLiquidWater(temp);

		return waterVisc;
	}

	TSwxReservoirData &getReservoirData() { return fReservoirData; }
	TSwxInjectionData &getInjectionData() { return fInjectionData; }
	TSwxConfinementData &getConfinementData() { return fConfinementData; }
	WaterDataInStateOfSaturation &getWaterInSaturationState() { return fWaterInSaturationState; }
	TSWXMaxSigmaTheta * getMaxSigmaThetaData() { return fMaxSigmaThetaData; }

	double getRho_C_Estrela(int i, double So);  // Sw, Sv, So sao as saturacoes de agua, vapor e oleo em uma celula. De 0 ate 1.

	private:

	TSwxReservoirData fReservoirData;
	TSwxInjectionData fInjectionData;
	TSwxConfinementData fConfinementData;
	WaterDataInStateOfSaturation fWaterInSaturationState;
	TSWXMaxSigmaTheta * fMaxSigmaThetaData;
 
};

#endif

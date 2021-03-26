//---------------------------------------------------------------------------


#pragma hdrstop

#include "TSWXConfinementData.h"
#include "ConvertionMethods.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

TSwxConfinementData::TSwxConfinementData(TSwxConfinementData &cpy) {
	fThermalConductivity = cpy.fThermalConductivity;
	fDensity = cpy.fDensity;
	fSpecificHeat = cpy.fSpecificHeat;
}

double TSwxConfinementData::getProductOfTheProperties() {
		double val = fThermalConductivity * fDensity * fSpecificHeat;
	return val;
}

void TSwxConfinementData::SetData(double condtermica, double densidade, double calorespecif)
{
	fThermalConductivity = ConvertWTokJperSecond(condtermica);
	fDensity = densidade;
	fSpecificHeat = calorespecif;
}

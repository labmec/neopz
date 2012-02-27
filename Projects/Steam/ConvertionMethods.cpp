/**
 * @file
 * @brief Contains the implementation of the convertion methods.
 */

#include "ConvertionMethods.h"

/// Length units
REAL ConvertFeetToMeter(REAL lenght) {
	return (lenght*0.30480061L);
}

/// Energy units
REAL ConvertBtuperFtCuboTokJpermCubo(REAL value) {
	return (37.259L*value);
}

REAL ConvertBtuTokJ(REAL value) {
	return value*1.055056L; 
}

REAL ConvertBtuPerHourTokJPerSecond(REAL value) {
	return value*(1.055056L/3600.0L);
}

//double ConvertBtuperFtCuboFTokJpermCuboK(double value) {
//	return (67.066*value);
//}

//double ConvertWTokJperH(double value) {
//	return (3.6*value);
//}

/// Temperature units
REAL ConvertCToF(REAL temperature) {
	return (1.8L*temperature + 32.L);
}

//double ConvertCToK(double temperature) {
//	return (temperature + 273.15);
//}

/// Time units
REAL ConvertHourToSecond(REAL horas) {
	return horas*3600L;
}

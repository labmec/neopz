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

REAL ConvertWTokJperH(REAL value) {
	return (3.6L*value);
}
REAL ConvertWTokJperSecond(REAL value) {
	return (0.001L*value);
}

/// Temperature units
REAL ConvertCToF(REAL temperature) {
	return (1.8L*temperature + 32.L);
}

REAL ConvertCToK(REAL temperature) {
	return (temperature + 273.15L);
}

/// Time units
REAL ConvertHourToSecond(REAL horas) {
	return horas*3600L;
}

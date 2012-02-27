/**
 * @file
 * @brief Contains some units conversions.
 */

#ifndef CONVERTIONMETHODSHH
#define CONVERTIONMETHODSHH

#include "pzreal.h"

/// Space measures
REAL ConvertFeetToMeter(REAL lenght);

/// Energy measures
REAL ConvertBtuperFtCuboTokJpermCubo(REAL value);
REAL ConvertBtuTokJ(REAL value);
REAL ConvertBtuPerHourTokJPerSecond(REAL value);

/// Temperature measures
REAL ConvertCToF(REAL temperature);

REAL ConvertCToK(REAL temperature);

REAL ConvertWTokJperH(REAL value);

inline REAL ConvertWTokJperH(REAL value) {
	return (3.6*value);
}

inline REAL ConvertCToK(REAL temperature) {
	return (temperature + 273.15);
}

/// Time measures
REAL ConvertHourToSecond(REAL horas);

#endif

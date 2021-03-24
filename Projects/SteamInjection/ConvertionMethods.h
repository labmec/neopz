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

REAL ConvertWTokJperSecond(REAL value);

/// Time measures
REAL ConvertHourToSecond(REAL horas);

#endif

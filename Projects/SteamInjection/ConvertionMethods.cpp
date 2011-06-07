/*
 *  ConvertionMethods.cpp
 *  FrenteDoVapor
 *
 *  Created by Jorge Calle on 11/04/10.
 *  Copyright 2010 Labmec. All rights reserved.
 *
 */

#include "ConvertionMethods.h"

double ConvertWTokJperSecond(double value) {
	return (0.001*value);
}

double ConvertCToK(double temperature) {
	return (temperature + 273.15);
}

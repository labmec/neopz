/*
 *  ConvertionMethods.cpp
 *  FrenteDoVapor
 *
 *  Created by Jorge Calle on 11/04/10.
 *  Copyright 2010 Labmec. All rights reserved.
 *
 */

#include "ConvertionMethods.h"

//--------- Unidades de comprimento --------------------
double ConvertFeetToMeter(double lenght) {
	return (lenght*0.30480061);
}

//------------ Unidades de energia -----------------------------
double ConvertBtuperFtCuboTokJpermCubo(double value) {
	return (37.259*value);
}

double ConvertBtuTokJ(double value) {
	return value*1.055056; 
}

double ConvertBtuPerHourTokJPerSecond(double value) {
	return value*(1.055056/3600.0);
}

//double ConvertBtuperFtCuboFTokJpermCuboK(double value) {
//	return (67.066*value);
//}

double ConvertWTokJperH(double value) {
	return (3.6*value);
}

//----------- Unidades de temperatura --------------------
double ConvertCToF(double temperature) {
	return (1.8*temperature + 32.);
}

double ConvertCToK(double temperature) {
	return (temperature + 273.15);
}

//----------- Unidades em tempo -----------------------
double ConvertHourToSecond(double horas) {
	return horas*3600;
}

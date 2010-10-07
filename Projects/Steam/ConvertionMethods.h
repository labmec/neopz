/*
 *  ConvertionMethods.h
 *  FrenteDoVapor
 *
 *  Created by Jorge Calle on 11/04/10.
 *  Copyright 2010 Labmec. All rights reserved.
 *
 */

#ifndef CONVERTIONMETHODSHH
#define CONVERTIONMETHODSHH

double ConvertFeetToMeter(double lenght);

double ConvertBtuperFtCuboTokJpermCubo(double value);
double ConvertBtuTokJ(double value);
double ConvertBtuPerHourTokJPerSecond(double value);

//double ConvertBtuperFtCuboFTokJpermCuboK(double value);

double ConvertCToF(double temperature);

double ConvertCToK(double temperature);

double ConvertWTokJperH(double value);


// Unidades em tempo
double ConvertHourToSecond(double horas);

#endif

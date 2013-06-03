//
//  pznlfluidstructureMaterials.h
//  PZ
//
//  Created by Cesar Lucci on 28/05/13.
//
//

#ifndef PZ_pznlfluidstructureMaterials_h
#define PZ_pznlfluidstructureMaterials_h

int const globReservMatId   = 1; //elastic
int const globPressureMatId = 2; //pressure
int const globMultiFisicMatId = 1;//multiphisics

int const globDirichletElastMatId = -1;
int const globBlockedXElastMatId  = -3;

int const globBCfluxIn  = -10; //bc pressure
int const globBCfluxOut = -20; //bc pressure

int const dirichlet = 0;
int const neumann   = 1;
int const mixed     = 2;

int const globDir_elast     = 11;
int const globMix_elast     = 20;
int const globNeum_pressure = 21;

#endif

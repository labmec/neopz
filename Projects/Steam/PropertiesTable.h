/*
 *  PropertiesTables.h
 *  FrenteDoVapor
 *
 *  Created by Jorge Calle on 5/05/10.
 *  Copyright 2010 Labmec. All rights reserved.
 *
 */

#ifndef PROPERTIESTABLESFORSOLIDSFLUIDSANDOILSHH
#define PROPERTIESTABLESFORSOLIDSFLUIDSANDOILSHH

//#include "string.h"

   struct PhysicalProperties {

	char MaterialName[64];
	
	double density;
	double specificheat;
	double thermalconductivity;
	
	double porosity;
	
	PhysicalProperties(char *name,double rho,double cp,double lamda,double poros);
	PhysicalProperties(int type,int rockid);

	void PutData(const char *name,double rho,double cp,double lamda,double poros);
};

struct TPBrScales
{
    static const double fPressureScale;
    static const double fEnergyScale;
    static const double fPermeability;
    static double fReferencePressure;
    
};

// Funções para inicializar as tabelas abaixo
//extern void InitializePhysicalPropertiesRocks();

// Variaveis externas para as tabelas das rochas da formação e confinantes
//extern PhysicalProperties FormationRock[20];
//extern PhysicalProperties ConfinantRock[10];

#endif
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
	
	double viscosity;
	double porosity;
	double PrandtlNumber;
	
	PhysicalProperties(char *name,double rho,double cp,double lamda,double poros,double mu,double prandtl);
	PhysicalProperties(int type,int rockid);

	void PutData(const char *name,double rho,double cp,double lamda,double poros,double mu,double prandtl);
};

// Funções para inicializar as tabelas abaixo
//extern void InitializePhysicalPropertiesRocks();

// Variaveis externas para as tabelas das rochas da formação e confinantes
//extern PhysicalProperties FormationRock[20];
//extern PhysicalProperties ConfinantRock[10];

#endif
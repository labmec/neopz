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
                //permeable      //impermeable
enum ERocktype {EPermeableLayer, EConfinantLayer};

struct PhysicalProperties {
	
	double density;
	double specificheat;
	double porosity;
	double thermalconductivity;

	PhysicalProperties(ERocktype type,int rockid);
	void PutData(double rho, double cp, double poros, double thermconduct);
};

// Funções para inicializar as tabelas abaixo
//extern void InitializePhysicalPropertiesRocks();

//Variaveis externas para as tabelas das rochas da formação e confinantes
//extern PhysicalProperties FormationRock[20];
//extern PhysicalProperties ConfinantRock[10];

#endif
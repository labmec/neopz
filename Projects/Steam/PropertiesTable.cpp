/*
 *  PropertiesTable.cpp
 *  FrenteDoVapor
 *
 *  Created by Jorge Calle on 5/05/10.
 *  Copyright 2010 Labmec. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "PropertiesTable.h"

// Inicializando as tabelas das rochas de formação e confinantes
//PhysicalProperties FormationRock[20];
//PhysicalProperties ConfinantRock[10];

// Métodos para a estrutura PhysicalProperties
PhysicalProperties::PhysicalProperties(char *name,double rho,double cp,double lambda,double poros) {
	strncpy(MaterialName,name,32);
	//MaterialName = name;
	density = rho;// Kg/m^3
	specificheat = cp; // KJ/(Kg C)
	thermalconductivity = lambda; // W/(m C)
	porosity = poros;
}
PhysicalProperties::PhysicalProperties(int type,int rockId) {
	switch(type) {
		case 0:
		{
			switch(rockId) {
				case 0:
					PutData("SandstoneSaturatedWithAir",2080.,0.766,0.877,0.196);
					break;
				case 1:
					PutData("SandstoneSaturatedWithWater",2275.,1.055,2.75,0.196);
					break;
				case 2:
					PutData("SilstoneSaturatedWithAir",1920.,0.854,0.685,0.199);
					break;
				case 3:
					PutData("ShaleSaturatedWithAir",2320.,0.804,1.04,0.071);
					break;
				case 4:
					PutData("ShaleSaturatedWithWater",2390.,0.892,1.69,0.071);
					break;
				case 5:
					PutData("LimestoneSaturatedWithAir",2195.,0.846,1.7,0.186);
					break;
				case 6:
					PutData("LimestoneSaturatedWithWater",2390.,1.114,3.55,0.186);
					break;
	
				case 11:
					PutData("DisaggregatedSandstoneSaturatedWithAir",1440.,0.837,0.493,0.40);
					break;
				case 12:
					PutData("DisaggregatedSandstoneSaturatedWithWater",1840.,1.566,1.82,0.40);
					break;
				case 13:
					PutData("DisaggregatedSandstoneSaturatedMix",1600.,1.0485,1.,0.40);
					break;
				case 14:
					PutData("DisaggregatedSandstoneSaturatedMixture",1945,1.7627,1.,0.40);
					break;
					
				case 15:
					PutData("DisaggregatedSilstoneSaturatedWithAir",1540.,0.846,0.585,0.36);
					break;
				case 16:
					PutData("DisaggregatedSilstoneSaturatedWithWater",1890.,1.47,1.79,0.36);
					break;
				default:
					PutData("UNKNOW",0.,0.,0.,0.);
					break;
			}
		}
			break;
		case 1:
		{
			switch(rockId) {
				case 0:
					PutData("Aluminum(pure)",2707.,0.896,204.,0.0);
					break;
				case 1:
					PutData("Copper(pure)",8954.,0.383,386.,0.0);
					break;
				case 2:
					PutData("Granite",2700.,0.82,2.,0.0);
					break;
				case 3:
					PutData("MARBLE_LOW",2500.,0.8,2.1,0.0);
					break;
				case 4:
					PutData("MARBLE_HIGH",2700.,0.8,2.9,0.0);
					break;
				case 5:
					PutData("BASALT_LOW",2650.,0.84,2.2,0.0);
					break;
				case 6:
					PutData("BASALT_HIGH",2700.,0.84,2.5,0.0);
					break;
				case 7:
					PutData("BASALT_TEST",2682.,0.84,2.41,0.0);
					break;
				default:
					PutData("UNKNOW",0.,0.,0.,0.);
					break;
			}
		}
			break;
		default:
			PutData("UNKNOW",0.,0.,0.,0.);
			break;
	}
}

void PhysicalProperties::PutData(const char *name,double rho,double cp,double lamda,double poros) {
	memset(MaterialName,'\0',64);
	strncpy(MaterialName,name,63);
	density = rho;
	specificheat = cp;
	thermalconductivity = lamda;
	porosity = poros;
}

const double TPBrScales::fPressureScale = 100.;
const double TPBrScales::fEnergyScale = 1.e4;
const double TPBrScales::fPermeability = 1.e-12;
double TPBrScales::fReferencePressure = 1.e6;



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
PhysicalProperties::PhysicalProperties(ERocktype type, int rockId) {
	switch(type) {
		case (EPermeableLayer):
		{
			switch(rockId) {
				case 0://"SandstoneSaturatedWithAir"
					PutData(2080.,0.766,0.196,0.877);
					break;
				case 1://"SandstoneSaturatedWithWater"
					PutData(2275.,1.055,0.196,2.75);
					break;
				case 2://"SilstoneSaturatedWithAir"
					PutData(1920.,0.854,0.199,0.685);
					break;
				case 3://"ShaleSaturatedWithAir"
					PutData(2320.,0.804,0.071,1.04);
					break;
				case 4://"ShaleSaturatedWithWater"
					PutData(2390.,0.892,0.071,1.69);
					break;
				case 5://"LimestoneSaturatedWithAir"
					PutData(2195.,0.846,0.186,1.7);
					break;
				case 6://"LimestoneSaturatedWithWater"
					PutData(2390.,1.114,0.186,3.55);
					break;
				case 11://"DisaggregatedSandstoneSaturatedWithAir"
					PutData(1440.,0.837,0.40,0.493);
					break;
				case 12://"DisaggregatedSandstoneSaturatedWithWater"
					PutData(1840.,1.566,0.40,1.82);
					break;
				case 13://"DisaggregatedSandstoneSaturatedMix"
					PutData(1600.,1.0485,0.40,1.);
					break;
				case 14://"DisaggregatedSandstoneSaturatedMixture"
					PutData(1945,1.7627,0.40,1.);
					break;
				case 15://"DisaggregatedSilstoneSaturatedWithAir"
					PutData(1540.,0.846,0.36,0.585);
					break;
				case 16://"DisaggregatedSilstoneSaturatedWithWater"
					PutData(1890.,1.47,0.36,1.79);
					break;
				default://"UNKNOW"
					PutData(0., 0., 0.,0.);
					break;
			}
		}
			break;

		case EConfinantLayer:
		{
			switch(rockId) {
				case 0://"Aluminum(pure)"
					PutData(2707.,0.896,0.,204.);
					break;
				case 1://"Copper(pure)"
					PutData(8954.,0.383,0.,386.);
					break;
				case 2://"Granite"
					PutData(2700.,0.82,0.,2.);
					break;
				case 3://"MARBLE_LOW"
					PutData(2500.,0.8,0.,2.1);
					break;
				case 4://"MARBLE_HIGH"
					PutData(2700.,0.8,0.,2.9);
					break;
				case 5://"BASALT_LOW"
					PutData(2650.,0.84,0.,2.2);
					break;
				case 6://"BASALT_HIGH"
					PutData(2700.,0.84,0.,2.5);
					break;
				case 7://"BASALT_TEST"
					PutData(2682.,0.84,0.,2.41);
					break;
				default://"UNKNOW"
					PutData(0.,0.,0.,0.);
					break;
			}
		}
			break;
		default://"UNKNOW"
			PutData(0.,0.,0.,0.);
			break;
	}
}

void PhysicalProperties::PutData(double rho, double cp, double poros, double thermconduct) {
	density = rho;
	specificheat = cp;
	porosity = poros;
	thermalconductivity = thermconduct;
}


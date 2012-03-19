//
//  TPBrSteamSimulation.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/15/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//
#include <fstream>

#include "TPBrSteamSimulation.h"
#include "tpbrsteammesh.h"
#include "ThermalMethodsTables.h"
#include "pzgengrid.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("br.steamsimulation"));
#endif

WaterDataInStateOfSaturation waterdata;
//OilData oildata;

using namespace std;

void MobilidadeRelativa();



int main()
{
#ifdef LOG4CXX
	InitializePZLOG("../log4cxx.cfg");
#endif
    MobilidadeRelativa();
#ifdef _AUTODIFF
    //    TPBrSteamMesh::TPBrSteamMesh(int numcells, REAL temperature, REAL pressure, REAL WellRadius, REAL ReservoirRadius, REAL oilsaturation)
    REAL wellradius = 0.15;
    REAL ReservoirRadius = 300;
    int ncells = 300;
    REAL reservoirheight = 30.;
    REAL reservoirtemperature = 30.;
    REAL farfieldpressure = 1.e5;
    TPBrSteamMesh mesh(ncells,reservoirtemperature,farfieldpressure,wellradius,ReservoirRadius,reservoirheight,0.);
    mesh.SetWaterInjection(0.7*reservoirheight, reservoirtemperature);
    if(ncells > 1)
    {
        REAL progression = TPZGenGrid::GeometricProgression(0.04, ReservoirRadius-0.15, ncells);
        mesh.SetGeometricProgression(progression);        
    }
    std::cout << "First cell size = " << mesh.CellSize(0) << std::endl;
    mesh.Print();
    int neq = mesh.NumEquations();
    
    TPZFMatrix<REAL> tangent(neq,neq),residual(neq,1);
/*    TPZStack<REAL> scales;
    mesh.StateScales(scales);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Statescales " << scales;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
 */
    mesh.SetSteamQuality(0.1,3.e5);
    mesh.SetTimeStep(0.);
//    mesh.ComputeTangent(tangent, residual,scales);
//    {
//        std::ofstream tfile("step.nb");
//        tangent.Print("tangent = ", tfile,EMathematicaInput);
//        residual.Print("rhs = ",tfile,EMathematicaInput);
//    }
//    mesh.Iterate();
//    mesh.fPrevState = mesh.fNextState;
    mesh.TimeStep(0.);
    REAL loctime=0.;
    mesh.PrintAll(loctime);
    REAL PrintInterval = 1000.;
    REAL delt = 50.;
    REAL nextime = loctime+PrintInterval;
    while(loctime <= 10000.)
    {
        int again = 1;
        while(again)
        {
            try {
                REAL step = delt;
                if(loctime+delt > nextime) 
                {
                    step = nextime-loctime;
                }
                mesh.TimeStep(step);
                if(loctime+delt > nextime) 
                {
                    mesh.PrintAll(loctime+step);
                    nextime += PrintInterval;
                }
                else
                {
                    delt *= 1.5;
                }
                loctime += step;
                std::cout << "Time " << loctime << " delt = " << delt << std::endl;
                again = 0;
            } catch (...) {
                mesh.fNextState = mesh.fPrevState;
                delt /= 10.;
            }
        }
    }
//    mesh.Print();
#endif
     
    return 0;
}

void MobilidadeRelativa()
{
    ofstream sout("$\\frac{\\rho_steam \\nu_water}{\\nu_steam \\rho_water.txt");
    REAL Temp = 100;
    for (Temp = 100.; Temp <= 370; Temp += 20) {
        REAL viscwater = waterdata.getSaturationStateViscosityToLiquidWater(Temp);
        REAL viscsteam = waterdata.getSaturationStateViscosityToSteam(Temp);
        REAL denswater = waterdata.getSaturationStateDensityToLiquidWater(Temp);
        REAL denssteam = waterdata.getSaturationStateDensityToSteam(Temp);
        REAL relmob = (denssteam/viscsteam)/(denswater/viscwater);
        sout << Temp << " " << relmob <<std::endl;
    }
}


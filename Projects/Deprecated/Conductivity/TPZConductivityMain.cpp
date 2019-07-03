//
//  TPZConductivityMain.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/15/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "TPZConductivityMain.h"
#include "TPZConductivityProblem.h"
#include <map>
#include <iostream>

int main()
{
    /// evolucao da conductividade com a relacao area/vazio
    TPZConductivityProblem problem;
    REAL bridgevoidratio = problem.GetBridgeVoidRatio();
    std::map<REAL,REAL> bridgeflux;
    int i;
    for (i=1; i<10; i++) {
        problem.SetBridgeVoidRatio(bridgevoidratio*i);
        bridgeflux[bridgevoidratio*i] = problem.ComputeFlux();
    }
    std::cout << "Influence of BridgeVoidRatio ";
    std::map<REAL,REAL>::iterator it = bridgeflux.begin();
    while (it != bridgeflux.end()) {
        std::cout << it->first << ' ' << it->second << ' ';
        it++;
    }
    std::cout << std::endl;
}
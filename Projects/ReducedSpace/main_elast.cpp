//
//  File.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/31/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//


#include "pzfstrmatrix.h"
#include "toolstransienttime.h"

#include "TPZVTKGeoMesh.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
//#endif



int main(int argc, char *argv[])
{
    //Propagation criterion
    REAL Lx = 400.;
    REAL Ly = 400.;
    REAL Lf = 50.;
    REAL Hf = 1.;
    REAL Young = 3.9E4;
    REAL Poiss = 0.25;
    REAL Fx = 0.;
    REAL Fy = 0.;
    REAL Visc = 0.001E-6;
    REAL SigN = 61.5;
    REAL QinjTot  = -0.2;//Lembre-se que a divisao por 2 (1 asa) e por Hf (na secao de 1 asa) eh feita no kernel. Aqui vai Qinj total mesmo!!!
    REAL Ttot = 50.;
    REAL deltaT = Ttot/20.;
    REAL Cl = 0.005;
    REAL Pe = 10.;
    REAL SigmaConf = 11.;
    REAL Pref = 60000.;
    REAL vsp = 0.001;
    REAL KIc = 300.;
    int p = 2;
    
    ToolsTransient ToolTrans(p, Lx, Ly, Lf, Hf, Young, Poiss, Fx, Fy, Visc, SigN, QinjTot, Ttot, deltaT, Cl, Pe, SigmaConf, Pref, vsp, KIc);
	
	ToolTrans.Run();
    
    return 0;
}



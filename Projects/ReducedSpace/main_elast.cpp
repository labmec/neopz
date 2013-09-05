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
#include "TPZRefPatternDataBase.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
//#endif



int main(int argc, char *argv[])
{
    //Propagation criterion
    REAL Lx = 1000.;
    REAL Ly = 600.;
    REAL Lf = 50.;
    REAL Hf = 1.;
    REAL Young = 3.9E4;
    REAL Poiss = 0.25;
    REAL Fx = 0.;
    REAL Fy = 0.;
    REAL preStressXX = 0.;
    REAL preStressXY = 0.;
    REAL preStressYY = -50.; //(negativo : estado de compressao)
    
    int NStripes = 1;
    REAL Visc = 0.001E-6;
    
    REAL SigN = 6.15;
    
    /**
     * Lembre-se que a divisao por 2 (1 asa) e por Hf (na secao de 1 asa) eh feita no kernel.
     * Aqui vai Qinj total mesmo!!!
     */
    REAL QinjTot  = -0.5;

    REAL Ttot = 20.; /** em segundos */
    REAL maxDeltaT = 4.; /** em segundos */
    int nTimes = 2; /**  */
    
    REAL Cl = 0.005;
    REAL Pe = 10.;
    REAL SigmaConf = 11.;
    REAL Pref = 60000.;
    REAL vsp = 0.001;
    REAL KIc = 1200.;
    REAL Jradius = 0.5;
    
    int p = 2;
    
    globFractInputData.SetData(Lx, Ly, Lf, Hf, Young, Poiss, Fx, Fy, preStressXX, preStressXY, preStressYY, NStripes, Visc, SigN,
                               QinjTot, Ttot, maxDeltaT, nTimes, Cl, Pe, SigmaConf, Pref, vsp, KIc, Jradius);
    ToolsTransient ToolTrans(p);
    
    std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
    std::cout << "****\n";
    std::cout << "*****\n";
    std::cout << "******\n";
    std::cout << "*******\n";
    std::cout << "********\n";
    std::cout << "**********\n";
    std::cout << "*************\n";
    std::cout << "*******************\n";
    std::cout << "*************************\n";
    std::cout << "*************************************\n";
    std::cout << "*************************************************\n";
    std::cout << "*******************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    std::cout << "Lembre-se de deixar os minT, maxT, actT etc como inteiros!!!\n";
    std::cout << "*******************************************************\n";
    std::cout << "*******************************************************\n";
    std::cout << "*******************************************************\n";
    std::cout << "*******************************************************\n";
    std::cout << "*******************************************************\n";
//    LEIA ACIMA!!!
	
	ToolTrans.Run();
    
    return 0;
}



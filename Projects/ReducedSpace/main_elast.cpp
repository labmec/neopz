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
//static PZLogger logger("pz.reducedspace.data");
//#endif



int main(int argc, char *argv[])
{
    //Propagation criterion
    REAL Lx = 80.;
    REAL Ly = 60.;
    REAL Lf = 5.;
    REAL Hf = 20.;
    REAL Lmax_edge = 1.;
    REAL Young1 = 4.1368543680E3;
    REAL Poiss1 = 0.25;
    REAL Young2 = 4.1368543680E3;
    REAL Poiss2 = 0.25;
    REAL Xinterface = 75.;
    REAL Fx = 0.;
    REAL Fy = 0.;
    REAL preStressXX = 0.;
    REAL preStressXY = 0.;
    REAL preStressYY = -3.4528254983; //(negativo : estado de compressao)
    
    int NStripes = 1;
    REAL Visc = 200.02E-10;
    
    REAL SigN = 1.;
    
    /**
     * Lembre-se que a divisao por 2 (1 asa) e por Hf (na secao de 1 asa) eh feita no kernel.
     * Aqui vai Qinj total mesmo (no poco)!!!
     */
    REAL QinjTot  = -0.05333333333333;

    REAL Ttot = 20.; /** em segundos */
    REAL maxDeltaT = 1.; /** em segundos */
    int nTimes = 10; /** quantidade de divisao do maxDeltaT para definir minDeltaT (minDeltaT = maxDeltaT/nTimes) */
    
    REAL Cl = 0.00019674755398733676;
    REAL Pe = 0.;
    REAL SigmaConf = -preStressYY;
    REAL Pref = 1.;
    REAL vsp = 0;
    REAL KIc = 0.109;
    REAL Jradius = 1.;
    
    int p = 2;
    
    globFractInputData.SetData(Lx, Ly, Lf, Hf, Lmax_edge, Young1, Poiss1, Young2, Poiss2, Xinterface,
                               Fx, Fy, preStressXX, preStressXY, preStressYY, NStripes, Visc, SigN,
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



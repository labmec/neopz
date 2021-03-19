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
#include "CohesiveTests.h"

#include "pzlog.h"

#ifdef PZ_LOG
static PZLogger logger("pz.reducedspace.data");
#endif


void CohesiveTest();
int mainFrac(int argc, char *argv[]);


int main(int argc, char *argv[])
{	
	if(1){
		mainFrac(argc,argv);
	}
	else {
		CohesiveTest();
	}
}
int mainFrac(int argc, char *argv[])
{
  //Propagation criterion
  
  REAL Lx = 4;
  REAL Ly = 2.5;
  REAL Lf = 0.2; //1
  REAL Hf = 1.;
  int ndivV = 10; // division in x for PGmesh
  int ndivH = 25; // division in y for PGmesh
  REAL q = 1.01; // PG order
  REAL Lmax_edge = 0.025;

  REAL Young1 = 25.e3;//9
  REAL Poiss1 = 0.2;
  
  
  REAL Young2 = 3.9E3;
  REAL Poiss2 = 0.25;
  REAL Xinterface = Lx + 10; // BECAUSE I WILL NOT USE INTERFACE
  
  REAL Fx = 0.;
  REAL Fy = 0.;
	// Now im using net pressure. So PreStress must be zero!!!
  REAL preStressXX = 0.; // era -50
  REAL preStressXY = 0.;
  REAL preStressYY = 0.; // era 25 //(positivo : estado de compressao)
  
  int NStripes = 1;
  REAL Visc = 5.e-8;//-2
  
  REAL SigN = 1.;
  
  //TPZMaterial::gBigNumber = 1.e15;
   
   // Lembre-se que a divisao por 2 (1 asa) e por Hf (na secao de 1 asa) eh feita no kernel.
   //Aqui vai Qinj total mesmo (no poco)!!!
  REAL QinjTot  = -0.0004; //-0.002;
  
  REAL Ttot = 6; // em segundos
  REAL maxDeltaT = 0.2; // em segundos
  int nTimes = 14; // quantidade de divisao do maxDeltaT para definir minDeltaT (minDeltaT = maxDeltaT/nTimes)
  
  // LeakOff Param
  REAL Cl = 0.001; //0.005
  REAL Pe = 0.; // 10.
  REAL SigmaConf = -preStressYY;
  REAL Pref = 60000.;
  REAL vsp = 0.*0.000001; //0.001
  REAL KIc = 25.;
  REAL Jradius = 0.5;
  bool usingLeakOff = false;
  
  // porder
  int p = 1;
  
	//cohesive param
	REAL DeltaC = 1.75*0.0001024; //0.0001024
	REAL DeltaT = 0.2 * DeltaC;
	REAL SigmaT = 3;//30 //3
  
  //MohrCoulomb parameters
  REAL cohesion = 0.13 * 5.77; //article is 5.77
  REAL phiMC = 30.*M_PI/180.; // article is 30 degres
	
	int NThreadsForAssemble = 8; // if set 0 it will be serial

  globFractInputData.SetData(Lx, Ly, Lf, Hf, Lmax_edge, Young1, Poiss1, Young2, Poiss2, Xinterface,
                             Fx, Fy, preStressXX, preStressXY, preStressYY, NStripes, Visc, SigN,
                             QinjTot, Ttot, maxDeltaT, nTimes, Cl, Pe, SigmaConf, Pref, vsp, KIc, Jradius,ndivV,ndivH,q,DeltaC,DeltaT,SigmaT,NThreadsForAssemble);
  
  globFractInputData.SetMohrCoulombData(cohesion,phiMC);  //Plastic Model is set here
  //globFractInputData.SetSandlerData();  //Plastic Model is set here
  
  globFractInputData.SetUsingLeakOff(usingLeakOff);
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

void CohesiveTest()
{
	ElastNLTestWithCohesive();
}

/* Dados do caju
 // APAGAR
 int ndivV = 10; // division in x for PGmesh
 int ndivH = 25; // division in y for PGm
 esh
 REAL q = 1.01; // PG order
 //Propagation criterion
 REAL Lx = 1000.;
 REAL Ly = 600.;
 REAL Lf = 51.;
 REAL Hf = 1.;
 REAL Lmax_edge = 10.;
 REAL Young1 = 3.9E4;
 REAL Poiss1 = 0.25;
 REAL Young2 = 3.9E5;
 REAL Poiss2 = 0.25;
 REAL Xinterface = Lx+10;
 REAL Fx = 0.;
 REAL Fy = 0.;
 REAL preStressXX = 0.;
 REAL preStressXY = 0.;
 REAL preStressYY = -50.; //(negativo : estado de compressao)
 
 int NStripes = 1;
 REAL Visc = 0.001E-6;
 
 REAL SigN = 1.;
 
 
 REAL QinjTot  = -0.5;
 
 REAL Ttot = 20.;
 REAL maxDeltaT = 4.;
 int nTimes = 1;
 
 REAL Cl = 0.005;
 REAL Pe = 10.;
 REAL SigmaConf = -preStressYY;
 REAL Pref = 60000.;
 REAL vsp = 0.001;
 REAL KIc = 25.;
 REAL Jradius = 0.5;
 
 int p = 1;
 // APAGAR

 */

/* Dados que funcionam
 REAL Lx = 10;
 REAL Ly = 5.;
 REAL Lf = 1;
 REAL Hf = 1.;
 int ndivV = 10; // division in x for PGmesh
 int ndivH = 25; // division in y for PGmesh
 REAL q = 1.01; // PG order
 
 REAL Lmax_edge = 0.2;
 REAL Young1 = 25.e3;//9
 REAL Poiss1 = 0.2;
 
 
 REAL Young2 = 3.9E3;
 REAL Poiss2 = 0.25;
 REAL Xinterface = Lx + 10; // BECAUSE I WILL NOT USE INTERFACE
 
 REAL Fx = 0.;
 REAL Fy = 0.;
 // Now im using net pressure. So PreStress must be zero!!!
 REAL preStressXX = 0.; // era -50
 REAL preStressXY = 0.;
 REAL preStressYY = 0.; // era 25 //(positivo : estado de compressao)
 
 int NStripes = 1;
 REAL Visc = 5.e-8;//-2
 
 REAL SigN = 1.;
 
 //TPZMaterial::gBigNumber = 1.e15;
 
 // Lembre-se que a divisao por 2 (1 asa) e por Hf (na secao de 1 asa) eh feita no kernel.
 //Aqui vai Qinj total mesmo (no poco)!!!
 REAL QinjTot  = -0.002;
 
 REAL Ttot = 5.; // em segundos
 REAL maxDeltaT = 0.5; // em segundos
 int nTimes = 1; // quantidade de divisao do maxDeltaT para definir minDeltaT (minDeltaT = maxDeltaT/nTimes)
 
 // LeakOff Param
 REAL Cl = 0.001; //0.005
 REAL Pe = 10.;
 REAL SigmaConf = -preStressYY;
 REAL Pref = 60000.;
 REAL vsp = 0.00001; //0.001
 REAL KIc = 25.;
 REAL Jradius = 0.5;
 bool usingLeakOff = true;
 
 // porder
 int p = 1;
 
 //cohesive param
 REAL DeltaC = 0.001024; //0.0001024
 REAL DeltaT = 0.2 * DeltaC;
 REAL SigmaT = 60;//30 //3
 
 //MohrCoulomb parameters
 REAL cohesion = 0.5 * 5.77; //article is 5.77
 REAL phiMC = 30.*M_PI/180.; // article is 30 degres
 
 int NThreadsForAssemble = 0; // if set 0 it will be serial
 */

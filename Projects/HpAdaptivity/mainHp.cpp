/*
 * @file
 * @author Denise de Siqueira
 * @since 6/9/11.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzanalysiserror.h"
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "pzmatrix.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "pzinterpolationspace.h"
#include "pzsubcmesh.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzbndcond.h"
#include "pzmaterial.h"
#include "tpzquadrilateral.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include "Tools.h"

#include "pzskylstrmatrix.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("HpAdaptivity.main"));

#endif

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads);

using namespace std;

/** Resolver o problema do tipo 
 * -Laplac(u) = 0
 * du/dn = lambda u em todo contorno
 */


using namespace std;
int main()
{
	
//#ifdef LOG4CXX
//	{
//		InitializePZLOG();
//		std::stringstream sout;
//		sout<< "Testando p adaptatividade"<<endl;
//		LOGPZ_DEBUG(logger, sout.str().c_str());
//	}
//#endif
	std::ofstream erro("TaxaArcTanTriangUni.txt");
    //std::ofstream erro("TaxaArcTanQuadUni.txt");
 //   std::ofstream erro("TaxaArcTanTriangNaoUni.txt");
    	
	TPZVec<REAL> calcErro;
	for (int porder=1; porder<5; porder++) {
		
		erro<<"ordem "<<porder <<std::endl;
			for(int h=1;h<8;h++){
			erro<<std::endl;
			
			//1. Criacao da malha geom. e computacional
					bool hrefine=false;//true nao uniforme
					bool prefine=false;
                TPZGeoMesh *gmesh = MalhaGeoT(h,hrefine);

                std::ofstream filemesh("MalhaGeoArcTan.vtk");
               PrintGMeshVTK( gmesh, filemesh);
//
//              RefiningNearCircunference(2,gmesh,h,1);
//                std::ofstream filemesh2("MalhaGeoQArcTanRefeineNearCirc.vtk");
//                PrintGMeshVTK( gmesh, filemesh2);
       

        TPZCompMesh *cmesh = CompMeshPAdap(*gmesh,porder,prefine);
                int nDofs;
                nDofs=cmesh->NEquations();
                
                erro<< "\n NRefinamento "<<h<< "   NDofs "<<nDofs<<std::endl;
            
        //2. Resolve o problema
		
		 TPZAnalysis analysis(cmesh);
         SolveSyst(analysis, cmesh, 8);
        // SolveLU ( analysis );
			
			
			
			//3. Calcula o erro 
					
//			TPZVec<REAL> calcErro;
//            analysis.SetExact(*SolArcTan);
//			analysis.PostProcess(calcErro,erro);
                
               ErrorHDiv(cmesh, erro);

        
			
			//4. visualizacao grafica usando vtk
        
//			 TPZVec<std::string> scalnames(4), vecnames(2);
//			 scalnames[3] = "Divergence";
//			 scalnames[2] = "ExactDiv";
//			 scalnames[0] = "Pressure";//"Solution";//
//			 scalnames[1] = "ExactPressure";
//			 //scalnames[1] = "ExactDiv";
//			 
//			 
//			 vecnames[0] = "ExactFlux";
//			 vecnames[1] = "Flux";//"Derivative";//
//                 
//			
//			 
//			 
//			 //vecnames[0] = "Derivate";
//               //  std::string plotfile("GraficArcTanHpAdaptivity.vtk");
//                 const int dim = 2;
//                 int div = 2;
//                 
//                 char buf[256] ;
//                 sprintf(buf,"ArcTanPhilUniMeshT_porder%d_h%d.vtk",porder,h);
//                 analysis.DefineGraphMesh(dim,scalnames,vecnames,buf);
//                 analysis.PostProcess(div);
               
			 
			 
			 /*std::string plotfile("GraficArcTanHpAdaptivity.vtk");
			 const int dim = 2;
			 int div = 2;
			 analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
			 analysis.PostProcess(div);
           */  
        
		}}
	
	
	return 0;
}

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads)
{
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
    full.SetNumThreads(numthreads);
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
    //	ofstream file("Solutout");
    //	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}



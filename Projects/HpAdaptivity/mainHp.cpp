/*
 * @file
 * @author Denise de Siqueira
 * @since 6/9/11.
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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
#include "TPZMaterial.h"
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
int main2()
{
	
#ifdef LOG4CXX
	{
		InitializePZLOG();
		std::stringstream sout;
		sout<< "Testando p adaptatividade"<<endl;
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif
    
    //std::ofstream erro("TaxaArcTanQuadUni.txt");
   // std::ofstream erro("TaxaArcTanTriangNaoUni.txt");
    //std::ofstream erro("TaxaProbJuan.txt");
    //std::ofstream erro("TaxaCNMACQuadUni.txt");
    //    std::ofstream erro("TaxaCNMACQuadNaoUni.txt");
    std::ofstream erro("TaxaCNMACTriangUni.txt");
    
   //bool ftriang=false;
    bool prefine=false;
//    bool nouniform=false;//false para ref uniforme
//    REAL Lx=1;
//    REAL Ly=1;
    
    
    
    //TPZGeoMesh *gmesh = GMesh(ftriang, Lx,  Ly);
    //UniformRefine(gmesh, 1);
    //RefiningNearCircunference(2,gmesh,1,1);
    //std::ofstream filemesh2("MalhaGeoIni.vtk");
    //PrintGMeshVTK( gmesh, filemesh2);
   // TPZGeoMesh *gmesh =GMesh(ftriang, Lx,  Ly);//MalhaGeoT(h,hrefine);
    TPZGeoMesh *gmesh =MalhaGeo2(1);
    ofstream arg("gmeshInicial.txt");
    gmesh->Print(arg);
    
	TPZVec<REAL> calcErro;
	for (int porder=1; porder<2; porder++) {
        
//        TPZCompMesh *cmesh = CompMeshPAdap(*gmesh,porder,prefine);
//        ofstream arg2("CmeshInicial.txt");
//        cmesh->Print(arg2);
//        
//        SetDifferentOrderP(cmesh,porder);
//        ofstream arg22("CmeshDifOrder.txt");
//        cmesh->Print(arg22);

        
        erro<<"ordem "<<porder <<std::endl;
        for(int h=1;h<4;h++){
        //    erro << "NRef " << h << std::endl;
            
            erro<<std::endl;
            
            //1. Criacao da malha geom. e computacional
 //           bool hrefine=true;//true nao uniforme
           // TPZGeoMesh *gmesh =MalhaGeo2(h);//MalhaGeoT(h,nouniform);//MalhaGeo(h,nouniform);//GMesh(ftriang, Lx,  Ly);//MalhaGeo(h,nonuniform);
            

//            UniformRefine(gmesh, h);
//            ofstream arg3("gmeshRefinada.txt");
//            gmesh->Print(arg3);
//         NoUniformRefine(gmesh, h);
//            std::ofstream filemesh("MalhaGeoNaoUniT.vtk");
//            std::ofstream filemesh("MalhaGeoUniQCNMAC.vtk");
            //std::ofstream filemesh("MalhaGeoNaoUniQCNMAC.vtk");
           std::ofstream filemesh("MalhaGeo4ElCNMAC.vtk");
            PrintGMeshVTK( gmesh, filemesh);
            

//            RefiningNearCircunference(2,gmesh,h,1);
//            std::ofstream filemesh2("MalhaGeoQArcTanRefeineNearCirc.vtk");
//            PrintGMeshVTK( gmesh, filemesh2);
            
            TPZCompMesh *cmesh = CompMeshPAdap(*gmesh,porder,prefine);//CompMeshPAdapJuan(*gmesh,porder,prefine);
            ofstream arg2("CmeshInicial.txt");
            cmesh->Print(arg2);
            
            SetDifferentOrderPMesh4Elem(cmesh,porder);
            ofstream arg22("CmeshDifOrder.txt");
            cmesh->Print(arg22);
            
            
            int nDofs;
            nDofs=cmesh->NEquations();
            
            erro<< "\nNRefinamento "<<h<< "   NDofs "<<nDofs<<std::endl;
            
            std::cout << "Neq = " << nDofs << std::endl;
            
            //2. Resolve o problema
            
            TPZAnalysis analysis(cmesh);
            int numthreads = 1;
            SolveSyst(analysis, cmesh, numthreads);
            // SolveLU ( analysis );
            
            
            
            //3. Calcula o erro

//			TPZVec<REAL> calcErro;
//            analysis.SetExact(*SolExata);
//			analysis.PostProcess(calcErro,erro);
            
            ErrorHDiv(cmesh, erro);
            
            UniformRefine(gmesh, h);
            ofstream arg3("gmeshRefinada.txt");
            gmesh->Print(arg3);
            
      
            /*4. visualizacao grafica usando vtk
            
            TPZVec<std::string> scalnames(5), vecnames(2);
            scalnames[0] = "Divergence";
            scalnames[1] = "ExactDiv";
            scalnames[2] = "Pressure";//"Solution";//
            scalnames[3] = "ExactPressure";
            scalnames[4] = "ExactDiv";
            
            
            vecnames[0] = "ExactFlux";
            vecnames[1] = "Flux";//"Derivative";//

            const int dim = 2;
            int div = 2;
            
            char buf[256] ;
            //sprintf(buf,"ArcTanPhilUniMeshQ_porder%d_h%d.vtk",porder,h);
           // sprintf(buf,"ProblemaJuan_porder%d_h%d.vtk",porder,h);
            sprintf(buf,"ProblemaCNAMC_porder%d_h%d.vtk",porder,h);
            analysis.DefineGraphMesh(dim,scalnames,vecnames,buf);
            analysis.PostProcess(div);
 
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                cmesh->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            */
            
            
            //              std::string plotfile("GraficArcTanHpAdaptivity.vtk");
            //              const int dim = 2;
            //              int div = 2;
            //              analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
            //              analysis.PostProcess(div);
            //
            //        UniformRefine(gmesh, 1);
            //         std::ofstream filemesh2("MalhaGeoDepois.vtk");
            //        PrintGMeshVTK( gmesh, filemesh2);
           // delete cmesh;
           // delete gmesh;
            
        }
    }
	
	
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



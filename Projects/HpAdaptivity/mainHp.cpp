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
	
#ifdef LOG4CXX
	{
		InitializePZLOG();
		std::stringstream sout;
		sout<< "Testando p adaptatividade"<<endl;
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif
    
    //std::ofstream erro("TaxaArcTanQuadUni.txt");
    //std::ofstream erro("TaxaArcTanTriangNaoUni.txt",ios::app);
    std::ofstream erro("TaxaProbJuan.txt",ios::app);
    
    bool ftriang=false;
    REAL Lx=1;
    REAL Ly=0.5;
    
    
    
    //TPZGeoMesh *gmesh = GMesh(ftriang, Lx,  Ly);
    //UniformRefine(gmesh, 1);
    //RefiningNearCircunference(2,gmesh,1,1);
    //std::ofstream filemesh2("MalhaGeoIni.vtk");
    //PrintGMeshVTK( gmesh, filemesh2);
    //TPZGeoMesh *gmesh = GMesh(ftriang, Lx,  Ly);//MalhaGeoT(h,hrefine);
    
    
	TPZVec<REAL> calcErro;
	for (int porder=1; porder<2; porder++) {
        
        erro<<"ordem "<<porder <<std::endl;
        for(int h=1;h<5;h++){
            erro << "NRef " << h << std::endl;
            
            erro<<std::endl;
            
            //1. Criacao da malha geom. e computacional
            bool hrefine=true;//true nao uniforme
            bool prefine=false;
            
            TPZGeoMesh *gmesh = GMesh(ftriang, Lx,  Ly);//MalhaGeoT(h,hrefine);
            UniformRefine(gmesh, h);
//            NoUniformRefine(gmesh, 1);
//            std::ofstream filemesh("MalhaGeoNaoUniT.vtk");
//            PrintGMeshVTK( gmesh, filemesh);
            

//            RefiningNearCircunference(2,gmesh,h,1);
//            std::ofstream filemesh2("MalhaGeoQArcTanRefeineNearCirc.vtk");
//            PrintGMeshVTK( gmesh, filemesh2);
            
            TPZCompMesh *cmesh = CompMeshPAdapJuan(*gmesh,porder,prefine);
            int nDofs;
            nDofs=cmesh->NEquations();
            
            erro<< "\nNRefinamento "<<h<< "   NDofs "<<nDofs<<std::endl;
            
            std::cout << "Neq = " << nDofs << std::endl;
            
            //2. Resolve o problema
            
            TPZAnalysis analysis(cmesh);
            int numthreads = 4;
            SolveSyst(analysis, cmesh, numthreads);
            // SolveLU ( analysis );
            
            
            
            //3. Calcula o erro
            
            //			TPZVec<REAL> calcErro;
            //            analysis.SetExact(*SolArcTan);
            //			analysis.PostProcess(calcErro,erro);
            
            ErrorHDiv(cmesh, erro);
            
            
            
            //4. visualizacao grafica usando vtk
            
            TPZVec<std::string> scalnames(5), vecnames(2);
            scalnames[0] = "Divergence";
            scalnames[1] = "ExactDiv";
            scalnames[2] = "Pressure";//"Solution";//
            scalnames[3] = "ExactPressure";
            scalnames[4] = "ExactDiv";
            
            
            vecnames[0] = "ExactFlux";
            vecnames[1] = "Flux";//"Derivative";//
            //
            //
            //
            //
            //			 //vecnames[0] = "Derivate";
            //               //  std::string plotfile("GraficArcTanHpAdaptivity.vtk");
            const int dim = 2;
            int div = 2;
            
            char buf[256] ;
            sprintf(buf,"ArcTanPhilUniMeshQ_porder%d_h%d.vtk",porder,h);
            analysis.DefineGraphMesh(dim,scalnames,vecnames,buf);
            analysis.PostProcess(div);
            
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                cmesh->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            
            //              std::string plotfile("GraficArcTanHpAdaptivity.vtk");
            //              const int dim = 2;
            //              int div = 2;
            //              analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
            //              analysis.PostProcess(div);
            //
            //        UniformRefine(gmesh, 1);
            //         std::ofstream filemesh2("MalhaGeoDepois.vtk");
            //        PrintGMeshVTK( gmesh, filemesh2);
            delete cmesh;
            delete gmesh;
            
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



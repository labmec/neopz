//
//  mainHpJan2015.cpp
//  PZ
//
//  Created by Denise de Siqueira on 1/16/15.
//
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"
#include "TPZVTKGeoMesh.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "mixedpoisson.h"
#include "PZMatPoissonControl.h"

#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <cmath>
#include <set>
#include "Tools.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;

const int MatId = 1;
const int dirichlet = 0;
const int neumann = 1;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;
void Forcing(const TPZVec<REAL> &pt, TPZVec<REAL> &res,TPZFMatrix<STATE> &disp);
TPZCompMesh *MeshH1(TPZGeoMesh *gmesh, int pOrder, int dim);
void SolutionControl(TPZAnalysis &an, std::string plotfile);
void EstadoAd(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
const int eps=100000;
const REAL Pi=M_PI;
void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads);



///Metodos para o ArcTangente


//solucoes exatas

bool Triang = false;
int main5(int argc, char *argv[])
{
    //InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    TPZVec<REAL> erros;
    ofstream saidaerros("../ErrosH1.txt",ios::app);
    saidaerros << "\nCalculo do Erro\n";
	
	
	
	
    for(int p=2; p<3; p++){
         
        if (Triang==true) {
              saidaerros << "\n ===========Malha Triang ===============\n";
        }

       
        
        for(int h=4; h<5; h++){
         //   saidaerros << "\nRefinamento: Ndiv = "<< h <<"\n";

            
            //------------ Etapa 1 ------------
            
           
			TPZGeoMesh * gmesh = GMesh(Triang,1,1);
			ofstream arg2("MalhaGeo.txt");
			gmesh->Print(arg2);
		
			UniformRefine(gmesh,h);
			std::ofstream filemesh1("MalhaGeoControle.vtk");
			TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh1, true);
            


			// REsolvendo o problema em H1
              TPZCompMesh *cmeshH1 = MeshH1(gmesh, p, 2);
				//Prefinamento(cmeshH1, h, pp);

                TPZAnalysis anh1(cmeshH1, true);
				
      
				int numthreads=1;
				SolveSyst(anh1, cmeshH1,numthreads);
            
                stringstream refh1,grauh1;
                grauh1 << p;
                refh1 << h;
                string strgh1 = grauh1.str();
                string strrh1 = refh1.str();
                std::string plotnameh1("SolutionControle");
                std::string Grauh1("P");
                std::string Refh1("h_");
                std::string VTKh1(".vtk");
                std::string plotDatah1;
                plotDatah1 = plotnameh1+Grauh1+strgh1+Refh1+strrh1+VTKh1;
                std::string plotfileh1(plotDatah1);
                
                SolutionControl(anh1, plotfileh1);
				anh1.SetExact(*EstadoAd);
				 saidaerros<< "\nNRefinamento "<<h<< "   NDofs "<<cmeshH1->NEquations()<<" N Comp. Elements " <<cmeshH1->NElements()<< std::endl;
				anh1.PostProcess(erros,saidaerros);
                
 
		}
    }
    
	return EXIT_SUCCESS;
}


void Forcing(const TPZVec<REAL> &pt, TPZVec<REAL> &res,TPZFMatrix<STATE> &disp){
	disp.Redim(2,1);
	double x = pt[0];
    double y = pt[1];
    res[0]= 0.;
}

void EstadoAd(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    u.Resize(1, 0.);
    du.Resize(3, 1.);
    du(0,0)=du(1,0)=du(2,0)=0.;

	const REAL alpha=0.001;
    
    const REAL sol = 10.*x*y*(1-x)*(1-y);
    u[0] = sol;
    

}
TPZCompMesh *MeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
{
	 /// criar materiais
    dim = 2;
    TPZMatPoissonControl *material = new TPZMatPoissonControl( MatId,  dim);
    material->NStateVariables();
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
	
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);

	    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(EstadoAd);
    material->SetForcingFunctionExact(solexata);
	
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    
    dum = new TPZDummyFunction<STATE>(Forcing);
    dum->SetPolynomialOrder(20);
    force = dum;
    material->SetForcingFunction(force);
	
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    return cmesh;
    
}

void SolutionControl(TPZAnalysis &an, std::string plotfile){
    
	TPZManVector<std::string,10> scalnames(2), vecnames(1);
    
	scalnames[0] = "Solution";
	scalnames[1]= "OrdemP";
	vecnames[0] = "Derivative";

	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
//	std::ofstream out("malhaH1.txt");
//	an.Print("nothing",out);
}
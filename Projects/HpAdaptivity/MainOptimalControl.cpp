//
//  mainHpJan2015.cpp
//  PZ
//
//  Created by Denise de Siqueira on 1/16/15.
//
//

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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
#include "ToolsOptimal.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;

const int MatId = 1;
const int dirichlet = 0;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;


void Forcing(const TPZVec<REAL> &pt, TPZVec<REAL> &res,TPZFMatrix<STATE> &disp);
TPZCompMesh *MeshH1(TPZGeoMesh *gmesh, int pOrder, int dim, bool hasbc);
TPZCompMesh *MeshL2(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *MalhaMultifisicaOpt(TPZVec<TPZCompMesh *> meshvec, TPZGeoMesh *gmesh);
void StateAd(const TPZVec<REAL>&pt,TPZVec<STATE> &res, TPZFMatrix<STATE> & disp);


void SolutionControl(TPZAnalysis &an, std::string plotfile);
void EstadoAd(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);


const int bcdirichlet = 0;



void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads);


bool Triang = false;
int main(/*int argc, char *argv[]*/)
{
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    
    int p =1;
    int ndiv = 1;
    TPZGeoMesh *gmesh = GMesh(false, 1, 1);
    UniformRefine( gmesh,ndiv);
    
    ofstream filegmesh1("gmes.txt");
    gmesh->Print(filegmesh1);
    
    
    TPZCompMesh *cmesh1 = MeshH1(gmesh, p, 2,true);//malha para y(estado)
        ofstream filemesh1("malhaY.txt");
        cmesh1->Print(filemesh1);
    TPZCompMesh *cmesh2 = MeshL2(gmesh, p, 2);//malha para controle(u)
    ofstream filemesh2("malhaU.txt");
    cmesh2->Print(filemesh2);
    TPZCompMesh *cmesh3 = MeshH1(gmesh, p, 2,false);//malha para mult.lagrange(p)
    ofstream filemesh3("malhaP.txt");
    cmesh3->Print(filemesh3);
    
    TPZManVector<TPZCompMesh *,3> meshvec(3);
    meshvec[0] = cmesh1;
    meshvec[1] = cmesh2;
    meshvec[2] = cmesh3;
    TPZCompMesh * mphysics = MalhaMultifisicaOpt(meshvec, gmesh);
    ofstream filemesh4("malhaMulti.txt");
    mphysics->Print(filemesh4);
    
    
    
    //Resolver problema
    TPZAnalysis analysis(mphysics,false);
    //TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(mphysics);
    //strmat.SetDecomposeType(ELDLt);
    TPZSkylineStructMatrix strmat(mphysics);
    
    //strmat.SetNumThreads(6);
    analysis.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); //caso simetrico
    analysis.SetSolver(step);
    analysis.Assemble();
    analysis.Solve();
    
    std::ofstream out("SolOpt.txt");
    analysis.Solution().Print("Solucao",out);

    
    
    //Post-Process
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    
    TPZManVector<std::string,10> scalnames(4), vecnames(0);
    scalnames[0] = "State";
    scalnames[1] = "Control";
    scalnames[2]="LagrangeMult";
    scalnames[3]="ExactState";
    
    std::stringstream name;
    name << "Solution_Opt" <<ndiv<< ".vtk";
    std::string paraviewfile(name.str());
    analysis.DefineGraphMesh(2,scalnames,vecnames,paraviewfile);
    analysis.PostProcess(0);
    
    //visualizar matriz no vtk
    //    TPZFMatrix<REAL> vismat(100,100);
    //    mphysics->ComputeFillIn(100,vismat);
    //    VisualMatrixVTK(vismat,"matrixstruct.vtk");
    
    
	return EXIT_SUCCESS;
}


void Forcing(const TPZVec<REAL> &pt, TPZVec<REAL> &res,TPZFMatrix<STATE> &disp){
	disp.Redim(2,1);
//	double x = pt[0];
  //  double y = pt[1];
    res[0]= 0.;
}

void EstadoAd(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    u.Resize(1, 0.);
    du.Resize(3, 1.);
    du(0,0)=du(1,0)=du(2,0)=0.;
    
	//const REAL alpha=0.001;
    
    const REAL sol = 10.*x*y*(1-x)*(1-y);
    u[0] = sol;
    
    
}

TPZCompMesh *MeshH1(TPZGeoMesh *gmesh, int pOrder, int dim,bool hasbc)
{
    /// criar materiais
    dim = 2;
    TPZMatPoisson3d *material = new TPZMatPoisson3d( MatId,  dim);
    material->NStateVariables();
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
	
    ///Inserir condicao de contorno
    if(hasbc){
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	
        TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	
        cmesh->InsertMaterialObject(BCond0);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
    }
    
    //solucao exata
//    TPZAutoPointer<TPZFunction<STATE> > solexata;
//    solexata = new TPZDummyFunction<STATE>(EstadoAd);
//    material->SetForcingFunctionExact(solexata);
//
//    //funcao do lado direito da equacao do problema
//    TPZAutoPointer<TPZFunction<STATE> > force;
//    TPZDummyFunction<STATE> *dum;
//    
//    dum = new TPZDummyFunction<STATE>(ForcingOpt);
//    dum->SetPolynomialOrder(20);
//    force = dum;
//    material->SetForcingFunction(force);
	
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    return cmesh;
    
}

TPZCompMesh *MeshL2(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    dim = 2;
    TPZMatPoisson3d *material = new TPZMatPoisson3d( MatId,  dim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
//    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//    
//    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
//    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
//    
//    cmesh->InsertMaterialObject(BCond0);
//    cmesh->InsertMaterialObject(BCond1);
//    cmesh->InsertMaterialObject(BCond2);
//    cmesh->InsertMaterialObject(BCond3);
    
//    //solucao exata
//    TPZAutoPointer<TPZFunction<STATE> > solexata;
//    solexata = new TPZDummyFunction<STATE>(EstadoAd);
//    material->SetForcingFunctionExact(solexata);
//    
//    //funcao do lado direito da equacao do problema
//    TPZAutoPointer<TPZFunction<STATE> > force;
//    TPZDummyFunction<STATE> *dum;
//    
//    dum = new TPZDummyFunction<STATE>(OptForcing);
//    dum->SetPolynomialOrder(20);
//    force = dum;
//    material->SetForcingFunction(force);
    
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    //cmesh->SetAllCreateFunctionsContinuous();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    return cmesh;
    
}

TPZCompMesh *MalhaMultifisicaOpt(TPZVec<TPZCompMesh *> meshvec, TPZGeoMesh *gmesh){
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =2;
    
    TPZMatPoissonControl *material = new TPZMatPoissonControl(MatId,dim);
    
    
    //incluindo os dados do problema
    REAL k=1;
    REAL alpha=1;
    material-> SetParameters( k, alpha);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(StateAd, 5);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(OptForcing, 5);
    dum->SetPolynomialOrder(20);
    force = dum;
    material->SetForcingFunction(force);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    mphysics->SetDimModel(dim);
    
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0, bcdirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1, bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2, bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3, bcdirichlet, val1, val2);
    
    ///Inserir condicoes de contorno
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
    
}

void StateAd(const TPZVec<REAL>&pt,TPZVec<STATE> &res, TPZFMatrix<STATE> & disp){
    disp.Redim(2,1);
    res[0]=0.;
    double x=pt[0];
    double y=pt[1];
    res[0]=10.*x*y*(1-x)*(1-y);
    
}
//
//void OptForcing(const TPZVec<REAL>&pt,TPZVec<REAL> &res, TPZFMatrix<STATE> & disp){
//    disp.Redim(2,1);
//    res[0]=0.;
//    res[1]=0.;
//    double x=pt[0];
//    double y=pt[1];
//    res[1]=10.*x*y*(1-x)*(1-y);
//    
//}


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

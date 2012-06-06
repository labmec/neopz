//
//  main_mixed.cpp
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>
using namespace std;

int const matId =1;
int const bc0=-1;
int const bc1=-2;
int const bc2=-3;
int const bc3=-4;

int const dirichlet =0;
int const neumann = 1;

REAL const Pi = 4.*atan(1.);

TPZGeoMesh *GMesh(bool triang_elements);
TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson* &mymaterial);
void Forcing1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
TPZCompMesh *CMeshHDivPressure(TPZGeoMesh *gmesh, int pOrder);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
#endif


int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	std::string logs("../logmixedproblem.cfg");
	InitializePZLOG("../logmixedproblem.cfg");
#endif
    
    int p =1;
	//primeira malha
	
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = GMesh(true);
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
        std::stringstream sout;
        sout<<"\n\n Malha Geometrica Inicial\n ";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
    
    // First computational mesh
	TPZCompMesh * cmesh1= CMeshFlux(gmesh,  p);
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
        std::stringstream sout;
        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
        cmesh1->Print(sout);
        LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
    
	// Second computational mesh
	TPZCompMesh * cmesh2 = CMeshPressure(gmesh, 0);
	
#ifdef DEBUG   
    int ncel = cmesh2->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh2->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
        std::stringstream sout;
        sout<<"\n\n Malha Computacional_2 pressure\n ";
        cmesh2->Print(sout);
        LOGPZ_DEBUG(logdata,sout.str());
	}
#endif
    
    //malha multifisica
    TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh1;
	meshvec[1] = cmesh2;
    TPZMixedPoisson * mymaterial;
    TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,mymaterial);
    //ofstream arg("mphysic.txt");
	//mphysics->Print(arg);
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
        std::stringstream sout;
        sout<<"\n\n Malha Computacional Multiphysic\n ";
        mphysics->Print(sout);
        LOGPZ_DEBUG(logdata,sout.str());
	}
#endif
    
//    TPZAnalysis an(mphysics);
//	SolveSyst(an, mphysics);
    
  //  std::string plotfile("saidaSolution_mphysics.vtk");
   // PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
    
    
    TPZCompMesh * cmesh3= CMeshHDivPressure(gmesh,  p);
    TPZAnalysis an2(cmesh3);
	SolveSyst(an2, cmesh3);
    
    
    
  
    return 0;
}

TPZGeoMesh *GMesh(bool triang_elements){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
    TPZVec <int> TopolTriang(3);
	TPZVec <int> TopolLine(2);
	
	//indice dos nos
	int id = 0;
	REAL valx, dx=1.;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = 1. - xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,1. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    if(triang_elements==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
	//ofstream arg("gmesh.txt");
    //	gmesh->Print(arg);
    
	return gmesh;
}

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	TPZMaterial * mat(material);
	material->NStateVariables();
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
	cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	material->NStateVariables();
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}


TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson * &mymaterial){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=1;
    int dim =2;
    
    REAL coefk=1.;
    mymaterial = new TPZMixedPoisson(MatId,dim);
    mymaterial->SetParameters(coefk);
    
    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
	mymaterial->SetForcingFunction(force1);
    
    //REAL fxy=1.;
   // mymaterial->SetInternalFlux(fxy);
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = mymaterial->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = mymaterial->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

void Forcing1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VTK 
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
		
	
	const int dim = 2;
	int div =0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
    
}

TPZCompMesh *CMeshHDivPressure(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
//	int dim = 2;
//	TPZMatPoisson3d *material;
//	material = new TPZMatPoisson3d(matId,dim); 
//	material->NStateVariables();
//	
//    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//    TPZMaterial * mat(material);
//    cmesh->InsertMaterialObject(mat);
//    
//    ///Inserir condicao de contorno
//	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
//	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
//    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
//    
//    
//	cmesh->SetAllCreateFunctionsHDivPressure();
//    cmesh->InsertMaterialObject(BCond0);
//    cmesh->InsertMaterialObject(BCond1);
//    cmesh->InsertMaterialObject(BCond2);
//    cmesh->InsertMaterialObject(BCond3);
//    
//	cmesh->SetDefaultOrder(pOrder);
//    cmesh->SetDimModel(dim);
//	
//	//Ajuste da estrutura de dados computacional
//	cmesh->AutoBuild();
    
    
    /// criar materiais
    int MatId=1;
    int dim =2;
	TPZMatPoisson3d *material = new TPZMatPoisson3d(MatId,dim);
	material->NStateVariables();
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    
    
    REAL diff = 1.;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
    material->SetParameters(diff, conv, convdir);
   
    TPZAutoPointer<TPZFunction<REAL> > force1 = new TPZDummyFunction<REAL>(Forcing1);
	material->SetForcingFunction(force1);
//    REAL fxy=1.;
//    material->SetInternalFlux(fxy);
//    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    cmesh->SetAllCreateFunctionsHDivPressure();
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    cmesh->AutoBuild();
	//cmesh->AdjustBoundaryElements();
	//cmesh->CleanUpUnconnectedNodes();	

	return cmesh;
}



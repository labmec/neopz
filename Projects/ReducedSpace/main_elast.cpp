//
//  File.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/31/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>


#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzreferredcompel.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"

#include "pzelasmat.h"
#include "pzmat1dlin.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "pzreducedspace.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzbuildmultiphysicsmesh.h"

#include <iostream>
#include <math.h>
using namespace std;

int const matId1 =1; //elastic
int const matId2 =2; //pressure
int const bc0=-1;
int const bc1=-2;
int const bc2=-3;
int const bc3=-4; //bc pressure


int const dirichlet =0;
int const neumann = 1;
int const mixed =2;

REAL const Pi = 4.*atan(1.);

TPZGeoMesh *GMesh(int nh,REAL w, REAL L);
TPZCompMesh *CMeshElastic(TPZGeoMesh *gmesh, int pOrder);
TPZCompMeshReferred *CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);

void MySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
void PosProcessamento1(TPZAnalysis &an, std::string plotfile);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.reducedspace.data"));
#endif


int main(int argc, char *argv[])
{
    int p = 1;
	//primeira malha
	
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = GMesh(0,1.,1.);
    ofstream arg1("gmesh_inicial.txt");
    gmesh->Print(arg1);
    
    //computational mesh elastic
	TPZCompMesh * cmesh_elast= CMeshElastic(gmesh, p);
    ofstream arg2("cmesh_inicial.txt");
    cmesh_elast->Print(arg2);
        
    TPZAnalysis an(cmesh_elast);
    MySolve(an, cmesh_elast);
    string plotfile("saidaSolution_mesh1.vtk");
    PosProcessamento1(an, plotfile);
    TPZFMatrix<REAL> solucao;
    solucao=cmesh_elast->Solution();
    solucao.Print();
    
    //computational mesh of reduced space
    TPZCompMeshReferred *cmesh_referred = CMeshReduced(gmesh, cmesh_elast, p);
    cmesh_referred->ComputeNodElCon();
    TPZFStructMatrix fstr(cmesh_referred);
    TPZFMatrix<STATE> rhs(1);
    TPZAutoPointer<TPZMatrix<STATE> > strmat = fstr.CreateAssemble(rhs,NULL);
    strmat->Print("rigidez");
    rhs.Print("forca");
    
    ofstream arg3("cmeshreferred_inicial.txt");
    cmesh_referred->Print(arg3);
	
    //computational mesh of pressure
    TPZCompMesh *cmesh_pressure = CMeshPressure(gmesh, p);
    ofstream arg4("cmeshpressure_inicial.txt");
    cmesh_pressure->Print(arg4);
    
    
    //------- computational mesh multiphysic ----------//
    
    // Cleaning reference of the geometric mesh to cmesh_referred
	gmesh->ResetReference();
	cmesh_referred->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh_referred,0);
	cmesh_referred->AdjustBoundaryElements();
	cmesh_referred->CleanUpUnconnectedNodes();
    ofstream arg5("cmeshreferred_final.txt");
    cmesh_referred->Print(arg5);
    
    // Cleaning reference of the geometric mesh to cmesh_pressure
	gmesh->ResetReference();
	cmesh_pressure->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh_pressure,0);
	cmesh_pressure->AdjustBoundaryElements();
	cmesh_pressure->CleanUpUnconnectedNodes();
    ofstream arg6("cmeshpressure_final.txt");
    cmesh_pressure->Print(arg6);
    
    //multiphysic mesh
    TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh_referred;
	meshvec[1] = cmesh_pressure;
    
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->AutoBuild();
	mphysics->CleanUpUnconnectedNodes();

	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    ofstream arg7("mphysics.txt");
    mphysics->Print(arg7);
    
    return 0;
}

TPZGeoMesh *GMesh(int nh,REAL w, REAL L){
    
    int Qnodes = 6;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
    TPZVec <int> TopolTriang(3);
	TPZVec <int> TopolLine(2);
    TPZVec <int> TopolPoint(1);
	
	//indice dos nos
	int id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*L/2.;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = L - xi*L/2.;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,w);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	//indice dos elementos
	id = 0;
    
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 4;
        TopolQuad[3] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId1,*gmesh);
        id++;
    
        TopolQuad[0] = 1;
        TopolQuad[1] = 2;
        TopolQuad[2] = 3;
        TopolQuad[3] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId1,*gmesh);
        id++;
    
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId2,*gmesh);
        id++;
    
        TopolPoint[0] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolLine,bc3,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 5;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
    
    
    
	gmesh->BuildConnectivity();
    
    //refinamento uniforme
    for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			//if (gel->Dimension() == 2) gel->Divide (filhos);
            gel->Divide (filhos);
		}//for i
	}//ref

	return gmesh;
}


TPZCompMesh *CMeshElastic(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	
    TPZVec<REAL> force(dim,0.);
    REAL E = 100.;
	REAL poisson = 0.35;
    int planestress = -1;
    
    TPZElasticityMaterial *material;
	material = new TPZElasticityMaterial(matId1, E, poisson, force[0], force[1], planestress);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    
    ///Inserir condicao de contorno
    REAL big = material->gBigNumber;
    REAL sign = -1.;
    
    TPZFMatrix<REAL> val11(2,2,0.), val21(2,1,0.);
    val21(0,0)=sign;
    val21(1,0)=sign;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc0,neumann, val11, val21);
    
    TPZFMatrix<REAL> val12(2,2,0.), val22(2,1,0.);
    val12(1,1) = big;
    TPZMaterial * BCond2 = material->CreateBC(mat, bc1,mixed, val12, val22);
    
    TPZFMatrix<REAL> val13(2,2,0.), val23(2,1,0.);
    val13(0,0) = big;
    TPZMaterial * BCond3 = material->CreateBC(mat, bc2,dirichlet, val13, val23);
    
    cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);
		
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
        
	return cmesh;
}

TPZCompMeshReferred *CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder){
    
    /// criar materiais
	int dim = 2;
    
    TPZVec<REAL> force(dim,0.);
    REAL E = 100.;
	REAL poisson = 0.35;
    int planestress = -1;
    
	TPZElasticityMaterial *material;
	material = new TPZElasticityMaterial(matId1, E, poisson, force[0], force[1], planestress); 
	material->NStateVariables();
    
    TPZVec<REAL> convdir;
    convdir.Resize(dim,0.);
    
    TPZCompMeshReferred *cmeshreferred = new TPZCompMeshReferred(gmesh);
    
    cmeshreferred->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmeshreferred->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    REAL big = material->gBigNumber;
    REAL sign = -1.;
    
    TPZFMatrix<REAL> val11(2,2,0.), val21(2,1,0.);
    val21(0,0)=sign;
    val21(1,0)=sign;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc0,neumann, val11, val21);
    
    TPZFMatrix<REAL> val12(2,2,0.), val22(2,1,0.);
    val12(1,1) = big;
    TPZMaterial * BCond2 = material->CreateBC(mat, bc1,mixed, val12, val22);
    
    TPZFMatrix<REAL> val13(2,2,0.), val23(2,1,0.);
    val13(0,0) = big;
    TPZMaterial * BCond3 = material->CreateBC(mat, bc2,dirichlet, val13, val23);

    
    int numsol = cmesh->Solution().Cols();
    cmeshreferred->AllocateNewConnect(numsol, 1, 1);
    
	TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmeshreferred);
    
    cmeshreferred->InsertMaterialObject(BCond1);
    cmeshreferred->InsertMaterialObject(BCond2);
    cmeshreferred->InsertMaterialObject(BCond3);
    
	cmeshreferred->SetDefaultOrder(pOrder);
    cmeshreferred->SetDimModel(dim);
	
    gmesh->ResetReference();
	//Ajuste da estrutura de dados computacional
	cmeshreferred->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    cmeshreferred->LoadReferred(cmesh);
   
    return cmeshreferred;
}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder){
    
    /// criar materiais
	int dim = 1;
	
    TPZMat1dLin *material;
	material = new TPZMat1dLin(matId2);
    
    TPZFMatrix<REAL> xk(1,1,1.);
    TPZFMatrix<REAL> xc(1,1,0.);
    TPZFMatrix<REAL> xb(1,1,0.);
    TPZFMatrix<REAL> xf(1,1,0.);
    material->SetMaterial(xk,xc,xb,xf);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    REAL vazao = 10.;
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(0,0)=vazao;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc3,neumann, val1, val2);
    
    
    cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
    cmesh->InsertMaterialObject(BCond1);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
    
}

void MySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void PosProcessamento1(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(1), vecnames(0);
	scalnames[0] = "SigmaX";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

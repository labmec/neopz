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
#include "pzreducedspace.h"
#include "pzfstrmatrix.h"
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

TPZGeoMesh *GMesh(int triang_elements, int nh);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);
TPZCompMeshReferred *CMeshReferred(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder);
void Forcingf(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);

void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
void PosProcessamento(TPZAnalysis &an, std::string plotfile);

#ifdef PZ_LOG
static TPZLogger logdata("pz.mixedpoisson.data");
#endif


int main3(int argc, char *argv[])
{
    int p = 2;
	//primeira malha
	gRefDBase.InitializeAllUniformRefPatterns();
	
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = GMesh(-1,2);
    ofstream arg1("gmesh_inicial.txt");
    gmesh->Print(arg1);
    
    // First computational mesh
	TPZCompMesh * cmesh= CMesh(gmesh, p);
    ofstream arg2("cmesh_inicial.txt");
    cmesh->Print(arg2);
        
    TPZAnalysis an(cmesh);
    mySolve(an, cmesh);
    string plotfile("saidaSolution_mesh1.vtk");
    PosProcessamento(an, plotfile);

    TPZFMatrix<STATE> solucao;
    solucao=cmesh->Solution();
    solucao.Print("solucao");
    
    TPZCompMeshReferred *cmeshreferred = CMeshReferred(gmesh, cmesh, p);
    cmeshreferred->ComputeNodElCon();
    TPZFStructMatrix fstr(cmeshreferred);
    TPZFMatrix<STATE> rhs(1);
    TPZAutoPointer<TPZMatrix<STATE> > strmat = fstr.CreateAssemble(rhs,NULL);
    
    strmat->Print("rigidez");
    rhs.Print("forca");
    ofstream arg3("cmeshreferred_inicial.txt");
    cmeshreferred->Print(arg3);
	    
    
    return 0;
}

TPZGeoMesh *GMesh(int triang_elements, int nh){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
	int64_t id = 0;
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
    
    if(triang_elements==1)
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
    
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Geometrica Inicial\n ";
    //        gmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int64_t n = gmesh->NElements();
		for ( int64_t i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			//if (gel->Dimension() == 2) gel->Divide (filhos);
            gel->Divide (filhos);
		}//for i
	}//ref

	return gmesh;
}


TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	material->NStateVariables();
    
    REAL diff=1.;
    REAL conv =0.;
    TPZVec<REAL> convdir;
    convdir.Resize(dim,0.);
    material-> SetParameters(diff, conv, convdir);
    
    REAL ff=-8.;
    material->SetInternalFlux(ff);
    
//    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcingf);
//    material->SetForcingFunction(force1);
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    
	cmesh->SetAllCreateFunctionsContinuous();
    
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
        
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_2 pressure\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
	
	return cmesh;
}

TPZCompMeshReferred *CMeshReferred(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder){
    
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	material->NStateVariables();
    
    REAL diff=1.;
    REAL conv =0.;
    TPZVec<REAL> convdir;
    convdir.Resize(dim,0.);
    material-> SetParameters(diff, conv, convdir);
    REAL ff=-8.;
    material->SetInternalFlux(ff);
    
    //TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcingf);
    //material->SetForcingFunction(force1);

    
    TPZCompMeshReferred *cmeshreferred = new TPZCompMeshReferred(gmesh);
    
    cmeshreferred->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmeshreferred->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    int64_t numsol = cmesh->Solution().Cols();
    cmeshreferred->AllocateNewConnect(numsol, 1, 1);
    
	TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmeshreferred);
    
    cmeshreferred->InsertMaterialObject(BCond0);
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


void Forcingf(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
    TPZBaseMatrix &anSol = an.Solution();
	anSol.Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void PosProcessamento(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(1), vecnames(0);
	scalnames[0] = "Solution";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

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

#include "pzpoisson3d.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzlog.h"

#include <iostream>
using namespace std;

int const matId =1;
int const bc0=-1;
int const bc1=-2;
int const bc2=-3;
int const bc3=-4;

int const dirichlet =0;
int const neumann = 1;

TPZGeoMesh *GMesh(bool triang_elements);
TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);


int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	std::string logs("../log_mixedproblem.cfg");
	InitializePZLOG();
#endif
    
    int p =1;
	//primeira malha
	
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = GMesh(true);
	ofstream arg1("gmeshZero.txt");
	gmesh->Print(arg1);
    
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
    id++;
    
    if(triang_elements==true)
    {
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolQuad[0] = 2;
        TopolQuad[1] = 1;
        TopolQuad[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
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
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsHDiv();
	cmesh->InsertMaterialObject(mat);
	
    
	///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond0);
    
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond1);
    
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond2);
    
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond3);
	
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
	TPZMaterial * mat(material);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->InsertMaterialObject(mat);
	
    
	///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond0);
    
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond1);
    
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond2);
    
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond3);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}



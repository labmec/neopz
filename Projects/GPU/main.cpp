/*{{{ includes*/
#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZRefPatternTools.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPattern.h"
#include "pzgraphmesh.h"

#include "tpzchangeel.h"
#include "TPZGeoElement.h"
#include "pzreftriangle.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPattern.h"
#include "TPZMatLaplacian.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"

#include "GPUMatrix.h"

using namespace pzgeom;
/*}}}*/

/*{{{ declaration*/
TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh);
/*}}}*/

int main(int argc, char *argv[]){/*{{{*/

	TPZGeoMesh	*gmesh = CreateGeoMesh();
    //std::ofstream file("/Users/santos/Desktop/CUDAgeomesh.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file);
    
    TPZCompMesh *cmesh = CreateCompMesh(gmesh);

    int64_t n = cmesh->NElements();
    for(int64_t i = 0; i < n; i++){
        
        TPZCompEl *compel = cmesh->Element(i);
        if(!compel) continue;
        
        TPZGeoEl *geoel = compel->Reference();
        if(!geoel) DebugStop();
        
        TPZVec<int64_t> SonsIndex;
        cmesh->Divide(i, SonsIndex, false);
    }
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->Reference()->BuildConnectivity();
    
    std::ofstream file("/Users/santos/Desktop/CUDAcompmesh.txt");
    cmesh->Print(file);
    
    
    TPZSkylineStructMatrix<STATE> skylstruct(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    TPZAnalysis an(cmesh);
    an.SetStructuralMatrix(skylstruct);
    an.SetSolver(step);
    an.Run();
    
    //post process
    int dimension = 2, resolution = 0;
    std::string plotfile("/Users/santos/Desktop/CUDAsolution.vtk");
    TPZVec <std::string> scalnames(1), vecnames(0);
    scalnames[0] = "Pressure";
    
    an.DefineGraphMesh(dimension, scalnames, vecnames, plotfile);
    an.PostProcess(resolution);

	return 0;
}
/*}}}*/

TPZGeoMesh *CreateGeoMesh(){/*{{{*/

    int nnodes 			= 4;
    int64_t id				= 0;
    int mat 			= 1;
	int reftype 		= 0;
	TPZGeoMesh *gmesh   = NULL;
	TPZManVector<REAL,3> coord(3,0);
    TPZManVector<int64_t,4> quad(4,0);
	TPZManVector<int64_t,2> boundary(2,0);

	gmesh = new TPZGeoMesh();
    gmesh->NodeVec().Resize( nnodes );

	//vertex 0
	coord[0]=0.; coord[1]=0.;
	gmesh->NodeVec()[0].SetCoord(coord);
	gmesh->NodeVec()[0].SetNodeId(0);
	//vertex 1
	coord[0]=10.; coord[1]=0.;
	gmesh->NodeVec()[1].SetCoord(coord);
	gmesh->NodeVec()[1].SetNodeId(1);
	//vertex 2
	coord[0]=10.; coord[1]=10.;
	gmesh->NodeVec()[2].SetCoord(coord);
	gmesh->NodeVec()[2].SetNodeId(2);
	//vertex 3
	coord[0]=0.; coord[1]=10.;
	gmesh->NodeVec()[3].SetCoord(coord);
	gmesh->NodeVec()[3].SetNodeId(3);

	//element quad 0
    quad[0]=0; quad[1]=1; quad[2]=2; quad[3]=3;
    gmesh->CreateGeoElement(EQuadrilateral,quad,mat,id,reftype);
    gmesh->ElementVec()[id]->SetId(0);

	//element line
	for(int64_t i=1;i<5;i++){
		boundary[0]=i-1; boundary[1]=(i<4)?i:0;
		gmesh->CreateGeoElement(EOned,boundary,mat+i,id,reftype);
		gmesh->ElementVec()[id]->SetId(i);	
	}

	gmesh->BuildConnectivity();

	return gmesh;
}
/*}}}*/

TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh){/*{{{*/
	
	const int dim                   = 2;
	const int pOrder                = 1;
	const int mat                   = 1;
    const int dirichlet             = 0;
    const int matface_in            = 5;
    const int matface_out           = 3;
	TPZMatLaplacian* matlaplacian 	= NULL;
	TPZCompMesh* cmesh              = NULL;

    //Creating material
    matlaplacian = new TPZMatLaplacian(mat);
    matlaplacian->SetDimension(dim);
    matlaplacian->SetParameters(1,0);
    matlaplacian->SetNoPenalty();
    matlaplacian->SetSymmetric();

    //Creating compmesh
	cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
    
    //Creating boundary conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond1 = matlaplacian->CreateBC(matlaplacian,matface_out,dirichlet, val1, val2);
    
    val2(0,0) = 1.0*1000.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = matlaplacian->CreateBC(matlaplacian,matface_in,dirichlet, val1, val2);

    //Inserting materials
    cmesh->InsertMaterialObject(matlaplacian);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    
	/// adjust the data structure
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();

	return cmesh;
}
/*}}}*/




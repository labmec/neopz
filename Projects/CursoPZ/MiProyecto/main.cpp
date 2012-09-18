
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"

#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

std::string Archivo = PZSOURCEDIR;

TPZGeoMesh *CreateGeoMesh(std::string &nome);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);

// bi-dimensional problem for elasticity
int main() {

	Archivo += "/Projects/CursoPZ/MiProyecto/";
	Archivo += "MiPlaca.dump";

	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMesh(Archivo);
	UniformRefine(gmesh,1);

	// Creating computational mesh (approximation space and materials)
	int p = 4;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMesh(gmesh);
	
	// Solving linear equations
	// Initial steps
	TPZAnalysis an (cmesh);
	TPZSkylineStructMatrix strskyl(cmesh);
	an.SetStructuralMatrix(strskyl);
	// Solver (is your choose) 
	TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
	direct->SetDirect(ECholesky);
	an.SetSolver(*direct);
	delete direct;
	direct = 0;
/*
	// Caso no simetrico
//	TPZFStructMatrix full(cmesh);
	TPZBandStructMatrix full(cmesh);
	an.SetStructuralMatrix(full);
	an.Solution().Zero();
	TPZStepSolver<REAL> step;
	step.SetDirect(ELU);
	an.SetSolver(step);
	*/
	an.Run();
	
	// Post processing
	TPZManVector<std::string> scalarnames(5), vecnames(3);
	scalarnames[0] = "SigmaX";
	scalarnames[1] = "SigmaY";
	scalarnames[2] = "Pressure";
	scalarnames[3] = "MaxStress";
	scalarnames[4] = "TauXY";
	vecnames[0] = "displacement";
	vecnames[1] = "PrincipalStress1";
	vecnames[2] = "PrincipalStress2";
	an.DefineGraphMesh(2,scalarnames,vecnames,"ElasticitySolutions.vtk");

	an.PostProcess(2);
}

//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh(std::string &archivo) {

	// Ejemplo uni-dimensional para la generacion de una malla para un reservatorio 
	TPZReadGIDGrid grid;
	TPZGeoMesh *meshgrid = grid.GeometricGIDMesh(archivo);
	if(!meshgrid->NElements())
		return 0;

	return meshgrid;
}
//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();

    // Creating elasticity material
    TPZMaterial * mat = new TPZElasticityMaterial(4,2000.,0.3,0.,0.);
    cmesh->InsertMaterialObject(mat);

	// Creating four boundary condition
    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
	TPZMaterial *bcBottom, *bcRigth, *bcTop, *bcLeft;

	// Condicion livre - nada para hacer
    bcLeft = mat->CreateBC(mat,5,3,val1,val2);
    cmesh->InsertMaterialObject(bcLeft);
	// Condicion de Dirichlet fijando la posicion de la placa
	val1(1,1) = 1000000.;
    bcBottom = mat->CreateBC(mat,1,0,val1,val2);
	cmesh->InsertMaterialObject(bcBottom);
	val1(1,1) = 0.;
	val2(1,0) = 10.;
    bcTop = mat->CreateBC(mat,2,1,val1,val2);
	cmesh->InsertMaterialObject(bcTop);
	val2(1,0) = 0.;
    bcRigth = mat->CreateBC(mat,3,1,val1,val2);
	cmesh->InsertMaterialObject(bcRigth);

	// Inserting boundary conditions into computational mesh
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
// #ifdef LOG4CXX
//     if (logger->isDebugEnabled())
//     {
//         std::stringstream sout;
//         cmesh->Print(sout);
//         LOGPZ_DEBUG(logger, sout.str())
//     }
// #endif
    return cmesh;
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {    
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
}

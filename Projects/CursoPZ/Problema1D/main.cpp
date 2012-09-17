
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzmatrix.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "MiViga1D.h"

#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

std::string Archivo = PZSOURCEDIR;

TPZGeoMesh *CreateGeoMesh(std::string &nome);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);

// uni-dimensional problem for elasticity
int main() {

	Archivo += "/Projects/CursoPZ/Problema1D/";
	Archivo += "Viga1D.dump";

	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMesh(Archivo);
//	UniformRefine(gmesh,1);

	// Creating computational mesh (approximation space and materials)
	int p = 3;
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
	TPZManVector<std::string> scalarnames(0), vecnames(2);
	vecnames[0] = "desplazamiento";
	vecnames[1] = "giro";
	an.DefineGraphMesh(1,scalarnames,vecnames,"MiSolucion1D.vtk");

	an.PostProcess(2);
}

//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh(std::string &archivo) {

	// Ejemplo uni-dimensional para la generacion de una malla para un reservatorio 
	TPZReadGIDGrid grid;
	TPZGeoMesh *meshgrid = grid.GeometricGIDMesh(archivo);

	int nelem = meshgrid->NElements();
	if(!nelem)
		return 0;

	// Buscando 

	// Criando los elementos puntuales en los extremos
	TPZGeoEl *gel = meshgrid->ElementVec()[0];
	TPZGeoElBC gbc1(gel,0,3);
	gel = meshgrid->ElementVec()[nelem-1];
	TPZGeoElBC gbc2(gel,1,2);

	// Re construindo la conectividad en la malha geometrica
	meshgrid->BuildConnectivity();

	return meshgrid;
}
//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());

    // Creating elasticity material
    TPZMaterial * mat = new TPZMiViga1D(1,2000.,0.3,0.4);
    cmesh->InsertMaterialObject(mat);
	((TPZMiViga1D *)mat)->SetLoad(0.);

	// Creating four boundary condition
    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
	val2(0,0) = 0.;
	val2(1,0) = 1.0;
	TPZMaterial *bcRight, *bcLeft;
	// Condicion Dirichlet desplazamiento y giro en el extremo inicial de la viga (v,phi) = (0,0)
    bcLeft = mat->CreateBC(mat,2,0,val1,val2);
    cmesh->InsertMaterialObject(bcLeft);
	// Condicion de Neumann colocando un peso en el extremo
	val2(0,0) = 1.;
	val2(1,0) = 0.2;
    //bcRight = mat->CreateBC(mat,3,1,val1,val2);
	//cmesh->InsertMaterialObject(bcRight);

	// Inserting boundary conditions into computational mesh
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
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
//	gmesh->BuildConnectivity();
}

/*
 * @file
 * @brief Contains the main function for solving a 2D linear partial differential equation.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzmat2dlin.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"

#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#include "PZ_Process.h"

#include <iostream>
#include <string>
#include <math.h>

using namespace std;

//
//	This program solve  a system of linear PDE - 1D I hope! 
//	- k.[(ddu/dxdx) + (ddu/dydy)] = f(x),   0<= x <=1    --->  Now k=1 and b=1 and  f(x)=x
//	using uDirichlet = u(x) = 0  on x = 0  and  x = 1.
//

/* RELATED TO MODEL */

// Right handside term of our Linear PDE
void ForcingFunction(const TPZVec<REAL> &pto, TPZVec<REAL> &xfloat) 
{
	int n = xfloat.NElements();
	for(int i=0;i<n;i++)
		xfloat[i] = 0.0;
}

// Exact Solution u(x)
void SolExata(TPZVec<REAL> &pto, TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact)
{
}

int main_2D()
//int main()
{
	// Initializing generation log comments as log4cxx done
	std::string logs("log4cxx.doubleprojection1d");
	InitializePZLOG();
	
	// Files for information (geometric mesh and computational mesh)
	std::ofstream infoGeoMesh("gmesh.txt");
	std::ofstream infoCompMesh("cmesh.txt");	
	// Files for output solution as multiplier coefficients
	std::ofstream FileSolution("Solution.txt");	
	std::ofstream FileError("erros.txt");

	// File for Mathematica output format
	std::ofstream outMath("Mathematica.txt");

	// p -> interpolation order
	int p = 21;
	// h -> level of uniform refinement of the initial mesh
	int h = 3;
	
	// Creating main extremes and material for current project
	TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the mesh.
	x1[2] = 0.0;
	// id of material and boundary conditions
	TPZVec<int> matId(1,1);
	TPZVec<int> bc(2,-1);
	bc[1] = -2;
	// type of boundary conditions
	TPZVec<int> bcType(2,0);   // 0 - dirichlet
	
	// Coefficients of the differential equation (Linear 2D). They will be required by material constructor
	TPZMat2dLin *material = new TPZMat2dLin(matId[0]);
	TPZFMatrix<REAL> xkin (1,1,1.0);
	TPZFMatrix<REAL> xcin (1,1,0.0);
	TPZFMatrix<REAL> xfin (1,1,0.0);
	material->SetMaterial(xkin,xcin,xfin);
	// inserting function force
	//material->SetForcingFunction(new TPZDummyFunction(ForcingFunction));
	
	// FEM PROCESS
	// Creating a geometric mesh and printing geometric information. Note: The coordinates of the nodes are 3D
	TPZGeoMesh *gmesh = GeomMesh2D(h,matId,bc,x0,x1,1);   // 1 -> triangular   0 -> quadrilateral
	gmesh->Print(infoGeoMesh);	

	// Creating a computational mesh and printing computational information.
	TPZCompMesh * cmesh = CompMesh(gmesh,p,material,bc,bcType);
	cmesh->Print(infoCompMesh);
	
	// Assembling and Solving linear system
	TPZAnalysis an(cmesh);
	SolveSist(an,cmesh);
	
	// Print Solution as multiplier coefficients
	TPZFMatrix<REAL> toprint = an.Solution();
	toprint.Print("Solution", FileSolution);
	
	// Solution output for Mathematica
	OutputMathematica(outMath,1,8,cmesh);

	// Plot erro (norms) 
	//an.SetExact(SolExata);
	//TPZVec<REAL> posproc;
	//FileError << "\n\n Solving with P = 1\n";
	//an.PostProcess(posproc,FileError); // Compute the errors
	//FileError<<endl;
	
	// Close output files
	infoGeoMesh.close();
	infoCompMesh.close();
	FileSolution.close();
	FileError.close();
	outMath.close();
	
	return EXIT_SUCCESS;
}

#include "TPZExtendGridDimension.h"
#include "pzmattest3d.h"

//int main3D()
int main()
{
	// Initializing generation log comments as log4cxx done
	std::string logs("log4cxx.doubleprojection1d");
	InitializePZLOG();
	
	// Files for information (geometric mesh and computational mesh)
	std::ofstream infoGeoMesh("gmesh.txt");
	std::ofstream infoCompMesh("cmesh.txt");	
	// Files for output solution as multiplier coefficients
	std::ofstream FileSolution("Solution.txt");	
	std::ofstream FileError("erros.txt");
	
	// File for Mathematica output format
	std::ofstream outMath("Mathematica.txt");
	
	// p -> interpolation order
	int p = 2;
	// h -> level of uniform refinement of the initial mesh
	int h = 3;
	
	// Creating main extremes and material for current project
	TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the mesh.
	x1[2] = 0.0;
	// id of material and boundary conditions
	TPZVec<int> matId(1,1);
	TPZVec<int> bc(2,-1);
	bc[1] = -2;
	// type of boundary conditions
	TPZVec<int> bcType(2,0);   // 0 - dirichlet
	
	// Coefficients of the differential equation (Linear 2D). They will be required by material constructor
	TPZMaterialTest3D *material = new TPZMaterialTest3D(matId[0]);
	TPZFMatrix<REAL> xfin (1,1,0.0);
	material->SetMaterial(xfin);
	// inserting function force
	//material->SetForcingFunction(new TPZDummyFunction(ForcingFunction));
	
	// FEM PROCESS
	// Creating a geometric mesh and printing geometric information. Note: The coordinates of the nodes are 3D
	TPZGeoMesh *gmesh2d = GeomMesh2D(h,matId,bc,x0,x1,1);

	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	TPZAutoPointer<TPZGeoMesh> gmesh = gmesh2d;
	TPZExtendGridDimension gmeshextend(gmesh,0.5);
	
	TPZGeoMesh *gmesh3D = gmeshextend.ExtendedMesh(2,-1,-1);
	gmesh3D->Print(infoGeoMesh);
	
	// Creating a computational mesh and printing computational information.
	TPZCompMesh *cmesh = CompMesh(gmesh3D,p,material,bc,bcType);
	cmesh->Print(infoCompMesh);
	
	// Assembling and Solving linear system
	TPZAnalysis an(cmesh);
	SolveSist(an,cmesh);
	
	return 0;
}

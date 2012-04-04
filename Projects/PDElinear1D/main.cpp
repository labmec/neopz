#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzmat1dlin.h"

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
//	- k.u''(x) + b.u'(x) = f(x),   0<= x <=1    --->  Now k=1 and b=1 and  f(x)=x
//	using uDirichlet = u(x) = 0  on  x = 0  and  x = 1.

/* RELATED TO MODEL */

// Right handside term of our Linear PDE
void ForcingFunction(const TPZVec<REAL> &pto, TPZVec<REAL> &xfloat) 
{
	double x=pto[0];
	xfloat[0] = x;
}

// Exact Solution u(x)
void SolExata(TPZVec<REAL> &pto, TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact)
{
	double x=pto[0];
	//	double y=pto[1];
	double xexact;
	xexact = 1.-sinh(x)/sinh(1.);
	u_exact[0] = xexact;
	//	u_exact[1] =0.;
	du_exact(0,0) = 1.-cosh(x)/sinh(1.);//dx
	//	du_exact(1,0) =0.;//dy
}

int main()
//int main_firs(int argc, char *argv[])
{
#ifdef LOG4CXX
	// Initializing generation log comments as log4cxx done
	std::string logs("log4cxx.doubleprojection1d");
	InitializePZLOG();
#endif

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
	int h = 4;
	
	// Creating main extremes and material for current project
	TPZManVector<REAL> x0(3,0.), x1(3,0.);  // Corners of the mesh. Coordinates are zeros.
	x1[0] = 1.0;
	// id of material and boundary conditions
	TPZVec<int> matId(1,1);
	TPZVec<int> bc(2,-1);
	bc[1] = -2;
	// type of boundary conditions
	TPZVec<int> bcType(2,0);   // 0 - dirichlet
	
	// Coefficients of the differential equation (Linear 1D). They will be required by material constructor
	TPZFMatrix<REAL> xkin(1,1,1.0);   //  k = 1
	TPZFMatrix<REAL> xbin(1,1,1.0);   //  b = 1
	TPZFMatrix<REAL> xcin(1,1,0.0);   //  c = 0
	TPZFMatrix<REAL> xfin(1,1,0.0);
	// Material data
	TPZMat1dLin *material;
	material = new TPZMat1dLin(matId[0]); 
	material->SetMaterial(xkin,xcin,xbin,xfin);
	// inserting function force
	material->SetForcingFunction(new TPZDummyFunction(ForcingFunction));
	
	// FEM PROCESS
	// Creating a geometric mesh and printing geometric information. Note: The coordinates of the nodes are 3D
	TPZGeoMesh *gmesh = GeomMesh(h,matId,bc,x0,x1);
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
	OutputMathematica(outMath,1,10, cmesh);

	// Plot erro (norms) 
	an.SetExact(SolExata);
	TPZVec<REAL> posproc;
	FileError << "\n\n Solving with P = 1\n";
	an.PostProcess(posproc,FileError); // Compute the errors
	FileError<<endl;
	
	// Close output files
	infoGeoMesh.close();
	infoCompMesh.close();
	FileSolution.close();
	FileError.close();
	outMath.close();

	return EXIT_SUCCESS;
}

//int main(int argc, char *argv[])
int main_refinement(int argc, char *argv[])
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
	// File for VTK visualization
	std::string outVTK("saidaT.vtk");

	// Creating nodes (extremes) and material for current project
	TPZManVector<REAL> x0(3,0.), x1(3,0.);  // Corners of the mesh. Coordinates are zeros.
	x1[0] = 1.0;
	// id of material and boundary conditions
	TPZVec<int> matId(1,1);
	TPZVec<int> bc(2,-1);
	bc[1] = -2;
	// type of boundary conditions
	TPZVec<int> bcType(2,0);   // 0 - dirichlet

	// Coefficients of the differential equation (Linear 1D). They will be required by material constructor
	TPZFMatrix<REAL> xkin(1,1,1.0);   //  k = 1
	TPZFMatrix<REAL> xbin(1,1,1.0);   //  b = 1
	TPZFMatrix<REAL> xcin(1,1,0.0);   //  c = 0
	TPZFMatrix<REAL> xfin(1,1,0.0);
	// Material data
	TPZMat1dLin *material;
	material = new TPZMat1dLin(matId[0]); 
	material->SetMaterial(xkin,xcin,xbin,xfin);
	// inserting function force
	material->SetForcingFunction(new TPZDummyFunction(ForcingFunction));
	
	// Loop over p -> interpolation order
	for(int p = 1; p < 7; p++)
	{
		// Loop over h -> level of uniform refinement
		for(int h = 2; h < 10; h+=2)
		{
			// Creating a geometric mesh and printing geometric information.  Note: The coordinates of the nodes are 3D
			TPZGeoMesh *gmesh = GeomMesh(h,matId,bc,x0,x1);
			gmesh->Print(infoGeoMesh);

			// Creating a geometric mesh and printing computational information.
			TPZCompMesh *cmesh = CompMesh(gmesh,p,material,bc,bcType);
			cmesh->Print(infoCompMesh);

			// Assembling and Solving linear system
			TPZAnalysis an(cmesh);
			SolveSist(an,cmesh);

			// Print Solution
			TPZFMatrix<REAL> toprint = an.Solution();
			toprint.Print("Solution", FileSolution);

			// Solution output for Mathematica
			OutputMathematica(outMath,1,5,cmesh);

			// Computing the errors (and plot)
			an.SetExact(SolExata);
			FileError << "\n\n Solving with P = " << p << " E H = " << pow(2.,h) << "\n";
			TPZVec<REAL> posproc;
			an.PostProcess(posproc,FileError); // Computing the errors.(pzanalysis.cpp)
			FileError << endl;

			// Solution output for VTK
		//	OutputVTK(outVTK,cmesh,an);
		}
		outMath << endl;
	}
	
	// Close output files
	infoGeoMesh.close();
	infoCompMesh.close();
	FileSolution.close();
	FileError.close();
	outMath.close();
	
	return EXIT_SUCCESS;
}


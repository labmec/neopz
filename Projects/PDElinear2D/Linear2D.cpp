/*
 * @file
 * @brief Contains the main function for solving a 2D linear partial differential equation.
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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

#ifndef USING_BOOST
#define USING_BOOST
#endif

//
//	This program solve  a system of linear PDE - 1D 
//	\f$ - k.[(ddu/dxdx) + (ddu/dydy)] = f(x)\f$,   \f$ 0<= x <=1\f$    --->  Now \f$k=1\f$ and \f$b=1\f$ and  \f$f(x)=x\f$
//	using uDirichlet \f$= u(x) = 0\f$  on \f$x = 0\f$  and  \f$x = 1\f$.
// WEAK FORMULATION
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

//int main_2D()
int main()
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
	int pp = 17, p;
	// h -> level of uniform refinement of the initial mesh
	int h = 2;
	
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
	TPZGeoMesh *gmesh = GeomMesh2D(h,matId,bc,x0,x1,0);   // 1 -> triangular   0 -> quadrilateral
	gmesh->Print(infoGeoMesh);	

for(p = 10;p < pp;p++) {
	// Creating a computational mesh and printing computational information.
	TPZCompMesh * cmesh = CompMesh(gmesh,p,material,bc,bcType);
	if(!p) cmesh->Print(infoCompMesh);
	
	// Assembling and Solving linear system
	TPZAnalysis an(cmesh);
	SolveSist(an,cmesh);
	
	// Print Solution as multiplier coefficients
	TPZFMatrix<REAL> toprint = an.Solution();
	if(!p) toprint.Print("Solution", FileSolution);
	
	// Plot erro (norms) 
	//an.SetExact(SolExata);
	//TPZVec<REAL> posproc;
	//FileError << "\n\n Solving with P = 1\n";
	//an.PostProcess(posproc,FileError); // Compute the errors
	//FileError<<endl;
	
	// Solution output for Mathematica
	OutputMathematica(outMath,1,8,cmesh);
		
	// Solution output for Paraview (VTK)
	std::string vtkString("SolVTK.vtk");
	OutputVTK(vtkString,cmesh,an);
}	
	// Close output files
	infoGeoMesh.close();
	infoCompMesh.close();
	FileSolution.close();
	FileError.close();
	outMath.close();

	return EXIT_SUCCESS;
}

// Testing for 3D problem 

#include "TPZExtendGridDimension.h"
#include "pzmattest3d.h"

int main3D()
//int main()
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
	int p = 14;
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
	
	// Testing
	
	
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

/**
RESULTADOS

Using triangles with h = 1 and p = 1, 2, 3, ... 20
we have 
Solution (right accuracy) to p = 2, 3, ..., 14, 16 and 17
Poor solution to p = 1
Solution with big oscillation (but between limits) p = 15
Solution out of range to p = 18, 19 and 20.

Using quadrilaterals with h = 2 and p = 1, 2, 3, ..., 20
IT IS IMPORTANT STORE THE FOLLOWING INFORMATION TO PROJECT RUN TIME OF THE ANY TEST:
It runs p = 1, 2, ..., 16
 break in p = 17
The time of runs were:
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:00.037513 Time for solving 00:00:00.000219
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:00.118877 Time for solving 00:00:00.003716
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:00.348824 Time for solving 00:00:00.049339
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:01.576123 Time for solving 00:00:00.174542
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:01.515140 Time for solving 00:00:00.619218
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:03.549315 Time for solving 00:00:01.994280
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:07.624906 Time for solving 00:00:05.370113
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:14.996559 Time for solving 00:00:12.268905
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:26.960360 Time for solving 00:00:25.649507
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:00:47.112816 Time for solving 00:00:50.689021
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:01:21.001119 Time for solving 00:01:28.554897
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:02:06.799302 Time for solving 00:02:30.355490
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:03:15.008663 Time for solving 00:04:06.016281
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:04:55.586688 Time for solving 00:06:26.601016
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 ****
 Time for assembly 00:07:23.095688 Time for solving 00:10:04.850374
 TPZBoostGraph::Resequence started 
 NNodes 289
 NElements 96
 Original Bandwidth: 173
 Original profile: 11067
 Reverse Expensive Cuthill-McKee ordering:
 TPZBoostGraph::Resequence finished 
 Warning: the current language does not match this frame.
 ****
 Time for assembly 00:10:25.432045 Time for solving 00:14:32.188744
*/

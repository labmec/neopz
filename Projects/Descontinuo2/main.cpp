
#include "oneElementMesh.cpp"
#include "reflectedShock.cpp"
#include "reflectedShock_nonalignedmesh.cpp"
#include "SimpleShock.cpp"
#include "ShockTube2d.cpp"
#include "SubsonicRadialShock.cpp"
#include "NACA4digit.cpp"
#include "pzeuleranalysis.h"
#include "pzconslaw.h"
#include "pzmaterial.h"
#include "pzeulerconslaw.h"
#include "pzartdiff.h"
#include "pzreal.h"
#include "pzvec.h"
#include "pzflowcmesh.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include <iostream>
#include <fstream>
#include "TPZGeoElement.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzrefquad.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzbstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzblock.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontNonSym.h"
#include "TPBSpStructMatrix.h"
#include "pzstring.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "pzrenumbering.h"
#include "pzvisualmatrix.h"


int gDebug;

void error(char *) {}

void InitialGuess(TPZVec<REAL> &x,TPZVec<REAL> &result){
// makes nothing
   result.Resize(4 /*(2 dimensões + rho + energy)*/);
   result.Fill(0.);
}

int main()
{
   gDebug = 0;

   //Creating the computational and geometric meshes.

   TPZString filename, file;
   int ProblemType;
   int temp;
   char number[32];
   TPZArtDiffType DiffType;
   TPZTimeDiscr Diff_TD,
                ConvVol_TD,
		ConvFace_TD;
   REAL delta = 0.;
   REAL CFL = 0.;
   int EvolCFL = 0;
   int MaxIter = 100;
   int p = 0;
   int nSubdiv = 0;
   TPZFlowCompMesh * cmesh;

   cout << "\nProblem type:\n\t0: OneElement\n\t1: SimpleShock\n\t2: ReflectedShock\n\t3: ReflectedShock - NonAlignedMesh\n\t4: ShockTube\n\t5: RadialShock\n\t6: NACA\n";

   cin >> ProblemType;

   cout << "\nDiffusion type:\n\t0: None\n\t1: LS\n\t2: SUPG\n\t3: Bornhaus\n\t4: Transposed LS\n";

   cin >> temp;

   if(temp == 0)
   {
      DiffType = None_AD;
      Diff_TD  = None_TD;
      filename += "NoDiff_";
   }else{
      switch(temp)
      {
         case 1:
	 DiffType = LeastSquares_AD;
	 filename += "LeastSqr";
	 break;
	 case 2:
	 DiffType = SUPG_AD;
	 filename += "SUPG";
	 break;
	 case 3:
	 DiffType = Bornhaus_AD;
	 filename += "Bornhaus";
	 break;
	 case 4:
	 DiffType = TrnLeastSquares_AD;
	 filename += "TrnLeastSqr";
	 break;
      }

      cout << "\nDiffusion Time Discr:\n\t0: ApproxImplicit\n\t1: Implicit\n\t2: Explicit\n";


      cin >> temp;

      switch(temp)
      {
         case 0:
         Diff_TD = ApproxImplicit_TD;
         filename += "ApproxImpl=";
         break;
         case 1:
         Diff_TD = Implicit_TD;
         filename += "Impl=";
         break;
         case 2:
         Diff_TD = Explicit_TD;
         filename += "Expl=";
         break;
      }

      cout << "\nDelta\n";
      cin >> delta;
      sprintf(number, "%lf_", delta);
      filename += number;
   }

   cout << "\nVolume convective Time Discr:\n\t0: None\n\t1: Implicit\n\t2: Explicit\n";

   cin >> temp;

   switch(temp)
   {
      case 0:
      ConvVol_TD = None_TD;
      filename += "NoConvVol_";
      break;
      case 1:
      ConvVol_TD = Implicit_TD;
      filename += "ImplConvVol_";
      break;
      case 2:
      ConvVol_TD = Explicit_TD;
      filename += "ExplConvVol_";
      break;
   }


   cout << "\nFace convective Time Discr:\n\t0: None\n\t1: Implicit\n\t2: Explicit\n";

   cin >> temp;

   switch(temp)
   {
      case 0:
      ConvFace_TD = None_TD;
      filename += "NoConvFace_";
      break;
      case 1:
      ConvFace_TD = Implicit_TD;
      filename += "ImplConvFace_";
      break;
      case 2:
      ConvFace_TD = Explicit_TD;
      filename += "ExplConvFace_";
      break;
   }

   cout << "\nCFL\n";
   cin >> CFL;
   sprintf(number, "CFL%lf_", CFL);
   filename += number;

   cout << "\nEvolute CFL? 0:[no] 1:[yes] 2:[super]\n";
   cin >> EvolCFL;
   if(EvolCFL == 1)
   {
      filename += "EvolCFL_";
   }


   cout << "\nInterpolation degree\n";
   cin >> p;
   sprintf(number, "P%d_", p);
   filename += number;

   cout << "\nNumber of Subdivisions\n";
   cin >> nSubdiv;
   sprintf(number, "N%d", nSubdiv);
   filename += number;

   switch(ProblemType)
   {
      case 0:
      file = "OE_";
      cmesh =
        OneElCompMesh(CFL, delta, p, nSubdiv, DiffType,
	              Diff_TD, ConvVol_TD, ConvFace_TD);
      break;
      case 1:
      file = "SS_";
      cmesh =
         SSCompMesh(CFL, delta, p, nSubdiv, DiffType,
	            Diff_TD, ConvVol_TD, ConvFace_TD);
      break;
      case 2:
      file = "RS_";
      cmesh =
         RSCompMesh(CFL, delta, p, nSubdiv, DiffType,
	            Diff_TD, ConvVol_TD, ConvFace_TD);
      break;
      case 3:
      file = "RSNA_";
      cmesh =
         RSNACompMesh(CFL, delta, p, nSubdiv, DiffType,
	              Diff_TD, ConvVol_TD, ConvFace_TD);
      break;
      case 4:
      file = "ST_";
      cmesh =
         STCompMesh(CFL, delta, p, nSubdiv, DiffType,
	            Diff_TD, ConvVol_TD, ConvFace_TD);
      break;
      case 5:
      file = "SRS_";
      cmesh =
         SRSCompMesh(CFL, delta, p, nSubdiv, DiffType,
	            Diff_TD, ConvVol_TD, ConvFace_TD);
      case 6:
      file = "NACA_";
      cmesh =
         NACACompMesh(CFL, delta, p, nSubdiv, DiffType,
	            Diff_TD, ConvVol_TD, ConvFace_TD);
      break;
   }
   file += filename;

   cout << "\nMaxIter\n";
   cin >> MaxIter;


// Creating the analysis object

   ofstream anFile("analysis.out");
   TPZManVector<REAL,3> normal(3,0.);
   normal[0] = 1.;
   normal[1] = 1.e-5;
   ResequenceByGeometry(cmesh,normal);
   TPZFMatrix fillin(100,100);
   cmesh->ComputeFillIn(100,fillin);
   VisualMatrix(fillin,"matrix.dx");
   TPZEulerAnalysis An(cmesh, anFile);
   An.SetEvolCFL(EvolCFL);

//Solver attributes

/*
   { // LU
   TPZFStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   TPZMatrix * mat = StrMatrix.Create();

   An.SetLinSysCriteria(1e-10, 10);
   An.SetNewtonCriteria(1e-10, 10);
   An.SetTimeIntCriteria(1e-10, MaxIter);

   TPZStepSolver Solver;
   Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
   Solver.SetMatrix(mat);
   An.SetSolver(Solver);
   }
//

   { // GMRES
   TPZSpStructMatrix StrMatrix(cmesh);
      //TPZFStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   TPZMatrix * mat = StrMatrix.Create();

   An.SetLinSysCriteria(1e-8, 100);
   An.SetNewtonCriteria(1e-8, 8);
   An.SetTimeIntCriteria(1e-8,MaxIter);

   //Preconditioner
   TPZStepSolver Pre;
   //Main Solver
//   Pre.SetSSOR(100, 1.1,
//		1e-10,
//		0);

   Pre.SetJacobi(1,//numiterations,
		 1e-8,//tol
		  0);//From Current
   Pre.SetMatrix(mat);

   //Main Solver
   TPZStepSolver Solver;
   Solver.SetGMRES(10000,
		1000,
		Pre,
		1e-9,
		0);
   Solver.SetMatrix(mat);
   An.SetSolver(Solver);
   }
//



*/
/*
   { // GMRES with block preconditioning
   TPZSpStructMatrix StrMatrix(cmesh);
      //TPZFStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   TPZMatrix * mat = StrMatrix.Create();

   An.SetLinSysCriteria(1e-8, 100);
   An.SetNewtonCriteria(1e-8, 4);
   An.SetTimeIntCriteria(1e-8,MaxIter);

   //Preconditioner
   TPZStepSolver Pre;
   //Main Solver
//   Pre.SetSSOR(100, 1.1,
//		1e-10,
//		0);

   TPZBlockDiagonalStructMatrix strBlockDiag(cmesh);
   //just to retrieve blocksizes
//   TPZVec<int> blocksizes;
//   blockDiag.BlockSizes(blocksizes);
//   TPZBlockDiagonal * block = new TPZBlockDiagonal(blocksizes);
   TPZBlockDiagonal * block = new TPZBlockDiagonal();//blockDiag.Create();
   strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure

   Pre.SetMatrix(block);
   Pre.SetDirect(ELU);

   //Main Solver
   TPZStepSolver Solver;
   Solver.SetGMRES(1000,
		300,
		Pre,
		1e-9,
		0);
   Solver.SetMatrix(mat);
   An.SetSolver(Solver);
   An.SetBlockDiagonalPrecond(block);
   }
//
*/

/*
   { // SSOR
   TPZFStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   TPZMatrix * mat = StrMatrix.Create();

   An.SetLinSysCriteria(1e-10, 1000);
   An.SetNewtonCriteria(1e-9, 10);
   An.SetTimeIntCriteria(1e-8,100);

   //Main Solver
   TPZStepSolver Solver;
   Solver.SetSSOR(1000, 1.1,
		1e-10,
		0);
   Solver.SetMatrix(mat);
   An.SetSolver(Solver);
   }
//*/


  TPZParFrontStructMatrix <TPZFrontNonSym> StrMatrix(cmesh);
  An.SetStructuralMatrix(StrMatrix);

  TPZMatrix * mat = NULL;//StrMatrix.CreateAssemble(An.Rhs());

   An.SetLinSysCriteria(1e-8, 100);
   An.SetNewtonCriteria(1e-8, 4);
   An.SetTimeIntCriteria(1e-8,MaxIter);

  TPZStepSolver Solver;
  Solver.SetDirect(ELU);
  Solver.SetMatrix(mat);
  An.SetSolver(Solver);

   cout << "Generating File:" << file.Str() << endl;

   ofstream * dxout = new ofstream((file+".dx" ).Str());
   ofstream *   out = new ofstream((file+".csv").Str());

   An.Run(* out, * dxout, max(0, p-1));

   return 0;
}

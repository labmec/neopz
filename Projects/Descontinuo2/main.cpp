
#include "oneElementMesh.cpp"
#include "reflectedShock.cpp"
#include "reflectedShock_nonalignedmesh.cpp"
#include "SimpleShock.cpp"
#include "ShockTube2d.cpp"
#include "SubsonicRadialShock.cpp"
#include "NACA4digit.cpp"
#include "sphere3D.cpp"
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
#include "pzsave.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "tpzoutofrange.h"
#include "pzlog.h"

 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;

int gDebug;

//void error(char *) {}

void InitialGuess(TPZVec<REAL> &x,TPZVec<REAL> &result){
// makes nothing
   result.Resize(4 /*(2 dimensões + rho + energy)*/);
   result.Fill(0.);
}


void GetSolutionGraph (int bc_id, ostream &arq, TPZFlowCompMesh *cmesh){
int nelem = cmesh->NElements();
int i,j;
TPZVec<int> mark (nelem,0);
TPZStack<TPZGeoEl *> bcStack;
//int nmark = 0;

//Get the selected bcs
 for (i=0;i<nelem;i++){
  TPZCompEl *cel = cmesh->ElementVec()[i];
  if (!cel) continue;
  TPZGeoEl * pGEl;
  pGEl = cel->Reference();
  int matid = cel->Material()->Id();//pGEl->MaterialId();
  if (matid== bc_id){
    bcStack.Push(cel->Reference());
  }
 }

 int nBcs = bcStack.NElements();

 //Get the neighbors of the select bcs
 for (i=0;i<nBcs;i++){
   TPZGeoEl *gel = bcStack[i];
   int nsides = gel->NSides();
   for (j=0;j<nsides;j++){
     TPZGeoElSide thisside (gel,j);
     TPZGeoElSide neighbour = thisside.Neighbour();
     while (neighbour.Exists() && neighbour != thisside){
       TPZCompEl * celneig = neighbour.Element()->Reference();
       if (celneig) {
         int index = celneig->Index();
         mark[index] = 1;
       }
       neighbour = neighbour.Neighbour();
     }
   }
 }

 int nstate = cmesh->MaterialVec()[0]->NStateVariables();

  arq << "State Variables\n";

 for (i=0;i<nelem;i++){
   if (!mark[i]) continue;
   TPZCompEl *cel = cmesh->ElementVec()[i];
   if (!cel || cel->Type() != EDiscontinuous) continue;
   TPZGeoEl *gel = cel->Reference();
   TPZVec<REAL> coord (3,0.), coordX(3, 0.);
   gel->CenterPoint(gel->NSides()-1,coord);
   gel->X(coord, coordX);

   int nconnects = cel->NConnects();
   int connectindex = cel->ConnectIndex(nconnects-1);
   if(connectindex < 0 )continue;

   for(j = 0; j < coordX.NElements(); j++)
   {
      arq << "\t" << coordX[j];
   }

   int seqnum = cmesh->ConnectVec()[connectindex].SequenceNumber();

   int pos = cmesh->Block().Position(seqnum);
   int size = cmesh->Block().Size(seqnum);

   //for (int sz = 0; sz < size; sz++)
   for(j = 0; j < nstate; j++)
   {
      //arq << "\t" << cmesh->Solution()(pos + size - nstate + j);
      arq << "\t" << cmesh->Solution()(pos + size - nstate + j , 0);
   }
   arq << "\n";
 }
}

// saveable test
int main1()
{
   RegisterMeshClasses();
   RegisterMatrixClasses();
   RegisterMaterialClasses();

   TPZEulerConsLaw2 euler(3, 10, 1.5, 3, SUPG_AD), * peuler2;
   euler.SetTimeDiscr(Implicit_TD, ApproxImplicit_TD, None_TD);

   {
      TPZFileStream fstr;
      fstr.OpenWrite("dump.dat");
      euler.Write(fstr,1);
   }


   {
      TPZFileStream fstr;
      fstr.OpenRead("dump.dat");
      TPZSaveable *sv = TPZSaveable::Restore(fstr,NULL);
      peuler2 = dynamic_cast<TPZEulerConsLaw2*>(sv);
   }

   return 0;
}

int run(istream & input, ostream & output)
{

   RegisterMeshClasses();
   RegisterMatrixClasses();
   RegisterMaterialClasses();

   gDebug = 0;

   //Creating the computational and geometric meshes.

   TPZString filename, file, startFileName;
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
   TPZGeoMesh * gmesh;

   output << "\nProblem type:\n\t0: OneElement\n\t1: SimpleShock\n\t2: ReflectedShock\n\t3: ReflectedShock - NonAlignedMesh\n\t4: ShockTube\n\t5: RadialShock\n\t6: NACA\n\t7: GenerateNACAProfile\n\t8: From File\n\t9: Sphere3D\n";

   input >> ProblemType;

   if(ProblemType != 7)
   {

   output << "\nDiffusion type:\n\t0: None\n\t1: LS\n\t2: SUPG\n\t3: Bornhaus\n\t4: Transposed LS\n";


      input >> temp;

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

         output << "\nDiffusion Time Discr:\n\t0: ApproxImplicit\n\t1: Implicit\n\t2: Explicit\n";


         input >> temp;

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

         output << "\nDelta\n";
         input >> delta;
         sprintf(number, "%lf_", delta);
         filename += number;
      }

      output << "\nVolume convective Time Discr:\n\t0: None\n\t1: Implicit\n\t2: Explicit\n";

      input >> temp;

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


      output << "\nFace convective Time Discr:\n\t0: None\n\t1: Implicit\n\t2: Explicit\n";

      input >> temp;

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

      output << "\nCFL\n";
      input >> CFL;
      sprintf(number, "CFL%lf_", CFL);
      filename += number;

      output << "\nEvolute CFL? 0:[no] 1:[yes] 2:[super]\n";
      input >> EvolCFL;
      if(EvolCFL == 1)
      {
         filename += "EvolCFL_";
      }

   }

   if(ProblemType<7 || ProblemType > 8)
   {
      output << "\nInterpolation degree\n";
      input >> p;
      sprintf(number, "P%d_", p);
      filename += number;

      output << "\nNumber of Subdivisions\n";
      input >> nSubdiv;
      sprintf(number, "N%d", nSubdiv);
      filename += number;
   }


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

      case 7:
      case 8:
      {
         file = "FromFile_";

         output << "\nEnter filename to restart from [without extension]:\n";
	 char inputChar[1024];
         input >> inputChar;
	 startFileName = inputChar;

	 TPZMaterial * pmat;
	 TPZEulerConsLaw2 * pEuler;

	 startFileName += ".pzf";
	 TPZFileStream fstr;
         fstr.OpenRead(startFileName.Str());
         TPZSaveable *sv = TPZSaveable::Restore(fstr,NULL);
	 gmesh = dynamic_cast<TPZGeoMesh *>(sv);
	 sv = TPZSaveable::Restore(fstr, gmesh);
         cmesh = dynamic_cast<TPZFlowCompMesh *>(sv);
         cmesh->SetCFL(CFL);
	 pmat = cmesh->GetFlowMaterial(0);
	 pEuler = dynamic_cast<TPZEulerConsLaw2 *>(pmat);
	 pEuler->SetTimeDiscr
	            (Diff_TD,
                     ConvVol_TD,
		     ConvFace_TD);
	 pEuler->ArtDiff().SetArtDiffType(DiffType);
	 pEuler->ArtDiff().SetDelta(delta);


	 if(ProblemType == 8)
	 {
	   output << "Interpolation Degree:\n";
	   input >> p;
	   sprintf(number, "P%d_", p);
	   filename += number;
	   filename += startFileName;
	 }

	 int i_el, nEl = cmesh->NElements();
	 TPZCompElDisc * pCEl;
	 TPZCompEl * pEl;
	 TPZGeoEl * pGEl;
	 TPZInterfaceElement * pIEl;


         for(i_el = 0; i_el < nEl; i_el++)
	 {
	    pEl = cmesh->ElementVec()[i_el];
	    pCEl = dynamic_cast<TPZCompElDisc *>(pEl);
	    pIEl = dynamic_cast<TPZInterfaceElement *>(pEl);
	    if(!pIEl && pCEl && ProblemType == 8)pCEl -> SetDegree(p);
	    pGEl = pEl->Reference();
	    pGEl->SetReference(pEl); // building cross references
	 }

         if(ProblemType == 8)
	 {
	    cmesh->ExpandSolution2();
	 }else
	 {
           ofstream arq("profile.csv");
           GetSolutionGraph (-1, arq, cmesh);
	   return 0;
	 }

      }
      break;
      case 9:
      file = "Sph3D_";
      cmesh =
         SphereCompMesh(CFL, delta, p, nSubdiv, DiffType,
	            Diff_TD, ConvVol_TD, ConvFace_TD);
      break;
   }
   file += filename;

   // computing number of generated elements and faces
   int nEl = cmesh->NElements();
   int nfaces= 0, nelements = 0;
   for(int i = 0; i < nEl; i++)
   {
      TPZCompElDisc * pCEl;
      TPZCompEl * pEl;
      TPZInterfaceElement * pIEl;
      pEl = cmesh->ElementVec()[i];
      pCEl = dynamic_cast<TPZCompElDisc *>(pEl);
      pIEl = dynamic_cast<TPZInterfaceElement *>(pEl);
      if(pCEl)nelements++;
      if(pIEl)nfaces++;
   }
   output << "Number of elements:" << nelements << " faces:" << nfaces << "\n";


   output << "\nMaxIter\n";
   input >> MaxIter;


// Creating the analysis object

   ofstream anFile("analysis.out");

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
   Solver.SetGMRES(3000,
		3000,
		Pre,
		1e-9,
		0);
   Solver.SetMatrix(mat);
   An.SetSolver(Solver);
   An.SetBlockDiagonalPrecond(block);
   }
//


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
/*
   TPZManVector<REAL,3> normal(3,0.);
   normal[0] = 1.;
   normal[1] = 1.e-5;
   ResequenceByGeometry(cmesh,normal);
   TPZFMatrix fillin(100,100);
   cmesh->ComputeFillIn(100,fillin);
   VisualMatrix(fillin,"matrix.dx");

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
*/
   output << "Generating File:" << file.Str() << endl;

   ofstream * dxout = new ofstream((file+".dx" ).Str());
   ofstream *   out = new ofstream((file+".csv").Str());


   An.Run(* out, * dxout, max(0, p-1));

   An.WriteCMesh((file+".pzf").Str());

   return 0;
}

int main()
{
  InitializePZLOG();
  //TPZOutofRange obj;
   try
   {
      run(cin, cout);
   } catch(TPZOutofRange obj)
   {
      cout << "main programa nao terminou normalmente\n";
      return -1;
   }
   return 0;
}

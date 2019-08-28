
#include "oneElementMesh.h"
#include "reflectedShock.h"
#include "reflectedShock_nonalignedmesh.h"
#include "SimpleShock.h"
#include "ShockTube2d.h"
#include "SubsonicRadialShock.h"
#include "NACA4digit.h"
#include "sphere3D.h"
#include "pzeuleranalysis.h"
#include "pzconslaw.h"
#include "TPZMaterial.h"
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
#include "TPZRenumbering.h"
#include "pzvisualmatrix.h"
#include "pzsave.h"
#include "TPZCompElDisc.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "tpzoutofrange.h"
#include "pzlog.h"

#include "pztransfer.h"

using namespace std;

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.converge"));
static LoggerPtr logtr(Logger::getLogger("pz.mesh.transfer"));
#endif

int gDebug;

void InitialGuess(TPZVec<REAL> &x,TPZVec<REAL> &result){
	// makes nothing
	result.Resize(4 /*(2 dimensıes + rho + energy)*/);
	result.Fill(0.);
}


void GetSolutionGraph (int bc_id, std::ostream &arq, TPZFlowCompMesh *cmesh){
	int nelem = cmesh->NElements();
	int i,j;
	TPZVec<int> mark (nelem,0);
	TPZStack<TPZGeoEl *> bcStack;
	
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
	
    TPZMaterial * mat = cmesh->MaterialVec().rbegin()->second;
	int nstate = mat->NStateVariables();
	
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
		
		for(j = 0; j < nstate; j++)
		{
			arq << "\t" << cmesh->Solution()(pos + size - nstate + j , 0);
		}
		arq << "\n";
	}
}

// saveable test
int main1()
{
	
	TPZEulerConsLaw euler(3, 10, 1.5, 3, SUPG_AD), * peuler2;
	euler.SetTimeDiscr(Implicit_TD, ApproxImplicit_TD, None_TD);
	
	{
		TPZFileStream fstr;
		fstr.OpenWrite("dump.dat");
		euler.Write(fstr,1);
	}
	
	
	{
		TPZFileStream fstr;
		fstr.OpenRead("dump.dat");
		TPZSavable *sv = TPZSavable::CreateInstance(fstr,NULL);
		peuler2 = dynamic_cast<TPZEulerConsLaw*>(sv);
	}
	
	return 0;
}

int run(std::istream & input, std::ostream & output)
{
	gDebug = 0;
	
	//Creating the computational and geometric meshes.
	
	TPZString filename, file, startFileName;
	int ProblemType;
	int temp;
	char number[32];
	TPZArtDiffType DiffType;
	TPZTimeDiscr Diff_TD, ConvVol_TD, ConvFace_TD;
	REAL delta = 0.;
	REAL CFL = 0.;
	int EvolCFL = 0;
	int MaxIter = 100;
	int p = 0;
	int nSubdiv = 0;
	TPZFlowCompMesh * cmesh;
	TPZGeoMesh * gmesh = new TPZGeoMesh;
    cmesh = new TPZFlowCompMesh(gmesh);
    
	std::ofstream options("options.txt");
	
	output << "\nProblem type:\n\t0: OneElement\n\t1: SimpleShock\n\t2: ReflectedShock\n\t3: ReflectedShock - NonAlignedMesh\n\t4: ShockTube\n\t5: RadialShock\n\t6: NACA\n\t7: GenerateNACAProfile\n\t8: From File\n\t9: Sphere3D\n";
	
	input >> ProblemType;
	options << ProblemType << std::endl;
	
	if(ProblemType != 7)
	{
		output << "\nDiffusion type:\n\t0: None\n\t1: LS\n\t2: SUPG\n\t3: Bornhaus\n\t4: Transposed LS\n";
		
		
		input >> temp;
		options << temp << std::endl;
		
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
			
			output << "\nInterpolation Space:\n\t0: Discontinuous\n\t1: Continuous\n";
			
			input >> temp;
			options << temp << std::endl;
			
			switch(temp)
			{
				case 0:
					cmesh->SetAllCreateFunctionsDiscontinuous();
                    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
					filename += "_disc";
					break;
				case 1:
					cmesh->SetAllCreateFunctionsContinuous();
					filename += "_cont";
					break;
				default:
					output << "Wrong parameter, setting discontinuous\n";
					cmesh->SetAllCreateFunctionsDiscontinuous();
                    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
					filename += "_disc";
					break;
			}
			
			output << "\nDiffusion Time Discr:\n\t0: ApproxImplicit\n\t1: Implicit\n\t2: Explicit\n";
			
			
			input >> temp;
			options << temp << std::endl;
			
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
			options << delta << std::endl;
			sprintf(number, "%lf_", ((double)delta));
			filename += number;
		}
		
		output << "\nVolume convective Time Discr:\n\t0: None\n\t1: Implicit\n\t2: Explicit\n";
		
		input >> temp;
		options << temp << std::endl;
		
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
		
		
		output << "\nFace convective Time Discr:\n\t0: None\n\t1: Implicit\n\t2: Explicit\n\t3: Approx Implicit\n";
		
		input >> temp;
		options << temp << std::endl;
		
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
			case 3:
				ConvFace_TD = ApproxImplicit_TD;
				filename += "ApproxImplConvFace_";
		}
		
		output << "\nCFL\n";
		input >> CFL;
		options << CFL << std::endl;
		sprintf(number, "CFL%lf_", (double)CFL);
		filename += number;
		
		output << "\nEvolute CFL? 0:[no] 1:[yes] 2:[super]\n";
		input >> EvolCFL;
		options << EvolCFL << std::endl;
		if(EvolCFL == 1)
		{
			filename += "EvolCFL_";
		}
		
	}
	
	if(ProblemType<7 || ProblemType > 8)
	{
		output << "\nInterpolation degree\n";
		input >> p;
		options << p << std::endl;
		sprintf(number, "P%d_", p);
		filename += number;
		
		output << "\nNumber of Subdivisions\n";
		input >> nSubdiv;
		options << nSubdiv << std::endl;
		sprintf(number, "N%d", nSubdiv);
		filename += number;
	}
	
	
	switch(ProblemType)
	{
		case 0:
			file = "OE_";
			cmesh =
			OneElCompMesh(cmesh,CFL, delta, p, nSubdiv, DiffType,
						  Diff_TD, ConvVol_TD, ConvFace_TD);
			break;
		case 1:
			file = "SS_";
			cmesh =
			SSCompMesh(cmesh,CFL, delta, p, nSubdiv, DiffType,
					   Diff_TD, ConvVol_TD, ConvFace_TD);
			break;
		case 2:
			file = "RS_";
			cmesh =
			RSCompMesh(cmesh,CFL, delta, p, nSubdiv, DiffType,
					   Diff_TD, ConvVol_TD, ConvFace_TD);
			break;
		case 3:
			file = "RSNA_";
			cmesh =
			RSNACompMesh(cmesh,CFL, delta, p, nSubdiv, DiffType,
						 Diff_TD, ConvVol_TD, ConvFace_TD);
			break;
		case 4:
			file = "ST_";
			cmesh =
			STCompMesh(cmesh,CFL, delta, p, nSubdiv, DiffType,
					   Diff_TD, ConvVol_TD, ConvFace_TD);
			break;
		case 5:
			file = "SRS_";
			cmesh =
			SRSCompMesh(cmesh, CFL, delta, p, nSubdiv, DiffType,
						Diff_TD, ConvVol_TD, ConvFace_TD);
		case 6:
			file = "NACA_";
			cmesh =
			NACACompMesh(cmesh, CFL, delta, p, nSubdiv, DiffType,
						 Diff_TD, ConvVol_TD, ConvFace_TD,options);
			break;
			
		case 7:
		case 8:
		{
			file = "FromFile_";
			
			output << "\nEnter filename to restart from [without extension]:\n";
			char inputChar[1024];
			input >> inputChar;
			options << inputChar << std::endl;
			startFileName = inputChar;
			
			TPZMaterial * pmat;
			TPZEulerConsLaw * pEuler;
			
			startFileName += ".pzf";
			TPZFileStream fstr;
			fstr.OpenRead(startFileName.Str());
			TPZSavable *sv = TPZSavable::CreateInstance(fstr,NULL);
			gmesh = dynamic_cast<TPZGeoMesh *>(sv);
			std::ofstream gout("geomesh.txt");
			gmesh->Print(gout);
			sv = TPZSavable::CreateInstance(fstr, gmesh);
			cmesh = dynamic_cast<TPZFlowCompMesh *>(sv);
			cmesh->SetCFL(CFL);
			pmat = cmesh->GetFlowMaterial();
			pEuler = dynamic_cast<TPZEulerConsLaw *>(pmat);
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
				options << p << std::endl;
				sprintf(number, "P%d_", p);
				filename += number;
				filename += startFileName;
			}
			
			int i_el, nEl = cmesh->NElements();
			TPZCompElDisc * pCEl;
			TPZCompEl * pEl;
			TPZGeoEl * pGEl;
			TPZInterfaceElement * pIEl;
			TPZCompMesh tempmesh(*cmesh);
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
			TPZTransfer<STATE> transfer;
			cmesh->BuildTransferMatrix(tempmesh,transfer);
			TPZFMatrix<STATE> coarsesol = tempmesh.Solution();
			coarsesol.Remodel(4,coarsesol.Rows()/4);
			coarsesol.Transpose();
			TPZFMatrix<STATE> finesol;
			transfer.Multiply(coarsesol,finesol,0);
			finesol.Transpose();
			finesol.Remodel(finesol.Rows()*finesol.Cols(),1);
			cmesh->LoadSolution(finesol);
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
	int numintel = 0;
	for(int i = 0; i < nEl; i++)
	{
		TPZCompElDisc * pCEl;
		TPZCompEl * pEl;
		TPZInterfaceElement * pIEl;
		TPZInterpolatedElement *pIntel;
		pEl = cmesh->ElementVec()[i];
		pCEl = dynamic_cast<TPZCompElDisc *>(pEl);
		pIEl = dynamic_cast<TPZInterfaceElement *>(pEl);
		pIntel = dynamic_cast<TPZInterpolatedElement *> (pEl);
		if(pCEl)nelements++;
		if(pIEl)nfaces++;
		if(pIntel) numintel++;
	}
	output << "Number of disc elements:" << nelements << " faces:" << nfaces << " number of cont elements " << numintel << "\n";

	output << "\nMaxIter\n";
	input >> MaxIter;
	options << MaxIter << std::endl;
	
	options.close();
	std::string optionname(&file[0]);
	optionname += ".input";
	rename("options.txt",optionname.c_str());
	
	
	// Creating the analysis object
	
	ofstream anFile("analysis.out");
	
	TPZEulerAnalysis An(cmesh, anFile);
	An.SetEvolCFL(EvolCFL);
	
	//Solver attributes
	
	/*
	 { // LU
	 TPZFStructMatrix StrMatrix(cmesh);
	 An.SetStructuralMatrix(StrMatrix);
	 
	 TPZMatrix<REAL> * mat = StrMatrix.Create();
	 
	 An.SetNewtonCriteria(1e-10, 10);
	 An.SetTimeIntCriteria(1e-10, MaxIter);
	 
	 TPZStepSolver Solver;
	 Solver.SetDirect(ELU);// ECholesky -> simÈtrica e positiva definida
	 Solver.SetMatrix(mat);
	 
	 An.SetSolver(Solver);
	 }
	 //
	 
	 { // GMRES
	 TPZSpStructMatrix StrMatrix(cmesh);
	 //TPZFStructMatrix StrMatrix(cmesh);
	 An.SetStructuralMatrix(StrMatrix);
	 
	 TPZMatrix<REAL> * mat = StrMatrix.Create();
	 
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
	 
	 TPZMatrix<REAL> * mat = StrMatrix.Create();
	 
	 int numnewton = 4;
	 REAL NewtonTol = 1.e-8;
	 #ifdef LOG4CXX
	 if(logger->isDebugEnabled())
	 {
	 std::stringstream sout;
	 sout << "Linear system criterium tol " << 1.e-8 << "\nMax number of GMRES iterations " << 100 << endl;
	 sout << "Newton criterium tol " << NewtonTol << "\nMax number of Newton iterations " << numnewton << endl;
	 sout << "Maximum number of steps " << MaxIter;
	 LOGPZ_DEBUG(logger,sout.str().c_str());
	 }
	 
	 #endif   
	 An.SetNewtonCriteria(NewtonTol, numnewton);
	 An.SetTimeIntCriteria(1e-8,MaxIter);
	 An.ComputeTimeStep();
	 
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
	 Solver.SetGMRES(300,
	 20,
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
	 
	 TPZMatrix<REAL> * mat = StrMatrix.Create();
	 
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
	 //  ResequenceByGeometry(cmesh,normal);
	 TPZFMatrix<REAL> fillin(100,100);
	 cmesh->ComputeFillIn(100,fillin);
	 VisualMatrix(fillin,"matrix.dx");
	 
	 TPZFrontStructMatrix <TPZFrontNonSym> StrMatrix(cmesh);
	 StrMatrix.SetQuiet(1);
	 An.SetStructuralMatrix(StrMatrix);
	 
	 TPZMatrix<REAL> * mat = NULL;//StrMatrix.CreateAssemble(An.Rhs());
	 
	 int numnewton = 4;
	 REAL NewtonTol = 1.e-8;
	 #ifdef LOG4CXX
	 if(logger->isDebugEnabled())
	 {
	 std::stringstream sout;
	 sout << "Linear system criterium tol " << 1.e-8 << "\nMax number of GMRES iterations " << 100 << endl;
	 sout << "Newton criterium tol " << NewtonTol << "\nMax number of Newton iterations " << numnewton << endl;
	 sout << "Maximum number of steps " << MaxIter;
	 LOGPZ_DEBUG(logger,sout.str().c_str());
	 }
	 
	 #endif   
	 An.SetNewtonCriteria(NewtonTol, numnewton);
	 An.SetTimeIntCriteria(1e-8,MaxIter);
	 
	 TPZStepSolver Solver;
	 Solver.SetDirect(ELU);
	 Solver.SetMatrix(ma2t);
	 An.SetSolver(Solver);
	 */   
	
	REAL linsystol = 1.e-8;
	int maxiter = 400;
	int numvec = 400;
	An.ComputeTimeStep();
	An.SetGMResBlock(linsystol,maxiter,numvec);
	int numnewton = 10;
	REAL NewtonTol = 1.e-8;
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Linear system criterium tol " << linsystol << "\nMax number of GMRES iterations " << maxiter << endl;
		sout << "Newton criterium tol " << NewtonTol << "\nMax number of Newton iterations " << numnewton << endl;
		sout << "Maximum number of steps " << MaxIter;
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
	
#endif
	An.SetNewtonCriteria(NewtonTol, numnewton);
	An.SetTimeIntCriteria(1e-8,MaxIter);
	
	
#ifdef LOG4CXX
	LOGPZ_DEBUG(logger,file.Str());
#endif
	output << "Generating File:" << file.Str() << endl;
	
	//ofstream * dxout = new ofstream((file+".dx" ).Str());
    std::string outfile(file+".vtk");
	ofstream *   out = new ofstream((file+".csv").Str());
	
	
	An.Run(* out, outfile, max(0, p-1));
	
	An.WriteCMesh((file+".pzf").Str());
	
	return 0;
}

int main()
{
	std::string path;
	std::string configfile;
#ifdef HAVE_CONFIG_H
	path = PZSOURCEDIR;
	path += "/Projects/Descontinuo2/";
#else
	path = "";
#endif
	configfile = path;
	configfile += "log4cxx.cfg";
	//   InitializePZLOG(configfile);
	InitializePZLOG();
	//  TPZInterfaceElement::SetCalcStiffPenalty();
	
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

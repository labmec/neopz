		//#ifdef HAVE_CONFIG_H
//		#include <config.h>
//		#endif
//
//		#include "pzvec.h"
//		#include "pzstack.h"
//		#include "pzfmatrix.h"
//		#include "pzfstrmatrix.h"
//		#include "pzlog.h"
//
//		#include "pzgmesh.h"
//		#include "pzcmesh.h"
//		#include "pzcompel.h"
//		#include "pzgeoelside.h"
//		#include "TPZGeoLinear.h"
//		#include "pzgeopoint.h"
//		#include "tpzgeoblend.h"
//
//		#include <pzgengrid.h>
//		#include "../MeshGeneration.h"
//		#include "../AnalyticalFunctions.h"
//		#include "../ElasticMatInterface2D.h"
//		#include "TPZInterfaceEl.h"
//		#include "pzdiscgal.h"
//
//		#include "TPZRefPattern.h"
//		#include "tpzgeoelrefpattern.h"
//		#include "tpzcompmeshreferred.h"
//		#include "tpzautopointer.h"
//		#include "pzbndcond.h"
//		#include "pzanalysis.h"
//		#include <tpzarc3d.h>
//
//		#include "TPZParSkylineStructMatrix.h"
//		#include "pzstepsolver.h"
//		#include "pzstrmatrix.h"
//		#include "TPZFrontNonSym.h"
//		#include "TPZFrontSym.h"
//		#include "TPBSpStructMatrix.h"
//		#include "TPZSpStructMatrix.h"
//		#include "pzbstrmatrix.h"
//		#include "pzl2projection.h"
//
//		#include "pzpoisson3d.h"
//		#include "pzpoisson3dreferred.h"
//
//		#include "pzelasmat.h"
//		#include "pzmultiphysicselement.h"
//		#include "pzmultiphysicscompel.h"
//		#include "pzbuildmultiphysicsmesh.h"
//		#include "TPZSpStructMatrix.h"
//		#include "../pzporoelastic2d.h"
//		#include "pzlog.h"
//		#include <iostream>
//		#include <string>
//		#include "TPZVTKGeoMesh.h"
//		#include "pzfunction.h"
//		#include "TPZReadGIDGrid.h"
//
//
//		#include <cmath>
//		#include <set>
//
//
//		// Task to do
//		// Document all the code
//		// Include all BC combinations
//		// Include Finite elasticity
//		// Conservation mass mini benchmark (tiago's suggestion)
//		// How I can divide the code in cases for the realistic analysis?
//		// Write Dimensionless Poroelastic Formulation
//
//
//
//		// Using Log4cXX as logging tool
//		//
//		//#ifdef LOG4CXX
//		//static LoggerPtr logger(Logger::getLogger("pz.poroelastic2d"));
//		//#endif
//
//		#ifdef LOG4CXX
//		static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic.data"));
//		#endif
//		//
//		// End Using Log4cXX as logging tool
//
//		// Defintions of names space
//		using namespace std;
//		using namespace pzgeom;	
//
//		// Defitions of Material index
//		// Internal Materials
//		// In this case just one material is used
//		const int TotalMats = 1; 
//		const int matId = 1;	
//
//		// Boundary Conditions Materials
//		// 2D Rectangular Domain
//		const int bcBottom = 2;
//		const int bcRight = 5;
//		const int bcTop = 4;
//		const int bcLeft = 3;
//		const int yfixedPoints = 7;
//		const int xfixedPoints = 8;
//		const int WellLine = 6;
//
//		// Boundary Conditions Definitions
//		// This Definitions are related with the Boudary conditions convention for multphysics simulation (in Poroelastic U/P formulation problems 3 boundary conditions are included) 
//		// not defined
//		const int dirichlet = 0;
//		const int neumann = 1;
//
//
//		// Defintions of Implemented Methods
//
//		//	This Create Computational Meshes
//		TPZCompMesh *ComputationalDiffusionMesh(TPZGeoMesh * gmesh,int pOrder);
//		TPZCompMesh *ComputationalElasticityMesh(TPZGeoMesh * gmesh,int pOrder, bool Initialstress);
//		TPZCompMesh *ComputationalPoroelasticityMesh(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial);
//		TPZCompMesh *SetInitialConditionsbyL2Projection(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &InitialCondition, bool &CustomFunction);
//
//		//	This Set intial Conditions
//		void SetInitialConditions(TPZAnalysis &AnalysisToInitialize, TPZCompMesh * ComputationalMesh, TPZFMatrix<REAL> &ComputationalMeshInitialSolution, void ( *PointerToFunction(REAL &,TPZVec <REAL> &)));
//		void *InitalPressureCalculations(REAL &PressureValue, TPZVec <REAL> &PointCoordinates);
//		void *InitialDisplacementCalculations();
//
//		//	This Solve Different analysis
//		void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);
//		TPZAutoPointer <TPZMatrix<REAL> > MassMatrix(TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics);
//		void SolveSistTransient(TPZVec <REAL> &PrintEachtimeStep, std::string FileName, REAL deltaT,REAL maxTime, TPZFMatrix<REAL> &InitialSolution,TPZVec <TPZPoroElastic2d  * > &mymaterial, 
//								TPZAnalysis &an, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
//		void StiffMatrixLoadVec(TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec);
//
//		//	This are tools for spatial and polinomyal refinement and Postprocess of solutions 
//		void PosProcessElasticity(TPZAnalysis &an, std::string plotfile);
//		void PostProcessDiffusion(TPZAnalysis &an, std::string plotfile);
//		void PostProcessPoroeasticity(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
//		void RefinamentoUniforme(TPZGeoMesh  *gMesh, int nh);
//		void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh, int MatId, int indexEl);
//		void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
//		void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
//		void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);
//		void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file);
//
//		int main(int argc, char *argv[])
//		{
//		#ifdef LOG4CXX
//			InitializePZLOG("../mylog4cxx.cfg");
//		#endif		
//			
//			
//			// GEOMETRICAL MESH CREATION	
//			
//			//	polynomial base degree approach
//			int p=2;
//			
//			//  Used this for Directional Refinement	
//			
//			gRefDBase.InitializeRefPatterns();
//			gRefDBase.InitializeAllUniformRefPatterns();
//			ofstream RefinElPatterns("RefElPatterns.txt");
//			gRefDBase.WriteRefPatternDBase(RefinElPatterns); 	
//			gRefDBase.ReadRefPatternDBase("RefElPatterns.txt");
//			
//			
//			// Begin Circular geometric mesh using tpzarc3d
//			
//			MeshGeneration mygenerator;
//			TPZVec <int> matIdlist(5,0);
//			matIdlist[0]= matId;
//			matIdlist[1]= bcTop;
//			matIdlist[2]= yfixedPoints;
//			matIdlist[3]= xfixedPoints;	
//			matIdlist[4]= WellLine;	
//			REAL ModelRadius = 100.0;
//			REAL Modelthickness = 1.0;
//			mygenerator.Setdimensions(ModelRadius, Modelthickness);
//			TPZGeoMesh * gmesh = mygenerator.GeometricMesh2DValidationQ(matIdlist);	
//			
//			
//			// End Circular mesh	
//			
//			
//			// END GEOMETRICAL MESH CREATION	
//			
//			
//			{
//				//	Print Geometrical Base Mesh
//				ofstream argument("BaseGeoMesh.txt");
//				gmesh->Print(argument);
//				ofstream Dummyfile("BaseGeoMesh.vtk");	
//				TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);	
//			}
//			
//			int Href = 0;
//			RefinamentoUniforme(gmesh, Href);
//			
//			// Directional Point source Refinement 
//			int ndirectdivp = 17;
//			set<int> SETmatPointRefDir2;
//			SETmatPointRefDir2.insert(WellLine);
//			//	SETmatPointRefDir2.insert(bcBottom);
//			//	SETmatPointRefDir2.insert(bcTop);
//			//	SETmatPointRefDir2.insert(bcRight);
//			//	SETmatPointRefDir2.insert(bcLeft);	
//			//	SETmatPointRefDir2.insert(xfixedPoints);
//			//	SETmatPointRefDir2.insert(yfixedPoints);	
//			for(int j = 0; j < ndirectdivp; j++)
//			{
//				int nel = gmesh->NElements();
//				for (int iref = 0; iref < nel; iref++)
//				{
//					TPZVec<TPZGeoEl*> filhos;
//					TPZGeoEl * gelP2 = gmesh->ElementVec()[iref];
//					if(!gelP2 || gelP2->HasSubElement()) continue;
//					TPZRefPatternTools::RefineDirectional(gelP2, SETmatPointRefDir2);
//				}		
//			}	
//			
//			{
//				//	Print Geometrical refined Base Mesh
//				ofstream argument("RefineGeoMesh.txt");
//				gmesh->Print(argument);
//				ofstream Dummyfile("RefineGeoMesh.vtk");
//				TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);	
//			}	
//			
//			
//			//	First computational mesh
//			
//			// initialization problem ...
//			bool Initialstress = false;
//			TPZCompMesh * ComputationalMeshElasticity = ComputationalElasticityMesh(gmesh, p+1, Initialstress);
//			//	Print First computational mesh
//			ofstream ArgumentElasticity("cmeshElasticity.txt");
//			ComputationalMeshElasticity->Print(ArgumentElasticity);
//			
//			//	Second computational mesh
//			TPZCompMesh * ComputationalMeshDiffusion = ComputationalDiffusionMesh(gmesh, p);
//			//	Print Second computational mesh
//			ofstream ArgumentDiffusion("cmeshDiffusion.txt");
//			ComputationalMeshDiffusion->Print(ArgumentDiffusion);	
//			
//			//	Cleaning reference of the geometric mesh for cmesh1
//			gmesh->ResetReference();
//			ComputationalMeshElasticity->LoadReferences();
//			
//			//	Using Uniform Refinement for first computational problem
//			RefinUniformElemComp(ComputationalMeshElasticity,0);
//			ComputationalMeshElasticity->AdjustBoundaryElements();
//			ComputationalMeshElasticity->CleanUpUnconnectedNodes();
//			
//			ofstream ArgumentElasticityRef("cmeshElasticityRef.txt");
//			ComputationalMeshElasticity->Print(ArgumentElasticityRef);
//			ofstream ArgumentElasticityGeoRef("GeoMeshElasticity.txt");
//			gmesh->Print(ArgumentElasticityGeoRef);
//			
//			// Vtk visualization Elasticity Geometric Mesh
//			ofstream DummyfileElasticity("GeoMeshElasticity.vtk");
//			PrintGMeshVTK(gmesh, DummyfileElasticity);
//			
//			// Clear reference of the geometric mesh to cmesh2
//			gmesh->ResetReference();
//			ComputationalMeshDiffusion->LoadReferences();
//			
//			//	Using Uniform Refinement for second computational problem	
//			RefinUniformElemComp(ComputationalMeshDiffusion,0);
//			ComputationalMeshDiffusion->AdjustBoundaryElements();
//			ComputationalMeshDiffusion->CleanUpUnconnectedNodes();
//			//	ComputationalMeshDiffusion->ExpandSolution();
//			
//			ofstream ArgumentDiffusionRef("cmeshDiffusionRef.txt");
//			ComputationalMeshDiffusion->Print(ArgumentDiffusionRef);
//			ofstream ArgumentDiffusionGeoRef("GeoMeshDiffusion.txt");
//			gmesh->Print(ArgumentDiffusionGeoRef);
//			
//			// Vtk visualization Diffusion Geometric Mesh	
//			ofstream DummyfileDiffusion("GeoMeshDiffusion.vtk");
//			PrintGMeshVTK(gmesh, DummyfileDiffusion);
//			
//			
//			//	Set initial conditions for pressure
//			TPZFMatrix<REAL> DiffusionSolution;
//			TPZAnalysis DiffusionAnalysis(ComputationalMeshDiffusion);
//			//	Without L2 projection
//			//	SetInitialConditions(DiffusionAnalysis, ComputationalMeshDiffusion, DiffusionSolution, *InitalPressureCalculations);
//			//	ComputationalMeshDiffusion->Solution() = DiffusionSolution;
//			
//			//	Set initial conditions for pressure by using L2 projection
//			int nrs = DiffusionAnalysis.Solution().Rows();
//			bool CustomFunction = true;
//			//	Using  constant initial pressure
//			TPZVec<REAL> InitialCondition(nrs,0.0);
//			TPZCompMesh  * ComputationalL2 = SetInitialConditionsbyL2Projection(gmesh, p, InitialCondition, CustomFunction);
//			TPZAnalysis L2Analysis(ComputationalL2);
//			SolveSist(L2Analysis, ComputationalL2);
//			DiffusionAnalysis.LoadSolution(L2Analysis.Solution());	
//			
//			
//			//#ifdef LOG4CXX
//			//	if(logdata->isDebugEnabled())
//			//	{		
//			//		std::stringstream sout;
//			//		DiffusionAnalysis.Solution().Print("Intial conditions = ", sout,EMathematicaInput);
//			//		LOGPZ_DEBUG(logdata,sout.str())
//			//	}
//			//#endif	
//			
//			TPZVec<TPZCompMesh *> ComputationalMeshVector(2);
//			ComputationalMeshVector[0] = ComputationalMeshElasticity;
//			ComputationalMeshVector[1] = ComputationalMeshDiffusion;
//			TPZVec <TPZPoroElastic2d *>  materialist(TotalMats,0) ;
//			TPZCompMesh * ComputationalMeshPoroelasticity = ComputationalPoroelasticityMesh(gmesh,ComputationalMeshVector,materialist);
//			
//			
//			// time control
//			REAL hour = 3600;
//			REAL day = 86400;
//			REAL month = 30*day;
//			REAL year = 365*day;
//			
//			REAL InitialTime = 0.0;
//			//	REAL Delta		= 0.0001;
//			//	REAL MaxTime	= 0.1;	
//			//	REAL Delta		= 1000.0;
//			//	REAL MaxTime	= 1.0e7;
//			REAL Delta		= 0.00001;
//			REAL MaxTime	= 1.0;	
//			//	REAL Delta = 0.0164*day;
//			//	REAL MaxTime =  164*day;
//			std::string	OutFileName("PoroelasticLineSourcetheta06Pb2");
//			
//			//	REAL Theta = 0.0;
//			//	REAL Theta = 0.3;
//			
//			//REAL Theta = 0.5;	
//			REAL Theta = 0.6667;
//			//REAL Theta = 1.0;
//			//	REAL Theta = 1+((1/Delta)-(1/(log(1+Delta))));	
//			//	REAL Theta = 1+((1/Delta)-(1/(log(1+Delta))));
//			//	int PrintEachtimeStep = 999;
//			TPZVec <REAL> PrintEachtimeStep(10,0);
//			//	PrintEachtimeStep[0] = 1.0e5;	
//			//	PrintEachtimeStep[1] = 8.0e5;
//			//	PrintEachtimeStep[2] = 2.0e6;	
//			//	PrintEachtimeStep[3] = 3.0e6;
//			PrintEachtimeStep[0] = 1*Delta;	
//			PrintEachtimeStep[1] = 5*Delta;	
//			PrintEachtimeStep[2] = 10*Delta;
//			PrintEachtimeStep[3] = 50*Delta;	
//			PrintEachtimeStep[4] = 100*Delta;	
//			PrintEachtimeStep[5] = 500*Delta;	
//			PrintEachtimeStep[6] = 1000*Delta;	
//			PrintEachtimeStep[7] = 10000*Delta;
//			PrintEachtimeStep[8] = 50000*Delta;	
//			PrintEachtimeStep[9] = 100000*Delta;	
//			//	PrintEachtimeStep[3] = 0.5;
//			//	PrintEachtimeStep[3] = 1.0;	
//			
//			//	// time control
//			//	REAL hour = 3600;
//			//	REAL day = 86400;
//			//	REAL month = 30*day;
//			//	REAL year = 365*day;
//			//
//			//	REAL InitialTime = 0.0;
//			//	//	REAL Delta		= 0.0001;
//			//	//	REAL MaxTime	= 0.1;	
//			//	//	REAL Delta		= 1000.0;
//			//	//	REAL MaxTime	= 1.0e7;
//			//	REAL Delta		= 0.0000000000000001;
//			//	REAL MaxTime	= 0.00000000000001;	
//			//	//	REAL Delta = 0.0164*day;
//			//	//	REAL MaxTime =  164*day;
//			//	std::string	OutFileName("PoroelasticLineSourceAtpb2");
//			//
//			//	//	REAL Theta = 0.0;
//			//	//	REAL Theta = 0.3;
//			//
//			//	//	REAL Theta = 0.5;	
//			//	//	REAL Theta = 0.6667;
//			//	REAL Theta = 1.0;
//			//	//	REAL Theta = 1+((1/Delta)-(1/(log(1+Delta))));	
//			//	//	REAL Theta = 1+((1/Delta)-(1/(log(1+Delta))));
//			//	//	int PrintEachtimeStep = 999;
//			//	TPZVec <REAL> PrintEachtimeStep(10,0);
//			//	//	PrintEachtimeStep[0] = 1.0e5;	
//			//	//	PrintEachtimeStep[1] = 8.0e5;
//			//	//	PrintEachtimeStep[2] = 2.0e6;	
//			//	//	PrintEachtimeStep[3] = 3.0e6;
//			//	PrintEachtimeStep[0] = 100*Delta;	
//			//	PrintEachtimeStep[1] = 50*Delta;	
//			//	PrintEachtimeStep[2] = 100*Delta;
//			//	PrintEachtimeStep[3] = 500*Delta;	
//			//	PrintEachtimeStep[4] = 1000*Delta;	
//			//	PrintEachtimeStep[5] = 5000*Delta;	
//			//	PrintEachtimeStep[6] = 10000*Delta;	
//			//	PrintEachtimeStep[7] = 100000*Delta;
//			//	PrintEachtimeStep[8] = 500000*Delta;	
//			//	PrintEachtimeStep[9] = 1000000*Delta;	
//			//	//	PrintEachtimeStep[3] = 0.5;
//			//	//	PrintEachtimeStep[3] = 1.0;	
//			
//			
//			
//			// Improve this
//			for(int imat = 0; imat < TotalMats; imat++)
//			{	
//				materialist[imat]->SetTimeStep(Delta,Theta);
//				materialist[imat]->SetTimeValue(InitialTime);			
//			}
//			
//			TPZAnalysis PoroelasticAnalysis(ComputationalMeshPoroelasticity);
//			
//			//Setting initial coditions 
//			TPZFMatrix<REAL> InitialSolution = PoroelasticAnalysis.Solution();
//			std::string output(OutFileName);
//			std::stringstream outputfiletemp;
//			outputfiletemp << output << ".vtk";
//			std::string plotfile = outputfiletemp.str();
//			//Print Initial conditions
//			PostProcessPoroeasticity(ComputationalMeshVector,ComputationalMeshPoroelasticity,PoroelasticAnalysis,plotfile);	
//			
//			//	TPZSpStructMatrix matsp(mphysics);
//			//	//TPZSkylineStructMatrix matsp(mphysics);	
//			//	std::set< int > materialid;
//			//	materialid.clear();	
//			//	//	int matid = 1;
//			//	//	materialid.insert(imat);
//			//	//	matsp.SetMaterialIds (materialid);
//			//	// Improve this
//			//	for(int imat = 0; imat < TotalMats; imat++)
//			//	{	
//			//		materialid.insert(imat+1);	
//			//	}
//			//	// Inserting Material ids for matrix K2 creation
//			//	matsp.SetMaterialIds(materialid);		
//			//	
//			//	TPZAutoPointer<TPZGuiInterface> guiInterface;
//			//	TPZFMatrix<REAL> Un;
//			//	TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
//			
//			//#ifdef LOG4CXX
//			//	if(logdata->isDebugEnabled())
//			//	{
//			//		//		std::stringstream sout;
//			//		//		matK2->Print("K2 = ", sout,EMathematicaInput);
//			//		//		LOGPZ_DEBUG(logdata,sout.str())
//			//	}
//			//#endif
//			
//			//	//Criando matriz K1
//			//	// Improve this
//			//	for(int imat = 0; imat < TotalMats; imat++)
//			//	{	
//			//		materialist[imat]->SetCurrentState();
//			//	}	
//			//	TPZSkylineStructMatrix matsk(mphysics);
//			//	TPZFMatrix<REAL> matK1;	
//			//	TPZFMatrix<REAL> fvec; //vetor de carga
//			//	an.SetStructuralMatrix(matsk);//	Set the structtural sky line matrix to our analysis
//			//	
//			//	
//			
//			//	an.Run(); //	Excecute an analysis to obtain the Rhs vector (In this case we start with zero initial values)
//			//	
//			//	matK1 = an.StructMatrix(); //Storage the global matrix and load vector
//			//	fvec = an.Rhs();
//			
//			// This code identify singular blocks
//			TPZStepSolver<REAL> step; //Create Solver object
//			step.SetDirect(ELDLt); //	Symmetric case
//			PoroelasticAnalysis.SetSolver(step); //	Set solver	
//			TPZStepSolver<REAL> & temp = dynamic_cast<TPZStepSolver<REAL> &> (PoroelasticAnalysis.Solver());
//			std::list <int> & zeropivot = temp.Singular(); 
//			if (zeropivot.size()) 
//			{
//				int eq = * zeropivot.begin();
//				PoroelasticAnalysis.Rhs().Zero();
//				PoroelasticAnalysis.Rhs()(eq,0) = -10000.0;
//				PoroelasticAnalysis.Solve();
//				TPZFMatrix<REAL> TempSolution = PoroelasticAnalysis.Solution();
//				
//		#ifdef LOG4CXX
//				// Print the temporal solution
//				if(logdata->isDebugEnabled())
//				{
//					std::stringstream sout;
//					TempSolution.Print("SingularNodes = ", sout,EMathematicaInput);
//					LOGPZ_DEBUG(logdata,sout.str())
//				}
//		#endif	
//				std::string output;
//				output = "SingularNodes";
//				std::stringstream outputfiletemp;
//				outputfiletemp << output << ".vtk";
//				std::string plotfile = outputfiletemp.str();
//				PostProcessPoroeasticity(ComputationalMeshVector,ComputationalMeshPoroelasticity,PoroelasticAnalysis,plotfile);
//				
//				// Probelem bad pose
//				DebugStop();
//			}
//			
//			//#ifdef LOG4CXX
//			//	if(logdata->isDebugEnabled())
//			//	{
//			//		// base system to invert
//			//		// just one for checking purpose
//			//		//		std::stringstream sout;
//			//		//		an.Solver().Matrix()->Print("K1 = ", sout,EMathematicaInput);
//			//		//		fvec.Print("fvec = ", sout,EMathematicaInput);		
//			//		//		//Print the temporal solution
//			//		//		Initialsolution.Print("Intial conditions = ", sout,EMathematicaInput);
//			//		//		TPZFMatrix<REAL> Temp;
//			//		//		TPZFMatrix<REAL> Temp2;
//			//		//		matK2->Multiply(Initialsolution,Temp);
//			//		//		Temp.Print("Temp K2 = ", sout,EMathematicaInput);	
//			//		//		LOGPZ_DEBUG(logdata,sout.str())
//			//	}
//			//#endif
//			
//			///start transient problem
//			SolveSistTransient(PrintEachtimeStep,OutFileName,Delta,MaxTime, InitialSolution, materialist, PoroelasticAnalysis, ComputationalMeshVector,ComputationalMeshPoroelasticity);	
//			
//			return EXIT_SUCCESS;
//		}
//
//
//		TPZCompMesh * ComputationalElasticityMesh(TPZGeoMesh * gmesh,int pOrder, bool Initialstress)
//		{
//			int dim = 2;
//			int planestress = 0;	
//			// All material with equal properties
//			// Needed to include diferents materials and properties via keyword file	
//			REAL gravity = 0.0; // SI system
//			REAL overburdendepth = 2000.0; // SI system
//			REAL layerthickness = 10.0;  // SI system	
//			TPZVec <REAL> E(TotalMats,100.0);
//			TPZVec <REAL> poisson(TotalMats,0.35);
//			TPZVec <REAL> rockrho(TotalMats,2330.0);	
//			TPZVec <REAL> force(2,0.);
//			TPZVec <TPZVec< REAL > > forceVec(TotalMats,force);	
//			TPZVec <TPZElasticityMaterial *> materialist(TotalMats,0);
//			
//			TPZCompEl::SetgOrder(pOrder);
//			TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//			cmesh->SetDimModel(dim);
//			cmesh->SetAllCreateFunctionsContinuous();	
//			
//			
//			
//			for(int imat = 0; imat < TotalMats; imat++)
//			{	
//				
//				// Non Gravity field
//				//		forceVec[imat][1]=-gravity*rockrho[imat];
//				// Creating and inserting material objec in to computational mesh
//				materialist[imat] = new TPZElasticityMaterial(imat+1, E[imat], poisson[imat], forceVec[imat][0], forceVec[imat][1], planestress); 
//				TPZMaterial * mat(materialist[imat]);
//				cmesh->InsertMaterialObject(mat);
//				
//				// This is a dummy boundary conditions
//				// Here is necessary to include the same set of computational information for all computational meshes -> Elastic, Pressure, Biot's Poroelasticity
//				TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
//				
//				// Top Boundary			
//				TPZMaterial * BCTOP = materialist[imat]->CreateBC(mat, bcTop,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCTOP);
//				
//				// Bottom Boundary			
//				TPZMaterial * BCBOTT = materialist[imat]->CreateBC(mat, bcBottom,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCBOTT);
//				
//				// Left Boundary			
//				TPZMaterial * BCLEFT = materialist[imat]->CreateBC(mat, bcLeft,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary			
//				TPZMaterial * BCRIGHT = materialist[imat]->CreateBC(mat, bcRight,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCRIGHT);
//				
//				// XfixedPoints Boundary Restriction			
//				TPZMaterial * BCXfIXED = materialist[imat]->CreateBC(mat, xfixedPoints,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCXfIXED);
//				
//				// YfixedPoints Boundary Restriction
//				TPZMaterial * BCYfIXED = materialist[imat]->CreateBC(mat, yfixedPoints,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCYfIXED);
//				
//				// Line Source
//				TPZMaterial * BCLINE = materialist[imat]->CreateBC(mat, WellLine,neumann, val1, val2);
//				cmesh->InsertMaterialObject(BCLINE);
//				
//			}
//			
//			if (Initialstress) 
//			{
//				
//				//		Inserting interface element
//				ElasticMatInterface2D *interfacemat = new ElasticMatInterface2D(41,100,0.35, 0.0, -gravity*2500.0, planestress); 
//				TPZMaterial * mat(interfacemat);
//				cmesh->InsertMaterialObject(mat);
//				
//				// setting Computational element structure
//				std::set<int> iset;			
//				for(int imat = 0; imat < TotalMats; imat++)
//				{	
//					iset.insert(imat+1);	
//					cmesh->AutoBuild(iset);	
//					cmesh->Reference()->ResetReference();			
//					iset.clear();
//					
//				}
//				cmesh->LoadReferences();
//				
//				cmesh->LoadReferences();
//				iset.insert(-1);	
//				cmesh->AutoBuild(iset);	
//				iset.clear();
//				
//				iset.insert(-2);	
//				cmesh->AutoBuild(iset);	
//				iset.clear();		
//				
//				iset.insert(-3);	
//				cmesh->AutoBuild(iset);	
//				iset.clear();
//				
//				iset.insert(-4);	
//				cmesh->AutoBuild(iset);	
//				iset.clear();
//				
//				iset.insert(-7);	
//				cmesh->AutoBuild(iset);	
//				iset.clear();
//				
//				iset.insert(-8);	
//				cmesh->AutoBuild(iset);	
//				iset.clear();	
//				
//			}
//			else 
//			{
//				cmesh->AutoBuild();
//			}
//			
//			
//			return cmesh;
//			
//		}
//
//		TPZCompMesh * ComputationalDiffusionMesh(TPZGeoMesh * gmesh, int pOrder)
//		{
//			/// criar materiais
//			int dim = 2;
//			// All material with equal properties
//			// Needed include diferents materials and properties via keyword file
//			TPZVec <REAL> convdir(3,0.);
//			TPZVec <REAL> diff(TotalMats,0.1);
//			TPZVec <REAL> conv(TotalMats,0.0);
//			TPZVec <TPZVec< REAL > > convdirVec(TotalMats,convdir);
//			TPZVec <REAL> flux(TotalMats,0.0);
//			
//			TPZVec <TPZMatPoisson3d *> materialist(TotalMats,0);	
//			TPZCompEl::SetgOrder(pOrder);
//			TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//			cmesh->SetDimModel(dim);
//			cmesh->SetAllCreateFunctionsContinuous();	
//			
//			for(int imat = 0; imat < TotalMats; imat++)
//			{
//				// Creating and inserting material objec in to computational mesh		
//				materialist[imat] = new TPZMatPoisson3d(imat+1,dim); 
//				TPZMaterial * mat(materialist[imat]);
//				materialist[imat]->SetParameters(diff[imat], conv[imat], convdirVec[imat]);
//				materialist[imat]->SetInternalFlux(flux[imat]);
//				materialist[imat]->NStateVariables();
//				cmesh->InsertMaterialObject(mat);
//				
//				// This is a dummy boundary conditions
//				// Here is necessary to include the same set of computational information for all computational meshes -> Elastic, Pressure, Biot's Poroelasticity
//				TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
//				
//				// Top Boundary			
//				TPZMaterial * BCTOP = materialist[imat]->CreateBC(mat, bcTop,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCTOP);
//				
//				// Bottom Boundary			
//				TPZMaterial * BCBOTT = materialist[imat]->CreateBC(mat, bcBottom,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCBOTT);
//				
//				// Left Boundary			
//				TPZMaterial * BCLEFT = materialist[imat]->CreateBC(mat, bcLeft,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary			
//				TPZMaterial * BCRIGHT = materialist[imat]->CreateBC(mat, bcRight,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCRIGHT);
//				
//				// XfixedPoints Boundary Restriction			
//				TPZMaterial * BCXfIXED = materialist[imat]->CreateBC(mat, xfixedPoints,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCXfIXED);
//				
//				// YfixedPoints Boundary Restriction
//				TPZMaterial * BCYfIXED = materialist[imat]->CreateBC(mat, yfixedPoints,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCYfIXED);
//				
//				// Line Source
//				TPZMaterial * BCLINE = materialist[imat]->CreateBC(mat, WellLine,neumann, val1, val2);
//				cmesh->InsertMaterialObject(BCLINE);		
//			}
//			
//			
//			//	// setting Computational element structure
//			//	std::set<int> iset;	
//			//	for(int imat = 0; imat < TotalMats; imat++)
//			//	{	
//			//		iset.insert(imat+1);	
//			//		cmesh->AutoBuild(iset);
//			//		cmesh->LoadReferences(); 		
//			//		cmesh->Reference()->ResetReference();			
//			//		iset.clear();
//			//		
//			//	}
//			//	cmesh->Reference()->ResetReference();		
//			cmesh->AutoBuild();
//			
//			return cmesh;
//		}
//
//		TPZCompMesh * ComputationalPoroelasticityMesh(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial)
//		{
//			//Creating computational mesh for multiphysic elements
//			int dim = 2;
//			int planestress = 0; // This is a Plain strain problem
//			
//			gmesh->ResetReference();
//			TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
//			mphysics->SetAllCreateFunctionsMultiphysicElem();
//			TPZAutoPointer<TPZFunction<STATE> > TimeDepForcingF;
//			TPZAutoPointer<TPZFunction<STATE> > TimeDepFExact;
//			TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolution2DLineSource);
//			
//			// Definitions
//			REAL lamb = 8.3e9;					//	[Pa]
//			REAL lambu = 13.4e9;						//	[Pa]
//			REAL alpha = 0.7;						//	[-]
//			REAL G= 5.5e9;							//	[Pa]
//			REAL rhof = 1000.0;						//	[kg/m3]
//			REAL poisson = (lamb)/(2*(lamb+G));						//	[-]
//			REAL Eyoung = (G*(3.0*lamb+2.0*G))/(lamb+G);					//	[Pa]			
//			REAL Ssig = ((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
//			REAL Se = ((pow(alpha,2))/((lambu-lamb)));		
//			REAL Ss = Se;//((pow(alpha,2))/((lambu-lamb)));
//			REAL Coeff = (alpha)/((lamb+2.0*G));
//			REAL rockrho = 2500.0;					//	[kg/m3]
//			REAL gravity = 0.0;					//	[m/s2]
//			REAL c = 0.083;							//	[m2/s]	
//			REAL viscd = 0.001;						//	[Pa.s]
//			REAL permd = c*Ssig*viscd;//7.8784288e-15;//1.109542e-14						//	[m2]
//			REAL visc = 1.0;						//	[Pa.s]
//			REAL perm = 1.0;//7.8784288e-15;//1.109542e-14						//	[m2]	
//			REAL qo = -20.0;	
//			REAL Bo = 1.0;
//			REAL PI = atan(1.)*4.;
//			TPZVec <REAL> force(2,0.);	
//			
//			
//			// Definitions for paramter vectors Al materials have the same porperties
//			TPZVec <REAL> lambVec(TotalMats,lamb*Ssig);
//			TPZVec <REAL> lambuVec(TotalMats,lambu);
//			TPZVec <REAL> alphaVec(TotalMats,alpha);
//			TPZVec <REAL> GVec(TotalMats,G*Ssig);
//			TPZVec <REAL> rhofVec(TotalMats,rhof);
//			TPZVec <REAL> poissonVec(TotalMats,poisson);
//			TPZVec <REAL> EyoungVec(TotalMats,Eyoung);
//			TPZVec <REAL> SsigVec(TotalMats,Ssig);
//			TPZVec <REAL> SeVec(TotalMats,Se);
//			TPZVec <REAL> rockrhoVec(TotalMats,rockrho);
//			TPZVec <REAL> cVec(TotalMats,c);
//			TPZVec <REAL> viscVec(TotalMats,visc);
//			TPZVec <REAL> permVec(TotalMats,perm);
//			TPZVec <REAL> qoVec(TotalMats,qo);
//			TPZVec <REAL> BoVec(TotalMats,Bo);
//			TPZVec <TPZVec< REAL > > forceVec(TotalMats,force);
//			for(int imat = 0; imat < TotalMats; imat++)
//			{	
//				
//				// No Gravity field
//				//		forceVec[imat][1]=gravity*rockrhoVec[imat];		
//				
//				mymaterial[imat] = new TPZPoroElastic2d (imat+1, dim);
//				mymaterial[imat]->SetParameters(lambVec[imat], GVec[imat],forceVec[imat][0], forceVec[imat][1]);
//				mymaterial[imat]->SetParameters(permVec[imat],viscVec[imat]);
//				mymaterial[imat]->SetfPlaneProblem(planestress);
//				mymaterial[imat]->SetBiotParameters(alpha,Se/Ssig);
//				mymaterial[imat]->SetTimeDependentFunctionExact(TimeDepFExact);		
//				
//				ofstream argm("mymaterial1.txt");
//				mymaterial[imat]->Print(argm);
//				TPZMaterial * mat(mymaterial[imat]);
//				mphysics->InsertMaterialObject(mat);
//				
//				///--- --- Inserting Boundary Conditions
//				// Squared Model for validation line source Validation give in ...........					
//				
//				// elastic problem -> 1 ,  Poisson problem -> 2
//				// Boundary condition N1D2 means: Elastic -> N, Pressure -> D Unbounded body with free flux contour
//				// Boundary condition N1N2 means: Elastic -> N, Pressure -> N Unbounded body with no flux contour		
//				TPZFMatrix<REAL> val1AllContour(3,2,0.), val2AllContour(3,1,0.);
//				int D1D2			= 000;		
//				int D1N2			= 001;
//				int N1N2			= 111;
//				int N1D2			= 110;
//				int PM1N2			= 3221;		
//				int WichCase		= N1N2;
//				REAL Nx				= 0.0;//7.48846e-5;
//				REAL Ny				= 0.0;//7.48846e-5;
//				REAL DiffusionCase	= 0.0;
//				//		val1AllContour(0,0)= 1.0;
//				//		val1AllContour(0,1)= 0.0;
//				//		val1AllContour(1,0)= 1.0;
//				//		val1AllContour(1,1)= 0.0;
//				//		val1AllContour(2,0)= 0.0;
//				//		val1AllContour(2,0)= 0.0;	
//				val2AllContour(0,0)= Nx;
//				val2AllContour(1,0)= Ny;
//				val2AllContour(2,0)= DiffusionCase;
//				
//				
//				//		// Top Boundary			
//				//		TPZMaterial * BCTOP = mymaterial[imat]->CreateBC(mat, bcTop, WichCase, val1AllContour, val2AllContour);
//				//		mphysics->InsertMaterialObject(BCTOP);
//				// Bottom Boundary		
//				TPZMaterial * BCBOTT = mymaterial[imat]->CreateBC(mat, bcTop, WichCase, val1AllContour, val2AllContour);
//				mphysics->InsertMaterialObject(BCBOTT);
//				//		// Left Boundary		
//				//		TPZMaterial * BCLEFT = mymaterial[imat]->CreateBC(mat, bcLeft, WichCase, val1AllContour, val2AllContour);
//				//		mphysics->InsertMaterialObject(BCLEFT);
//				//		// Right Boundary		
//				//		TPZMaterial * BCRIGHT = mymaterial[imat]->CreateBC(mat, bcRight, WichCase, val1AllContour, val2AllContour);
//				//		mphysics->InsertMaterialObject(BCRIGHT);
//				
//				/// Setting Constriction to avoid rotatinal and traslational displacement DYFY1N2 means: Elastic -> Dirichlet Y component, Pressure -> Neumman
//				TPZFMatrix<REAL> val1YFix(3,2,0.), val2YFix(3,1,0.0);
//				int DYFX1N2		= 101;
//				REAL uDFixedy	= 0.0;
//				REAL pNFixedy	= 0.0;
//				val2YFix(1,0)=uDFixedy;
//				val2YFix(2,0)=pNFixedy;
//				TPZMaterial * BCyfixedPoints = mymaterial[imat]->CreateBC(mat, yfixedPoints, DYFX1N2, val1YFix, val2YFix);
//				mphysics->InsertMaterialObject(BCyfixedPoints);
//				
//				/// Setting Constriction to avoid rotatinal and traslational displacement DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman
//				TPZFMatrix<REAL> val1XFix(3,2,0.), val2XFix(3,1,0.0);
//				int DXFY1N2		= 011;
//				REAL uDFixedx	= 0.0;
//				REAL pNFixedx	= 0.0;
//				val2XFix(0,0)=uDFixedx;			
//				val2XFix(2,0)=pNFixedx;
//				TPZMaterial * BCxfixedPoints = mymaterial[imat]->CreateBC(mat, xfixedPoints, DXFY1N2, val1XFix, val2XFix);
//				mphysics->InsertMaterialObject(BCxfixedPoints);			
//				
//				//  Source Production/injection Well line (Q is the injected/produced volume per time [kg/s])
//				int Well = 111;
//				TPZFMatrix<REAL> valP1(3,2,0.), valP2(3,1,0.);
//				//		valP2(2,0) = qo/(rhof);
//				//		valP2(2,0) = -4*PI;		
//				valP2(2,0) = -1.0/4;		
//				TPZMaterial * BCLINE = mymaterial[imat]->CreateBC(mat, WellLine, Well, valP1, valP2);
//				mphysics->InsertMaterialObject(BCLINE);	
//				
//				
//			}
//			
//			
//			//	// setting Computational element structure
//			//	std::set<int> iset;	
//			//	for(int imat = 0; imat < TotalMats; imat++)
//			//	{	
//			//		iset.insert(imat+1);	
//			//		mphysics->AutoBuild(iset);
//			//		mphysics->LoadReferences(); 		
//			//		mphysics->Reference()->ResetReference();			
//			//		iset.clear();
//			//		
//			//	}
//			////	
//			////	mphysics->Reference()->ResetReference();		
//			mphysics->AutoBuild();
//			mphysics->AdjustBoundaryElements();
//			mphysics->CleanUpUnconnectedNodes();
//			
//			// Creating multiphysic elements into mphysics computational mesh
//			TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//			TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//			TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
//			
//		#ifdef LOG4CXX
//			{
//				//        std::stringstream sout;
//				//        mphysics->Print(sout);
//				//        LOGPZ_DEBUG(logger, sout.str())
//			}
//		#endif
//			ofstream arg("Mphysic.txt");
//			mphysics->Print(arg);
//			
//			return mphysics;
//		}
//
//		#define VTK
//		void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
//		{			
//			//TPZBandStructMatrix full(fCmesh);
//			TPZSkylineStructMatrix full(fCmesh); //caso simetrico
//			an.SetStructuralMatrix(full);
//			TPZStepSolver<REAL> step;
//			step.SetDirect(ELDLt); //caso simetrico
//			//step.SetDirect(ELU);
//			an.SetSolver(step);
//			an.Run();
//			
//			//	ofstream file("Solution.out");
//			//	an.Solution().Print("solution", file); 
//		}
//
//
//
//		void PostProcessDiffusion(TPZAnalysis &an, std::string plotfile)
//		{
//			TPZManVector<std::string,10> scalnames(1), vecnames(1);
//			scalnames[0] = "Solution";
//			vecnames[0]= "MinusKGradU";
//			
//			const int dim = 2;
//			int div = 2;
//			an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//			an.PostProcess(div,dim);
//			std::ofstream out("malha.txt");
//			an.Print("nothing",out);
//		}
//
//		void PostProcessElasticity(TPZAnalysis &an, std::string plotfile)
//		{
//			TPZManVector<std::string,10> scalnames(3), vecnames(1);
//			scalnames[0] = "SigmaX";
//			scalnames[1] = "SigmaY";
//			scalnames[2] = "Pressure";	
//			vecnames[0]= "displacement";
//			
//			
//			const int dim = 2;
//			int div = 2;
//			an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//			an.PostProcess(div,dim);
//			std::ofstream out("malha.txt");
//			an.Print("nothing",out);
//		}
//
//		void PostProcessPoroeasticity(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
//		{
//			TPZBuildMultiphysicsMesh * Objectdumy;
//			Objectdumy->TransferFromMultiPhysics(meshvec, mphysics);
//			TPZManVector<std::string,10> scalnames(15), vecnames(3);
//			scalnames[0] = "SigmaX";
//			scalnames[1] = "SigmaY";
//			scalnames[2] = "ETAUXY";	
//			scalnames[3] = "TotalSigmaX";
//			scalnames[4] = "TotalSigmaY";	
//			scalnames[5] = "TauXY";	
//			scalnames[6] = "DisplacementX";
//			scalnames[7] = "DisplacementY";
//			scalnames[8] = "SolidPressure";
//			scalnames[9] = "FluidPressure";
//			scalnames[10] = "EDisplacementX";
//			scalnames[11] = "EDisplacementY";
//			scalnames[12] = "EPressure";
//			scalnames[13] = "ESIGX";
//			scalnames[14] = "ESIGY";
//			vecnames[0]  = "Displacement";
//			vecnames[1]  = "FluxVector";
//			vecnames[2]  = "EFluxVector";
//			
//			
//			
//			const int dim = 2;
//			int div = 2;
//			an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//			an.PostProcess(div,dim);
//			std::ofstream out("malha.txt");
//			an.Print("nothing",out);
//		}
//
//		void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh)
//		{
//			for ( int ref = 0; ref < nh; ref++ ){
//				TPZVec<TPZGeoEl *> filhos;
//				int n = gMesh->NElements();
//				for ( int i = 0; i < n; i++ ){
//					TPZGeoEl * gel = gMesh->ElementVec() [i];
//					if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
//				}//for i
//			}//ref
//		}
//
//		void RefinamentoUniforme(TPZGeoMesh * gMesh, int nh, int MatId, int indexEl)
//		{
//			for ( int ref = 0; ref < nh; ref++ ){
//				TPZVec<TPZGeoEl *> filhos;
//				int n = gMesh->NElements();
//				for ( int i = 0; i < n; i++ ){
//					TPZGeoEl * gel = gMesh->ElementVec() [i];
//					if (gel->Dimension() == 2 || gel->Dimension() == 1){
//						if (gel->MaterialId()== MatId && gel-> Index()==indexEl){
//							gel->Divide (filhos);
//						}
//					}
//				}//for i
//			}//ref
//		}
//
//		void RefinElemComp(TPZCompMesh  *cMesh, int indexEl)
//		{
//			
//			TPZVec<int > subindex; 
//			int nel = cMesh->ElementVec().NElements(); 
//			for(int el=0; el < nel; el++){
//				TPZCompEl * compEl = cMesh->ElementVec()[el];
//				if(!compEl) continue;
//				int ind = compEl->Index();
//				if(ind==indexEl){
//					compEl->Divide(indexEl, subindex, 1);
//				}
//			}	
//		}
//
//		void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
//		{
//			
//			TPZVec<int > subindex;
//			for (int iref = 0; iref < ndiv; iref++) {
//				TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
//				int nel = elvec.NElements(); 
//				for(int el=0; el < nel; el++){
//					TPZCompEl * compEl = elvec[el];
//					if(!compEl) continue;
//					int ind = compEl->Index();
//					compEl->Divide(ind, subindex, 0);
//				}
//			}
//		}
//
//		void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file)
//		{
//			file.clear();
//			int nelements = gmesh->NElements();
//			
//			std::stringstream node, connectivity, type;
//			
//			//Header
//			file << "# vtk DataFile Version 3.0" << std::endl;
//			file << "TPZGeoMesh VTK Visualization" << std::endl;
//			file << "ASCII" << std::endl << std::endl;
//			
//			file << "DATASET UNSTRUCTURED_GRID" << std::endl;
//			file << "POINTS ";
//			
//			int actualNode = -1, size = 0, nVALIDelements = 0;
//			
//			for(int el = 0; el < nelements; el++)
//			{        
//				if(gmesh->ElementVec()[el]->Type() == EPoint)//Exclude Lines and Arc3D
//				{
//					continue;
//				}
//				if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
//				{
//					continue;
//				}
//				if(gmesh->ElementVec()[el]->HasSubElement())
//				{
//					continue;
//				}
//				
//				int elNnodes = gmesh->ElementVec()[el]->NNodes();
//				size += (1+elNnodes);
//				connectivity << elNnodes;
//				
//				for(int t = 0; t < elNnodes; t++)
//				{
//					for(int c = 0; c < 3; c++)
//					{
//						double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
//						node << coord << " ";
//					}            
//					node << std::endl;
//					
//					actualNode++;
//					connectivity << " " << actualNode;
//				}
//				connectivity << std::endl;
//				
//				int elType = -1;
//				switch (gmesh->ElementVec()[el]->Type())
//				{
//					case (ETriangle):
//					{
//						elType = 5;
//						break;                
//					}
//					case (EQuadrilateral ):
//					{
//						elType = 9;
//						break;                
//					}
//					case (ETetraedro):
//					{
//						elType = 10;
//						break;                
//					}
//					case (EPiramide):
//					{
//						elType = 14;
//						break;                
//					}
//					case (EPrisma):
//					{
//						elType = 13;
//						break;                
//					}
//					case (ECube):
//					{
//						elType = 12;
//						break;                
//					}
//					default:
//					{
//						//ElementType NOT Found!!!
//						DebugStop();
//						break;    
//					}
//				}
//				
//				type << elType << std::endl;
//				nVALIDelements++;
//			}
//			node << std::endl;
//			actualNode++;
//			file << actualNode << " float" << std::endl << node.str();
//			
//			file << "CELLS " << nVALIDelements << " ";
//			
//			file << size << std::endl;
//			file << connectivity.str() << std::endl;
//			
//			file << "CELL_TYPES " << nVALIDelements << std::endl;
//			file << type.str();
//			
//			file.close();
//		}
//
//
//		void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file)
//		{
//			TPZGeoMesh *gmesh;
//			PrintGMeshVTK(gmesh, file);
//		}
//
//
//
//
//		TPZAutoPointer <TPZMatrix<REAL> > MassMatrix(TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics)
//		{
//			
//			TPZSpStructMatrix matsp(mphysics);
//			std::set< int > materialid;
//			
//			for(int imat = 0; imat < TotalMats; imat++)
//			{
//				mymaterial[imat]->SetLastState();
//				materialid.insert(imat+1);
//			}	
//			
//			matsp.SetMaterialIds(materialid);
//			TPZAutoPointer<TPZGuiInterface> guiInterface;
//			TPZFMatrix<REAL> Un;
//			TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
//			
//			return matK2;
//			
//		}
//
//
//		void StiffMatrixLoadVec(TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec)
//		{
//			
//			for(int imat = 0; imat < TotalMats; imat++)
//			{	
//				mymaterial[imat]->SetCurrentState();
//			}		
//			
//			//TPZFStructMatrix matsk(mphysics);
//			TPZSkylineStructMatrix matsk(mphysics);
//			an.SetStructuralMatrix(matsk); 
//			TPZStepSolver<REAL> step; 
//			step.SetDirect(ELDLt); 
//			an.SetSolver(step); 
//			an.Assemble(); 
//			matK1 = an.StructMatrix();
//			fvec = an.Rhs();
//			
//		}
//
//
//
//		void SolveSistTransient(TPZVec <REAL> &PrintEachtimeStep, std::string FileName, REAL deltaT,REAL maxTime, TPZFMatrix<REAL> &InitialSolution,TPZVec <TPZPoroElastic2d  * > &mymaterial, 
//								TPZAnalysis &an, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics){
//			
//			// Calculationg Mass Matrix    
//			TPZAutoPointer <TPZMatrix<REAL> > PoroelasticMassMatrix = MassMatrix(mymaterial, mphysics);
//			
//			// Calculationg Stiffness Matrix and Load Vector
//			TPZFMatrix<REAL> PoroelasticStiffnessMatrix;	
//			TPZFMatrix<REAL> PoroelasticLoadVector; 
//			StiffMatrixLoadVec(mymaterial, mphysics, an, PoroelasticStiffnessMatrix, PoroelasticLoadVector);
//			
//			// Printing Matrix
//		#ifdef LOG4CXX
//			if(logdata->isDebugEnabled())
//			{
//				std::stringstream sout;
//				PoroelasticStiffnessMatrix.Print("KStiffness = ", sout,EMathematicaInput);
//				PoroelasticLoadVector.Print("LoadVec = ", sout,EMathematicaInput);
//				PoroelasticMassMatrix->Print("MassMatrix = ", sout,EMathematicaInput);		
//				InitialSolution.Print("IntialConditions = ", sout,EMathematicaInput);
//				LOGPZ_DEBUG(logdata,sout.str())
//			}
//			
//		#endif
//			
//			int nrows;
//			nrows = PoroelasticMassMatrix->Rows();
//			TPZFMatrix<REAL> TotalRhs(nrows,1,0.0);
//			TPZFMatrix<REAL> RhsTemporal(nrows,1,0.0);
//			TPZFMatrix<REAL> LastSolution = InitialSolution;
//			
//			
//			
//			REAL	TimeValue	= 0.0;
//			int		cent		= 1;
//			int		control		= 0;
//			
//			TimeValue = cent*deltaT; 
//			
//			while (TimeValue <= maxTime)
//			{	
//				
//				for(int imat = 0; imat < TotalMats; imat++)
//				{
//					mymaterial[imat]->SetTimeValue(TimeValue);
//				}
//				PoroelasticMassMatrix->Multiply(LastSolution,RhsTemporal);
//				
//				//#ifdef LOG4CXX
//				//        if(logdata->isDebugEnabled())
//				//        {
//				//            std::stringstream sout;
//				//            sout<< " Current Time = " << TimeValue;
//				//            LastSolution.Print("CurrenState = ", sout,EMathematicaInput);
//				//            RhsTemporal.Print("Mat Mass x Last solution = ", sout,EMathematicaInput);	
//				//            LOGPZ_DEBUG(logdata,sout.str())
//				//        }
//				//#endif
//				
//				TotalRhs = PoroelasticLoadVector + RhsTemporal;
//				an.Rhs() = TotalRhs;
//				an.Solve(); 
//				LastSolution = an.Solution();
//				
//				
//				// Information Print Control
//				
//				if (control == PrintEachtimeStep.size()) 
//				{
//					control = 0;
//				}				
//				if (abs(PrintEachtimeStep[control] - TimeValue) <= 1.0e-10) 
//				{
//					std::stringstream outputfiletemp;
//					outputfiletemp << FileName << ".vtk";
//					std::string plotfile = outputfiletemp.str();
//					PostProcessPoroeasticity(meshvec,mphysics,an,plotfile);
//					control++;
//					
//					//		// Total mass calculation
//					//		
//					//		int totalel=mphysics->NElements();
//					//		int vartointegrate = 13;
//					//		REAL TotalMass=0;
//					//		
//					//		for ( int el = 0; el < totalel; el++ )
//					//		{
//					//			TPZVec <REAL> result(1,0);
//					//			TPZCompEl * celvar = mphysics->ElementVec()[el];
//					//			if (celvar->Reference()->MaterialId() > 0 )// avoid bc conditions
//					//			{
//					//				celvar->Integrate(vartointegrate, result);
//					//				TotalMass += result[0];	
//					//			}
//					//		}
//					//		cout << "Numerical total mass  " << TotalMass << "  at time " << TimeValue << endl;
//					//		cout << "Theorical total mass  " << 0.002*TimeValue << "  at time " << TimeValue << endl;		
//					//		
//					//		// End TotalMass Calculation			
//					
//				}
//				
//				
//				
//				cent++;
//				TimeValue = cent*deltaT;
//			}
//			
//		}
//
//		void *InitalPressureCalculations(REAL &PressureValue,TPZVec <REAL> &PointCoordinates)
//		{
//			PressureValue = 0.0;
//		}
//
//		void SetInitialConditions(TPZAnalysis &AnalysisToInitialize, TPZCompMesh * ComputationalMesh, TPZFMatrix<REAL> &ComputationalMeshInitialSolution, void ( *PointerToFunction(REAL &,TPZVec <REAL> &)))
//		{
//			
//			int RowsNumber = AnalysisToInitialize.Solution().Rows();
//			ComputationalMeshInitialSolution.Resize(RowsNumber,1);
//			
//			for (int i = 0; i < RowsNumber; i++) 
//			{
//				ComputationalMeshInitialSolution(i,0) = 0.0;
//			}
//			
//			int TotalComputationalElements	= ComputationalMesh->NElements();		
//			for ( int el = 0; el < TotalComputationalElements ; el++ )
//			{
//				
//				TPZCompEl * ComputationalElement = ComputationalMesh->ElementVec()[el];
//				TPZGeoEl * GeoemtricElement;
//				GeoemtricElement = ComputationalElement->Reference();
//				int NumberNodeCorners = ComputationalElement->Reference()->NCornerNodes();
//				
//				for (int ConnectIndex = 0; ConnectIndex < NumberNodeCorners ; ConnectIndex++) 
//				{
//					TPZConnect & ThisConnect = ComputationalElement->Connect(ConnectIndex);
//					TPZVec <REAL> Coordenates(3,0.0);
//					REAL InitialValue;			
//					Coordenates[0] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(0);
//					Coordenates[1] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(1);
//					Coordenates[2] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(2);
//					
//					PointerToFunction(InitialValue, Coordenates);			
//					
//					int Sequence = ThisConnect.SequenceNumber();
//					int BlockPos = ComputationalMesh->Block().Position(Sequence);
//					
//					ComputationalMeshInitialSolution(BlockPos,0)= InitialValue;
//				}
//			}	
//		}
//
//		TPZCompMesh *SetInitialConditionsbyL2Projection(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &InitialCondition, bool &CustomFunction)
//		{
//			
//			int Dimension = 2, NumberStateVariables = 1;
//			TPZL2Projection * MyProjection;
//			//	This create material for projection	
//			MyProjection = new TPZL2Projection(matId, Dimension, NumberStateVariables, InitialCondition, pOrder);
//			
//			TPZCompMesh * ComputationalMeshMyProjection = new TPZCompMesh(gmesh);
//			ComputationalMeshMyProjection->SetDimModel(Dimension);
//			ComputationalMeshMyProjection->SetAllCreateFunctionsContinuous();
//			//	ComputationalMeshMyProjection->SetAllCreateFunctionsDiscontinuous();	
//			ComputationalMeshMyProjection->SetDefaultOrder(pOrder);
//			
//			
//			TPZMaterial * mat(MyProjection);
//			if (CustomFunction) {
//				TPZAutoPointer < TPZFunction<STATE> > InitialPressureCalculations;
//				InitialPressureCalculations = new TPZDummyFunction<STATE> (InitialPressureDistribution);
//				MyProjection->SetForcingFunction(InitialPressureCalculations);	
//			}
//			
//			ComputationalMeshMyProjection->InsertMaterialObject(mat);
//			
//			
//			
//			//    ///inserir connect da pressao
//			//    int ncon = cmesh->NConnects();
//			//    for(int i=0; i<ncon; i++)
//			//    {
//			//        TPZConnect &newnod = cmesh->ConnectVec()[i];
//			//        newnod.SetPressure(true);
//			//    }
//			//    
//			//    ///set order total da shape
//			//    int nel = cmesh->NElements();
//			//    for(int i=0; i<nel; i++){
//			//        TPZCompEl *cel = cmesh->ElementVec()[i];
//			//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//			//        celdisc->SetConstC(1.);
//			//        celdisc->SetCenterPoint(0, 0.);
//			//        celdisc->SetCenterPoint(1, 0.);
//			//        celdisc->SetCenterPoint(2, 0.);
//			//        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//			//        {
//			//            if(triang==true) celdisc->SetTotalOrderShape();
//			//            else celdisc->SetTensorialShape();
//			//        }
//			//    }
//			//    
//			//    
//			//#ifdef DEBUG   
//			//    int ncel = cmesh->NElements();
//			//    for(int i =0; i<ncel; i++){
//			//        TPZCompEl * compEl = cmesh->ElementVec()[i];
//			//        if(!compEl) continue;
//			//        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
//			//        if(facel)DebugStop();
//			//        
//			//    }
//			//#endif
//			ComputationalMeshMyProjection->AutoBuild();    
//			return ComputationalMeshMyProjection;
//		}
//
//
//

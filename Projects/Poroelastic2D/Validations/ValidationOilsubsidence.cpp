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
//		const int TotalMats = 2; 
//		const int matId = 1;	
//
//		// Boundary Conditions Materials
//		// 2D Rectangular Domain for sideburden
//		const int bcBottomSB = 7;
//		const int bcRightSB = 10;
//		const int bcTopSB = 9;
//		const int bcLeftSB = 8;
//
//		// 2D Rectangular Domain for sideburden
//		const int bcBottomR = 3;
//		const int bcRightR = 6;
//		const int bcTopR = 5;
//		const int bcLeftR = 4;
//
//		// line source identity
//		const int WellLine = 20;
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
//		TPZCompMesh *ComputationalPoroelasticityInitialMesh(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial);
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
//		void SolveSistTransient(TPZVec <REAL> &PrintEachtimeStep, std::string FileName, REAL deltaT,REAL maxTime, TPZFMatrix<REAL> &InitialSolution, TPZVec <TPZPoroElastic2d  * > &mymaterial , TPZAnalysis &an, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);	
//		void StiffMatrixLoadVec(TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec);
//
//		//	This are tools for spatial and polinomyal refinement and Postprocess of solutions 
//		void PostProcessElasticity(TPZAnalysis &an, std::string plotfile);
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
//			//	gRefDBase.InitializeRefPatterns();
//			//	gRefDBase.InitializeAllUniformRefPatterns();
//			//	ofstream RefinElPatterns("RefElPatterns.txt");
//			//	gRefDBase.WriteRefPatternDBase(RefinElPatterns); 	
//			//	gRefDBase.ReadRefPatternDBase("RefElPatterns.txt");
//			
//			TPZReadGIDGrid mygenerator;
//			mygenerator.SetfDimensionlessL(1.0);
//			//	TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh("CirclePlanStrain.dump");
//			//	TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh("OilExtractionsubsidencemoreCoarse.dump");	
//			//	TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh("Segall1985GoodAproxx.dump");
//			//	TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh("Segall1985.dump");
//			//	TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh("Segall1985extended.dump");
//			TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh("Segall1985extendedCoarse.dump");	
//			//	TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh("Segall1985Well.dump");
//			//	TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh("Segall1985WellSimple.dump");	
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
//			//	// Directional Point source Refinement 
//			//	int ndirectdivp = 17;
//			//	set<int> SETmatPointRefDir2;
//			//	SETmatPointRefDir2.insert(WellLine);
//			//	//	SETmatPointRefDir2.insert(bcBottom);
//			//	//	SETmatPointRefDir2.insert(bcTop);
//			//	//	SETmatPointRefDir2.insert(bcRight);
//			//	//	SETmatPointRefDir2.insert(bcLeft);	
//			//	//	SETmatPointRefDir2.insert(xfixedPoints);
//			//	//	SETmatPointRefDir2.insert(yfixedPoints);	
//			//	for(int j = 0; j < ndirectdivp; j++)
//			//	{
//			//		int nel = gmesh->NElements();
//			//		for (int iref = 0; iref < nel; iref++)
//			//		{
//			//			TPZVec<TPZGeoEl*> filhos;
//			//			TPZGeoEl * gelP2 = gmesh->ElementVec()[iref];
//			//			if(!gelP2 || gelP2->HasSubElement()) continue;
//			//			TPZRefPatternTools::RefineDirectional(gelP2, SETmatPointRefDir2);
//			//		}		
//			//	}	
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
//			//--------------------------------------------------------------------------------------------------------------------------------//	
//			//--------------------------------------------------------------------------------------------------------------------------------//	
//			
//			// Setting intial conditions
//			
//			//	//	Set initial conditions for pressure
//			//	TPZFMatrix<REAL> DiffusionSolution;
//			//	TPZAnalysis ElasticAnalysis(ComputationalMeshElasticity);
//			//	TPZAnalysis DiffusionAnalysis(ComputationalMeshDiffusion);	
//			//	//	//	//	Solving Second Computational problem
//			//	//	SolveSist(DiffusionAnalysis, ComputationalMeshDiffusion);
//			//	//	std::string plotfile("SolutionDiffusion.vtk");		
//			//	//	PostProcessDiffusion(DiffusionAnalysis, plotfile);	
//			//	
//			//	SolveSist(ElasticAnalysis,ComputationalMeshElasticity);		
//			//	std::string OutTubeElasticity("SolutionElasticity.vtk");		
//			//	PostProcessElasticity(ElasticAnalysis, OutTubeElasticity);		
//			//	
//			//	SolveSist(DiffusionAnalysis,ComputationalMeshDiffusion);
//			//	std::string OutTubeDiffusion("SolutionDiffusion.vtk");	
//			//	PostProcessDiffusion(DiffusionAnalysis, OutTubeDiffusion);		
//			//	
//			//	
//			//	Without L2 projection
//			//	SetInitialConditions(DiffusionAnalysis, ComputationalMeshDiffusion, DiffusionSolution, *InitalPressureCalculations);
//			//	DiffusionAnalysis.LoadSolution(DiffusionSolution);
//			
//			//	Set initial conditions for pressure by using L2 projection
//			//	int nrs = DiffusionAnalysis.Solution().Rows();
//			//	bool CustomFunction = false;
//			//	//	Using  constant initial pressure
//			//    TPZVec<REAL> InitialCondition(nrs,1000.0);
//			//    TPZCompMesh  * ComputationalL2 = SetInitialConditionsbyL2Projection(gmesh, p, InitialCondition, CustomFunction);
//			//    TPZAnalysis L2Analysis(ComputationalL2);
//			//    SolveSist(L2Analysis, ComputationalL2);
//			//    DiffusionAnalysis.LoadSolution(L2Analysis.Solution());	
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
//			TPZVec<TPZCompMesh *> ComputationalMeshVectorInitial(2);
//			ComputationalMeshVectorInitial[0] = ComputationalMeshElasticity;
//			ComputationalMeshVectorInitial[1] = ComputationalMeshDiffusion;
//			
//			TPZVec<TPZCompMesh *> ComputationalMeshVector(2);
//			ComputationalMeshVector[0] = ComputationalMeshElasticity;
//			ComputationalMeshVector[1] = ComputationalMeshDiffusion;	
//			
//			TPZVec <TPZPoroElastic2d *>  materialistini(TotalMats,0) ;
//			TPZVec <TPZPoroElastic2d *>  materialist(TotalMats,0) ;	
//			
//			TPZCompMesh * ComputationalMeshPoroelasticityInitial = ComputationalPoroelasticityInitialMesh(gmesh,ComputationalMeshVectorInitial,materialistini);
//			TPZCompMesh * ComputationalMeshPoroelasticityReCurrent = ComputationalPoroelasticityMesh(gmesh,ComputationalMeshVector,materialist);	
//			
//			//--------------------------------------------------------------------------------------------------------------------------------//	
//			//--------------------------------------------------------------------------------------------------------------------------------//	
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
//			REAL Delta		= 0.1*year;
//			REAL MaxTime	= 100.0*year;	
//			//	REAL Delta = 20;
//			//	REAL MaxTime = 1000.0;
//			std::string	OutFileName("PoroelasticSolution");
//			//	REAL Theta = 0.6667;
//			REAL Theta = 1.0;	
//			//	REAL Theta = 1+((1/Delta)-(1/(log(1+Delta))));
//			//	int PrintEachtimeStep = 999;
//			
//			TPZVec <REAL> OneStepToInitial(5,0);
//			OneStepToInitial[0] = 20.0*year;	
//			OneStepToInitial[1] = 40.0*year;	
//			OneStepToInitial[2] = 60.0*year;
//			OneStepToInitial[3] = 80.0*year;	
//			OneStepToInitial[4] = 100.0*year;		
//			
//			TPZVec <REAL> PrintEachtimeStep(6,0);
//			PrintEachtimeStep[0] = 0.2*year;	
//			PrintEachtimeStep[1] = 20.0*year;	
//			PrintEachtimeStep[2] = 40.0*year;	
//			PrintEachtimeStep[3] = 60.0*year;
//			PrintEachtimeStep[4] = 80.0*year;	
//			PrintEachtimeStep[5] = 100.0*year;	
//			
//			// Improve this
//			for(int imat = 0; imat < TotalMats; imat++)
//			{	
//				materialistini[imat]->SetTimeStep(Delta,Theta);
//				materialistini[imat]->SetTimeValue(InitialTime);		
//				materialist[imat]->SetTimeStep(Delta,Theta);
//				materialist[imat]->SetTimeValue(InitialTime);			
//			}
//			
//			TPZAnalysis PoroelasticAnalysisInitial(ComputationalMeshPoroelasticityInitial);	
//			TPZAnalysis PoroelasticAnalysis(ComputationalMeshPoroelasticityReCurrent);
//			
//			//Setting initial coditions 
//			TPZFMatrix<REAL> TrialInitialSolution = PoroelasticAnalysisInitial.Solution();
//			std::string output("PoroelasticSolutionInitial");
//			std::stringstream outputfiletemp;
//			outputfiletemp << output << ".vtk";
//			std::string plotfile = outputfiletemp.str();
//			
//			//	Calculation of Initial conditions
//			SolveSistTransient(OneStepToInitial,output,Delta,MaxTime, TrialInitialSolution, materialistini, PoroelasticAnalysisInitial, ComputationalMeshVectorInitial,ComputationalMeshPoroelasticityInitial);
//			//Setting initial coditions 
//			PoroelasticAnalysis.LoadSolution(PoroelasticAnalysisInitial.Solution());	
//			TPZFMatrix<REAL> InitialSolution = PoroelasticAnalysis.Solution();	
//			
//			//	//	Print Initial conditions	
//			//	PostProcessPoroeasticity(ComputationalMeshVectorInitial,ComputationalMeshPoroelasticityInitial,PoroelasticAnalysisInitial,plotfile);	
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
//				PostProcessPoroeasticity(ComputationalMeshVector,ComputationalMeshPoroelasticityReCurrent,PoroelasticAnalysis,plotfile);
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
//			SolveSistTransient(PrintEachtimeStep,OutFileName,Delta,MaxTime, InitialSolution, materialist, PoroelasticAnalysis, ComputationalMeshVector,ComputationalMeshPoroelasticityReCurrent);	
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
//			REAL gravity = 9.81; // SI system
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
//			for(int imat = 0; imat < TotalMats-1; imat++)
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
//				TPZMaterial * BCTOP = materialist[imat]->CreateBC(mat, bcTopSB,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCTOP);
//				
//				// Bottom Boundary			
//				TPZMaterial * BCBOTT = materialist[imat]->CreateBC(mat, bcBottomSB,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCBOTT);
//				
//				// Left Boundary			
//				TPZMaterial * BCLEFT = materialist[imat]->CreateBC(mat, bcLeftSB,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary			
//				TPZMaterial * BCRIGHT = materialist[imat]->CreateBC(mat, bcRightSB,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCRIGHT);
//				
//				//		// Line Source
//				//		TPZMaterial * BCLINE = materialist[imat]->CreateBC(mat, WellLine,neumann, val1, val2);
//				//		cmesh->InsertMaterialObject(BCLINE);
//				
//			}
//			
//			for(int imat = 1; imat < TotalMats; imat++)
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
//				TPZMaterial * BCTOP = materialist[imat]->CreateBC(mat, bcTopR,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCTOP);
//				
//				// Bottom Boundary			
//				TPZMaterial * BCBOTT = materialist[imat]->CreateBC(mat, bcBottomR,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCBOTT);
//				
//				// Left Boundary			
//				TPZMaterial * BCLEFT = materialist[imat]->CreateBC(mat, bcLeftR,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary			
//				TPZMaterial * BCRIGHT = materialist[imat]->CreateBC(mat, bcRightR,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCRIGHT);
//				
//				// Line Source
//				TPZMaterial * BCLINE = materialist[imat]->CreateBC(mat, WellLine,neumann, val1, val2);
//				cmesh->InsertMaterialObject(BCLINE);
//				
//			}	
//			
//			
//			if (Initialstress) 
//			{
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
//			
//			for(int imat = 0; imat < TotalMats-1; imat++)
//			{
//				//		// Creating and inserting material objec in to computational mesh		
//				//		materialist[imat] = new TPZMatPoisson3d(imat+1,dim); 
//				//		TPZMaterial * mat(materialist[imat]);
//				//		materialist[imat]->SetParameters(diff[imat], conv[imat], convdirVec[imat]);
//				//		materialist[imat]->SetInternalFlux(flux[imat]);
//				//		materialist[imat]->NStateVariables();
//				//		cmesh->InsertMaterialObject(mat);
//				//		
//				//		// This is a dummy boundary conditions
//				//		// Here is necessary to include the same set of computational information for all computational meshes -> Elastic, Pressure, Biot's Poroelasticity
//				//		TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
//				//		
//				//		// Top Boundary			
//				//		TPZMaterial * BCTOP = materialist[imat]->CreateBC(mat, bcTopSB,dirichlet, val1, val2);
//				//		cmesh->InsertMaterialObject(BCTOP);
//				//		
//				//		// Bottom Boundary			
//				//		TPZMaterial * BCBOTT = materialist[imat]->CreateBC(mat, bcBottomSB,dirichlet, val1, val2);
//				//		cmesh->InsertMaterialObject(BCBOTT);
//				//		
//				//		// Left Boundary			
//				//		TPZMaterial * BCLEFT = materialist[imat]->CreateBC(mat, bcLeftSB,dirichlet, val1, val2);
//				//		cmesh->InsertMaterialObject(BCLEFT);
//				//		
//				//		// Right Boundary			
//				//		TPZMaterial * BCRIGHT = materialist[imat]->CreateBC(mat, bcRightSB,dirichlet, val1, val2);
//				//		cmesh->InsertMaterialObject(BCRIGHT);
//				
//			}	
//			
//			
//			for(int imat = 1; imat < TotalMats; imat++)
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
//				TPZMaterial * BCTOP = materialist[imat]->CreateBC(mat, bcTopR,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCTOP);
//				
//				// Bottom Boundary			
//				TPZMaterial * BCBOTT = materialist[imat]->CreateBC(mat, bcBottomR,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCBOTT);
//				
//				// Left Boundary			
//				TPZMaterial * BCLEFT = materialist[imat]->CreateBC(mat, bcLeftR,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary			
//				TPZMaterial * BCRIGHT = materialist[imat]->CreateBC(mat, bcRightR,dirichlet, val1, val2);
//				cmesh->InsertMaterialObject(BCRIGHT);
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
//			
//			cmesh->AutoBuild();
//			
//			return cmesh;
//		}
//
//		TPZCompMesh * ComputationalPoroelasticityInitialMesh(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial)
//		{
//			//Creating computational mesh
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
//			REAL lamb = 5.33333e9;							//	[Pa]
//			REAL lambu = 1.55294e10;					//	[Pa]
//			REAL alpha = 0.814536;						//	[-]
//			REAL G= 8.0e9;							//	[Pa]
//			REAL Phi = 0.3;
//			REAL rhof = 1000.0;						//	[kg/m3]	
//			REAL S = ((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
//			REAL Se = ((pow(alpha,2))/((lambu-lamb)));		
//			REAL Coeff = (alpha)/((lamb+2.0*G));
//			REAL rockrho = 2300.0;					//	[kg/m3]
//			REAL rho = (1-Phi)*rockrho+Phi*rhof;
//			REAL gravity = 9.81;					//	[m/s2]
//			REAL c = 0.1;							//	[m2/s]	
//			REAL visc = 0.001;						//	[Pa.s]
//			REAL perm = c*S*visc;//7.8784288e-15;//1.109542e-14						//	[m2]
//			REAL viscd = 1.0;						//	[Pa.s]
//			REAL permd = 1.0;//7.8784288e-15;//1.109542e-14						//	[m2]	
//			REAL PI = atan(1.)*4.;
//			REAL qo = 0.0;
//			TPZVec <REAL> force(2,0.);	
//			
//			
//			// Definitions for paramter vectors Al materials have the same porperties
//			TPZVec <REAL> lambVec(TotalMats,lamb);
//			TPZVec <REAL> lambuVec(TotalMats,lambu);
//			TPZVec <REAL> alphaVec(TotalMats,alpha);
//			TPZVec <REAL> GVec(TotalMats,G);
//			TPZVec <REAL> rhofVec(TotalMats,rhof);
//			//	TPZVec <REAL> poissonVec(TotalMats,poisson);
//			//	TPZVec <REAL> EyoungVec(TotalMats,Eyoung);
//			TPZVec <REAL> SVec(TotalMats,S);
//			TPZVec <REAL> SeVec(TotalMats,Se);
//			TPZVec <REAL> rockrhoVec(TotalMats,rockrho);
//			TPZVec <REAL> cVec(TotalMats,c);
//			TPZVec <REAL> viscVec(TotalMats,visc);
//			TPZVec <REAL> permVec(TotalMats,perm);
//			TPZVec <TPZVec< REAL > > forceVec(TotalMats,force);
//			
//			for(int imat = 0; imat < TotalMats-1; imat++)
//			{	
//				
//				// No Gravity field
//				forceVec[imat][1]= - gravity*rockrhoVec[imat];	
//				mymaterial[imat] = new TPZPoroElastic2d (imat+1, dim);
//				mymaterial[imat]->SetParameters(lambVec[imat], GVec[imat],forceVec[imat][0], forceVec[imat][1]);
//				mymaterial[imat]->SetParameters(0.0*permVec[imat],viscVec[imat]);
//				mymaterial[imat]->SetfPlaneProblem(planestress);
//				mymaterial[imat]->SetBiotParameters(0.0*alpha,0.0*Se);
//				mymaterial[imat]->SetTimeDependentFunctionExact(TimeDepFExact);		
//				
//				ofstream argm("mymaterial1.txt");
//				mymaterial[imat]->Print(argm);
//				TPZMaterial * mat(mymaterial[imat]);
//				mphysics->InsertMaterialObject(mat);
//				
//				///--- --- Inserting Boundary Conditions
//				// Squared Model for validation line source Validation give in ...........					
//				///--- --- Boundary Conditions
//				// Squared Model 1D compactation
//				
//				// Elastic problem -> 1 ,  Diffsion problem -> 2
//				
//				// Top Boundary	-3
//				// Setting free Boundary condition N1N2 means: Elastic -> Neumman, Diffusion -> Dirichlet
//				int N1D2 = 110;
//				TPZFMatrix<REAL> val13(3,2,0.), val23(3,1,0.);
//				REAL SigTopx	= 0.0;
//				REAL SigTopy	= 0.0;			
//				REAL pTop		= 0.0;
//				val23(0,0)=SigTopx;
//				val23(1,0)=SigTopy;			
//				val23(2,0)=pTop;
//				TPZMaterial * BCTOP = mymaterial[imat]->CreateBC(mat,bcTopSB, N1D2, val13, val23);
//				mphysics->InsertMaterialObject(BCTOP);			
//				
//				// Bottom Boundary -1					
//				/// Setting free displacement non permeable surface Boundary condition DYFX1N2 means: Elastic -> Dirichlet Y component, Diffusion -> Neumman
//				TPZFMatrix<REAL> val11(3,2,0.), val21(3,1,0.);
//				int DYFX1N2 = 100; 			
//				REAL uDBottx	=0.0;
//				REAL uDBotty	=0.0;
//				REAL fluxBott	=0.0;
//				val21(0,0)=uDBottx;
//				val21(1,0)=uDBotty;
//				val21(2,0)=fluxBott;
//				TPZMaterial * BCBOTT = mymaterial[imat]->CreateBC(mat,bcBottomSB, DYFX1N2, val11, val21);
//				mphysics->InsertMaterialObject(BCBOTT);
//				
//				// Letf Boundary -4		
//				/// Setting free displacement non permeable surface Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman
//				int DXFY1N2 = 010;
//				TPZFMatrix<REAL> val14(3,2,0.), val24(3,1,0.);
//				REAL uDLeftx	= 0.0;
//				REAL uDLefty	= 0.0;			
//				REAL fluxLeft	= 0.0;
//				val24(0,0)=uDLeftx;
//				val24(1,0)=uDLefty;			
//				val24(2,0)=fluxLeft;
//				TPZMaterial * BCLEFT = mymaterial[imat]->CreateBC(mat, bcLeftSB, DXFY1N2, val14, val24);
//				mphysics->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary -2			
//				/// Setting free displacement non permeable surface  Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman			
//				TPZFMatrix<REAL> val12(3,2,0.), val22(3,1,0.);
//				
//				REAL uDRightx	= 0.0;
//				REAL uDRighty	= 0.0;			
//				REAL fluxRight	= 0.0;
//				val22(0,0)=uDRightx;
//				val22(1,0)=uDRighty;			
//				val22(2,0)=fluxRight;
//				TPZMaterial * BCRIGHT = mymaterial[imat]->CreateBC(mat, bcRightSB, DXFY1N2, val12, val22);
//				mphysics->InsertMaterialObject(BCRIGHT);
//				
//			}
//			
//			
//			for(int imat = 1; imat < TotalMats; imat++)
//			{	
//				
//				// No Gravity field
//				forceVec[imat][1]= - gravity*rho;
//				mymaterial[imat] = new TPZPoroElastic2d (imat+1, dim);
//				mymaterial[imat]->SetParameters(lambVec[imat], GVec[imat],forceVec[imat][0], forceVec[imat][1]);
//				mymaterial[imat]->SetParameters(permVec[imat],viscVec[imat]);
//				mymaterial[imat]->SetfPlaneProblem(planestress);
//				mymaterial[imat]->SetBiotParameters(alpha,0.0*Se);
//				mymaterial[imat]->SetTimeDependentFunctionExact(TimeDepFExact);		
//				
//				ofstream argm("mymaterial1.txt");
//				mymaterial[imat]->Print(argm);
//				TPZMaterial * mat(mymaterial[imat]);
//				mphysics->InsertMaterialObject(mat);
//				
//				///--- --- Inserting Boundary Conditions
//				// Squared Model for validation line source Validation give in ...........					
//				///--- --- Boundary Conditions
//				// Squared Model 1D compactation
//				
//				// Elastic problem -> 1 ,  Diffsion problem -> 2
//				
//				// Top Boundary	-3
//				// Setting free Boundary condition N1N2 means: Elastic -> Neumman, Diffusion -> Dirichlet
//				int N1D2 = 111;
//				TPZFMatrix<REAL> val13(3,2,0.), val23(3,1,0.);
//				REAL SigTopx	= 0.0;
//				REAL SigTopy	= 0.0;			
//				REAL pTop		= 0.0;
//				val23(0,0)=SigTopx;
//				val23(1,0)=SigTopy;			
//				val23(2,0)=pTop;
//				TPZMaterial * BCTOP = mymaterial[imat]->CreateBC(mat,bcTopR, N1D2, val13, val23);
//				mphysics->InsertMaterialObject(BCTOP);			
//				
//				// Bottom Boundary -1					
//				/// Setting free displacement non permeable surface Boundary condition DYFX1N2 means: Elastic -> Dirichlet Y component, Diffusion -> Neumman
//				TPZFMatrix<REAL> val11(3,2,0.), val21(3,1,0.);
//				int DYFX1N2 = 111; 			
//				REAL uDBottx	=0.0;
//				REAL uDBotty	=0.0;
//				REAL fluxBott	=0.0;
//				val21(0,0)=uDBottx;
//				val21(1,0)=uDBotty;
//				val21(2,0)=fluxBott;
//				TPZMaterial * BCBOTT = mymaterial[imat]->CreateBC(mat,bcBottomR, DYFX1N2, val11, val21);
//				mphysics->InsertMaterialObject(BCBOTT);
//				
//				// Letf Boundary -4		
//				/// Setting free displacement non permeable surface Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman
//				int DXFY1N2 = 011;
//				TPZFMatrix<REAL> val14(3,2,0.), val24(3,1,0.);
//				REAL uDLeftx	= 0.0;
//				REAL uDLefty	= 0.0;			
//				REAL fluxLeft	= 0.0;
//				val24(0,0)=uDLeftx;
//				val24(1,0)=uDLefty;			
//				val24(2,0)=fluxLeft;
//				TPZMaterial * BCLEFT = mymaterial[imat]->CreateBC(mat, bcLeftR, DXFY1N2, val14, val24);
//				mphysics->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary -2			
//				/// Setting free displacement non permeable surface  Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman			
//				TPZFMatrix<REAL> val12(3,2,0.), val22(3,1,0.);
//				
//				REAL uDRightx	= 0.0;
//				REAL uDRighty	= 0.0;			
//				REAL fluxRight	= 0.0;
//				val22(0,0)=uDRightx;
//				val22(1,0)=uDRighty;			
//				val22(2,0)=fluxRight;
//				TPZMaterial * BCRIGHT = mymaterial[imat]->CreateBC(mat, bcRightR, DXFY1N2, val12, val22);
//				mphysics->InsertMaterialObject(BCRIGHT);
//				
//				//		//  Source Production/injection Well line (Q is the injected/produced volume per time [kg/s])
//				//		int Well = 011;
//				//		TPZFMatrix<REAL> valP1(3,2,0.), valP2(3,1,0.);
//				//		valP2(2,0) = qo/(rhof);
//				//		TPZMaterial * BCLINE = mymaterial[imat]->CreateBC(mat, WellLine, Well, valP1, valP2);
//				//		mphysics->InsertMaterialObject(BCLINE);			
//				
//				
//			}
//			
//			
//			
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
//			REAL lamb = 5.33333e9;							//	[Pa]
//			REAL lambu = 1.55294e10;					//	[Pa]
//			REAL alpha = 0.814536;						//	[-]
//			REAL G= 8.0e9;							//	[Pa]
//			REAL Phi = 0.3;
//			REAL rhof = 1000.0;						//	[kg/m3]	
//			REAL S = ((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
//			REAL Se = ((pow(alpha,2))/((lambu-lamb)));		
//			REAL Coeff = (alpha)/((lamb+2.0*G));
//			REAL rockrho = 2300.0;					//	[kg/m3]
//			REAL rho = (1-Phi)*rockrho+Phi*rhof;
//			REAL gravity = 9.81;					//	[m/s2]
//			REAL c = 0.1;							//	[m2/s]	
//			REAL visc = 0.001;						//	[Pa.s]
//			REAL perm = c*S*visc;//7.8784288e-15;//1.109542e-14						//	[m2]
//			REAL viscd = 1.0;						//	[Pa.s]
//			REAL permd = 1.0;//7.8784288e-15;//1.109542e-14						//	[m2]	
//			REAL PI = atan(1.)*4.;
//			REAL Bo = 1.0;
//			REAL qo = 0.00063419*(1+0.4);//0.00200063419;
//			TPZVec <REAL> force(2,0.);		
//			
//			
//			// Definitions for paramter vectors Al materials have the same porperties
//			TPZVec <REAL> lambVec(TotalMats,lamb);
//			TPZVec <REAL> lambuVec(TotalMats,lambu);
//			TPZVec <REAL> alphaVec(TotalMats,alpha);
//			TPZVec <REAL> GVec(TotalMats,G);
//			TPZVec <REAL> rhofVec(TotalMats,rhof);
//			//	TPZVec <REAL> poissonVec(TotalMats,poisson);
//			//	TPZVec <REAL> EyoungVec(TotalMats,Eyoung);
//			TPZVec <REAL> SVec(TotalMats,S);
//			TPZVec <REAL> SeVec(TotalMats,Se);
//			TPZVec <REAL> rockrhoVec(TotalMats,rockrho);
//			TPZVec <REAL> cVec(TotalMats,c);
//			TPZVec <REAL> viscVec(TotalMats,visc);
//			TPZVec <REAL> permVec(TotalMats,perm);
//			
//			TPZVec <TPZVec< REAL > > forceVec(TotalMats,force);
//			
//			for(int imat = 0; imat < TotalMats-1; imat++)
//			{	
//				
//				// No Gravity field
//				forceVec[imat][1]= - gravity*rockrhoVec[imat];	
//				mymaterial[imat] = new TPZPoroElastic2d (imat+1, dim);
//				mymaterial[imat]->SetParameters(lambVec[imat], GVec[imat],forceVec[imat][0], forceVec[imat][1]);
//				mymaterial[imat]->SetParameters(0.0*permVec[imat],viscVec[imat]);
//				mymaterial[imat]->SetfPlaneProblem(planestress);
//				mymaterial[imat]->SetBiotParameters(0.0*alpha,0.0*Se);
//				mymaterial[imat]->SetTimeDependentFunctionExact(TimeDepFExact);		
//				
//				ofstream argm("mymaterial1.txt");
//				mymaterial[imat]->Print(argm);
//				TPZMaterial * mat(mymaterial[imat]);
//				mphysics->InsertMaterialObject(mat);
//				
//				///--- --- Inserting Boundary Conditions
//				// Squared Model for validation line source Validation give in ...........					
//				///--- --- Boundary Conditions
//				// Squared Model 1D compactation
//				
//				// Elastic problem -> 1 ,  Diffsion problem -> 2
//				
//				// Top Boundary	-3
//				// Setting free Boundary condition N1N2 means: Elastic -> Neumman, Diffusion -> Dirichlet
//				int N1D2 = 110;
//				TPZFMatrix<REAL> val13(3,2,0.), val23(3,1,0.);
//				REAL SigTopx	= 0.0;
//				REAL SigTopy	= 0.0;			
//				REAL pTop		= 0.0;
//				val23(0,0)=SigTopx;
//				val23(1,0)=SigTopy;			
//				val23(2,0)=pTop;
//				TPZMaterial * BCTOP = mymaterial[imat]->CreateBC(mat,bcTopSB, N1D2, val13, val23);
//				mphysics->InsertMaterialObject(BCTOP);			
//				
//				// Bottom Boundary -1					
//				/// Setting free displacement non permeable surface Boundary condition DYFX1N2 means: Elastic -> Dirichlet Y component, Diffusion -> Neumman
//				TPZFMatrix<REAL> val11(3,2,0.), val21(3,1,0.);
//				int DYFX1N2 = 100; 			
//				REAL uDBottx	=0.0;
//				REAL uDBotty	=0.0;
//				REAL fluxBott	=0.0;
//				val21(0,0)=uDBottx;
//				val21(1,0)=uDBotty;
//				val21(2,0)=fluxBott;
//				TPZMaterial * BCBOTT = mymaterial[imat]->CreateBC(mat,bcBottomSB, DYFX1N2, val11, val21);
//				mphysics->InsertMaterialObject(BCBOTT);
//				
//				// Letf Boundary -4		
//				/// Setting free displacement non permeable surface Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman
//				int DXFY1N2 = 010;
//				TPZFMatrix<REAL> val14(3,2,0.), val24(3,1,0.);
//				REAL uDLeftx	= 0.0;
//				REAL uDLefty	= 0.0;			
//				REAL fluxLeft	= 0.0;
//				val24(0,0)=uDLeftx;
//				val24(1,0)=uDLefty;			
//				val24(2,0)=fluxLeft;
//				TPZMaterial * BCLEFT = mymaterial[imat]->CreateBC(mat, bcLeftSB, DXFY1N2, val14, val24);
//				mphysics->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary -2			
//				/// Setting free displacement non permeable surface  Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman			
//				TPZFMatrix<REAL> val12(3,2,0.), val22(3,1,0.);
//				
//				REAL uDRightx	= 0.0;
//				REAL uDRighty	= 0.0;			
//				REAL fluxRight	= 0.0;
//				val22(0,0)=uDRightx;
//				val22(1,0)=uDRighty;			
//				val22(2,0)=fluxRight;
//				TPZMaterial * BCRIGHT = mymaterial[imat]->CreateBC(mat, bcRightSB, DXFY1N2, val12, val22);
//				mphysics->InsertMaterialObject(BCRIGHT);
//				
//			}
//			
//			
//			for(int imat = 1; imat < TotalMats; imat++)
//			{	
//				
//				// No Gravity field
//				forceVec[imat][1]= - gravity*rho;
//				mymaterial[imat] = new TPZPoroElastic2d (imat+1, dim);
//				mymaterial[imat]->SetParameters(lambVec[imat], GVec[imat],forceVec[imat][0], forceVec[imat][1]);
//				mymaterial[imat]->SetParameters(permVec[imat],viscVec[imat]);
//				mymaterial[imat]->SetfPlaneProblem(planestress);
//				mymaterial[imat]->SetBiotParameters(alpha,Se);
//				mymaterial[imat]->SetTimeDependentFunctionExact(TimeDepFExact);		
//				
//				ofstream argm("mymaterial1.txt");
//				mymaterial[imat]->Print(argm);
//				TPZMaterial * mat(mymaterial[imat]);
//				mphysics->InsertMaterialObject(mat);
//				
//				///--- --- Inserting Boundary Conditions
//				// Squared Model for validation line source Validation give in ...........					
//				///--- --- Boundary Conditions
//				// Squared Model 1D compactation
//				
//				// Elastic problem -> 1 ,  Diffsion problem -> 2
//				
//				// Top Boundary	-3
//				// Setting free Boundary condition N1N2 means: Elastic -> Neumman, Diffusion -> Dirichlet
//				int N1D2 = 111;
//				TPZFMatrix<REAL> val13(3,2,0.), val23(3,1,0.);
//				REAL SigTopx	= 0.0;
//				REAL SigTopy	= 0.0;			
//				REAL pTop		= 0.0;
//				val23(0,0)=SigTopx;
//				val23(1,0)=SigTopy;			
//				val23(2,0)=pTop;
//				TPZMaterial * BCTOP = mymaterial[imat]->CreateBC(mat,bcTopR, N1D2, val13, val23);
//				mphysics->InsertMaterialObject(BCTOP);			
//				
//				// Bottom Boundary -1					
//				/// Setting free displacement non permeable surface Boundary condition DYFX1N2 means: Elastic -> Dirichlet Y component, Diffusion -> Neumman
//				TPZFMatrix<REAL> val11(3,2,0.), val21(3,1,0.);
//				int DYFX1N2 = 111; 			
//				REAL uDBottx	=0.0;
//				REAL uDBotty	=0.0;
//				REAL fluxBott	=0.0;
//				val21(0,0)=uDBottx;
//				val21(1,0)=uDBotty;
//				val21(2,0)=fluxBott;
//				TPZMaterial * BCBOTT = mymaterial[imat]->CreateBC(mat,bcBottomR, DYFX1N2, val11, val21);
//				mphysics->InsertMaterialObject(BCBOTT);
//				
//				// Letf Boundary -4		
//				/// Setting free displacement non permeable surface Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman
//				int DXFY1N2 = 011;
//				TPZFMatrix<REAL> val14(3,2,0.), val24(3,1,0.);
//				REAL uDLeftx	= 0.0;
//				REAL uDLefty	= 0.0;			
//				REAL fluxLeft	= 0.0;
//				val24(0,0)=uDLeftx;
//				val24(1,0)=uDLefty;			
//				val24(2,0)=fluxLeft;
//				TPZMaterial * BCLEFT = mymaterial[imat]->CreateBC(mat, bcLeftR, DXFY1N2, val14, val24);
//				mphysics->InsertMaterialObject(BCLEFT);
//				
//				// Right Boundary -2			
//				/// Setting free displacement non permeable surface  Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman			
//				TPZFMatrix<REAL> val12(3,2,0.), val22(3,1,0.);
//				
//				REAL uDRightx	= 0.0;
//				REAL uDRighty	= 0.0;			
//				REAL fluxRight	= 0.0;
//				val22(0,0)=uDRightx;
//				val22(1,0)=uDRighty;			
//				val22(2,0)=fluxRight;
//				TPZMaterial * BCRIGHT = mymaterial[imat]->CreateBC(mat, bcRightR, DXFY1N2, val12, val22);
//				mphysics->InsertMaterialObject(BCRIGHT);
//				
//				//  Source Production/injection Well line (Q is the injected/produced volume per time [kg/s])
//				int Well = 011;
//				TPZFMatrix<REAL> valP1(3,2,0.), valP2(3,1,0.);
//				valP2(2,0) = qo/(rhof);
//				TPZMaterial * BCLINE = mymaterial[imat]->CreateBC(mat, WellLine, Well, valP1, valP2);
//				mphysics->InsertMaterialObject(BCLINE);			
//				
//				
//			}
//			
//			
//			
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
//			TPZManVector<std::string,10> scalnames(12), vecnames(3);
//			scalnames[0] = "SigmaX";
//			scalnames[1] = "SigmaY";
//			scalnames[2] = "TauXY";	
//			scalnames[3] = "DisplacementX";
//			scalnames[4] = "DisplacementY";
//			scalnames[5] = "SolidPressure";
//			scalnames[6] = "FluidPressure";
//			scalnames[7] = "EDisplacementX";
//			scalnames[8] = "EDisplacementY";
//			scalnames[9] = "EPressure";
//			scalnames[10] = "ESIGX";
//			scalnames[11] = "ESIGY";
//			vecnames[0]  = "Displacement";
//			vecnames[1]  = "FluxVector";
//			vecnames[2]  = "EFluxVector";	
//			
//			
//			
//			const int dim = 2;
//			int div =0;
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
//		void SolveSistTransient(TPZVec <REAL> &PrintEachtimeStep, std::string FileName, REAL deltaT,REAL maxTime, TPZFMatrix<REAL> &InitialSolution, TPZVec <TPZPoroElastic2d  * > &mymaterial ,
//								TPZAnalysis &an, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics)
//		{
//			
//			
//			// Calculationg Mass Matrix    
//			TPZAutoPointer <TPZMatrix<REAL> > PoroelasticMassMatrix = MassMatrix(mymaterial, mphysics);
//			
//			// Calculationg Stiffness Matrix and Load Vector
//			TPZFMatrix<REAL> PoroelasticStiffnessMatrix, PoroelasticStiffnessMatrixInverse;	
//			TPZFMatrix<REAL> PoroelasticLoadVector; 
//			StiffMatrixLoadVec(mymaterial, mphysics, an, PoroelasticStiffnessMatrix, PoroelasticLoadVector);
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
//				if (PrintEachtimeStep[control] == TimeValue) 
//				{
//					std::stringstream outputfiletemp;
//					outputfiletemp << FileName << ".vtk";
//					std::string plotfile = outputfiletemp.str();
//					PostProcessPoroeasticity(meshvec,mphysics,an,plotfile);
//					control++;					
//				}
//				
//				cent++;
//				TimeValue = cent*deltaT;
//			}
//			
//		}
//
//
//		//	 TPZVec <REAL> ModelDimensions, REAL Freaticlevel, REAL FluidDensity 
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

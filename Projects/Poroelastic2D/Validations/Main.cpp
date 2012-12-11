#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"

#include <pzgengrid.h>
#include "../MeshGeneration.h"
#include "../AnalyticalFunctions.h"
#include "../ElasticMatInterface2D.h"
//#include "../PoroElasticMatInterface2D.h"
#include "TPZInterfaceEl.h"
#include "pzdiscgal.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include <tpzarc3d.h>

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzl2projection.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

#include "pzelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSpStructMatrix.h"
#include "../pzporoelastic2d.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "TPZReadGIDGrid.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "tinystr.h"
#include "tinyxml.h"


#include <cmath>
#include <set>


// Task to do
// Document all the code
// Include all BC combinations
// Include Finite elasticity
// Conservation mass mini benchmark (tiago's suggestion)
// How I can divide the code in cases for the realistic analysis?
// Write Dimensionless Poroelastic Formulation



// Using Log4cXX as logging tool
//
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.poroelastic2d"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic.data"));
#endif
//
// End Using Log4cXX as logging tool

// Defintions of names space
using namespace std;
using namespace pzgeom;	

// Defitions of Material index
// Internal Materials
// In this case just one material is used
const int TotalMats = 2; 
const int matId = 1;	

// Boundary Conditions Materials
// 2D Rectangular Domain for sideburden
const int bcBottomSB = 7;
const int bcRightSB = 10;
const int bcTopSB = 9;
const int bcLeftSB = 8;

// 2D Rectangular Domain for sideburden
const int bcBottomR = 3;
const int bcRightR = 6;
const int bcTopR = 5;
const int bcLeftR = 4;

// line source identity
const int WellLine = 20;

// 2D Rectangular Domain
const int bcBottom = -1;
const int bcRight = -2;
const int bcTop = -3;
const int bcLeft = -4;
const int yfixedPoints = -6;
const int xfixedPoints = -7;
const int pointsource = -5;

// Rock Blocks bondaries
// Right Block 
int RBright = 5;
int RBleft = 3;
int RBBot = 2;
int RBTop = 4;

//int RBright = 10;
//int RBleft = 8;
//int RBBot = 7;
//int RBTop = 9;

// Left Block 
int LBright = 6;
int LBleft = 8;
int LBBot = 5;
int LBTop = 9;

// Graven Block 
int GBright = 6;
int GBleft = 4;
int GBBot = 3;
int GBTop = 5;


// Boundary Conditions Definitions
// This Definitions are related with the Boudary conditions convention for multphysics simulation (in Poroelastic U/P formulation problems 3 boundary conditions are included) 
// not defined
const int dirichlet = 0;
const int neumann = 1;


// Defintions of Implemented Methods

//	This Create Computational Meshes
TPZCompMesh *ComputationalDiffusionMesh(TiXmlHandle ControlDoc,TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh,int pOrder,TPZVec < TPZFMatrix<REAL> > Events, int TStep);
TPZCompMesh *ComputationalElasticityMesh(TiXmlHandle ControlDoc,TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh,int pOrder, bool Initialstress,TPZVec < TPZFMatrix<REAL> > Events, int TStep);
TPZCompMesh *ComputationalPoroelasticityMesh(TiXmlHandle ControlDoc, TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,TPZVec <TPZPoroElastic2d  * > &mymaterial,TPZVec < TPZFMatrix<REAL> > Events, int TStep);
//TPZCompMesh *ComputationalPoroelasticityInitialMesh(TiXmlHandle ControlDoc, TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial);
TPZCompMesh *SetInitialConditionsbyL2Projection(TiXmlHandle ControlDoc, TPZReadGIDGrid GeometryInfo, TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &InitialCondition, bool &CustomFunction);

//	This Set intial Conditions
void SetInitialConditions(TPZAnalysis &AnalysisToInitialize, TPZCompMesh * ComputationalMesh, TPZFMatrix<REAL> &ComputationalMeshInitialSolution, void ( *PointerToFunction(REAL &,TPZVec <REAL> &)));
void *InitalPressureCalculations(REAL &PressureValue, TPZVec <REAL> &PointCoordinates);
void *InitialDisplacementCalculations();

//	This Solve Different analysis
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);
TPZAutoPointer <TPZMatrix<REAL> > MassMatrix(TPZReadGIDGrid GeometryInfo,TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics);
void SolveSistTransient(TPZVec < TPZFMatrix<REAL> > Events,TiXmlHandle ControlDoc,int Nthreads, TPZVec <REAL> &PrintEachtimeStep, std::string FileName, REAL deltaT,REAL maxTime, TPZFMatrix<REAL> &InitialSolution, 
						TPZVec <TPZPoroElastic2d  * > &mymaterial , TPZAnalysis &an, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZReadGIDGrid GeometryInfo);	
void StiffMatrixLoadVec(TPZReadGIDGrid GeometryInfo,int Nthreads, TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec);

//	This are tools for spatial and polinomyal refinement and Postprocess of solutions 
void PostProcessElasticity(TPZAnalysis &an, std::string plotfile);
void PostProcessDiffusion(TPZAnalysis &an, std::string plotfile);
void PostProcessPoroeasticity(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void RefinamentoUniforme(TPZGeoMesh  *gMesh, int nh);
void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh, int MatId, int indexEl);
void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);
void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file);

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	InitializePZLOG("../mylog4cxx.cfg");
#endif		
	
	
	int Nthreads = 1;
	int ref = 0;
	
	
	//	Using main arguments 
	{
		int argisize = argc;
		
		using namespace std;
		// If the user didn't provide a filename command line argument,
		// print an error and exit.
		if (argc != 3)
		{
			cout << "Size: " << argc << " Number of Arguments " << endl;
			cout << "Usage: " << argv[0] << " Nthreads  ref" << endl;
			cout << "Used " << " Default" << endl;
			//			DebugStop();
		}
		
		if (argc == 3)
		{
			cout << "Used int: " << argv[1] << " Nthreads" << endl;
			cout << "Used int: " << argv[2] << " ref" << endl;
			Nthreads		= atoi(argv[1]);
			ref				= atoi(argv[2]);			
		}
		
	}		
	
	
	
	std::string ControlFile("MyDataStructure.xml");	
	
	// control reading of xml file
	// load the named file and dump its structure to STDOUT
	TiXmlDocument doc("MyDataStructure.xml");
	bool loadOkay = false;
	loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		cout << "Th file is ok! \n" << ControlFile << endl;
	}
	else
	{
		cout << "Failed to load file \n" << ControlFile << endl;
	}
	
	
	// Accesing data
	
	TiXmlHandle docHandle( &doc );
	TiXmlElement* GridDumpName = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "GridDumpFile" ).ToElement();
	const char * GridFileName = GridDumpName->Attribute("GridFile");
	cout << "Find = \n" << GridFileName << endl;	
	TiXmlElement* LmaxDimension = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "MaxGeometricDim" ).ToElement();	
	const char * LDimension = LmaxDimension->Attribute("L");
	cout << "Find = \n" << LDimension << endl;	
	
	
	
	
	// GEOMETRICAL MESH CREATION	
	
	//	polynomial base degree approach
	int p=2;
	
	
	TPZReadGIDGrid GeometryInfo;
	GeometryInfo.SetfDimensionlessL(atof(LDimension));
	TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);	
	
	
	{
		//	Print Geometrical Base Mesh
		ofstream argument("BaseGeoMesh.txt");
		gmesh->Print(argument);
		ofstream Dummyfile("BaseGeoMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
	}
	
	
	int NumOfThreads=0;
	int Href = 0;
	//	
	//	if (Nthreads != 0) 
	//	{
	NumOfThreads = Nthreads;
	Href = ref;
	//	}	
	
	
	RefinamentoUniforme(gmesh, Href);
	
	//	// Directional Point source Refinement 
	//	int ndirectdivp = 17;
	//	set<int> SETmatPointRefDir2;
	//	SETmatPointRefDir2.insert(WellLine);
	//	//	SETmatPointRefDir2.insert(bcBottom);
	//	//	SETmatPointRefDir2.insert(bcTop);
	//	//	SETmatPointRefDir2.insert(bcRight);
	//	//	SETmatPointRefDir2.insert(bcLeft);	
	//	//	SETmatPointRefDir2.insert(xfixedPoints);
	//	//	SETmatPointRefDir2.insert(yfixedPoints);	
	//	for(int j = 0; j < ndirectdivp; j++)
	//	{
	//		int nel = gmesh->NElements();
	//		for (int iref = 0; iref < nel; iref++)
	//		{
	//			TPZVec<TPZGeoEl*> filhos;
	//			TPZGeoEl * gelP2 = gmesh->ElementVec()[iref];
	//			if(!gelP2 || gelP2->HasSubElement()) continue;
	//			TPZRefPatternTools::RefineDirectional(gelP2, SETmatPointRefDir2);
	//		}		
	//	}	
	
	{
		//	Print Geometrical refined Base Mesh
		ofstream argument("RefineGeoMesh.txt");
		gmesh->Print(argument);
		ofstream Dummyfile("RefineGeoMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);	
	}	
	
	
	// TiXmlElement dummy object
	TiXmlElement* Container;
	const char * CharContainer;	
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Events" ).FirstChild( "WellAlterations" ).ToElement();			
	CharContainer = Container->Attribute("EventsNumber");	
	int NEvents = atoi(CharContainer);	
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Events" ).FirstChild( "WellAlterations" ).ToElement();			
	CharContainer = Container->Attribute("WellsPerEvent");	
	int WellsPerEvent = atoi(CharContainer);
	
	// Reading well recurrent Data
	TPZFMatrix<REAL> ALter(WellsPerEvent,8);
	TPZVec < TPZFMatrix<REAL> > Events(NEvents);
	// Reading alterations
	int iEvent = 0;
	if (NEvents!=0) 
	{
		Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "WellLines" ).FirstChild( "Event" ).ToElement();
		for( Container; Container; Container=Container->NextSiblingElement())
		{
			
			int TemptimeStep = atoi(Container->Attribute( "TimeStep"));
			int iAlter = 0;
			TiXmlElement* LineSource_el = Container->FirstChild("LineSource")->ToElement();
			for( LineSource_el; LineSource_el; LineSource_el=LineSource_el->NextSiblingElement())
			{	
				ALter(iAlter,0)=TemptimeStep;
				CharContainer = LineSource_el->Attribute( "ID" );
				ALter(iAlter,1)=atoi(CharContainer);
				CharContainer = LineSource_el->Attribute( "WitchDomainMaterial" );			
				ALter(iAlter,2)=atoi(CharContainer);
				CharContainer = LineSource_el->Attribute( "ConstrictionType" );			
				ALter(iAlter,3)=atoi(CharContainer);
				CharContainer = LineSource_el->Attribute( "Val1" );			
				ALter(iAlter,4)=atof(CharContainer);
				CharContainer = LineSource_el->Attribute( "Val2First" );			
				ALter(iAlter,5)=atof(CharContainer);
				CharContainer = LineSource_el->Attribute( "Val2Second" );			
				ALter(iAlter,6)=atof(CharContainer);
				CharContainer = LineSource_el->Attribute( "Val2Third" );			
				ALter(iAlter,7)=atof(CharContainer);
				iAlter++;
			}
			//		cout << "Data readed" << ALter(iAlter,0) << ALter(iAlter,1) << ALter(iAlter,2) << ALter(iAlter,3) << ALter(iAlter,4) << ALter(iAlter,5) << ALter(iAlter,6) << ALter(iAlter,7) <<endl;
			Events[iEvent] = ALter; 
			iEvent++;
		}	
	}
	
	
	//	Intial Calculations
	int TimeStep = 0;	
	
	//	First computational mesh
	
	// initialization problem ...
	bool Initialstress = true;
	TPZCompMesh * ComputationalMeshElasticity = ComputationalElasticityMesh(docHandle,GeometryInfo, gmesh, p+1, Initialstress, Events, TimeStep);
	//	Print First computational mesh
	ofstream ArgumentElasticity("cmeshElasticity.txt");
	ComputationalMeshElasticity->Print(ArgumentElasticity);
	
	//	Second computational mesh
	TPZCompMesh * ComputationalMeshDiffusion = ComputationalDiffusionMesh(docHandle,GeometryInfo, gmesh, p, Events, TimeStep);
	//	Print Second computational mesh
	ofstream ArgumentDiffusion("cmeshDiffusion.txt");
	ComputationalMeshDiffusion->Print(ArgumentDiffusion);	
	
	//	Cleaning reference of the geometric mesh for cmesh1
	gmesh->ResetReference();
	ComputationalMeshElasticity->LoadReferences();
	
	//	Using Uniform Refinement for first computational problem
	RefinUniformElemComp(ComputationalMeshElasticity,0);
	ComputationalMeshElasticity->AdjustBoundaryElements();
	ComputationalMeshElasticity->CleanUpUnconnectedNodes();
	
	ofstream ArgumentElasticityRef("cmeshElasticityRef.txt");
	ComputationalMeshElasticity->Print(ArgumentElasticityRef);
	ofstream ArgumentElasticityGeoRef("GeoMeshElasticity.txt");
	gmesh->Print(ArgumentElasticityGeoRef);
	
	// Vtk visualization Elasticity Geometric Mesh
	ofstream DummyfileElasticity("GeoMeshElasticity.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, DummyfileElasticity,true);
	
	// Clear reference of the geometric mesh to cmesh2
	gmesh->ResetReference();
	ComputationalMeshDiffusion->LoadReferences();
	
	//	Using Uniform Refinement for second computational problem	
	RefinUniformElemComp(ComputationalMeshDiffusion,0);
	ComputationalMeshDiffusion->AdjustBoundaryElements();
	ComputationalMeshDiffusion->CleanUpUnconnectedNodes();
	//	ComputationalMeshDiffusion->ExpandSolution();
	
	ofstream ArgumentDiffusionRef("cmeshDiffusionRef.txt");
	ComputationalMeshDiffusion->Print(ArgumentDiffusionRef);
	ofstream ArgumentDiffusionGeoRef("GeoMeshDiffusion.txt");
	gmesh->Print(ArgumentDiffusionGeoRef);
	
	// Vtk visualization Diffusion Geometric Mesh	
	ofstream DummyfileDiffusion("GeoMeshDiffusion.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, DummyfileDiffusion,true);
	
	
	//--------------------------------------------------------------------------------------------------------------------------------//	
	//--------------------------------------------------------------------------------------------------------------------------------//	
	
	// Setting intial conditions
	
	//	//	Set initial conditions for pressure
	//	TPZFMatrix<REAL> DiffusionSolution;
	//	TPZAnalysis ElasticAnalysis(ComputationalMeshElasticity);
	//	TPZAnalysis DiffusionAnalysis(ComputationalMeshDiffusion);	
	//	//	Solving Second Computational problem
	//	SolveSist(DiffusionAnalysis, ComputationalMeshDiffusion);
	//	std::string plotfilediffusion("SolutionDiffusion.vtk");		
	//	PostProcessDiffusion(DiffusionAnalysis, plotfilediffusion);	
	//	
	//	SolveSist(ElasticAnalysis,ComputationalMeshElasticity);		
	//	std::string OutTubeElasticity("SolutionElasticity.vtk");		
	//	PostProcessElasticity(ElasticAnalysis, OutTubeElasticity);		
	//	
	//	SolveSist(DiffusionAnalysis,ComputationalMeshDiffusion);
	//	std::string OutTubeDiffusion("SolutionDiffusion.vtk");	
	//	PostProcessDiffusion(DiffusionAnalysis, OutTubeDiffusion);		
	//	
	//	
	//	Without L2 projection
	//	SetInitialConditions(DiffusionAnalysis, ComputationalMeshDiffusion, DiffusionSolution, *InitalPressureCalculations);
	//	DiffusionAnalysis.LoadSolution(DiffusionSolution);
	
	//	Set initial conditions for pressure by using L2 projection
	//	int nrs = DiffusionAnalysis.Solution().Rows();
	//	bool CustomFunction = false;
	//	//	Using  constant initial pressure
	//    TPZVec<REAL> InitialCondition(nrs,1000.0);
	//    TPZCompMesh  * ComputationalL2 = SetInitialConditionsbyL2Projection(gmesh, p, InitialCondition, CustomFunction);
	//    TPZAnalysis L2Analysis(ComputationalL2);
	//    SolveSist(L2Analysis, ComputationalL2);
	//    DiffusionAnalysis.LoadSolution(L2Analysis.Solution());	
	
	
	//#ifdef LOG4CXX
	//	if(logdata->isDebugEnabled())
	//	{		
	//		std::stringstream sout;
	//		DiffusionAnalysis.Solution().Print("Intial conditions = ", sout,EMathematicaInput);
	//		LOGPZ_DEBUG(logdata,sout.str())
	//	}
	//#endif	
	
	TPZVec<TPZCompMesh *> ComputationalMeshVectorInitial(2);
	ComputationalMeshVectorInitial[0] = ComputationalMeshElasticity;
	ComputationalMeshVectorInitial[1] = ComputationalMeshDiffusion;
	
	TPZVec<TPZCompMesh *> ComputationalMeshVector(2);
	ComputationalMeshVector[0] = ComputationalMeshElasticity;
	ComputationalMeshVector[1] = ComputationalMeshDiffusion;	
	
	TPZVec <TPZPoroElastic2d *>  materialist(GeometryInfo.MatNumber,0);	
	
	TPZCompMesh * ComputationalMeshPoroelasticityInitial = ComputationalPoroelasticityMesh(docHandle,GeometryInfo, gmesh,ComputationalMeshVectorInitial,materialist,Events,TimeStep);
	TPZCompMesh * ComputationalMeshPoroelasticityReCurrent = ComputationalPoroelasticityMesh(docHandle,GeometryInfo, gmesh,ComputationalMeshVector,materialist,Events,TimeStep);	
	
	//--------------------------------------------------------------------------------------------------------------------------------//	
	//--------------------------------------------------------------------------------------------------------------------------------//	
	
	// time control
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "FileName" ).ToElement();		
	CharContainer = Container->Attribute("ProblemName");
	std::string	OutFileName(CharContainer);
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeControls" ).ToElement();		
	CharContainer = Container->Attribute("Initial");	
	REAL InitialTime = atof(CharContainer);
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeControls" ).ToElement();			
	CharContainer = Container->Attribute("Final");	
	REAL MaxTime = atof(CharContainer);	
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeControls" ).ToElement();		
	CharContainer = Container->Attribute("DeltaTime");		
	REAL Delta = atof(CharContainer);	
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeControls" ).ToElement();			
	CharContainer = Container->Attribute("Theta");	
	REAL Theta = atof(CharContainer);
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeSchedule" ).ToElement();			
	CharContainer = Container->Attribute("Nvalues");	
	int NTimeValues = atoi(CharContainer);
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeSchedule" ).ToElement();			
	CharContainer = Container->Attribute("ScheduleTime");	
	std::ifstream Timefile (CharContainer);
	
	// Reading time schedule
	REAL TimeValue	=	0.0;
	TPZVec <REAL> PrintStep(NTimeValues,0);	
	for (int itime = 0; itime < NTimeValues; itime++) 
	{
		char buf[1024];
		Timefile.getline(buf, 1024);
		Timefile >> TimeValue;
		PrintStep[itime] =  TimeValue;
	}	
	
	// Improve this
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{		
		materialist[imat]->SetTimeStep(Delta,Theta);
		materialist[imat]->SetTimeValue(InitialTime);			
	}
	
	TPZAnalysis PoroelasticAnalysisInitial(ComputationalMeshPoroelasticityInitial);	
	TPZAnalysis PoroelasticAnalysis(ComputationalMeshPoroelasticityReCurrent);
	
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "FileName" ).ToElement();		
	CharContainer = Container->Attribute("InitialProblemName");
	//Setting initial coditions 
	TPZFMatrix<REAL> InitialSolution = PoroelasticAnalysisInitial.Solution();
	std::string output(CharContainer);
	std::stringstream outputfiletemp;
	outputfiletemp << output << ".vtk";
	std::string plotfile = outputfiletemp.str();
	
	//	Calculation of Initial conditions
	cout << "Calculation of initial conditions !" << endl;
	SolveSistTransient(Events , docHandle, NumOfThreads,PrintStep,output,Delta,MaxTime,InitialSolution, materialist, PoroelasticAnalysisInitial, ComputationalMeshVectorInitial,ComputationalMeshPoroelasticityInitial,GeometryInfo);
	//	//Setting initial coditions 
	//	PoroelasticAnalysis.LoadSolution(PoroelasticAnalysisInitial.Solution());	
	//	TPZFMatrix<REAL> InitialSolution = PoroelasticAnalysis.Solution();	
	//	
	//	//	//	Print Initial conditions	
	//	//	PostProcessPoroeasticity(ComputationalMeshVectorInitial,ComputationalMeshPoroelasticityInitial,PoroelasticAnalysisInitial,plotfile);	
	//	
	//	//	TPZSpStructMatrix matsp(mphysics);
	//	//	//TPZSkylineStructMatrix matsp(mphysics);	
	//	//	std::set< int > materialid;
	//	//	materialid.clear();	
	//	//	//	int matid = 1;
	//	//	//	materialid.insert(imat);
	//	//	//	matsp.SetMaterialIds (materialid);
	//	//	// Improve this
	//	//	for(int imat = 0; imat < TotalMats; imat++)
	//	//	{	
	//	//		materialid.insert(imat+1);	
	//	//	}
	//	//	// Inserting Material ids for matrix K2 creation
	//	//	matsp.SetMaterialIds(materialid);		
	//	//	
	//	//	TPZAutoPointer<TPZGuiInterface> guiInterface;
	//	//	TPZFMatrix<REAL> Un;
	//	//	TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
	//	
	//	//#ifdef LOG4CXX
	//	//	if(logdata->isDebugEnabled())
	//	//	{
	//	//		//		std::stringstream sout;
	//	//		//		matK2->Print("K2 = ", sout,EMathematicaInput);
	//	//		//		LOGPZ_DEBUG(logdata,sout.str())
	//	//	}
	//	//#endif
	//	
	//	//	//Criando matriz K1
	//	//	// Improve this
	//	//	for(int imat = 0; imat < TotalMats; imat++)
	//	//	{	
	//	//		materialist[imat]->SetCurrentState();
	//	//	}	
	//	//	TPZSkylineStructMatrix matsk(mphysics);
	//	//	TPZFMatrix<REAL> matK1;	
	//	//	TPZFMatrix<REAL> fvec; //vetor de carga
	//	//	an.SetStructuralMatrix(matsk);//	Set the structtural sky line matrix to our analysis
	//	//	
	//	//	
	//	
	//	//	an.Run(); //	Excecute an analysis to obtain the Rhs vector (In this case we start with zero initial values)
	//	//	
	//	//	matK1 = an.StructMatrix(); //Storage the global matrix and load vector
	//	//	fvec = an.Rhs();
	//	
	//	// This code identify singular blocks
	//	TPZStepSolver<REAL> step; //Create Solver object
	//	step.SetDirect(ELDLt); //	Symmetric case
	//	PoroelasticAnalysis.SetSolver(step); //	Set solver	
	//	TPZStepSolver<REAL> & temp = dynamic_cast<TPZStepSolver<REAL> &> (PoroelasticAnalysis.Solver());
	//	std::list <int> & zeropivot = temp.Singular(); 
	//	if (zeropivot.size()) 
	//	{
	//		int eq = * zeropivot.begin();
	//		PoroelasticAnalysis.Rhs().Zero();
	//		PoroelasticAnalysis.Rhs()(eq,0) = -10000.0;
	//		PoroelasticAnalysis.Solve();
	//		TPZFMatrix<REAL> TempSolution = PoroelasticAnalysis.Solution();
	//		
	//#ifdef LOG4CXX
	//		// Print the temporal solution
	//		if(logdata->isDebugEnabled())
	//		{
	//			std::stringstream sout;
	//			TempSolution.Print("SingularNodes = ", sout,EMathematicaInput);
	//			LOGPZ_DEBUG(logdata,sout.str())
	//		}
	//#endif	
	//		std::string output;
	//		output = "SingularNodes";
	//		std::stringstream outputfiletemp;
	//		outputfiletemp << output << ".vtk";
	//		std::string plotfile = outputfiletemp.str();
	//		PostProcessPoroeasticity(ComputationalMeshVector,ComputationalMeshPoroelasticityReCurrent,PoroelasticAnalysis,plotfile);
	//		
	//		// Probelem bad pose
	//		DebugStop();
	//	}
	//	
	//	//#ifdef LOG4CXX
	//	//	if(logdata->isDebugEnabled())
	//	//	{
	//	//		// base system to invert
	//	//		// just one for checking purpose
	//	//		//		std::stringstream sout;
	//	//		//		an.Solver().Matrix()->Print("K1 = ", sout,EMathematicaInput);
	//	//		//		fvec.Print("fvec = ", sout,EMathematicaInput);		
	//	//		//		//Print the temporal solution
	//	//		//		Initialsolution.Print("Intial conditions = ", sout,EMathematicaInput);
	//	//		//		TPZFMatrix<REAL> Temp;
	//	//		//		TPZFMatrix<REAL> Temp2;
	//	//		//		matK2->Multiply(Initialsolution,Temp);
	//	//		//		Temp.Print("Temp K2 = ", sout,EMathematicaInput);	
	//	//		//		LOGPZ_DEBUG(logdata,sout.str())
	//	//	}
	//	//#endif
	//	
	//	///start transient problem
	//	
	//	SolveSistTransient(NumOfThreads,PrintEachtimeStep,OutFileName,Delta,MaxTime, InitialSolution, materialist, PoroelasticAnalysis, ComputationalMeshVector,ComputationalMeshPoroelasticityReCurrent);	
	
	
	
	
	//	Fault Stability Factor Calculations
	
	//	int totalel = ComputationalMeshPoroelasticityInitial->NElements();
	//	int nsegments = 4;
	//	
	//	std::map< REAL,TPZVec<REAL> > sN1T1N2T2;
	//	TPZVec<REAL> vecDelta(3,0.);
	//	REAL s = 0.;
	//	
	//	for (int iel = 0; iel < totalel; iel++) 
	//	{
	//		// procurando elemento computacional ne MaterialId = 7 (elemento 1D da falha)
	//		TPZCompEl * CompEl = ComputationalMeshPoroelasticityInitial->ElementVec()[iel];
	//		
	//		#ifdef DEGUG
	//		if(CompEl->Dimension() != 1)
	//		{
	//			DebugStop();//nao eh elemento 1D!!!!!!!
	//		}
	//		#endif
	//		
	//		int CompElMatID = CompEl->Material()->Id();
	//		
	//		if (CompElMatID == 8) 
	//		{	
	//			// elemento geometrico assiciado ao CompEl
	//			TPZGeoEl * GeoEl = CompEl->Reference();
	//		
	//			TPZMultiphysicsInterfaceElement * InterFace =  dynamic_cast<TPZMultiphysicsInterfaceElement *> (CompEl);
	////			
	////			TPZManVector<TPZMaterialData,6> datavecleft,datavecright;
	////			TPZMaterialData data;
	////			InitMaterialData(datavecleft, leftel);
	////			InitMaterialData(datavecright, rightel);
	////			TPZManVector<TPZTransform> leftcomptr, rightcomptr;
	////			leftel->AffineTransform(leftcomptr);
	////			rightel->AffineTransform(rightcomptr);
	////			InitMaterialData(data);			
	//
	//			
	//			//Capturando as coordenadas dos nohs inicial e final
	//			TPZVec<REAL> node0Coord(3,0.), node1Coord(3,0.);
	//			GeoEl->NodePtr(0)->GetCoordinates(node0Coord);
	//			GeoEl->NodePtr(1)->GetCoordinates(node1Coord);
	//			
	////			TPZVec< TPZFMatrix<REAL> Sigma(?,?) > SigmaPt(nsegments+1);
	////					->loop para inicializar cada posicao do vetor SigmaPt
	//
	//			// percorrendo os pontos intermediarios entre os nohs inicial e final
	//			for(int seg = 0; seg <= nsegments; seg++)
	//			{
	//				TPZVec<REAL> ptX(3,0.);
	//				for(int c = 0; c < 3; c++)
	//				{
	//					//ponto EM X intermediario correspondente aa sec
	//					vecDelta[c] = seg*(node1Coord[c]-node0Coord[c])/nsegments;
	//					ptX[c] = node0Coord[c] + vecDelta[c];
	//					
	//					//capturando os vizinhos ao elemento 1D (GeoEl)
	//					TPZGeoElSide GeoElSide1D(GeoEl, 2);
	//					TPZGeoElSide NeighGeoElSide1D(GeoElSide1D.Neighbour());
	//					
	//					int nElsFounded = 0;
	//					while(GeoElSide1D != NeighGeoElSide1D)
	//					{
	//						TPZGeoEl * neighGeoEl = NeighGeoElSide1D.Element();
	//						if(neighGeoEl && neighGeoEl->Dimension() == 2 && 
	//						   (neighGeoEl->MaterialId() == 1 || neighGeoEl->MaterialId() == 2))
	//						{
	//							nElsFounded++;
	//							
	//							//Transferindo o ponto EM X para o espaco parametrico do elemento neighGeoEl
	//							TPZVec<REAL> qsi(neighGeoEl->Dimension(),0.);
	//							bool succeeded = neighGeoEl->ComputeXInverse(ptX,qsi);
	//							
	//							#ifdef DEBUG
	//							if(succeeded == false)
	//							{
	//								DebugStop();//nao conseguiu calcular o qsi associado ao ptX!!!!
	//							}
	//							#endif
	//							
	////							neighGeoEl->Reference()->Solution(qsi,var,Mydata)
	//							//		calcular normal
	//							
	//	//						neighGeoEl->Reference()->ComputeSolution(qsi,var,SigmaPt[seg]);
	////							TPZFMatrix<REAL> SigmaN(2,0.), SigmaTau(2,0.);
	//							
	//							int side = NeighGeoElSide1D.Side();
	//							
	//						}
	//					}
	//					
	//					#ifdef DEBUG
	//					if(nElsFounded != 2)
	//					{
	//						DebugStop();//deveria encontrar exatamente 2 elementos que satisfazem a exigencia de ser dim=2 e matId=(1||2)
	//					}
	//					#endif
	//				}
	//				
	//			}
	
	//>>>>>>
	//			TPZManVector<TPZMaterialData,6> datavecleft,datavecright;
	//			TPZMaterialData data;
	//			TPZMultiphysicsInterfaceElement * InterFace =  dynamic_cast<TPZMultiphysicsInterfaceElement *> (CompEl);
	//			
	//			
	//			TPZCompElSide LeftSideElement;	
	//			TPZCompElSide RightSideElement;			
	//			
	//			InterFace->GetLeftRightElement(LeftSideElement,RightSideElement);
	//			
	//			TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (LeftSideElement.Element());
	//			TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(RightSideElement.Element());	
	//					
	//			InterFace->InitMaterialData(datavecleft, leftel);
	//			InterFace->InitMaterialData(datavecright, rightel);
	//			
	//			
	//
	////			TPZVec<TPZMaterialData> &datavec
	//			int var = 6;
	//			TPZVec<REAL> SolPressure;
	////			
	//			TPZPoroElastic2d *Poroelasticleftel = dynamic_cast<TPZPoroElastic2d *> (leftel);
	//			TPZPoroElastic2d *Poroelasticrightel = dynamic_cast<TPZPoroElastic2d *>(rightel);				
	//			
	//			
	//			
	//			Poroelasticleftel->FillDataRequirements(datavecleft);
	//			Poroelasticrightel->FillDataRequirements(datavecright);		
	//			
	//			int Presusress = SolPressure[0];
	//			
	//			std::cout <<  "Pressures " << SolPressure[0] << std::endl;		
	
	
	
	//	for material on right 
	
	//			TPZInterpolationSpace *mspleft  = dynamic_cast <TPZInterpolationSpace *>(LeftSideElement.Element());
	//			if(!mspleft) continue;
	//			mspleft->InitMaterialData(datavec[iref]);
	//			TPZMaterialData::MShapeFunctionType shapetype = datavec[iref].fShapeType;
	//			if(shapetype==datavec[iref].EVecShape) continue;
	//			
	//			trvec[iref].Apply(qsi, myqsi);
	//			datavec[iref].p = msp->MaxOrder();
	//			msp->ComputeShape(qsi,datavec[iref]);
	//			msp->ComputeSolution(myqsi, datavec[iref]);
	//			
	//			datavec[iref].x.Resize(3);
	//			msp->Reference()->X(myqsi, datavec[iref].x);			
	
	//			// for material on left
	//			
	//			TPZInterpolationSpace *mspright  = dynamic_cast <TPZInterpolationSpace *>(RightSideElement);
	//			if(!msp) continue;
	//			msp->InitMaterialData(datavec[iref]);
	//			TPZMaterialData::MShapeFunctionType shapetype = datavec[iref].fShapeType;
	//			if(shapetype==datavec[iref].EVecShape) continue;
	//			
	//			trvec[iref].Apply(qsi, myqsi);
	//			datavec[iref].p = msp->MaxOrder();
	//			msp->ComputeShape(qsi,datavec[iref]);
	//			msp->ComputeSolution(myqsi, datavec[iref]);
	//			
	//			datavec[iref].x.Resize(3);
	//			msp->Reference()->X(myqsi, datavec[iref].x);
	//			
	//			
	//			void TPZMultiphysicsCompEl<TGeometry>::Solution(TPZVec<REAL> &qsi, int var,TPZVec<REAL> &sol) {
	//				
	//				if(var >= 100) {
	//					TPZCompEl::Solution(qsi,var,sol);
	//					return;
	//				}
	//				
	//				TPZMaterial * material = this->Material();
	//				if(!material){
	//					sol.Resize(0);
	//					return;
	//				}
	//				
	//				TPZManVector<TPZTransform> trvec;
	//				AffineTransform(trvec);
	//				
	//				TPZVec<REAL> myqsi;
	//				myqsi.resize(qsi.size());
	//				
	//				int nref = fElementVec.size();
	//				TPZVec<TPZMaterialData> datavec;
	//				datavec.resize(nref);
	//				
	//				for (int iref = 0; iref<nref; iref++)
	//				{		
	//
	//				}
	//				
	//				material->Solution(datavec, var, sol);
	//			}			
	//			
	//			TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
	//			TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
	//			TPZGeoEl *leftgel = leftel->Reference();
	//			TPZGeoEl *rightgel = rightel->Reference();			
	//			CompEl->InitMaterialData(datavecright, rightel);
	
	
	//		}
	//		
	//	}
	
	cout << "Check:: Calculation finished successfully" << endl;	
	return EXIT_SUCCESS;
}


TPZCompMesh * ComputationalElasticityMesh(TiXmlHandle ControlDoc, TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh,int pOrder, bool Initialstress,TPZVec < TPZFMatrix<REAL> > Events, int TStep)
{
	// TiXmlElement dummy object
	TiXmlElement* Container;
	const char * CharContainer;
	
	// Plane strain assumption
	int planestress = 0;
	
	// Getting mesh dimension
	int dim = GeometryInfo.fProblemDimension;
	
	// Aproximation Space of order -> pOrder
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravitationalConstant" ).ToElement();		
	CharContainer = Container->Attribute("Gravity");
	REAL Gravitationalconstant = atof(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravityDirection" ).ToElement();		
	CharContainer = Container->Attribute("x-direction");	
	REAL Gravity_Xdirection = atof(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravityDirection" ).ToElement();		
	CharContainer = Container->Attribute("y-direction");	
	REAL Gravity_Ydirection = atof(CharContainer);	
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "FaultInterface" ).FirstChild( "Contact" ).ToElement();		
	CharContainer = Container->Attribute("Interfaces");		
	int InterfacesCutoff = atof(CharContainer);	
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "FaultInterface" ).FirstChild( "Contact" ).ToElement();		
	CharContainer = Container->Attribute("ContacFlag");	
	int InterfaceBehavior = atoi(CharContainer);	
	
	// Data for elastic problem
	REAL Eyoung			=	0.0;
	REAL PoissonRatio	=	0.0;
	REAL Lambda			=	0.0;
	REAL G				=	0.0;	
	REAL RockDensity	=	0.0;	
	REAL BodyForceX		=	0.0;
	REAL BodyForceY		=	0.0;	
	TPZElasticityMaterial * MaterialElastic;
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{	
		
		// Using Gravity field
		BodyForceX		= Gravitationalconstant * Gravity_Xdirection * GeometryInfo.fMaterialDataVec[imat].Properties[0];		
		BodyForceY		= Gravitationalconstant * Gravity_Ydirection * GeometryInfo.fMaterialDataVec[imat].Properties[0];
		Lambda			= GeometryInfo.fMaterialDataVec[imat].Properties[6];
		G				= GeometryInfo.fMaterialDataVec[imat].Properties[7];		
		RockDensity		= GeometryInfo.fMaterialDataVec[imat].Properties[0]; 
		
		Eyoung = (G*(3*Lambda+2*G))/(Lambda+G);
		PoissonRatio = (Lambda)/(2*(Lambda+G));
		MaterialElastic = new TPZElasticityMaterial(GeometryInfo.fMaterialDataVec[imat].MatID, Eyoung, PoissonRatio, BodyForceX, BodyForceY, planestress); 
		TPZMaterial * Material(MaterialElastic);
		cmesh->InsertMaterialObject(Material);
		
		// Inserting boundary conditions
		for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
		{
			if(GeometryInfo.fMaterialDataVec[imat].MatID == int(GeometryInfo.fBCMaterialDataVec[ibc].Properties[0])) 
			{	
				TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);				
				TPZMaterial * BCond = MaterialElastic->CreateBC(Material,GeometryInfo.fBCMaterialDataVec[ibc].MatID,dirichlet, val1, val2);
				cmesh->InsertMaterialObject(BCond);				
			}
		}
		
		// Inseting mass line sources (Wells)
		int NumberOfEvents= Events.size();
		if (NumberOfEvents !=0 ) 
		{
			int WellsPerEvent = Events[0].Rows();
			for (int itime = 0; itime < NumberOfEvents; itime++) 
			{
				for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
				{
					if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].MatID) {
						if (TStep == int(Events[itime](iwell,0))) 
						{
							TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
							TPZMaterial * BCLINE = MaterialElastic->CreateBC(Material, int(Events[itime](iwell,1)),neumann, val1, val2);
							cmesh->InsertMaterialObject(BCLINE);
						}
					}
				}
				
			}
		}
		
		
	}
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "BuildGroups" ).ToElement();		
	CharContainer = Container->Attribute("NGroups");
	int BuildGroups = atoi(CharContainer);	
	
	if (InterfaceBehavior) 
	{	
		
		// Building computational objects by groups
		std::set<int> iset;
		for (int ngroup = 1; ngroup <= BuildGroups; ngroup++) 
		{
			for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
			{	
				
				if(ngroup == GeometryInfo.fMaterialDataVec[imat].Properties[9]) 
				{	
					iset.insert(GeometryInfo.fMaterialDataVec[imat].MatID);
					// Inseting mass line sources (Wells)
					int NumberOfEvents= Events.size();
					if (NumberOfEvents != 0) {
						int WellsPerEvent = Events[0].Rows();
						for (int itime = 0; itime < NumberOfEvents; itime++) 
						{
							for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
							{
								if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].MatID) {
									if (TStep == int(Events[itime](iwell,0))) 
									{
										iset.insert(int(Events[itime](iwell,1)));
									}
								}
							}
							
						}
					}
					
				}				
				
			}
			
			for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
			{
				if(ngroup == GeometryInfo.fBCMaterialDataVec[ibc].Properties[5]) 
				{	
					iset.insert(GeometryInfo.fBCMaterialDataVec[ibc].MatID);
				}			
				
			}			
			
			cmesh->AutoBuild(iset);
			iset.clear();		
			cmesh->Reference()->ResetReference();
			
		}			
	}
	else 
	{
		cmesh->AutoBuild();
	}
	
	
	// Inserting Interfaces
	for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
	{
		
		if (GeometryInfo.fBCMaterialDataVec[ibc].Properties[0] == InterfacesCutoff) 
		{
			if (InterfaceBehavior) 
			{
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[2];
				
				
				ElasticMatInterface2D *InterfaceMat = new ElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].MatID,Eyoung, PoissonRatio, BodyForceX, BodyForceY, planestress); 
				InterfaceMat->SetPenalty(KNU,KTU);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				cmesh->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				cmesh->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int ngel = cmesh->Reference()->NElements();
				
				for(int igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].MatID) {
						TPZStack < TPZCompElSide > neigh;
						int side = Gel->NSides()-1;
						TPZGeoElSide gelside(Gel,side);
						gelside.EqualLevelCompElementList(neigh, 1, 0);
						int temp = neigh.size();
						if (temp == 0) {
							continue;
						}				
						if (neigh.size() != 2) {
							DebugStop();
						}
						int gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}				
				
			}
			else 
			{
				
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[2];
				
				ElasticMatInterface2D *InterfaceMat = new ElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].MatID,Eyoung, PoissonRatio, BodyForceX, BodyForceY, planestress); 
				InterfaceMat->SetPenalty(KNU,KTU);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				cmesh->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				cmesh->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int ngel = cmesh->Reference()->NElements();
				
				for(int igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].MatID) {
						TPZStack < TPZCompElSide > neigh;
						int side = Gel->NSides()-1;
						TPZGeoElSide gelside(Gel,side);
						gelside.EqualLevelCompElementList(neigh, 1, 0);
						int temp = neigh.size();
						if (temp == 0) {
							continue;
						}				
						if (neigh.size() != 2) {
							DebugStop();
						}
						int gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}
				// Check the interfcae behaviour
				cmesh->Reference()->ResetReference();				
				cmesh->AutoBuild();
				
			}
			
			
		}
		
	}	
	
	
	return cmesh;
	
}

TPZCompMesh * ComputationalDiffusionMesh(TiXmlHandle ControlDoc,TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh, int pOrder,TPZVec < TPZFMatrix<REAL> > Events, int TStep)
{
	
	// TiXmlElement dummy object
	TiXmlElement* Container;
	const char * CharContainer;
	
	// Pure Diffusion problem
	TPZVec <REAL> convdir(3,0.);
	REAL diff	=	0.0;
	REAL conv	=	0.0;
	REAL flux	=	0.0;
	
	// Getting mesh dimension
	int dim = GeometryInfo.fProblemDimension;
	
	// Aproximation Space of order -> pOrder
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravitationalConstant" ).ToElement();		
	CharContainer = Container->Attribute("Gravity");
	REAL Gravitationalconstant = atof(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravityDirection" ).ToElement();		
	CharContainer = Container->Attribute("x-direction");	
	REAL Gravity_Xdirection = atof(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravityDirection" ).ToElement();		
	CharContainer = Container->Attribute("y-direction");	
	REAL Gravity_Ydirection = atof(CharContainer);	
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "FaultInterface" ).FirstChild( "Contact" ).ToElement();		
	CharContainer = Container->Attribute("Interfaces");		
	int InterfacesCutoff = atof(CharContainer);	
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "FaultInterface" ).FirstChild( "Contact" ).ToElement();		
	CharContainer = Container->Attribute("ContacFlag");	
	int InterfaceBehavior = atoi(CharContainer);
	
	
	// Data for elastic problem
	REAL Eyoung			=	0.0;
	REAL PoissonRatio	=	0.0;
	REAL Permeability	=	0.0;
	REAL Viscosity		=	0.0;	
	REAL FluidDensity	=	0.0;	
	REAL BodyForceX		=	0.0;
	REAL BodyForceY		=	0.0;	
	TPZMatPoisson3d * MaterialDiffusion;
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{	
		
		// Using Gravity field
		BodyForceX		= Gravitationalconstant * Gravity_Xdirection * GeometryInfo.fMaterialDataVec[imat].Properties[1];		
		BodyForceY		= Gravitationalconstant * Gravity_Ydirection * GeometryInfo.fMaterialDataVec[imat].Properties[1];
		Permeability	= GeometryInfo.fMaterialDataVec[imat].Properties[3];
		Viscosity		= GeometryInfo.fMaterialDataVec[imat].Properties[4];		
		FluidDensity	= GeometryInfo.fMaterialDataVec[imat].Properties[1]; 
		
		diff = (Permeability)/(Viscosity);
		
		MaterialDiffusion = new TPZMatPoisson3d(GeometryInfo.fMaterialDataVec[imat].MatID, dim);
		MaterialDiffusion->SetParameters(diff, conv, convdir);
		MaterialDiffusion->SetInternalFlux(flux);
		MaterialDiffusion->NStateVariables();		
		TPZMaterial * Material(MaterialDiffusion);
		cmesh->InsertMaterialObject(Material);
		
		// Inserting boundary conditions
		for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
		{
			if(GeometryInfo.fMaterialDataVec[imat].MatID == int(GeometryInfo.fBCMaterialDataVec[ibc].Properties[0])) 
			{	
				TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);				
				TPZMaterial * BCond = MaterialDiffusion->CreateBC(Material,GeometryInfo.fBCMaterialDataVec[ibc].MatID,dirichlet, val1, val2);
				cmesh->InsertMaterialObject(BCond);				
			}
		}
		
		// Inseting mass line sources (Wells)
		int NumberOfEvents= Events.size();
		if (NumberOfEvents !=0 ) 
		{
			int WellsPerEvent = Events[0].Rows();
			for (int itime = 0; itime < NumberOfEvents; itime++) 
			{
				for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
				{
					if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].MatID) {
						if (TStep == int(Events[itime](iwell,0))) 
						{
							TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
							TPZMaterial * BCLINE = MaterialDiffusion->CreateBC(Material, int(Events[itime](iwell,1)),neumann, val1, val2);
							cmesh->InsertMaterialObject(BCLINE);
						}
					}
				}
				
			}
		}		
		
	}
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "BuildGroups" ).ToElement();		
	CharContainer = Container->Attribute("NGroups");
	int BuildGroups = atoi(CharContainer);	
	
	if (InterfaceBehavior) 
	{	
		
		// Building computational objects by groups
		std::set<int> iset;
		for (int ngroup = 1; ngroup <= BuildGroups; ngroup++) 
		{
			for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
			{	
				
				if(ngroup == GeometryInfo.fMaterialDataVec[imat].Properties[9]) 
				{	
					iset.insert(GeometryInfo.fMaterialDataVec[imat].MatID);
					// Inseting mass line sources (Wells)
					int NumberOfEvents= Events.size();
					if (NumberOfEvents != 0) {
						int WellsPerEvent = Events[0].Rows();
						for (int itime = 0; itime < NumberOfEvents; itime++) 
						{
							for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
							{
								if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].MatID) {
									if (TStep == int(Events[itime](iwell,0))) 
									{
										iset.insert(int(Events[itime](iwell,1)));
										cout << "Inserting well with ID " << int(Events[itime](iwell,1)) << " With boundary contition " << ngroup << endl;										
									}
								}
							}
							
						}
					}
					
				}				
				
			}
			for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
			{
				if(ngroup == int(GeometryInfo.fBCMaterialDataVec[ibc].Properties[5])) 
				{	
					
					if (0 < GeometryInfo.fBCMaterialDataVec[ibc].Properties[6]) 
					{
						iset.insert(GeometryInfo.fBCMaterialDataVec[ibc].MatID);
					}
				}			
				
			}			
			
			cmesh->AutoBuild(iset);
			iset.clear();		
			cmesh->Reference()->ResetReference();
			
		}			
	}
	else 
	{
		std::set<int> iset;
		for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
		{	
			
			if (0 < GeometryInfo.fMaterialDataVec[imat].Properties[10]) 
			{
				iset.insert(GeometryInfo.fMaterialDataVec[imat].MatID);
				// Inseting mass line sources (Wells)
				int NumberOfEvents= Events.size();
				if (NumberOfEvents != 0) {
					int WellsPerEvent = Events[0].Rows();
					for (int itime = 0; itime < NumberOfEvents; itime++) 
					{
						for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
						{
							if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].MatID) {
								if (TStep == int(Events[itime](iwell,0))) 
								{
									iset.insert(int(Events[itime](iwell,1)));
									
								}
							}
						}
						
					}
				}				
			}
			
		}
		
		for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
		{
			
			if (0 < GeometryInfo.fBCMaterialDataVec[ibc].Properties[6]) 
			{
				iset.insert(GeometryInfo.fBCMaterialDataVec[ibc].MatID);				
			}
			
			
		}		
		
		cmesh->AutoBuild(iset);
		iset.clear();		
		cmesh->Reference()->ResetReference();	
		
		
		//		cmesh->AutoBuild();
	}	
	
	
	// Inserting Interfaces
	for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
	{
		
		if (GeometryInfo.fBCMaterialDataVec[ibc].Properties[0] == InterfacesCutoff) 
		{
			if (InterfaceBehavior) 
			{
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[2];
				
				
				ElasticMatInterface2D *InterfaceMat = new ElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].MatID,Eyoung, PoissonRatio, BodyForceX, BodyForceY, 0); 
				InterfaceMat->SetPenalty(KNU,KTU);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				cmesh->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				cmesh->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int ngel = cmesh->Reference()->NElements();
				
				for(int igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].MatID) {
						TPZStack < TPZCompElSide > neigh;
						int side = Gel->NSides()-1;
						TPZGeoElSide gelside(Gel,side);
						gelside.EqualLevelCompElementList(neigh, 1, 0);
						int temp = neigh.size();
						if (temp == 0) {
							continue;
						}				
						if (neigh.size() != 2) {
							DebugStop();
						}
						int gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}				
				
			}
			else 
			{
				
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[2];
				
				ElasticMatInterface2D *InterfaceMat = new ElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].MatID,Eyoung, PoissonRatio, BodyForceX, BodyForceY, 0); 
				InterfaceMat->SetPenalty(KNU,KTU);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				cmesh->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				cmesh->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int ngel = cmesh->Reference()->NElements();
				
				for(int igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].MatID) {
						TPZStack < TPZCompElSide > neigh;
						int side = Gel->NSides()-1;
						TPZGeoElSide gelside(Gel,side);
						gelside.EqualLevelCompElementList(neigh, 1, 0);
						int temp = neigh.size();
						if (temp == 0) {
							continue;
						}				
						if (neigh.size() != 2) {
							DebugStop();
						}
						int gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}
				// Check the interfcae behaviour
				cmesh->Reference()->ResetReference();				
				cmesh->AutoBuild();
				
			}
			
			
		}
		
	}	
	
	
	return cmesh;
	
}

TPZCompMesh * ComputationalPoroelasticityMesh(TiXmlHandle ControlDoc, TPZReadGIDGrid GeometryInfo, 
											  TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec ,TPZVec <TPZPoroElastic2d  * > &mymaterial,TPZVec < TPZFMatrix<REAL> > Events,  int TStep)
{	
	// TiXmlElement dummy object
	TiXmlElement* Container;
	const char * CharContainer;
	
	// Plane strain assumption
	int planestress = 0;
	
	// Getting mesh dimension
	int dim = GeometryInfo.fProblemDimension;	
	
	// Aproximation Space continuous	
	gmesh->ResetReference();
	TPZCompMesh * Multiphysics = new TPZCompMesh(gmesh);
	Multiphysics->SetAllCreateFunctionsMultiphysicElem();
	TPZAutoPointer<TPZFunction<STATE> > TimeDepForcingF;
	TPZAutoPointer<TPZFunction<STATE> > TimeDepFExact;
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Analyticfunction" ).FirstChild( "Function" ).ToElement();		
	CharContainer = Container->Attribute("FunctionName");
	//	if (CharContainer == "ExactSolution2DLineSource") {
	TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolutionfiniteColumn1D);
	//	}
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravitationalConstant" ).ToElement();		
	CharContainer = Container->Attribute("Gravity");
	REAL Gravitationalconstant = atof(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravityDirection" ).ToElement();		
	CharContainer = Container->Attribute("x-direction");	
	REAL Gravity_Xdirection = atof(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Gravity" ).FirstChild( "GravityDirection" ).ToElement();		
	CharContainer = Container->Attribute("y-direction");	
	REAL Gravity_Ydirection = atof(CharContainer);	
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "FaultInterface" ).FirstChild( "Contact" ).ToElement();		
	CharContainer = Container->Attribute("Interfaces");		
	int InterfacesCutoff = atof(CharContainer);	
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "FaultInterface" ).FirstChild( "Contact" ).ToElement();		
	CharContainer = Container->Attribute("ContacFlag");	
	int InterfaceBehavior = atoi(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "DimensionlessCalculation" ).FirstChild( "Calculation" ).ToElement();		
	CharContainer = Container->Attribute("CalcualtionFlag");	
	int Dimensionless = atoi(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeControls" ).ToElement();		
	CharContainer = Container->Attribute("DeltaTime");		
	REAL Delta = atof(CharContainer);	
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeControls" ).ToElement();			
	CharContainer = Container->Attribute("Theta");	
	REAL Theta = atof(CharContainer);		
	
	// Data for poroelastic problem
	REAL LambdaU		=	0.0;
	REAL Lambda			=	0.0;
	REAL G				=	0.0;
	REAL alpha			=	0.0;	
	REAL MixtureDensity	=	0.0;
	REAL RockDensity	=	0.0;
	REAL RockPorosity	=	0.0;	
	REAL FluidDensity	=	0.0;
	REAL FluidViscosity	=	0.0;
	REAL Permeability	=	0.0;
	REAL BodyForceX		=	0.0;
	REAL BodyForceY		=	0.0;
	REAL diff			=	0.0;
	REAL S				=	0.0;
	REAL Se				=	0.0;	
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{	
		TPZPoroElastic2d * MaterialPoroElastic;
		LambdaU			= GeometryInfo.fMaterialDataVec[imat].Properties[5];
		Lambda			= GeometryInfo.fMaterialDataVec[imat].Properties[6];		
		G				= GeometryInfo.fMaterialDataVec[imat].Properties[7];
		alpha			= GeometryInfo.fMaterialDataVec[imat].Properties[8];
		RockDensity		= GeometryInfo.fMaterialDataVec[imat].Properties[0];		
		Permeability	= GeometryInfo.fMaterialDataVec[imat].Properties[3];
		FluidViscosity	= GeometryInfo.fMaterialDataVec[imat].Properties[4];
		FluidDensity	= GeometryInfo.fMaterialDataVec[imat].Properties[1];		
		RockPorosity	= GeometryInfo.fMaterialDataVec[imat].Properties[2];		
		diff = (Permeability)/(FluidViscosity);
		MixtureDensity = (1 - RockPorosity) * RockDensity + RockPorosity * FluidDensity;
		// Using Gravity field
		BodyForceX		= Gravitationalconstant * Gravity_Xdirection * MixtureDensity;	
		BodyForceY		= Gravitationalconstant * Gravity_Ydirection * MixtureDensity;
		
		if (Lambda != LambdaU) 
		{
			S = ((pow(alpha,2))/((LambdaU-Lambda)))*((LambdaU+2.0*G)/(Lambda+2.0*G));
			Se = (pow(alpha,2))/((LambdaU-Lambda));
			
		}
		else
		{
			S = (pow(alpha,2)/(Lambda+2.0*G));
			Se = 0;
			
			
		}		
		
		if (Dimensionless) 
		{
			Lambda			= Lambda*S;
			G				= G*S;
			Se				= Se/S;
			Permeability	= 1.0;
			FluidViscosity	= 1.0;
			
		}	
		
		MaterialPoroElastic = new TPZPoroElastic2d (GeometryInfo.fMaterialDataVec[imat].MatID, dim);
		MaterialPoroElastic->SetParameters(Lambda, G,BodyForceX,BodyForceY);
		MaterialPoroElastic->SetParameters(Permeability,FluidViscosity);
		MaterialPoroElastic->SetfPlaneProblem(planestress);
		MaterialPoroElastic->SetBiotParameters(alpha,Se);
		MaterialPoroElastic->SetTimeStep(Delta,Theta);
		MaterialPoroElastic->SetTimeDependentFunctionExact(TimeDepFExact);			
		TPZMaterial * Material(MaterialPoroElastic);
		Multiphysics->InsertMaterialObject(Material);
		
		
		// Inserting boundary conditions
		for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
		{
			if(GeometryInfo.fMaterialDataVec[imat].MatID == int(GeometryInfo.fBCMaterialDataVec[ibc].Properties[0])) 
			{	
				TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
				val2(0,0)=GeometryInfo.fBCMaterialDataVec[ibc].Properties[2];
				val2(1,0)=GeometryInfo.fBCMaterialDataVec[ibc].Properties[3];			
				val2(2,0)=GeometryInfo.fBCMaterialDataVec[ibc].Properties[4];
				TPZMaterial * BCond = MaterialPoroElastic->CreateBC(Material,GeometryInfo.fBCMaterialDataVec[ibc].MatID,GeometryInfo.fBCMaterialDataVec[ibc].Properties[1], val1, val2);
				Multiphysics->InsertMaterialObject(BCond);				
			}
		}
		
		// Inseting mass line sources (Wells)
		int NumberOfEvents= Events.size();
		if (NumberOfEvents != 0) {
			int WellsPerEvent = Events[0].Rows();
			for (int itime = 0; itime < NumberOfEvents; itime++) 
			{
				for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
				{
					if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].MatID) {
						if (TStep == int(Events[itime](iwell,0))) 
						{
							TPZFMatrix<REAL> val1(3,2,Events[itime](iwell,4)), val2(3,1,0.);
							val2(0,0)=Events[itime](iwell,5);
							val2(1,0)=Events[itime](iwell,6);			
							val2(2,0)=Events[itime](iwell,7);
							TPZMaterial * BCLINE = MaterialPoroElastic->CreateBC(Material, int(Events[itime](iwell,1)),int(Events[itime](iwell,3)), val1, val2);
							cout << "Inserting well with ID " << int(Events[itime](iwell,1)) << " With boundary contition " << int(Events[itime](iwell,7)) << endl;
							Multiphysics->InsertMaterialObject(BCLINE);
						}
					}
				}
				
			}		
		}
		
		mymaterial[imat] = MaterialPoroElastic;
		
	}	
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "BuildGroups" ).ToElement();		
	CharContainer = Container->Attribute("NGroups");
	int BuildGroups = atoi(CharContainer);	
	
	if (InterfaceBehavior) 
	{	
		
		// Building computational objects by groups
		std::set<int> iset;
		for (int ngroup = 1; ngroup <= BuildGroups; ngroup++) 
		{
			for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
			{	
				
				if(ngroup == GeometryInfo.fMaterialDataVec[imat].Properties[9]) 
				{	
					iset.insert(GeometryInfo.fMaterialDataVec[imat].MatID);
					// Inseting mass line sources (Wells)
					int NumberOfEvents= Events.size();
					if (NumberOfEvents != 0) {
						int WellsPerEvent = Events[0].Rows();
						for (int itime = 0; itime < NumberOfEvents; itime++) 
						{
							for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
							{
								if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].MatID) {
									if (TStep == int(Events[itime](iwell,0))) 
									{
										iset.insert(int(Events[itime](iwell,1)));
									}
								}
							}
							
						}
					}
					
				}				
				
			}
			
			for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
			{
				if(ngroup == GeometryInfo.fBCMaterialDataVec[ibc].Properties[5]) 
				{	
					iset.insert(GeometryInfo.fBCMaterialDataVec[ibc].MatID);
				}			
				
			}			
			
			Multiphysics->AutoBuild(iset);
			iset.clear();		
			Multiphysics->Reference()->ResetReference();
			
		}			
	}
	else 
	{
		Multiphysics->AutoBuild();
	}
	
	
	// Inserting Interfaces
	for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
	{
		
		if (GeometryInfo.fBCMaterialDataVec[ibc].Properties[0] == InterfacesCutoff) 
		{
			if (InterfaceBehavior) 
			{
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[2];
				REAL KNP = GeometryInfo.fBCMaterialDataVec[ibc].Properties[3];
				REAL KTP = GeometryInfo.fBCMaterialDataVec[ibc].Properties[4];
				
				PoroElasticMatInterface2D *InterfaceMat = new PoroElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].MatID,1); 
				InterfaceMat->SetPenalty(KNU,KTU,KNP,KNP);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				Multiphysics->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				Multiphysics->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int ngel = Multiphysics->Reference()->NElements();
				
				for(int igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].MatID) {
						TPZStack < TPZCompElSide > neigh;
						int side = Gel->NSides()-1;
						TPZGeoElSide gelside(Gel,side);
						gelside.EqualLevelCompElementList(neigh, 0, 0);
						int temp = neigh.size();
						if (temp == 0) {
							continue;
						}				
						if (neigh.size() != 2) {
							DebugStop();
						}
						int gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZMultiphysicsInterfaceElement(* Multiphysics, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}			
				
			}
			else 
			{
				
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].Properties[2];
				REAL KNP = GeometryInfo.fBCMaterialDataVec[ibc].Properties[3];
				REAL KTP = GeometryInfo.fBCMaterialDataVec[ibc].Properties[4];
				
				PoroElasticMatInterface2D *InterfaceMat = new PoroElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].MatID,1); 
				InterfaceMat->SetPenalty(KNU,KTU,KNP,KNP);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				Multiphysics->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				Multiphysics->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int ngel = Multiphysics->Reference()->NElements();
				
				for(int igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].MatID) {
						TPZStack < TPZCompElSide > neigh;
						int side = Gel->NSides()-1;
						TPZGeoElSide gelside(Gel,side);
						gelside.EqualLevelCompElementList(neigh, 0, 0);
						int temp = neigh.size();
						if (temp == 0) {
							continue;
						}				
						if (neigh.size() != 2) {
							DebugStop();
						}
						int gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZMultiphysicsInterfaceElement(* Multiphysics, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}
				// Check the interface behaviour
				Multiphysics->Reference()->ResetReference();				
				Multiphysics->AutoBuild();
				
			}
			
			
		}
		
	}		
	
	Multiphysics->AdjustBoundaryElements();
	Multiphysics->CleanUpUnconnectedNodes();
	
	// Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, Multiphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,Multiphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, Multiphysics);
	
	ofstream arg("Mphysic.txt");
	Multiphysics->Print(arg);
	
	return Multiphysics;
}

#define VTK
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//	ofstream file("Solution.out");
	//	an.Solution().Print("solution", file); 
}



void PostProcessDiffusion(TPZAnalysis &an, std::string plotfile)
{
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
	vecnames[0]= "MinusKGradU";
	
	const int dim = 2;
	int div = 2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void PostProcessElasticity(TPZAnalysis &an, std::string plotfile)
{
	TPZManVector<std::string,10> scalnames(3), vecnames(1);
	scalnames[0] = "SigmaX";
	scalnames[1] = "SigmaY";
	scalnames[2] = "Pressure";	
	vecnames[0]= "displacement";
	
	
	const int dim = 2;
	int div = 2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void PostProcessPoroeasticity(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
{
	TPZBuildMultiphysicsMesh * Objectdumy;
	Objectdumy->TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(12), vecnames(3);
	scalnames[0] = "SigmaX";
	scalnames[1] = "SigmaY";
	scalnames[2] = "TauXY";	
	scalnames[3] = "DisplacementX";
	scalnames[4] = "DisplacementY";
	scalnames[5] = "SolidPressure";
	scalnames[6] = "FluidPressure";
	scalnames[7] = "EDisplacementX";
	scalnames[8] = "EDisplacementY";
	scalnames[9] = "EPressure";
	scalnames[10] = "ESIGX";
	scalnames[11] = "ESIGY";
	vecnames[0]  = "Displacement";
	vecnames[1]  = "FluxVector";
	vecnames[2]  = "EFluxVector";
	
	
	const int dim = 1;
	int div =0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}


void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh)
{
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gMesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
		}//for i
	}//ref
}

void RefinamentoUniforme(TPZGeoMesh * gMesh, int nh, int MatId, int indexEl)
{
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gMesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2 || gel->Dimension() == 1){
				if (gel->MaterialId()== MatId && gel-> Index()==indexEl){
					gel->Divide (filhos);
				}
			}
		}//for i
	}//ref
}

void RefinElemComp(TPZCompMesh  *cMesh, int indexEl)
{
	
	TPZVec<int > subindex; 
	int nel = cMesh->ElementVec().NElements(); 
	for(int el=0; el < nel; el++){
		TPZCompEl * compEl = cMesh->ElementVec()[el];
		if(!compEl) continue;
		int ind = compEl->Index();
		if(ind==indexEl){
			compEl->Divide(indexEl, subindex, 1);
		}
	}	
}

void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
{
	
	TPZVec<int > subindex;
	for (int iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		int nel = elvec.NElements(); 
		for(int el=0; el < nel; el++){
			TPZCompEl * compEl = elvec[el];
			if(!compEl) continue;
			int ind = compEl->Index();
			compEl->Divide(ind, subindex, 0);
		}
	}
}

void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	std::stringstream node, connectivity, type;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{        
		if(gmesh->ElementVec()[el]->Type() == EPoint)//Exclude Lines and Arc3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->HasSubElement())
		{
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}            
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType = -1;
		switch (gmesh->ElementVec()[el]->Type())
		{
			case (ETriangle):
			{
				elType = 5;
				break;                
			}
			case (EQuadrilateral ):
			{
				elType = 9;
				break;                
			}
			case (ETetraedro):
			{
				elType = 10;
				break;                
			}
			case (EPiramide):
			{
				elType = 14;
				break;                
			}
			case (EPrisma):
			{
				elType = 13;
				break;                
			}
			case (ECube):
			{
				elType = 12;
				break;                
			}
			default:
			{
				//ElementType NOT Found!!!
				DebugStop();
				break;    
			}
		}
		
		type << elType << std::endl;
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();
	
	file.close();
}


void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file)
{
	TPZGeoMesh *gmesh;
	PrintGMeshVTK(gmesh, file);
}




TPZAutoPointer <TPZMatrix<REAL> > MassMatrix(TPZReadGIDGrid GeometryInfo,TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics)
{
	
	
	TPZSpStructMatrix matsp(mphysics);
	std::set< int > materialid;
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{
		mymaterial[imat]->SetLastState();
		materialid.insert(GeometryInfo.fMaterialDataVec[imat].MatID );
	}	
	
	matsp.SetMaterialIds(materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<REAL> Un;
	TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
	
	return matK2;
	
}


void StiffMatrixLoadVec(TPZReadGIDGrid GeometryInfo,int Nthreads, TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec)
{
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{	
		mymaterial[imat]->SetCurrentState();
	}		
	
	//TPZFStructMatrix matsk(mphysics);
	TPZSkylineStructMatrix matsk(mphysics);
	an.SetStructuralMatrix(matsk); 
	//	an.StructMatrix()->SetNumThreads(Nthreads);
	TPZStepSolver<REAL> step; 
	step.SetDirect(ELDLt); 
	an.SetSolver(step); 
	an.Assemble(); 
	matK1 = an.StructMatrix();
	fvec = an.Rhs();
	
}



void SolveSistTransient(TPZVec < TPZFMatrix<REAL> > Events,TiXmlHandle ControlDoc,int Nthreads, TPZVec <REAL> &PrintStep, std::string FileName, REAL deltaT,REAL maxTime, TPZFMatrix<REAL> &InitialSolution, TPZVec <TPZPoroElastic2d  * > &mymaterial ,
						TPZAnalysis &an, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics , TPZReadGIDGrid GeometryInfo)
{
	
	
	// Calculation of Mass Matrix
	cout << "Calculate Mass Matrix " << endl;    
	TPZAutoPointer <TPZMatrix<REAL> > PoroelasticMassMatrix = MassMatrix(GeometryInfo,mymaterial, mphysics);
	cout << "Calculate Mass Matrix was done! " << endl;	
	
	// Calculation of Stiffness Matrix and Load Vector
	TPZFMatrix<REAL> PoroelasticStiffnessMatrix, PoroelasticStiffnessMatrixInverse;	
	TPZFMatrix<REAL> PoroelasticLoadVector;
	cout << "Calculate Stiffness Matrix and Load Vector " << endl;	
	StiffMatrixLoadVec(GeometryInfo,Nthreads,mymaterial, mphysics, an, PoroelasticStiffnessMatrix, PoroelasticLoadVector);
	cout << "Calculate Stiffness Matrix and Load Vector was done! " << endl;	
	
	// TiXmlElement dummy object
	TiXmlElement* Container;
	const char * CharContainer;
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeControls" ).ToElement();			
	CharContainer = Container->Attribute("Theta");	
	REAL Theta = atof(CharContainer);		
	int NEvents = Events.size();	
	
	int nrows;
	nrows = PoroelasticMassMatrix->Rows();
	TPZFMatrix<REAL> TotalRhs(nrows,1,0.0);
	TPZFMatrix<REAL> RhsTemporal(nrows,1,0.0);
	TPZFMatrix<REAL> LastSolution = InitialSolution;
	
	
	
	REAL	TimeValue	= 0.0;
	int		cent		= 1;
	int		control		= 0;
	
	TimeValue = cent*deltaT; 
	
	while (TimeValue <= maxTime)
	{	
		
		for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
		{
			mymaterial[imat]->SetTimeValue(TimeValue);
		}
		PoroelasticMassMatrix->Multiply(LastSolution,RhsTemporal);
		
		TotalRhs = PoroelasticLoadVector + RhsTemporal;
		an.Rhs() = TotalRhs;
		an.Solve(); 
		LastSolution = an.Solution();
		
		
		// Information Print Control
		
		if (control == PrintStep.size()) 
		{
			control = 0;
		}				
		if (abs(PrintStep[control] - TimeValue) < 1.0e-8) 
		{
			std::stringstream outputfiletemp;
			outputfiletemp << FileName << ".vtk";
			std::string plotfile = outputfiletemp.str();
			PostProcessPoroeasticity(meshvec,mphysics,an,plotfile);
			control++;
			
		}
		
		if (NEvents != 0) 
		{
			for (int ievent = 0; ievent < NEvents; ievent++) 
			{
				
				if (cent == Events[ievent](0,0) ) 
				{
					cout << "Modifiying Poroelastic load Vector ...." << endl;
					TPZCompMesh * DumpMphysics = ComputationalPoroelasticityMesh(ControlDoc, GeometryInfo, meshvec[0]->Reference(), meshvec, mymaterial, Events, cent);
					TPZAnalysis DumpAnalysis(DumpMphysics);
					for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
					{
						mymaterial[imat]->SetTimeStep(deltaT,Theta);
						mymaterial[imat]->SetTimeValue(TimeValue);
						mymaterial[imat]->SetCurrentState();
					}
					
					TPZSkylineStructMatrix matsk(DumpMphysics);
					DumpAnalysis.SetStructuralMatrix(matsk);
					TPZStepSolver<REAL> step;
					DumpAnalysis.SetSolver(step);
					DumpAnalysis.AssembleResidual();
					PoroelasticLoadVector = DumpAnalysis.Rhs();
					cout << "Modifiying Poroelastic done! at time step = " << cent << endl;				
				}
			}
		}
		

		
		
		
		cent++;
		TimeValue = cent*deltaT;
	}
	
}

void *InitalPressureCalculations(REAL &PressureValue,TPZVec <REAL> &PointCoordinates)
{
	PressureValue = 0.0;
}

void SetInitialConditions(TPZAnalysis &AnalysisToInitialize, TPZCompMesh * ComputationalMesh, TPZFMatrix<REAL> &ComputationalMeshInitialSolution, void ( *PointerToFunction(REAL &,TPZVec <REAL> &)))
{
	
	int RowsNumber = AnalysisToInitialize.Solution().Rows();
	ComputationalMeshInitialSolution.Resize(RowsNumber,1);
	
	for (int i = 0; i < RowsNumber; i++) 
	{
		ComputationalMeshInitialSolution(i,0) = 0.0;
	}
	
	int TotalComputationalElements	= ComputationalMesh->NElements();		
	for ( int el = 0; el < TotalComputationalElements ; el++ )
	{
		
		TPZCompEl * ComputationalElement = ComputationalMesh->ElementVec()[el];
		TPZGeoEl * GeoemtricElement;
		GeoemtricElement = ComputationalElement->Reference();
		int NumberNodeCorners = ComputationalElement->Reference()->NCornerNodes();
		
		for (int ConnectIndex = 0; ConnectIndex < NumberNodeCorners ; ConnectIndex++) 
		{
			TPZConnect & ThisConnect = ComputationalElement->Connect(ConnectIndex);
			TPZVec <REAL> Coordenates(3,0.0);
			REAL InitialValue;			
			Coordenates[0] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(0);
			Coordenates[1] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(1);
			Coordenates[2] = GeoemtricElement->NodePtr(ConnectIndex)->Coord(2);
			
			PointerToFunction(InitialValue, Coordenates);			
			
			int Sequence = ThisConnect.SequenceNumber();
			int BlockPos = ComputationalMesh->Block().Position(Sequence);
			
			ComputationalMeshInitialSolution(BlockPos,0)= InitialValue;
		}
	}	
}

TPZCompMesh *SetInitialConditionsbyL2Projection(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &InitialCondition, bool &CustomFunction)
{
	
	int Dimension = 2, NumberStateVariables = 1;
	TPZL2Projection * MyProjection;
	//	This create material for projection	
	MyProjection = new TPZL2Projection(matId, Dimension, NumberStateVariables, InitialCondition, pOrder);
	
	TPZCompMesh * ComputationalMeshMyProjection = new TPZCompMesh(gmesh);
	ComputationalMeshMyProjection->SetDimModel(Dimension);
	ComputationalMeshMyProjection->SetAllCreateFunctionsContinuous();
	//	ComputationalMeshMyProjection->SetAllCreateFunctionsDiscontinuous();	
	ComputationalMeshMyProjection->SetDefaultOrder(pOrder);
	
	
	TPZMaterial * mat(MyProjection);
	if (CustomFunction) {
		TPZAutoPointer < TPZFunction<STATE> > InitialPressureCalculations;
		InitialPressureCalculations = new TPZDummyFunction<STATE> (InitialPressureDistribution);
		MyProjection->SetForcingFunction(InitialPressureCalculations);	
	}
	
	ComputationalMeshMyProjection->InsertMaterialObject(mat);
	
	
	
	//    ///inserir connect da pressao
	//    int ncon = cmesh->NConnects();
	//    for(int i=0; i<ncon; i++)
	//    {
	//        TPZConnect &newnod = cmesh->ConnectVec()[i];
	//        newnod.SetPressure(true);
	//    }
	//    
	//    ///set order total da shape
	//    int nel = cmesh->NElements();
	//    for(int i=0; i<nel; i++){
	//        TPZCompEl *cel = cmesh->ElementVec()[i];
	//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
	//        celdisc->SetConstC(1.);
	//        celdisc->SetCenterPoint(0, 0.);
	//        celdisc->SetCenterPoint(1, 0.);
	//        celdisc->SetCenterPoint(2, 0.);
	//        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
	//        {
	//            if(triang==true) celdisc->SetTotalOrderShape();
	//            else celdisc->SetTensorialShape();
	//        }
	//    }
	//    
	//    
	//#ifdef DEBUG   
	//    int ncel = cmesh->NElements();
	//    for(int i =0; i<ncel; i++){
	//        TPZCompEl * compEl = cmesh->ElementVec()[i];
	//        if(!compEl) continue;
	//        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
	//        if(facel)DebugStop();
	//        
	//    }
	//#endif
	ComputationalMeshMyProjection->AutoBuild();    
	return ComputationalMeshMyProjection;
}




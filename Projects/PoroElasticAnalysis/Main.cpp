#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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
#include "AnalyticalFunctions.h"
#include "ElasticMatInterface2D.h"
#include "PoroElasticMatInterface2D.h"
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
#include "pzporoelastic2d.h"
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


//	// Using Log4cXX as logging tool
//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.poroelastic2d"));
//#endif
//
//#ifdef LOG4CXX
//static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic.data"));
//#endif

//#undef REFPATTERNDIR
//#define REFPATTERNDIR "RefPatterns/fgfchghghg"

using namespace std;
using namespace pzgeom;	

// Dymmy Boundary Conditions
const int dirichlet = 0;
const int neumann = 1;


// Defintions of Implemented Methods

//	This Create Computational Meshes
TPZCompMesh *ComputationalDiffusionMesh(TiXmlHandle ControlDoc,TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh,int pOrder,TPZVec < TPZFMatrix<REAL> > Events, int TStep);
TPZCompMesh *ComputationalElasticityMesh(TiXmlHandle ControlDoc,TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh,int pOrder,TPZVec < TPZFMatrix<REAL> > Events, int TStep);
TPZCompMesh *ComputationalPoroelasticityMesh(TiXmlHandle ControlDoc, TPZReadGIDGrid &GeometryInfo, TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,TPZVec <TPZPoroElastic2d  * > &mymaterial,TPZVec < TPZFMatrix<REAL> > Events, int TStep);
TPZCompMesh *SetInitialConditionsbyL2Projection(TiXmlHandle ControlDoc, TPZReadGIDGrid GeometryInfo, TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &InitialCondition, bool &CustomFunction);

//	This Set intial Conditions
void SetInitialConditions(TPZAnalysis &AnalysisToInitialize, TPZCompMesh * ComputationalMesh, TPZFMatrix<REAL> &ComputationalMeshInitialSolution, void ( *PointerToFunction(REAL &,TPZVec <REAL> &)));
void *InitalPressureCalculations(REAL &PressureValue, TPZVec <REAL> &PointCoordinates);
void *InitialDisplacementCalculations();

//	This Solve Different analysis
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);
TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZReadGIDGrid GeometryInfo,TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics);
void SolveSistTransient(bool IsInitial,TPZVec < TPZFMatrix<REAL> > Events,TiXmlHandle ControlDoc,int Nthreads, TPZVec <REAL> &PrintEachtimeStep, std::string FileName, REAL deltaT,REAL maxTime, TPZFMatrix<STATE> &InitialSolution, 
						TPZVec <TPZPoroElastic2d  * > &mymaterial , TPZAnalysis &an, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZReadGIDGrid & GeometryInfo);	
void StiffMatrixLoadVec(TPZReadGIDGrid GeometryInfo,int Nthreads, TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec);

//	These are tools for spatial and polynomial refinement and Postprocess of solutions
void PostProcessElasticity(TPZAnalysis &an, std::string plotfile);
void PostProcessDiffusion(TPZAnalysis &an, std::string plotfile);
void PostProcessPoroeasticity(TiXmlHandle ControlDoc,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile, int dim);
void UniformRefinement(TPZGeoMesh  *gMesh, int nh);
void UniformRefinement(TPZGeoMesh *gMesh, int nh, int MatId);
void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);

int main(int argc, char *argv[])
{

	//	Reading arguments 
	char *Filetoload = NULL;
	{
		using namespace std;
		if (argc != 2)
		{
			cout << "Size: " << argc << " Number of Arguments " << endl;
			cout << "Usage: " << argv[0] << " MyContorlFile.xml " << endl;
			cout <<	"Program stop: not xml file found \n" << endl;
			DebugStop();
		}
		
		if (argc == 2)
		{
			cout << "Control File used : " << argv[1] << "\n" << endl;
			Filetoload		= argv[1];
		}
	}

	TiXmlDocument doc(Filetoload);	
	bool loadOkay = false;
	loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		cout << "This Xml is ok! ->> " << Filetoload << endl;
	}
	else
	{
		cout << "Failed to load file \n" << Filetoload << endl;
		cout <<	"Check your xml structure. \n" << endl;
		DebugStop();
	}	
	
	// Accesing xml data
	
	// TiXmlElement dummy object
	TiXmlElement* Container;
	const char * CharContainer;	
	TiXmlHandle docHandle( &doc );
	
	int returnvalue;
		
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "InputDir" ).ToElement();
	CharContainer = Container->Attribute("Name");
	returnvalue = system(CharContainer);
	cout << "Creating directory using : "<< CharContainer << "\n" << endl;
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "OutputDir" ).ToElement();
	CharContainer = Container->Attribute("Name");	
	returnvalue = system(CharContainer);
	cout << "Creating directory using : "<< CharContainer << "\n" << endl;	
	
	//	Creation of Mesh info folder
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "DumpFolder" ).ToElement();
	CharContainer = Container->Attribute("Name");
	returnvalue = system(CharContainer);
	cout << "Creating directory using : "<< CharContainer << "\n" << endl;	

	TiXmlElement* GridDumpName = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "GridDumpFile" ).ToElement();
	const char * GridFileName = GridDumpName->Attribute("GridFile");
	
	TiXmlElement* GridDumpNameini = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "GridDumpFileini" ).ToElement();
	const char * GridFileNameini = GridDumpNameini->Attribute("GridFile");

	TiXmlElement* LmaxDimension = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "MaxGeometricDim" ).ToElement();
	const char * LDimension = LmaxDimension->Attribute("L");	
	
	// GEOMETRICAL MESH CREATION	
	
	//	Polynomial degree approach
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Porder" ).FirstChild( "Elasticity" ).ToElement();			
	CharContainer = Container->Attribute("order");	
	int PElasticity = atoi(CharContainer);
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Porder" ).FirstChild( "Diffusion" ).ToElement();			
	CharContainer = Container->Attribute("order");	
	int PDiffusion = atoi(CharContainer);		
	
	TPZReadGIDGrid GeometryInfoini;
	GeometryInfoini.SetfDimensionlessL(atof(LDimension));
	TPZGeoMesh * gmeshini = NULL;
    gmeshini = GeometryInfoini.GeometricGIDMesh(GridFileNameini);
	
	TPZReadGIDGrid GeometryInfo;
	GeometryInfo.SetfDimensionlessL(atof(LDimension));
	TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);	

	{
		//	Print Geometrical Base Mesh
		ofstream argument("DumpFolder/GeometicMesh.txt");
		gmesh->Print(argument);
		ofstream Dummyfile("DumpFolder/GeometricMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
	}

	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Calculations" ).FirstChild( "Threads" ).ToElement();			
	CharContainer = Container->Attribute("Nthreads");	
	int NumOfThreads = atoi(CharContainer);	
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "RefUniform" ).ToElement();			
	CharContainer = Container->Attribute("URef");	
	int Href = atoi(CharContainer);	
	
	if (NumOfThreads == 0 || Href < 0) 
	{
		NumOfThreads = 1;
		Href = 0;
		cout << "Using default number of trheads	:" << 1 << endl;
		cout << "Using default refinement			:" << 0 << endl;		
	}	
	
	UniformRefinement(gmesh, Href);
	
	//  Used this for Refinement at specific materials
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "RefMaterials" ).ToElement();			
	CharContainer = Container->Attribute("Use");	
	bool UseRefMat = atoi(CharContainer);	
	
	if (UseRefMat) 
	{
        Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "RefMaterials" ).FirstChild( "Mat" ).ToElement();
		for( ; Container; Container=Container->NextSiblingElement())
		{
			int HrefMat = atoi(Container->Attribute( "URef"));
			int Id = atoi(Container->Attribute( "Id"));	
			UniformRefinement(gmesh, HrefMat,Id);
		}
	}
	
	//  Used this for Directional Refinement

	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "GridDirectionalRefinement" ).ToElement();			
	CharContainer = Container->Attribute("Use");	
	bool UseDirectional = atoi(CharContainer);
	
	if (UseDirectional) 
	{

		std::string Path = "RefPatterns";
		gRefDBase.ImportRefPatterns(Path);
		ofstream RefinElPatterns("DumpFolder/RefElPatterns.txt");
		gRefDBase.WriteRefPatternDBase(RefinElPatterns); 	
		gRefDBase.ReadRefPatternDBase("DumpFolder/RefElPatterns.txt");		
		set<int> SetMatsRefDir;
		int ndirectdivp;
		Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "GridDirectionalRefinement" ).FirstChild( "List" ).ToElement();
		for( ; Container; Container=Container->NextSiblingElement())
		{
			ndirectdivp = atoi(Container->Attribute( "NDiv"));
			int Id = atoi(Container->Attribute( "Id"));			
			SetMatsRefDir.insert(Id);
			for(int j = 0; j < ndirectdivp; j++)
			{
				int nel = gmesh->NElements();
				for (int iref = 0; iref < nel; iref++)
				{
					TPZVec<TPZGeoEl*> filhos;
					TPZGeoEl * gelP2 = gmesh->ElementVec()[iref];
					if(!gelP2 || gelP2->HasSubElement()) continue;
					TPZRefPatternTools::RefineDirectional(gelP2, SetMatsRefDir);
					
				}
			}
			SetMatsRefDir.clear();
		}
	}
	
	
	{
		//	Print Geometrical refined Base Mesh
		ofstream argument("DumpFolder/RefinedGeometricMesh.txt");
		gmesh->Print(argument);
		ofstream Dummyfile("DumpFolder/RefinedGeometricMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);	
	}	
	
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
		for( ; Container; Container=Container->NextSiblingElement())
		{
			
			int TemptimeStep = atoi(Container->Attribute( "TimeStep"));
			int iAlter = 0;
			TiXmlElement* LineSource_el = Container->FirstChild("LineSource")->ToElement();
			for( ; LineSource_el; LineSource_el=LineSource_el->NextSiblingElement())
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
			Events[iEvent] = ALter; 
			iEvent++;
		}	
	}
	
	
	//	Intial Calculations
	int TimeStep = -1;	
	
	//	First computational mesh
	
	TPZCompMesh * ComputationalMeshElasticity = ComputationalElasticityMesh(docHandle,GeometryInfo, gmesh, PElasticity, Events, TimeStep);
	//	Print First computational mesh
	ofstream ArgumentElasticity("DumpFolder/ComputationalMeshForElasticity.txt");
	ComputationalMeshElasticity->Print(ArgumentElasticity);
	
	//	Second computational mesh
	TPZCompMesh * ComputationalMeshDiffusion = ComputationalDiffusionMesh(docHandle,GeometryInfo, gmesh, PDiffusion, Events, TimeStep);
	//	Print Second computational mesh
	ofstream ArgumentDiffusion("DumpFolder/ComputationalMeshForDiffusion.txt");
	ComputationalMeshDiffusion->Print(ArgumentDiffusion);	
	
	//	Cleaning reference of the geometric mesh for cmesh1
	gmesh->ResetReference();
	ComputationalMeshElasticity->LoadReferences();
	
	//	Using Uniform Refinement for first computational problem

	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "RefCMeshElasticity" ).ToElement();			
	CharContainer = Container->Attribute("URef");	
	int RefForElasticity = atoi(CharContainer);	
	
	RefinUniformElemComp(ComputationalMeshElasticity,RefForElasticity);
	ComputationalMeshElasticity->AdjustBoundaryElements();
	ComputationalMeshElasticity->CleanUpUnconnectedNodes();
	
	ofstream ArgumentElasticityRef("DumpFolder/ComputationalMeshForElasticityRef.txt");
	ComputationalMeshElasticity->Print(ArgumentElasticityRef);
	ofstream ArgumentElasticityGeoRef("DumpFolder/GeometricMeshForElasticity.txt");
	gmesh->Print(ArgumentElasticityGeoRef);
	
	// Vtk visualization Elasticity Geometric Mesh
	ofstream DummyfileElasticity("DumpFolder/GeometricMeshForElasticity.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, DummyfileElasticity,true);
	
	// Clear reference of the geometric mesh to cmesh2
	gmesh->ResetReference();
	ComputationalMeshDiffusion->LoadReferences();
	
	//	Using Uniform Refinement for second computational problem
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "Grid" ).FirstChild( "RefCMeshDiffusion" ).ToElement();			
	CharContainer = Container->Attribute("URef");	
	int RefForDiffusion = atoi(CharContainer);
	
	RefinUniformElemComp(ComputationalMeshDiffusion,RefForDiffusion);
	ComputationalMeshDiffusion->AdjustBoundaryElements();
	ComputationalMeshDiffusion->CleanUpUnconnectedNodes();
	//	ComputationalMeshDiffusion->ExpandSolution();
	
	ofstream ArgumentDiffusionRef("DumpFolder/ComputationalMeshForDiffusionRef.txt");
	ComputationalMeshDiffusion->Print(ArgumentDiffusionRef);
	ofstream ArgumentDiffusionGeoRef("DumpFolder/GeometricMeshForDiffusion.txt");
	gmesh->Print(ArgumentDiffusionGeoRef);
	
	// Vtk visualization Diffusion Geometric Mesh	
	ofstream DummyfileDiffusion("DumpFolder/GeometricMeshForDiffusion.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, DummyfileDiffusion,true);
	
	// Visualization of computational meshes
	TPZAnalysis ElasticAnalysis(ComputationalMeshElasticity);
	TPZAnalysis DiffusionAnalysis(ComputationalMeshDiffusion);
	
	std::string ElasticityOutput;
	ElasticityOutput = "DumpFolder/ComputationalMeshElasticity";
	std::stringstream ElasticityOutputfiletemp;
	ElasticityOutputfiletemp << ElasticityOutput << ".vtk";
	std::string ElasticityPlotfile = ElasticityOutputfiletemp.str();
	PostProcessElasticity(ElasticAnalysis, ElasticityPlotfile);
	
	std::string DiffusionOutput;
	DiffusionOutput = "DumpFolder/ComputationalMeshDiffusion";
	std::stringstream DiffusionOutputfiletemp;
	DiffusionOutputfiletemp << DiffusionOutput << ".vtk";
	std::string DiffusionPlotfile = DiffusionOutputfiletemp.str();	
	PostProcessDiffusion(DiffusionAnalysis, DiffusionPlotfile);	
	
	TPZVec<TPZCompMesh *> ComputationalMeshVectorInitial(2);
	ComputationalMeshVectorInitial[0] = ComputationalMeshElasticity;
	ComputationalMeshVectorInitial[1] = ComputationalMeshDiffusion;
	
	TPZVec<TPZCompMesh *> ComputationalMeshVector(2);
	ComputationalMeshVector[0] = ComputationalMeshElasticity;
	ComputationalMeshVector[1] = ComputationalMeshDiffusion;	
	
	TPZVec <TPZPoroElastic2d *>  materialist(GeometryInfo.MatNumber,0);	
	
	TPZCompMesh * ComputationalMeshPoroelasticityInitial = ComputationalPoroelasticityMesh(docHandle,GeometryInfoini, gmesh,ComputationalMeshVectorInitial,materialist,Events,TimeStep);
	TPZCompMesh * ComputationalMeshPoroelasticityReCurrent = ComputationalPoroelasticityMesh(docHandle,GeometryInfo, gmesh,ComputationalMeshVector,materialist,Events,TimeStep);	
	
	//--------------------------------------------------------------------------------------------------------------------------------//	
	//--------------------------------------------------------------------------------------------------------------------------------//	
	
	// time control
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "FileName" ).ToElement();		
	CharContainer = Container->Attribute("ProblemName");
	std::string	OutFileName(CharContainer);
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "TimeData" ).FirstChild( "TimeControls" ).ToElement();		
	CharContainer = Container->Attribute("Initial");	
//	REAL InitialTime = atof(CharContainer);
	
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
//	REAL TimeValue	=	0.0;
	TPZVec <REAL> PrintStep(NTimeValues,0);
	int itime = 0;
	
	TiXmlDocument Times(CharContainer);	
	loadOkay = false;
	loadOkay = Times.LoadFile();
	if (loadOkay)
	{
		cout << "This Xml is ok! ->> " << CharContainer << endl;
	}
	else
	{
		cout << "Failed to load file \n" << CharContainer << endl;
		cout <<	"Check your xml structure. \n" << endl;
		DebugStop();
	}	
	TiXmlHandle TimesHandle( &Times );
	
	Container = TimesHandle.FirstChild( "Times" ).FirstChild( "Time" ).ToElement();
	for( ; Container; Container=Container->NextSiblingElement())
	{
		CharContainer = Container->Attribute( "Value" );
		PrintStep[itime] =  atof(CharContainer);
		itime++;
	}
	
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{		
		materialist[imat]->SetTimeStep(1.0,Theta);
		materialist[imat]->SetTimeValue(0.0);			
	}	
	
	TPZAnalysis PoroelasticAnalysisInitial(ComputationalMeshPoroelasticityInitial);	
	TPZAnalysis PoroelasticAnalysis(ComputationalMeshPoroelasticityReCurrent);
	
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "FileName" ).ToElement();		
	CharContainer = Container->Attribute("InitialProblemName");
	
	// Calculation of initial conditions 
	bool Initial = true;
	TPZFMatrix<STATE> InitialSolution = PoroelasticAnalysisInitial.Solution();
	std::string output(CharContainer);
	std::stringstream outputfiletemp;
	outputfiletemp << output << ".vtk";
	
	returnvalue = system("clear");
	cout << "Calculation of initial conditions !" << endl;
	TPZVec <REAL> PrintInital(1,0);
	PrintInital[0]=0.0;
	SolveSistTransient(Initial,Events, docHandle, NumOfThreads,PrintInital,output,0.0,0.0,InitialSolution, materialist, PoroelasticAnalysisInitial, ComputationalMeshVector,ComputationalMeshPoroelasticityInitial,GeometryInfoini);

	
	// This code identify singular blocks
	TPZStepSolver<STATE> step;		// Create Solver object
	step.SetDirect(ELDLt);			//	Symmetric case
	PoroelasticAnalysis.SetSolver(step); //	Set solver	
	TPZStepSolver<STATE> & temp = dynamic_cast<TPZStepSolver<STATE> &> (PoroelasticAnalysis.Solver());
	std::list <int64_t> & zeropivot = temp.Singular(); 
	if (zeropivot.size()) 
	{
		int eq = * zeropivot.begin();
		PoroelasticAnalysis.Rhs().Zero();
		PoroelasticAnalysis.Rhs()(eq,0) = -10000.0;
		PoroelasticAnalysis.Solve();
		TPZFMatrix<STATE> TempSolution = PoroelasticAnalysis.Solution();
		std::string output;
		output = "DumpFolder/SingularNodes";
		std::stringstream outputfiletemp;
		outputfiletemp << output << ".vtk";
		std::string plotfile = outputfiletemp.str();
		PostProcessPoroeasticity(docHandle,ComputationalMeshVector,ComputationalMeshPoroelasticityReCurrent,PoroelasticAnalysis,plotfile,2);
			
		DebugStop();
	}
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{		
		materialist[imat]->SetTimeStep(Delta,Theta);			
	}	
	
//	PoroelasticAnalysis.LoadSolution(PoroelasticAnalysisInitial.Solution());
	Container = docHandle.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "FileName" ).ToElement();		
	CharContainer = Container->Attribute("ProblemName");
	std::string Outputre(CharContainer);
	std::stringstream outputfilere;
	outputfilere << Outputre << ".vtk";
	std::string plotfile = outputfilere.str();
	Initial	= false;
	SolveSistTransient(Initial,Events, docHandle, NumOfThreads,PrintStep,Outputre,Delta,MaxTime,InitialSolution, materialist, PoroelasticAnalysis, ComputationalMeshVector,ComputationalMeshPoroelasticityReCurrent,GeometryInfo);
	
	cout << "Check:: Calculation finished successfully" << endl;	
	return EXIT_SUCCESS;
}


TPZCompMesh * ComputationalElasticityMesh(TiXmlHandle ControlDoc, TPZReadGIDGrid GeometryInfo, TPZGeoMesh * gmesh,int pOrder,TPZVec < TPZFMatrix<REAL> > Events, int TStep)
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
	CharContainer = Container->Attribute("Use");	
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
		BodyForceX		= Gravitationalconstant * Gravity_Xdirection * GeometryInfo.fMaterialDataVec[imat].fProperties[0];		
		BodyForceY		= Gravitationalconstant * Gravity_Ydirection * GeometryInfo.fMaterialDataVec[imat].fProperties[0];
		Lambda			= GeometryInfo.fMaterialDataVec[imat].fProperties[5];
		G				= GeometryInfo.fMaterialDataVec[imat].fProperties[6];		
		RockDensity		= GeometryInfo.fMaterialDataVec[imat].fProperties[0]; 	
		
		Eyoung = (G*(3*Lambda+2*G))/(Lambda+G);
		PoissonRatio = (Lambda)/(2*(Lambda+G));
		MaterialElastic = new TPZElasticityMaterial(GeometryInfo.fMaterialDataVec[imat].fMatID, Eyoung, PoissonRatio, BodyForceX, BodyForceY, planestress); 
		TPZMaterial * Material(MaterialElastic);
		cmesh->InsertMaterialObject(Material);
		
		// Inserting boundary conditions
		for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
		{
			if(GeometryInfo.fMaterialDataVec[imat].fMatID == int(GeometryInfo.fBCMaterialDataVec[ibc].fProperties[0])) 
			{	
				TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);				
				TPZMaterial * BCond = MaterialElastic->CreateBC(Material,GeometryInfo.fBCMaterialDataVec[ibc].fMatID,dirichlet, val1, val2);
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
					if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].fMatID) {
						if (TStep == int(Events[itime](iwell,0))) 
						{
							TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
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
				
				if(ngroup == GeometryInfo.fMaterialDataVec[imat].fProperties[9]) 
				{	
					iset.insert(GeometryInfo.fMaterialDataVec[imat].fMatID);
					// Inseting mass line sources (Wells)
					int NumberOfEvents= Events.size();
					if (NumberOfEvents != 0) {
						int WellsPerEvent = Events[0].Rows();
						for (int itime = 0; itime < NumberOfEvents; itime++) 
						{
							for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
							{
								if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].fMatID) {
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
				if(ngroup == GeometryInfo.fBCMaterialDataVec[ibc].fProperties[5]) 
				{	
					iset.insert(GeometryInfo.fBCMaterialDataVec[ibc].fMatID);
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
		
		if (GeometryInfo.fBCMaterialDataVec[ibc].fProperties[0] == InterfacesCutoff) 
		{
			if (InterfaceBehavior) 
			{
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[2];
				
				
				ElasticMatInterface2D *InterfaceMat = new ElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].fMatID,Eyoung, PoissonRatio, BodyForceX, BodyForceY, planestress); 
				InterfaceMat->SetPenalty(KNU,KTU);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				cmesh->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				cmesh->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int64_t ngel = cmesh->Reference()->NElements();
				
				for(int64_t igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].fMatID) {
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
						int64_t gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}				
				
			}
			else 
			{
				
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[2];
				
				ElasticMatInterface2D *InterfaceMat = new ElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].fMatID,Eyoung, PoissonRatio, BodyForceX, BodyForceY, planestress); 
				InterfaceMat->SetPenalty(KNU,KTU);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				cmesh->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				cmesh->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int64_t ngel = cmesh->Reference()->NElements();
				
				for(int64_t igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].fMatID) {
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
						int64_t gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}
				
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
	CharContainer = Container->Attribute("Use");	
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
		BodyForceX		= Gravitationalconstant * Gravity_Xdirection * GeometryInfo.fMaterialDataVec[imat].fProperties[1];		
		BodyForceY		= Gravitationalconstant * Gravity_Ydirection * GeometryInfo.fMaterialDataVec[imat].fProperties[1];
		Permeability	= GeometryInfo.fMaterialDataVec[imat].fProperties[3];
		Viscosity		= GeometryInfo.fMaterialDataVec[imat].fProperties[4];		
		FluidDensity	= GeometryInfo.fMaterialDataVec[imat].fProperties[1]; 
		
		diff = (Permeability)/(Viscosity);
		
		MaterialDiffusion = new TPZMatPoisson3d(GeometryInfo.fMaterialDataVec[imat].fMatID, dim);
		MaterialDiffusion->SetParameters(diff, conv, convdir);
		MaterialDiffusion->SetInternalFlux(flux);
		MaterialDiffusion->NStateVariables();		
		TPZMaterial * Material(MaterialDiffusion);
		cmesh->InsertMaterialObject(Material);
		
		// Inserting boundary conditions
		for (int ibc = 0; ibc < GeometryInfo.BCNumber ; ibc++) 
		{
			if(GeometryInfo.fMaterialDataVec[imat].fMatID == int(GeometryInfo.fBCMaterialDataVec[ibc].fProperties[0])) 
			{	
				TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);				
				TPZMaterial * BCond = MaterialDiffusion->CreateBC(Material,GeometryInfo.fBCMaterialDataVec[ibc].fMatID,dirichlet, val1, val2);
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
					if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].fMatID) {
						if (TStep == int(Events[itime](iwell,0))) 
						{
							TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
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
				
				if(ngroup == GeometryInfo.fMaterialDataVec[imat].fProperties[9]) 
				{	
					iset.insert(GeometryInfo.fMaterialDataVec[imat].fMatID);
					// Inseting mass line sources (Wells)
					int NumberOfEvents= Events.size();
					if (NumberOfEvents != 0) {
						int WellsPerEvent = Events[0].Rows();
						for (int itime = 0; itime < NumberOfEvents; itime++) 
						{
							for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
							{
								if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].fMatID) {
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
				if(ngroup == int(GeometryInfo.fBCMaterialDataVec[ibc].fProperties[5])) 
				{	
					
					if (0 < GeometryInfo.fBCMaterialDataVec[ibc].fProperties[6]) 
					{
						iset.insert(GeometryInfo.fBCMaterialDataVec[ibc].fMatID);
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
			
			if (0 < GeometryInfo.fMaterialDataVec[imat].fProperties[10]) 
			{
				iset.insert(GeometryInfo.fMaterialDataVec[imat].fMatID);
				// Inseting mass line sources (Wells)
				int NumberOfEvents= Events.size();
				if (NumberOfEvents != 0) {
					int WellsPerEvent = Events[0].Rows();
					for (int itime = 0; itime < NumberOfEvents; itime++) 
					{
						for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
						{
							if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].fMatID) {
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
			
			if (0 < GeometryInfo.fBCMaterialDataVec[ibc].fProperties[6]) 
			{
				iset.insert(GeometryInfo.fBCMaterialDataVec[ibc].fMatID);				
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
		
		if (GeometryInfo.fBCMaterialDataVec[ibc].fProperties[0] == InterfacesCutoff) 
		{
			if (InterfaceBehavior) 
			{
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[2];
				
				
				ElasticMatInterface2D *InterfaceMat = new ElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].fMatID,Eyoung, PoissonRatio, BodyForceX, BodyForceY, 0); 
				InterfaceMat->SetPenalty(KNU,KTU);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				cmesh->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				cmesh->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int64_t ngel = cmesh->Reference()->NElements();
				
				for(int64_t igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].fMatID) {
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
						int64_t gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}				
				
			}
			else 
			{
				
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[2];
				
				ElasticMatInterface2D *InterfaceMat = new ElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].fMatID,Eyoung, PoissonRatio, BodyForceX, BodyForceY, 0); 
				InterfaceMat->SetPenalty(KNU,KTU);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				cmesh->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				cmesh->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int64_t ngel = cmesh->Reference()->NElements();
				
				for(int64_t igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].fMatID) {
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
						int64_t gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}
				
			}
			
			
		}
		
	}	
	
	
	return cmesh;
	
}

TPZCompMesh * ComputationalPoroelasticityMesh(TiXmlHandle ControlDoc, TPZReadGIDGrid &GeometryInfo, 
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
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "AnalyticFunction" ).FirstChild( "Function" ).ToElement();		
	CharContainer = Container->Attribute("FunctionName");
	
	// Analitical functions for validation
	switch (atoi(CharContainer)) {
		case 1:
		{
				TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolutionfiniteColumn1D, 5);
		}
			break;
		case 2:
		{
			TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolutionSemiInfiniteColumn1D, 5);
		}
			break;
		case 3:
		{
			TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolution2DLineSource, 5);
		}
			break;			
		default:
		{
			TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolutionfiniteColumn1D, 5);
		}			
			break;
	}
	
	if (atoi(CharContainer) == 1) 
	{

	}
	
	if (atoi(CharContainer) == 2) 
	{
		TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolutionSemiInfiniteColumn1D, 5);
	}
	
	if (atoi(CharContainer) == 3) 
	{
		TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolution2DLineSource, 5);
	}
	
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
	CharContainer = Container->Attribute("Use");	
	int InterfaceBehavior = atoi(CharContainer);
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "DimensionlessCalculation" ).FirstChild( "Calculation" ).ToElement();		
	CharContainer = Container->Attribute("Use");	
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
	REAL SFluid			=	0.0;	
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{	
		TPZPoroElastic2d * MaterialPoroElastic;
		LambdaU			= GeometryInfo.fMaterialDataVec[imat].fProperties[7];
		Lambda			= GeometryInfo.fMaterialDataVec[imat].fProperties[5];		
		G				= GeometryInfo.fMaterialDataVec[imat].fProperties[6];
		alpha			= GeometryInfo.fMaterialDataVec[imat].fProperties[8];
		RockDensity		= GeometryInfo.fMaterialDataVec[imat].fProperties[0];		
		Permeability	= GeometryInfo.fMaterialDataVec[imat].fProperties[3];
		FluidViscosity	= GeometryInfo.fMaterialDataVec[imat].fProperties[4];
		FluidDensity	= GeometryInfo.fMaterialDataVec[imat].fProperties[1];		
		RockPorosity	= GeometryInfo.fMaterialDataVec[imat].fProperties[2];		
		diff = (Permeability)/(FluidViscosity);
		MixtureDensity = (1 - RockPorosity) * RockDensity + RockPorosity * FluidDensity;
		SFluid = LambdaU;
		// Using Gravity field
		BodyForceX		= Gravitationalconstant * Gravity_Xdirection * MixtureDensity;	
		BodyForceY		= Gravitationalconstant * Gravity_Ydirection * MixtureDensity;
		
		
		
		if (Lambda != LambdaU) 
		{
			if (alpha != 0.0) {
				S = ((pow(alpha,2))/((LambdaU-Lambda)))*((LambdaU+2.0*G)/(Lambda+2.0*G));
				Se = (pow(alpha,2))/((LambdaU-Lambda));
			}
			else 
			{
				S = (1)/(3*Lambda+2.0*G);
				Se = RockPorosity*(SFluid + S);
			}
			
		}
		else
		{
			
			if (alpha != 0.0) {
				S = (pow(alpha,2)/(Lambda+2.0*G));
				Se = 0;
			}
			else 
			{
				S = (1.0)/(3*Lambda+2.0*G);
				Se = 0;
			}			
			

			
		}		
		
		
		if (Dimensionless) 
		{
			if (alpha != 0.0) {
				Se				= Se/S;
			}
			else 
			{
				if (Se != 0.0) {
					Se				= 1.0;
				}
			}			
			Lambda			= Lambda*S;
			LambdaU			= LambdaU*S;			
			G				= G*S;
			Permeability	= 1.0;
			FluidViscosity	= 1.0;
			
		}	

		MaterialPoroElastic = new TPZPoroElastic2d (GeometryInfo.fMaterialDataVec[imat].fMatID, dim);
		MaterialPoroElastic->SetParameters(Lambda,G,LambdaU,BodyForceX,BodyForceY);
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
			if(GeometryInfo.fMaterialDataVec[imat].fMatID == int(GeometryInfo.fBCMaterialDataVec[ibc].fProperties[0])) 
			{	
				TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
				val2(0,0)=GeometryInfo.fBCMaterialDataVec[ibc].fProperties[2];
				val2(1,0)=GeometryInfo.fBCMaterialDataVec[ibc].fProperties[3];			
				val2(2,0)=GeometryInfo.fBCMaterialDataVec[ibc].fProperties[4];
				TPZMaterial * BCond = MaterialPoroElastic->CreateBC(Material,GeometryInfo.fBCMaterialDataVec[ibc].fMatID,GeometryInfo.fBCMaterialDataVec[ibc].fProperties[1], val1, val2);
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
					if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].fMatID) {
						if (TStep == int(Events[itime](iwell,0))) 
						{
							TPZFMatrix<STATE> val1(3,2,Events[itime](iwell,4)), val2(3,1,0.);
							val2(0,0)=Events[itime](iwell,5);
							val2(1,0)=Events[itime](iwell,6);			
							val2(2,0)=Events[itime](iwell,7);
							TPZMaterial * BCLINE = MaterialPoroElastic->CreateBC(Material, int(Events[itime](iwell,1)),int(Events[itime](iwell,3)), val1, val2);
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
				
				if(ngroup == GeometryInfo.fMaterialDataVec[imat].fProperties[9]) 
				{	
					iset.insert(GeometryInfo.fMaterialDataVec[imat].fMatID);
					// Inseting mass line sources (Wells)
					int NumberOfEvents= Events.size();
					if (NumberOfEvents != 0) {
						int WellsPerEvent = Events[0].Rows();
						for (int itime = 0; itime < NumberOfEvents; itime++) 
						{
							for (int iwell = 0; iwell < WellsPerEvent; iwell++) 
							{
								if (int(Events[itime](iwell,2)) == GeometryInfo.fMaterialDataVec[imat].fMatID) {
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
				if(ngroup == GeometryInfo.fBCMaterialDataVec[ibc].fProperties[5]) 
				{	
					iset.insert(GeometryInfo.fBCMaterialDataVec[ibc].fMatID);
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
		
		if (GeometryInfo.fBCMaterialDataVec[ibc].fProperties[0] == InterfacesCutoff) 
		{
			if (InterfaceBehavior) 
			{
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[2];
				REAL KNP = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[3];
				REAL KTP;
                KTP = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[4];
				REAL friction = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[5];				
				bool docontribute = true;
				
				PoroElasticMatInterface2D *InterfaceMat = new PoroElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].fMatID,1,docontribute,friction); 
				InterfaceMat->SetPenalty(KNU,KTU,KNP,KNP);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				Multiphysics->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				Multiphysics->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int64_t ngel = Multiphysics->Reference()->NElements();
				
				for(int64_t igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].fMatID) {
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
						int64_t gelindex;
						
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZMultiphysicsInterfaceElement(* Multiphysics, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}			
				
			}
			else 
			{
				
				//		Inserting interface element					
				REAL KNU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[1];
				REAL KTU = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[2];
				REAL KNP = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[3];
				REAL KTP;
                KTP = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[4];
				REAL friction = GeometryInfo.fBCMaterialDataVec[ibc].fProperties[5];	
				bool docontribute = false;				
				
				PoroElasticMatInterface2D *InterfaceMat = new PoroElasticMatInterface2D(GeometryInfo.fBCMaterialDataVec[ibc].fMatID,1,docontribute,friction); 
				InterfaceMat->SetPenalty(KNU,KTU,KNP,KNP);		
				TPZMaterial * MaterialInterface(InterfaceMat);
				Multiphysics->InsertMaterialObject(MaterialInterface);			
				
				// Load references for search interface locations
				Multiphysics->LoadReferences();
				//				ofstream Argument("CMeshBefore.txt");		
				//				cmesh->Print(Argument);
				int64_t ngel = Multiphysics->Reference()->NElements();
				
				for(int64_t igel = 0; igel < ngel ; igel++ )
				{
					TPZGeoEl * Gel = gmesh->ElementVec()[igel];
					if (!Gel) {
						continue;
					}
					int MatId = Gel->MaterialId();
					if (MatId == GeometryInfo.fBCMaterialDataVec[ibc].fMatID) {
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
						int64_t gelindex;
						// TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);				
						new TPZMultiphysicsInterfaceElement(* Multiphysics, Gel, gelindex, neigh[0], neigh[1] );		
						
					}
					
				}
				
			}
			
			
		}
		
	}		
	
	Multiphysics->AdjustBoundaryElements();
	Multiphysics->CleanUpUnconnectedNodes();
	
	// Creating multiphysics elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, Multiphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,Multiphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, Multiphysics);
	
	ofstream arg("DumpFolder/ComputationalMultiphysicMesh.txt");
	Multiphysics->Print(arg);
	
	return Multiphysics;
}

#define VTK
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{			
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
	an.SetSolver(step);
	an.Run();
}

void PostProcessDiffusion(TPZAnalysis &an, std::string plotfile)
{
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
	vecnames[0]= "MinusKGradU";
	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
}

void PostProcessElasticity(TPZAnalysis &an, std::string plotfile)
{
	TPZManVector<std::string,10> scalnames(2), vecnames(1);
	scalnames[0] = "SigmaX";
	scalnames[1] = "SigmaY";	
	vecnames[0]= "displacement";
	
	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
}

void PostProcessPoroeasticity(TiXmlHandle ControlDoc, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile, int dim)
{
	TPZBuildMultiphysicsMesh * Objectdumy = NULL;
	Objectdumy->TransferFromMultiPhysics(meshvec, mphysics);
	
	// TiXmlElement dummy object
	TiXmlElement* Container;
	const char * CharContainer;	
	
	if (dim == 2) 
	{
		Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "FileName" ).ToElement();		
		CharContainer = Container->Attribute("NDiv");
		int div = atoi(CharContainer);
		CharContainer = Container->Attribute("NVectorials");
		int NVectorials = atoi(CharContainer);
		CharContainer = Container->Attribute("NScalars");
		int NScalars = atoi(CharContainer);	
		
		TPZManVector<std::string,10> scalnames(NScalars), vecnames(NVectorials);
		
		int iscalar = 0;
		Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "Scalars" ).FirstChild( "Var" ).ToElement();
		for( ; Container; Container=Container->NextSiblingElement())
		{
			CharContainer = Container->Attribute("Name");
			scalnames[iscalar] = CharContainer;
			iscalar++;
		}
		
		int ivectorial = 0;
		Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "OutputControls" ).FirstChild( "Vectorials" ).FirstChild( "Var" ).ToElement();
		for( ; Container; Container=Container->NextSiblingElement())
		{
			CharContainer = Container->Attribute("Name");
			vecnames[ivectorial] = CharContainer;
			ivectorial++;
		}
		
		const int dimension = dim;
		an.DefineGraphMesh(dimension,scalnames,vecnames,plotfile);
		an.PostProcess(div,dim);		
	}
	else 
	{
		Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "OutputControlsInterface" ).FirstChild( "FileName" ).ToElement();		
		CharContainer = Container->Attribute("NDiv");
		int div = atoi(CharContainer);
		CharContainer = Container->Attribute("NVectorials");
		int NVectorials = atoi(CharContainer);
		CharContainer = Container->Attribute("NScalars");
		int NScalars = atoi(CharContainer);	
		
		TPZManVector<std::string,10> scalnames(NScalars), vecnames(NVectorials);
		
		int iscalar = 0;
		Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "OutputControlsInterface" ).FirstChild( "Scalars" ).FirstChild( "Var" ).ToElement();
		for( ; Container; Container=Container->NextSiblingElement())
		{
			CharContainer = Container->Attribute("Name");
			scalnames[iscalar] = CharContainer;
			iscalar++;
		}
		
		int ivectorial = 0;
		Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "OutputControlsInterface" ).FirstChild( "Vectorials" ).FirstChild( "Var" ).ToElement();
		for( ; Container; Container=Container->NextSiblingElement())
		{
			CharContainer = Container->Attribute("Name");
			vecnames[ivectorial] = CharContainer;
			ivectorial++;
		}
		
		const int dimension = dim;
		an.DefineGraphMesh(dimension,scalnames,vecnames,plotfile);
		an.PostProcess(div,dim);		
	}
	
}


void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int64_t n = gMesh->NElements();
		for ( int64_t i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
		}//for i
	}//ref
}

void UniformRefinement(TPZGeoMesh * gMesh, int nh, int MatId)
{
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int64_t n = gMesh->NElements();
		for ( int64_t i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2 || gel->Dimension() == 1){
				if (gel->MaterialId()== MatId){
					gel->Divide (filhos);
				}
			}
		}//for i
	}//ref
}

void RefinElemComp(TPZCompMesh  *cMesh, int indexEl)
{
	
	TPZVec<int64_t > subindex; 
	int64_t nel = cMesh->ElementVec().NElements(); 
	for(int64_t el=0; el < nel; el++){
		TPZCompEl * compEl = cMesh->ElementVec()[el];
		if(!compEl) continue;
		int64_t ind = compEl->Index();
		if(ind==indexEl){
			compEl->Divide(indexEl, subindex, 1);
		}
	}	
}

void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
{
	
	TPZVec<int64_t > subindex;
	for (int64_t iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		int64_t nel = elvec.NElements(); 
		for(int64_t el=0; el < nel; el++){
			TPZCompEl * compEl = elvec[el];
			if(!compEl) continue;
			int64_t ind = compEl->Index();
			compEl->Divide(ind, subindex, 0);
		}
	}
}

TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZReadGIDGrid GeometryInfo,TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics)
{
	
	
	TPZSpStructMatrix matsp(mphysics);
//	TPZSkylineStructMatrix matsp(mphysics);
	std::set< int > materialid;
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{
		mymaterial[imat]->SetLastState();
		materialid.insert(GeometryInfo.fMaterialDataVec[imat].fMatID );
	}	
	
	matsp.SetMaterialIds(materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<STATE> Un;
	TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
	
	return matK2;
	
}


void StiffMatrixLoadVec(TPZReadGIDGrid GeometryInfo,int Nthreads, TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<STATE> &fvec)
{
	
	for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
	{	
		mymaterial[imat]->SetCurrentState();
	}		
	
	//TPZFStructMatrix matsk(mphysics);
	TPZSkylineStructMatrix matsk(mphysics);
	an.SetStructuralMatrix(matsk); 
	an.StructMatrix()->SetNumThreads(Nthreads);
	TPZStepSolver<STATE> step; 
	step.SetDirect(ELDLt); 
	an.SetSolver(step); 
	an.Assemble(); 
	matK1 = an.StructMatrix();
	fvec = an.Rhs();
	
}



void SolveSistTransient(bool IsInitial, TPZVec < TPZFMatrix<REAL> > Events,TiXmlHandle ControlDoc,int Nthreads, TPZVec <REAL> &PrintStep, std::string FileName, REAL deltaT,REAL maxTime, TPZFMatrix<STATE> &InitialSolution, TPZVec <TPZPoroElastic2d  * > &mymaterial ,
						TPZAnalysis &an, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics , TPZReadGIDGrid & GeometryInfo)
{
	
	
	// Calculation of Mass Matrix
	cout << "Calculate Mass Matrix " << endl;    
	TPZAutoPointer <TPZMatrix<STATE> > PoroelasticMassMatrix = MassMatrix(GeometryInfo,mymaterial, mphysics);
	cout << "Calculate Mass Matrix was done! " << endl;	
	
	// Calculation of Stiffness Matrix and Load Vector
	TPZFMatrix<REAL> PoroelasticStiffnessMatrix, PoroelasticStiffnessMatrixInverse;	
	TPZFMatrix<STATE> PoroelasticLoadVector;
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
	
	Container = ControlDoc.FirstChild( "ProblemData" ).FirstChild( "FaultInterface" ).FirstChild( "Contact" ).ToElement();			
	CharContainer = Container->Attribute("UseNormal");		
	int NormalCalculations = atoi(CharContainer);	
	
	int nrows;
	nrows = PoroelasticMassMatrix->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> RhsTemporal(nrows,1,0.0);
	TPZFMatrix<STATE> LastSolution = InitialSolution;
	TPZVec <TPZPoroElastic2d *>  Dummymaterialist(GeometryInfo.MatNumber,0);	
	
	
	REAL	TimeValue	= 0.0;
	int		cent		= 0;
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
		
		if (IsInitial) 
		{
			TotalRhs=PoroelasticLoadVector;
			an.Rhs() = TotalRhs;
			an.Solve();
			InitialSolution = an.Solution();
			std::stringstream outputfiletemp;
			outputfiletemp << FileName << ".vtk";
			std::string plotfile = outputfiletemp.str();
			PostProcessPoroeasticity(ControlDoc,meshvec,mphysics,an,plotfile,2);
			if (NormalCalculations) 
			{
				std::stringstream filetemp;
				filetemp << "Output/FaultGeometry" << ".vtk";
				std::string FaultFile = filetemp.str();					
				PostProcessPoroeasticity(ControlDoc,meshvec,mphysics,an,FaultFile,1);				
			}			
			break;
		}
		
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
			PostProcessPoroeasticity(ControlDoc,meshvec,mphysics,an,plotfile,2);
			if (NormalCalculations) 
			{
				std::stringstream filetemp;
				filetemp << "Output/Fault" << ".vtk";
				std::string FaultFile = filetemp.str();								
				PostProcessPoroeasticity(ControlDoc,meshvec,mphysics,an,FaultFile,1);
				cout << "Created file with normals at time step = " << cent << " (time = " << TimeValue << ")" << endl;					
			}
			cout << "Created Vtk file at time step = " << cent << " (time = " << TimeValue << ")" << endl;			
			control++;
			
		}
		
		if (NEvents != 0) 
		{
			for (int ievent = 0; ievent < NEvents; ievent++) 
			{
				
				if (cent == Events[ievent](0,0) ) 
				{
					cout << "Modifiying Poroelastic load Vector ...." << endl;
					TPZCompMesh * DumpMphysics = ComputationalPoroelasticityMesh(ControlDoc, GeometryInfo, meshvec[0]->Reference(), meshvec, Dummymaterialist, Events, cent);
					TPZAnalysis DumpAnalysis(DumpMphysics);
					for(int imat = 0; imat < GeometryInfo.MatNumber; imat++)
					{
						Dummymaterialist[imat]->SetTimeStep(deltaT,Theta);
						Dummymaterialist[imat]->SetTimeValue(TimeValue);
						Dummymaterialist[imat]->SetCurrentState();
					}
					
					TPZSkylineStructMatrix matsk(DumpMphysics);
					DumpAnalysis.SetStructuralMatrix(matsk);
					DumpAnalysis.StructMatrix()->SetNumThreads(Nthreads);
					TPZStepSolver<STATE> step;
					step.SetDirect(ELDLt); 					
					DumpAnalysis.SetSolver(step);
					DumpAnalysis.AssembleResidual();
					PoroelasticLoadVector = DumpAnalysis.Rhs();
					cout << "Modifiying Poroelastic load Vector -> done! at time step = " << cent << " ,  or time value = " << TimeValue << endl;	
				}
			}
		}		
		
		
		cent++;
		TimeValue = cent*deltaT;
	}
	int returnvalue;
	returnvalue = system("clear");
}

void *InitalPressureCalculations(REAL &PressureValue,TPZVec <REAL> &PointCoordinates)
{
	PressureValue = 0.0;
	return NULL;
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





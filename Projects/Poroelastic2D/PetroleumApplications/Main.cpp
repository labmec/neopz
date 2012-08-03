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
#include "MeshGeneration.h"
#include "AnalyticalFunctions.h"
#include "ElasticMatInterface2D.h"
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


#include <cmath>
#include <set>


//	//	Use this tools for Exponential Integral Fuction (Implementation defined at end of this main file) 
//	#include <math.h>		// required for fabsl(), expl() and logl()        
//	#include <float.h>		// required for LDBL_EPSILON, DBL_MAX
//	//	Internally Defined Routines
//	double      Exponential_Integral_Ei( double x );
//	long double xExponential_Integral_Ei( long double x );
//	static long double Continued_Fraction_Ei( long double x );
//	static long double Power_Series_Ei( long double x );
//	static long double Argument_Addition_Series_Ei( long double x);
//	//	Internally Defined Constants
//	static const long double epsilon = 10.0 * LDBL_EPSILON;
//	//	End Use this tools for Exponential Integral Fuction (Implementation defined at end of this main file)

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
const int TotalMats = 3; 
const int matId = 1;// 1 Rock Formation
const int matIdl2 = 2;// 2 Rock Formation
const int matIdl3 = 3;// 3 Rock Formation
const int matIdl4 = 4;// 4 Rock Formation
// const int matIdl5 = n;// N Rock formation	

// Boundary Conditions Materials
// 2D Rectangular Domain
const int bcBottom = -1;
const int bcRight = -2;
const int bcTop = -3;
const int bcLeft = -4;
const int yfixedPoints = -6;
const int xfixedPoints = -7;
const int pointsource = -5;

// 2D Cylindrical Domain
const int arc1 = -1;
const int arc2 = -2;
const int arc3 = -3;
const int arc4 = -4;
const int WellPoint = -5;

// Case 1 2D Geological Fault Reactivation stess state = gravity field
// Case 2 2D Geological Fault Reactivation, initialization with fault consistent stress state Geological Graven Structure

// Rock Blocks bondaries
// Right Block 
//const int RBright = -12;
//const int RBleft = -22;
//const int RBBot = -32;
//const int RBTop = -42;
//
// Left Block 
//const int LBright = -12;
//const int LBleft = -22;
//const int LBBot = -32;
//const int LBTop = -42;
//
// Graven Block 
//const int GBright = -12;
//const int GBleft = -22;
//const int GBBot = -32;
//const int GBTop = -42;

// Rock Blocks bondaries
// Right Block 
int RBright = -23;
int RBleft = -43;
int RBBot = -13;
int RBTop = -33;

// Left Block 
int LBright = -22;
int LBleft = -42;
int LBBot = -12;
int LBTop = -32;

// Graven Block 
int GBright = -21;
int GBleft = -41;
int GBBot = -11;
int GBTop = -31;

const int WellLine = -5;

// 2D SPE comparative Study
// not implemented

// Boundary Conditions Definitions
// This Definitions are related with the Boudary conditions convention for multphysics simulation (in Poroelastic U/P formulation problems 3 boundary conditions are included) 
// not defined
const int dirichlet = 0;
const int neumann = 1;


// Defintions of Implemented Methods

//	This Create Computational Meshes
TPZCompMesh *MalhaCompPressao(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh *MalhaCompElast(TPZGeoMesh * gmesh,int pOrder, bool Initialstress);
TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial);

//	This Solve Different analysis
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);
void SolveSistTransient(REAL delt,REAL maxtime,TPZFMatrix<REAL> matK1, TPZAutoPointer < TPZMatrix<REAL> > matK2, TPZFMatrix<REAL> fvec, TPZFMatrix<REAL> &Initialsolution, TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZAnalysis &an,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);

//	This are tools for spatial and polinomyal refinement and Postprocess of solutions 
void PosProcess(TPZAnalysis &an, std::string plotfile);
void PosProcess2(TPZAnalysis &an, std::string plotfile);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
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
	
	
	// GEOMETRICAL MESH CREATION	
	
	//	polynomial base degree approach
	int p=2;
	
	//  Used this for Directional Refinement	
	
	//	gRefDBase.InitializeRefPatterns();
	//	gRefDBase.InitializeAllUniformRefPatterns();
	//	ofstream RefinElPatterns("RefElPatterns.txt");
	//	gRefDBase.WriteRefPatternDBase(RefinElPatterns); 	
	//	gRefDBase.ReadRefPatternDBase("RefElPatterns.txt");
	
	//	Rectangular geometric mesh for graven analysis
	
	int nLayers			= 5; 
	bool InterfaceEl	= true; 
	REAL LlengthFootFault	= 1600.0;
	REAL DipFaultAngleleft		= 45.0;
	REAL DipFaultAngleright		= 45.0;
	REAL WellFaultlength		= 800.0;
	TPZVec <bool> wichProductionlayer(nLayers+1,false);
	//	wichProductionlayer[0] = true; // Horizons -> Top
	//	wichProductionlayer[1] = true; // Horizons -> Top	
	//	wichProductionlayer[2] = true; // Horizons -> Top
	//	wichProductionlayer[3] = true; // Horizons -> Top
	//	wichProductionlayer[3] = true; // Horizons -> Top
	
	
	std::string FiletoRead;
//	FiletoRead = "Anticline3wells.dump"; 
	FiletoRead = "labyrinth.dump";
//	FiletoRead = "laberinto.dump";	
//	FiletoRead = "AnticlineColoradoTet.dump"; 	
//	FiletoRead = "AnticlineTet.dump"; 
//	FiletoRead = "AnticlineQuad.dump"; 
//	FiletoRead = "AnticlineQuadfine.dump"; 
//	FiletoRead = "cube.dump"; 	
	
//	MeshGeneration mygenerator;
//	TPZGeoMesh * gmesh = mygenerator.GeometricGIDMesh(FiletoRead);
	
	TPZReadGIDGrid myreader;
	TPZGeoMesh * gmesh = myreader.GeometricGIDMesh(FiletoRead);
	
//	TPZGeoMesh * gmesh = mygenerator.MalhaGeoGravenobj(nLayers,LlengthFootFault,DipFaultAngleleft,DipFaultAngleright,WellFaultlength,wichProductionlayer,InterfaceEl);		
	
	// End Rectangular geometric mesh for graven analysis
	
	
	// END GEOMETRICAL MESH CREATION	
	
	{
	//	Print Geometrical Base Mesh
	ofstream arg1("BaseGeoMesh.txt");
	gmesh->Print(arg1);
	ofstream file1("BaseGeoMesh.vtk");	
	//	In this option true -> let you use shrink paraview filter
	//	shrink filter in paraview just use for "amazing" visualization
	//	PrintGMeshVTK(gmesh,file1);	
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file1, true);	
	}
		
//	// Directional Point source Refinement 
//	int ndirectdivp = 17;
//	set<int> SETmatPointRefDir2;
//	SETmatPointRefDir2.insert(WellPoint);	
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
//	
	//	Modifying Geometric Mesh
	int Href = 1;
	RefinamentoUniforme(gmesh, Href);
	
	{
	//	Print Geometrical refined Base Mesh
	ofstream arg1("RefineGeoMesh.txt");
	gmesh->Print(arg1);
	ofstream file1("RefineGeoMesh.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file1, true);	
	}
		
	//	First computational mesh
	
	// initialization problem ...
	bool Initialstress = true;
	TPZCompMesh * cmesh1 = MalhaCompElast(gmesh, p+1, Initialstress);
	//	Print First computational mesh
	ofstream arg2("cmeshElasticity.txt");
	cmesh1->Print(arg2);
	
	//	Second computational mesh
	//	TPZCompMesh * cmesh2= MalhaCompPressao(gmesh, p);
	////	Print Second computational mesh
	//	ofstream arg3("cmesh2.txt");
	//	cmesh2->Print(arg3);	
	
	//	Clear reference of the geometric mesh to cmesh1
//	gmesh->ResetReference();
//	cmesh1->LoadReferences();
//	
//	//	Using Uniform Refinement
//	RefinUniformElemComp(cmesh1,4);
//	cmesh1->AdjustBoundaryElements();
//	cmesh1->CleanUpUnconnectedNodes();
	
	//	ofstream arg4("cmesh12.txt");
	//	cmesh1->Print(arg4);
	//	ofstream arg5("gmesh2.txt");
	//	gmesh->Print(arg5);
	//	ofstream file3("malhageo1.vtk");
	//	PrintGMeshVTK(gmesh, file3);
	
	// Clear reference of the geometric mesh to cmesh2
	//	gmesh->ResetReference();
	//	cmesh2->LoadReferences();
	
	//	Using Uniform Refinement
	//	RefinUniformElemComp(cmesh2,1);
	//	cmesh2->AdjustBoundaryElements();
	//	cmesh2->CleanUpUnconnectedNodes();
	//	cmesh2->ExpandSolution();
	
	//	ofstream arg6("cmesh22.txt");
	//	cmesh2->Print(arg6);
	//	ofstream arg7("gmesh3.txt");
	//	gmesh->Print(arg7);
	//	ofstream file5("malhageo2.vtk");
	//	PrintGMeshVTK(gmesh, file5);
	
	//	Solving First Computational problem
	TPZAnalysis an1(cmesh1);
	SolveSist(an1, cmesh1);
	std::string plotfile("SolutionElasticity.vtk");		
	PosProcess2(an1, plotfile);	
	
	// This code identify singular blocks
	TPZStepSolver<REAL> & temp = dynamic_cast<TPZStepSolver<REAL> &> (an1.Solver());
	std::list <int> & zeropivot = temp.Singular(); 
	if (zeropivot.size()) 
	{
		int eq = * zeropivot.begin();
		an1.Rhs().Zero();
		an1.Rhs()(eq,0) = -10000.0;
		an1.Solve();
		TPZFMatrix<REAL> TempSolution = an1.Solution();
		
#ifdef LOG4CXX
		// Print the temporal solution
		if(logdata->isDebugEnabled())
		{
			std::stringstream sout;
			TempSolution.Print("Singularnodes = ", sout,EMathematicaInput);
			LOGPZ_DEBUG(logdata,sout.str())
		}
#endif	
		std::string output;
		output = "Singularnodes";
		std::stringstream outputfiletemp;
		outputfiletemp << output << ".vtk";
		std::string plotfile = outputfiletemp.str();
		PosProcess2(an1,plotfile);
		
	}		
	
	
	
	//	Solving Second Computational problem
	//	TPZAnalysis an2(cmesh2);
	//	SolveSist(an2, cmesh2);
	//	std::string plotfile2("saidaSolution_cmesh2.vtk");
	//	PosProcess(an2, plotfile2);
	
	
	//	Solving Multiphysic Computational problem
	//	TPZVec<TPZCompMesh *> meshvec(2);
	//	meshvec[0] = cmesh1;
	//	meshvec[1] = cmesh2;
	//	TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,1);
	//	 
	//	ofstream file6("mphysics.vtk");
	//	PrintGMeshVTK(gmesh, file6);
	//	 
	//	TPZAnalysis an(mphysics);
	//	SolveSist(an, mphysics);
	//	std::string plotfile3("saidaMultphysics.vtk");
	//	PosProcessMultphysics(meshvec,mphysics, an, plotfile3);
	//
	
	// This part is for Poroelastic analysis	
	
	//	Set initial conditions for pressure 
	//	int nrs = an2.Solution().Rows();
	//	TPZFMatrix<REAL> solucao1(nrs,1,0.0);
	//	cmesh2->Solution() = solucao1;
	//	
	//	TPZVec<TPZCompMesh *> meshvec(2);
	//	meshvec[0] = cmesh1;
	//	meshvec[1] = cmesh2;
	//	TPZVec <TPZPoroElastic2d *>  materialist(TotalMats,0) ;
	//	TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,materialist);
	//	
	//	
	//	// time control
	//	REAL hour = 3600;
	//	REAL day = 86400;
	//	REAL month = 30*day;
	//	REAL year = 365*day;
	//	REAL delta = 81.999999999*day;
	//	REAL MaxTime = 164*day;
	//	
	//	// Improve this
	//	for(int imat = 0; imat < TotalMats; imat++)
	//	{	
	//	materialist[imat]->SetTimeStep(delta);		
	//	materialist[imat]->SetTimeValue(0.0);
	//	//Criando matriz K2
	//	materialist[imat]->SetLastState();
	//	}
	//	
	//	TPZAnalysis an(mphysics);
	//	
	//	//Setting initial coditions 
	//	TPZFMatrix<REAL> Initialsolution = an.Solution();
	//	std::string output;
	//	output = "TransientSolution";
	//	std::stringstream outputfiletemp;
	//	outputfiletemp << output << ".vtk";
	//	std::string plotfile = outputfiletemp.str();
	//	//Print Initial conditions
	//	PosProcessMultphysics(meshvec,mphysics,an,plotfile);	
	//	
	//	TPZSpStructMatrix matsp(mphysics);
	//	//TPZSkylineStructMatrix matsp(mphysics);	
	//	std::set< int > materialid;
	//	materialid.clear();	
	//	int matid = 1;
	//	materialid.insert(imat);
	//	matsp.SetMaterialIds (materialid);
	//	// Improve this
	//	for(int imat = 0; imat < TotalMats; imat++)
	//	{	
	//		materialid.insert(imat+1);	
	//	}
	//	// Inserting Material ids for matrix K2 creation
	//	matsp.SetMaterialIds(materialid);		
	//	
	//	TPZAutoPointer<TPZGuiInterface> guiInterface;
	//	TPZFMatrix<REAL> Un;
	//	TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
	//	
	//#ifdef LOG4CXX
	//	if(logdata->isDebugEnabled())
	//	{
	//		std::stringstream sout;
	//		matK2->Print("K2 = ", sout,EMathematicaInput);
	//		LOGPZ_DEBUG(logdata,sout.str())
	//	}
	//#endif
	//	
	//	//Criando matriz K1
	//	// Improve this
	//	for(int imat = 0; imat < TotalMats; imat++)
	//	{	
	//	materialist[imat]->SetCurrentState();
	//	}	
	//	TPZSkylineStructMatrix matsk(mphysics);
	//	TPZFMatrix<REAL> matK1;	
	//	TPZFMatrix<REAL> fvec; //vetor de carga
	//	an.SetStructuralMatrix(matsk);//	Set the structtural sky line matrix to our analysis
	//	
	//	
	//	TPZStepSolver<REAL> step; //Create Solver object
	//	step.SetDirect(ELDLt); //	Symmetric case
	//	//step.SetDirect(ELU);
	//	an.SetSolver(step); //	Set solver
	//	an.Run(); //	Excecute an analysis to obtain the Rhs vector (In this case we start with zero initial values)
	//	
	//	matK1 = an.StructMatrix(); //Storage the global matrix and load vector
	//	fvec = an.Rhs();
	//	
	//	// This code identify singular blocks
	//	TPZStepSolver<REAL> & temp = dynamic_cast<TPZStepSolver<REAL> &> (an.Solver());
	//	std::list <int> & zeropivot = temp.Singular(); 
	//	if (zeropivot.size()) 
	//	{
	//		int eq = * zeropivot.begin();
	//		an.Rhs().Zero();
	//		an.Rhs()(eq,0) = -10000.0;
	//		an.Solve();
	//		TPZFMatrix<REAL> TempSolution = an.Solution();
	//		
	//#ifdef LOG4CXX
	//		// Print the temporal solution
	//		if(logdata->isDebugEnabled())
	//		{
	//			std::stringstream sout;
	//			TempSolution.Print("Singularnodes = ", sout,EMathematicaInput);
	//			LOGPZ_DEBUG(logdata,sout.str())
	//		}
	//#endif	
	//		std::string output;
	//		output = "Singularnodes";
	//		std::stringstream outputfiletemp;
	//		outputfiletemp << output << ".vtk";
	//		std::string plotfile = outputfiletemp.str();
	//		PosProcessMultphysics(meshvec,mphysics,an,plotfile);
	//		
	//	}
	//	
	//#ifdef LOG4CXX
	//	if(logdata->isDebugEnabled())
	//	{
	//		// base system to invert
	//		// just one for checking purpose
	//		std::stringstream sout;
	//		an.Solver().Matrix()->Print("K1 = ", sout,EMathematicaInput);
	//		fvec.Print("fvec = ", sout,EMathematicaInput);		
	//		//Print the temporal solution
	//		Initialsolution.Print("Intial conditions = ", sout,EMathematicaInput);
	//		TPZFMatrix<REAL> Temp;
	//		TPZFMatrix<REAL> Temp2;
	//		matK2->Multiply(Initialsolution,Temp);
	//		Temp.Print("Temp K2 = ", sout,EMathematicaInput);	
	//		LOGPZ_DEBUG(logdata,sout.str())
	//	}
	//#endif
	//	
	//	///start transient problem
	//	SolveSistTransient(delta,MaxTime,matK1, matK2, fvec, Initialsolution, materialist, an, meshvec,  mphysics);
	
	
	return EXIT_SUCCESS;
}


TPZCompMesh*MalhaCompElast(TPZGeoMesh * gmesh,int pOrder, bool Initialstress)
{
	int dim = 2;
	int planestress = 0;	
	// All material with equal properties
	// Needed include diferents materials and properties via keyword file	
	REAL gravity = 9.81; // SI system
	REAL overburdendepth = 2000.0; // SI system
	REAL layerthickness = 10.0;  // SI system	
	TPZVec <REAL> E(TotalMats,10.0e9);
	TPZVec <REAL> poisson(TotalMats,0.35);
	TPZVec <REAL> rockrho(TotalMats,2330.0);	
	TPZVec <REAL> force(2,0.);
	TPZVec <TPZVec< REAL > > forceVec(TotalMats,force);	
	TPZVec <TPZElasticityMaterial *> materialist(TotalMats,0);
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();	
	
	
	
	for(int imat = 0; imat < TotalMats; imat++)
	{	
		
		// Using Gravity field
		forceVec[imat][1]=-gravity*rockrho[imat];
		materialist[imat] = new TPZElasticityMaterial(imat+1, E[imat], poisson[imat], forceVec[imat][0], forceVec[imat][1], planestress); 
		
		
		TPZMaterial * mat(materialist[imat]);
		cmesh->InsertMaterialObject(mat);
		
		// Here is necessary to include the same set of computational information for all computational meshes -> Elastic, Pressure, Biot's Poroelasticity
		TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
		REAL uNUy=0.0;
		val2(1,0)=uNUy;
		
		// 10  displacement	
		TPZFMatrix<REAL> unitdis(2,1,0.);
		unitdis(0,0) = 0.0;
		unitdis(1,0) = 0.0;		
		
		TPZMaterial * BCondNRBTop = materialist[imat]->CreateBC(mat, RBTop,neumann, val1, unitdis);
		cmesh->InsertMaterialObject(BCondNRBTop);
		
		TPZMaterial * BCondNGBTop = materialist[imat]->CreateBC(mat, GBTop,neumann, val1, unitdis);
		cmesh->InsertMaterialObject(BCondNGBTop);
		
		TPZMaterial * BCondNLBTop = materialist[imat]->CreateBC(mat, LBTop,neumann, val1, unitdis);
		cmesh->InsertMaterialObject(BCondNLBTop);		
		
		REAL uNLy=0.0;
		val2(1,0)=uNLy;	
		
		TPZMaterial * BCondNL = materialist[imat]->CreateBC(mat, LBleft,4, val1, val2);
		cmesh->InsertMaterialObject(BCondNL);
		
		
		TPZFMatrix<REAL> val12(2,2,0.), val22(2,1,0.);
		REAL uDL=0.0;
		val22(0,0)=uDL;
		TPZMaterial * BCondDLBBot = materialist[imat]->CreateBC(mat, LBBot,3, val1, val2);
		cmesh->InsertMaterialObject(BCondDLBBot);
		
		TPZFMatrix<REAL> valbot(3,1,0.);
		REAL uDx=0.0;
		REAL uDy=-10.0;
		valbot(0,0)=uDx;
		valbot(1,0)=uDy;		
		
		TPZMaterial * BCondDGBBot = materialist[imat]->CreateBC(mat, GBBot,dirichlet, val1, valbot);
		cmesh->InsertMaterialObject(BCondDGBBot);
		
		TPZMaterial * BCondDRBBot = materialist[imat]->CreateBC(mat, RBBot,3, val1, val2);
		cmesh->InsertMaterialObject(BCondDRBBot);		
		
		REAL uDR=0.0;
		val22(0,0)=uDR;
		TPZMaterial * BCondDR = materialist[imat]->CreateBC(mat, RBright,4, val12, val2);
		cmesh->InsertMaterialObject(BCondDR);
		
//		TPZMaterial * BCondGDR = materialist[imat]->CreateBC(mat, -7,dirichlet, val1, val2);
//		cmesh->InsertMaterialObject(BCondGDR);
//		
//		TPZMaterial * BCondRBDR = materialist[imat]->CreateBC(mat, -8,dirichlet, val1, val2);
//		cmesh->InsertMaterialObject(BCondRBDR);			
		
		
		//	// Point source
		//	TPZMaterial * BCondNpoint = materialist[imat]->CreateBC(mat, pointsource,neumann, val12, val22);
		//	cmesh->InsertMaterialObject(BCondNpoint);
		//	
		//	TPZMaterial * BCyfixedPoints = materialist[imat]->CreateBC(mat, yfixedPoints,neumann, val12, val22);
		//	cmesh->InsertMaterialObject(BCyfixedPoints);
		//	
		//	TPZMaterial * BCxfixedPoints = materialist[imat]->CreateBC(mat, xfixedPoints,neumann, val12, val22);
		//	cmesh->InsertMaterialObject(BCxfixedPoints);	
		
	}
	
	if (Initialstress) 
	{
		
		//		Inserting interface element
		ElasticMatInterface2D *interfacematr = new ElasticMatInterface2D(43,100,0.35, 0.0, -gravity*2500.0, planestress); 
		TPZMaterial * matINteR(interfacematr);
		cmesh->InsertMaterialObject(matINteR);
		
		//		Inserting interface element
		ElasticMatInterface2D *interfacematl = new ElasticMatInterface2D(22,100,0.35, 0.0, -gravity*2500.0, planestress); 
		TPZMaterial * matINteL(interfacematl);
		cmesh->InsertMaterialObject(matINteL);		
		
		
		// setting Computational element structure
		// Definitions of first group pf materials
		
		// Definitions of Second group pf materials		
		std::set<int> iset;
				
		iset.insert(2);
		iset.insert(-12);
		iset.insert(-22);		
		iset.insert(-32);
		iset.insert(-42);		
		cmesh->AutoBuild(iset);
		iset.clear();		
		cmesh->Reference()->ResetReference();
		
		
		iset.insert(3);
		iset.insert(-13);		
		iset.insert(-23);
		iset.insert(-33);
		iset.insert(-43);		
		cmesh->AutoBuild(iset);
		cmesh->Reference()->ResetReference();			
		iset.clear();		
		
		iset.insert(1);
		iset.insert(-11);
		iset.insert(-31);		
		cmesh->AutoBuild(iset);
		cmesh->Reference()->ResetReference();			
		iset.clear();
				
		
//		for(int imat = 0; imat < TotalMats; imat++)
//		{	
//			iset.insert(imat+1);	
//		}
//		cmesh->LoadReferences();
//		
//		cmesh->LoadReferences();
//		iset.insert(-1);	
//		cmesh->AutoBuild(iset);	
//		iset.clear();
//		
//		iset.insert(-2);	
//		cmesh->AutoBuild(iset);	
//		iset.clear();		
//		
//		iset.insert(-3);	
//		cmesh->AutoBuild(iset);	
//		iset.clear();
//		
//		iset.insert(-4);	
//		cmesh->AutoBuild(iset);	
//		iset.clear();
//		
//		iset.insert(-7);	
//		cmesh->AutoBuild(iset);	
//		iset.clear();
//		
//		iset.insert(-8);	
//		cmesh->AutoBuild(iset);	
//		iset.clear();	
		
		
		// Load references for search interface locations
		cmesh->LoadReferences();		
		int ngel = cmesh->Reference()->NElements();
		
		for(int igel = 0; igel < ngel ; igel++ )
		{
			TPZGeoEl * Gel = gmesh->ElementVec()[igel];
			if (!Gel) {
				continue;
			}
			int MatId = Gel->MaterialId();
			if (MatId == 43) {
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
				// Inserting Interace via previus method				
				new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
				
			}
			
			if (MatId == 22) {
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
				// Inserting Interace via previus method				
				new TPZInterfaceElement(* cmesh, Gel, gelindex, neigh[0], neigh[1] );		
				
			}			
			
		}
		iset.insert(43);
		iset.insert(22);		
		cmesh->AutoBuild(iset);
		cmesh->Reference()->ResetReference();			
		iset.clear();		
		
		
	}
	else 
	{
//		cmesh->AutoBuild();
	}
	
	
	return cmesh;
	
}

TPZCompMesh*MalhaCompPressao(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	// All material with equal properties
	// Needed include diferents materials and properties via keyword file
	TPZVec <REAL> convdir(3,0.);
	TPZVec <REAL> diff(TotalMats,0.1);
	TPZVec <REAL> conv(TotalMats,0.0);
	TPZVec <TPZVec< REAL > > convdirVec(TotalMats,convdir);
	TPZVec <REAL> flux(TotalMats,0.0);
	
	TPZVec <TPZMatPoisson3d *> materialist(TotalMats,0);	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();	
	
	for(int imat = 0; imat < TotalMats; imat++)
	{
		materialist[imat] = new TPZMatPoisson3d(imat+1,dim); 
		TPZMaterial * mat(materialist[imat]);
		
		materialist[imat]->SetParameters(diff[imat], conv[imat], convdirVec[imat]);
		materialist[imat]->SetInternalFlux(flux[imat]);
		materialist[imat]->NStateVariables();
		
		cmesh->InsertMaterialObject(mat);
		
		// Here is necessary to include the same set of computational information for all computational meshes -> Elastic, Pressure, Biot's Poroelasticity
		TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
		REAL pN=0.;
		val2(0,0)=pN;
		TPZMaterial * BCondNL = materialist[imat]->CreateBC(mat, bcBottom,neumann, val1, val2);
		cmesh->InsertMaterialObject(BCondNL);
		
		TPZMaterial * BCondNU = materialist[imat]->CreateBC(mat, bcTop,neumann, val1, val2);
		cmesh->InsertMaterialObject(BCondNU);
		
		TPZFMatrix<REAL> val12(2,2,0.), val22(2,1,0.);
		REAL uDL=4000.;
		val22(0,0)=uDL;
		TPZMaterial * BCondDL = materialist[imat]->CreateBC(mat, bcLeft,dirichlet, val12, val22);
		cmesh->InsertMaterialObject(BCondDL);
		
		REAL uDR=3000.;
		val22(0,0)=uDR;
		TPZMaterial * BCondDR = materialist[imat]->CreateBC(mat, bcRight,dirichlet, val12, val22);
		cmesh->InsertMaterialObject(BCondDR);
		
		// Point source
		TPZMaterial * BCondNpoint = materialist[imat]->CreateBC(mat, pointsource,neumann, val12, val22);
		cmesh->InsertMaterialObject(BCondNpoint);
		
		TPZMaterial * BCyfixedPoints = materialist[imat]->CreateBC(mat, yfixedPoints,neumann, val12, val22);
		cmesh->InsertMaterialObject(BCyfixedPoints);	
		
		TPZMaterial * BCxfixedPoints = materialist[imat]->CreateBC(mat, xfixedPoints,neumann, val12, val22);
		cmesh->InsertMaterialObject(BCxfixedPoints);		
	}
	
	
	//	// setting Computational element structure
	//	std::set<int> iset;	
	//	for(int imat = 0; imat < TotalMats; imat++)
	//	{	
	//		iset.insert(imat+1);	
	//		cmesh->AutoBuild(iset);
	//		cmesh->LoadReferences(); 		
	//		cmesh->Reference()->ResetReference();			
	//		iset.clear();
	//		
	//	}
	//	cmesh->Reference()->ResetReference();		
	cmesh->AutoBuild();	
	
	
	return cmesh;
}

TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial)
{
	//Creating computational mesh for multiphysic elements
	int dim = 2;
	int planestress = 0; // This is a Plain strain problem
	
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
	mphysics->SetAllCreateFunctionsMultiphysicElem();
	TPZAutoPointer<TPZFunction<STATE> > TimeDepForcingF;
	TPZAutoPointer<TPZFunction<STATE> > TimeDepFExact;
	//	TimeDepForcingF = new TPZDummyFunction<STATE>(ForcingTimeDependFunction);
//	TimeDepFExact = new TPZDummyFunction<STATE>(SolucaoExata2DLinesource);
	
	// Definitions
	REAL lamb = 8333.33;					//	[Pa]
	REAL lambu = 13.4e9;						//	[Pa]
	REAL alpha = 0.7;						//	[-]
	REAL G= 12500.0;							//	[Pa]
	REAL rhof = 1000.0;						//	[kg/m3]
	REAL poisson = (lamb)/(2*(lamb+G));						//	[-]
	REAL Eyoung = (G*(3.0*lamb+2.0*G))/(lamb+G);					//	[Pa]			
	REAL Ssig = ((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
	REAL Se =  ((pow(alpha,2))/((lambu-lamb)));
	REAL rockrho = 2500.0;					//	[kg/m3]
	REAL gravity = 9.81;					//	[m/s2]
	REAL c = 0.083;							//	[m2/s]	
	REAL visc = 0.001;						//	[Pa.s]
	REAL perm = c*Ssig*visc;//7.8784288e-15;//1.109542e-14						//	[m2]
	REAL qo = -0.001;
	REAL Bo = 1.0;
	REAL PI = atan(1.)*4.;
	REAL Yforce = -1000.0;
	TPZVec <REAL> force(2,0.);	
	REAL fx=0.0;
	REAL fy=-rockrho*gravity;
	
	
	// Definitions for paramter vectors Al materials have the same porperties
	TPZVec <REAL> lambVec(TotalMats,lamb);
	TPZVec <REAL> lambuVec(TotalMats,lambu);
	TPZVec <REAL> alphaVec(TotalMats,alpha);
	TPZVec <REAL> GVec(TotalMats,G);
	TPZVec <REAL> rhofVec(TotalMats,rhof);
	TPZVec <REAL> poissonVec(TotalMats,poisson);
	TPZVec <REAL> EyoungVec(TotalMats,Eyoung);
	TPZVec <REAL> SsigVec(TotalMats,Ssig);
	TPZVec <REAL> SeVec(TotalMats,Se);
	TPZVec <REAL> rockrhoVec(TotalMats,rockrho);
	TPZVec <REAL> cVec(TotalMats,c);
	TPZVec <REAL> viscVec(TotalMats,visc);
	TPZVec <REAL> permVec(TotalMats,perm);
	TPZVec <REAL> qoVec(TotalMats,qo);
	TPZVec <REAL> BoVec(TotalMats,Bo);
	TPZVec <REAL> YforceVec(TotalMats,Yforce);
	TPZVec <TPZVec< REAL > > forceVec(TotalMats,force);
	for(int imat = 0; imat < TotalMats; imat++)
	{	
		
		// Using Gravity field
		forceVec[imat][1]=gravity*rockrhoVec[imat];		
		
		mymaterial[imat] = new TPZPoroElastic2d (imat+1, dim);
		mymaterial[imat]->SetParameters(EyoungVec[imat], poissonVec[imat],forceVec[imat][0], forceVec[imat][1]);
		mymaterial[imat]->SetParameters(permVec[imat],viscVec[imat]);
		mymaterial[imat]->SetfPlaneProblem(planestress);
		mymaterial[imat]->SetBiotParameters(alpha,Se);
		//	mymaterial[imat]->SetTimeDependentForcingFunction(TimeDepForcingF);
		mymaterial[imat]->SetTimeDependentFunctionExact(TimeDepFExact);		
		
		ofstream argm("mymaterial1.txt");
		mymaterial[imat]->Print(argm);
		TPZMaterial * mat(mymaterial[imat]);
		mphysics->InsertMaterialObject(mat);
		
		///--- --- Inserir condicoes de contorno
		// Radial  Model (1)
		// Squared Model (2) for consistency
		// Squared Model (3) for validation
		// SemiCircular HalfSpace Model (4) for validation	
		// Graven test (5)		
		
		int model = 5;
		
		switch (model) {
			case 1:
			{
				
				//	Parameters
				//	REAL Eyoung = 3.e4;
				//	REAL poisson = 0.2;
				//	REAL alpha=0.0;
				//	REAL Se=1.0;
				//	REAL rockrho = 2330.0; // SI system
				//	REAL gravity = 0.0;//-9.8; // SI system
				//	REAL fx=0.0;
				//	REAL fy=gravity*rockrho;
				//	REAL overburdendepth = 2000.0; // SI system
				//	REAL layerthickness = 10.0;  // SI system
				//	REAL perm = 1.e-10;
				//	REAL visc = 1.e-3;
				
				// elastic problem -> 1 ,  Poisson problem -> 2
				// Boundary condition D1N2 means: Elastic -> Dirichlet, Pressure -> Neumman
				//			TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
				//			int D1D2 = 1;
				//			REAL uDx=0.0;
				//			REAL uDy=0.0;
				//			REAL Pressure=0.0;
				//			val2(0,0)=uDx;
				//			val2(1,0)=uDy;
				//			val2(2,0)=Pressure;
				//			TPZMaterial * BonArc1 = mymaterial[imat]->CreateBC(mat, arc1,D1D2, val1, val2);
				//			mphysics->InsertMaterialObject(BonArc1);
				//			
				//			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				//			TPZMaterial * BonArc2 = mymaterial[imat]->CreateBC(mat, arc2,D1D2, val1, val2);
				//			mphysics->InsertMaterialObject(BonArc2);
				//			
				//			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				//			TPZMaterial * BonArc3 = mymaterial[imat]->CreateBC(mat, arc3,D1D2, val1, val2);
				//			mphysics->InsertMaterialObject(BonArc3);
				//			
				//			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				//			TPZMaterial * BonArc4 = mymaterial[imat]->CreateBC(mat, arc4,D1D2, val1, val2);
				//			mphysics->InsertMaterialObject(BonArc4);
				
				
				//		Production/Injection Well Just Neumman for Pressure
				int WellQ = 400;
				TPZFMatrix<REAL> valP13(3,2,0.), valP23(3,1,0.);
				REAL MassQ = qo/rhof;
				valP23(2,0)=MassQ;
				TPZMaterial * BCPoint = mymaterial[imat]->CreateBC(mat, WellPoint, WellQ, valP13, valP23);
				mphysics->InsertMaterialObject(BCPoint);	
				//		-----------
				
			}
				
				break;
			case 2:
			{
				
				// elastic problem -> 1 ,  Poisson problem -> 2
				// Boundary condition D1D2 means: Elastic -> Dirichlet, Pressure -> Dirichlet
				TPZFMatrix<REAL> val12(3,2,0.), val22(3,1,0.);
				int D1D2 = 100;
				REAL uDx=0.0;
				REAL uDy=0.0;
				REAL Pressure=1000.0;
				val22(0,0)=uDx;
				val22(1,0)=uDy;
				val22(2,0)=Pressure;
				
				TPZMaterial * BCondT = mymaterial[imat]->CreateBC(mat, bcTop,D1D2, val12, val22);
				mphysics->InsertMaterialObject(BCondT);	
				
				TPZMaterial * BCondBt = mymaterial[imat]->CreateBC(mat, bcBottom,D1D2, val12, val22);
				mphysics->InsertMaterialObject(BCondBt);
				
				TPZMaterial * BCondL = mymaterial[imat]->CreateBC(mat, bcLeft, D1D2, val12, val22);
				mphysics->InsertMaterialObject(BCondL);
				
				TPZMaterial * BCondR = mymaterial[imat]->CreateBC(mat, bcRight, D1D2, val12, val22);
				mphysics->InsertMaterialObject(BCondR);
				
				//		Production/Injection Well Just Neumman for Pressure
				int WellQ = 400;
				TPZFMatrix<REAL> valPP12(3,2,0.), valPP22(3,1,0.);
				REAL MassQ=10.0;
				valPP22(2,0)=MassQ;
				TPZMaterial * BCPoint = mymaterial[imat]->CreateBC(mat, pointsource, WellQ, valPP12, valPP22);
				mphysics->InsertMaterialObject(BCPoint);	
				//		-----------
				
			}
				
				break;
			case 3:
			{
				
				//	a poroelastic contribution to the reservoir stress path
				//	Parameters
				//	REAL Eyoung = 1.43e10;					//	[Pa]
				//	REAL poisson = 0.3;						//	[-]
				//	REAL alpha=0.7;							//	[-]
				//	REAL Se=0.0;//9.60784e-11;				//	[-]
				//	REAL rockrho = 2330.0;					//	[kg/m3]
				//	REAL gravity = 0.0;//-9.8;				//	[m/s2]
				//	REAL perm = 8.2e-8;						//	[m2]
				//	REAL visc = 1.e-3;						//	[Pa.s]
				//	REAL c = 0.082;							//	[m2/s]
				
				//	REAL fx=0.0;
				//	REAL fy=gravity*rockrho;			
				
				// elastic problem -> 1 ,  Poisson problem -> 2
				// Boundary condition D1D2 means: Elastic -> N, Pressure -> N
				TPZFMatrix<REAL> valall12(3,2,0.), valall22(3,1,0.);
				int N1N2 = 11;
				REAL uDx=0.0;
				REAL uDy=0.0;
				REAL Pressure=0.0;
				valall22(0,0)=uDx;
				valall22(1,0)=uDy;
				valall22(2,0)=Pressure;
				
				TPZMaterial * BCondT = mymaterial[imat]->CreateBC(mat, bcTop,N1N2, valall12, valall22);
				mphysics->InsertMaterialObject(BCondT);
				TPZMaterial * BCondB = mymaterial[imat]->CreateBC(mat, bcBottom,N1N2, valall12, valall22);
				mphysics->InsertMaterialObject(BCondB);
				TPZMaterial * BCondL = mymaterial[imat]->CreateBC(mat, bcLeft,N1N2, valall12, valall22);
				mphysics->InsertMaterialObject(BCondL);
				TPZMaterial * BCondR = mymaterial[imat]->CreateBC(mat, bcRight,N1N2, valall12, valall22);
				mphysics->InsertMaterialObject(BCondR);
				
				/// Setting drained free surface  Boundary condition N1D2 means: Elastic -> Neumman, Pressure -> Dirichlet
				TPZFMatrix<REAL> val11(3,2,0.), val21(3,1,0.);
				int DYFX1N2 = 200;
				REAL sigmax = 0.0;
				REAL sigmay = 0.0;
				REAL ptop=0.0;
				val21(0,0)=sigmax;
				val21(1,0)=sigmay;
				val21(2,0)=ptop;
				TPZMaterial * BCyfixedPoints = mymaterial[imat]->CreateBC(mat, yfixedPoints, DYFX1N2, val11, val21);
				mphysics->InsertMaterialObject(BCyfixedPoints);
				
				/// Setting free displacement non permeable surface  Boundary condition DXFY1N2 means: Elastic -> Dirichlet x component, Pressure -> Neumman
				int DXFY1N2 = 300;
				TPZFMatrix<REAL> val13(3,2,0.), val23(3,1,0.);
				REAL uDleftx = 0.0;
				REAL uDlefty = 0.0;			
				REAL pNleft  = 0.0;
				val23(0,0)=uDleftx;
				val23(1,0)=uDlefty;			
				val23(2,0)=pNleft;
				TPZMaterial * BCxfixedPoints = mymaterial[imat]->CreateBC(mat, xfixedPoints, DXFY1N2, val13, val23);
				mphysics->InsertMaterialObject(BCxfixedPoints);			
				//			
				//			/// Setting free displacement non permeable surface  Boundary condition DYFX1N2 means: Elastic -> Dirichlet y component, Pressure -> Neumman
				//			TPZFMatrix<REAL> val12(3,2,0.), val22(3,1,0.);
				//			int DYFX1N2 = 200; 			
				//			REAL uDbotx=0.0;
				//			REAL uDboty=0.0;
				//			REAL pNbot=0.0;
				//			val22(0,0)=uDbotx;
				//			val22(1,0)=uDboty;
				//			val22(2,0)=pNbot;
				//			TPZMaterial * BCondBt = mymaterial[imat]->CreateBC(mat, bcBottom, DYFX1N2, val12, val22);
				//			mphysics->InsertMaterialObject(BCondBt);
				
				//			/// Setting free displacement non permeable surface  Boundary condition DXFY1N2 means: Elastic -> Dirichlet x component, Pressure -> Neumman
				//			int DXFY1N2 = 300;
				//			TPZFMatrix<REAL> val13(3,2,0.), val23(3,1,0.);
				//			REAL uDleftx = 0.0;
				//			REAL uDlefty = 0.0;			
				//			REAL pNleft  = 0.0;
				//			val23(0,0)=uDleftx;
				//			val23(1,0)=uDlefty;			
				//			val23(2,0)=pNleft;
				//			TPZMaterial * BCondL = mymaterial[imat]->CreateBC(mat, bcLeft, DXFY1N2, val13, val23);
				//			mphysics->InsertMaterialObject(BCondL);
				
				//			/// Setting free displacement non permeable surface  Boundary condition DXFY1N2 means: Elastic -> Dirichlet x component, Pressure -> Neumman			
				//			TPZFMatrix<REAL> val14(3,2,0.), val24(3,1,0.);
				//			int D1D2 = 1;
				//			REAL uDrightx = 0.0;
				//			REAL uDrighty = 0.0;			
				//			REAL pNright  = 0.0;
				//			val24(0,0)=uDrightx;
				//			val24(1,0)=uDrighty;			
				//			val24(2,0)=pNright;
				//			TPZMaterial * BCondR = mymaterial[imat]->CreateBC(mat, bcRight, DXFY1N2, val14, val24);
				//			mphysics->InsertMaterialObject(BCondR);
				
				//	Production or injection Well point (Q is the injected/produced volume per time [kg/s])
				int Well = 400;
				TPZFMatrix<REAL> valP13(3,2,0.), valP23(3,1,0.);
				valP23(2,0)=qo/rhof;
				TPZMaterial * BCPoint = mymaterial[imat]->CreateBC(mat, WellPoint, Well, valP13, valP23);
				mphysics->InsertMaterialObject(BCPoint);	
				//-----------
				
			}
				break;
			case 4:
			{
				
				//			// Definitions Flamant problem
				//			REAL lamb = 1.0e9;					//	[Pa]
				//			REAL alpha = 0.0;						//	[-]
				//			REAL G= 1.0e9;							//	[Pa]
				//			REAL rhof = 1.0;						//	[kg/m3]
				//			REAL poisson3D = 0.3;						//	[-]
				//			REAL Eyoung3D = (G*(3.0*lamb+2.0*G))/(lamb+G);					//	[Pa]	
				//			REAL poisson = (poisson3D)/(1-poisson3D);						//	[-]
				//			REAL Eyoung = (Eyoung3D)/(1-pow(poisson3D,2.0));					//	[Pa]		
				//			REAL Se = 0.0;//((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
				//			REAL rockrho = 0.0;					//	[kg/m3]
				//			REAL gravity = 0.0;					//	[m/s2]
				//			REAL c = 1.0;							//	[m2/s]	
				//			REAL visc = 1.0;						//	[Pa.s]
				//			REAL perm = c*Se;//7.8784288e-15;//1.109542e-14						//	[m2]
				//			REAL qo = 1.0;			
				//			REAL Yforce = 1.0;
				//			REAL Bo = 1.0;
				//			REAL PI = atan(1.)*4.;		
				//			
				//			REAL fx=0.0;
				//			REAL fy=0.0;			
				
				// elastic problem -> 1 ,  Poisson problem -> 2
				// Boundary condition D1N2 means: Elastic -> Dirichlet, Pressure -> Neumman
				TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
				int D1D2 = 0;
				int N1D2 = 10;			
				REAL uDx=0.0;
				REAL uDy=0.0;
				REAL Pressure=0.0;
				val2(0,0)=uDx;
				val2(1,0)=uDy;
				val2(2,0)=Pressure;
				TPZMaterial * BonArc1 = mymaterial[imat]->CreateBC(mat, bcLeft,D1D2, val1, val2);
				mphysics->InsertMaterialObject(BonArc1);
				
				///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				TPZMaterial * BonArc2 = mymaterial[imat]->CreateBC(mat, bcRight,D1D2, val1, val2);
				mphysics->InsertMaterialObject(BonArc2);
				
				///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				TPZMaterial * BonArc3 = mymaterial[imat]->CreateBC(mat, bcBottom,D1D2, val1, val2);
				mphysics->InsertMaterialObject(BonArc3);
				
				///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				TPZMaterial * BonArc4 = mymaterial[imat]->CreateBC(mat, bcTop,N1D2, val1, val2);
				mphysics->InsertMaterialObject(BonArc4);			
				
				//		loaded vertical point force Y
				int N1N2 = 11;
				TPZFMatrix<REAL> valP13(3,2,0.), valP23(3,1,0.);
				// 170 Elasticity Theory
				valP23(1,0)=Yforce;
				TPZMaterial * BCPoint = mymaterial[imat]->CreateBC(mat, WellPoint, N1N2, valP13, valP23);
				mphysics->InsertMaterialObject(BCPoint);	
				//		-----------
				
			}			
				break;
			case 5:
			{
				
				// elastic problem -> 1 ,  Poisson problem -> 2
				// Boundary condition D1N2 means: Elastic -> Dirichlet, Pressure -> Neumman
				TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
				int D1D2 = 10;
				int N1N2 = 200;
				int D1FN2 = 300;		
				REAL uDx=0.0;
				REAL uDy=0.0;
				REAL Pressure=0.0;
				val2(0,0)=uDx;
				val2(1,0)=uDy;
				val2(2,0)=Pressure;
				TPZMaterial * BonArc1 = mymaterial[imat]->CreateBC(mat, GBTop,D1D2, val1, val2);
				mphysics->InsertMaterialObject(BonArc1);
				
				///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				TPZMaterial * BonArc2 = mymaterial[imat]->CreateBC(mat, GBBot,N1N2, val1, val2);
				mphysics->InsertMaterialObject(BonArc2);
				
				///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				TPZMaterial * BonArc3 = mymaterial[imat]->CreateBC(mat, LBleft,D1FN2, val1, val2);
				mphysics->InsertMaterialObject(BonArc3);
				
				///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
				TPZMaterial * BonArc4 = mymaterial[imat]->CreateBC(mat, RBright,D1FN2, val1, val2);
				mphysics->InsertMaterialObject(BonArc4);			
				
				//			//		Production/Injection Well Just Neumman for Pressure
				//			int WellQ = 400;
				//			TPZFMatrix<REAL> valP13(3,2,0.), valP23(3,1,0.);
				//			REAL MassQ = qo/rhof;
				//			valP23(2,0)=MassQ;
				//			TPZMaterial * BCPoint = mymaterial[imat]->CreateBC(mat, WellPoint, WellQ, valP13, valP23);
				//			mphysics->InsertMaterialObject(BCPoint);
				//		-----------
				
			}			
				break;			
			default:
			{
#ifdef LOG4CXX
				{
					LOGPZ_DEBUG(logger, "Nothing to say, nothig to do")
				}
#endif			
				
			}
				
				break;
		}
		
	}
	
	
	//	// setting Computational element structure
	//	std::set<int> iset;	
	//	for(int imat = 0; imat < TotalMats; imat++)
	//	{	
	//		iset.insert(imat+1);	
	//		mphysics->AutoBuild(iset);
	//		mphysics->LoadReferences(); 		
	//		mphysics->Reference()->ResetReference();			
	//		iset.clear();
	//		
	//	}
	////	
	////	mphysics->Reference()->ResetReference();		
	mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
	
	// Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
	
#ifdef LOG4CXX
	{
		//        std::stringstream sout;
		//        mphysics->Print(sout);
		//        LOGPZ_DEBUG(logger, sout.str())
	}
#endif
	ofstream arg("mphysic.txt");
	mphysics->Print(arg);
	
	return mphysics;
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

void SolveSistTransient(REAL deltime,REAL maxtime,TPZFMatrix<REAL> matK1, TPZAutoPointer < TPZMatrix<REAL> > matK2, TPZFMatrix<REAL> fvec, TPZFMatrix<REAL> &Initialsolution, TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZAnalysis &an,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics){
	
	int nrows;
	nrows = matK2->Rows();
	TPZFMatrix<REAL> TotalRhs(nrows,1,0.0);
	TPZFMatrix<REAL> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<REAL> Lastsolution = Initialsolution;
	
	std::string outputfile;
	outputfile = "TransientSolution";
	
	REAL delt = deltime;
	REAL Maxtime = maxtime;
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*delt; 
	while (TimeValue < Maxtime)
	{	
		// This time solution i for Transient Analytic Solution
		for(int imat = 0; imat < TotalMats; imat++)
		{
			mymaterial[imat]->SetTimeValue(TimeValue);
		}
		matK2->Multiply(Lastsolution,TotalRhstemp);
		TotalRhs = fvec + TotalRhstemp;
		an.Rhs() = TotalRhs;
		
		
#ifdef LOG4CXX
		if(logdata->isDebugEnabled())
		{
			//	Print the temporal solution
			//			std::stringstream sout;
			//			matK1.Print("Temporal Solution = ", sout,EMathematicaInput);
			//			TotalRhstemp.Print("Temporal Solution = ", sout,EMathematicaInput);
		}
#endif	
		
		an.Solve(); //	Solve the current Linear system
		Lastsolution = an.Solution(); //	Save the current solution
		//Lastsolution.Print("Sol = ");
		
#ifdef LOG4CXX
		//Print the temporal solution
		//		if(logdata->isDebugEnabled()){
		//			std::stringstream sout;
		//			Lastsolution.Print("Temporal Solution = ", sout,EMathematicaInput);
		//			LOGPZ_DEBUG(logdata,sout.str())
		//		}
#endif
		
		//General post-processing
		//TPZBuildMultiphysicsMesh * Objectdumy;
		//Objectdumy->TransferFromMultiPhysics(meshvec, mphysics);
		std::stringstream outputfiletemp;
		outputfiletemp << outputfile << ".vtk";
		std::string plotfile = outputfiletemp.str();
		PosProcessMultphysics(meshvec,mphysics,an,plotfile);		
		
		// Total mass calculation
		//
		
		//		int totalel=mphysics->NElements();
		//		int vartointegrate = 13;
		//		REAL TotalMass=0;
		//		
		//		
		//		
		//		for ( int el = 0; el < totalel; el++ )
		//		{
		//		
		//			TPZVec <REAL> result(1,0);
		//			TPZCompEl * celvar = mphysics->ElementVec()[el];
		//			if (celvar->Reference()->MaterialId() > 0 )// avoid bc conditions
		//			{
		//				celvar->Integrate(vartointegrate, result);
		//				TotalMass += result[0];
		//				
		//			}
		//			
		//		}
		//		
		//		cout << "Numerical total mass  " << TotalMass << "  at time " << TimeValue << endl;
		//		cout << "Theorical total mass  " << 0.000001*TimeValue << "  at time " << TimeValue << endl;		
		
		//
		// End TotalMass Calculation
		
		// Next Calculation
		cent++;
		TimeValue = cent*delt;
	}
	
	
}

void PosProcess(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
	vecnames[0]= "MinusKGradU";
	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void PosProcess2(TPZAnalysis &an, std::string plotfile){
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

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
{
	TPZBuildMultiphysicsMesh * Objectdumy;
	Objectdumy->TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(13), vecnames(3);
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
	scalnames[12] = "ETAUXY";	
	vecnames[0]  = "Displacement";
	vecnames[1]  = "FluxVector";
	vecnames[2]  = "EDisplacement";	
	
	
	
	const int dim = 2;
	int div =0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh){
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gMesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			int GelDimension = gel->Dimension();
			if (GelDimension == 2 || GelDimension == 1) 
			{
				gel->Divide (filhos);
			}
		}//for i
	}//ref
}

void RefinamentoUniforme(TPZGeoMesh * gMesh, int nh, int MatId, int indexEl){
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

void RefinElemComp(TPZCompMesh  *cMesh, int indexEl){
	
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

void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv){
	
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
	
	//    RefPatternMesh();
	//TPZGeoMesh * gmesh = refp->Mesh();
	PrintGMeshVTK(gmesh, file);
}

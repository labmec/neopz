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


#include <cmath>
#include <set>


//	Use this tools for Exponential Integral Fuction (Implementation defined at end of this main file) 
#include <math.h>		// required for fabsl(), expl() and logl()        
#ifdef USING_BOOST
	#include <boost/math/special_functions/erf.hpp> //Required for erfc function on windows
#endif
#include <float.h>		// required for LDBL_EPSILON, DBL_MAX
//	Internally Defined Routines
double      Exponential_Integral_Ei( double x );
long double xExponential_Integral_Ei( long double x );
static long double Continued_Fraction_Ei( long double x );
static long double Power_Series_Ei( long double x );
static long double Argument_Addition_Series_Ei( long double x);
//	Internally Defined Constants
static const long double epsilon = 10.0 * LDBL_EPSILON;
//	End Use this tools for Exponential Integral Fuction (Implementation defined at end of this main file)

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
const int RBright = -2;
const int RBleft = -4;
const int RBBot = -1;
const int RBTop = -3;

// Left Block 
const int LBright = -2;
const int LBleft = -4;
const int LBBot = -1;
const int LBTop = -3;

// Graven Block 
const int GBright = -2;
const int GBleft = -4;
const int GBBot = -1;
const int GBTop = -3;

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
TPZCompMesh *MalhaCompElast(TPZGeoMesh * gmesh,int pOrder);
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

// This Methods are Related whit Analytic Solutions of direferent Validations Problems (in this Main are created 5 problems for validation)
// Analytical solution for paper -> X
void SolucaoExata1D(TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
// Analytical solution for paper -> X
void DeslocamentoYExata(TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
// Analytical solution for paper -> X
void SigmaYExata(TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
// Analytical solution for paper -> Y
void SolucaoExata2D(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
// Analytical solution for paper -> Z
void SolucaoExata2DLinesource(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
// Analytical solution for book -> Rosa seila
void SolucaoExataRosa1D(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
// Analytical solution for book -> Rosa seila
void SolucaoExataRosa1DPseudo(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
// Analytical solution for flamant problem Elasticity Theory: Applications and Numerics 
void ExactSolutionFlamantProblem(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);


//	Forcing transient function definition
void ForcingTimeDependFunction(TPZVec<REAL> &loc, REAL TimeValue,int WhichStateVariable,double &StateVariable);



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


//	Rectangular geometric mesh
	
	
	MeshGeneration mygenerator;
	TPZVec <int> matIdlist(8,0);
	matIdlist[0]= matId;
	matIdlist[1]= bcBottom;
	matIdlist[2]= bcRight;
	matIdlist[3]= bcTop;
	matIdlist[4]= bcLeft;
	matIdlist[5]= WellLine;
	matIdlist[6]= xfixedPoints;
	matIdlist[7]= yfixedPoints;	
	REAL ModelLenght = 100000.0;
	REAL ModelHight = 100000.0;
	REAL Modelthickness = 1.0;
	mygenerator.Setdimensions(ModelLenght, ModelHight, Modelthickness);		
	TPZGeoMesh * gmesh = mygenerator.MalhaGeom(matIdlist);
	

//	End Rectangular geometric mesh	
	
// Begin Circular geometric mesh using tpzarc3d

//	MeshGeneration mygenerator;
//	TPZVec <int> matIdlist(6,0);
//	matIdlist[0]= matId;
//	matIdlist[1]= arc1;
//	matIdlist[2]= arc1;
//	matIdlist[3]= arc1;
//	matIdlist[4]= arc1;
//	matIdlist[5]= WellLine;	
//	REAL ModelRadius = 100000.0;
//	REAL Modelthickness = 1.0;
//	mygenerator.Setdimensions(ModelRadius, Modelthickness);
//	TPZGeoMesh * gmesh = mygenerator.GeometricMesh2DValidation(matIdlist);	

// End Circular mesh	
	
//	Rectangular geometric mesh for graven analysis
	
//	int nLayers			= 5; 
//	bool InterfaceEl	= false; 
//	REAL LlengthFootFault	= 1000.0;
//	REAL DipFaultAngleleft		= 45.0;
//	REAL DipFaultAngleright		= 45.0;
//	REAL WellFaultlength		= 500.0;
//	TPZVec <bool> wichProductionlayer(nLayers+1,false);
////	wichProductionlayer[0] = true; // Horizons -> Top
////	wichProductionlayer[1] = true; // Horizons -> Top	
////	wichProductionlayer[2] = true; // Horizons -> Top
////	wichProductionlayer[2] = true; // Horizons -> Top
////	wichProductionlayer[3] = true; // Horizons -> Top
//
//	MeshGeneration mygenerator;
//	TPZGeoMesh * gmesh = mygenerator.MalhaGeoGravenobj(nLayers,LlengthFootFault,DipFaultAngleleft,DipFaultAngleright,WellFaultlength,wichProductionlayer,InterfaceEl);		

// End Rectangular geometric mesh for graven analysis

	
	
// Rectangular geometric mesh using TPZGenGrid 
	
//	TPZVec < int > nx(2);
//	TPZVec < REAL > corx0(2);
//	TPZVec < REAL > corx1(2);
//	int  	numlayer = 1; // Layers Numbers
//	REAL  	rotation = 0.5; // For testing purpose 
//	
//	// refinement level
//	nx[0] = 2;
//	nx[1] = 2;
//	//	x0	lower left coordinate
//	corx0[0] = 0.0;
//	corx0[1] = 0.0;	
//	//	x1	upper right coordinate 
//	corx1[0] = 50000.0;
//	corx1[1] = 50000.0;	
//	
//	TPZGenGrid geomesh(nx,corx0,corx1,numlayer);
//	TPZGeoMesh * gmesh = new TPZGeoMesh;
//	geomesh.Read(*gmesh);
//	
//	// Setting BC conditions
//	geomesh.SetBC(gmesh,0,bcBottom);
//	geomesh.SetBC(gmesh,1,bcRight);
//	geomesh.SetBC(gmesh,2,bcTop);
//	geomesh.SetBC(gmesh,3,bcLeft);	
//	
//	TPZVec < REAL > PointSourceCor(3);
//	// Injection point in the center of the model
//	PointSourceCor[0]=25000.0;
//	PointSourceCor[1]=25000.0;
//	PointSourceCor[2]=0.0;	
//	geomesh.SetPointBC(gmesh,PointSourceCor, pointsource);	
//	
//	gmesh->BuildConnectivity();	

// End Rectangular geometric mesh using TPZGenGrid
	

	
// Begin Half Sapce for Flamant Problem Semi-Circular geometric mesh using tpzarc3d
	
//	int nodenumber = 6;
//	REAL ModelRadius = 30000.0;
//	TPZGeoMesh * gmesh = new TPZGeoMesh;
//	gmesh->NodeVec().Resize(nodenumber);
//	
//	// Setting node coordantes for Arc3D 1
//	int id = 0;
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
//	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
//	id++;
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
//	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
//	id++;
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
//	gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y	
//	id++;	
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
//	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
//	id++;	
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,-sqrt(2)*ModelRadius/2.);//coord X
//	gmesh->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y	
//	id++;
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,sqrt(2)*ModelRadius/2.);//coord X
//	gmesh->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y	
//	id++;
//	
//	int elementid = 0;
//	// Create Geometrical Arc #1
//	// Definition of Arc coordenates
//	TPZVec < int > nodeindex(3,0.0);
//	nodeindex[0] = 1;	
//	nodeindex[1] = 4;
//	nodeindex[2] = 2;
//	TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc1 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1,*gmesh);
//	elementid++;
//	
//	// Create Geometrical Arc #2
//	nodeindex[0] = 2;	
//	nodeindex[1] = 5;
//	nodeindex[2] = 3;			
//	TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc2 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2,*gmesh);
//	elementid++;
//
//	// Create Geometrical triangle #1	
//	nodeindex.resize(3);
//	nodeindex[0] = 0;
//	nodeindex[1] = 1;
//	nodeindex[2] = 2;	
//	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle> > *triangle1 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle> > (elementid,nodeindex, matId,*gmesh);
//	elementid++;
//	
//	// Create Geometrical triangle #2		
//	nodeindex[0] = 0;
//	nodeindex[1] = 2;
//	nodeindex[2] = 3;		
//	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle> > *triangle2 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle> > (elementid,nodeindex, matId,*gmesh);
//	elementid++;
//	
//	// Create Geometrical Line #1	
//	nodeindex.resize(2);
//	nodeindex[0] = 3;
//	nodeindex[0] = 0;	
//	TPZGeoElRefPattern < pzgeom::TPZGeoLinear > * elline1 = new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc3,*gmesh); 
//	elementid++;
//	
//	// Create Geometrical Line #2	
//	nodeindex.resize(2);
//	nodeindex[0] = 0;
//	nodeindex[0] = 1;	
//	TPZGeoElRefPattern < pzgeom::TPZGeoLinear > * elline2 = new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc4,*gmesh); 
//	elementid++;	
//	
//	// Create Geometrical Point for fluid injection or Production #1	
//	nodeindex.resize(1);
//	nodeindex[0] = 0;	
//	TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, WellPoint,*gmesh); 
//	elementid++;
	
//	gmesh->BuildConnectivity();

//	TPZVec < int > nx(2);
//	TPZVec < REAL > corx0(2);
//	TPZVec < REAL > corx1(2);
//	int  	numlayer = 1; // Layers Numbers
//	REAL  	rotation = 0.5; // For testing purpose 
//	
//	// refinement level
//	nx[0] = 2;
//	nx[1] = 2;
//	//	x0	lower left coordinate
//	corx0[0] = 0.0;
//	corx0[1] = 0.0;	
//	//	x1	upper right coordinate 
//	corx1[0] = 100000.0;
//	corx1[1] = 50000.0;	
//	
//	TPZGenGrid geomesh(nx,corx0,corx1,numlayer);
//	TPZGeoMesh * gmesh = new TPZGeoMesh;
//	geomesh.Read(*gmesh);
//	
//	// Setting BC conditions
//	geomesh.SetBC(gmesh,0,bcBottom);
//	geomesh.SetBC(gmesh,1,bcRight);
//	geomesh.SetBC(gmesh,2,bcTop);
//	geomesh.SetBC(gmesh,3,bcLeft);	
//	
//	TPZVec < REAL > PointSourceCor(3);
//	// Injection point in the center of the model
//	PointSourceCor[0]=50000.0;
//	PointSourceCor[1]=50000.0;
//	PointSourceCor[2]=0.0;	
//	geomesh.SetPointBC(gmesh,PointSourceCor, pointsource);	
//
//	gmesh->BuildConnectivity();		
	
//	
// End Half Sapce for Flamant Problem Semi-Circular geometric mesh using tpzarc3d	
	
// Use this for irregular mesh created with GID format
// Not implemented	
// End Use this for irregular mesh created with GID format	

// END GEOMETRICAL MESH CREATION	


//	Print Geometrical Base Mesh
	ofstream arg1("gmesh1.txt");
	gmesh->Print(arg1);
	ofstream file1("malhageoInicial.vtk");	
	//	In this option true -> let you use shrink paraview filter
	//	shrink filter in paraview just use for "amazing" visualization
//	PrintGMeshVTK(gmesh,file1);	
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file1, true);	
	
	
// Directional Point source Refinement 
	int ndirectdivp = 17;
	set<int> SETmatPointRefDir2;
	SETmatPointRefDir2.insert(WellPoint);	
	for(int j = 0; j < ndirectdivp; j++)
	{
		int nel = gmesh->NElements();
		for (int iref = 0; iref < nel; iref++)
		{
			TPZVec<TPZGeoEl*> filhos;
			TPZGeoEl * gelP2 = gmesh->ElementVec()[iref];
			if(!gelP2 || gelP2->HasSubElement()) continue;
			TPZRefPatternTools::RefineDirectional(gelP2, SETmatPointRefDir2);
		}		
	}

//	Need test this
//	int ndirectdivL = 2;	
//	set<int> SETmatRefDir12;
//	for(int j = 0; j < ndirectdivL; j++)
//	{
//		int nel = gmesh->NElements();
//		for (int iref = 0; iref < nel; iref++)
//		{
//			TPZVec<TPZGeoEl*> filhos;
//			TPZGeoEl * gelP2 = gmesh->ElementVec()[iref];
//			if(!gelP2 || gelP2->HasSubElement()) continue;
//			SETmatRefDir12.insert(bcLeft);
//			TPZRefPatternTools::RefineDirectional(gelP2, SETmatRefDir12);
//		}		
//	}


	
//	First computational mesh
	TPZCompMesh * cmesh1 = MalhaCompElast(gmesh, p+1);
//	Print First computational mesh
	ofstream arg2("cmesh1.txt");
	cmesh1->Print(arg2);
	
//	Second computational mesh
	TPZCompMesh * cmesh2= MalhaCompPressao(gmesh, p);
//	Print Second computational mesh
	ofstream arg3("cmesh2.txt");
	cmesh2->Print(arg3);	
	
//	Clear reference of the geometric mesh to cmesh1
	gmesh->ResetReference();
	cmesh1->LoadReferences();
	
//	Using Uniform Refinement
	RefinUniformElemComp(cmesh1,1);
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();
	
	ofstream arg4("cmesh12.txt");
	cmesh1->Print(arg4);
	ofstream arg5("gmesh2.txt");
	gmesh->Print(arg5);
	ofstream file3("malhageo1.vtk");
	PrintGMeshVTK(gmesh, file3);
	
	// Clear reference of the geometric mesh to cmesh2
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	
//	Using Uniform Refinement
	RefinUniformElemComp(cmesh2,1);
	cmesh2->AdjustBoundaryElements();
	cmesh2->CleanUpUnconnectedNodes();
	cmesh2->ExpandSolution();
	
	ofstream arg6("cmesh22.txt");
	cmesh2->Print(arg6);
	ofstream arg7("gmesh3.txt");
	gmesh->Print(arg7);
	ofstream file5("malhageo2.vtk");
	PrintGMeshVTK(gmesh, file5);
	
	
//	Solving First Computational problem
//	TPZAnalysis an1(cmesh1);
//	SolveSist(an1, cmesh1);
//	std::string plotfile("saidaSolution_cmesh1.vtk");
//	PosProcess2(an1, plotfile);
	
//	Solving Second Computational problem
	TPZAnalysis an2(cmesh2);
	SolveSist(an2, cmesh2);
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
	
//	Set initial conditions for pressure 
	int nrs = an2.Solution().Rows();
	TPZFMatrix<REAL> solucao1(nrs,1,0.0);
	cmesh2->Solution() = solucao1;
	
	TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh1;
	meshvec[1] = cmesh2;
	TPZVec <TPZPoroElastic2d *>  materialist(1,0) ;
	TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,materialist);
	
	
	// time control
	REAL hour = 3600;
	REAL day = 86400;
	REAL month = 30*day;
	REAL year = 365*day;
	REAL delta = 81.999999999*day;
	REAL MaxTime = 164*day;
	
	// Improve this
	materialist[0]->SetTimeStep(delta);		
	materialist[0]->SetTimeValue(0.0);
	//Criando matriz K2
	materialist[0]->SetLastState();
	TPZAnalysis an(mphysics);
	
	//Setting initial coditions 
	TPZFMatrix<REAL> Initialsolution = an.Solution();
	std::string output;
	output = "TransientSolution";
	std::stringstream outputfiletemp;
	outputfiletemp << output << ".vtk";
	std::string plotfile = outputfiletemp.str();
	//Print Initial conditions
	PosProcessMultphysics(meshvec,mphysics,an,plotfile);	
	
	TPZSpStructMatrix matsp(mphysics);
	//TPZSkylineStructMatrix matsp(mphysics);	
	std::set< int > materialid;
	int matid = 1;
	materialid.insert(matid);
	matsp.SetMaterialIds (materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<REAL> Un;
	TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
	
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
//		std::stringstream sout;
//		matK2->Print("K2 = ", sout,EMathematicaInput);
//		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
	
	//Criando matriz K1
	materialist[0]->SetCurrentState();
	TPZSkylineStructMatrix matsk(mphysics);
	TPZFMatrix<REAL> matK1;	
	TPZFMatrix<REAL> fvec; //vetor de carga
	an.SetStructuralMatrix(matsk);//	Set the structtural sky line matrix to our analysis
	
	
	TPZStepSolver<REAL> step; //Create Solver object
	step.SetDirect(ELDLt); //	Symmetric case
	//step.SetDirect(ELU);
	an.SetSolver(step); //	Set solver
	an.Run(); //	Excecute an analysis to obtain the Rhs vector (In this case we start with zero initial values)
	
	matK1 = an.StructMatrix(); //Storage the global matrix and load vector
	fvec = an.Rhs();
	
	// This code identify singular blocks
	TPZStepSolver<REAL> & temp = dynamic_cast<TPZStepSolver<REAL> &> (an.Solver());
	std::list <int> & zeropivot = temp.Singular(); 
	if (zeropivot.size()) 
	{
		int eq = * zeropivot.begin();
		an.Rhs().Zero();
		an.Rhs()(eq,0) = -10000.0;
		an.Solve();
		TPZFMatrix<REAL> TempSolution = an.Solution();
		
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
		PosProcessMultphysics(meshvec,mphysics,an,plotfile);
		
	}
	
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		// base system to invert
		// just one for checking purpose
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
	}
#endif
	
	///start transient problem
	SolveSistTransient(delta,MaxTime,matK1, matK2, fvec, Initialsolution, materialist, an, meshvec,  mphysics);

	
	return EXIT_SUCCESS;
}


void SolucaoExata1D(TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
	//REAL x = ptx[0];
	REAL x = ptx[1];
	
	REAL pini = 1000.;
	REAL lamb = 8333.33;
	REAL mi = 12500.;
	REAL visc =0.001; 
	REAL perm =  1.e-10;
	REAL H=1.;
	REAL tp = 10.;
	int in;
	REAL pD, uD, sigD;
	REAL PI = atan(1.)*4.;
	
	sol[0]=0.;
	sol[1]=0.;
	sol[2]=0.;
	
	REAL tD = (lamb+2.*mi)*perm*tp/(visc*H);
	REAL xD = fabs(x-1.)/H;
	for (in =0; in<1000; in++) {
		
		REAL M = PI*(2.*in+1.)/2.;
		pD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
		uD += (2./(M*M))*cos(M*xD)*exp(-1.*M*M*tD);
		sigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
	}
	
	sol[0] = pD*pini;
	sol[1] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
	sol[2] = (-1.+ sigD)*pini;
}


void SolucaoExata2D(TPZVec<REAL> &ptx,REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
	//REAL x = ptx[0];
	REAL x = ptx[0];//-25000;
	REAL y = ptx[1];//-25000;
	REAL z = 0.0;
	REAL r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
	if ( r == 0.0) {
		x=0.0001;
		y=0.0001;
		r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
	}
	
	
	//	Parameters
	REAL lamb = 8.3e9;					//	[Pa]
	REAL lambu = 13.4e9;						//	[Pa]
	REAL alphaMod=0.7;						//	[-]
	REAL G= 5.5e9;							//	[Pa]
	REAL rhof = 1000.0;						//	[kg/m3]
	REAL Pressure = 0.0;					//	[Pa]
	REAL segtime = 0.0;					//	[s]
	REAL qMod = 20.0;						//	[kg/s]
	REAL cMod = 0.082;						//	[m2/s]
	REAL PI = atan(1.)*4.;
	
	segtime = timestep;
	
	if (segtime == 0.0) {
		segtime = 1.0e-12;
	}
	
	sol[0]=0.;
	sol[1]=0.;
	sol[2]=0.;
	
#ifdef USING_BOOST
	REAL ERFCC = boost::math::erfc((0.5)*(r/sqrt(cMod*segtime)));
#else
	REAL ERFCC = erfc((0.5)*(r/sqrt(cMod*segtime)));
#endif
	
	Pressure = (qMod/(4.0*PI*rhof*cMod*r))*(((lambu-lamb)*(lamb+2.0*G))/(pow(alphaMod,2.0)*(lambu+2.0*G)))*ERFCC;
	
	sol[0] = Pressure;
//	sol[1] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
//	sol[2] = (-1.+ sigD)*pini;
}


void SolucaoExata2DLinesource(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
	
	// Defition of varibles
	REAL x = ptx[0];//-50000.0;
	REAL y = ptx[1];//-50000.0;
	REAL z = 0.0;
	REAL r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
	if ( r == 0.0) {
		x=1.0e-20;
		y=1.0e-20;
		r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
	}
	
	
	// Definitions of Parameters
	REAL lamb = 8.3e9;						//	[Pa]
	REAL lambu = 13.4e9;						//	[Pa]
	REAL alpha = 0.7;						//	[-]
	REAL G= 5.5e9;							//	[Pa]
	REAL rhof = 1000.0;						//	[kg/m3]	
	REAL Ssig = ((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
	REAL segtime = 0.0;						//	[s]
	REAL qMod = -0.001;						//	[kg/s]
	REAL cMod = 0.083;						//	[m2/s]
	REAL kdif = cMod*Ssig;
	REAL PI = atan(1.)*4.;
	segtime = timestep;
	
	if (segtime == 0.0) {
		segtime = 1.0e-20;
	}
	
	
	sol[0]=0.; // Pressure
	sol[1]=0.; // Ux
	sol[2]=0.; // Uy
	flux(0,0)=0.0; // SigX
	flux(1,0)=0.0; // SigY
	flux(2,0)=0.0; // TauXY
	
	
	REAL Pressure = 0.0;					//	[Pa]
	REAL Ux = 0.0;							//	[m]
	REAL Uy = 0.0;							//	[m]
	REAL Sigx = 0.0;						//	[Pa]
	REAL Sigy = 0.0;						//	[Pa]
	REAL Tauxy = 0.0;						//	[Pa]
	
	REAL Zz = (pow(r, 2)/(4.0*cMod*segtime));
	REAL Ei = Exponential_Integral_Ei(-Zz);
	REAL Den = (8*PI*rhof*kdif*(lamb+2.0*G));
	REAL Nem = qMod*alpha*x;

	Pressure = (qMod/(4*PI*rhof*kdif))*Exponential_Integral_Ei(-Zz);
	Ux = ((-qMod*alpha*x)/(8*PI*rhof*kdif*(lamb+2.0*G)))*(((1/Zz)*(1-exp(-Zz)))-Exponential_Integral_Ei(-Zz));
	Uy = ((-qMod*alpha*y)/(8*PI*rhof*kdif*(lamb+2.0*G)))*(((1/Zz)*(1-exp(-Zz)))-Exponential_Integral_Ei(-Zz));
	Sigx = (qMod*alpha*G/(4*PI*rhof*kdif*(2.0*G+lamb)))*(((1/Zz)*(1-exp(-Zz))*(1-(2*pow(x,2)/pow(r,2))))+Exponential_Integral_Ei(-Zz));
	Sigy = (qMod*alpha*G/(4*PI*rhof*kdif*(2.0*G+lamb)))*(((1/Zz)*(1-exp(-Zz))*(1-(2*pow(y,2)/pow(r,2))))+Exponential_Integral_Ei(-Zz));
	Tauxy = (2.0*qMod*alpha*G*x*y/(4*PI*rhof*kdif*(2.0*G+lamb)*pow(r,2)))*((1/Zz)*(1-exp(-Zz)));
	
	sol[0] = Pressure;
	sol[1] = Ux;
	sol[2] = Uy;
	
	flux(0,0)=Sigx; // SigX
	flux(1,0)=Sigy; // SigY
	flux(2,0)=Exponential_Integral_Ei(-Zz); // TauXY
}


// Right handside term of our Linear PDE
void ForcingTimeDependFunction(TPZVec<REAL> &ptx, REAL TimeValue,int WhichStateVariable,double &StateVariable) 
{
	
	// Define the relations for each variable in the right hand side of the StateVariable at the current PDE.
	
	REAL x = ptx[0];
	REAL y = ptx[1];
	REAL z = 0.0;
	
	REAL hour = 3600;
	REAL day = 86400;
	REAL month = 30*day;
	REAL year = 365*day;
	REAL delta = 99.9999*hour;
	REAL MaxTime = 100.0*hour;
	
	
	switch (WhichStateVariable) 
	{
		case 0:
		{
			//	Ux
			StateVariable = 0.0;
			break;
		}
		case 1:
		{
			//	Uy
			StateVariable = 0.0;
			break;			
		}
		case 2:
		{
			//	Pressure
//			if ((abs(x-347.15922101486848) < 1.0e-4) && (abs(y-347.15922101486848) < 1.0e-4)) 
//			{
//			StateVariable = -0.25*5.0e-5;
//			}
//			else 
//			{
				StateVariable = 0.0;
//			}
			break;
		}
		default:
			break;
	}
}

void SolucaoExataRosa1D(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
	//REAL x = ptx[0];
	REAL x = ptx[0];//-12500;
	
	
	//	Parameters
	REAL Eyoung = 1.43e10;					//	[Pa]
	REAL poisson = 0.3;						//	[-]
	REAL alpha=0.0;///((1.19667e10)*3);//1.70667e10		1.19667e10					//	[-]
	
	REAL Phi= 1.0;//9.60784e-11;//1.35695e-10				//	[-]
	REAL Ct = 1.0;					//	[kg/m3]
	REAL Se = Ct*Phi;							//	[m2/s]	
	REAL Visc = 1.0;						//	[Pa.s]
	REAL Kmed = 1.0;//7.8784288e-15;//1.109542e-14						//	[m2]
	REAL Qo = 2.0;
	REAL Bo = 1.0;
	REAL PI = atan(1.)*4.;
	REAL segtime = 0.0;					//	[s]	

	
//	REAL Phi= 0.18;//9.60784e-11;//1.35695e-10				//	[-]
//	REAL Ct = (150.0e-6)*(1/(98066.50));					//	[kg/m3]
//	REAL Se = Ct*Phi;							//	[m2/s]	
//	REAL Visc = 0.8*(1.e-3);						//	[Pa.s]
//	REAL Kmed = 20*(9.86923e-16);//7.8784288e-15;//1.109542e-14						//	[m2]
//	REAL Qo = 400.0/(86400);
//	REAL Bo = 1.2;
//	REAL PI = atan(1.)*4.;
//	REAL segtime = 0.0;					//	[s]		
	
	segtime = timestep;
	
	if (segtime == 0.0) {
		segtime = 1.0e-12;
	}
	
	x = abs(x);
	sol[0]=0.;
	sol[1]=0.;
	sol[2]=0.;
	
	REAL Eta = (Kmed)/(Phi*Visc*Ct);
	REAL Pressure = ((0.5*Qo*Bo*Visc)/(Kmed*1.0))*(sqrt((4*Eta*segtime)/PI)*(exp(-(pow(x,2)/(4*Eta*segtime))))-(x*erfc(x/sqrt(4*Eta*segtime))));
	
	sol[0] = Pressure;
	//	sol[1] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
	//	sol[2] = (-1.+ sigD)*pini;
}

// Analytical solution for flamant problem Elasticity Theory: Applications and Numerics 
void ExactSolutionFlamantProblem(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux)
{
	
	// Defition of variables
	REAL x = ptx[0]-50000.0;
	REAL y = ptx[1]-50000.0;
	REAL z = 0.0;
	REAL Yforce = 1000.0;
	REAL PI = atan(1.)*4.;	
	REAL r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
	if ( r == 0.0) {
		x=1.0e-10;
		y=1.0e-10;
		r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
	}
	
	
	// Definitions of Parameters
//	REAL lamb = 1.0e9;						//	[Pa]
//	REAL G= 1.0e9;							//	[Pa]
//	REAL rhof = 1.0;						//	[kg/m3]	
//	REAL qMod = -1.0;						//	[kg/s]
//	REAL cMod = 1.0;						//	[m2/s]
//	REAL kdif = cMod*Se;
//	REAL PI = atan(1.)*4.;
	REAL sigXX = -2.0*((Yforce*pow(x,2)*y)/(PI*(pow(r,2))));
	REAL sigYY = -2.0*((Yforce*pow(y,3))/(PI*(pow(r,2))));
	REAL tauXY = -2.0*((Yforce*pow(y,2)*x)/(PI*(pow(r,2))));

	sol[0] = sigXX;
	sol[1] = sigYY;
	sol[2] = tauXY;
}


void SolucaoExataRosa1DPseudo(TPZVec<REAL> &ptx, REAL timestep, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
	//REAL x = ptx[0];
	REAL x = ptx[0];//-12500;
	REAL L = 20000.0;
//	x = abs(x);
	sol[0]=0.;
	sol[1]=0.;
	sol[2]=0.;

	REAL Pressure = 1000 + L*((x/L)-0.5*(pow((x/L),2)));
	
	sol[0] = Pressure;
	//	sol[1] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
	//	sol[2] = (-1.+ sigD)*pini;
}


TPZCompMesh*MalhaCompPressao(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	TPZMaterial * mat(material);
	
	REAL diff = 0.1;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
	REAL flux = 0.;
	
	material->SetParameters(diff, conv, convdir);
	material->SetInternalFlux( flux);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	
	
	///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	REAL pN=0.;
	val2(0,0)=pN;
	TPZMaterial * BCondNL = material->CreateBC(mat, bcBottom,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondNL);
	
	TPZMaterial * BCondNU = material->CreateBC(mat, bcTop,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondNU);
	
	TPZFMatrix<REAL> val12(2,2,0.), val22(2,1,0.);
	REAL uDL=4000.;
	val22(0,0)=uDL;
	TPZMaterial * BCondDL = material->CreateBC(mat, bcLeft,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDL);
	
	REAL uDR=3000.;
	val22(0,0)=uDR;
	TPZMaterial * BCondDR = material->CreateBC(mat, bcRight,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDR);
	
	// Point source
	TPZMaterial * BCondNpoint = material->CreateBC(mat, pointsource,neumann, val12, val22);
	cmesh->InsertMaterialObject(BCondNpoint);
	
	TPZMaterial * BCyfixedPoints = material->CreateBC(mat, yfixedPoints,neumann, val12, val22);
	cmesh->InsertMaterialObject(BCyfixedPoints);	
	
	TPZMaterial * BCxfixedPoints = material->CreateBC(mat, xfixedPoints,neumann, val12, val22);
	cmesh->InsertMaterialObject(BCxfixedPoints);	

	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}

TPZCompMesh*MalhaCompElast(TPZGeoMesh * gmesh,int pOrder)
{
	REAL rockrho = 2330.0; // SI system
	REAL gravity = 9.8; // SI system
	REAL overburdendepth = 2000.0; // SI system
	REAL layerthickness = 10.0;  // SI system
	
	REAL E = 100.;
	REAL poisson = 0.35;
	int dim = 2;
	TPZVec<REAL> force(2,0.);
	force[1]=gravity*rockrho;
	//force[0] = 0.;
	TPZElasticityMaterial *material;
	int planestress = -1;
	material = new TPZElasticityMaterial(matId, E, poisson, force[0], force[1], planestress); 
	
	
	TPZMaterial * mat(material);
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	
	
	///Inserir condicao de contorno
	
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	REAL uNUy=rockrho*gravity*overburdendepth;
	val2(1,0)=uNUy;
	
	TPZMaterial * BCondNU = material->CreateBC(mat, bcTop,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondNU);
	
	REAL uNLy=rockrho*gravity*(overburdendepth+layerthickness);
	val2(1,0)=uNLy;	
	
	TPZMaterial * BCondNL = material->CreateBC(mat, bcBottom,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondNL);
	
	
	TPZFMatrix<REAL> val12(2,2,0.), val22(2,1,0.);
	REAL uDL=0.0;
	val22(0,0)=uDL;
	TPZMaterial * BCondDL = material->CreateBC(mat, bcLeft,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDL);
	
	REAL uDR=0.0;
	val22(0,0)=uDR;
	TPZMaterial * BCondDR = material->CreateBC(mat, bcRight,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDR);
	
	
	// Point source
	TPZMaterial * BCondNpoint = material->CreateBC(mat, pointsource,neumann, val12, val22);
	cmesh->InsertMaterialObject(BCondNpoint);
	
	TPZMaterial * BCyfixedPoints = material->CreateBC(mat, yfixedPoints,neumann, val12, val22);
	cmesh->InsertMaterialObject(BCyfixedPoints);
	
	TPZMaterial * BCxfixedPoints = material->CreateBC(mat, xfixedPoints,neumann, val12, val22);
	cmesh->InsertMaterialObject(BCxfixedPoints);	
	
//	//Ajuste da estrutura de dados computacional
//	std::set<int> set1;
//	set1.insert(1);
//	cmesh->AutoBuild(set1);
//	set1.clear();
//	set1.insert(2);
//	cmesh->Reference()->ResetReference();
//	cmesh->AutoBuild(set1);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();	
	
	return cmesh;
	
}

TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial)
{
	//Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
	mphysics->SetAllCreateFunctionsMultiphysicElem();
	
	
	
	int MatId = 1;
	int dim = 2;
//	Parameters
//	
//	// Definitions
//	REAL lamb = 8.3e9;					//	[Pa]
//	REAL lambu = 13.4e9;						//	[Pa]
//	REAL alpha = 0.7;						//	[-]
//	REAL G= 5.5e9;							//	[Pa]
//	REAL rhof = 1000.0;						//	[kg/m3]
//	REAL poisson = 0.3;						//	[-]
//	REAL Eyoung = 2*G*(1+poisson);					//	[Pa]	
//	REAL Se = ((pow(alpha,2))/(lambu-lamb));//1.35695e-10				//	[-]
//	REAL rockrho = 0.0;					//	[kg/m3]
//	REAL gravity = 0.0;					//	[m/s2]
//	REAL c = 0.082;							//	[m2/s]	
//	REAL visc = 1.0;						//	[Pa.s]
//	REAL perm = c*Se;//7.8784288e-15;//1.109542e-14						//	[m2]
//	REAL qo = -20.0;
//	REAL Bo = 1.0;
//	REAL PI = atan(1.)*4.;	

	// Definitions
	REAL lamb = 8.3e9;					//	[Pa]
	REAL lambu = 13.4e9;						//	[Pa]
	REAL alpha = 0.7;						//	[-]
	REAL G= 5.5e9;							//	[Pa]
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
	
	
//	// Definitions Flamant problem
//	REAL lamb = 1.0e9;					//	[Pa]
//	REAL alpha = 0.0;						//	[-]
//	REAL G= 1.0e9;							//	[Pa]
//	REAL rhof = 1.0;						//	[kg/m3]
//	REAL poisson3D = 0.3;						//	[-]
//	REAL Eyoung3D = (G*(3.0*lamb+2.0*G))/(lamb+G);					//	[Pa]	
//	REAL poisson = (poisson3D)/(1-poisson3D);						//	[-]
//	REAL Eyoung = (Eyoung3D)/(1-pow(poisson3D,2.0));					//	[Pa]		
//	REAL Se = 1.0;//((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
//	REAL rockrho = 2500.0;					//	[kg/m3]
//	REAL gravity = -9.81;					//	[m/s2]
//	REAL c = 1.0;							//	[m2/s]	
//	REAL visc = 1.0;						//	[Pa.s]
//	REAL perm = c*Se;//7.8784288e-15;//1.109542e-14						//	[m2]
//	REAL qo = 1.0;
	REAL Yforce = -1000.0;
//	REAL Bo = 1.0;
//	REAL PI = atan(1.)*4.;		
	
	REAL fx=0.0;
	REAL fy=0.0;//-rockrho*gravity;	
	
	
	int planestress = 0; // This is a Plain strain problem
	mymaterial[0] = new TPZPoroElastic2d (MatId, dim);
	mymaterial[0]->SetParameters(Eyoung, poisson, fx, fy);
	mymaterial[0]->SetParameters(perm,visc);
	mymaterial[0]->SetfPlaneProblem(planestress);
	mymaterial[0]->SetBiotParameters(alpha,Se);
	
	TPZAutoPointer<TPZFunction<STATE> > TimeDepForcingF;
	TPZAutoPointer<TPZFunction<STATE> > TimeDepFExact;
	
//	TimeDepForcingF = new TPZDummyFunction<STATE>(ForcingTimeDependFunction);
	TimeDepFExact = new TPZDummyFunction<STATE>(SolucaoExata2DLinesource);
//	mymaterial->SetTimeDependentForcingFunction(TimeDepForcingF);
	mymaterial[0]->SetTimeDependentFunctionExact(TimeDepFExact);

	
	ofstream argm("mymaterial1.txt");
	mymaterial[0]->Print(argm);
	TPZMaterial * mat(mymaterial[0]);
	mphysics->InsertMaterialObject(mat);
	
	///--- --- Inserir condicoes de contorno
	// Radial  Model (1)
	// Squared Model (2) for consistency
	// Squared Model (3) for validation
	// SemiCircular HalfSpace Model (4) for validation	
	// Graven test (5)		
	
	int model = 1;
	
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
//			TPZMaterial * BonArc1 = mymaterial[0]->CreateBC(mat, arc1,D1D2, val1, val2);
//			mphysics->InsertMaterialObject(BonArc1);
//			
//			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
//			TPZMaterial * BonArc2 = mymaterial[0]->CreateBC(mat, arc2,D1D2, val1, val2);
//			mphysics->InsertMaterialObject(BonArc2);
//			
//			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
//			TPZMaterial * BonArc3 = mymaterial[0]->CreateBC(mat, arc3,D1D2, val1, val2);
//			mphysics->InsertMaterialObject(BonArc3);
//			
//			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
//			TPZMaterial * BonArc4 = mymaterial[0]->CreateBC(mat, arc4,D1D2, val1, val2);
//			mphysics->InsertMaterialObject(BonArc4);
						
			
//		Production/Injection Well Just Neumman for Pressure
		int WellQ = 400;
		TPZFMatrix<REAL> valP13(3,2,0.), valP23(3,1,0.);
			REAL MassQ = qo/rhof;
		valP23(2,0)=MassQ;
		TPZMaterial * BCPoint = mymaterial[0]->CreateBC(mat, WellPoint, WellQ, valP13, valP23);
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
			
			TPZMaterial * BCondT = mymaterial[0]->CreateBC(mat, bcTop,D1D2, val12, val22);
			mphysics->InsertMaterialObject(BCondT);	

			TPZMaterial * BCondBt = mymaterial[0]->CreateBC(mat, bcBottom,D1D2, val12, val22);
			mphysics->InsertMaterialObject(BCondBt);
			
			TPZMaterial * BCondL = mymaterial[0]->CreateBC(mat, bcLeft, D1D2, val12, val22);
			mphysics->InsertMaterialObject(BCondL);

			TPZMaterial * BCondR = mymaterial[0]->CreateBC(mat, bcRight, D1D2, val12, val22);
			mphysics->InsertMaterialObject(BCondR);
			
//		Production/Injection Well Just Neumman for Pressure
		int WellQ = 400;
		TPZFMatrix<REAL> valPP12(3,2,0.), valPP22(3,1,0.);
		REAL MassQ=10.0;
		valPP22(2,0)=MassQ;
		TPZMaterial * BCPoint = mymaterial[0]->CreateBC(mat, pointsource, WellQ, valPP12, valPP22);
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
			
			TPZMaterial * BCondT = mymaterial[0]->CreateBC(mat, bcTop,N1N2, valall12, valall22);
			mphysics->InsertMaterialObject(BCondT);
			TPZMaterial * BCondB = mymaterial[0]->CreateBC(mat, bcBottom,N1N2, valall12, valall22);
			mphysics->InsertMaterialObject(BCondB);
			TPZMaterial * BCondL = mymaterial[0]->CreateBC(mat, bcLeft,N1N2, valall12, valall22);
			mphysics->InsertMaterialObject(BCondL);
			TPZMaterial * BCondR = mymaterial[0]->CreateBC(mat, bcRight,N1N2, valall12, valall22);
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
			TPZMaterial * BCyfixedPoints = mymaterial[0]->CreateBC(mat, yfixedPoints, DYFX1N2, val11, val21);
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
			TPZMaterial * BCxfixedPoints = mymaterial[0]->CreateBC(mat, xfixedPoints, DXFY1N2, val13, val23);
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
//			TPZMaterial * BCondBt = mymaterial[0]->CreateBC(mat, bcBottom, DYFX1N2, val12, val22);
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
//			TPZMaterial * BCondL = mymaterial[0]->CreateBC(mat, bcLeft, DXFY1N2, val13, val23);
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
//			TPZMaterial * BCondR = mymaterial[0]->CreateBC(mat, bcRight, DXFY1N2, val14, val24);
//			mphysics->InsertMaterialObject(BCondR);
			
			//	Production or injection Well point (Q is the injected/produced volume per time [kg/s])
			int Well = 400;
			TPZFMatrix<REAL> valP13(3,2,0.), valP23(3,1,0.);
			valP23(2,0)=qo/rhof;
			TPZMaterial * BCPoint = mymaterial[0]->CreateBC(mat, WellPoint, Well, valP13, valP23);
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
			TPZMaterial * BonArc1 = mymaterial[0]->CreateBC(mat, bcLeft,D1D2, val1, val2);
			mphysics->InsertMaterialObject(BonArc1);
			
			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
			TPZMaterial * BonArc2 = mymaterial[0]->CreateBC(mat, bcRight,D1D2, val1, val2);
			mphysics->InsertMaterialObject(BonArc2);
			
			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
			TPZMaterial * BonArc3 = mymaterial[0]->CreateBC(mat, bcBottom,D1D2, val1, val2);
			mphysics->InsertMaterialObject(BonArc3);
			
			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
			TPZMaterial * BonArc4 = mymaterial[0]->CreateBC(mat, bcTop,N1D2, val1, val2);
			mphysics->InsertMaterialObject(BonArc4);			
			
			//		loaded vertical point force Y
			int N1N2 = 11;
			TPZFMatrix<REAL> valP13(3,2,0.), valP23(3,1,0.);
			// 170 Elasticity Theory
			valP23(1,0)=Yforce;
			TPZMaterial * BCPoint = mymaterial[0]->CreateBC(mat, WellPoint, N1N2, valP13, valP23);
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
			TPZMaterial * BonArc1 = mymaterial[0]->CreateBC(mat, GBTop,D1D2, val1, val2);
			mphysics->InsertMaterialObject(BonArc1);
			
			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
			TPZMaterial * BonArc2 = mymaterial[0]->CreateBC(mat, GBBot,N1N2, val1, val2);
			mphysics->InsertMaterialObject(BonArc2);
			
			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
			TPZMaterial * BonArc3 = mymaterial[0]->CreateBC(mat, GBleft,D1FN2, val1, val2);
			mphysics->InsertMaterialObject(BonArc3);
			
			///Inserir cc de Dirichlet para elasticidade e Neumann para pressao
			TPZMaterial * BonArc4 = mymaterial[0]->CreateBC(mat, GBright,D1FN2, val1, val2);
			mphysics->InsertMaterialObject(BonArc4);			
			
//			//		Production/Injection Well Just Neumman for Pressure
//			int WellQ = 400;
//			TPZFMatrix<REAL> valP13(3,2,0.), valP23(3,1,0.);
//			REAL MassQ = qo/rhof;
//			valP23(2,0)=MassQ;
//			TPZMaterial * BCPoint = mymaterial[0]->CreateBC(mat, WellPoint, WellQ, valP13, valP23);
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
		mymaterial[0]->SetTimeValue(TimeValue);
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
	int div = 2;
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
			if (gel->Dimension() == 2) gel->Divide (filhos);
		}//for i
	}//ref
}

void RefinamentoUniforme(TPZGeoMesh * gMesh, int nh, int MatId, int indexEl){
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gMesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2){
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


////////////////////////////////////////////////////////////////////////////////
// File: exponential_integral_Ei.c                                            //
// Routine(s):                                                                //
//    Exponential_Integral_Ei                                                 //
//    xExponential_Integral_Ei                                                //
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_Ei( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Exponential_Integral_Ei( x );                                      //
////////////////////////////////////////////////////////////////////////////////
double Exponential_Integral_Ei( double x )
{
	return (double) xExponential_Integral_Ei( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xExponential_Integral_Ei( long double x )                      //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the exponential integral Ei().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xExponential_Integral_Ei( x );                                     //
////////////////////////////////////////////////////////////////////////////////

long double xExponential_Integral_Ei( long double x )
{
	if ( x < -5.0L ) return Continued_Fraction_Ei(x);
	if ( x == 0.0L ) return -DBL_MAX;
	if ( x < 6.8L )  return Power_Series_Ei(x);
	if ( x < 50.0L ) return Argument_Addition_Series_Ei(x);
	return Continued_Fraction_Ei(x);
}

////////////////////////////////////////////////////////////////////////////////
// static long double Continued_Fraction_Ei( long double x )                  //
//                                                                            //
//  Description:                                                              //
//     For x < -5 or x > 50, the continued fraction representation of Ei      //
//     converges fairly rapidly.                                              //
//                                                                            //
//     The continued fraction expansion of Ei(x) is:                          //
//        Ei(x) = -exp(x) { 1/(-x+1-) 1/(-x+3-) 4/(-x+5-) 9/(-x+7-) ... }.    //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Continued_Fraction_Ei( long double x )
{
	long double Am1 = 1.0L;
	long double A0 = 0.0L;
	long double Bm1 = 0.0L;
	long double B0 = 1.0L;
	long double a = expl(x);
	long double b = -x + 1.0L;
	long double Ap1 = b * A0 + a * Am1;
	long double Bp1 = b * B0 + a * Bm1;
	int j = 1;
	
	a = 1.0L;
	while ( fabsl(Ap1 * B0 - A0 * Bp1) > epsilon * fabsl(A0 * Bp1) ) {
		if ( fabsl(Bp1) > 1.0L) {
			Am1 = A0 / Bp1;
			A0 = Ap1 / Bp1;
			Bm1 = B0 / Bp1;
			B0 = 1.0L;
		} else {
			Am1 = A0;
			A0 = Ap1;
			Bm1 = B0;
			B0 = Bp1;
		}
		a = -j * j;
		b += 2.0L;
		Ap1 = b * A0 + a * Am1;
		Bp1 = b * B0 + a * Bm1;
		j += 1;
	}
	return (-Ap1 / Bp1);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Power_Series_Ei( long double x )                        //
//                                                                            //
//  Description:                                                              //
//     For -5 < x < 6.8, the power series representation for                  //
//     (Ei(x) - gamma - ln|x|)/exp(x) is used, where gamma is Euler's gamma   //
//     constant.                                                              //
//     Note that for x = 0.0, Ei is -inf.  In which case -DBL_MAX is          //
//     returned.                                                              //
//                                                                            //
//     The power series expansion of (Ei(x) - gamma - ln|x|) / exp(x) is      //
//        - Sum(1 + 1/2 + ... + 1/j) (-x)^j / j!, where the Sum extends       //
//        from j = 1 to inf.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Power_Series_Ei( long double x )
{ 
	long double xn = -x;
	long double Sn = -x;
	long double Sm1 = 0.0L;
	long double hsum = 1.0L;
	long double g = 0.5772156649015328606065121L;
	long double y = 1.0L;
	long double factorial = 1.0L;
	
	if ( x == 0.0L ) return (long double) -DBL_MAX;
	
	while ( fabsl(Sn - Sm1) > epsilon * fabsl(Sm1) ) {
		Sm1 = Sn;
		y += 1.0L;
		xn *= (-x);
		factorial *= y;
		hsum += (1.0 / y);
		Sn += hsum * xn / factorial;
	}
	return (g + logl(fabsl(x)) - expl(x) * Sn);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Argument_Addition_Series_Ei(long double x)              //
//                                                                            //
//  Description:                                                              //
//     For 6.8 < x < 50.0, the argument addition series is used to calculate  //
//     Ei.                                                                    //
//                                                                            //
//     The argument addition series for Ei(x) is:                             //
//     Ei(x+dx) = Ei(x) + exp(x) Sum j! [exp(j) expj(-dx) - 1] / x^(j+1),     //
//     where the Sum extends from j = 0 to inf, |x| > |dx| and expj(y) is     //
//     the exponential polynomial expj(y) = Sum y^k / k!, the Sum extending   //
//     from k = 0 to k = j.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////
static long double Argument_Addition_Series_Ei(long double x)
{
	static long double ei[] = {
		1.915047433355013959531e2L,  4.403798995348382689974e2L,
		1.037878290717089587658e3L,  2.492228976241877759138e3L,
		6.071406374098611507965e3L,  1.495953266639752885229e4L,
		3.719768849068903560439e4L,  9.319251363396537129882e4L,
		2.349558524907683035782e5L,  5.955609986708370018502e5L,
		1.516637894042516884433e6L,  3.877904330597443502996e6L,
		9.950907251046844760026e6L,  2.561565266405658882048e7L,
		6.612718635548492136250e7L,  1.711446713003636684975e8L,
		4.439663698302712208698e8L,  1.154115391849182948287e9L,
		3.005950906525548689841e9L,  7.842940991898186370453e9L,
		2.049649711988081236484e10L, 5.364511859231469415605e10L,
		1.405991957584069047340e11L, 3.689732094072741970640e11L,
		9.694555759683939661662e11L, 2.550043566357786926147e12L,
		6.714640184076497558707e12L, 1.769803724411626854310e13L,
		4.669055014466159544500e13L, 1.232852079912097685431e14L,
		3.257988998672263996790e14L, 8.616388199965786544948e14L,
		2.280446200301902595341e15L, 6.039718263611241578359e15L,
		1.600664914324504111070e16L, 4.244796092136850759368e16L,
		1.126348290166966760275e17L, 2.990444718632336675058e17L,
		7.943916035704453771510e17L, 2.111342388647824195000e18L,
		5.614329680810343111535e18L, 1.493630213112993142255e19L,
		3.975442747903744836007e19L, 1.058563689713169096306e20L
	};
	int  k = (int) (x + 0.5);
	int  j = 0;
	long double xx = (long double) k;
	long double dx = x - xx;
	long double xxj = xx;
	long double edx = expl(dx);
	long double Sm = 1.0L;
	long double Sn = (edx - 1.0L) / xxj;
	long double term = DBL_MAX;
	long double factorial = 1.0L;
	long double dxj = 1.0L;
	
	while (fabsl(term) > epsilon * fabsl(Sn) ) {
		j++;
		factorial *= (long double) j;
		xxj *= xx;
		dxj *= (-dx);
		Sm += (dxj / factorial);
		term = ( factorial * (edx * Sm - 1.0L) ) / xxj;
		Sn += term;
	}
	
	return ei[k-7] + Sn * expl(xx); 
}

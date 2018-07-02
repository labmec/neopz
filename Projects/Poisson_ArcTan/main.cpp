/**
 * @file
 * @brief This file contains the tests for validation auto-adaptive algorithms
 * @note The solution is a function (Arc Tangent) with high gradient on points with radius r=0.25 from a center C with coordinates equal to 0.5 . The gradient depends on parameter e, we are using e = 10^6 to high gradient.
 */

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"
#include "pzgeotetrahedra.h"
#include "TPZGeoCube.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"

#include "pzgeoelbc.h"

#include "pzcclonemesh.h"
#include "pzonedref.h"
#include "pzadaptmesh.h"

#include "pzlog.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzsbstrmatrix.h"
#include "pzfstrmatrix.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzpoisson3d.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzshapelinear.h"

#include "TPZRefPatternTools.h"

#include <time.h>
#include <stdio.h>

#include <fstream>
#include <cmath>

using namespace std;
using namespace pzshape;
using namespace pzgeom;

/** VARIABLES */
/** Printing level */
int gPrintLevel = 0;
int printing = 1;

int materialId = 1;
int id_bc0 = -1;
int id_bc1 = -2;
int materialBC1 = 2;
int anothertests = 0;

char saida[512];

ofstream out("OutPoissonArcTan.txt");             // To store output of the console
ofstream outLaplace("OutLaplace.txt");

int gDebug = 0;


/** Rotation data to construct deformed meshes */
bool rotating = false;
TPZFNMatrix<16,REAL> Rot(4,4,0.),RotInv(4,4,0.);
/** angle ??? */
REAL alfa = M_PI/6.;
REAL transx = 3.;
REAL transy = 0.;


/** PROBLEM WITH HIGH GRADIENT ON CIRCUNFERENCE  ---  DATA */
STATE ValueK = 100000;


/** To identify localization of PZ resources */
std::string Archivo = PZSOURCEDIR;


/** Functions to construction of geometry of problems */
TPZGeoMesh *CreateLShapeGeoMesh(MElementType typeel);
TPZGeoMesh *CreateGeoMesh(MElementType typeel);
TPZGeoMesh *CreateGeoMesh(std::string &nome);
TPZGeoMesh *CreateGeoMeshWithClassesPre(MElementType typeel);
// Crea malla computacional sem forcingfunction quando hasforcingfunction = 0, ou toma diferentes forcingfuncition para diferentes
// valores de hasforcingfunction
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);
TPZGeoMesh *ConstructingPositiveCube(REAL InitialL,MElementType typeel);
TPZGeoMesh *ConstructingTetrahedraInCube(REAL InitialL);
TPZGeoMesh *ConstructingPrismsInCube(REAL InitialL);
TPZGeoMesh *ConstructingPyramidsInCube(REAL InitialL);
TPZGeoMesh *ConstructingSeveral3DElementsInCube(REAL InitialL,MElementType typeel);

/** Fucntions to apply refinement. */
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=true, const int matidtodivided=1);
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<TPZVec<REAL> > &points,REAL &distance,bool &isdefined);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &points,REAL r,REAL &distance,bool &isdefined);
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,REAL &radius,int ntyperefs);

void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points);

void DeterminingPOrderOnLevelHRefinement(TPZCompMesh *cmesh,int p);

// Printing in VTK format the geometric mesh but taking geometric elements as reference or computational elements as reference
void PrintGeoMeshInVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);
void PrintGeoMeshAsCompMeshInVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);
void PrintGeoMeshAsCompMeshInVTKWithElementIDAsData(TPZGeoMesh *gmesh,char *filename);
void PrintGeoMeshAsCompMeshInVTKWithElementIndexAsData(TPZGeoMesh *gmesh,char *filename);
/** Print the elements of the computational mesh with associated element data */
void PrintGeoMeshAsCompMeshInVTKWithElementData(TPZGeoMesh *gmesh,char *filename,TPZVec<REAL> &elData);


/** Functions for Differential equation - Right terms, solutions and derivatives */
void RightTermCircle(const TPZVec<REAL> &x, TPZVec<STATE> &force, TPZFMatrix<STATE> &dforce);

void ExactSolCircle(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
void ExactSolLaplace(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
void ExactSolLaplaceBC(const TPZVec<REAL> &x, TPZVec<STATE> &sol);


/** Utilitaries */
void formatTimeInSec(char *strtime,int lenstrtime,int timeinsec);

void GradientReconstructionByLeastSquares(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var=0,bool continuous=false);

int DefineDimensionOverElementType(MElementType typeel);
void GetFilenameFromGID(MElementType typeel, std::string &name);
/** Detects the bigger dimension of the computational elements into cmesh to set the Model Dimension */
bool DefineModelDimension(TPZCompMesh *cmesh);

// To save meshes in disk
void SaveCompMesh(TPZCompMesh *cmesh, int timessave,TPZCompMesh *cmeshmodified=NULL,bool check=false);
// Save information of the current mesh to compare with cloned mesh (geometric mesh plus computational mesh)
//void SaveCompMesh(TPZCompCloneMesh *cmesh, int timessave,TPZCompCloneMesh *cmeshmodified=NULL,bool check=false);


/** PROBLEMS */
bool SolveSymmetricPoissonProblemOnCubeMesh();
bool SolveLaplaceProblemOnLShapeMesh();

// Generic data for problems to solve
bool usethreads = false;
int MaxPOrder = TPZOneDRef::gMaxP;

TPZManVector<REAL,3> CCircle(3,0.5);
REAL RCircle = 0.25;


// MAIN FUNCTION TO NUMERICAL SOLVE WITH AUTO ADAPTIVE HP REFINEMENTS
/** Laplace equation on square 1D 2D 3D - Volker John article 2000 */
int main() {

#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();
//    gRefDBase.InitializeRefPatterns();

    // Solving laplace problema on LShape domain in 2D.
    if(!SolveLaplaceProblemOnLShapeMesh())
        return 2;
    
    // Solving symmetricPoissonProblem on [0,1]^d with d=1, d=2 and d=3
    if(!SolveSymmetricPoissonProblemOnCubeMesh())
        return 1;
    
    return 0;
}


bool SolveSymmetricPoissonProblemOnCubeMesh() {
	// To compute processing times
	time_t sttime;
	time_t endtime;
	int time_elapsed;
	char time_formated[256];
	memset(time_formated,0,256);
	
	// Output files
    std::ofstream convergence("convergence.txt");
	std::ofstream fileerrors("ErrorsHP_ArcTan.txt");   // To store all errors calculated by TPZAnalysis (PosProcess)
	// Initial message to print computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	int nref = 1, NRefs = 12;
    int ninitialrefs = 1;
	int nthread = 2, NThreads = 4;
    int dim;
	
    //Working on regular meshes
    for(int regular=1; regular>0; regular--) {
		fileerrors << "Type of mesh: " << regular << " Level. " << endl;
		MElementType typeel;
//		for(int itypeel=(int)ECube;itypeel<(int)EPolygonal;itypeel++)
		for(int itypeel=(int)ETetraedro;itypeel<(int)EPiramide;itypeel++)
		{
			typeel = (MElementType)itypeel;
			fileerrors << "Type of element: " << typeel << endl;
			TPZGeoMesh *gmesh;
			if(!regular) {
				std::string nombre;
				GetFilenameFromGID(typeel,nombre);
				// Generating geometric mesh 
				gmesh = CreateGeoMesh(nombre);
			}
			else {
				gmesh = CreateGeoMesh(typeel);
			}
			dim = DefineDimensionOverElementType(typeel);
			
			// Defining initial refinements and total refinements depends on dimension of the model
			if(dim==3) {
                MaxPOrder = 3;
                NRefs = 3;
            }
            else if(dim==2) {
                MaxPOrder = 10;
                NRefs = 15;
            }
            else {
                NRefs = 25;
            }
			// Printing geometric mesh to validate
			if(gDebug) {
				sprintf(saida,"gmesh_%02dD_H%dTR%dE%d.vtk",dim,nref,regular,typeel);
				PrintGeoMeshInVTKWithDimensionAsData(gmesh,saida);
			}
            UniformRefinement(ninitialrefs,gmesh,dim);

			// Creating computational mesh (approximation space and materials)
			int p = 2, pinit;
			pinit = p;
			TPZCompEl::SetgOrder(p);
			TPZCompMesh *cmesh = CreateMesh(gmesh,dim,1);               // Forcing function is out 2013_07_25
			gmesh->SetName("Malha Geometrica original");
			cmesh->SetName("Malha Computacional Original");
			if(gDebug) {
				sprintf(saida,"gmesh_%02dD_H%dTR%dE%dIndex.vtk",dim,nref,regular,typeel);
				PrintGeoMeshAsCompMeshInVTKWithElementIndexAsData(gmesh,saida);
			}

			// Selecting orthogonal polynomial family to construct shape functions
			if(anothertests)
				TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;  // Setting Chebyshev polynomials as orthogonal sequence generating shape functions
			
			/** Variable names for post processing */
			TPZStack<std::string> scalnames, vecnames;
			scalnames.Push("POrder");
			scalnames.Push("Solution");

#ifdef LOG4CXX
            // Cleaning log4cxx files
            InitializePZLOG();
#endif
            
			// Solving adaptive process
			for(nref=0;nref<NRefs;nref++) {
				out << "\nConstructing Poisson problem " << dim << "D. Refinement: " << nref << " Threads: " << nthread << " Regular: " << regular << " TypeElement: " << typeel << endl;
                std::cout << "\nConstructing Poisson problem. Type element: " << typeel << std::endl;
				if(nref > 5) nthread = 2*NThreads;
				else nthread = NThreads;
				
				// Initializing the generation mesh process
				time(& sttime);
				
				// Introduzing exact solution depending on the case
				TPZAnalysis an(cmesh);
				an.SetExact(ExactSolCircle);
				{
					std::stringstream sout;
					sout << "Poisson" << dim << "D_MESH" << regular << "E" << typeel << "Thr" << nthread << "H" << std::setprecision(2) << nref << "P" << pinit << ".vtk";
					an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
				}
				std::string MeshFileName;
				{
					std::stringstream sout;
					sout << "meshAngle" << dim << "D_MESH" << regular << "E" << typeel << "Thr" << nthread << "H" << std::setprecision(2) << nref << "P" << pinit << ".vtk";
					MeshFileName = sout.str();
				}
				
				cmesh->SetName("Malha computacional adaptada");
				// Printing geometric and computational mesh
				if(gDebug) {
					cmesh->Reference()->Print(std::cout);
					cmesh->Print(std::cout);
				}
				
				// Solve using symmetric matrix then using Cholesky (direct method)
				TPZSkylineStructMatrix strskyl(cmesh);
                if(usethreads) strskyl.SetNumThreads(nthread);
				an.SetStructuralMatrix(strskyl);
				
				TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
				direct->SetDirect(ECholesky);
				an.SetSolver(*direct);
				delete direct;
				direct = 0;
				
				an.Run();
				
				// Post processing
				an.PostProcess(0,dim);
				if(gDebug) {
					std::ofstream out(MeshFileName.c_str());
					cmesh->LoadReferences();
					TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(),out,false);
				}
                
				// generation mesh process finished
				time(&endtime);
				time_elapsed = endtime - sttime;
				formatTimeInSec(time_formated,256,time_elapsed);
				out << "\tRefinement: " << nref+1 << " Regular Mesh: " << regular << " TypeElement: " << typeel << " Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n\n";
				
				// Initializing the auto adaptive process
				REAL valerror =0.;
				REAL valtruerror=0.;
				TPZVec<REAL> ervec,truervec,effect;
				
				TPZAdaptMesh adapt(MaxPOrder);
				adapt.SetCompMesh(cmesh);
				
				out << "\n\nEntering Auto Adaptive Methods... step " << nref << "\n";
                std::cout << "\n\nEntering Auto Adaptive Methods... step " << nref << "\n";
				fileerrors << "\n\nEntering Auto Adaptive Methods... step " << nref << "\n";

				TPZCompMesh *adaptmesh = NULL;
				if(NRefs>1) {
					time(&sttime);
					adaptmesh = adapt.GetAdaptedMesh(valerror,valtruerror,ervec,ExactSolCircle,truervec,effect,fileerrors,0,typeel,printing);
                    // Saving computational mesh before adaptive process
                    if(!adaptmesh)
                        return false;

					time_t endtime;
					time(&endtime);
					int time_elapsed = endtime - sttime;
					out << "\n\nExiting Auto Adaptive Methods....step " << nref << "time elapsed " << time_elapsed << "\n\n\n\n";
					fileerrors << "\n\nExiting Auto Adaptive Methods....step " << nref << "time elapsed " << time_elapsed << "\n\n\n\n";
					
					int prt;
					out << "neq = " << cmesh->NEquations() << " error estimate = " << valerror << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
					fileerrors << "neq = " << cmesh->NEquations() << " error estimate = " << valerror << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;

					convergence  << cmesh->NEquations() << "\t" << valerror << "\t" << valtruerror << "\t" << ( valtruerror / valerror ) <<  "\t" << sttime <<std::endl;
					for (prt=0;prt<ervec.NElements();prt++) {
						out <<"error " << ervec[prt] << "  truerror = " << truervec[prt] << "  Effect " << effect[prt] << std::endl;
					}
				}
				
				out.flush();

				cmesh->Reference()->ResetReference();
				cmesh->LoadReferences();
				adapt.DeleteElements(cmesh);
				delete cmesh;
				cmesh = 0;
				if(NRefs>1) {
					cmesh = adaptmesh;
					cmesh->CleanUpUnconnectedNodes();
				}
			}
			if(gmesh)
				delete gmesh;
		}
	}
	
	fileerrors << std::endl << "Finished running.\n" << std::endl << std::endl;
	fileerrors.close();
    std::cout << std::endl << "Finished running.\n" << std::endl << std::endl;
	out.close();
    return true;
}

bool SolveLaplaceProblemOnLShapeMesh() {
    // To compute processing times
	time_t sttime;
	time_t endtime;
	int time_elapsed;
	char time_formated[256];
	memset(time_formated,0,256);
	
	// Output files
    std::ofstream convergence("convergenceLP.txt");
	std::ofstream fileerrors("ErrorsHP_Laplace.txt");   // To store all errors calculated by TPZAnalysis (PosProcess)
	// Initial message to print computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	int nref = 1, NRefs = 12;
    int ninitialrefs = 2;
	int nthread = 1, NThreads = 4;
	
    //Working on regular meshes
    for(int regular=1; regular>0; regular--) {
		fileerrors << "Type of mesh: " << regular << " Level. " << endl;
		MElementType typeel;
		for(int itypeel=(int)ETriangle;itypeel<(int)ETetraedro;itypeel++)
		{
            MaxPOrder = 7;
			typeel = (MElementType)itypeel;
			fileerrors << "Type of element: " << typeel << endl;
            std::cout << "\nConstructing Poisson problem. Type element: " << typeel << std::endl;
			TPZGeoMesh *gmesh;
			if(!regular) {
				std::string nombre;
				GetFilenameFromGID(typeel,nombre);
				// Generating geometric mesh
				gmesh = CreateGeoMesh(nombre);
			}
			else {
				gmesh = CreateLShapeGeoMesh(typeel);
			}
			
			// Defining initial refinements and total refinements depends on dimension of the model
            UniformRefinement(ninitialrefs,gmesh,2);
            
			// Creating computational mesh (approximation space and materials)
			int p = 2, pinit;
			pinit = p;
			TPZCompEl::SetgOrder(p);
			TPZCompMesh *cmesh = CreateMesh(gmesh,2,1);
			gmesh->SetName("Malha Geometrica original");
			cmesh->SetName("Malha Computacional Original");
			
			// Printing geometric mesh to validate
			if(gDebug) {
				sprintf(saida,"gmeshL_H%dTR%dE%d.vtk",nref,regular,typeel);
				PrintGeoMeshAsCompMeshInVTKWithDimensionAsData(gmesh,saida);
			}
            
			// Selecting orthogonal polynomial family to construct shape functions
			if(!anothertests)
				TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;  // Setting Chebyshev polynomials as orthogonal sequence generating shape functions
			
			/** Variable names for post processing */
			TPZStack<std::string> scalnames, vecnames;
			scalnames.Push("POrder");
			scalnames.Push("Solution");
			
			// Solving adaptive process
			for(nref=1;nref<NRefs;nref++) {
				outLaplace << "\nConstructing Laplace problem. Refinement: " << nref << " Threads: " << nthread << " Regular: " << regular << " TypeElement: " << typeel << endl;
				if(nref > 5) nthread = 2*NThreads;
				else nthread = NThreads;
				
				// Initializing the generation mesh process
				time(& sttime);
				
				// Introduzing exact solution depending on the case
				TPZAnalysis an(cmesh);
				an.SetExact(ExactSolLaplace);
				{
					std::stringstream sout;
					sout << "Laplace_MESH" << regular << "E" << typeel << "Thr" << nthread << "H" << std::setprecision(2) << nref << "P" << pinit << ".vtk";
					an.DefineGraphMesh(2,scalnames,vecnames,sout.str());
				}
				std::string MeshFileName;
				{
					std::stringstream sout;
					sout << "Laplace_mesh_" << regular << "E" << typeel << "Thr" << nthread << "H" << std::setprecision(2) << nref << "P" << pinit << ".vtk";
					MeshFileName = sout.str();
				}
				
				cmesh->SetName("Malha computacional adaptada");
				// Printing geometric and computational mesh
				if(gDebug) {
					cmesh->Reference()->Print(std::cout);
					cmesh->Print(std::cout);
				}
				
				// Solve using symmetric matrix then using Cholesky (direct method)
				TPZSkylineStructMatrix strskyl(cmesh);
                if(usethreads) strskyl.SetNumThreads(nthread);
				an.SetStructuralMatrix(strskyl);
				
				TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
				direct->SetDirect(ECholesky);
				an.SetSolver(*direct);
				delete direct;
				direct = 0;
				
				an.Run();
				
				// Post processing
				an.PostProcess(0,2);
				if(gDebug) {
					std::ofstream out(MeshFileName.c_str());
					cmesh->LoadReferences();
					TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(),out,false);
				}
				// generation mesh process finished
				
				time(&endtime);
				time_elapsed = endtime - sttime;
				formatTimeInSec(time_formated,256,time_elapsed);
				outLaplace << "\tRefinement: " << nref+1 << " Regular Mesh: " << regular << " TypeElement: " << typeel << " Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n\n";
				
				// Initializing the auto adaptive process
				REAL valerror =0.;
				REAL valtruerror=0.;
				TPZVec<REAL> ervec,truervec,effect;
				
				TPZAdaptMesh adapt(MaxPOrder);
				adapt.SetCompMesh(cmesh);
				
				outLaplace << "\n\nEntering Auto Adaptive Methods... step " << nref << "\n";
				fileerrors << "\n\nEntering Auto Adaptive Methods... step " << nref << "\n";
                std::cout << "\n\nEntering Auto Adaptive Methods... step " << nref << "\n";
                
				TPZCompMesh *adaptmesh;
				if(NRefs>1) {
					time(&sttime);
					adaptmesh = adapt.GetAdaptedMesh(valerror,valtruerror,ervec,ExactSolLaplace,truervec,effect,fileerrors,0,typeel);
                    if(!adaptmesh) return false;
					
					time_t endtime;
					time(&endtime);
					
					int time_elapsed = endtime - sttime;
					outLaplace << "\n\nExiting Auto Adaptive Methods....step " << nref << "time elapsed " << time_elapsed << "\n\n\n\n";
					fileerrors << "\n\nExiting Auto Adaptive Methods....step " << nref << "time elapsed " << time_elapsed << "\n\n\n\n";
					
					int prt;
					outLaplace << "neq = " << cmesh->NEquations() << " error estimate = " << valerror << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
					fileerrors << "neq = " << cmesh->NEquations() << " error estimate = " << valerror << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
                    std::cout << "neq = " << cmesh->NEquations() << " error estimate = " << valerror << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
                    
					convergence  << cmesh->NEquations() << "\t" << valerror << "\t" << valtruerror << "\t" << ( valtruerror / valerror ) <<  "\t" << sttime <<std::endl;
					for (prt=0;prt<ervec.NElements();prt++) {
						outLaplace <<"error " << ervec[prt] << "  truerror = " << truervec[prt] << "  Effect " << effect[prt] << std::endl;
					}
				}
				
				outLaplace.flush();
				cmesh->Reference()->ResetReference();
				cmesh->LoadReferences();
				adapt.DeleteElements(cmesh);
				delete cmesh;
				cmesh = 0;
				if(NRefs>1) {
					cmesh = adaptmesh;
					cmesh->CleanUpUnconnectedNodes();
				}
			}
			if(gmesh)
				delete gmesh;
		}
	}
	
	fileerrors << std::endl << std::endl;
	fileerrors.close();
	outLaplace.close();
    return true;
}

TPZGeoMesh *CreateLShapeGeoMesh(MElementType typeel) {
	TPZGeoMesh* gmesh = new TPZGeoMesh;
	
	switch(typeel) {
		case ETriangle:
		{
            const int nelem = 12;
            const int nnodes = 11;
            REAL co[nnodes][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.},{0.5,-0.5},{0.5,0.5},{-0.5,0.5}};
            int indices[nelem][3] = {{0,3,8},{0,1,8},{1,2,8},{2,3,8},{0,3,9},{3,4,9},{4,5,9},{0,5,9},{0,7,10},{0,5,10},{5,6,10},{6,7,10}};
            TPZGeoEl *elvec[nelem];
            int nod;
            for(nod=0; nod<nnodes; nod++) {
                int nodind = gmesh->NodeVec().AllocateNewElement();
                TPZVec<REAL> coord(3,0.0);
                coord[0] = co[nod][0];
                coord[1] = co[nod][1];
                gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
            }
            
            int el;
            for(el=0; el<nelem; el++) {
                TPZVec<int64_t> nodind(3);
                for(nod=0; nod<3; nod++) nodind[nod]=indices[el][nod];
                int64_t index;
                elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
            }
            gmesh->BuildConnectivity();
            
            // bc -1 -> Dirichlet (Fixed - 0.0 in this side)
            TPZGeoElBC gbc1(elvec[1],3,id_bc0);
            // bc -2 -> Dirichlet with value of exact solution on this side
            TPZGeoElBC gbc2(elvec[2],3,id_bc1);
            TPZGeoElBC gbc3(elvec[3],3,id_bc1);
            TPZGeoElBC gbc4(elvec[5],3,id_bc1);
            TPZGeoElBC gbc5(elvec[6],3,id_bc1);
            TPZGeoElBC gbc6(elvec[8],3,id_bc1);
            TPZGeoElBC gbc7(elvec[10],3,id_bc1);
            TPZGeoElBC gbc8(elvec[11],3,id_bc1);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case EQuadrilateral:
		{
            REAL co[8][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.}};
            int indices[3][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7}};
            TPZGeoEl *elvec[3];
            int nnode = 8;
            int nod;
            for(nod=0; nod<nnode; nod++) {
                int nodind = gmesh->NodeVec().AllocateNewElement();
                TPZVec<REAL> coord(3,0.0);
                coord[0] = co[nod][0];
                coord[1] = co[nod][1];
                gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
            }
            
            int el;
            int nelem = 3;
            for(el=0; el<nelem; el++) {
                TPZVec<int64_t> nodind(4);
                for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
                int64_t index;
                elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
            }
            
            gmesh->BuildConnectivity();
            // bc -1 -> Dirichlet
            TPZGeoElBC gbc1(elvec[0],4,id_bc0);
            // bc -2 -> Dirichlet with value of exact solution on this side
            TPZGeoElBC gbc2(elvec[0],5,id_bc1);
            TPZGeoElBC gbc3(elvec[0],6,id_bc1);
            TPZGeoElBC gbc4(elvec[1],5,id_bc1);
            TPZGeoElBC gbc5(elvec[1],6,id_bc1);
            TPZGeoElBC gbc6(elvec[2],5,id_bc1);
            TPZGeoElBC gbc7(elvec[2],6,id_bc1);
            TPZGeoElBC gbc8(elvec[2],7,id_bc1);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
        default:
            break;
    }
    return gmesh;
}

//**** Creating Geometric Mesh as square */
TPZGeoMesh *CreateGeoMesh(MElementType typeel) {
	TPZManVector<REAL> point(3,0.), pointlast(3,0.);
	TPZGeoMesh* gmesh;
	switch (typeel) {
		case EOned:
		{
			pointlast[0] = 1.;
			gmesh = new TPZGeoMesh;
			int Qnodes = 2;
			
			gmesh->SetMaxNodeId(Qnodes-1);
			gmesh->NodeVec().Resize(Qnodes);
			TPZVec<TPZGeoNode> Node(Qnodes);
			
			TPZVec <int64_t> TopolLine(2);
			TPZVec <int64_t> TopolPoint(1);
			
			int64_t id = 0;
			for (int j=0; j<2;j++) {
				Node[id].SetNodeId(id);
				if(!j) Node[id].SetCoord(point);//coord x
				else Node[id].SetCoord(pointlast);
				gmesh->NodeVec()[id] = Node[id];
				id++;
			}
			
			//indice dos elementos
			id = 0;
			
			TopolLine[0] = 0;
			TopolLine[1] = 1;
			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,materialId,*gmesh);
			id++;
			
			TopolPoint[0] = 0;
			new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,id_bc0,*gmesh);
			id++;
			TopolPoint[0] = 1;
			new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,id_bc0,*gmesh);

			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case EQuadrilateral:
		{
			const int nnode = 4;
			REAL co[nnode][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};
			
			TPZGeoEl *elvec[1];
			gmesh = new TPZGeoMesh();
			
			int nod;
			for(nod=0; nod<nnode; nod++) {
				int nodind = gmesh->NodeVec().AllocateNewElement();
				TPZVec<REAL> coord(3);
				coord[0] = co[nod][0];
				coord[1] = co[nod][1];
				coord[2] = co[nod][2];
				gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
			}
			
			int el = 0;
			TPZVec<int64_t> nodind(nnode);
			for(nod=0; nod<nnode; nod++) nodind[nod]=nod;
			int64_t index;
			elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,materialId,index);
			
			gmesh->BuildConnectivity();
			
			// bc -1 -> Dirichlet
			TPZGeoElBC gbc1(elvec[0],4,id_bc0);
			TPZGeoElBC gbc2(elvec[0],5,id_bc0);
			TPZGeoElBC gbc3(elvec[0],6,id_bc0);
			TPZGeoElBC gbc4(elvec[0],7,id_bc0);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case ETriangle:
		{
			const int nelem = 4;
			const int nnode = 5;
			
			REAL co[nnode][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.},{0.5,0.5}};
            int indices[nelem][nnode] = {{0,1,4},{1,2,4},{2,3,4},{0,3,4}};
			
			TPZGeoEl *elvec[nelem];
			gmesh = new TPZGeoMesh();
			
			int nod;
			for(nod=0; nod<nnode; nod++) {
				int64_t nodind = gmesh->NodeVec().AllocateNewElement();
				TPZVec<REAL> coord(3,0.0);
				coord[0] = co[nod][0];
				coord[1] = co[nod][1];
				gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
			}
			
			int el;
			for(el=0; el<nelem; el++) {
				TPZVec<int64_t> nodind(3);
				for(nod=0; nod<3; nod++) nodind[nod]=indices[el][nod];
				int64_t index;
				elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
			}
			
			gmesh->BuildConnectivity();
			
			// bc -1 -> Dirichlet
			TPZGeoElBC gbc1(elvec[0],3,id_bc0);
			TPZGeoElBC gbc2(elvec[1],3,id_bc0);
			TPZGeoElBC gbc3(elvec[2],3,id_bc0);
			TPZGeoElBC gbc4(elvec[3],3,id_bc0);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case ETetraedro:
			gmesh = ConstructingTetrahedraInCube(1.);
			break;
		case EPrisma:
			gmesh = ConstructingPrismsInCube(1.);
			break;
		case EPiramide:
			gmesh = ConstructingPyramidsInCube(1.);
			break;
		case ECube:
			gmesh = ConstructingPositiveCube(1.,typeel);
			break;
		default:
            gmesh = 0;
			break;
	}

	return gmesh;
}

#include "TPZRefPatternDataBase.h"
TPZGeoMesh *ConstructingPositiveCube(REAL InitialL,MElementType typeel) {
	// CREATING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // Dependig on dimension of the typeel
	const int nelem = 1;
	const int nnode = 8;
	REAL co[nnode][3] = {
		{0.,0.,0.},
		{InitialL,0.,0.},
		{InitialL,InitialL,0.},
		{0.,InitialL,0.},
		{0.,0.,InitialL},
		{InitialL,0.,InitialL},
		{InitialL,InitialL,InitialL},
		{0.,InitialL,InitialL}
	};
	TPZVec<TPZVec<int64_t> > indices(nelem);
	indices[0].Resize(nnode);
	int nod;
	for(nod=0;nod<nnode;nod++)
		indices[0][nod] = nod;
	
	TPZGeoEl *elvec[nelem];
	TPZGeoMesh *gmesh = new TPZGeoMesh();

	for(nod=0; nod<nnode; nod++) {
		int64_t nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord(3);
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		coord[2] = co[nod][2];
		gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
	}
	
	int el;
	for(el=0; el<nelem; el++) {
		int64_t index;
		elvec[el] = gmesh->CreateGeoElement(ECube,indices[el],1,index);
	}
    gmesh->BuildConnectivity();
	
	// Introduzing boundary condition for cube - ALL DIRICHLET
    // Boundary condition on one dimensional sides
	TPZGeoElBC gbc00(gmesh->ElementVec()[0],8,id_bc0);
	TPZGeoElBC gbc01(gmesh->ElementVec()[0],9,id_bc0);
	TPZGeoElBC gbc02(gmesh->ElementVec()[0],10,id_bc0);
	TPZGeoElBC gbc03(gmesh->ElementVec()[0],11,id_bc0);
	TPZGeoElBC gbc04(gmesh->ElementVec()[0],12,id_bc0);
	TPZGeoElBC gbc05(gmesh->ElementVec()[0],13,id_bc0);
	TPZGeoElBC gbc06(gmesh->ElementVec()[0],14,id_bc0);
	TPZGeoElBC gbc07(gmesh->ElementVec()[0],15,id_bc0);
	TPZGeoElBC gbc08(gmesh->ElementVec()[0],16,id_bc0);
	TPZGeoElBC gbc09(gmesh->ElementVec()[0],17,id_bc0);
	TPZGeoElBC gbc10(gmesh->ElementVec()[0],18,id_bc0);
	TPZGeoElBC gbc11(gmesh->ElementVec()[0],19,id_bc0);
	// face 0 (20) bottom XY - face 1 (21) lateral left XZ - face 4 (24) lateral back YZ : Dirichlet
	TPZGeoElBC gbc20(gmesh->ElementVec()[0],20,id_bc0);
	TPZGeoElBC gbc21(gmesh->ElementVec()[0],21,id_bc0);
	TPZGeoElBC gbc22(gmesh->ElementVec()[0],24,id_bc0);
	// face 2 (22) Neumann - Partial derivative (du/dx) - lateral front
	TPZGeoElBC gbc23(gmesh->ElementVec()[0],22,id_bc0);
	// face 3 (23) Neumann - Partial derivative (du/dy) - lateral right
	TPZGeoElBC gbc24(gmesh->ElementVec()[0],23,id_bc0);
	// face 5 (25) Neumann - Partial derivative (du/dz) - top
	TPZGeoElBC gbc25(gmesh->ElementVec()[0],25,id_bc0);

	TPZVec<TPZGeoEl *> sub;
	std::string filename = REFPATTERNDIR;
    switch (typeel) {
        case ENoType:
        {
            char buf[1024];
			std::istringstream str(buf);
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
			TPZAutoPointer<TPZRefPattern> refpatFound = gRefDBase.FindRefPattern(refpat);
			if(!refpatFound)
			{
				gRefDBase.InsertRefPattern(refpat);
			}
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
        }
    break;

        case EPrisma:   // hexahedron -> four prisms
		{
            // Dividing hexahedron in four prisms (anymore)
            filename += "/3D_Hexa_directional_2faces.rpt";
            
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
			TPZAutoPointer<TPZRefPattern> refpatFound = gRefDBase.FindRefPattern(refpat);
            if(!refpatFound)
            {
                gRefDBase.InsertRefPattern(refpat);
            }
			else
			{
				refpatFound->SetName(refpat->Name());
			}
			refpat->InsertPermuted();
            TPZGeoEl *gel = gmesh->ElementVec()[0];
            TPZGeoElRefPattern <TPZGeoCube> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
            gelrp->SetRefPattern(refpat);
            gel->Divide(sub);
		}   
            break;
        case EPiramide:
		{
            // Dividing hexahedron in four pyramids (anymore)
            filename += "/3D_Hexa_Rib_Side_08.rpt";
            
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
            if(!gRefDBase.FindRefPattern(refpat))
            {
                gRefDBase.InsertRefPattern(refpat);
            }
            TPZGeoEl *gel = gmesh->ElementVec()[0];
            TPZGeoElRefPattern <TPZGeoCube> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
            gelrp->SetRefPattern(refpat);
            gel->Divide(sub);
		}
			break;
        case ETetraedro:
		{
            // Dividing hexahedron in two tetrahedras, two prisms and one pyramid
            filename += "/3D_Hexa_Rib_Side_16_17_18.rpt";
            
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
            if(!gRefDBase.FindRefPattern(refpat))
            {
                gRefDBase.InsertRefPattern(refpat);
            }
            TPZGeoEl *gel = gmesh->ElementVec()[0];
            TPZGeoElRefPattern <TPZGeoCube> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
            gelrp->SetRefPattern(refpat);
            gel->Divide(sub);
		}
			break;
		case ECube:
			break;
        default:
		{
            // hexahedron -> three prisms (anymore)
            // Dividing hexahedron in prisms
            filename += "/3D_Hexa_Rib_Side_16_18.rpt";
            
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
            if(!gRefDBase.FindRefPattern(refpat))
            {
                gRefDBase.InsertRefPattern(refpat);
            }
            TPZGeoEl *gel = gmesh->ElementVec()[0];
            TPZGeoElRefPattern <TPZGeoCube> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
            gelrp->SetRefPattern(refpat);
            gel->Divide(sub);
		}
            break;
    }
	
//	gmesh->BuildConnectivity();
	
/*	switch(typeel) {
		case ECube:
		{
			// face 0 (20) bottom XY - face 1 (21) lateral left XZ - face 4 (24) lateral back YZ : Dirichlet
			TPZGeoElBC gbc10(gmesh->ElementVec()[0],20,-1);
			TPZGeoElBC gbc11(gmesh->ElementVec()[0],21,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[0],24,-1);
			
			// face 2 (22) Neumann - Partial derivative (du/dx) - lateral front
			TPZGeoElBC gbc13(gmesh->ElementVec()[0],22,-1);
			// face 3 (23) Neumann - Partial derivative (du/dy) - lateral right
			TPZGeoElBC gbc14(gmesh->ElementVec()[0],23,-1);
			// face 5 (25) Neumann - Partial derivative (du/dz) - top
			TPZGeoElBC gbc15(gmesh->ElementVec()[0],25,-1);
		}
			break;
		case EPrisma:
		{
			// First sub element - faces: 15, 16, 18 and 19
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],15,-1);
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],16,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],18,-1);
			TPZGeoElBC gbc13(gmesh->ElementVec()[1],19,-1);
			// Second sub element - faces: 15 and 19
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],15,-1);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],19,-1);
			// Third sub element - faces: 15, 16, 17 and 19
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],15,-1);
			TPZGeoElBC gbc31(gmesh->ElementVec()[3],17,-1);
			TPZGeoElBC gbc32(gmesh->ElementVec()[3],18,-1);
			TPZGeoElBC gbc33(gmesh->ElementVec()[3],19,-1);
			// Fouthrd sub element - faces: 15, 17, 18 and 19
			TPZGeoElBC gbc40(gmesh->ElementVec()[4],15,-1);
			TPZGeoElBC gbc41(gmesh->ElementVec()[4],17,-1);
			TPZGeoElBC gbc42(gmesh->ElementVec()[4],18,-1);
			TPZGeoElBC gbc43(gmesh->ElementVec()[4],19,-1);
			gmesh->ElementVec()[1]->Divide(sub);
			gmesh->ElementVec()[3]->Divide(sub);
		}
			break;
		case EPiramide:
		{
			// First sub element - faces: 13, 14 and 17
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],13,-1);
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],14,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],17,-1);
			// Second sub element - faces: 13 and 14
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],13,-1);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],14,-1);
			// Third sub element - faces: 13, 14 and 17
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],13,-1);
			TPZGeoElBC gbc31(gmesh->ElementVec()[3],14,-1);
			TPZGeoElBC gbc32(gmesh->ElementVec()[3],17,-1);
			// Fouthrd sub element - faces: 13 and 14
			TPZGeoElBC gbc40(gmesh->ElementVec()[4],13,-1);
			TPZGeoElBC gbc41(gmesh->ElementVec()[4],14,-1);
		}
			break;
		case ETetraedro:
		{
			// First sub element - faces: 10, 11 and 13
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],10,-1);
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],11,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],13,-1);
			// Second sub element - faces: 10, 11 and 13
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],10,-1);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],11,-1);
			TPZGeoElBC gbc22(gmesh->ElementVec()[2],13,-1);
			// Third sub element - faces: 15, 16, 18 and 19
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],15,-1);
			TPZGeoElBC gbc31(gmesh->ElementVec()[3],16,-1);
			TPZGeoElBC gbc32(gmesh->ElementVec()[3],18,-1);
			TPZGeoElBC gbc33(gmesh->ElementVec()[3],19,-1);
			// Fouthrd sub element - faces: 15, 18 and 19
			TPZGeoElBC gbc40(gmesh->ElementVec()[4],15,-1);
			TPZGeoElBC gbc41(gmesh->ElementVec()[4],18,-1);
			TPZGeoElBC gbc42(gmesh->ElementVec()[4],19,-1);
			// Fifth sub element - faces: 13 and 15
			TPZGeoElBC gbc50(gmesh->ElementVec()[5],14,-1);
			TPZGeoElBC gbc51(gmesh->ElementVec()[5],16,-1);
			gmesh->ElementVec()[3]->Divide(sub);
			gmesh->ElementVec()[4]->Divide(sub);
		}
			break;
		default:
		{
			// First sub element - faces: 15, 16, 18 and 19
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],15,-1);
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],16,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],18,-1);
			TPZGeoElBC gbc13(gmesh->ElementVec()[1],19,-1);
			// Second sub element - faces: 15, 16, 18 and 19
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],15,-1);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],16,-1);
			TPZGeoElBC gbc22(gmesh->ElementVec()[2],18,-1);
			TPZGeoElBC gbc23(gmesh->ElementVec()[2],19,-1);
			// Third sub element - faces: 15, 17 and 19
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],15,-1);
			TPZGeoElBC gbc31(gmesh->ElementVec()[3],17,-1);
			TPZGeoElBC gbc32(gmesh->ElementVec()[3],19,-1);
			gmesh->ElementVec()[1]->Divide(sub);
			gmesh->ElementVec()[2]->Divide(sub);
			gmesh->ElementVec()[3]->Divide(sub);
		}
			break;
	}*/
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	return gmesh;
}

TPZGeoMesh *ConstructingTetrahedraInCube(REAL InitialL) {
	// CONSIDERING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // And dividing into five tetrahedras
	TPZGeoMesh *gmesh = new TPZGeoMesh();
    
	const int nelem = 5;
	const int nnode = 8;
	REAL co[nnode][3] = {
		{0.,0.,0.},
		{InitialL,0.,0.},
		{InitialL,InitialL,0.},
		{0.,InitialL,0.},
		{0.,0.,InitialL},
		{InitialL,0.,InitialL},
		{InitialL,InitialL,InitialL},
		{0.,InitialL,InitialL},
	};
	int nod;
	for(nod=0; nod<nnode; nod++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord(3);
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		coord[2] = co[nod][2];
		gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
	}
    
	TPZVec<TPZVec<int64_t> > indices(nelem);
	int nnodebyelement = 4;
	int el;
	for(el=0;el<nelem;el++)
		indices[el].Resize(nnodebyelement);
	// nodes to first element
	indices[0][0] = 0;
	indices[0][1] = 1;
	indices[0][2] = 3;
	indices[0][3] = 4;
	// nodes to second element
	indices[1][0] = 1;
	indices[1][1] = 2;
	indices[1][2] = 3;
	indices[1][3] = 6;
	// nodes to third element
	indices[2][0] = 4;
	indices[2][1] = 5;
	indices[2][2] = 6;
	indices[2][3] = 1;
	// nodes to fourth element
	indices[3][0] = 6;
	indices[3][1] = 7;
	indices[3][2] = 4;
	indices[3][3] = 3;
	// nodes to fifth element
	indices[4][0] = 1;
	indices[4][1] = 4;
	indices[4][2] = 6;
	indices[4][3] = 3;
    
	TPZGeoEl *elvec[nelem];
	for(el=0; el<nelem; el++) {
		int64_t index;
		elvec[el] = gmesh->CreateGeoElement(ETetraedro,indices[el],materialId,index);
	}
    gmesh->BuildConnectivity();
	
	// Introduzing boundary condition for cube - ALL DIRICHLET
    // boundary condition over one dimensional sides
	TPZGeoElBC gbc01(gmesh->ElementVec()[0],4,id_bc0);
	TPZGeoElBC gbc02(gmesh->ElementVec()[0],5,id_bc0);
	TPZGeoElBC gbc03(gmesh->ElementVec()[0],6,id_bc0);
	TPZGeoElBC gbc04(gmesh->ElementVec()[0],7,id_bc0);
	TPZGeoElBC gbc05(gmesh->ElementVec()[0],8,id_bc0);
	TPZGeoElBC gbc06(gmesh->ElementVec()[0],9,id_bc0);
	TPZGeoElBC gbc07(gmesh->ElementVec()[1],5,id_bc0);
	TPZGeoElBC gbc08(gmesh->ElementVec()[1],6,id_bc0);
	TPZGeoElBC gbc09(gmesh->ElementVec()[1],7,id_bc0);
	TPZGeoElBC gbc10(gmesh->ElementVec()[1],8,id_bc0);
	TPZGeoElBC gbc11(gmesh->ElementVec()[1],9,id_bc0);
	TPZGeoElBC gbc12(gmesh->ElementVec()[2],4,id_bc0);
	TPZGeoElBC gbc13(gmesh->ElementVec()[2],5,id_bc0);
	TPZGeoElBC gbc14(gmesh->ElementVec()[2],6,id_bc0);
	TPZGeoElBC gbc15(gmesh->ElementVec()[2],8,id_bc0);
	TPZGeoElBC gbc16(gmesh->ElementVec()[3],4,id_bc0);
	TPZGeoElBC gbc17(gmesh->ElementVec()[3],5,id_bc0);
	TPZGeoElBC gbc18(gmesh->ElementVec()[3],8,id_bc0);
    
	// face 0 (20) bottom XY
	TPZGeoElBC gbc20(gmesh->ElementVec()[0],10,id_bc0);
	TPZGeoElBC gbc21(gmesh->ElementVec()[0],11,id_bc0);
	TPZGeoElBC gbc22(gmesh->ElementVec()[0],13,id_bc0);
	TPZGeoElBC gbc23(gmesh->ElementVec()[1],10,id_bc0);
	TPZGeoElBC gbc24(gmesh->ElementVec()[1],11,id_bc0);
	TPZGeoElBC gbc25(gmesh->ElementVec()[1],12,id_bc0);
	TPZGeoElBC gbc26(gmesh->ElementVec()[2],10,id_bc0);
	TPZGeoElBC gbc27(gmesh->ElementVec()[2],11,id_bc0);
	TPZGeoElBC gbc28(gmesh->ElementVec()[2],12,id_bc0);
	TPZGeoElBC gbc29(gmesh->ElementVec()[3],10,id_bc0);
	TPZGeoElBC gbc30(gmesh->ElementVec()[3],11,id_bc0);
	TPZGeoElBC gbc31(gmesh->ElementVec()[3],12,id_bc0);
/*	std::string filename = REFPATTERNDIR;
    filename += "/2D_Triang_Rib_3.rpt";
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
    TPZAutoPointer<TPZRefPattern> refpatFound = gRefDBase.FindRefPattern(refpat);
    if(!refpatFound)
        gRefDBase.InsertRefPattern(refpat);
    else
        refpatFound->SetName(refpat->Name());
    refpat->InsertPermuted();

    TPZGeoEl* gel = gbc20.CreatedElement();
    TPZGeoElRefPattern <TPZGeoTriangle> *gelrp;
    if(gel) {
        gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
        gelrp->SetRefPattern(refpat);
    }
    gel = gbc21.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc22.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc23.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc24.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc25.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc26.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc27.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc28.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc29.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc30.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
    gel = gbc31.CreatedElement();
    gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle> *> (gel);
    gelrp->SetRefPattern(refpat);
*/
	return gmesh;
}

TPZGeoMesh *ConstructingPyramidsInCube(REAL InitialL) {
	// CONSIDERING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // And dividing into six pyramids
	TPZGeoMesh *gmesh = new TPZGeoMesh();

	const int nelem = 6;
	const int nnode = 9;
	REAL co[nnode][3] = {
		{0.,0.,0.},
		{InitialL,0.,0.},
		{InitialL,InitialL,0.},
		{0.,InitialL,0.},
		{0.,0.,InitialL},
		{InitialL,0.,InitialL},
		{InitialL,InitialL,InitialL},
		{0.,InitialL,InitialL},
		{0.5*InitialL,0.5*InitialL,0.5*InitialL}
	};
	int nod;
	for(nod=0; nod<nnode; nod++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord(3);
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		coord[2] = co[nod][2];
		gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
	}

	TPZVec<TPZVec<int64_t> > indices(nelem);
	int nnodebyelement = 5;
	int el;
	for(el=0;el<nelem;el++)
		indices[el].Resize(nnodebyelement);
	// nodes to first element
	indices[0][0] = 0;
	indices[0][1] = 1;
	indices[0][2] = 2;
	indices[0][3] = 3;
	indices[0][4] = 8;
	// nodes to second element
	indices[1][0] = 0;
	indices[1][1] = 1;
	indices[1][2] = 5;
	indices[1][3] = 4;
	indices[1][4] = 8;
	// nodes to third element
	indices[2][0] = 1;
	indices[2][1] = 2;
	indices[2][2] = 6;
	indices[2][3] = 5;
	indices[2][4] = 8;
	// nodes to fourth element
	indices[3][0] = 2;
	indices[3][1] = 3;
	indices[3][2] = 7;
	indices[3][3] = 6;
	indices[3][4] = 8;
	// nodes to fifth element
	indices[4][0] = 0;
	indices[4][1] = 3;
	indices[4][2] = 7;
	indices[4][3] = 4;
	indices[4][4] = 8;
	// nodes to sixth element
	indices[5][0] = 4;
	indices[5][1] = 5;
	indices[5][2] = 6;
	indices[5][3] = 7;
	indices[5][4] = 8;

	TPZGeoEl *elvec[nelem];
	for(el=0; el<nelem; el++) {
		int64_t index;
		elvec[el] = gmesh->CreateGeoElement(EPiramide,indices[el],materialId,index);
	}
    gmesh->BuildConnectivity();
	
	// Introduzing boundary condition for cube - ALL DIRICHLET
    // boundary condition on one dimensional sides
	TPZGeoElBC gbc01(gmesh->ElementVec()[0],5,id_bc0);
	TPZGeoElBC gbc02(gmesh->ElementVec()[0],6,id_bc0);
	TPZGeoElBC gbc03(gmesh->ElementVec()[0],7,id_bc0);
	TPZGeoElBC gbc04(gmesh->ElementVec()[0],8,id_bc0);
	TPZGeoElBC gbc05(gmesh->ElementVec()[1],6,id_bc0);
	TPZGeoElBC gbc06(gmesh->ElementVec()[1],7,id_bc0);
	TPZGeoElBC gbc07(gmesh->ElementVec()[1],8,id_bc0);
	TPZGeoElBC gbc08(gmesh->ElementVec()[2],6,id_bc0);
	TPZGeoElBC gbc09(gmesh->ElementVec()[2],7,id_bc0);
	TPZGeoElBC gbc10(gmesh->ElementVec()[3],6,id_bc0);
	TPZGeoElBC gbc11(gmesh->ElementVec()[3],7,id_bc0);
	TPZGeoElBC gbc12(gmesh->ElementVec()[4],7,id_bc0);
    
	// face 0 (20) bottom XY 
	TPZGeoElBC gbc21(gmesh->ElementVec()[0],13,id_bc0);
	TPZGeoElBC gbc22(gmesh->ElementVec()[1],13,id_bc0);
	TPZGeoElBC gbc23(gmesh->ElementVec()[2],13,id_bc0);
	TPZGeoElBC gbc24(gmesh->ElementVec()[3],13,id_bc0);
	TPZGeoElBC gbc25(gmesh->ElementVec()[4],13,id_bc0);
	TPZGeoElBC gbc26(gmesh->ElementVec()[5],13,id_bc0);

	return gmesh;
}

TPZGeoMesh *ConstructingPrismsInCube(REAL InitialL) {
	// CREATING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // And dividing into four prisms
	TPZGeoMesh *gmesh = new TPZGeoMesh();

	const int nelem = 4;
	const int nnode = 12;
	REAL co[nnode][3] = {
		{0.,0.,0.},
		{InitialL,0.,0.},
		{InitialL,InitialL,0.},
		{0.,InitialL,0.},
		{0.,0.,InitialL},
		{InitialL,0.,InitialL},
		{InitialL,InitialL,InitialL},
		{0.,InitialL,InitialL},
		{0.,0.,0.5*InitialL},
		{InitialL,0.,0.5*InitialL},
		{InitialL,InitialL,0.5*InitialL},
		{0.,InitialL,0.5*InitialL},
	};
	for(int nod=0; nod<nnode; nod++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord(3);
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		coord[2] = co[nod][2];
		gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
	}

	TPZVec<TPZVec<int64_t> > indices(nelem);
	int nnodebyelement = 6;
	int el;
	for(el=0;el<nelem;el++)
		indices[el].Resize(nnodebyelement);
	// nodes to first element
	indices[0][0] = 0;
	indices[0][1] = 1;
	indices[0][2] = 2;
	indices[0][3] = 8;
	indices[0][4] = 9;
	indices[0][5] = 10;
	// nodes to second element
	indices[1][0] = 0;
	indices[1][1] = 2;
	indices[1][2] = 3;
	indices[1][3] = 8;
	indices[1][4] = 10;
	indices[1][5] = 11;
	// nodes to third element
	indices[2][0] = 8;
	indices[2][1] = 9;
	indices[2][2] = 10;
	indices[2][3] = 4;
	indices[2][4] = 5;
	indices[2][5] = 6;
	// nodes to fourth element
	indices[3][0] = 8;
	indices[3][1] = 10;
	indices[3][2] = 11;
	indices[3][3] = 4;
	indices[3][4] = 6;
	indices[3][5] = 7;

	TPZGeoEl *elvec[nelem];
	for(el=0; el<nelem; el++) {
		int64_t index;
		elvec[el] = gmesh->CreateGeoElement(EPrisma,indices[el],materialId,index);
	}
    gmesh->BuildConnectivity();
	
	// Introduzing boundary condition for cube - ALL DIRICHLET
    // boundary condition on one dimensional sides
	TPZGeoElBC gbc00(gmesh->ElementVec()[0],6,id_bc0);
	TPZGeoElBC gbc01(gmesh->ElementVec()[0],7,id_bc0);
	TPZGeoElBC gbc02(gmesh->ElementVec()[0],8,id_bc0);
	TPZGeoElBC gbc03(gmesh->ElementVec()[0],9,id_bc0);
	TPZGeoElBC gbc04(gmesh->ElementVec()[0],10,id_bc0);
	TPZGeoElBC gbc05(gmesh->ElementVec()[0],11,id_bc0);
	TPZGeoElBC gbc06(gmesh->ElementVec()[0],12,id_bc0);
	TPZGeoElBC gbc07(gmesh->ElementVec()[0],13,id_bc0);
	TPZGeoElBC gbc08(gmesh->ElementVec()[1],7,id_bc0);
	TPZGeoElBC gbc09(gmesh->ElementVec()[1],8,id_bc0);
	TPZGeoElBC gbc10(gmesh->ElementVec()[1],11,id_bc0);
	TPZGeoElBC gbc11(gmesh->ElementVec()[1],13,id_bc0);
	TPZGeoElBC gbc12(gmesh->ElementVec()[1],14,id_bc0);
	TPZGeoElBC gbc13(gmesh->ElementVec()[2],9,id_bc0);
	TPZGeoElBC gbc14(gmesh->ElementVec()[2],10,id_bc0);
	TPZGeoElBC gbc15(gmesh->ElementVec()[2],11,id_bc0);
	TPZGeoElBC gbc16(gmesh->ElementVec()[2],12,id_bc0);
	TPZGeoElBC gbc17(gmesh->ElementVec()[2],13,id_bc0);
	TPZGeoElBC gbc18(gmesh->ElementVec()[2],14,id_bc0);
	TPZGeoElBC gbc19(gmesh->ElementVec()[3],11,id_bc0);
	TPZGeoElBC gbc20(gmesh->ElementVec()[3],13,id_bc0);
	TPZGeoElBC gbc21(gmesh->ElementVec()[3],14,id_bc0);
    
	// face 0 (20) bottom XY
	TPZGeoElBC gbc30(gmesh->ElementVec()[0],15,id_bc0);
	TPZGeoElBC gbc31(gmesh->ElementVec()[0],16,id_bc0);
	TPZGeoElBC gbc32(gmesh->ElementVec()[0],17,id_bc0);
	TPZGeoElBC gbc33(gmesh->ElementVec()[1],15,id_bc0);
	TPZGeoElBC gbc34(gmesh->ElementVec()[1],17,id_bc0);
	TPZGeoElBC gbc35(gmesh->ElementVec()[1],18,id_bc0);
	TPZGeoElBC gbc36(gmesh->ElementVec()[2],16,id_bc0);
	TPZGeoElBC gbc37(gmesh->ElementVec()[2],17,id_bc0);
	TPZGeoElBC gbc38(gmesh->ElementVec()[2],19,id_bc0);
	TPZGeoElBC gbc39(gmesh->ElementVec()[3],17,id_bc0);
	TPZGeoElBC gbc40(gmesh->ElementVec()[3],18,id_bc0);
	TPZGeoElBC gbc41(gmesh->ElementVec()[3],19,id_bc0);

	return gmesh;
}
TPZGeoMesh *ConstructingSeveral3DElementsInCube(REAL InitialL,MElementType typeel) {
	// CREATING A CUBE WITH MASS CENTER (0.5*INITIALL, 0.5*INITIALL, 0.5*INITIALL) AND VOLUME = INITIALL*INITIALL*INITIALL
    // Dependig on dimension of the typeel
	const int nelem = 1;
	const int nnode = 8;
	REAL co[nnode][3] = {
		{0.,0.,0.},
		{InitialL,0.,0.},
		{InitialL,InitialL,0.},
		{0.,InitialL,0.},
		{0.,0.,InitialL},
		{InitialL,0.,InitialL},
		{InitialL,InitialL,InitialL},
		{0.,InitialL,InitialL}
	};
	TPZVec<TPZVec<int64_t> > indices(nelem);
	indices[0].Resize(nnode);
	int nod;
	for(nod=0;nod<nnode;nod++)
		indices[0][nod] = nod;
	
	TPZGeoEl *elvec[nelem];
	TPZGeoMesh *gmesh = new TPZGeoMesh();

	for(nod=0; nod<nnode; nod++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord(3);
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		coord[2] = co[nod][2];
		gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
	}
	
	int el;
	for(el=0; el<nelem; el++) {
		int64_t index;
		elvec[el] = gmesh->CreateGeoElement(ECube,indices[el],1,index);
	}
    gmesh->BuildConnectivity();
	
	// Introduzing boundary condition for cube - ALL DIRICHLET
	// face 0 (20) bottom XY - face 1 (21) lateral left XZ - face 4 (24) lateral back YZ : Dirichlet
	TPZGeoElBC gbc10(gmesh->ElementVec()[0],20,id_bc0);
	TPZGeoElBC gbc11(gmesh->ElementVec()[0],21,id_bc0);
	TPZGeoElBC gbc12(gmesh->ElementVec()[0],24,id_bc0);
	// face 2 (22) Neumann - Partial derivative (du/dx) - lateral front
	TPZGeoElBC gbc13(gmesh->ElementVec()[0],22,id_bc0);
	// face 3 (23) Neumann - Partial derivative (du/dy) - lateral right
	TPZGeoElBC gbc14(gmesh->ElementVec()[0],23,id_bc0);
	// face 5 (25) Neumann - Partial derivative (du/dz) - top
	TPZGeoElBC gbc15(gmesh->ElementVec()[0],25,id_bc0);

	return gmesh;
}

//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh(std::string &archivo) {
	
	// Generacion de una malla utilizando o GID previamente 
	TPZReadGIDGrid grid;
	TPZGeoMesh *meshgrid = grid.GeometricGIDMesh(archivo);
	if(!meshgrid->NElements())
		return 0;
	
	return meshgrid;
}


//************************************************************************
//**********   Creating computational mesh with materials    *************
//************************************************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
    // Creating Poisson material
	TPZMaterial *mat = new TPZMatPoisson3d(materialId,dim);
	switch(hasforcingfunction) {
		case 1:
        {
            TPZVec<REAL> convd(3,0.);
            ((TPZMatPoisson3d *)mat)->SetParameters(ValueK,0.,convd);
			mat->SetForcingFunction(new TPZDummyFunction<STATE>(RightTermCircle));
        }
			break;
		case 2:
			break;
		default:
			break;
	}
    cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	
	// Creating four boundary condition
    TPZFMatrix<STATE> val1(dim,dim,0.),val2(dim,1,0.);
    for(int i=0;i<dim;i++)
        val1.PutVal(i,i,1.);
	TPZMaterial *bc = 0, *bc1 = 0;
    switch(hasforcingfunction) {
        case 0:
        case 1:
            // Condicion de Dirichlet fijando la posicion de la placa
            bc = mat->CreateBC(mat,id_bc0,0,val1,val2);
            break;
        case 2:
            // Condicion de Dirichlet fijando la posicion de la placa
            bc = mat->CreateBC(mat,id_bc0,0,val1,val2);
            bc1 = mat->CreateBC(mat,id_bc1,0,val1,val2);
            bc1->SetForcingFunction(new TPZDummyFunction<STATE>(ExactSolLaplaceBC));
            break;
        default:
            break;
    }

	if(bc) cmesh->InsertMaterialObject(bc);
	if(bc1) cmesh->InsertMaterialObject(bc1);

    cmesh->AutoBuild();
	cmesh->ExpandSolution();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}


/** PROBLEMS TO SOLVING */

/** ARC TANGENT PROBLEM */

/** We are considering - f, because is as TPZMatPoisson3d was implemented in Contribute method */
void RightTermCircle(const TPZVec<REAL> &x, TPZVec<STATE> &force, TPZFMatrix<STATE> &dforce) {
	int dim = dforce.Rows();
	
	REAL B = ValueK/M_PI;
	REAL F = 2*sqrt(ValueK);
	REAL Coeff;
	if(dim==1)
		 Coeff = -2.;
	else if(dim==2)
		Coeff = 8.;
	else
		Coeff = -32.;
	B *= (2.*Coeff);

	REAL Prod, prodx, prody, prodz;
    REAL Soma, arc, den;
    prodx = x[0]*(x[0] - 1.);
	if(dim==1) {
		arc = F*(RCircle*RCircle - (x[0]-CCircle[0])*(x[0]-CCircle[0]));
        prody = prodz = 1.;
        Soma = 1.;
	}
	else if(dim == 2) {
		arc = F*(RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1])));
        prody = x[1]*(x[1] - 1.);
        prodz = 1.;
        Soma = prodx + prody;
	}
	else {
		arc = F*(RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1]) + (x[2]-CCircle[2])*(x[2]-CCircle[2])));
		prody = x[1]*(x[1] - 1.);
		prodz = x[2]*(x[2] - 1.);
        Soma = prodx*prody + prodx*prodz + prody*prodz;
	}
    Prod = prodx*prody*prodz;
	den = 1. + arc*arc;
	force[0] = Soma*(M_PI + 2.*atan(arc));
	force[0] += (-2.)*(F/den)*(5*dim*Prod+Soma);
	force[0] += (F*Prod*arc*(8.*arc-0.5*F)/(den*den));
	force[0] *= B;
}

void ExactSolCircle(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
	int dim = dsol.Rows();
	REAL F = 2*sqrt(ValueK);
	REAL Coeff, B;
	if(dim==1)
		Coeff = -2.;
	else if(dim==2)
		Coeff = 8.;
	else
		Coeff = -32.;
	B = Coeff/M_PI;
    REAL arc = F;
    REAL Prod, temp;
    REAL prody = 1.;
    REAL prodz = 1.;
    REAL prodx = x[0]*(x[0]-1.);
	if(dim == 1) {
		arc *= (RCircle*RCircle - (x[0] - CCircle[0])*(x[0] - CCircle[0]));
	}
	else if(dim == 2) {
		arc *= (RCircle*RCircle - ((x[0] - CCircle[0])*(x[0] - CCircle[0]) + (x[1] - CCircle[1])*(x[1] - CCircle[1])));
		prody = x[1]*(x[1]-1.);
	}
	else if(dim == 3) {
		arc *= (RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1]) + (x[2]-CCircle[2])*(x[2]-CCircle[2])));
		prody = x[1]*(x[1]-1.);
		prodz = x[2]*(x[2]-1.);
	}
	else {
		DebugStop();
	}
    Prod = prodx*prody*prodz;
    temp = M_PI + 2.*atan(arc);
    sol[0] = B*Prod*temp;
    dsol(0,0) = B*prody*prodz*(2*x[0]-1.)*(temp - ((2*F*prodx)/(1+arc*arc)));
    if(dim==2) {
        dsol(1,0) = B*prodx*prodz*(2*x[1]-1.)*(temp - ((2*F*prody)/(1+arc*arc)));
    }
    else if(dim==3) {
        dsol(2,0) = B*prodx*prody*(2*x[2]-1.)*(temp - ((2*F*prodz)/(1+arc*arc)));
    }
}
REAL PartialDerivateX(int dim,const TPZVec<REAL> &x) {
	REAL F = 2*sqrt(ValueK);
	REAL Coeff, B;
	if(dim==1)
		Coeff = -2.;
	else if(dim==2)
		Coeff = 8.;
	else
		Coeff = -32.;
	B = Coeff/M_PI;
    REAL arc=F, prody=1., prodz=1., Prod, temp;
    REAL prodx = x[0]*(x[0]-1.);
	if(dim == 1) {
		arc *= (RCircle*RCircle - (x[0]-CCircle[0])*(x[0]-CCircle[0]));
	}
	else if(dim == 2) {
		arc *= (RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1])));
		prody = x[1]*(x[1]-1.);
	}
	else if(dim == 3) {
		arc *= (RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1]) + (x[2]-CCircle[2])*(x[2]-CCircle[2])));
		prody = x[1]*(x[1]-1.);
		prodz = x[2]*(x[2]-1.);
	}
	else {
		DebugStop();
        return 0.;
	}
    Prod = prodx*prody*prodz;
    temp = M_PI + 2.*atan(arc);
    return (B*prody*prodz*(2*x[0]-1.)*(temp - ((2*F*prodx)/(1+arc*arc))));
}

REAL PartialDerivateY(int dim,const TPZVec<REAL> &x) {
	REAL F = 2*sqrt(ValueK);
	REAL Coeff, B;
	if(dim==1)
		Coeff = -2.;
	else if(dim==2)
		Coeff = 8.;
	else
		Coeff = -32.;
	B = Coeff/M_PI;
    REAL arc=F, prody=1., prodz=1., Prod, temp;
    REAL prodx = x[0]*(x[0]-1.);
	if(dim == 1) {
		arc *= (RCircle*RCircle - (x[0]-CCircle[0])*(x[0]-CCircle[0]));
	}
	else if(dim == 2) {
		arc *= (RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1])));
		prody = x[1]*(x[1]-1.);
	}
	else if(dim == 3) {
		arc *= (RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1]) + (x[2]-CCircle[2])*(x[2]-CCircle[2])));
		prody = x[1]*(x[1]-1.);
		prodz = x[2]*(x[2]-1.);
	}
	else {
		DebugStop();
        return 0.;
	}
    Prod = prodx*prody*prodz;
    temp = M_PI + 2.*atan(arc);
    return (B*prodx*prodz*(2*x[1]-1.)*(temp - ((2*F*prody)/(1+arc*arc))));
}

REAL PartialDerivateZ(int dim,const TPZVec<REAL> &x) {
	REAL F = 2*sqrt(ValueK);
	REAL Coeff, B;
	if(dim==1)
		Coeff = -2.;
	else if(dim==2)
		Coeff = 8.;
	else
		Coeff = -32.;
	B = Coeff/M_PI;
    REAL arc=F, prody=1., prodz=1., Prod, temp;
    REAL prodx = x[0]*(x[0]-1.);
	if(dim == 1) {
		arc *= (RCircle*RCircle - (x[0]-CCircle[0])*(x[0]-CCircle[0]));
	}
	else if(dim == 2) {
		arc *= (RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1])));
		prody = x[1]*(x[1]-1.);
	}
	else if(dim == 3) {
		arc *= (RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1]) + (x[2]-CCircle[2])*(x[2]-CCircle[2])));
		prody = x[1]*(x[1]-1.);
		prodz = x[2]*(x[2]-1.);
	}
	else {
		DebugStop();
        return 0.;
	}
    Prod = prodx*prody*prodz;
    temp = M_PI + 2.*atan(arc);
    return (B*prodx*prody*(2*x[2]-1.)*(temp - ((2*F*prodz)/(1+arc*arc))));
}

/** LAPLACE PROBLEM ON L-SHAPE DOMAIN (2D) */

void ExactSolLaplaceBC(const TPZVec<REAL> &x, TPZVec<STATE> &sol) {
	REAL radius = sqrt(x[0]*x[0] + x[1]*x[1]);
    REAL angle = atan2(x[1],x[0]);;
    /*    if(IsZero(x[0])) {
     if(x[1]>0)
     angle = 0.5*M_PI;
     else {
     sol[0] = 0.;
     return;
     }
     }
     else
     angle = atan(x[1]/x[0]);
     if(angle < -0.5*M_PI || angle > M_PI)
     DebugStop();*/
    sol[0] = 0.5*pow(radius,(REAL(1./3.)))*(sqrt(3.)*sin(angle/3.)+cos(angle/3.));
}
void ExactSolLaplace(const TPZVec<REAL> &x, TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol) {
	REAL radius = sqrt(x[0]*x[0] + x[1]*x[1]);
    REAL angle = 0;
    if(IsZero(x[0])) {
        if(x[1]>0)
            angle = 0.5*M_PI;
        else {
            sol[0] = 0.;
            return;
        }
    }
    else
        angle = atan(x[1]/x[0]);
    sol[0] = 0.5*pow(radius,(REAL(1./3.)))*(sqrt(3.)*sin(angle/3.)+cos(angle/3.));
    dsol.Zero();
}


/** UTILITARY TOOL */

void UniformRefine(TPZGeoMesh* gmesh,int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
	// Re-constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}

/** Detects the bigger dimension of the computational elements into cmesh to set the Model Dimension */
bool DefineModelDimension(TPZCompMesh *cmesh) {
	if(!cmesh || !cmesh->NElements()) return false;
	TPZCompEl *cel;
	int dim = -1;
	// Run over all computational elements and check its type to define the dimension of the model
	for(int i=0;i<cmesh->NElements();i++) {
		cel = cmesh->ElementVec()[i];
		if(!cel) continue;
		int type = cel->Type();
		if(!type) dim = (dim > -1) ? dim : 0;
		else if(type==1) dim = (dim > 0) ? dim : 1;
		else if(type > 1 && type < 4)
			dim = (dim > 1) ? dim : 2;
		else if(type > 3 && type < 8)
			dim = 3;
		// If exist a three dimensional element, finish
		if(dim == 3) break;
	}
	// Whether the dimension is invalid return false
	if(dim == -1) return false;
	// If dimension is valid set into the computational mesh
	cmesh->SetDimModel(dim);
	return true;
}
int DefineDimensionOverElementType(MElementType typeel) {
	int dim = 0;
	switch(typeel) {
		case EOned:
			dim = 1;
			break;
		case EQuadrilateral:
		case ETriangle:
			dim = 2;
			break;
		case ETetraedro:
		case EPrisma:
		case EPiramide:
		case ECube:
			dim = 3;
			break;
		default:
			break;
	}
	return dim;
}

void GetFilenameFromGID(MElementType typeel, std::string &nombre) {
	switch (typeel) {
		case EOned:
			nombre = "LinhaReta.dump";
			break;
		case EQuadrilateral:
			nombre = "RegionQuadrada.dump";
			break;
		case ETriangle:
			nombre = "RegionQuadradaT.dump";
			break;
		case ETetraedro:
			nombre = "RegionCuboEnTetrahedros.dump";
			break;
		case EPrisma:
			nombre = "RegionCuboEnPrismas.dump";
			break;
		case EPiramide:
			nombre = "RegionCuboEnPiramides.dump";
			break;
		case ECube:
			nombre = "RegionCuboEnHexahedros.dump";
			break;
		default:
			DebugStop();
			break;
	}
}


// Save information of the current mesh in disk
void SaveCompMesh(TPZCompMesh *cmesh, int timessave,TPZCompMesh *cmeshmodified,bool check) {
    if(!cmesh || timessave < 0) {
        std::cout << "SaveCompMesh - Bad argument: " << (void *)cmesh << " " << timessave << std::endl;
        return;
    }
#ifdef LOG4CXX
    {
        TPZFileStream fstrthis;
        std::stringstream soutthis;
        if(cmeshmodified) soutthis << (void*)cmeshmodified;
        else soutthis << (void*)cmesh;
        // Rename the computational mesh
        cmesh->SetName(soutthis.str());
        soutthis << "_" << timessave;
        std::string filenamethis("LOG/");
        filenamethis.append(soutthis.str());
        filenamethis.append(".txt");
        fstrthis.OpenWrite(filenamethis);
        
        // Renaming the geometric mesh
        std::stringstream gout;
        gout << (void*)cmesh->Reference();
        cmesh->Reference()->SetName(gout.str());
        
        // Save geometric mesh data
        int classid = cmesh->Reference()->ClassId();
        fstrthis.Write(&classid,1);   // this first data is necessary to use TPZSavable::Restore
        cmesh->Reference()->Write(fstrthis,0);
        // Save computational mesh data
        classid = cmesh->ClassId();
        fstrthis.Write(&classid,1);   // this first data is necessary to use TPZSavable::Restore
        cmesh->Write(fstrthis,0);
        // To check printing computational mesh data in file
        if(check) {
            std::string filename("Mesh_");
            filename.append(soutthis.str());
            filename.append(".txt");
            std::ofstream arq(filename.c_str());
            cmesh->Print(arq);
        }
    }
#endif
}

/////

/////   ANOTHER TESTS

// bi-dimensional problem for elasticity on square domain
int main_GID() {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for quadrilaterals and triangles
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	
	Archivo += "/Projects/CursoPZ/MiProyecto/";
	Archivo += "MiPlaca.dump";
	TPZGeoMesh *gmesh;
	time_t sttime;
	time_t endtime;
    int dim = 2;
	
	int nelem;
	REAL distance = 0.;
	bool isdefined = false;
	for(int ii=0;ii<2;ii++) {
		time (& sttime);
		// Creating geometric mesh
		gmesh = CreateGeoMesh(Archivo);
		if(!ii) {			
			// Refinando nas esquinas desejadas
			nelem=0;
			int nrefs = 5;
			TPZVec<REAL> point(3,0.);
			TPZVec<TPZVec<REAL> > points(3);
			points[0] = point;
			point[1] = -1.;
			points[1] = point;
			point[0] = 1.;
			points[2] = point;
			
			for(int i=0;i<nrefs;i++) {
				distance = 1./((i+1)*13);
				RefineGeoElements(2,gmesh,points,distance,isdefined);
			}
			// Constructing connectivities
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
		else
			// Refinamento uniforme para toda a malla
			UniformRefine(gmesh,2);
		
		// Creating computational mesh (approximation space and materials)
		int p;
		if(!ii) p = 5;
		else p = 2;
		TPZCompEl::SetgOrder(p);
		TPZCompMesh *cmesh = CreateMesh(gmesh,dim,false);
		// Colocando a menor ordem para elementos subdivididos
		nelem = 0;
		while(!ii && nelem < cmesh->NElements()-1) {
			TPZCompEl *cel = cmesh->ElementVec()[nelem++];
			if(cel && cel->Reference()->Father()) {
				if(cel->Reference()->Father()->Father())
					((TPZInterpolatedElement*)cel)->PRefine(1);
				((TPZInterpolatedElement*)cel)->PRefine(3);
			}
		}
		cmesh->AutoBuild();
		cmesh->AdjustBoundaryElements();
		cmesh->CleanUpUnconnectedNodes();
		
		// Solving linear equations
		// Initial steps
		TPZAnalysis an (cmesh);
		TPZSkylineStructMatrix strskyl(cmesh);
		an.SetStructuralMatrix(strskyl);
		// Solver (is your choose) 
		TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
		direct->SetDirect(ECholesky);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;
		
		/*
		 // Caso no simetrico
		 //	TPZFStructMatrix full(cmesh);
		 TPZBandStructMatrix full(cmesh);
		 an.SetStructuralMatrix(full);
		 an.Solution().Zero();
		 TPZStepSolver<STATE> step;
		 step.SetDirect(ELU);
		 an.SetSolver(step);
		 */
		an.Run();
		
		// Calculando o tempo que demorou para calcular em cada cenario 
		time (& endtime);
		int time_elapsed = endtime - sttime;
		std::cout << "\n\n\tHP-Adaptive Methods....step: " << ii+1 << " time elapsed " << time_elapsed << "\n\n\n";
		
		// Post processing
		TPZStack<std::string> scalarnames, vecnames;
		std::string filename;
		if(!ii) filename = "ElasticitySolutions.vtk";
		else filename += "ElasticitySolutions1.vtk";
		scalarnames.Push("POrder");
		scalarnames.Push("SigmaX");
		scalarnames.Push("SigmaY");
		scalarnames.Push("Pressure");
		scalarnames.Push("MaxStress");
		scalarnames.Push("TauXY");
		vecnames.Push("displacement");
		vecnames.Push("PrincipalStress1");
		vecnames.Push("PrincipalStress2");
		//vecnames.Push("POrder");
		an.DefineGraphMesh(2,scalarnames,vecnames,filename);
		
		an.PostProcess(0);
	}
	return 0;
}

/** Laplace equation on L-domain */
int main_LDomain() {
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for quadrilaterals and triangles
	//gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	
	time_t sttime;
	time_t endtime;
	
    // First rectangular mesh:
	// The rectangular mesh has four corners: (0,-1,0), (1,-1,0), (1,0,0) and (0,0,0)
	// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
	// Has 4 elements, 9 connects
	cout << "Generating geometric mesh bi-dimensional ...\n";
	TPZManVector<REAL> point(3,0.), pointlast(3,0.);
	TPZGeoMesh* gmesh1 = new TPZGeoMesh;
	TPZManVector<REAL> x0(3,0.), x1(3,0.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
	x0[1] = -1.; x1[0] = 1.;
	TPZManVector<int> nx(2,2);   // subdivisions in X and in Y. 
	TPZGenGrid gen1(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
	gen1.SetElementType(EQuadrilateral);       // type = 0 means rectangular elements
	gen1.Read(gmesh1,materialId);             // generating grid in gmesh
	
	// Selecting base functions on vertices
	if(anothertests) {
		// Setting Chebyshev polynomials as orthogonal sequence generating shape functions
		TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;
		sprintf(saida,"meshextrudedLeg.vtk");
		
	}
	else {
		sprintf(saida,"meshextrudedTChe.vtk");
	}
	
	int nelem;
	REAL radius = 0.;
	bool isdefined = false;
	
	for(int ii=0;ii<2;ii++) {
		// Constructing a geometric mesh
		TPZGeoMesh* gmesh = new TPZGeoMesh;
		x0[0] = -1.; x0[1] = 0.;
		x1[0] = 1.; x1[1] = 1.;
		nx[0] = 4; //nx[1] *= 2;
		TPZGenGrid gen(nx,x0,x1);
		gen.SetElementType(EQuadrilateral);
		gen.ReadAndMergeGeoMesh(gmesh,gmesh1,materialId);
		// Inserting boundary elements with associated material
		// Bottom is fixed
		point[0] = 0.; point[1] = -1;
		pointlast[0] = 1.; pointlast[1] = -1.;
		gen.SetBC(gmesh,point,pointlast,1);
		// Top boundary has vertical force applied
		point[0] = -1; point[1] = 1.;
		pointlast[0] = 1.; pointlast[1] = 1.;
		gen.SetBC(gmesh,point,pointlast,2);
		// Vertical right boundary has horizontal force applied to left
		point[0] = 1; point[1] = -1.;
		pointlast[0] = 1.; pointlast[1] = 1.;
		gen.SetBC(gmesh,point,pointlast,3);
		
		// Initializing the process
		time (& sttime);
		if(!ii) {			
			// Refinando nas esquinas desejadas
			nelem=0;
			int nrefs = 3;
			point[0] = point[1] = point[2] = 0.;
			TPZVec<TPZVec<REAL> > points(3);
			points[0] = point;
			point[1] = -1.;
			points[1] = point;
			point[0] = 1.;
			points[2] = point;
			
			for(int i=0;i<nrefs;i++) {
				RefineGeoElements(2,gmesh,points,radius,isdefined);
				radius *= .5;
			}
			// Constructing connectivities
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
		else {
			// Refinamento uniforme para toda a malla
			UniformRefine(gmesh,1);
		}
		
		// Creating computational mesh (approximation space and materials)
		int p;
		if(!ii) p = 7;
		else p = 3;
		TPZCompEl::SetgOrder(p);
		TPZCompMesh *cmesh = CreateMesh(gmesh,2,false);
		// Disminuindo a ordem p dos elementos subdivididos
		// A cada nivel disminue em uma unidade o p, mas não será menor de 1.
		nelem = 0;
		TPZGeoEl *gelem;
		while(!ii && nelem < cmesh->NElements()-1) {
			REAL pCopy = p;
			TPZCompEl *cel = cmesh->ElementVec()[nelem++];
			if(cel) {
				gelem = cel->Reference();
				while(gelem) {
					gelem = gelem->Father();
					if(gelem) {
						if(pCopy != 1) pCopy--;
						((TPZInterpolatedElement*)cel)->PRefine(pCopy);
					}
				}
			}
		}
		cmesh->AutoBuild();
		cmesh->AdjustBoundaryElements();
		cmesh->CleanUpUnconnectedNodes();
		
		// Solving linear equations
		// Initial steps
		TPZAnalysis an (cmesh);
		TPZSkylineStructMatrix strskyl(cmesh);
		an.SetStructuralMatrix(strskyl);
		// Solver (is your choose) 
		TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
		direct->SetDirect(ECholesky);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;
		
		an.Run();
		
		// Calculando o tempo que demorou para calcular em cada cenario 
		time (& endtime);
		int time_elapsed = endtime - sttime;
		std::cout << "\n\n\tHP-Adaptive Methods....step: " << ii+1 << " time elapsed " << time_elapsed << "\n\n\n";
		
		// Post processing
		TPZStack<std::string> scalarnames, vecnames;
		std::string filename;
		if(!ii) filename = "ElasticitySolutions.vtk";
		else filename += "ElasticitySolutions1.vtk";
		scalarnames.Push("POrder");
		scalarnames.Push("SigmaX");
		scalarnames.Push("SigmaY");
		scalarnames.Push("Pressure");
		scalarnames.Push("MaxStress");
		scalarnames.Push("TauXY");
		vecnames.Push("displacement");
		vecnames.Push("PrincipalStress1");
		vecnames.Push("PrincipalStress2");
		//vecnames.Push("POrder");
		an.DefineGraphMesh(2,scalarnames,vecnames,filename);
		
		an.PostProcess(1);
		
		delete cmesh;
		delete gmesh;
	}
	return 0;
}

/** Reconstrucción del gradiente utilizando la linearizacion (Taylor) de la solución para los centros de todos los elementos vecinos */
/** Formula: u(xbi,ybi,zbi) = u(xa,ya,za) + a*(xbi-xa) + b*(ybi-ya) + c*(zbi-za)  ->  donde Grad(u) ~= (a,b,c) */
/** (xa,ya,za) es el centro del elemento donde queremos aproximar o gradiente de u */
/** (xbi,ybi,zbi) son los centros de los elementos vecinos al elemento corriente por alguno de sus lados, e enumerados por i */
void GradientReconstructionByLeastSquares(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var,bool continuous) {
	int i, nstates=0;
	TPZCompEl *cel;
	int dim = cmesh->Dimension();
	for(i=0;i<cmesh->NElements();i++) {
		cel = cmesh->ElementVec()[i];
		if(cel && cel->Dimension() == dim) {
			nstates = cel->Material()->NSolutionVariables(var);
			break;
		}
	}
	
	int nelem = cmesh->NElements();
    gradients.Redim(nelem,4*dim);
	
	int k, side;
	int counter = 0;
	
	TPZStack<TPZCompElSide> neighs;
	int nneighs = 0;
	
	TPZManVector<REAL> normal(3,0.0);
	TPZManVector<REAL> centerpsi(3,0.0);
	TPZManVector<REAL> center(3,0.0), centerbeta(3,0.0);
	TPZManVector<STATE> solalfa(nstates,0.0), solbeta(nstates,0.0);
	
	TPZFMatrix<REAL> A(dim,dim);    // Linear System matrix
	TPZFMatrix<REAL> B(dim,1,0.);   // Linear System vector
	
	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> DeltaH(nneighs,dim,0.);
	TPZFMatrix<REAL> DeltaHTranspose(dim,nneighs,0.);
	TPZFMatrix<REAL> DifSol(nneighs,1,0.);
	REAL Grad;
	
	// Calculando el gradiente por elemento computacional
	for(i=0;i<nelem;i++) {
		cel = cmesh->ElementVec()[i];
		// Nada sera realizado para elementos con dimension diferente de la dimension del problema
		if(!cel || cel->Dimension()!=dim) continue;
		
		// Limpiando las matrizes
		A.Zero(); B.Zero();
		// Encontramos el centro del elemento corriente cel
		TPZGeoEl* gelalfa = cel->Reference();
		gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
		center.Fill(0.);
		gelalfa->X(centerpsi,center);
		cel->Solution(centerpsi,var,solalfa);
		
		// PREFERENCIAL PARA CASOS DE CALCULO CON FUNCIONES DISCONTINUAS - Pues utiliza los valores de la solución en los elementos vecinos
		if(!continuous) {
			neighs.Resize(0);
			// Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
            //			for(side = cel->Reference()->NCornerNodes(); side < cel->NConnects(); side++) {
			for(side = 0; side < cel->NConnects(); side++) {
				TPZCompElSide celside(cel,side);
				celside.ConnectedElementList(neighs,1,0);
			}
			nneighs = neighs.NElements();
			// si no hay vecinos continuamos con el siguiente elemento
			if(!nneighs) continue;
			// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente			
			// Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
			// y el valor de la solucion en su centro solbeta
			DeltaH.Redim(nneighs,dim);
			DeltaHTranspose.Redim(dim,nneighs);
			DifSol.Redim(nneighs,1);
			// Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
			for(int ineighs=0;ineighs<nneighs;ineighs++) {
				TPZGeoEl* gelbeta = neighs[ineighs].Element()->Reference();
				if(!gelbeta)
					DebugStop();
				centerpsi.Fill(0.0);
				centerbeta.Fill(0.0);
				gelbeta->CenterPoint(gelbeta->NSides()-1,centerpsi);
				gelbeta->X(centerpsi,centerbeta);
				gelbeta->Reference()->Solution(centerpsi,var,solbeta);
				for(k=0;k<dim;k++)
					DeltaH(ineighs,k) = centerbeta[k] - center[k];
				DifSol(ineighs,0) = solbeta[n_var] - solalfa[n_var];
			}
		}
		else {
			int nsides = cel->NConnects()-1;
			// Para cada lado calculamos los deltaH (desde el centro del elemento al centro del lado de dimension menor a él
			// y el valor de la solucion en su centro solbeta
			DeltaH.Redim(nsides,dim);
			DeltaHTranspose.Redim(dim,nsides);
			DifSol.Redim(nsides,1);
			// Procuramos todos los puntos medios de cada lado del elemento y calculamos baseados en los valores de la solucion sobre ellos
			for(side = 0; side < nsides; side++) {
				centerpsi.Fill(0.0);
				centerbeta.Fill(0.0);
				cel->Reference()->CenterPoint(side,centerpsi);
				cel->Reference()->X(centerpsi,centerbeta);
				cel->Solution(centerpsi,var,solbeta);
				for(k=0;k<dim;k++)
					DeltaH(side,k) = centerbeta[k] - center[k];
				DifSol(side,0) = solbeta[n_var] - solalfa[n_var];
				
			}
		}
		// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u) 
		DeltaH.Transpose(&DeltaHTranspose);
		B = DeltaHTranspose*DifSol;
		A = DeltaHTranspose*DeltaH;
		A.SolveDirect(B,ELU);
		
		// Normalizando el vector gradiente
		Grad = 0.0;
		for(k=0;k<dim;k++)
			Grad += (B(k,0)*B(k,0));
		// Almacenando los gradientes encontrados
		for(k=0;k<dim;k++) {
			if(!IsZero(B(k))) {
				gradients(counter,k) = B(k,0)/sqrt(Grad);
			}
			gradients(counter,dim+k) = center[k];
			if(!k) {
				REAL dudx = PartialDerivateX(dim,center);
				REAL dudy = PartialDerivateY(dim,center);
				REAL dist = sqrt(dudx*dudx + dudy*dudy);
				if(!IsZero(dist)) {
					gradients(counter,2*dim) = dudx/dist;
					gradients(counter,2*dim+1) = dudy/dist;
				}
				dist = sqrt((center[0]-0.5)*(center[0]-0.5)+(center[1]-0.5)*(center[1]-0.5));
				if(!IsZero(dist)) {
					gradients(counter,3*dim) = (0.5-center[0])/dist;
					gradients(counter,3*dim+1) = (0.5-center[1])/dist;
				}
			}
		}
		counter++;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
//  //  FROM ADAPT_HP_JORGE
////////////////////////////////////////////////////////////////////////////////////////////////

#include "pzelast3d.h"

// Global variable
int gLMax;
int NUniformRefs = 2;
//REAL alfa = M_PI/6.;
//bool anothertests = false;
int nstate = 2;

/** Printing level */
//int gPrintLevel = 0;
//bool gDebug = false;

//TPZFNMatrix<16,REAL> Rot(4,4,0.),RotInv(4,4,0.);

void Exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);

void InitializeSolver(TPZAnalysis &an);
void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh *cmesh);


/**
 * @brief This project shows the creation of a rectangular mesh (two-dimensional) and the creation of a three-dimensional cube mesh using extrude method (ExtendMesh).
 */
int main_AdaptHP(int argc, char *argv[]) {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	const int L = 4;
	gLMax = L-1;
	char saida[260];
	
	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	
	int r, dim;
	
	// Initializing a ref patterns
	//gRefDBase.InitializeAllUniformRefPatterns();
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	
	// output files
    std::ofstream convergence("conv3d.txt");
    std::ofstream out("output.txt");
	
	/** Set polynomial order */
	int p, pmax = 2;
	for(p=1;p<pmax;p++) {
		// First rectangular mesh:
		// The rectangular mesh has four corners: (0,0,0), (1,0,0), (1,1,0) and (0,1,0)
		// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
		// Has 4 elements, 9 connects and 8 bc elements
		cout << "Generating geometric mesh bi-dimensional ...\n";
		TPZGeoMesh* gmesh = new TPZGeoMesh;
		TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
		x1[2] = 0.;
		TPZManVector<int> nx(3,3);   // subdivisions in X and in Y. 
		TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
		gen.SetElementType(EQuadrilateral);       // type = 0 means rectangular elements
		gen.Read(gmesh);             // generating grid in gmesh
		
		// Applying hp adaptive techniques 2012/10/01
		if(anothertests) {
			// Setting Chebyshev polynomials as orthogonal sequence generating shape functions
			TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;
			sprintf(saida,"meshextrudedLeg.vtk");
			
		}
		else {
			sprintf(saida,"meshextrudedTChe.vtk");
		}	
		
		// Refinement of the some element	
		TPZGeoEl *gel;   //  *gel1, *gel2, *gel3;
		TPZVec<TPZGeoEl *> sub;
		TPZVec<TPZGeoEl *> subsub;
		gel = gmesh->ElementVec()[4];
		//	gel1 = gmesh->ElementVec()[1];
		//	gel2 = gmesh->ElementVec()[2];
		//	gel3 = gmesh->ElementVec()[3];
		gel->Divide(sub);
		//	sub[0]->Divide(subsub);
		//		sub[1]->Divide(subsub);
		sub[2]->Divide(subsub);
		//	sub[3]->Divide(subsub);
		/*	gel1->Divide(sub);
		 sub[0]->Divide(subsub);
		 sub[1]->Divide(subsub);
		 sub[2]->Divide(subsub);
		 sub[3]->Divide(subsub);
		 gel2->Divide(sub);
		 sub[0]->Divide(subsub);
		 sub[1]->Divide(subsub);
		 sub[2]->Divide(subsub);
		 sub[3]->Divide(subsub);
		 gel3->Divide(sub);
		 sub[0]->Divide(subsub);
		 sub[1]->Divide(subsub);
		 sub[2]->Divide(subsub);
		 sub[3]->Divide(subsub);
		 */	
		// Constructing connectivities
		gmesh->ResetConnectivities();
		gmesh->BuildConnectivity();
		gmesh->Print();
		// Printing COMPLETE initial geometric mesh 
		PrintGeoMeshInVTKWithDimensionAsData(gmesh,saida);
		
		TPZCompEl::SetgOrder(p);
		// Creating computational mesh
		TPZCompMesh *comp = new TPZCompMesh(gmesh);
		
		// Creating and inserting materials into computational mesh
		TPZMaterial * mat = new TPZElasticityMaterial(1,1.e5,0.2,.5,0);   // two-dimensional
		comp->InsertMaterialObject(mat);
		dim = mat->Dimension();
		nstate = mat->NStateVariables();
		
		// Boundary conditions
		// Dirichlet
		TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,5.);
		val1(0,0) = 1.;
		TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
		comp->InsertMaterialObject(bnd);
		// Neumann
		val2(0,0)=30.; val2(1,0) = 10.;
		bnd = mat->CreateBC(mat,-2,1,val1,val2);
		comp->InsertMaterialObject(bnd);
		
		// Constructing and adjusting computational mesh
		comp->AutoBuild();
		comp->AdjustBoundaryElements();   // Adjust boundary elements and higher level of refinement, clean elements but not connects into them
		comp->CleanUpUnconnectedNodes();  // Clean connects not connected at least one element enabled.
		comp->Print();
		
		
		//--- END construction of the meshes
		/** Variable names for post processing */
		TPZStack<std::string> scalnames, vecnames;
		if(mat->NSolutionVariables(mat->VariableIndex("POrder")) == 1)
			scalnames.Push("POrder");
		else
			vecnames.Push("POrder");
		if(mat->NSolutionVariables(mat->VariableIndex("Error")) == 1)
			scalnames.Push("Error");
		else
			vecnames.Push("Error");
		if(mat->NSolutionVariables(mat->VariableIndex("state")) == 1)
			scalnames.Push("state");
		else
			vecnames.Push("state");
		
		if(nstate == 1) {
			scalnames.Push("TrueError");
			scalnames.Push("EffectivityIndex");
		}else if(nstate == 2) {
			scalnames.Push("sig_x");
			scalnames.Push("sig_y");
			scalnames.Push("tau_xy");
		}
		if(nstate == 3) {
			scalnames.Push("StressX");
			scalnames.Push("StressY");
			scalnames.Push("StressZ");
			vecnames.Push("PrincipalStress");
			vecnames.Push("PrincipalStrain");
		}
		// END Determining the name of the variables
		
		// INITIAL POINT FOR SOLVING AND APPLYING REFINEMENT
		for(r=0;r<NUniformRefs;r++) {
			// Printing computational mesh to information
			if(comp->NElements() < 200)
				comp->Print(std::cout);
			else {
				std::cout << "Computacional mesh : NElements = " << comp->NElements() << "\t NConnects = " << comp->NConnects() << std::endl;
			}
			
			// Introduzing exact solution depending on the case
			TPZAnalysis an (comp);
			an.SetExact(Exact);		
			{   // To print solution
				std::stringstream sout;
				int angle = (int) (alfa*180./M_PI + 0.5);
				if(anothertests) sout << "Leg_";
				sout << "hptestAngo" << angle << "." << r << ".vtk";
				an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
			}
			std::string MeshFileName;
			{   // To print computational mesh
				std::stringstream sout;
				int angle = (int) (alfa*180./M_PI + 0.5);
				if(anothertests) sout << "Leg_";
				sout << "meshAngle" << angle << "." << r << ".vtk";
				MeshFileName = sout.str();
			}
			comp->SetName("Malha computacional adaptada");
			
			// Solve using symmetric matrix then using Cholesky (direct method)
			TPZSkylineStructMatrix strskyl(comp);
			an.SetStructuralMatrix(strskyl);
			
			TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
			direct->SetDirect(ECholesky);
			an.SetSolver(*direct);
			delete direct;
			direct = 0;
			
            //			int neq = comp->NEquations();
			//	an.NEquations();
			//		an.Solution().Print();
			an.Run();
			
			// Computing approximation of gradient
			/** 
			 * @brief Method to reconstruct a gradient after run Solve of the analysis
			 * @param cmesh Computational mesh with solution */
			TPZFMatrix<REAL> gradients;
			GradientReconstructionByLeastSquares(gradients,comp,0,0,true);
			gradients.Print("gradients");
			
			// Post processing
			an.PostProcess(1,dim);
			{
				std::ofstream out(MeshFileName.c_str());
				comp->LoadReferences();
				TPZVTKGeoMesh::PrintGMeshVTK(comp->Reference(), out, false);
			}
		}
	}
	return 0;
}

int main_AdaptHP_3D(int argc, char *argv[]) {
	
#ifdef LOG4CXX
	if (argc > 1) {
		std::string logpath ( argv[1] );
		cout << "initializing LOG usign the following configuration file " << logpath << endl;
		InitializePZLOG ( logpath );	
	} else {
		cout << "initializing LOG\n";
		InitializePZLOG();
	}
#endif
	const int L = 4;
	gLMax = L-1;
	char saida[260];
	
	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	
	int r, dim;
	
	// Initializing a ref patterns
	gRefDBase.InitializeAllUniformRefPatterns();
	//gRefDBase.InitializeRefPatterns();
	//gRefDBase.InitializeUniformRefPattern(EOned);
	//gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	//gRefDBase.InitializeUniformRefPattern(ETriangle);
	// Inserting a special file with refinement pattern 
	std::string filename = REFPATTERNDIR;
	filename += "/3D_Hexa_Rib_Side_16_16_18_18.rpt";
	//filename += "/3D_Hexa_Face_20.rpt";
	
	TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
	if(!gRefDBase.FindRefPattern(refpat))
	{
		gRefDBase.InsertRefPattern(refpat);
	}
	refpat->InsertPermuted();
	
	// output files
    std::ofstream convergence("conv3d.txt");
    std::ofstream out("output.txt");
	
    // First rectangular mesh:
	// The rectangular mesh has four corners: (0,0,0), (1,0,0), (1,1,0) and (0,1,0)
	// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
	// Has 4 elements, 9 connects and 8 bc elements
	cout << "Generating geometric mesh bi-dimensional ...\n";
    TPZGeoMesh* gmesh = new TPZGeoMesh;
	TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
	TPZManVector<int> nx(2,2);   // subdivisions in X and in Y. 
	TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
	gen.SetElementType(EQuadrilateral);       // type = 0 means rectangular elements
	gen.Read(gmesh);             // generating grid in gmesh
	
	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	cout << "Generating geometric mesh three-dimensional (extruding) ...\n";
	TPZExtendGridDimension gmeshextend(gmesh,0.5);
	TPZGeoMesh *gmesh3D = gmeshextend.ExtendedMesh(2,2,2);
	
	// Applying hp adaptive techniques 2012/10/01
	if(anothertests) {
		// Setting Chebyshev polynomials as orthogonal sequence generating shape functions
		TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;
		sprintf(saida,"meshextrudedLeg.vtk");
		
	}
	else {
		sprintf(saida,"meshextrudedTChe.vtk");
	}	
	// Uniform Refinement - Some times for three dimensional elements
	//	UniformRefinement(2,gmesh3D,3);
	// Refinement of the some element	
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	TPZVec<TPZGeoEl *> subsub;
	int nele;
	//	for(int ii=0;ii<3;ii++) {
	//		int ngelem = gmesh->NElements()-1;
	nele = 0;
	//		for(;nele<ngelem;nele++) {
	gel = gmesh3D->ElementVec()[nele];
	
	//			if(gel->Dimension() != 3) continue;
	//	gel->SetRefPattern(refpat);
	
	gel->Divide(sub);
	//			int jj = 0;
	//			for(jj=0;jj<4;jj++) {
	//				gel = sub[jj];
	//				gel->Divide(subsub);
	//			}
	//		}
	//			TPZVec<REAL> coord(3,0.);
	//			TPZVec<REAL> point(3,-1.);
	//			gel->X(point,coord);
	//			if(!IsZero(coord[0]) || !IsZero(coord[1]) || !IsZero(coord[2])) continue;
	//			if(gel->Dimension() != 3) continue;
	//			gel->SetRefPattern(refpat);
	gel->Divide(sub);
	//			for(int jj=0;jj<4;jj++) {
	//				gel = sub[jj];
	//				if(gel->Dimension() != 3) continue;
	//				gel->SetRefPattern(refpat);
	//				gel->Divide(subsub);
	//				for(int kk=0;kk<gel->NSubElements();kk++)
	//					subsub[kk]->SetRefPattern(refpat);
	//			}
	//			gel = subsub[0];
	//		gel->Divide(sub);
	//		gel = sub[0];
	//		gel->Divide(subsub);
	//		gel = subsub[0];
	//		gel->Divide(sub);
	//			break;
	//		}
	//	}
	// Constructing connectivities
	gmesh3D->ResetConnectivities();
	gmesh3D->BuildConnectivity();
	gmesh3D->Print();
	// Printing COMPLETE initial geometric mesh 
	PrintGeoMeshInVTKWithDimensionAsData(gmesh3D,saida);
	
    // Creating computational mesh
    TPZCompMesh *comp = new TPZCompMesh(gmesh3D);
	/** Set polynomial order */
	int p = 2;
    TPZCompEl::SetgOrder(p);
	
	TPZVec<STATE> forces(3,0.);
	// Creating and inserting materials into computational mesh
    //TPZMaterial * mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);   // two-dimensional
	TPZMaterial *mat = new TPZElasticity3D(1,1.e5,0.2,forces);          // three-dimensional
    comp->InsertMaterialObject(mat);
	dim = mat->Dimension();
	nstate = mat->NStateVariables();
	
    // Boundary conditions
    // Dirichlet
    TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,5.);
	val1(0,0) = val1(1,1) = 1.;
    TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
    comp->InsertMaterialObject(bnd);
	// Neumann
    val2(0,0)=30.; val2(1,0) = 10.;
    bnd = mat->CreateBC(mat,-2,1,val1,val2);
    comp->InsertMaterialObject(bnd);
    
    // Constructing and adjusting computational mesh
    comp->AutoBuild();
    comp->AdjustBoundaryElements();   // Adjust boundary elements and higher level of refinement, clean elements but not connects into them
    comp->CleanUpUnconnectedNodes();  // Clean connects not connected at least one element enabled.
	
	//--- END construction of the meshes
	
	/** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
	if(mat->NSolutionVariables(mat->VariableIndex("POrder")) == 1)
		scalnames.Push("POrder");
	else
		vecnames.Push("POrder");
	if(mat->NSolutionVariables(mat->VariableIndex("Error")) == 1)
		scalnames.Push("Error");
	else
		vecnames.Push("Error");
	if(mat->NSolutionVariables(mat->VariableIndex("state")) == 1)
		scalnames.Push("state");
	else
		vecnames.Push("state");
    
    if(nstate == 1) {
        scalnames.Push("TrueError");
        scalnames.Push("EffectivityIndex");
    }else if(nstate == 2) {
        scalnames.Push("sig_x");
        scalnames.Push("sig_y");
        scalnames.Push("tau_xy");
    }
    if(nstate == 3) {
        scalnames.Push("StressX");
        scalnames.Push("StressY");
        scalnames.Push("StressZ");
		vecnames.Push("PrincipalStress");
		vecnames.Push("PrincipalStrain");
    }
	
	// END Determining the name of the variables
	
	// TO MAKE MERGE ANOTHER DOMAIN
	if(anothertests) {
		// Second rectangular domain - subdivisions and corners of the second rectangular mesh
		TPZAutoPointer<TPZGeoMesh> gmesh2 = new TPZGeoMesh;
		x0[1] = 0.2;                 // left and right extremes of the new geo mesh. Coordinates: (0.,0.2,0.0) (3.,1.,0.) 
		x1[0] = 3.; x1[1] = 1.;
		nx[0] = 15; nx[1] = 8;       // subdivision in X and Y. hx = 0.2 and hy = 0.1
		TPZGenGrid gen2(nx,x0,x1);   // second mesh generator
		gen2.SetElementType(EQuadrilateral);      // type = 0 means rectangular elements, type = 1 means triangular elements
		
		// Generating gmesh2 with last data and after this the gmesh is merged into the gmesh2. But gmesh is unmodified
		// if exist boundary elements into the mesh merged it will be deleted
		gen2.ReadAndMergeGeoMesh(gmesh2,gmesh);
		
		// setting bc condition -1 [no flux - is wall] from (0.0, 0.0) until (3.0, 0.2)
		x0[1] = 0.0; x1[1] = 0.2;
		gen2.SetBC(gmesh2,x0,x1,-1);
		// setting bc condition -1 from (3.0, 1.0) until (0.0, 1.0)
		x0[0] = 3.; x0[1] = 1.;
		x1[0] = 0.;	x1[1] = 1.;
		gen2.SetBC(gmesh2,x0,x1,-1);
		// setting bc condition -2 [free flux] from (3.0, 1.0) until (3.0, 0.2)
		x1[0] = 3.;	x1[1] = .2;
		gen2.SetBC(gmesh2,x1,x0,-2);
		// setting bc condition -2 [free flux] from (0.0, 1.0) until (0.0, 0.0)
		x0[0] = 0.;
		x1[0] = x1[1] = 0.;
		gen2.SetBC(gmesh2, x0, x1, -2);
		
#ifdef PZDEBUG
		sprintf(saida,"original.vtk");
		PrintGeoMeshInVTKWithDimensionAsData(gmesh,saida);
		sprintf(saida,"meshes.vtk");
		PrintGeoMeshInVTKWithDimensionAsData(gmesh2.operator->(),saida);
#endif
		// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
		// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
		TPZExtendGridDimension gmeshextend(gmesh2,0.3);
		TPZGeoMesh *gmesh3D = gmeshextend.ExtendedMesh(2,-5,-6);
		if(gmesh3D && gmesh3D->NElements() < 400) gmesh3D->Print();
#ifdef PZDEBUG
		sprintf(saida,"meshextrudedend.vtk");
		PrintGeoMeshInVTKWithDimensionAsData(gmesh3D,saida);
#endif
		dim = 3;
	}
	
	// INITIAL POINT FOR SOLVING AND APPLYING REFINEMENT
	for(r=0;r<NUniformRefs;r++) {
		// Printing computational mesh to information
		if(comp->NElements() < 200)
			comp->Print(std::cout);
		else {
			std::cout << "Computacional mesh : NElements = " << comp->NElements() << "\t NConnects = " << comp->NConnects() << std::endl;
		}
		
		// Introduzing exact solution depending on the case
		TPZAnalysis an (comp);
		an.SetExact(Exact);		
		{   // To print solution
			std::stringstream sout;
			int angle = (int) (alfa*180./M_PI + 0.5);
			if(anothertests) sout << "Leg_";
			sout << "hptestAngo" << angle << "." << r << ".vtk";
			an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
		}
		std::string MeshFileName;
		{   // To print computational mesh
			std::stringstream sout;
			int angle = (int) (alfa*180./M_PI + 0.5);
			if(anothertests) sout << "Leg_";
			sout << "meshAngle" << angle << "." << r << ".vtk";
			MeshFileName = sout.str();
		}
		comp->SetName("Malha computacional adaptada");
		
		// Solve using symmetric matrix then using Cholesky (direct method)
		TPZSkylineStructMatrix strskyl(comp);
		an.SetStructuralMatrix(strskyl);
		
		TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
		direct->SetDirect(ECholesky);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;
		
		an.Run();
		
		// Post processing
		an.PostProcess(1,dim);
		{
			std::ofstream out(MeshFileName.c_str());
			comp->LoadReferences();
			TPZVTKGeoMesh::PrintGMeshVTK(comp->Reference(), out, false);
		}
		
		/*
		 REAL valerror =0.;
		 REAL valtruerror=0.;
		 TPZVec<REAL> ervec,truervec,effect;
		 
		 TPZAdaptMesh adapt;
		 adapt.SetCompMesh (comp);
		 
		 std::cout << "\n\n\n\nEntering Auto Adaptive Methods... step " << r << "\n\n\n\n";
		 
		 time_t sttime;
		 time (& sttime);
		 TPZCompMesh *adptmesh;
		 
		 adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror,ervec,Exact,truervec,effect,0);
		 
		 time_t endtime;
		 time (& endtime);
		 
		 int time_elapsed = endtime - sttime;
		 std::cout << "\n\n\n\nExiting Auto Adaptive Methods....step " << r
		 << "time elapsed " << time_elapsed << "\n\n\n\n";
		 
		 int prt;
		 std::cout << "neq = " << comp->NEquations() << " error estimate = " << valerror
		 << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
		 
		 #ifdef LOG4CXX
		 if (loggerconv->isDebugEnabled())
		 {
		 std::stringstream sout;
		 sout << "neq = " << comp->NEquations() << " error estimate = " << valerror
		 << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
		 LOGPZ_DEBUG(loggerconv, sout.str())
		 }
		 #endif
		 
		 convergence  << comp->NEquations() << "\t"
		 << valerror << "\t"
		 << valtruerror << "\t"
		 << ( valtruerror / valerror ) <<  "\t"
		 << sttime <<std::endl;
		 for (prt=0;prt<ervec.NElements();prt++){
		 std::cout <<"error " << ervec[prt] << "  truerror = " << truervec[prt] << "  Effect " << effect[prt] << std::endl;
		 // convergence << '\t' << ervec[prt] << '\t' << truervec[prt] << "  Effect " << effect[prt] <<  std::endl;
		 //  adptmesh->Print(cout);
		 }
		 
		 std::cout.flush();
		 comp->Reference()->ResetReference();
		 comp->LoadReferences();
		 adapt.DeleteElements(comp);
		 delete comp;
		 comp = adptmesh;
		 
		 comp->CleanUpUnconnectedNodes();
		 */
	}
	/* Uniform refinement. Two times
	 UniformRefinement(2,gmesh3D,3);
	 
	 #ifdef PZDEBUG
	 sprintf(saida,"meshrefined.vtk");
	 PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
	 #endif
	 
	 // Creating a computational mesh (interpolation space)
	 TPZCompMesh *cmesh = CreateMeshMultires(gmesh3D);
	 #ifdef PZDEBUG
	 sprintf(saida,"aftercmesh.vtk");
	 PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
	 #endif
	 
	 // To work with the temporal variable.
	 REAL timeStep;
	 timeStep = ComputeTimeStep(1.,L,cmesh->Reference());
	 
	 #ifdef PZDEBUG
	 {
	 ofstream malhas("malhas.vtk");
	 cmesh->Print(malhas);
	 }
	 #endif
	 
	 TPZExplFinVolAnal an(cmesh, cout);
	 
	 InitializeSolver(an);
	 const double PhysicalTime = 0.1;
	 int niter = PhysicalTime/timeStep+1;
	 cout << "L = " << L << endl;
	 cout << "\nnequations = " << cmesh->NEquations();
	 cout << "\nNiter = " << niter << "\n";
	 
	 TPZFMatrix<REAL> InitialSol;
	 InitialSolutionLinearConvection(InitialSol,cmesh);
	 
	 an.SetInitialSolution(InitialSol);
	 
	 an.Set(timeStep,niter,1e-10);
	 an.SetSaveFrequency(niter/6,0);
	 
	 double Epsl = 1.e12;
	 an.MultiResolution( Epsl );
	 #ifdef PZDEBUG
	 sprintf(saida,"meshInitialSol.vtk");
	 TPZVTKGeoMesh::PrintGMeshVTK(gmesh3D,saida,0);
	 #endif
	 
	 return EXIT_SUCCESS;
	 */
	return 0;
}

/** Exact solutions to calculate the rate of convergence */

static REAL onethird = 0.33333333333333333;
static REAL PI = 3.141592654;

void Exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
    TPZManVector<REAL,3> x2(x);
	//    TransformInvX(x2,RotInv);
  	REAL r = sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
  	REAL theta = atan2(x2[1],x2[0]);
  	REAL rexp = pow(r,(REAL)onethird);
  	sol[0] = rexp*sin(onethird*(theta+PI/2));
    TPZFNMatrix<3,REAL> grad(4,1,0.),grad2(4,1,0.);
  	grad(0,0) = onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
  	grad(1,0) = onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
	//    Rot.Multiply(grad, grad2);
    dsol(0,0) = grad2(0,0);
    dsol(1,0) = grad2(1,0);
}

// Printing in VTK format the geometric mesh but taking geometric elements as reference or computational elements as reference
void PrintGeoMeshInVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename) {
	int i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		if(gmesh->ElementVec()[i])
			DataElement[i] = (gmesh->ElementVec()[i])->Dimension();
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}

void PrintGeoMeshAsCompMeshInVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename) {
	int i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		TPZGeoEl *gel = gmesh->ElementVec()[i];
		if(gel && gel->Reference())
			DataElement[i] = gel->Dimension();
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}
void PrintGeoMeshAsCompMeshInVTKWithElementIDAsData(TPZGeoMesh *gmesh,char *filename) {
	int i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		TPZGeoEl *gel = gmesh->ElementVec()[i];
		if(gel && gel->Reference())
			DataElement[i] = gel->Id();
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}
void PrintGeoMeshAsCompMeshInVTKWithElementIndexAsData(TPZGeoMesh *gmesh,char *filename) {
	int i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		TPZGeoEl *gel = gmesh->ElementVec()[i];
		if(gel && gel->Reference())
			DataElement[i] = gel->Id();
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}

/*
void PrintGeoMeshAsCompMeshInVTKWithElementData(TPZGeoMesh *gmesh,char *filename,TPZVec<REAL> &elData) {
	int i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		TPZGeoEl *gel = gmesh->ElementVec()[i];
		if(gel && gel->Reference())
			DataElement[i] = elData[gel->Reference()->Index()];
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}
*/
void InitializeSolver(TPZAnalysis &an) {
	TPZStepSolver<STATE> step;
	TPZBandStructMatrix matrix(an.Mesh());
	an.SetStructuralMatrix(matrix);
	step.SetDirect(ELU);
	an.SetSolver(step);
}

void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh * cmesh){
	InitialSol.Redim(cmesh->NEquations(),1);
	InitialSol.Zero();
	for(int iel = 0; iel < cmesh->NElements(); iel++){
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
		if(!disc) continue;
		if(disc->NConnects() == 0) continue;
		int bl = disc->Connect(0).SequenceNumber();
		int blpos = cmesh->Block().Position(bl);
		int blocksize = cmesh->Block().Size(bl);
		
		TPZGeoEl * gel = cel->Reference();
		TPZVec<REAL> xi(3), xVec(3);
		gel->CenterPoint(gel->NSides()-1,xi);
		gel->X(xi,xVec);
		double x = xVec[0];
		double y = xVec[1];
		double u = 0.125;
		
		double xCircle = 0.25;
		double yCircle = 0.5;
		double R = 0.1;
		if( (x-xCircle)*(x-xCircle)+(y-yCircle)*(y-yCircle) <= R*R ) u = 1.;
		
		InitialSol(blpos+blocksize-20+0,0) = u;
		InitialSol(blpos+blocksize-20+1,0) = 0.;
		InitialSol(blpos+blocksize-20+2,0) = 0.;
		InitialSol(blpos+blocksize-20+3,0) = 0.;
		InitialSol(blpos+blocksize-20+4,0) = 0.;
		
	}//for iel
	
	TPZVec<REAL> celerity(3,0.);
	celerity[0] = 1.;
#ifdef LinearConvection
	TPZEulerEquation::SetLinearConvection(cmesh, celerity);
#endif
	
}//method

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial, const int matidtodivided) {
	TPZManVector<TPZGeoEl*> filhos;
    for(int D=0; D<nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {    
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
			if(!gel || gel->HasSubElement())
				continue;
			if(dim > 0 && gel->Dimension() != dim) continue;
			if(!allmaterial){
				if(gel->MaterialId() == matidtodivided) {
					gel->Divide(filhos);
				}
			}
			else{
				gel->Divide(filhos);
			}
        }
    }
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
void formatTimeInSec(char *strtime,int lenstrtime,int timeinsec) {
	if(!strtime) return;
	memset(strtime,0,strlen(strtime));

	int anos=0, meses=0, dias=0, horas=0, minutos=0, segundos=0;
	while(1) {
		if(timeinsec < 60) {
			segundos = timeinsec;
			break;
		}
		else {
			timeinsec -= 60;
			minutos++;
			if(minutos > 59) {
				minutos -= 60;
				horas++;
				if(horas > 23) {
					horas -= 24;
					dias++;
					if(dias > 29) {
						dias -= 30;
						meses++;
						if(meses > 11) {
							meses -= 12;
							anos++;
						}
					}
				}
			}
		}
	}
	// Formating
	if(anos)
#ifdef WIN32
		sprintf_s(strtime,lenstrtime,"%d a, %d m, %d d, %02d:%02d:%02d",anos,meses,dias,horas,minutos,segundos);
#else
		sprintf(strtime,"%d a, %d m, %d d, %02d:%02d:%02d",anos,meses,dias,horas,minutos,segundos);
#endif
	else {
		if(meses) 
#ifdef WIN32
			sprintf_s(strtime,lenstrtime,"%d m, %d d, %02d:%02d:%02d",meses,dias,horas,minutos,segundos);
#else
			sprintf(strtime,"%d m, %d d, %02d:%02d:%02d",meses,dias,horas,minutos,segundos);
#endif
		else {
			if(dias)
#ifdef WIN32
				sprintf_s(strtime,lenstrtime,"%d d, %02d:%02d:%02d",dias,horas,minutos,segundos);
#else
				sprintf(strtime,"%d d, %02d:%02d:%02d",dias,horas,minutos,segundos);
#endif
			else
#ifdef WIN32
				sprintf_s(strtime,lenstrtime,"%02d:%02d:%02d",horas,minutos,segundos);
#else
				sprintf(strtime,"%02d:%02d:%02d",horas,minutos,segundos);
#endif
		}
	}
}

///// TO MANUAL REFINEMENT OVER RADIUS VALUE  //////////////////////////
////////////////////////////////////////////////////////////////////////

void DeterminingPOrderOnLevelHRefinement(TPZCompMesh *cmesh,int p) {
	int level = 0, highlevel = 0;
	int pinit;
	int nelem = 0;
	while(nelem < cmesh->NElements()) {
		TPZCompEl *cel = cmesh->ElementVec()[nelem++];
		if(cel) {
			level = cel->Reference()->Level();
		}
		if(level > highlevel)
			highlevel = level;
	}
	// Identifying maxime interpolation order
	if(highlevel>p-1) pinit = p;
	else pinit = highlevel+1;
	// Put order 1 for more refined element and (highlevel - level)+1 for others, but order not is greater than initial p
	nelem = 0;
	while(highlevel && nelem < cmesh->NElements()) {
		TPZCompEl *cel = cmesh->ElementVec()[nelem++];
		if(!cel) continue;
		level = cel->Reference()->Level();
		p = (highlevel-level);
		if(!p) p = 1;     // Fazendo os dois maiores niveis de refinamento devem ter ordem 1
		if(p > pinit) p = pinit;
		((TPZInterpolatedElement*)cel)->PRefine(p);
	}
	cmesh->ExpandSolution();
	cmesh->CleanUpUnconnectedNodes();
}

void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points) {
	Points.Resize(npoints);
	TPZManVector<REAL> point(3,0.);
	REAL angle = (2*M_PI)/npoints;
	for(int i=0;i<npoints;i++) {
		point[0] = center[0]+radius*cos(i*angle);
		point[1] = center[1]+radius*sin(i*angle);
		Points[i] = point;
	}
}
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,REAL &radius,int ntyperefs) {
	TPZVec<REAL> point(3);
	point[0] = point[1] = 0.5; point[2] = 0.0;
	REAL r = 0.25;
	bool isdefined = true;
	
	if(ntyperefs==2) {
		// To refine elements with center near to points than radius
		RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
		RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
	}
	else {
		// To refine elements with center near to points than radius
		RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
	}
	// Constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
    // reducing the distance from circunference for next iteration
    if(ntyperefs == 2) radius *= 0.5;
    else radius *= 0.7;
}

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			REAL distancia = TPZGeoEl::Distance(center,FirstNode);
			if(distancia > distance) distance = distancia;
			isdefined = true;
		}
		REAL centerdist = TPZGeoEl::Distance(center,point);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
}

void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
	int i;
	bool isdefined = false;
	
	// Refinando no local desejado
	TPZVec<REAL> point(3);
	point[0] = point[1] = 0.5; point[2] = 0.0;
	REAL r = 0.25;
	
	if(ntyperefs==2) {
		REAL radius = 0.19;
		for(i=0;i<nref;i+=2) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			if(nref < 5) radius *= 0.35;
			else if(nref < 7) radius *= 0.2;
			else radius *= 0.1;
		}
		if(i==nref) {
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
		}
	}
	else {
		REAL radius = 0.2;
		for(i=0;i<nref+1;i++) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			if(nref < 5) radius *= 0.6;
			else if(nref < 7) radius *= 0.3;
			else radius *= 0.15;
		}
	}
	// Constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<TPZVec<REAL> > &points,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	int i, npoints = points.NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			distance = 1.1*TPZGeoEl::Distance(center,FirstNode);
			isdefined = true;
		}
		for(i=0;i<npoints;i++) {
			REAL semidiag = TPZGeoEl::Distance(center,points[i]);
			if(semidiag < distance) {
				gel->Divide(sub);
				break;
			}
		}
	}
}

// MAIN FUNCTION TO NUMERICAL SOLVE WITH HP ADAPTIVE REFINEMENTS BUT IT HASN´T AUTOMATIC STRATEGY TO ADAPT
/** Laplace equation on square 1D 2D 3D - Volker John article 2000 */
/*
int main_NoAutoHP() {

#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();
 //   gRefDBase.InitializeRefPatterns();

	// To compute processing times
	time_t sttime;
	time_t endtime;
	int time_elapsed;
	const int lenstrtime = 512;
	char time_formated[lenstrtime];
	
	// Output files
    std::ofstream convergence("convergence.txt");
	std::ofstream fileerrors("ErrorsHP_ArcTan.txt");   // To store all errors calculated by TPZAnalysis (PosProcess)
	
	// To compute the errors
	TPZManVector<REAL> ervec(100,0.0);
	// Printing computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	int nref = 1, NRefs = 5;
	int nthread = 1, NThreads = 2;
    int dim;
	
    //Working on regular meshes
    for(int regular=1; regular>0; regular--) {
		fileerrors << "Type of mesh: " << regular << " Level. " << endl;
		MElementType typeel;
		for(int itypeel=(int)EOned;itypeel<(int)EPolygonal;itypeel++)
//		for(int itypeel=(int)EPrisma;itypeel<(int)EPolygonal;itypeel++)
		{
			typeel = (MElementType)itypeel;
			fileerrors << "Type of element: " << typeel << endl;
			TPZGeoMesh *gmesh;
			if(!regular) {
				std::string nombre;
				GetFilenameFromGID(typeel,nombre);
				// Generating geometric mesh 
				gmesh = CreateGeoMesh(nombre);
			}
			else {
				gmesh = CreateGeoMesh(typeel);
			}
			dim = DefineDimensionOverElementType(typeel);
			
			// Some refinements as initial step
			UniformRefinement(1,gmesh,dim);

			// Creating computational mesh (approximation space and materials)
			int p = 3, pinit;
			pinit = p;
			TPZCompEl::SetgOrder(p);
			TPZCompMesh *cmesh = CreateMesh(gmesh,dim,1);
			gmesh->SetName("Malha Geometrica original");
			cmesh->SetName("Malha Computacional Original");
			
			// Selecting orthogonal polynomial family to construct shape functions
			if(anothertests)
				TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;  // Setting Chebyshev polynomials as orthogonal sequence generating shape functions
			
			// Variable names for post processing 
			TPZStack<std::string> scalnames, vecnames;
			scalnames.Push("POrder");
			scalnames.Push("Solution");
			
			// Solving adaptive process
			for(nref=1;nref<NRefs;nref++) {
				cout << "\nConstructing Poisson problem " << dim << "D. Refinement: " << nref << " Threads: " << nthread << " Regular: " << regular << " TypeElement: " << typeel << endl;
				if(nref > 5) nthread = 2*NThreads;
				else nthread = NThreads;
				// Initializing the generation mesh process
				time (& sttime);
				
				// Introduzing exact solution depending on the case
				TPZAnalysis an(cmesh);
				an.SetExact(ExactSolCircle);
				{
					std::stringstream sout;
					int angle = (int) (alfa*180./M_PI + 0.5);
					sout << "Poisson" << dim << "D_MESH" << regular << "E" << typeel << "Thr" << nthread << "H" << nref << "P" << pinit << "_Ang" << angle << ".vtk";
					an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
				}
				std::string MeshFileName;
				{
					std::stringstream sout;
					int angle = (int) (alfa*180./M_PI + 0.5);
					sout << "meshAngle" << dim << "D_MESH" << regular << "E" << typeel << "Thr" << nthread << "H" << nref << "P" << pinit << "_Ang" << angle << ".vtk";
					MeshFileName = sout.str();
				}
				
				cmesh->SetName("Malha computacional adaptada");
				
				// Solve using symmetric matrix then using Cholesky (direct method)
			//	TPZParSkylineStructMatrix strskyl(cmesh,nthread);
				TPZSkylineStructMatrix strskyl(cmesh);
				an.SetStructuralMatrix(strskyl);
				
				TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
				direct->SetDirect(ECholesky);
				an.SetSolver(*direct);
				delete direct;
				direct = 0;
				
				an.Run();
				// Post processing
				an.PostProcess(0,dim);

				// generation mesh process finished
				
				time(&endtime);
				time_elapsed = endtime - sttime;
				formatTimeInSec(time_formated,lenstrtime,time_elapsed);
				out << "\tRefinement: " << nref+1 << " Regular Mesh: " << regular << " TypeElement: " << typeel << " Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n\n";
								
				std::cout << "\n\n\nEntering Auto Adaptive Methods... step " << nref << "\n\n";
				
				if(NRefs>1) {
					time (& sttime);

					REAL radius = 0.2;
			
					// h_refinement
					// Refining near the points belong a circunference with radio r - maxime distance radius
					RefiningNearCircunference(dim,gmesh,radius,1);
					radius *= 0.6;
					
					time_t endtime;
					time (& endtime);
					
					int time_elapsed = endtime - sttime;
					std::cout << "\n\nExiting Auto Adaptive Methods....step " << nref << "time elapsed " << time_elapsed << "\n\n\n\n";
					
//					int prt;
	//				std::cout << "neq = " << cmesh->NEquations() << " error estimate = " << valerror << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
					
		//			convergence  << cmesh->NEquations() << "\t" << valerror << "\t" << valtruerror << "\t" << ( valtruerror / valerror ) <<  "\t" << sttime <<std::endl;
			//		for (prt=0;prt<ervec.NElements();prt++) {
				//		std::cout <<"error " << ervec[prt] << "  truerror = " << truervec[prt] << "  Effect " << effect[prt] << std::endl;
					//}
				}
				
//				std::cout.flush();
	//			cmesh->Reference()->ResetReference();
		//		cmesh->LoadReferences();
			//	adapt.DeleteElements(cmesh);
				delete cmesh;
				cmesh = 0;
				if(NRefs>1) {
					CreateMesh(gmesh,dim,1);
					cmesh->CleanUpUnconnectedNodes();
				}
			}
			if(gmesh)
				delete gmesh;
		}
	}
	
	fileerrors << std::endl << std::endl;
	fileerrors.close();
	out.close();
	return 0;
}
 */
 
//int main_NoAutoHP() {
/*
int main_NoAutoHP_old() {
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();
	// To compute processing times
	time_t sttime;
	time_t endtime;
	int time_elapsed;
	const int lenstrtime = 256;
	char time_formated[lenstrtime];
	
	ofstream fileerrors("ErrorsHP2D_ArcTan.txt");   // To store all errors calculated by TPZAnalysis (PosProcess)
	
	// To compute the errors
	TPZManVector<REAL> ervec(100,0.0);
	// Printing computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	int nref, NRefs = 6;
	int nthread, NThreads = 3;
	int dim = 2;
	
	
	for(int ntyperefs=2;ntyperefs>0;ntyperefs--) {
		fileerrors << "Type of refinement: " << ntyperefs << " Level. " << endl;
		for(int itypeel=0;itypeel<2;itypeel++) {
			MElementType typeel;
			if(!itypeel) typeel = EOned;
			else typeel = EQuadrilateral;
			fileerrors << "Type of element: " << typeel << " (0-quadrilateral, 1-triangle." << endl;
			// Generating geometric mesh 2D
			cout << "\nConstructing Poisson problem. Refinement: " << nref+1 << " Threads: " << nthread << " TypeRef: " << ntyperefs << " TypeElement: " << typeel << endl;
			
			TPZGeoMesh *gmesh;
			// geometric mesh (initial)
			std::string nombre = PZSOURCEDIR;
//			if(!typeel)
//				nombre += "/Projects/Poisson_ArcTan/RegionQuadrada.dump";
//			else 
//				nombre += "/Projects/Poisson_ArcTan/RegionQuadradaT.dump";
//			gmesh = CreateGeoMesh(nombre);
			gmesh = CreateGeoMesh(typeel);
			REAL radius = 0.2;
			
			for(nref=3;nref<NRefs;nref++) {
				if(nref > 5) nthread = 2*NThreads;
				else nthread = NThreads;
				
				// Initializing the generation mesh process
				time (& sttime);
				
				// h_refinement
				// Refining near the points belong a circunference with radio r - maxime distance radius
				RefiningNearCircunference(dim,gmesh,radius,ntyperefs);
				if(ntyperefs==2) {
					nref++;
					radius *= 0.35;
				}
				else
					radius *= 0.6;
				
				//		if(nref == NRefs-1) {
				//			sprintf(saida,"gmesh_2DArcTan_H%dTR%dE%d.vtk",nref,ntyperefs,typeel);
				//			PrintGeoMeshVTKWithDimensionAsData(gmesh,saida);
				//		}
				
				// Creating computational mesh (approximation space and materials)
				int p = 8, pinit;
				TPZCompEl::SetgOrder(1);
				TPZCompMesh *cmesh = CreateMesh(gmesh,dim,1);
				dim = cmesh->Dimension();
				
				// Selecting orthogonal polynomial family to construct shape functions
				if(anothertests)
					TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;  // Setting Chebyshev polynomials as orthogonal sequence generating shape functions
				
				// Primeiro sera calculado o mayor nivel de refinamento. Remenber, the first level is zero level.
				// A cada nivel disminue em uma unidade o p, mas não será menor de 1.
				int level = 0, highlevel = 0;
				int nelem = 0;
				while(nelem < cmesh->NElements()) {
					TPZCompEl *cel = cmesh->ElementVec()[nelem++];
					if(cel) {
						level = cel->Reference()->Level();
					}
					if(level > highlevel)
						highlevel = level;
				}
				// Identifying maxime interpolation order
				if(highlevel>p-1) pinit = p;
				else pinit = highlevel+1;
				// Put order 1 for more refined element and (highlevel - level)+1 for others, but order not is greater than initial p
				nelem = 0;
				while(highlevel && nelem < cmesh->NElements()) {
					TPZCompEl *cel = cmesh->ElementVec()[nelem++];
					if(!cel) continue;
					level = cel->Reference()->Level();
					p = (highlevel-level);
					if(!p) p = 1;     // Fazendo os dois maiores niveis de refinamento devem ter ordem 1
					if(p > pinit) p = pinit;
					((TPZInterpolatedElement*)cel)->PRefine(p);
				}
				cmesh->ExpandSolution();
				cmesh->CleanUpUnconnectedNodes();
				
				// closed generation mesh process
				time (& endtime);
				time_elapsed = endtime - sttime;
				time_elapsed = endtime - sttime;
				formatTimeInSec(time_formated,lenstrtime,time_elapsed);
				out << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n";
				
				// SOLVING PROCESS
				// Initial steps
				TPZAnalysis an(cmesh);
				
				TPZParSkylineStructMatrix strskyl(cmesh,nthread);
				an.SetStructuralMatrix(strskyl);
				out << "Solving HP-Adaptive Methods...\n";
				
				TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
				direct->SetDirect(ECholesky);
				an.SetSolver(*direct);
				delete direct;
				direct = 0;
				
				// Initializing the solving process
				time (& sttime);
				// Solving
				an.Run();
				
				// Calculando o tempo que demorou para calcular em cada cenario 
				time (& endtime);
				time_elapsed = endtime - sttime;
				formatTimeInSec(time_formated,lenstrtime,time_elapsed);
				
				out << "\tRefinement: " << nref+1 << " TypeRef: " << ntyperefs << " TypeElement: " << typeel << " Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n\n";
				
				// Post processing
				std::string filename = "PoissonSolution";
				char pp[256];
				sprintf(pp,"TR%1dE%1dT%02dH%02dP%02d",ntyperefs,typeel,nthread,(nref+1),pinit);
				filename += pp;
				filename += ".vtk";
				
				TPZStack<std::string> scalarnames, vecnames;
				scalarnames.Push("Solution");
				scalarnames.Push("POrder");
				scalarnames.Push("KDuDx");
				scalarnames.Push("KDuDy");
				scalarnames.Push("KDuDz");
				scalarnames.Push("NormKDu");
				scalarnames.Push("Pressure");
				
				vecnames.Push("Derivative");
				vecnames.Push("Flux");
				vecnames.Push("MinusKGradU");
				an.DefineGraphMesh(dim,scalarnames,vecnames,filename);
				
				an.PostProcess(0,dim);
				
				// Computing error
                an.SetExact(ExactSolCircle);
				
				fileerrors << "Refinement: " << nref+1 << "  Threads: " << nthread << "  NEquations: " << cmesh->NEquations();
				an.PostProcessError(ervec,out);
				for(int rr=0;rr<ervec.NElements();rr++)
					fileerrors << "  Error_" << rr+1 << ": " << ervec[rr]; 
				fileerrors << "  TimeElapsed: " << time_elapsed << " <-> " << time_formated << std::endl;
				
				delete cmesh;
			}
			delete gmesh;
		}
	}
	
	fileerrors << std::endl << std::endl;
	fileerrors.close();
	out.close();
	return 0;
}
*/
TPZGeoMesh *CreateGeoMeshWithClassesPre(MElementType typeel) {
	TPZManVector<REAL> point(3,0.), pointlast(3,0.);
	TPZGeoMesh* gmesh;
	switch (typeel) {
		case EOned:
		{
			pointlast[0] = 1.;
			gmesh = new TPZGeoMesh;
			int Qnodes = 2;
			
			gmesh->SetMaxNodeId(Qnodes-1);
			gmesh->NodeVec().Resize(Qnodes);
			TPZVec<TPZGeoNode> Node(Qnodes);
			
			TPZVec <int64_t> TopolLine(2);
			TPZVec <int64_t> TopolPoint(1);
			
			int64_t id = 0;
			for (int j=0; j<2;j++) {
				Node[id].SetNodeId(id);
				if(!j) Node[id].SetCoord(point);//coord x
				else Node[id].SetCoord(pointlast);
				gmesh->NodeVec()[id] = Node[id];
				id++;
			}
			
			//indice dos elementos
			id = 0;
			TopolLine[0] = 0;
			TopolLine[1] = 1;
			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,materialId,*gmesh);
			id++;
			
			TopolPoint[0] = 0;
			new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,id_bc0,*gmesh);
			id++;
			TopolPoint[0] = 1;
			new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,id_bc0,*gmesh);

			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case EQuadrilateral:
		{
			gmesh = new TPZGeoMesh;
			TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
			x1[2] = 0.;
			TPZManVector<int> nx(2,1);   // subdivisions in X and Y. 
			TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
			gen.SetElementType(EQuadrilateral);       // typeel = 0 means rectangular elements, typeel = 1 means triangular elements
			gen.Read(gmesh,materialId);  // generating grid in gmesh
			gmesh->BuildConnectivity();
			TPZGeoElBC gbc10(gmesh->ElementVec()[0],4,id_bc0);
			TPZGeoElBC gbc11(gmesh->ElementVec()[0],5,id_bc0);
			TPZGeoElBC gbc12(gmesh->ElementVec()[0],6,id_bc0);
			TPZGeoElBC gbc13(gmesh->ElementVec()[0],7,id_bc0);

			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case ETriangle:
		{
			gmesh = new TPZGeoMesh;
			TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
			x1[2] = 0.;
			TPZManVector<int> nx(2,1);   // subdivisions in X and Y. 
			TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
			gen.SetElementType(ETriangle);       // typeel = 0 means rectangular elements, typeel = 1 means triangular elements
			gen.Read(gmesh,materialId);             // generating grid in gmesh
			gmesh->BuildConnectivity();
			TPZGeoElBC gbc10(gmesh->ElementVec()[0],3,id_bc0);
			TPZGeoElBC gbc11(gmesh->ElementVec()[0],4,id_bc0);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],4,id_bc0);
			TPZGeoElBC gbc13(gmesh->ElementVec()[1],5,id_bc0);

			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
			break;
		case ETetraedro:
		case EPrisma:
		case EPiramide:
		case ECube:
			gmesh = ConstructingPositiveCube(1.,typeel);
			break;
		default:
            gmesh = 0;
			break;
	}

	return gmesh;
}

int main_Failed() {

#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();
 //   gRefDBase.InitializeRefPatterns();

	// To compute processing times
	time_t sttime;
	time_t endtime;
	int time_elapsed;
	const int lenstrtime = 512;
	char time_formated[lenstrtime];
	
	// Output files
    std::ofstream convergence("convergence.txt");
	std::ofstream fileerrors("ErrorsHP_ArcTan.txt");   // To store all errors calculated by TPZAnalysis (PosProcess)
	
	// To compute the errors
	TPZManVector<REAL> ervec(100,0.0);
	// Printing computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	int nref, NRefs = 15;
    int dim;
	int nthread = 1;
	
    //Working on regular meshes
    for(int regular=1; regular>0; regular--) {
		fileerrors << "Type of mesh: " << regular << " Level. " << endl;
		MElementType typeel;
//		for(int itypeel=(int)EOned;itypeel<(int)EPolygonal;itypeel++)
//		for(int itypeel=(int)EOned;itypeel<(int)ETriangle;itypeel++)
//		for(int itypeel=(int)ETriangle;itypeel<(int)ETetraedro;itypeel++)
		for(int itypeel=(int)ECube;itypeel<(int)EPolygonal;itypeel++)
		{
			typeel = (MElementType)itypeel;
			fileerrors << "Type of element: " << typeel << endl;
			TPZGeoMesh *gmesh;
			if(!regular) {
				std::string nombre;
				GetFilenameFromGID(typeel,nombre);
				// Generating geometric mesh 
				gmesh = CreateGeoMesh(nombre);
			}
			else {
//				gmesh = CreateGeoMesh(typeel);                  // Run it First - comments the next line
				gmesh = CreateGeoMeshWithClassesPre(typeel);	// Run it Second - Comments the previous line - The mesh is the same but we are using a Pre object
			}
			dim = DefineDimensionOverElementType(typeel);
			
			// Some refinements as initial step
			UniformRefinement(1,gmesh,dim);

			// Creating computational mesh (approximation space and materials)
			int p = 1, pinit;
			pinit = p;
			TPZCompEl::SetgOrder(p);
			TPZCompMesh *cmesh = CreateMesh(gmesh,dim,1);
			gmesh->SetName("Malha Geometrica original");
			cmesh->SetName("Malha Computacional Original");
			
			// Printing geometric mesh to validate
			if(gDebug) {
				sprintf(saida,"gmesh_%dD_H%dTR%dE%d.vtk",dim,nref,regular,typeel);
				PrintGeoMeshAsCompMeshInVTKWithDimensionAsData(gmesh,saida);
			}

			// Selecting orthogonal polynomial family to construct shape functions
			if(anothertests)
				TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;  // Setting Chebyshev polynomials as orthogonal sequence generating shape functions
			
			/** Variable names for post processing */
			TPZStack<std::string> scalnames, vecnames;
			scalnames.Push("POrder");
			scalnames.Push("Solution");
			
			// Solving adaptive process
			for(nref=1;nref<NRefs;nref++) {
				cout << "\nConstructing Poisson problem " << dim << "D. Refinement: " << nref << " Threads: " << nthread << " Regular: " << regular << " TypeElement: " << typeel << endl;
				
				// Initializing the generation mesh process
				time (& sttime);
				
				// Introduzing exact solution depending on the case
				TPZAnalysis an(cmesh);
				an.SetExact(ExactSolCircle);
				{
					std::stringstream sout;
					int angle = (int) (alfa*180./M_PI + 0.5);
					sout << "Poisson" << dim << "D_MESH" << regular << "E" << typeel << "Thr" << nthread << "H" << nref << "P" << pinit << "_Ang" << angle << ".vtk";
					an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
				}
				std::string MeshFileName;
				{
					std::stringstream sout;
					int angle = (int) (alfa*180./M_PI + 0.5);
					sout << "meshAngle" << dim << "D_MESH" << regular << "E" << typeel << "Thr" << nthread << "H" << nref << "P" << pinit << "_Ang" << angle << ".vtk";
					MeshFileName = sout.str();
				}
				
				cmesh->SetName("Malha computacional adaptada");
				// Printing geometric and computational mesh
				if(gDebug) {
					cmesh->Reference()->Print(std::cout);
					cmesh->Print(std::cout);
				}
				
				// Solve using symmetric matrix then using Cholesky (direct method)
				TPZSkylineStructMatrix strskyl(cmesh);
				an.SetStructuralMatrix(strskyl);
				
				TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
				direct->SetDirect(ECholesky);
				an.SetSolver(*direct);
				delete direct;
				direct = 0;
				
				an.Run();
				
				// Post processing
				an.PostProcess(0,dim);
				if(gDebug) {
					std::ofstream out(MeshFileName.c_str());
					cmesh->LoadReferences();
					TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(),out,false);
				}
				
				// closed generation mesh process
				time(&endtime);
				time_elapsed = endtime - sttime;
				formatTimeInSec(time_formated,lenstrtime,time_elapsed);
				out << "\tRefinement: " << nref+1 << " Regular Mesh: " << regular << " TypeElement: " << typeel << " Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n\n";
				
				// Initializing the auto adaptive process
				REAL valerror =0.;
				REAL valtruerror=0.;
				TPZVec<REAL> ervec,truervec,effect;
				
				TPZAdaptMesh adapt;
				adapt.SetCompMesh(cmesh);
				std::cout << "\n\n\nEntering Auto Adaptive Methods... step " << nref << "\n\n";
				
				TPZCompMesh *adptmesh;
				if(NRefs>1) {
					time (& sttime);
					adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror,ervec,ExactSolCircle,truervec,effect,fileerrors,0,typeel);
					
					time_t endtime;
					time (& endtime);
					
					int time_elapsed = endtime - sttime;
					std::cout << "\n\nExiting Auto Adaptive Methods....step " << nref << "time elapsed " << time_elapsed << "\n\n\n\n";
					
					int prt;
					std::cout << "neq = " << cmesh->NEquations() << " error estimate = " << valerror << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;

					convergence  << cmesh->NEquations() << "\t"
					<< valerror << "\t"
					<< valtruerror << "\t"
					<< ( valtruerror / valerror ) <<  "\t"
					<< sttime <<std::endl;
					for (prt=0;prt<ervec.NElements();prt++) {
						std::cout <<"error " << ervec[prt] << "  truerror = " << truervec[prt] << "  Effect " << effect[prt] << std::endl;
					}
				}
				
				std::cout.flush();
				cmesh->Reference()->ResetReference();
				cmesh->LoadReferences();
				adapt.DeleteElements(cmesh);
				delete cmesh;
				cmesh = 0;
				if(NRefs>1) {
					cmesh = adptmesh;
					cmesh->CleanUpUnconnectedNodes();
				}
			}
			if(gmesh)
				delete gmesh;
		}
	}
	
	fileerrors << std::endl << std::endl;
	fileerrors.close();
	out.close();
	return 0;	
}

/**
 * @file Poisson 3D in hexahedra with shock problem
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
#include "TPZFrontStructMatrix.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzsbstrmatrix.h"
#include "pzfstrmatrix.h"

#include "pzmaterial.h"
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
int printingsol = 0;

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
STATE ValueK = 1000000;
REAL GlobalMaxError = 0.0;

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

void ExactSolutionSphere(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
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
int MaxPOrder = 8;

TPZManVector<REAL,3> CCircle(3,0.5);
REAL RCircle = 0.25;

int ninitialrefs = 2;


REAL ProcessingError(TPZAnalysis &analysis,TPZVec<REAL> &ervec,TPZVec<REAL> &ervecbyel);
void LoadSolutionFirstOrder(TPZCompMesh *cmesh, void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv));
void ApplyingStrategyHPAdaptiveBasedOnErrors(TPZAnalysis &analysis,REAL GlobalL2Error,TPZVec<REAL> &ervecbyel);
void ApplyingUpStrategyHPAdaptiveBasedOnGradient(TPZAnalysis &analysis,REAL &GlobalNormGradient);
void ApplyingDownStrategyHPAdaptiveBasedOnGradient(TPZAnalysis &analysis,REAL &GlobalNormGradient);
void ApplyingStrategyHPAdaptiveBasedOnExactSolution(TPZAnalysis &analysis,TPZVec<REAL> &ervecbyel,REAL MaxError,int ref);

REAL GradientNorm(TPZInterpolatedElement *el);
REAL Laplacian(TPZInterpolatedElement *el);
REAL GradientNormOnCorners(TPZInterpolatedElement *el);
REAL LaplacianOnCorners(TPZInterpolatedElement *el);

// Least Squares Method
bool AdjustingWithElipse(int dim,TPZVec<REAL> &points);
void FillingPoints2D(TPZVec<REAL> &Points);
void FillingPoints3D(TPZVec<REAL> &Points);

// MAIN FUNCTION TO NUMERICAL SOLVE WITH AUTO ADAPTIVE HP REFINEMENTS
/** Laplace equation on square 1D 2D 3D - Volker John article 2000 */
int main() {
    
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	int dim = 2;
	TPZVec<REAL> Points;
	if(dim==2)
		FillingPoints2D(Points);
	else
		FillingPoints3D(Points);
	AdjustingWithElipse(dim,Points);

	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//	gRefDBase.InitializeAllUniformRefPatterns();
    //    gRefDBase.InitializeRefPatterns();
    
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
	std::ofstream fileerrors("ErrorsHP_Poisson.txt",ios::app);   // To store all errors calculated by TPZAnalysis (PosProcess)
	// Initial message to print computed errors
	time(&sttime);
	formatTimeInSec(time_formated,256,sttime);
	fileerrors << "Approximation Error in " << time_formated << std::endl;
	
	int nref = 1, NRefs = 4;
	int nthread = 2, NThreads = 4;
    int dim;
	
    //Working on regular meshes
    for(int regular=1; regular>0; regular--) {
		// Initializing the auto adaptive process
		TPZVec<REAL> ervec,ErrorVec;
		TPZVec<long> NEquations;

		fileerrors << "Type of mesh: " << regular << " Level. " << endl;
		MElementType typeel;
        //		for(int itypeel=(int)ECube;itypeel<(int)EPolygonal;itypeel++)
		for(int itypeel=(int)ETriangle;itypeel<(int)EQuadrilateral;itypeel++)
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
                MaxPOrder = 5;
                NRefs = 4;
            }
            else if(dim==2) {
                MaxPOrder = 7;
                NRefs = 6;
            }
            else {
				MaxPOrder = 20;
                NRefs = 12;
            }
			ErrorVec.Resize(NRefs,0.0);
			NEquations.Resize(NRefs,0);
			// Printing geometric mesh to validate
			if(gDebug) {
				sprintf(saida,"gmesh_%02dD_H%dTR%dE%d.vtk",dim,nref,regular,typeel);
				PrintGeoMeshInVTKWithDimensionAsData(gmesh,saida);
			}

			/** Variable names for post processing */
			TPZStack<std::string> scalnames, vecnames;
			scalnames.Push("POrder");
			scalnames.Push("Solution");

			if(printingsol) {
				TPZGeoMesh *gmeshfirst = CreateGeoMesh(typeel);
				UniformRefinement(8,gmesh,dim);
				TPZCompEl::SetgOrder(1);
				TPZCompMesh *cmeshfirst = CreateMesh(gmesh,dim,1);
				TPZAnalysis ann(cmeshfirst,false);
				LoadSolutionFirstOrder(cmeshfirst,ExactSolutionSphere);
				{
					std::stringstream sut;
					sut << "Poisson" << dim << "D_MESHINIT_E" << typeel << "H" << std::setprecision(2) << nref << ".vtk";
					ann.DefineGraphMesh(dim,scalnames,vecnames,sut.str());
				}
				ann.PostProcess(3,dim);
				delete cmeshfirst;
				delete gmeshfirst;
			}
            UniformRefinement(ninitialrefs,gmesh,dim);
            
			// Creating computational mesh (approximation space and materials)
			int p = 1, pinit;
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

			for(nref=0;nref<NRefs;nref++) {
				out << "\nConstructing Poisson problem " << dim << "D. Refinement: " << nref << " Threads: " << nthread << " Regular: " << regular << " TypeElement: " << typeel << endl;
                std::cout << "\nConstructing Poisson problem. Type element: " << typeel << std::endl;
				if(usethreads) {
					if(nref > 5) nthread = 2*NThreads;
					else nthread = NThreads;
				}
				
				// Initializing the generation mesh process
				time(& sttime);
				
				// Solving adaptive process
				TPZAnalysis an(cmesh,true);
				an.SetExact(ExactSolutionSphere);
				// Introduzing exact solution depending on the case
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
//				TPZSBandStructMatrix strmat(cmesh);
				TPZSkylineStructMatrix strmat(cmesh);
//				TPZBandStructMatrix strmat(cmesh);
			//	TPZFrontStructMatrix<TPZFrontNonSym<REAL> > strmat(cmesh);
			//	strmat.SetQuiet(1);
				an.SetStructuralMatrix(strmat);

				TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
				direct->SetDirect(ECholesky);
//				direct->SetDirect(ELU);
				an.SetSolver(*direct);
				delete direct;
				direct = 0;
				
				an.Run();
				
				// Post processing
				an.PostProcess(0,dim);
				if(gDebug) {
					std::ofstream out(MeshFileName.c_str());
					cmesh->LoadReferences();
					TPZVTKGeoMesh::PrintCMeshVTK(cmesh->Reference(),out,false);
				}
                
				// generation mesh process finished
				time(&endtime);
				time_elapsed = endtime - sttime;
				formatTimeInSec(time_formated,256,time_elapsed);
				out << "\tRefinement: " << nref+1 << " Regular Mesh: " << regular << " TypeElement: " << typeel << " Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n\n";
				
				out << "\n\nEntering Adaptive Methods... step " << nref << "\n";
                std::cout << "\n\nEntering Adaptive Methods... step " << nref << "\n";
				fileerrors << "\n\nEntering Adaptive Methods... step " << nref << "\n";
				// Printing degree of freedom (number of equations)
				fileerrors << "Refinement: " << nref << "  Dimension: " << dim << "  NEquations: " << cmesh->NEquations();
				an.PostProcessError(ervec,out);
				ErrorVec[nref] = ervec[1];
				NEquations[nref] = cmesh->NEquations();
				// Printing obtained errors
				for(int rr=0;rr<ervec.NElements();rr++)
					fileerrors << "  Error_" << rr+1 << ": " << ervec[rr]; 
				fileerrors << "  TimeElapsed: " << time_elapsed << " <-> " << time_formated << std::endl;
				REAL Tol;
				ZeroTolerance(Tol);
				if(NRefs > 1 && nref < (NRefs-1)) {
					TPZManVector<REAL> ervecbyel;
                    REAL MaxError = 0.;
                    MaxError = ProcessingError(an,ervec,ervecbyel);
					if(MaxError > ervec[1])
						std::cout << "Local error is bigger than Global error, Ref " << nref << "." << std::endl;
					if(ervec[1] < 100*Tol) {
						std::cout << "Tolerance reached but no maxime refinements, Ref " << nref << "." << std::endl;
						fileerrors.flush();
						out.flush();
						break;
					}
//					ApplyingStrategyHPAdaptiveBasedOnErrors(an,MaxError,ervecbyel);
					ApplyingStrategyHPAdaptiveBasedOnExactSolution(an,ervecbyel,MaxError,nref);
                }
				fileerrors.flush();
				out.flush();
				// Cleaning computational mesh and creating a new computational mesh
//				if(cmesh) {
//					delete cmesh;
//					cmesh = CreateMesh(gmesh,dim,1);
//				}
			}
			if(cmesh)
				delete cmesh;
			cmesh = NULL;
			if(gmesh)
				delete gmesh;
			gmesh = NULL;
			// Writing a relation between number of degree of freedom and L2 error.
			fileerrors << "NEquations = {";
			for(nref=0;nref<NRefs-1;nref++) {
				fileerrors << NEquations[nref] << ", ";
			}
			fileerrors << NEquations[nref] << "};" << std::endl << "L2Error = {";
			for(nref=0;nref<NRefs-1;nref++) {
				fileerrors << ErrorVec[nref] << ", ";
			}
			fileerrors << ErrorVec[nref] << "};" << std::endl << "LogNEquations = Table[Log[NEquations[[i]]],{i,1,Length[NEquations]}];" << std::endl;
			fileerrors << "LogL2Errors = Table[Log[L2Error[[i]]],{i,1,Length[L2Error]}];" << std::endl;
			fileerrors << "ListPlot[Table[{LogNEquations[[i]],LogL2Errors[[i]]},{i,1,Length[LogNEquations]}],Joined->True]" << std::endl;
		}
	}
	
	fileerrors << std::endl << "Finished running.\n" << std::endl << std::endl;
	fileerrors.close();
    std::cout << std::endl << "Finished running.\n" << std::endl << std::endl;
	out.close();
    return true;
}

void ApplyingStrategyHPAdaptiveBasedOnExactSolution(TPZAnalysis &analysis,TPZVec<REAL> &ervecbyel,REAL MaxError,int nref) {

	TPZCompMesh *cmesh = analysis.Mesh();
	if(!cmesh) return;
	long nels = cmesh->NElements();
	TPZVec<long> subels;
	int j, k, pelement, dp;
	dp = 1;
	TPZVec<long> subsubels;
	TPZInterpolatedElement *el;
	float Tol;
	ZeroTolerance(Tol);
	int level;

	for(long i=0L;i<nels;i++) {
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el) continue;
		// If error is small and laplacian value is very little then the order will be minimized
        STATE GradNorm, LaplacianValue;
//		if(nref < 5)
		GradNorm = GradientNormOnCorners(el);
		if(ervecbyel[i] < 100*Tol && nref > 3) 
			continue;
        if(GradNorm > 2.) {
            // Dividing element one level
            el->Divide(el->Index(),subels,0);
            // Dividing sub elements one level more
	        for(j=0;j<subels.NElements();j++) {
				TPZInterpolatedElement* scel = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[j]]);
				level = scel->Reference()->Level();
				LaplacianValue = Laplacian(scel);
				if(LaplacianValue > 1.) {
					pelement = scel->PreferredSideOrder(scel->NConnects() - 1);
					// Applying p+1 order for all subelements
					if(pelement+1 < MaxPOrder)
						scel->PRefine(pelement+1);
				}
				if(level < 7)
					scel->Divide(subels[j],subsubels,0);
			}
		}
		else {
			if(ervecbyel[i] > 0.2*MaxError) {
				LaplacianValue = Laplacian(el);
				level = el->Reference()->Level();
	            // Dividing element one level
				if(level < 5) {
//					pelement = el->PreferredSideOrder(el->NConnects() - 1);
			        el->Divide(el->Index(),subels,0);
					el = 0;
//					if((LaplacianValue > 2.) && (pelement+dp < MaxPOrder-1)) {
	//					for(j=0;j<subels.NElements();j++) {
		//					TPZInterpolatedElement* scel = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[j]]);
			//				scel->PRefine(pelement+dp);
				//		}
					//}
				}
				if(el && LaplacianValue > 1.) {
					pelement = el->PreferredSideOrder(el->NConnects() - 1);
					el->PRefine(pelement+dp);
				}
			}
/*			else if(GradNorm < 0.1) {
				LaplacianValue = Laplacian(el);
				if(LaplacianValue < 0.1) {
					pelement = el->PreferredSideOrder(el->NConnects() - 1);
					// Applying p+1 order for all subelements
					if(pelement > 1)
						el->PRefine(pelement-1);
				}
			}*/
		}
	}
}
/**
 * Get Global L2 Error for solution and the L2 error for each element.
 * Return the maxime L2 error by elements.
 */
REAL ProcessingError(TPZAnalysis &analysis,TPZVec<REAL> &ervec,TPZVec<REAL> &ervecbyel) {
    long neq = analysis.Mesh()->NEquations();
    TPZVec<REAL> ux(neq);
    TPZVec<REAL> sigx(neq);
    TPZManVector<REAL,10> values(10,0.);
    analysis.Mesh()->LoadSolution(analysis.Solution());

	TPZAdmChunkVector<TPZCompEl *> elvec = analysis.Mesh()->ElementVec();
    TPZManVector<REAL,10> errors(10);
    errors.Fill(0.0);
    long i, nel = elvec.NElements();
	ervecbyel.Resize(nel,0.0);
	REAL maxError = 0.0;
    for(i=0L;i<nel;i++) {
        TPZCompEl *el = (TPZCompEl *) elvec[i];
        if(el) {
            errors.Fill(0.0);
            el->EvaluateError(analysis.fExact, errors, 0);
            int nerrors = errors.NElements();
            values.Resize(nerrors, 0.);
            for(int ier = 0; ier < nerrors; ier++)
                values[ier] += errors[ier] * errors[ier];
			// L2 error for each element
			ervecbyel[i] = sqrt(errors[1]*errors[1]);
			if(ervecbyel[i] > maxError)
				maxError = ervecbyel[i];
        }
    }
    
    int nerrors = errors.NElements();
	ervec.Resize(nerrors);
	ervec.Fill(-1.0);
    
	// Returns the square of the calculated errors.
	for(i=0;i<nerrors;i++)
		ervec[i] = sqrt(values[i]);
    return maxError;
}
/**
 * Criteria: Given GlobalL2Error = GE and MaxError (ME) over all the elements
 * If ElementError(EE) > 0.75*ME => twice h-refinement and p-2
 * Else if EE > 0.5*ME  =>  h-refinement and p--
 * Else if EE > 0.25*ME  =>  h-refinement
 * in the other hand => p++
 */
void ApplyingStrategyHPAdaptiveBasedOnErrors(TPZAnalysis &analysis,REAL GlobalL2Error,TPZVec<REAL> &ervecbyel) {

	TPZCompMesh *cmesh = analysis.Mesh();
	if(!cmesh) return;
	long nels = cmesh->NElements();
	TPZVec<long> subels;
	int j, k, pelement, dp = 1;
	TPZVec<long> subsubels;
	TPZInterpolatedElement *el;
	REAL errorcel = 0.0;
	for(long i=0L;i<nels;i++) {
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el) continue;
		errorcel = ervecbyel[i];
		if(errorcel < 1000*ZeroTolerance()) continue;
		// If error is small and laplacian value is very little then the order will be minimized
		if(errorcel < 0.01*GlobalL2Error) {
			REAL LaplacianValue = Laplacian(el);
			if(LaplacianValue < 0.1) {
				pelement = el->PreferredSideOrder(el->NConnects() - 1);
				// Applying p+1 order for all subelements
				if(pelement > 1)
					el->PRefine(pelement-1);
			}
		}
		else if(errorcel > 0.3*GlobalL2Error) {
            REAL GradNorm = GradientNorm(el);
            if(GradNorm > 3) {
				int level = el->Reference()->Level();
                // Dividing element one level
                el->Divide(el->Index(),subels,0);
                // Dividing sub elements one level more
				if(level < ninitialrefs+3) {
	                for(j=0;j<subels.NElements();j++) {
						TPZInterpolatedElement* scel = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[j]]);
					    scel->Divide(subels[j],subsubels,0);
				        for(k=0;k<subsubels.NElements();k++) {
							scel = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subsubels[k]]);
							REAL LaplacianValue = Laplacian(scel);
							if(LaplacianValue > 0.25) {
								pelement = scel->PreferredSideOrder(scel->NConnects() - 1);
						        // Applying p+1 order for all subelements
					            if(pelement+dp < MaxPOrder-1)
					                scel->PRefine(pelement+dp);
							}
						}
	                }
				}
            }
            else {
                el->Divide(el->Index(),subels,0);
                // Dividing sub elements one level more
                for(j=0;j<subels.NElements();j++) {
					TPZInterpolatedElement* scel = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[j]]);
					REAL LaplacianValue = Laplacian(scel);
					if(LaplacianValue > 0.25) {
						pelement = scel->PreferredSideOrder(scel->NConnects() - 1);
	                    // Applying p+1 order for all subelements
						scel->PRefine(pelement+1);
					}
					else if(IsZero(LaplacianValue)) {
						pelement = scel->PreferredSideOrder(scel->NConnects() - 1);
	                    // Applying p-1 order for all subelements
						if(pelement > 1) scel->PRefine(pelement-1);
					}
				}
            }
		}
/*		else {
			el->Divide(el->Index(),subels,0);
            // Dividing sub elements one level more
            for(j=0;j<subels.NElements();j++) {
				TPZInterpolatedElement* scel = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[j]]);
				REAL LaplacianValue = Laplacian(scel);
				if(LaplacianValue > 0.25) {
					pelement = scel->PreferredSideOrder(scel->NConnects() - 1);
	                // Applying p+1 order for all subelements
					scel->PRefine(pelement+1);
				}
			}
		}
		*/
	}
}
/*
REAL GradientNorm(TPZInterpolatedElement *el) {
    TPZSolVec sol;
    TPZGradSolVec dsol;
	TPZVec<REAL> qsi(3,0.0);
    int nshape = el->NShapeF();
    int dim = el->Dimension();
    REAL temp = 0.0;
    TPZFMatrix<REAL> phi(nshape,1);
    TPZFMatrix<REAL> dphi(dim,nshape);
    TPZFMatrix<REAL> axes(dim,dim,0.0);
    el->Reference()->CenterPoint(el->NConnects() - 1, qsi);
    el->Shape(qsi,phi,dphi);
    el->ComputeSolution(qsi,phi,dphi,axes,sol,dsol);
    for(int idsol=0;idsol<dim;idsol++)
        temp += dsol[0](idsol,0)*dsol[0](idsol,0);
    return (sqrt(temp)/el->VolumeOfEl());
}
*/
REAL GradientNorm(TPZInterpolatedElement *el) {
	int nstates = el->Material()->NStateVariables();
    int dim = el->Dimension();
    TPZVec<STATE> sol(nstates,(STATE)0.0);
    TPZFMatrix<STATE> dsol(dim,nstates);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<REAL> point(3,0.0);
    REAL temp = 0.0;

	el->Reference()->CenterPoint(el->NConnects() - 1, qsi);
	el->Reference()->X(qsi,point);
	ExactSolutionSphere(point,sol,dsol);
    for(int idsol=0;idsol<dim;idsol++)
        temp += dsol(idsol,0)*dsol(idsol,0);
    return sqrt(temp);
}
REAL Laplacian(TPZInterpolatedElement *el) {
	int nvar = el->Material()->NStateVariables();
    int dim = el->Dimension();
	TPZVec<STATE> sol(nvar,0.);
    TPZFMatrix<STATE> dsol(nvar,dim,0.0);
	TPZVec<REAL> qsi(3,0.0);
	TPZVec<REAL> point(3,0.0);
    REAL temp = 0.0;
    el->Reference()->CenterPoint(el->NConnects() - 1, qsi);
	el->Reference()->X(qsi,point);
	RightTermCircle(point,sol,dsol);
	return fabs(sol[0]/ValueK);
}
REAL GradientNormOnCorners(TPZInterpolatedElement *el) {
	int nstates = el->Material()->NStateVariables();
    int dim = el->Dimension();
    TPZVec<STATE> sol(nstates,(STATE)0.0);
    TPZFMatrix<STATE> dsol(dim,nstates);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<REAL> point(3,0.0);
	REAL MaxGrad = 0.0;
	int ncorners = el->NCornerConnects();;
	for(int i=0;i<ncorners;i++) {
	    REAL temp = 0.0;
		el->Reference()->CenterPoint(i, qsi);
		el->Reference()->X(qsi,point);
		ExactSolutionSphere(point,sol,dsol);
		for(int idsol=0;idsol<dim;idsol++)
			temp += dsol(idsol,0)*dsol(idsol,0);
		MaxGrad = (MaxGrad > sqrt(temp)) ? MaxGrad : sqrt(temp);
	}
    return MaxGrad;
}
REAL LaplacianOnCorners(TPZInterpolatedElement *el) {
	int nvar = el->Material()->NStateVariables();
    int dim = el->Dimension();
	TPZVec<STATE> sol(nvar,0.);
    TPZFMatrix<STATE> dsol(nvar,dim,0.0);
	TPZVec<REAL> qsi(3,0.0);
	TPZVec<REAL> point(3,0.0);
	REAL MaxLaplacian = 0.0;
	int ncorners = el->NCornerConnects();;
	for(int i=0;i<ncorners;i++) {
	    el->Reference()->CenterPoint(i, qsi);
		el->Reference()->X(qsi,point);
		RightTermCircle(point,sol,dsol);
		MaxLaplacian = (MaxLaplacian > fabs(sol[0])) ? MaxLaplacian : fabs(sol[0]);
	}
	return (MaxLaplacian/ValueK);
}
void ApplyingUpStrategyHPAdaptiveBasedOnGradient(TPZAnalysis &analysis,REAL &GlobalNormGradient) {
	TPZVec<REAL> ervecbyel;
	TPZCompMesh *cmesh = analysis.Mesh();
	if(!cmesh) return;
	long nels = cmesh->NElements();
	ervecbyel.Resize(nels,0.0);
	TPZVec<long> subels;
	TPZVec<long> subsubels;
	int j, k, pelement;
	TPZInterpolatedElement *el;
	// Computing norm of gradient on center of each element
	TPZVec<REAL> qsi(3,0.0);
	REAL temp;
	int dim = analysis.Mesh()->Dimension();

    // Computing reference gradient (global)
	for(long ii=0L;ii<nels;ii++) {
		temp = 0.0;
		TPZSolVec sol;
		TPZGradSolVec dsol;
		int nshape;
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[ii]);
		if(!el || el->Dimension() != dim) continue;
		nshape = el->NShapeF();
		TPZFMatrix<REAL> phi(nshape,1);
		TPZFMatrix<REAL> dphi(dim,nshape);
		TPZFMatrix<REAL> axes(dim,dim,0.0);
		el->Reference()->CenterPoint(el->NConnects() - 1, qsi);
		el->Shape(qsi,phi,dphi);
		el->ComputeSolution(qsi,phi,dphi,axes,sol,dsol);
		for(int idsol=0;idsol<dim;idsol++)
			temp += dsol[0](idsol,0)*dsol[0](idsol,0);
		ervecbyel[ii] = sqrt(temp);
		GlobalNormGradient = (GlobalNormGradient > ervecbyel[ii]) ? GlobalNormGradient : ervecbyel[ii];
	}
    
    // Applying hp refinement based on local gradient comparative with global gradient (in norm)
	for(long i=0L;i<nels;i++) {
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el) continue;
		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		if(el->Reference()->Level() < 4 && ervecbyel[i] > 0.6*GlobalNormGradient) {
			// Dividing element one level
			el->Divide(el->Index(),subels,0);
			// Dividing sub elements one level more
			for(j=0;j<subels.NElements();j++) {
				if(cmesh->ElementVec()[subels[j]]->Reference()->Level() < 5) {
                    cmesh->ElementVec()[subels[j]]->Divide(subels[j],subsubels,0);
//                    if(cmesh->ElementVec()[subsubels[0]]->Reference()->Level() < 3)
                    for(k=0;k<subsubels.NElements();k++) {
                        if(pelement < MaxPOrder-1) pelement++;
						// Applying p-2 order for all subelements
                        ((TPZInterpolatedElement*)cmesh->ElementVec()[subsubels[k]])->PRefine(pelement);
                    }
                }
			}
		}
		else if(ervecbyel[i] > 0.2*GlobalNormGradient) {
			// Dividing element one level
			el->Divide(el->Index(),subels,0);
//			if(ervecbyel[i] > 0.4*GlobalNormGradient) {
  //              if(pelement < MaxPOrder) pelement++;
	//			for(j=0;j<subels.NElements();j++)
					// Applying p-2 order for all subelements
	//				((TPZInterpolatedElement*)cmesh->ElementVec()[subels[j]])->PRefine(pelement);
	//		}
		}
//		else if(ervecbyel[i] < 0.1*GlobalNormGradient) {
//			pelement--;
//			if(pelement > 1)
//				el->PRefine(pelement);
//		}
	}
}
void ApplyingDownStrategyHPAdaptiveBasedOnGradient(TPZAnalysis &analysis,REAL &GlobalNormGradient) {
	TPZVec<REAL> ervecbyel;
	TPZCompMesh *cmesh = analysis.Mesh();
	if(!cmesh) return;
	long nels = cmesh->NElements();
	ervecbyel.Resize(nels,0.0);
	TPZVec<long> subels;
	TPZVec<long> subsubels;
	int j, k, pelement;
	TPZInterpolatedElement *el;
	// Computing norm of gradient on center of each element
	TPZVec<REAL> qsi(3,0.0);
	REAL temp;
	int dim = analysis.Mesh()->Dimension();
    
    // Computing reference gradient (global)
	for(long ii=0L;ii<nels;ii++) {
		temp = 0.0;
		TPZSolVec sol;
		TPZGradSolVec dsol;
		int nshape;
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[ii]);
		if(!el || el->Dimension() != dim) continue;
		nshape = el->NShapeF();
		TPZFMatrix<REAL> phi(nshape,1);
		TPZFMatrix<REAL> dphi(dim,nshape);
		TPZFMatrix<REAL> axes(dim,dim,0.0);
		el->Reference()->CenterPoint(el->NConnects() - 1, qsi);
		el->Shape(qsi,phi,dphi);
		el->ComputeSolution(qsi,phi,dphi,axes,sol,dsol);
		for(int idsol=0;idsol<dim;idsol++)
			temp += dsol[0](idsol,0)*dsol[0](idsol,0);
		ervecbyel[ii] = sqrt(temp);
		GlobalNormGradient = (GlobalNormGradient > ervecbyel[ii]) ? GlobalNormGradient : ervecbyel[ii];
	}
    
    // Applying hp refinement based on local gradient comparative with global gradient (in norm)
	for(long i=0L;i<nels;i++) {
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el) continue;
		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		if(ervecbyel[i] > 0.6*GlobalNormGradient) {
			// Dividing element one level
			el->Divide(el->Index(),subels,0);
			// Dividing sub elements one level more
			for(j=0;j<subels.NElements();j++) {
				if(cmesh->ElementVec()[subels[j]]->Reference()->Level() < 5) {
                    cmesh->ElementVec()[subels[j]]->Divide(subels[j],subsubels,0);
                    //                    if(cmesh->ElementVec()[subsubels[0]]->Reference()->Level() < 3)
                    for(k=0;k<subsubels.NElements();k++) {
                        if(pelement > 1) pelement--;
						// Applying p-2 order for all subelements
                        ((TPZInterpolatedElement*)cmesh->ElementVec()[subsubels[k]])->PRefine(pelement);
                    }
                }
			}
		}
		else if(ervecbyel[i] > 0.2*GlobalNormGradient) {
			// Dividing element one level
			el->Divide(el->Index(),subels,0);
			if(ervecbyel[i] > 0.5*GlobalNormGradient) {
                if(pelement < MaxPOrder) pelement++;
				for(j=0;j<subels.NElements();j++)
					// Applying p-2 order for all subelements
					((TPZInterpolatedElement*)cmesh->ElementVec()[subels[j]])->PRefine(pelement);
			}
		}
		else if(ervecbyel[i] < 0.1*GlobalNormGradient) {
			pelement++;
			if(pelement < MaxPOrder)
				el->PRefine(pelement);
		}
	}
}

void LoadSolutionFirstOrder(TPZCompMesh *cmesh, void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)) {
	TPZFMatrix<STATE> solution(cmesh->Solution());
	solution.Zero();
	TPZVec<STATE> sol(1);
	TPZFMatrix<STATE> dsol(cmesh->Dimension(),1);
	int nels = cmesh->NElements();
	for(int i=0;i<nels;i++) {
		TPZCompEl *el=cmesh->ElementVec()[i];
		if(!el) continue;
		int nconn = el->Reference()->NCornerNodes();
		for(int j=0;j<nconn;j++) {
			TPZConnect &df = el->Connect(j);
			int seqnum = df.SequenceNumber();
			int pos = cmesh->Block().Position(seqnum);
			TPZVec<REAL> coord(3,0.), coordinates(3,0.);
			el->Reference()->CenterPoint(j,coord);
			el->Reference()->X(coord,coordinates);
			f(coordinates,sol,dsol);
			cmesh->Solution()(pos,0) = sol[0];
		}
	}
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


void ExactSolutionSphere(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
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
			
			TPZVec <long> TopolLine(2);
			TPZVec <long> TopolPoint(1);
			
			long id = 0;
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
			TPZVec<long> nodind(nnode);
			for(nod=0; nod<nnode; nod++) nodind[nod]=nod;
			long index;
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
				int nodind = gmesh->NodeVec().AllocateNewElement();
				TPZVec<REAL> coord(3,0.0);
				coord[0] = co[nod][0];
				coord[1] = co[nod][1];
				gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
			}
			
			int el;
			for(el=0; el<nelem; el++) {
				TPZVec<long> nodind(3);
				for(nod=0; nod<3; nod++) nodind[nod]=indices[el][nod];
				long index;
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
	TPZVec<TPZVec<long> > indices(nelem);
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
		long index;
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
    
	TPZVec<TPZVec<long> > indices(nelem);
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
		long index;
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

	TPZVec<TPZVec<long> > indices(nelem);
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
		long index;
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

	TPZVec<TPZVec<long> > indices(nelem);
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
		long index;
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
	TPZVec<TPZVec<long> > indices(nelem);
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
		long index;
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


bool LeastSquaresToGetElipse(int dim,TPZVec<REAL> &points,TPZFMatrix<REAL> &Coefficients) {
	if(dim < 2 || dim > 3) return false;
	long npoints = points.NElements()/dim;
	int nincog, i;
	if(dim == 2) nincog = 4;
	else if(dim == 3) nincog = 6;

	if(npoints<nincog) return false;

	// Dimensioning vector of coefficients
	Coefficients.Redim(nincog,1);
	Coefficients.Zero();

	// Will be solved y^2 = p*x^2 + q*x + r*y + s
	// Constructing matrix H and Transpose of H to compute by least squares method
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	TPZFMatrix<REAL> A;

	// Redimensioning
	A.Redim(nincog,nincog);
    DeltaH.Redim(npoints,nincog);
    DeltaHTranspose.Redim(nincog,npoints);
    DifSol.Redim(npoints,1);
//	DifSol.Zero();

	// Filling y^2 into Coefficients
	for(i=0;i<npoints;i++)
		DifSol.PutVal(i,0,points[2*i+1]*points[2*i+1]);

	// Filling elements for H matrix
	for(int i=0;i<npoints;i++) {
		DeltaH.PutVal(i,0,points[2*i]*points[2*i]);
		DeltaH.PutVal(i,1,points[2*i]);
		DeltaH.PutVal(i,2,points[2*i+1]);
		DeltaH.PutVal(i,3,1.);
	}
	DeltaH.Print(std::cout);

    // Solving by least squares using product of matrix: DeltaH_t * DifSol = DeltaH_t * DeltaH * Coeffs(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    A = DeltaHTranspose*DeltaH;
	A.Print(std::cout);
	Coefficients = DeltaHTranspose*DifSol;
//	A.SolveDirect(DifSol,ELU);
	Coefficients.Print(std::cout);
	A.SolveDirect(Coefficients,ELU);
//	Coefficients = DifSol;
}

bool StandardFormatForElipse(TPZFMatrix<REAL> &Coeffs,TPZVec<REAL> &Center,TPZVec<REAL> &Ratios) {
	int dim = Center.NElements();
	int ncoeffs = Coeffs.Rows();
	REAL temp;
	if(ncoeffs != 2*dim || Coeffs.Cols()!=1)
		return false;
	if(dim ==2) {
		Center[0] = -(Coeffs(1,0)/(2.*Coeffs(0,0)));
		Center[1] = 0.5*Coeffs(2,0);
		// Computing Ratios[1] in temp
		temp = Coeffs(3,0)-(Center[0]*Center[0]*Coeffs(0,0))+(Center[1]*Center[1]);
		// Computing Ratios[0] in Ratios[1]
		Ratios[1] = -temp/Coeffs(0,0);
		if(temp < 0. || Ratios[1] < 0.)
			return false;
		Ratios[0] = sqrt(Ratios[1]);
		Ratios[1] = sqrt(temp);
	}
	else {
	}
	return true;
}
bool AdjustingWithElipse(int dim,TPZVec<REAL> &Points) {

	TPZFMatrix<REAL> Coeffs;
	TPZVec<REAL> Center(dim,0.);
	TPZVec<REAL> Ratios(dim,0.);

	// Applying least squares for these five points
	Points.Print(std::cout);
	LeastSquaresToGetElipse(dim,Points,Coeffs);
	Coeffs.Print(std::cout);
	std::cout << "\n\nSolution:";

	// Making zero depending on Tolerance
	float Tol;
	ZeroTolerance(Tol);
	for(int i=0;i<Coeffs.Rows();i++)
		if(fabs(Coeffs(i)) < 100*Tol)
			Coeffs.PutVal(i,0,0.);
	if(dim == 2) {
		std::cout << std::endl << "y*y = " << Coeffs(0,0) << "x*x + " << Coeffs(1,0) << "x + " << Coeffs(2,0) << "y + " << Coeffs(3,0) << "\n";

		if(!StandardFormatForElipse(Coeffs,Center,Ratios))
			return false;
		std::cout << "\nElipse: (x - " << Center[0] << ")/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")/" << Ratios[1]*Ratios[1] << " = 1.\n" << std::endl;
	}
	else {
		std::cout << std::endl << "z*z = " << Coeffs(0,0) << "x*x + " << Coeffs(1,0) << "x + " << Coeffs(2,0) << "y*y + " << Coeffs(3,0) << "y +";
		std::cout << Coeffs(4,0) << "z + " << Coeffs(5,0) << std::endl;
		if(!StandardFormatForElipse(Coeffs,Center,Ratios))
			return false;
		std::cout << "\nElipse: (x - " << Center[0] << ")/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")/" << Ratios[1]*Ratios[1];
		std::cout << " + (z - " << Center[2] << ")/" << Ratios[2]*Ratios[2] << " = 1.\n" << std::endl;
	}
	return true;
}
void FillingPoints2D(TPZVec<REAL> &Points) {
	Points.Resize(16);

	// Filling points coordinates
	Points[0] = 1.;
	Points[1] = 3;
	Points[2] = 1.;
	Points[3] = -3;
	Points[4] = -1.;
	Points[5] = 0.;
	Points[6] = 3.;
	Points[7] = 0.;
	Points[8] = 0.;
	Points[9] = (3./2.)*sqrt(3.);
	Points[10] = 0.;
	Points[11] = -(3./2.)*sqrt(3.);
	Points[12] = 2.97203;
	Points[13] = 0.5;
	Points[14] = 2.49071;
	Points[15] = -2.;
}
void FillingPoints3D(TPZVec<REAL> &Points) {
	Points.Resize(18);

	// Filling points coordinates
	Points[0] = 1.;
	Points[1] = 3;
	Points[2] = 0.;

	Points[3] = 1.;
	Points[4] = -3;
	Points[5] = 0.;

	Points[6] = -1.;
	Points[7] = 0.;
	Points[8] = 0.;

	Points[9] = 3.;
	Points[10] = 0.;
	Points[11] = 0.;

	Points[12] = 0.;
	Points[13] = (3./2.)*sqrt(3.);
	Points[14] = 0.;

	Points[15] = 1.;
	Points[16] = 0.;
	Points[17] = 1.;
}
/******   *****   ****   ******************** /////////////////////////////

using namespace std;
using namespace pzshape;
using namespace pzgeom;
*/
/**
 * @addtogroup Tutorials
 * @{
 */

/** VARIABLES */
/** Printing level */
/*
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

// Global variable
int gLMax;
int NUniformRefs = 2;
int nstate = 2;

// Alfa -> Coefficient of the arctang argument
REAL ALFA = 50.;

char saida[512];
ofstream out("ConsolePoisson3D.txt");   // output file from console  -> Because it has many energy faults


void ExactShock(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);

void BCDirichletShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol);
void BCNeumannLateralFrontShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol);
void BCNeumannLateralRightShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol);
void BCNeumannTopShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol);
void FforcingShock(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix<REAL> &df);

void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=true, const int matidtodivided=1);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL,3> &points,REAL r,REAL &distance,int ntyperefs);
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,REAL radius,int ntyperefs);

TPZGeoMesh *ConstructingPositiveCube(REAL L,MElementType typeel);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);

void formatTimeInSec(char *strtime,int timeinsec);

/ ** Detects the bigger dimension of the computational elements into cmesh to set the Model Dimension * /
bool DefineModelDimension(TPZCompMesh *cmesh);

REAL Radius[3][10] = {{0.3,0.25,0.2,0.15,0.075,0.025,0.008,0.004,0.002,0.001},{0.3,0.2,0.075,0.008,0.002,0.0005,0.0,0.0,0.0,0.0},{0.3,0.2,0.025,0.002,0.0005,0.0,0.0,0.0,0.0,0.0}};

int main_old(int argc, char *argv[]) {

#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing a ref patterns
	gRefDBase.InitializeAllUniformRefPatterns();
    //	gRefDBase.InitializeRefPatterns();
    
	// To compute processing times
	time_t sttime;
	time_t endtime;
	int time_elapsed;
	char tempo[512];
    
	// To print computed errors
    char errorname[512];
    ofstream fileerrors;
	TPZVec<REAL> ervec(100,0.0);
	fileerrors << "Approximation Error - Shock problem: " << std::endl;
	
	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	REAL InitialL = 1.0;
	int nref, NRefs = 3;
    int NMaxTypeRefs = 4;
	int nthread, NThreads = 2;
    for(int dim=3;dim<4;dim++) {
        sprintf(errorname,"ErrorsHP_%dD.txt",dim);
        fileerrors.open(errorname);
		MElementType typeel;
        REAL radius = 0.0;
		for(int itypeel=(int)ECube;itypeel<(int)EPolygonal;itypeel++) {
			typeel = (MElementType)itypeel;
			for(int ntyperefs=1;ntyperefs<NMaxTypeRefs;ntyperefs++) {
                if(ntyperefs==1)
                    NRefs = 9;
                else if(ntyperefs==2)
                    NRefs = 4;
                else
                    NRefs = 3;
				fileerrors << "Type of refinement: " << ntyperefs << " Level. " << endl;
                fileerrors << "Type of element: " << typeel << endl;
                // Constructing geometric mesh as hexahedra
                cout << "\nConstructing Shock problem in cube [0,1]^" << dim << ". TypeRef: " << ntyperefs << " TypeElement: " << typeel << std::endl;
                TPZGeoMesh *gmesh = ConstructingPositiveCube(InitialL,typeel);
                UniformRefinement(1,gmesh,3);
                for(nref=0;nref<NRefs;nref++) {
                    radius = Radius[ntyperefs-1][nref];
                    if(nref > 4) nthread = 2*NThreads;
                    else nthread = NThreads;
                    
                    // Initializing the generation mesh process
                    time(&sttime);
                    cout << "Refinement: " << nref+1 << " Threads: " << nthread << std::endl << std::endl;
                    // h_refinement
                    // Refining near to the origin
//                    if(!nref) RefiningNearCircunference(dim,gmesh,radius,1);
//                    else
                    RefiningNearCircunference(dim,gmesh,radius,ntyperefs);
                    
                    // Creating computational mesh
                    / ** Set polynomial order * /
                    int p = 5, pinit;
                    pinit = p;
                    TPZCompEl::SetgOrder(1);
                    TPZCompMesh *cmesh = CreateMesh(gmesh,dim,1);
                    cmesh->SetName("Computational mesh for Fichera problem");
                    dim = cmesh->Dimension();
                    
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
                        if(!p || p<0) p = 1;     // Fazendo os dois maiores niveis de refinamento devem ter ordem 1
                        if(p > pinit) p = pinit;
                        ((TPZInterpolatedElement*)cel)->PRefine(p);
                    }
                    cmesh->ExpandSolution();
                    cmesh->CleanUpUnconnectedNodes();
                    
                    // closed generation mesh process
                    time (& endtime);
                    time_elapsed = endtime - sttime;
                    formatTimeInSec(tempo, time_elapsed);
                    out << "  Time elapsed " << time_elapsed << " <-> " << tempo << "\n\n";
                    
                    //--- END construction of the meshes
                    
                    // Solving linear equations
                    // Initial steps
                    out << "Solving HP-Adaptive Methods...\n";
                    
                    TPZAnalysis an (cmesh);
                    
                    // Solve using symmetric matrix then using Cholesky (direct method)
                    TPZParSkylineStructMatrix strskyl(cmesh,nthread);
                    an.SetStructuralMatrix(strskyl);
                    
                    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
                    direct->SetDirect(ECholesky);
                    an.SetSolver(*direct);
                    delete direct;
                    direct = 0;
                    
                    // Initializing the solving process
                    time (& sttime);
                    an.Run();
                    time(&endtime);
                    time_elapsed = endtime - sttime;
                    formatTimeInSec(tempo,time_elapsed);
                    
                    out << "\tRefinement: " << nref << " TypeRef: " << ntyperefs << " TypeElement: " << typeel << " Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << tempo << "\n\n\n";
                    
                    // Post processing
                    char pp[64];
                    std::string filename = "Poisson3DSol_";
                    sprintf(pp,"TR%1dE%1dT%02dH%02dP%02d",ntyperefs,typeel,nthread,nref,pinit);
                    filename += pp;
                    filename += ".vtk";
                    
                    / ** Variable names for post processing * /
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
                    // END Determining the name of the variables
                    
                    an.DefineGraphMesh(dim,scalarnames,vecnames,filename);
                    
                    an.PostProcess(0,dim);
                    ervec.Fill(0.0);
                    
                    // Computing error
                    an.SetExact(ExactShock);
                    fileerrors << "Refinement: " << nref << "  Dimension: " << dim << "  NEquations: " << cmesh->NEquations();
                    an.PostProcessError(ervec,out);
                    for(int rr=0;rr<ervec.NElements();rr++)
                        fileerrors << "  Error_" << rr+1 << ": " << ervec[rr]; 
                    fileerrors << "  TimeElapsed: " << time_elapsed << " <-> " << tempo << std::endl;
                    
                    delete cmesh;
                }
                delete gmesh;
            }
        }
        fileerrors.close();
    }
	out.close();
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////   SHOCK PROBLEM    ///////////////////
////////////////////////////////////////////////////////////////////////////////////////
void ExactShock(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
	TPZVec<REAL> C0(3,0.5);
    REAL Radio = 0.25;
	REAL Product = 8.;
	REAL R0 = 0.;
    int dim = dsol.Cols();
    for(int i=0;i<dim;i++) {
        Product *= x[i]*(x[i] - 1.);
        R0 += (x[i]-C0[i])*(x[i]-C0[i]);
    }
	sol[0] = Product * (1. + (2./M_PI)*atan( 2*sqrt(ALFA) * ( Radio*Radio - R0 )));
    
	REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	if(IsZero(den))
		DebugStop();
	dsol(0,0) = (ALFA*(x[0]-C0[0]))/den;
	dsol(1,0) = (ALFA*(x[1]-C0[1]))/den;
	dsol(2,0) = (ALFA*(x[2]-C0[2]))/den;
}

void BCDirichletShock(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	bcsol[0] = atan(ALFA * ( R0 - sqrt(3.)) );
}
void BCNeumannLateralFrontShock(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	if(IsZero(den))
		DebugStop();
	bcsol[0] = (ALFA*(x[0]-C0[0]))/den;
}
void BCNeumannLateralRightShock(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	if(IsZero(den))
		DebugStop();
	bcsol[0] = (ALFA*(x[1]-C0[1]))/den;
}
void BCNeumannTopShock(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	if(IsZero(den))
		DebugStop();
	bcsol[0] = (ALFA*(x[2]-C0[2]))/den;
}

/ ** NOTE: Forcing function in TPZMatPoisson3d is negative * /
void FforcingShock(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	REAL temp =  (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	REAL den = R0*temp*temp;
	if(IsZero(den))
		DebugStop();
	f[0] = (2*ALFA*(1.+(3.*ALFA*ALFA)-(ALFA*ALFA*sqrt(3.)*R0)))/den;
    df.Zero();
}

TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction) {
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
	// Creating Poisson material
	TPZMaterial *mat = new TPZMatPoisson3d(1,dim);
	TPZVec<REAL> convd(3,0.);
	((TPZMatPoisson3d *)mat)->SetParameters(1.,0.,convd);
	if(hasforcingfunction) {
		mat->SetForcingFunction(new TPZDummyFunction<STATE>(FforcingShock));
	}
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	
	// Boundary conditions
	// Dirichlet on face into the coordinate planes XY YZ XZ
	TPZAutoPointer<TPZFunction<STATE> > FunctionBC = new TPZDummyFunction<STATE>(BCDirichletShock);
	TPZFMatrix<STATE> val1(dim,dim,0.),val2(dim,1,0.);
	val1.PutVal(0,0,1.); val1.PutVal(1,1,1.); val1.PutVal(2,2,1.);
	TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
	bnd->SetForcingFunction(FunctionBC);
	cmesh->InsertMaterialObject(bnd);
	// Neuman on face lateral front
	TPZAutoPointer<TPZFunction<STATE> > FunctionBCF = new TPZDummyFunction<STATE>(BCNeumannLateralFrontShock);
	TPZMaterial *bndF = mat->CreateBC(mat,-2,1,val1,val2);
	bndF->SetForcingFunction(FunctionBCF);
	cmesh->InsertMaterialObject(bndF);
	// Neuman on face lateral right
	TPZAutoPointer<TPZFunction<STATE> > FunctionBCR = new TPZDummyFunction<STATE>(BCNeumannLateralRightShock);
	TPZMaterial *bndR = mat->CreateBC(mat,-3,1,val1,val2);
	bndR->SetForcingFunction(FunctionBCR);
	cmesh->InsertMaterialObject(bndR);
	// Neuman on face top
	TPZAutoPointer<TPZFunction<STATE> > FunctionBCT = new TPZDummyFunction<STATE>(BCNeumannTopShock);
	TPZMaterial *bndT = mat->CreateBC(mat,-4,1,val1,val2);
	bndT->SetForcingFunction(FunctionBCT);
	cmesh->InsertMaterialObject(bndT);
	
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->ExpandSolution();
	return cmesh;
}
/ ** Detects the bigger dimension of the computational elements into cmesh to set the Model Dimension * /
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
	TPZVec<TPZVec<int> > indices;
    indices.Resize(nelem);
    indices[0].Resize(nnode);
    for(int i=0;i<nnode;i++) {
        indices[0][i] = i;
    }
    switch(typeel) {
        case EOned:
		case ETriangle:
		case EQuadrilateral:
			return 0;
        default:
            break;
    }
    

	TPZGeoEl *elvec[nelem];
	TPZGeoMesh *gmesh = new TPZGeoMesh();

	int nod;
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
		int index;
		elvec[el] = gmesh->CreateGeoElement(typeel,indices[el],1,index);
	}
	TPZVec<TPZGeoEl *> sub;
    switch (typeel) {
        case EPrisma:   // hexahedron -> four prisms
		{
            // Dividing hexahedron in prisms
            std::string filename = REFPATTERNDIR;
            filename += "/3D_Hexa_directional_2faces.rpt";
            
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
        case EPiramide:
		{
            // Dividing hexahedron in four pyramids
            std::string filename = REFPATTERNDIR;
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
        case ETetraedro:
		{
            // Dividing hexahedron in two tetrahedras, two prisms and one pyramid
            std::string filename = REFPATTERNDIR;
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
        case ECube:
            break;
        default:
		{
            // hexahedron -> three prisms
            // Dividing hexahedron in prisms
            std::string filename = REFPATTERNDIR;
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

	gmesh->BuildConnectivity();
	
	switch(typeel) {
		case ECube:
		{
			// face 0 (20) bottom XY - face 1 (21) lateral left XZ - face 4 (24) lateral back YZ : Dirichlet
			TPZGeoElBC gbc10(gmesh->ElementVec()[0],20,-1);
			TPZGeoElBC gbc11(gmesh->ElementVec()[0],21,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[0],24,-1);
			
			// face 2 (22) Neumann - Partial derivative (du/dx) - lateral front
			TPZGeoElBC gbc13(gmesh->ElementVec()[0],22,-2);
			// face 3 (23) Neumann - Partial derivative (du/dy) - lateral right
			TPZGeoElBC gbc14(gmesh->ElementVec()[0],23,-3);
			// face 5 (25) Neumann - Partial derivative (du/dz) - top
			TPZGeoElBC gbc15(gmesh->ElementVec()[0],25,-4);
		}
			break;
		case 1:
		{
			// First sub element - faces: 15, 16, 18 and 19
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],15,-1);
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],16,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],18,-1);
			TPZGeoElBC gbc13(gmesh->ElementVec()[1],19,-3);
			// Second sub element - faces: 15 and 19
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],15,-1);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],19,-3);
			// Third sub element - faces: 15, 16, 17 and 19
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],15,-1);
			TPZGeoElBC gbc31(gmesh->ElementVec()[3],17,-2);
			TPZGeoElBC gbc32(gmesh->ElementVec()[3],18,-4);
			TPZGeoElBC gbc33(gmesh->ElementVec()[3],19,-3);
			// Fouthrd sub element - faces: 15, 17, 18 and 19
			TPZGeoElBC gbc40(gmesh->ElementVec()[4],15,-1);
			TPZGeoElBC gbc41(gmesh->ElementVec()[4],17,-4);
			TPZGeoElBC gbc42(gmesh->ElementVec()[4],18,-1);
			TPZGeoElBC gbc43(gmesh->ElementVec()[4],19,-3);
			gmesh->ElementVec()[1]->Divide(sub);
			gmesh->ElementVec()[3]->Divide(sub);
		}
			break;
		case 2:
		{
			// First sub element - faces: 13, 14 and 17
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],13,-2);
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],14,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],17,-1);
			// Second sub element - faces: 13 and 14
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],13,-3);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],14,-1);
			// Third sub element - faces: 13, 14 and 17
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],13,-1);
			TPZGeoElBC gbc31(gmesh->ElementVec()[3],14,-1);
			TPZGeoElBC gbc32(gmesh->ElementVec()[3],17,-1);
			// Fouthrd sub element - faces: 13 and 14
			TPZGeoElBC gbc40(gmesh->ElementVec()[4],13,-4);
			TPZGeoElBC gbc41(gmesh->ElementVec()[4],14,-1);
		}
			break;
		case 3:
		{
			// First sub element - faces: 10, 11 and 13
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],10,-4);
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],11,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],13,-2);
			// Second sub element - faces: 10, 11 and 13
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],10,-4);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],11,-2);
			TPZGeoElBC gbc22(gmesh->ElementVec()[2],13,-3);
			// Third sub element - faces: 15, 16, 18 and 19
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],15,-3);
			TPZGeoElBC gbc31(gmesh->ElementVec()[3],16,-1);
			TPZGeoElBC gbc32(gmesh->ElementVec()[3],18,-4);
			TPZGeoElBC gbc33(gmesh->ElementVec()[3],19,-1);
			// Fouthrd sub element - faces: 15, 18 and 19
			TPZGeoElBC gbc40(gmesh->ElementVec()[4],15,-1);
			TPZGeoElBC gbc41(gmesh->ElementVec()[4],18,-1);
			TPZGeoElBC gbc42(gmesh->ElementVec()[4],19,-3);
			// Fifth sub element - faces: 13 and 15
			TPZGeoElBC gbc50(gmesh->ElementVec()[5],14,-2);
			TPZGeoElBC gbc51(gmesh->ElementVec()[5],16,-4);
			gmesh->ElementVec()[3]->Divide(sub);
			gmesh->ElementVec()[4]->Divide(sub);
		}
			break;
		case 4:
		{
			// First sub element - faces: 15, 16, 18 and 19
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],15,-3);
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],16,-1);
			TPZGeoElBC gbc12(gmesh->ElementVec()[1],18,-4);
			TPZGeoElBC gbc13(gmesh->ElementVec()[1],19,-1);
			// Second sub element - faces: 15, 16, 18 and 19
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],15,-1);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],16,-2);
			TPZGeoElBC gbc22(gmesh->ElementVec()[2],18,-4);
			TPZGeoElBC gbc23(gmesh->ElementVec()[2],19,-3);
			// Third sub element - faces: 15, 17 and 19
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],15,-1);
			TPZGeoElBC gbc31(gmesh->ElementVec()[3],17,-1);
			TPZGeoElBC gbc32(gmesh->ElementVec()[3],19,-3);
			gmesh->ElementVec()[1]->Divide(sub);
			gmesh->ElementVec()[2]->Divide(sub);
			gmesh->ElementVec()[3]->Divide(sub);
		}
			break;
		default:
			return 0;
	}
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	return gmesh;
}

void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,REAL radius,int ntyperefs) {
	TPZManVector<REAL,3> point(3,-0.25);
	REAL r = sqrt(3.0);

    // To refine elements with center near to points than radius, some times as ntyperefs
    RefineGeoElements(dim,gmesh,point,r,radius,ntyperefs);

	// Constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL,3> &point,REAL r,REAL &distance,int ntyperefs) {
	TPZManVector<REAL,3> centerpsi(3), center(3);
    int nsubs, nsubacum, p, k, q;
    REAL centerdist;
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
    TPZVec<TPZGeoEl *> subacum;
	TPZVec<TPZGeoEl *> subsub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
        subacum.Resize(0);
		gel = gmesh->ElementVec()[nelem++];
		if(!gel || gel->Dimension()!=dim || gel->HasSubElement()) continue;
        // element will be divided if any of their nodes is near to circunference
/ *        for(int i=0;i<gel->NCornerNodes();i++) {
            TPZGeoNode* node = gel->NodePtr(i);
            node->GetCoordinates(center);
            centerdist = TPZGeoEl::Distance(center,point);
            if(fabs(r-centerdist) < distance && !gel->NSubElements()) {
                gel->Divide(sub);
                nsubs = gel->NSubElements();
                nsubacum = subacum.NElements();
                subacum.Resize(nsubacum+nsubs,NULL);
                for(p=0;p<nsubs;p++)
                    subacum[nsubacum+p] = sub[p];
            }
        }* /
//        if(!gel->NSubElements()) {
            gel->CenterPoint(gel->NSides()-1,centerpsi);
            gel->X(centerpsi,center);
            centerdist = TPZGeoEl::Distance(center,point);
            if(fabs(r-centerdist) < distance && !gel->NSubElements()) {
                gel->Divide(sub);
                nsubs = gel->NSubElements();
                nsubacum = subacum.NElements();
                subacum.Resize(nsubacum+nsubs,NULL);
                for(p=0;p<nsubs;p++)
                    subacum[nsubacum+p] = sub[p];
            }
  //      }
        int nsubacumtot = subacum.NElements();
        if(ntyperefs>1) {
            for(k=0;k<nsubacumtot;k++) {
                if(!subacum[k]) continue;
                subacum[k]->Divide(subsub);
                if(ntyperefs>2) {
                    TPZVec<TPZGeoEl *> subsubsub;
                    for(q=0;q<subsub.NElements();q++)
                        subsub[q]->Divide(subsubsub);
                }
            }
        }
	}
}

void formatTimeInSec(char *strtime,int timeinsec) {
	if(!strtime) return;
	memset(strtime,0,strlen(strtime));
	//	strtime[0] = '\0';
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
		sprintf(strtime,"%d a, %d m, %d d, %02d:%02d:%02d",anos,meses,dias,horas,minutos,segundos);
	else {
		if(meses) 
			sprintf(strtime,"%d m, %d d, %02d:%02d:%02d",meses,dias,horas,minutos,segundos);
		else {
			if(dias)
				sprintf(strtime,"%d d, %02d:%02d:%02d",dias,horas,minutos,segundos);
			else
				sprintf(strtime,"%02d:%02d:%02d",horas,minutos,segundos);
		}
	}
}


void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename) {
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

void InitializeSolver(TPZAnalysis &an) {
	TPZStepSolver<STATE> step;
	TPZBandStructMatrix matrix(an.Mesh());
	an.SetStructuralMatrix(matrix);
	step.SetDirect(ELU);
	an.SetSolver(step);
}


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
				if(gel->MaterialId() == matidtodivided){
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

*/

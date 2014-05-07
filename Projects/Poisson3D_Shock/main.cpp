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
#include <math.h>

#include <fstream>
#include <cmath>


#include "problem.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;


/**  Global variables  */
int gDebug = 0;
bool usethreads = false;
// Maximum number of equations allowed
long MaxEquations = 1500000;
// Input - output
ofstream out("OutPoissonArcTan.txt",ios::app);             // To store output of the console
// ABOUT H P ADAPTIVE
int MaxPOrder = 13;     // Maximum order for p refinement allowed
int MaxHLevel = 6;      // Maximum level for h refinement allowed
int MaxHUsed = 0;
int MaxPUsed = 0;
// Poisson problem
STATE ValueK = 100000;
STATE F = sqrt(ValueK);
int ModelDimension;
// Circunference with high gradient - data
TPZManVector<REAL,3> CCircle(3,0.5);
REAL RCircle = 0.25;


/** Functions to construction of geometry of problems */
TPZGeoMesh *CreateGeomMesh(MElementType typeel,int mat,int bc0,int bc1=0,int bc2=0);
TPZGeoMesh *CreateGeomMesh(std::string &nome);
// Crea malla computacional sem forcingfunction quando hasforcingfunction = 0, ou toma diferentes forcingfuncition para diferentes
// valores de hasforcingfunction
TPZCompMesh *CreateComputationalMesh(TPZGeoMesh *gmesh,int dim,int matId,int hasforcingfunction,int bc0,int bc1=0,int bc2=0);
TPZGeoMesh *ConstructingPositiveCube(REAL InitialL,MElementType typeel,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingTetrahedraInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingPrismsInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingPyramidsInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingSeveral3DElementsInCube(REAL InitialL,MElementType typeel,int id_bc0,int id_bc1=0,int id_bc2=0);

/** Fucntions to apply refinement. */
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=true, const int matidtodivided=1);


void PrintGeoMeshInVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);
void PrintGeoMeshAsCompMeshInVTKWithElementIndexAsData(TPZGeoMesh *gmesh,char *filename);

// To save meshes in disk
void SaveCompMesh(TPZCompMesh *cmesh, int timessave,TPZCompMesh *cmeshmodified=NULL,bool check=false);

/** Utilitaries */
void formatTimeInSec(char *strtime,int lenstrtime,int timeinsec);

int DefineDimensionOverElementType(MElementType typeel);
void GetFilenameFromGID(MElementType typeel, std::string &name);


/** PROBLEMS */
bool SolveSymmetricPoissonProblemOnCubeMesh(int itypeel);
bool SolveLaplaceProblemOnLShapeMesh();


/**
 * Get Global L2 Error for solution and the L2 error for each element.
 * Return the maxime L2 error by elements. Also return in MinErrorByElement argument the minime L2 error for all elements of the mesh.
 */
REAL ProcessingError(TPZAnalysis &analysis,TPZVec<REAL> &ervec,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL &MinErrorByElement,REAL &);
void LoadSolutionFirstOrder(TPZCompMesh *cmesh, void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv,TPZVec<STATE> &ddsol));
bool ApplyingStrategyHPAdaptiveBasedOnErrorOfSolution(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int ref,int itypeel,REAL &factorError);
bool ApplyingStrategyHPAdaptiveBasedOnErrorOfSolutionAndGradient(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int ref,int itypeel,REAL &factorError);
void ApplyingStrategyPAdaptiveBasedOnExactSphereSolution(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,int ref);

// Writing a relation between number of degree of freedom and L2 error.
bool PrintResultsInMathematicaFormat(TPZVec<REAL> &ErrrVec,TPZVec<long> &NEquations,std::ostream &fileerrors);

void PrintNRefinementsByType(int nref, long nels,long newnels,TPZVec<long> &counter,ostream &out = std::cout);

void AdjustingOrder(TPZCompMesh *cmesh);
int MaxLevelReached(TPZCompMesh *cmesh);

// MAIN FUNCTION TO NUMERICAL SOLVE WITH AUTO ADAPTIVE HP REFINEMENTS
/** Laplace equation on square 1D 2D 3D - Volker John article 2000 */
int main(int argc,char *argv[]) {
    
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();
//    gRefDBase.InitializeRefPatterns();

	// Getting input data
	int itypeel = 4;
	int count = 0;
	do {
		if(argc > 1)
			itypeel = atoi(argv[count+1]);
		if(itypeel > 7 || itypeel < 2)
			itypeel = 7;
		count++;
		// Solving symmetricPoissonProblem on [0,1]^d with d=1, d=2 and d=3
	    if(!SolveSymmetricPoissonProblemOnCubeMesh(itypeel))
		    return 1;
	} while(count < argc-1);
    
    return 0;
}

bool SolveSymmetricPoissonProblemOnCubeMesh(int itypeel) {
	// Variables

	/** Printing level */
	int gPrintLevel = 0;
	int printingsol = 0;
	int printsave = 0;

	int materialId = 1;
	int id_bc0 = -1;
	int id_bc1 = -2;
	// Generic data for problems to solve
	int NRefs = 50;
	int ninitialrefs = 3;
	// Percent of error permited
	REAL factorError = .3;

	// auxiliar string
	char saida[512];

	// To compute processing times
	time_t sttime;
	time_t endtime;
	int time_elapsed;
	char time_formated[256];
	char * ptime = time_formated;
	memset(time_formated,0,256);
	
	// Output files
	std::ofstream fileerrors("ErrorsHP_Poisson.txt",ios::app);   // To store all errors calculated by TPZAnalysis (PosProcess)
	// Initial message to print computed errors
	time(&sttime);
	ptime = ctime(&sttime);
	fileerrors << "Approximation Error in " << time_formated << std::endl;
	
	int nref = 1;
	int nthread = 2, NThreads = 4;

	// Initializing the auto adaptive process
	TPZVec<REAL> ervec, ErrorVec(100,0.0);
	TPZVec<long> NEquations(100,0L);
	TPZVec<REAL> ervecbyel;
	TPZVec<REAL> gradervecbyel;

	MElementType typeel;

	/** Solving for type of geometric elements */
	typeel = (MElementType)itypeel;
	fileerrors << "\nType of element: " << typeel << endl;
	TPZGeoMesh *gmesh;
	gmesh = CreateGeomMesh(typeel,materialId,id_bc0,id_bc1);
	ModelDimension = DefineDimensionOverElementType(typeel);
			
	// Printing geometric mesh to validate
	if(gDebug) {
		sprintf(saida,"gmesh_%02dD_H%dE%d.vtk",ModelDimension,nref,typeel);
		PrintGeoMeshInVTKWithDimensionAsData(gmesh,saida);
	}

	/** Variable names for post processing */
	TPZStack<std::string> scalnames, vecnames;
	scalnames.Push("POrder");
	scalnames.Push("Solution");

	fileerrors.flush();
	out.flush();

	// Adjusting parameters
    NRefs = 10;

	// Initial uniform refinement or printing solution on mesh with 7-h refinements
	if(printingsol) {
		TPZGeoMesh *gmeshfirst = CreateGeomMesh(typeel,materialId,id_bc0,id_bc1);
		UniformRefinement(7,gmesh,ModelDimension);
		TPZCompEl::SetgOrder(1);
		TPZCompMesh *cmeshfirst = CreateComputationalMesh(gmesh,ModelDimension,materialId,1,id_bc0,id_bc1);
		TPZAnalysis ann(cmeshfirst,false);
		LoadSolutionFirstOrder(cmeshfirst,ExactSolutionArcTangent);
		{
			std::stringstream sut;
			sut << "Poisson" << ModelDimension << "D_MESHINIT_E" << typeel << "H" << std::setprecision(2) << nref << ".vtk";
			ann.DefineGraphMesh(ModelDimension,scalnames,vecnames,sut.str());
		}
		ann.PostProcess(3,ModelDimension);
		long countels = 0;
		for(int ii=0;ii<cmeshfirst->NElements();ii++) {
			if(!cmeshfirst->ElementVec()[ii] || cmeshfirst->ElementVec()[ii]->Dimension()!=ModelDimension) continue;
			countels++;
		}
		out << std::endl << "Number of elements 2D: " << countels << std::endl << std::endl ;
		delete cmeshfirst;
		delete gmeshfirst;
		printingsol = false;
	}
	else
		UniformRefinement(ninitialrefs,gmesh,ModelDimension);
            
	// Creating computational mesh (approximation space and materials)
	int p = 1, pinit;
	MaxPUsed = pinit = p;
	MaxHUsed = 1;
	TPZCompEl::SetgOrder(p);
	TPZCompMesh *cmesh;
	gmesh->SetName("Malha Geometrica original");
	if(gDebug) {
		sprintf(saida,"gmesh_%02dD_H%dE%dIndex.vtk",ModelDimension,nref,typeel);
		PrintGeoMeshAsCompMeshInVTKWithElementIndexAsData(gmesh,saida);
	}
	cmesh = CreateComputationalMesh(gmesh,ModelDimension,materialId,1,id_bc0,id_bc1);     // Forcing function is out 2013_07_25

	// To storing number of equations and errors obtained for all iterations
	ErrorVec.Resize(NRefs);
	ErrorVec.Fill(0.0L);
	NEquations.Resize(NRefs);
	NEquations.Fill(0L);

	int countermesh=0;
	// loop solving iteratively
	for(nref=0;nref<NRefs;nref++) {
		if(printsave > 0) {
#ifdef LOG4CXX
			InitializePZLOG();
#endif
			SaveCompMesh(cmesh,countermesh++);
		}
		out << "\nSolving Poisson problem " << ModelDimension << "D. Refinement: " << nref << " Threads: " << nthread << " TypeElement: " << typeel << endl;
		std::cout << "\nSolving Poisson problem. Type of element: " << typeel << std::endl;
		if(usethreads) {
			if(nref > 5) nthread = 2*NThreads;
			else nthread = NThreads;
		}
				
		// Initializing the generation mesh process
		time(& sttime);
				
		// Introduzing exact solution depending on the case
		// Solving adaptive process
		TPZAnalysis an(cmesh,true);
		an.SetExact(ExactSolutionArcTangent);
		{
			std::stringstream sout;
			sout << "Poisson" << ModelDimension << "D_E" << typeel << "Thr" << nthread << "H" << std::setprecision(2) << nref << "P" << pinit << ".vtk";
			an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
		}
		std::string MeshFileName;
		{
			std::stringstream sout;
			sout << "meshAngle" << ModelDimension << "D_E" << typeel << "Thr" << nthread << "H" << std::setprecision(2) << nref << "P" << pinit << ".vtk";
			MeshFileName = sout.str();
		}
				
		cmesh->SetName("Malha computacional adaptada");
		// Printing geometric and computational mesh
		if(gDebug) {
			cmesh->Reference()->Print(std::cout);
			cmesh->Print(std::cout);
		}
				
		// Solve using symmetric matrix then using Cholesky (direct method)
		TPZParSkylineStructMatrix strmat(cmesh,4);
		an.SetStructuralMatrix(strmat);
		
		TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
		direct->SetDirect(ECholesky);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;
				
		out << "\tRefinement: " << nref << " TypeElement: " << typeel << "NEquations " << cmesh->NEquations() << "\n";
		an.Run();
				
		// Post processing
		if(nref > 8 || !(nref%3) || nref==NRefs-1)
			an.PostProcess(2,ModelDimension);
		if(gDebug) {
			std::ofstream out(MeshFileName.c_str());
			cmesh->LoadReferences();
			TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(),out,false);
		}
                
		// generation mesh process finished
		time(&endtime);
		time_elapsed = endtime - sttime;
		formatTimeInSec(time_formated,256,time_elapsed);
		out << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n";
		fileerrors << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n";
		std::cout << "  Time elapsed " << time_elapsed << " <-> " << time_formated << "\n\n";
				
		REAL MinErrorByElement, MinGradErrorByElement;
		ervecbyel.Resize(0);
		gradervecbyel.Resize(0);
		REAL MaxErrorByElement = ProcessingError(an,ervec,ervecbyel,gradervecbyel,MinErrorByElement,MinGradErrorByElement);
		// Printing obtained errors
		if(ervec[1] > 10. || ervec[1] < 0.) {
			std::cout << "L2 Error is wrong (BIG?!) By! \n\n";
			break;
		}
		ErrorVec[nref] = ervec[1];
		NEquations[nref] = cmesh->NEquations();

		std::cout << "\n NRef " << nref << "\tL2 Error " << ervec[1] << "  NEquations: " << NEquations[nref] << " PUsed " << MaxPUsed << " HMax " << MaxHUsed << std::endl << std::endl;
		out << "\n NRef " << nref << "\tL2 Error " << ervec[1] << "  NEquations: " << NEquations[nref] << " PUsed " << MaxPUsed << " HMax " << MaxHUsed << std::endl << std::endl;
		if(cmesh->NEquations() > MaxEquations) {
			NRefs = nref+1;							// final iteration
			ErrorVec.Resize(NRefs);
			NEquations.Resize(NRefs);
			continue;
		}
		fileerrors.flush();
		out.flush();
		if(NRefs > 1 && nref < (NRefs-1)) {
			out << "\n\nApplying Adaptive Methods... step " << nref << "\n";
			std::cout << "\n\nApplying Adaptive Methods... step " << nref << "\n";
			while(!ApplyingStrategyHPAdaptiveBasedOnErrorOfSolutionAndGradient(cmesh,ervecbyel,gradervecbyel,MaxErrorByElement,MinErrorByElement,MinGradErrorByElement,nref,itypeel,factorError)) {
				factorError -= 0.05;
				out << "\nFactorError\nFactorError\nFactorError\n " << factorError << std::endl;
				if(factorError < 0.05) {
					nref = NRefs;
					break;
				}
			}
		}
		fileerrors.flush();
		out.flush();
		// Sometimes Writing a relation between number of degree of freedom and L2 error.
		if(nref && !(nref%4))
			PrintResultsInMathematicaFormat(ErrorVec,NEquations,fileerrors);
	}
	if(cmesh)
		delete cmesh;
	cmesh = NULL;
	if(gmesh)
		delete gmesh;
	gmesh = NULL;
	// Writing a relation between number of degree of freedom and L2 error.
	if(!PrintResultsInMathematicaFormat(ErrorVec,NEquations,fileerrors))
		std::cout << "\nThe errors and nequations values in Mathematica format was not done.\n";
	
	fileerrors << std::endl << "Finished running for element " << itypeel << std::endl << std::endl;
	fileerrors.close();
	std::cout << std::endl << "\tFinished running for element " << itypeel << std::endl << std::endl;
	out.close();
	return true;
}

bool ApplyingStrategyHPAdaptiveBasedOnErrorOfSolutionAndGradient(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int nref,int itypeel,REAL &factorError) {
	if(!cmesh) return false;
	bool result = true;
	long nels = cmesh->NElements();
    
	TPZVec<long> subels;
	TPZVec<long> subsubels;
	int pelement;
	int level;
	TPZInterpolatedElement *el;
	// To see where the computation is doing
	long index = -1;
	TPZVec<long> counterreftype(50,0);
	long i, ii;
    
	REAL factorGrad = .5;
	REAL factorSmall = .1;
	REAL factorErrorBig = 0.8;
    
	REAL BigError = factorErrorBig*MaxErrorByElement + (1.-factorErrorBig)*MinErrorByElement;
	REAL SmallError = factorSmall*MaxErrorByElement + (1. - factorSmall)*MinErrorByElement;
	REAL MaxGrad = factorGrad*gradervecbyel[nels] + (1.-factorGrad)*MinGrad;
	REAL SmallGrad = factorSmall*gradervecbyel[nels] + (1.-factorSmall)*MinGrad;

    
	REAL LaplacianVal;
	REAL MaxLaplacianVal, MinLaplacianVal;
    
	REAL factorLap = 0.7;
	ComputingMaxLaplacian(cmesh,MaxLaplacianVal,MinLaplacianVal);

    REAL LimitLaplace = factorLap*MaxLaplacianVal + (1.-factorLap)*MinLaplacianVal;
    REAL MediumError = factorError*MaxErrorByElement + (1.-factorError)*MinErrorByElement;

	/* Printing maximum and minimun values of the errors */
	out << "\nErro ->   Max " << MaxErrorByElement << "    Min " << MinErrorByElement << "\nGrad ->   Max " << gradervecbyel[nels] << "   Min " << MinGrad;
	out << "\nMaxGrad " << MaxGrad << "  SmallGrad " << SmallGrad << "    BigError " << BigError << "  SError " << SmallError << "  FactorError " << factorError;
	cout << "\nErro ->   Max " << MaxErrorByElement << "    Min " << MinErrorByElement << "\nGrad ->   Max " << gradervecbyel[nels] << "   Min " << MinGrad;
	cout << "\nMaxGrad " << MaxGrad << "  SmallGrad " << SmallGrad << "    BigError " << BigError << "  SError " << SmallError << "  FactorError " << factorError;
    
	// Applying hp refinement only for elements with dimension as model dimension
	for(i=0L;i<nels;i++) {
		bool hused = false, pused = false;
		subels.Resize(0);
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) {
			counterreftype[0]++;
			continue;
		}
        
		// element data
		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		pelement++;
		index = el->Index();
		level = el->Reference()->Level();
        
		if(!LaplacianValue(el,LaplacianVal))
			DebugStop();

        if(ervecbyel[index] > BigError && level < MaxHLevel) {
			if(gradervecbyel[index] > MaxGrad) {
				bool flag;
				flag = false;
				counterreftype[1]++;
				if(LaplacianVal > LimitLaplace && pelement<MaxPOrder) {
					el->PRefine(pelement);
					pused = true;
					counterreftype[2]++;
					flag = true;
				}
				el->Divide(index,subels);
                level++;
                hused = true;
				if(!flag && level < MaxHLevel) {
					counterreftype[3]++;
					for(ii=0;ii<subels.NElements();ii++) {
						el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[ii]]);
						el->Divide(subels[ii],subsubels);
					}
                    level++;
				}
			}
			else {
                level++;
                hused = true;
				el->Divide(index,subels);
				counterreftype[4]++;
			}
			counterreftype[7]++;
		}
		else if(ervecbyel[index] > MediumError) {
			counterreftype[10]++;
			if((gradervecbyel[index] > MaxGrad) && level < MaxHLevel) {
                level++;
                hused = true;
				counterreftype[11]++;
				el->Divide(index,subels);
			}
			else if(pelement<MaxPOrder) {
				el->PRefine(pelement);
				pused = true;
				counterreftype[12]++;
			}
            else {
                counterreftype[13]++;
            }
		}
		else if(gradervecbyel[index] > SmallGrad || ervecbyel[index] > SmallError) {
			counterreftype[20]++;
			if(pelement < MaxPOrder) {
				el->PRefine(pelement);
				pused = true;
				counterreftype[21]++;
			}
            else
				counterreftype[22]++;
		}

		if(!pused && !hused) {
			counterreftype[40]++;
		}
		if(pused)
			MaxPUsed = (pelement > MaxPUsed) ? pelement : MaxPUsed;
		if(hused)
			MaxHUsed = (level > MaxHUsed) ? level : MaxHUsed;
	}
	cmesh->ExpandSolution();
	if(!counterreftype[10] && !counterreftype[20]) {
		result = false;
	}
    
	// Printing information stored
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype,out);
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype);
	return result;
}

bool ApplyingStrategyHPAdaptiveBasedOnErrorOfSolutionAndGradientMorePOnGrad(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int nref,int itypeel,REAL &factorError) {
	if(!cmesh) return false;
	bool result = true;
	long nels = cmesh->NElements();

	TPZVec<long> subels;
	TPZVec<long> subsubels;
	int pelement;
	int level;
	TPZInterpolatedElement *el;
	// To see where the computation is doing
	long index = -1;
	TPZVec<long> counterreftype(50,0);
	long i, ii;

	REAL factorGrad = .6;
	REAL factorSGrad = .15;
	REAL factorErrorBig = 0.8;

	REAL BigError = factorErrorBig*MaxErrorByElement + (1.-factorErrorBig)*MinErrorByElement;
	REAL SmallError = factorError*MaxErrorByElement + (1.-factorError)*MinErrorByElement;
	REAL MaxGrad = factorGrad*gradervecbyel[nels] + (1.-factorGrad)*MinGrad;
	REAL SmallGrad = factorSGrad*gradervecbyel[nels] + (1.-factorSGrad)*MinGrad;

	/* Printing maximum and minimun values of the errors */
	out << "\nErro ->   Max " << MaxErrorByElement << "    Min " << MinErrorByElement << "\nGrad ->   Max " << gradervecbyel[nels] << "   Min " << MinGrad;
	out << "\nMaxGrad " << MaxGrad << "  SmallGrad " << SmallGrad << "    BigError " << BigError << "  SError " << SmallError << "  FactorError " << factorError;
	cout << "\nErro ->   Max " << MaxErrorByElement << "    Min " << MinErrorByElement << "\nGrad ->   Max " << gradervecbyel[nels] << "   Min " << MinGrad;
	cout << "\nMaxGrad " << MaxGrad << "  SmallGrad " << SmallGrad << "    BigError " << BigError << "  SError " << SmallError << "  FactorError " << factorError;

	// Applying hp refinement only for elements with dimension as model dimension
	for(i=0L;i<nels;i++) {
		bool hused = false, pused = false;
		subels.Resize(0);
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) {
			counterreftype[0]++;
			continue;
		}

		// element data
		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		pelement++;
		index = el->Index();
		level = el->Reference()->Level();

		if(nref < 6 && (gradervecbyel[i] > MaxGrad || ervecbyel[i] > BigError) && level < MaxHLevel) {
			counterreftype[10]++;
			el->Divide(index,subels);
			el = NULL;
			level++;
			hused = true;
			if(nref && pelement < MaxPOrder) {
				counterreftype[12]++;
				for(ii=0;ii<subels.NElements();ii++) {
					dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[ii]])->PRefine(pelement);
				}
				pused = true;
			}
			if(gradervecbyel[i] > MaxGrad && level < MaxHLevel) {
				counterreftype[11]++;
				for(ii=0;ii<subels.NElements();ii++) {
					cmesh->ElementVec()[subels[ii]]->Divide(subels[ii],subsubels);
					subels[ii]=0;
				}
				level++;
			}
		}
		else if((gradervecbyel[i] > SmallGrad || ervecbyel[i] > SmallError) && pelement < MaxPOrder && nref) {
			counterreftype[20]++;
			el->PRefine(pelement);
			pused = true;
		}
		else {
			counterreftype[30]++;
		}

		if(!pused && !hused) {
			counterreftype[40]++;
		}
		if(pused)
			MaxPUsed = (pelement > MaxPUsed) ? pelement : MaxPUsed;
		if(hused)
			MaxHUsed = (level > MaxHUsed) ? level : MaxHUsed;
	}
	cmesh->ExpandSolution();
	if(!counterreftype[10] && !counterreftype[20]) {
		result = false;
	}

	// Printing information stored
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype,out);
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype);
	return result;
}

bool ApplyingStrategyHPAdaptiveBasedOnErrorOfSolution(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int nref,int itypeel,REAL &factorError) {
	if(!cmesh) return false;
	bool result = true;
	long nels = cmesh->NElements();
	int dim = cmesh->Dimension();

	TPZVec<long> subels;
	TPZVec<long> subsubels;
	int pelement;
	int level;
	TPZInterpolatedElement *el;
	// To see where the computation is doing
	long index = -1;
	TPZVec<long> counterreftype(50,0);
	long i, ii;

	REAL factorGrad = .2;
	if(dim == 2) factorGrad = 1./3.;
	
	REAL SmallError = factorError*MaxErrorByElement + (1.-factorError)*MinErrorByElement;
	REAL MaxGrad = factorGrad*gradervecbyel[nels] + (1.-factorGrad)*MinGrad;

	TPZVec<REAL> Laplacian(1);
	TPZFMatrix<REAL> dLap(3);
	TPZVec<REAL> psi(3,0.);
	TPZVec<REAL> center(3,0.);

	/* Printing maximum and minimun values of the errors */
	out << "\nErro ->   Min " << MinErrorByElement << "    Max " << MaxErrorByElement << std::endl << "Grad ->   Min " << MinGrad << "   Max " << gradervecbyel[nels] << "\t";
	out << "\nMaxGrad " << MaxGrad << " Factor " << factorGrad << "     SError " << SmallError << " Factor " << factorError;

	// Applying hp refinement only for elements with dimension as model dimension
	for(i=0L;i<nels;i++) {
		bool hused = false, pused = false;
		subels.Resize(0);
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) {
			counterreftype[0]++;
			continue;
		}

		// Getting center of the element
		el->Reference()->CenterPoint(el->Reference()->NSides()-1,psi);
		el->Reference()->X(psi,center);
		RightTermArcTangent(center,Laplacian,dLap);
		Laplacian[0] /= ValueK;

		// element data
		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		pelement++;
		index = el->Index();
		level = el->Reference()->Level();

		// Applying hp refinement depends on high gradient and high laplacian value, and depends on computed error by element
		if(gradervecbyel[i] > MaxGrad && level < MaxHLevel) {
			counterreftype[10]++;
			el->Divide(index,subels);
			el = NULL;
			level++;
			hused = true;
		}
		if((Laplacian[0] > 10. || ervecbyel[i] > SmallError) && pelement < MaxPOrder && nref) {
			counterreftype[20]++;
			if(el)
				el->PRefine(pelement);
			else {
				for(ii=0;ii<subels.NElements();ii++) {
					dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[ii]])->PRefine(pelement);
				}
			}
			pused = true;
		}
		if(!pused && !hused) {
			counterreftype[30]++;
		}
		if(pused)
			MaxPUsed = (pelement > MaxPUsed) ? pelement : MaxPUsed;
		if(hused)
			MaxHUsed = (level > MaxHUsed) ? level : MaxHUsed;
	}
	cmesh->ExpandSolution();
	if(!counterreftype[10] && !counterreftype[20]) {
		result = false;
	}

	// Printing information stored
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype,out);
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype);
	return result;
}



int MaxLevelReached(TPZCompMesh *cmesh) {
	int levelreached = 0, level;
	TPZCompEl *el = 0;
	for(int i=0;i<cmesh->NElements();i++) {
		el = cmesh->ElementVec()[i];
		if(!el || el->Dimension() != cmesh->Dimension()) continue;
		level = el->Reference()->Level();
		levelreached = (level > levelreached) ? level : levelreached; 
	}
	return levelreached;
}

void PrintNRefinementsByType(int nref, long nels,long newnels,TPZVec<long> &counter,ostream &out) {
	out << "\nHP Refinement done, on  " << nels << " elements, given " << newnels << " elements. "<< std::endl;
	out << "NRef = " << nref << std::endl;
	for(int j=0;j<counter.NElements();j++)
		if(counter[j]) {
			out << "Refinement type " << j << " : " << counter[j] << std::endl;
		}
	out << "Processed elements " << (nels-counter[0]);
}
/**
 * Get Global L2 Error for solution and the L2 error for each element.
 * Return the maxime L2 error by elements. Also return in MinErrorByElement argument the minime L2 error for all elements of the mesh.
 */

REAL ProcessingError(TPZAnalysis &analysis,TPZVec<REAL> &ervec,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL &MinErrorByElement,REAL &MinGradErrorByElement) {
    long neq = analysis.Mesh()->NEquations();
	if(ModelDimension != analysis.Mesh()->Dimension())
		DebugStop();
    TPZVec<REAL> ux(neq);
    TPZVec<REAL> sigx(neq);
    TPZManVector<REAL,10> values(10,0.);
    analysis.Mesh()->LoadSolution(analysis.Solution());

	TPZAdmChunkVector<TPZCompEl *> elvec = analysis.Mesh()->ElementVec();
    TPZManVector<REAL,10> errors(10);
    errors.Fill(0.0);
    long i, nel = elvec.NElements();
	ervecbyel.Resize(nel,0.0);
	// The last position will be store the maxime value of the gradient errors
	gradervecbyel.Resize(nel+1,0.0);
	REAL maxError = 0.0;
	MinErrorByElement = 1000.0;
	MinGradErrorByElement = 10000.0;

	/** Computing error for all elements with same dimension of the model */
    for(i=0L;i<nel;i++) {
        TPZCompEl *el = (TPZCompEl *) elvec[i];
		if(!el || el->Dimension() != ModelDimension) continue;
        if(el) {
            errors.Fill(0.0);
            el->EvaluateError(analysis.fExact, errors, 0);
            int nerrors = errors.NElements();
            values.Resize(nerrors, 0.);
            for(int ier = 0; ier < nerrors; ier++)
                values[ier] += errors[ier] * errors[ier];
			// L2 error for each element
			ervecbyel[i] = sqrt(errors[1]*errors[1]);
			gradervecbyel[i] = sqrt(errors[2]*errors[2]);
			if(gradervecbyel[i] > gradervecbyel[nel])
				gradervecbyel[nel] = gradervecbyel[i];
			if(gradervecbyel[i] < MinGradErrorByElement)
				MinGradErrorByElement = gradervecbyel[i];
			// The computed error by current element is compared with max and min values to return
			if(ervecbyel[i] > maxError)
				maxError = ervecbyel[i];
			else if(ervecbyel[i] < MinErrorByElement)
				MinErrorByElement = ervecbyel[i];
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
// Writing a relation between number of degree of freedom and L2 error.
bool PrintResultsInMathematicaFormat(TPZVec<REAL> &ErrorVec,TPZVec<long> &NEquations,std::ostream &fileerrors) {
	int nref;
	long NRefs = ErrorVec.NElements();
	// setting format for ostream
	fileerrors << setprecision(13);
	fileerrors.setf(std::ios::fixed, std::ios::floatfield);
	fileerrors << "\n\nNEquations = {";

	// printing number of equations into a list
	for(nref=0;nref<NRefs-1;nref++) {
		fileerrors << NEquations[nref] << ", ";
	}
	fileerrors << NEquations[nref] << "};" << std::endl << "L2Error = {";
	// printing error values into a list
	for(nref=0;nref<NRefs-1;nref++) {
		fileerrors << ErrorVec[nref] << ", ";
	}
	fileerrors << ErrorVec[nref] << "};";
	// printing lines to create lists of logarithms
	fileerrors << std::endl << "LogNEquations = Table[Log[10,NEquations[[i]]],{i,1,Length[NEquations]}];" << std::endl;
	fileerrors << "LogL2Errors = Table[Log[10,L2Error[[i]]],{i,1,Length[L2Error]}];" << std::endl;
	fileerrors << "ListPlot[{Table[{LogNEquations[[i]],LogL2Errors[[i]]},{i,1,Length[LogNEquations]}]";
	fileerrors << "},Joined->True,PlotRange->All]\n" << std::endl;
	return true;
}

/*************************************************/
/**** Creating Geometric Mesh in cube or square **/
/*************************************************/
TPZGeoMesh *CreateGeomMesh(MElementType typeel,int materialId,int id_bc0,int id_bc1,int id_bc2) {
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
			gmesh = ConstructingTetrahedraInCube(1.,materialId,id_bc0,id_bc1);
			break;
		case EPrisma:
			gmesh = ConstructingPrismsInCube(1.,materialId,id_bc0,id_bc1);
			break;
		case EPiramide:
			gmesh = ConstructingPyramidsInCube(1.,materialId,id_bc0,id_bc1);
			break;
		case ECube:
			gmesh = ConstructingPositiveCube(1.,typeel,materialId,id_bc0,id_bc1);
			break;
		default:
            gmesh = 0;
			break;
	}

	return gmesh;
}

#include "TPZRefPatternDataBase.h"
TPZGeoMesh *ConstructingPositiveCube(REAL InitialL,MElementType typeel,int materialId,int id_bc0,int id_bc1,int id_bc2) {
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
	std::string filename = "D:\\";
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
	
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	return gmesh;
}

TPZGeoMesh *ConstructingTetrahedraInCube(REAL InitialL,int materialId,int id_bc0,int id_bc1,int id_bc2) {
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
	return gmesh;
}

TPZGeoMesh *ConstructingPyramidsInCube(REAL InitialL,int materialId,int id_bc0,int id_bc1,int id_bc2) {
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

TPZGeoMesh *ConstructingPrismsInCube(REAL InitialL,int materialId,int id_bc0,int id_bc1,int id_bc2) {
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
TPZGeoMesh *ConstructingSeveral3DElementsInCube(REAL InitialL,MElementType typeel,int id_bc0,int id_bc1,int id_bc2) {
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
TPZGeoMesh *CreateGeomMesh(std::string &archivo) {
	
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
TPZCompMesh *CreateComputationalMesh(TPZGeoMesh *gmesh,int dim,int materialId,int hasforcingfunction,int id_bc0,int id_bc1,int id_bc2) {
    
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
			mat->SetForcingFunction(new TPZDummyFunction<STATE>(RightTermArcTangent));
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





/** Auxiliar functions  */

void ApplyingStrategyPAdaptiveBasedOnExactSphereSolution(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,int nref) {
    
	if(!cmesh) return;
	long nels = cmesh->NElements();
	TPZVec<long> subels;
	int pelement, dp;
	dp = 1;
	TPZVec<long> subsubels;
	TPZInterpolatedElement *el;
	STATE Tol;
	ZeroTolerance(Tol);
    
	// To see where the computation is doing
	long index = -1;
	TPZVec<long> counterreftype(30,0);
	REAL GradError, SolError;
	long i;
	//	REAL IncrementError = MaxErrorByElement-MinErrorByElement;
	REAL factorGrad = 0.5;
	REAL factorErrorLower = 0.1;
	REAL LaplacianValue, GradNorm;
	
	REAL MaxGradErrorByElement = gradervecbyel[nels];
	std::cout << "\nErroMax " << MaxErrorByElement << "   GradError " << MaxGradErrorByElement << "\n";
	
	// Applying hp refinement only for elements with dimension as model dimension
	for(i=0L;i<nels;i++) {
		bool pused = false;
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) continue;
		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		index = el->Index();
		if(index != i)
			DebugStop();
		// If the element error is little enough do nothing
		if(ervecbyel[index] < (0.1*Tol)) {
			counterreftype[0]++;
			continue;
		}
		
		// If error is small and laplacian value is very little then the order will be minimized
		if(!GradientAndLaplacianOnCorners(el,GradNorm,LaplacianValue))
			DebugStop();
		// Applying hp refinement depends on high gradient and high laplacian value, and depends on computed error by element
        pelement++;
		GradError = gradervecbyel[i];
		SolError = ervecbyel[i];
		if(SolError > factorErrorLower*MaxErrorByElement) {
			if(pelement<MaxPOrder) {
				el->PRefine(pelement);
				pused = true;
				counterreftype[1]++;
			}
		}
		if(pused)
			MaxPUsed = (pelement > MaxPUsed) ? pelement : MaxPUsed;
	}
	cmesh->ExpandSolution();
	// Printing information stored
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype,out);
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype);
}


void AdjustingOrder(TPZCompMesh *cmesh) {
	if(!cmesh) return;
	long nels = cmesh->NElements();
	TPZInterpolatedElement *el;
	STATE Tol;
	ZeroTolerance(Tol);

	// To see where the computation is doing
	long i;
	int level, level0 = 100;
	// Searching highest level of computational elements
	for(i=0L;i<nels;i++) {
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) continue;
		level = el->Reference()->Level();
		level0 = (level < level0) ? level : level0;
	}
	// Applying hp refinement only for elements with dimension as model dimension
	for(i=0L;i<nels;i++) {
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) continue;
		level = el->Reference()->Level();
//		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		if(level-level0+1 < MaxPOrder)
			el->PRefine(level-level0+1);
		else
			el->PRefine(MaxPOrder);
	}
	// Adjusting boundary elements
	cmesh->AdjustBoundaryElements();
}


/**
 * Criteria: Given GlobalL2Error = GE and MaxError (ME) over all the elements
 * If ElementError(EE) > 0.75*ME => twice h-refinement and p-2
 * Else if EE > 0.5*ME  =>  h-refinement and p--
 * Else if EE > 0.25*ME  =>  h-refinement
 * in the other hand => p++
 */
void ApplyingStrategyHPAdaptiveBasedOnErrors(TPZAnalysis &analysis,REAL GlobalL2Error,TPZVec<REAL> &ervecbyel) {

	int ninitialrefs = 2;
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
            STATE GradNorm = GradientNorm(el);
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

void LoadSolutionFirstOrder(TPZCompMesh *cmesh, void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv,TPZVec<STATE> &ddsol)) {
	TPZFMatrix<STATE> solution(cmesh->Solution());
	solution.Zero();
	TPZVec<STATE> sol(1);
	TPZFMatrix<STATE> dsol(cmesh->Dimension(),1);
	TPZVec<STATE> ddsol(9,0.0);
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
			f(coordinates,sol,dsol,ddsol);
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
//        filenamethis.append(".txt");
        fstrthis.OpenWrite(filenamethis);
        
        // Renaming the geometric mesh
        std::stringstream gout;
        gout << (void*)cmesh->Reference();
        cmesh->Reference()->SetName(gout.str());
        
        // Save geometric mesh data
        int classid = cmesh->Reference()->ClassId();
        fstrthis.Write(&classid,1);   // this first data is necessary to use TPZSaveable::Restore
        cmesh->Reference()->Write(fstrthis,0);
        // Save computational mesh data
        classid = cmesh->ClassId();
        fstrthis.Write(&classid,1);   // this first data is necessary to use TPZSaveable::Restore
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

/*
 void ApplyingStrategyHPAdaptiveBasedOnExactCircleSolution(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,int nref) {
 
 if(!cmesh) return;
 long nels = cmesh->NElements();
 TPZVec<long> subels;
 int pelement, dp;
 dp = 1;
 TPZVec<long> subsubels;
 TPZInterpolatedElement *el;
 STATE Tol;
 ZeroTolerance(Tol);
 
 // To see where the computation is doing
 long index = -1;
 TPZVec<long> counterreftype(50,0);
 REAL GradNorm, LaplacianValue;
 REAL MaxGrad, MaxLaplacian;
 long i;
 REAL IncrementError = MaxErrorByElement-MinErrorByElement;
 REAL factorGrad= 0.3;
 REAL factorLap = 0.7;
 REAL factorError = 0.2;
 REAL factorErrorM = 0.8;
 if(nref>1) {
 factorError += (nref-2)*0.08;
 factorErrorM += (nref-2)*0.02;
 }
 if(2<nref)
 factorGrad += (nref-2)*0.1;
 if(factorGrad>0.9)
 factorGrad = 0.95;
 
 ComputingMaxGradientAndLaplacian(cmesh,MaxGrad,MaxLaplacian);
 MaxGrad = gradervecbyel[nels];
 
 // Applying hp refinement only for elements with dimension as model dimension
 for(i=0L;i<nels;i++) {
 bool pused = false;
 el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
 if(!el || el->Dimension()!=cmesh->Dimension()) continue;
 pelement = el->PreferredSideOrder(el->NConnects() - 1);
 index = el->Index();
 // If the element error is little enough do nothing
 if(ervecbyel[i] < (0.1*Tol)) {
 counterreftype[0]++;
 continue;
 }
 
 // If error is small and laplacian value is very little then the order will be minimized
 if(!GradientAndLaplacianOnCorners(el,GradNorm,LaplacianValue))
 DebugStop();
 
 // Applying hp refinement depends on high gradient and high laplacian value, and depends on computed error by element
 pelement++;
 if(ervecbyel[index] > factorErrorM*MaxErrorByElement && IncrementError > 10*Tol) {
 if(gradervecbyel[i] > factorGrad*MaxGrad) {
 bool flag;
 flag = false;
 counterreftype[1]++;
 if(LaplacianValue > factorLap*MaxLaplacian && pelement<MaxPOrder) {
 el->PRefine(pelement);
 pused = true;
 counterreftype[2]++;
 flag = true;
 }
 el->Divide(index,subels);
 if(!flag) {
 counterreftype[3]++;
 for(int ii=0;ii<subels.NElements();ii++) {
 el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[ii]]);
 el->Divide(subels[ii],subsubels);
 }
 }
 }
 else {
 el->Divide(index,subels);
 counterreftype[4]++;
 }
 counterreftype[7]++;
 }
 else if(ervecbyel[index] > factorError*MaxErrorByElement) {
 counterreftype[10]++;
 if((gradervecbyel[i] > factorGrad*MaxGrad)) {
 counterreftype[11]++;
 el->Divide(index,subels);
 }
 else if(pelement<MaxPOrder) {
 el->PRefine(pelement);
 pused = true;
 counterreftype[12]++;
 }
 else {
 counterreftype[13]++;
 el->Divide(index,subels);
 }
 }
 else {
 counterreftype[20]++;
 if(pelement < MaxPOrder) {
 el->PRefine(pelement);
 pused = true;
 counterreftype[21]++;
 }
 //			else if(nref<8) {
 //			el->Divide(index,subels);
 //		counterreftype[9]++;
 //}
 else
 counterreftype[22]++;
 }
 if(pused)
 MaxPUsed = (pelement > MaxPUsed) ? pelement : MaxPUsed;
 }
 cmesh->ExpandSolution();
 
 std::cout << "\nMaxLaplacian " << MaxLaplacian << "\n";
 
 // Printing information stored
 PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype,out);
 PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype);
 }
 
 /*void ApplyingStrategyHPAdaptiveBasedOnExactCircleSolution(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,int nref) {
 
 if(!cmesh) return;
 long nels = cmesh->NElements();
 TPZVec<long> subels;
 int j, pelement, dp;
 dp = 1;
 TPZVec<long> subsubels;
 TPZInterpolatedElement *el;
 STATE Tol;
 ZeroTolerance(Tol);
 
 // To see where the computation is doing
 long index = -1;
 TPZVec<long> counterreftype(50,0);
 REAL GradNorm, LaplacianValue;
 REAL MaxGrad, MaxLaplacian, MaxGradErVecByEl = gradervecbyel[nels];
 long i;
 //	REAL IncrementError = MaxErrorByElement-MinErrorByElement;
 REAL factorGrad = 0.7;
 REAL factorErrorLower = 0.25;
 if(nref==1) {
 factorErrorLower = 0.2;
 }
 if(nref>2) {
 factorGrad = 0.7-(nref-2)*0.1;
 factorErrorLower = 0.1;
 }
 if(factorGrad < 0.3)
 factorGrad = 0.3;
 //	REAL factorLap = 1. + nref;
 REAL factorErrorHigh = 0.75 + nref*0.02;
 
 ComputingMaxGradientAndLaplacian(cmesh,MaxGrad,MaxLaplacian);
 
 // Applying hp refinement only for elements with dimension as model dimension
 for(i=0L;i<nels;i++) {
 bool pused = false;
 el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
 if(!el || el->Dimension()!=cmesh->Dimension()) continue;
 pelement = el->PreferredSideOrder(el->NConnects() - 1);
 index = el->Index();
 // If the element error is little enough do nothing
 if(ervecbyel[i] < (0.1*Tol)) {
 counterreftype[0]++;
 continue;
 }
 
 // If error is small and laplacian value is very little then the order will be minimized
 if(!GradientAndLaplacianOnCorners(el,GradNorm,LaplacianValue))
 DebugStop();
 
 // Applying hp refinement depends on high gradient and high laplacian value, and depends on computed error by element
 pelement++;
 
 if(ervecbyel[i] > factorErrorHigh*MaxErrorByElement) {
 bool flag = false;
 if(LaplacianValue > 1. && pelement < MaxPOrder) {
 el->PRefine(pelement);
 pused = true;
 flag = true;
 counterreftype[12]++;
 }
 counterreftype[10]++;
 el->Divide(index,subels);
 if(GradNorm > factorGrad*MaxGrad || gradervecbyel[i] > factorGrad*MaxGradErVecByEl) {
 counterreftype[11]++;
 for(int ii=0;ii<subels.NElements();ii++) {
 el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[ii]]);
 el->Divide(subels[ii],subsubels);
 }
 }
 }
 else if(ervecbyel[i] > factorErrorLower*MaxErrorByElement) {
 if(LaplacianValue > 1. && pelement < MaxPOrder) {
 el->PRefine(pelement);
 pused = true;
 counterreftype[22]++;
 }
 else if(gradervecbyel[i] > factorGrad*MaxGradErVecByEl) {
 counterreftype[20]++;
 el->Divide(index,subels);
 }
 }
 else
 counterreftype[30]++;
 
 if(pused)
 MaxPUsed = (pelement > MaxPUsed) ? pelement : MaxPUsed;
 }
 cmesh->ExpandSolution();
 
 std::cout << "\nMaxLaplacian " << MaxLaplacian << "\n";
 
 // Printing information stored
 PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype,out);
 PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype);
 }

 void ApplyingStrategyPAdaptiveBasedOnExactSphereSolution_old(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,int nref) {
    
	if(!cmesh) return;
	long nels = cmesh->NElements();
	TPZVec<long> subels;
	int j, pelement, dp;
	dp = 1;
	TPZVec<long> subsubels;
	TPZInterpolatedElement *el;
	STATE Tol;
	ZeroTolerance(Tol);
    
	// To see where the computation is doing
	long index = -1;
	TPZVec<long> counterreftype(30,0);
	REAL GradError, SolError;
	long i;
	//	REAL IncrementError = MaxErrorByElement-MinErrorByElement;
	REAL factorGrad = 0.5;
	REAL factorErrorLower = 0.1;
	REAL LaplacianValue, GradNorm;
	
	//	if(2<nref)
	//	factorGrad += 0.1;
	
	REAL MaxGradErrorByElement = gradervecbyel[nels];
	std::cout << "\nErroMax " << MaxErrorByElement << "   GradError " << MaxGradErrorByElement << "\n";
	
	// Applying hp refinement only for elements with dimension as model dimension
	for(i=0L;i<nels;i++) {
		bool pused = false;
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) continue;
		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		index = el->Index();
		if(index != i)
			DebugStop();
		// If the element error is little enough do nothing
		if(ervecbyel[index] < (0.1*Tol)) {
			counterreftype[0]++;
			continue;
		}
		
		// If error is small and laplacian value is very little then the order will be minimized
		if(!GradientAndLaplacianOnCorners(el,GradNorm,LaplacianValue))
			DebugStop();
		// Applying hp refinement depends on high gradient and high laplacian value, and depends on computed error by element
        pelement++;
		GradError = gradervecbyel[i];
		SolError = ervecbyel[i];
		if(SolError > factorErrorLower*MaxErrorByElement) {
			if(pelement<MaxPOrder) {
				el->PRefine(pelement);
				pused = true;
				counterreftype[1]++;
			}
		}
		if(pused)
			MaxPUsed = (pelement > MaxPUsed) ? pelement : MaxPUsed;
	}
	cmesh->ExpandSolution();
	// Printing information stored
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype,out);
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype);
}
bool ApplyingStrategyHPAdaptiveBasedOnErrorOfSolutionAndGradient(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int nref,int itypeel,REAL &factorError) {
	if(!cmesh) return false;
	bool result = true;
	long nels = cmesh->NElements();
	int dim = cmesh->Dimension();

	TPZVec<long> subels;
	TPZVec<long> subsubels;
	int pelement;
	int level;
	TPZInterpolatedElement *el;
	// To see where the computation is doing
	long index = -1;
	TPZVec<long> counterreftype(50,0);
	long i, ii;

	REAL factorGrad = .6;
	REAL factorErrorBig = 0.85;
	if(MinGrad < 0 || MinErrorByElement < 0) DebugStop();
	REAL BigError = factorErrorBig*MaxErrorByElement;
	REAL SmallError = factorError*MaxErrorByElement;
	REAL MaxGrad = factorGrad*gradervecbyel[nels];
	REAL LaplacianValue;

	TPZVec<REAL> Laplacian(1);
	TPZFMatrix<REAL> dLap(3);
	TPZVec<REAL> psi(3,0.);
	TPZVec<REAL> center(3,0.);

	// Printing maximum and minimun values of the errors 
	out << "\nErro ->   Min " << MinErrorByElement << "    Max " << MaxErrorByElement << std::endl << "Grad ->   Min " << MinGrad << "   Max " << gradervecbyel[nels] << "\t";
	out << "\nMaxGrad " << MaxGrad << " Factor " << factorGrad << "     SError " << SmallError << " Factor " << factorError;

	// Applying hp refinement only for elements with dimension as model dimension
	for(i=0L;i<nels;i++) {
		bool hused = false, pused = false;
		subels.Resize(0);
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) {
			counterreftype[0]++;
			continue;
		}

		// Getting center of the element
		el->Reference()->CenterPoint(el->Reference()->NSides()-1,psi);
		el->Reference()->X(psi,center);
		RightTermArcTangent(center,Laplacian,dLap);
		LaplacianValue = fabs(Laplacian[0])/ValueK;

		// element data
		pelement = el->PreferredSideOrder(el->NConnects() - 1);
		pelement++;
		index = el->Index();
		level = el->Reference()->Level();

		if(ervecbyel[i] > BigError && level < MaxHLevel) {
			counterreftype[10]++;
			el->Divide(index,subels);
			el = NULL;
			level++;
			hused = true;
			if(gradervecbyel[i] > MaxGrad && level < MaxHLevel) {
				counterreftype[11]++;
				for(ii=0;ii<subels.NElements();ii++) {
					cmesh->ElementVec()[subels[ii]]->Divide(subels[ii],subsubels);
					subels[ii]=0;
				}
				level++;
			}
			else {
				counterreftype[12]++;
				for(ii=0;ii<subels.NElements();ii++) {
					dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[subels[ii]])->PRefine(pelement);
				}
				pused = true;
			}
		}
		else if(ervecbyel[i] > SmallError && pelement < MaxPOrder && nref) {
			counterreftype[20]++;
			if(gradervecbyel[i] > MaxGrad && level < MaxHLevel) {
				counterreftype[21]++;
				el->Divide(index,subels);
				hused = true;
				level++;
			}
			else {
				counterreftype[22]++;
				el->PRefine(pelement);
				pused = true;
			}
		}
		else if(LaplacianValue > 9. && pelement < MaxPOrder && nref) {
			counterreftype[30]++;
			el->PRefine(pelement);
			pused = true;
		}

		if(!pused && !hused) {
			counterreftype[40]++;
		}
		if(pused)
			MaxPUsed = (pelement > MaxPUsed) ? pelement : MaxPUsed;
		if(hused)
			MaxHUsed = (level > MaxHUsed) ? level : MaxHUsed;
	}
	cmesh->ExpandSolution();
	if(!counterreftype[10] && !counterreftype[20] && !counterreftype[30]) {
		result = false;
	}

	// Printing information stored
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype,out);
	PrintNRefinementsByType(nref,nels,cmesh->NElements(),counterreftype);
	return result;
}

*/
 /**
 * @file Poisson 3D in hexahedra with shock problem
 */

#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "pzgnode.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"
#include "tpzgeoelrefpattern.h"

#include "pzintel.h"
#include "pzcompel.h"
#include "pzcmesh.h"

#include "pzmatrix.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzvec_extras.h"

#include "pzlog.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"
#include "pzcheckmesh.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzsbstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"

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

#include "CreateAndRefineMeshes.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;


/*** New Data Type to make simulation case as unity ***/
// Simulation Case
struct SimulationCase {
    int hpcase;     // hpcase = 1 Table On U, hpcase = 2 Table On U e DU as in article, hpcase = 3 Table 2 simplified,
                    // hpcase = 4 Table 2 simplified and improved, hpcase = 5 Table defined on 3x3 regions on error measure,
                    // hpcase = 6 Table 5 improved and uniformed
    bool  IsHdivQ;
    int   n_acc_terms;
    int   eltype;
    int   nthreads;
    std::string  dir_name;
    
    SimulationCase() : hpcase(3), IsHdivQ(false), n_acc_terms(0), eltype(7), nthreads(0) {
    }
    SimulationCase(const SimulationCase &other) : hpcase(other.hpcase), IsHdivQ(other.IsHdivQ), n_acc_terms(other.n_acc_terms), eltype(other.eltype), nthreads(other.nthreads), dir_name(other.dir_name) {
    }
    void SetDirName() {
        std::stringstream sout;
        sout << "E" << eltype << "HPCase" << hpcase;
        dir_name = sout.str();

    }
};


/**  Global variables  */
REAL GlobScale = 1.;
// Maximum number of equations allowed
long MaxEquations = 700000;
// Input - output
ofstream out("OutPoissonArcTan.txt",ios::app);             // To store output of the console
// ABOUT H P ADAPTIVE
int MaxPOrder = 10;     // Maximum order for p refinement allowed
int MaxHLevel = 8;      // Maximum level for h refinement allowed
int MaxHUsed = 0;
int MaxPUsed = 0;

int ninitialrefs = 3;

// Poisson problem
STATE ValueK = 100000;
STATE F = sqrt(ValueK);
int ModelDimension;

// Circunference with high gradient - data
TPZManVector<REAL,3> CCircle(3,0.5);
REAL RCircle = 0.25;

// To run one time
bool Once = false;

//**********   Creating computational mesh with materials    *************
TPZCompMesh *CreateComputationalMesh(TPZGeoMesh *gmesh,int dim,int materialId,int hasforcingfunction,int id_bc0,int id_bc1=0,int id_bc2=0);
//void UnwrapMesh(TPZCompMesh *cmesh);

int DefineDimensionOverElementType(int typeel);
void GetFilenameFromGID(MElementType typeel, std::string &name);

/** PROBLEMS */
bool SolveSymmetricPoissonProblemOnCubeMesh(SimulationCase &sim_case);

/**
 * Get Global L2 Error for solution and the L2 error for each element.
 * Return the maxime L2 error by elements. Also return in MinErrorByElement argument the minime L2 error for all elements of the mesh.
 */
bool ProcessingErrorUAndDUKnowingExactSol(TPZAnalysis &an,TPZVec<REAL> &ErrorVecByIteration,int nref,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU);
//REAL ProcessingError(TPZAnalysis &analysis,TPZVec<REAL> &ervec,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL &MinErrorByElement,REAL &);
void LoadSolutionFirstOrder(TPZCompMesh *cmesh, void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv,TPZVec<STATE> &ddsol));

// HP adaptive for strategies in specific tables
bool ApplyingHPAdaptiveStrategyBasedOnU_I(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol);
bool ApplyingHPAdaptiveStrategyBasedOnUAndDUAsArticle_II(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol);
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_III(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol);
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_IV(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol);
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_V(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol);
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_VI(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol);

// Writing a relation between number of degree of freedom and L2 error.
bool PrintResultsInMathematicaFormat(int typeel,int hpcase,int nref,TPZVec<REAL> &ErrrVec,TPZVec<long> &NEquations,std::ostream &fileerrors);

void AdjustingOrder(TPZCompMesh *cmesh);
int MaxLevelReached(TPZCompMesh *cmesh);


/** Utilitaries Over Date And Time */
void formatTimeInSec(char *strtime,int lenstrtime,int timeinsec);
bool CreateCurrentResultDirectory(SimulationCase &sim);

#ifdef LOG4CXX
static LoggerPtr  logger(Logger::getLogger("pz.refine"));
#endif



// MAIN FUNCTION TO NUMERICAL SOLVE WITH AUTO ADAPTIVE HP REFINEMENTS
/** Laplace equation on square 1D 2D 3D - Volker John article 2000 */

int main(int argc,char *argv[]) {
    
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();

	// Getting input data
	// 1 -> line	2 -> triangles	3 -> quadrilateral  4 -> tetraedro	5 -> pyramid	6 -> prisma     7 -> cubo
    struct SimulationCase dummied;

	// Type of elements
	int itypeel = 3;

    // loop over all element types
    do {
        dummied.hpcase = 1;
        // loop over use of specific strategy hp-adaptive table
        do {
            dummied.eltype = itypeel;
            dummied.SetDirName();
            // Solving symmetricPoissonProblem on [0,1]^d with d=1, d=2 and d=3
            if(!SolveSymmetricPoissonProblemOnCubeMesh(dummied))
                return 1;

            dummied.hpcase++;
        }while(dummied.hpcase < 7);

		itypeel++;
	} while(itypeel < 8 && !Once);
    out.close();
	return 0;
}

bool SolveSymmetricPoissonProblemOnCubeMesh(SimulationCase &sim_case) {
    if(CreateCurrentResultDirectory(sim_case))
        return false;
    
    // To compute processing times
    time_t sttime;
    time_t endtime;
    int time_elapsed;
    char * ptime; // = time_formated;

	// Tolerance for applying hp adaptivity
	TPZManVector<REAL,3> Tol(3, 1.e-8);
    Tol[1] = sqrt(Tol[0]); Tol[2] = sqrt(sqrt(Tol[1]));

	int materialId = 1;
	int id_bc0 = -1;
	int id_bc1 = -2;
    
	// Generic data for problems to solve
	int NRefs = 8;

	// Output files
	std::stringstream sout;
	sout << sim_case.dir_name.c_str() << "/ErrorsHP_Poisson.txt";
	std::ofstream fileerrors(sout.str().c_str());   // To store all errors calculated by TPZAnalysis (PosProcess)
    // Initial message to print computed errors
	time(&sttime);
	ptime = ctime(&sttime);
	fileerrors << "\nApproximation Error in " << ptime << std::endl << "\nType of element: " << sim_case.eltype << endl;

	/** Solving for type of geometric elements */
	TPZGeoMesh *gmesh;
	gmesh = CreateGeomMesh(sim_case.eltype,materialId,id_bc0,id_bc1);
	ModelDimension = DefineDimensionOverElementType(sim_case.eltype);
	UniformRefinement(ninitialrefs, gmesh);

	// Printing initial geometric mesh
    std::stringstream sout1;
    sout1 << sim_case.dir_name.c_str() << "/InitialGMesh_" << ModelDimension << "D_E" << sim_case.eltype << ".vtk";
	ofstream arg2(sout1.str().c_str());
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, arg2);
    arg2.close();

	// Creating computational mesh (approximation space and materials)
	int p = 1, pinit;
	MaxPUsed = pinit = p;
	TPZCompEl::SetgOrder(p);
	TPZCompMesh *cmesh;
	gmesh->SetName("Malha Geometrica original");
    
    cmesh = CreateComputationalMesh(gmesh,ModelDimension,materialId,1,id_bc0,id_bc1);     // Forcing function is out 2013_07_25
	MaxHUsed = MaxLevelReached(cmesh);
	// Initializing the vectors of errors to store the errors for any iteration
	TPZMaterial *mater = cmesh->FindMaterial(materialId);
	int nerros = mater->NEvalErrors();
	TPZVec<REAL> ErrorVecByIteration(nerros*NRefs, 0.0);
	TPZVec<long> NEquations(NRefs, 0L);
	TPZVec<STATE> ErrorU, ErrorDU;
	
    //AdjustFluxPolynomialOrders(cmesh, 0);
    //if(HDiv) ReconstructHDivMesh(cmesh, meshvec, hdivplusplus);

	int nref = 0;
    bool tolachieved = false;

	// loop solving iteratively
	do {
        out << "\n\nSOLVING POISSON PROBLEM " << ModelDimension << "D." << " ELEMENT Type: " << sim_case.eltype << " STRATEGY: " << sim_case.hpcase << " Iteration: " << nref << std::endl;
        std::cout << "\n\nSOLVING POISSON PROBLEM " << ModelDimension << "D." << " ELEMENT Type: " << sim_case.eltype << " STRATEGY: " << sim_case.hpcase << " Iteration: " << nref << std::endl;
		
		// Initializing the generation mesh process
		time(& sttime);
				
		// Introduzing exact solution depending on the case
		// Solving adaptive process
        cmesh->CleanUpUnconnectedNodes();
        ErrorU.Resize(0);
        ErrorDU.Resize(0);
        
		TPZAnalysis an(cmesh,true);
		an.SetExact(ExactSolutionArcTangent);
        
		// Solve using symmetric matrix then using Cholesky (direct method)
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh);
        strmat.SetNumThreads(1);
        strmat.SetDecomposeType(ELDLt);
		
		TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
		direct->SetDirect(ELU);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;
				
		out << "\tRefinement: " << nref << " TypeElement: " << sim_case.eltype << "NEquations " << cmesh->NEquations() << "\n";
		an.Assemble();
        an.Solve();
		
		// Post processing
		/** Variable names for post processing */
		TPZStack<std::string> scalnames, vecnames;
		scalnames.Push("POrder");
		scalnames.Push("Pressure");

		std::stringstream sout3;
        sout3 << sim_case.dir_name.c_str() << "/" << "Poisson" << ModelDimension << "D_E" << sim_case.eltype << "H" << std::setprecision(2) << nref << "P" << pinit << ".vtk";
        an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout3.str().c_str());
		an.PostProcess(1,ModelDimension);
		cmesh->LoadReferences();
        
		// generation mesh process finished
		time(&endtime);
		time_elapsed = endtime - sttime;
		formatTimeInSec(ptime,256,time_elapsed);
		out << " Time elapsed " << time_elapsed << " <-> " << ptime << "\n";
		fileerrors << " Time elapsed " << time_elapsed << " <-> " << ptime << "\n";
		std::cout << " Time elapsed " << time_elapsed << " <-> " << ptime << "\n";

		if(!ProcessingErrorUAndDUKnowingExactSol(an,ErrorVecByIteration,nref,ErrorU,ErrorDU))
            DebugStop();

		NEquations[nref] = cmesh->NEquations();

		std::cout << " NRef " << nref << "\tL2 Error " << ErrorVecByIteration[nerros*nref+1] << "  NEquations: " << NEquations[nref] << " PUsed " << MaxPUsed << " HMax " << MaxHUsed << std::endl << std::endl;
		out << " NRef " << nref << "\tL2 Error " << ErrorVecByIteration[nerros*nref+1] << "  NEquations: " << NEquations[nref] << " PUsed " << MaxPUsed << " HMax " << MaxHUsed << std::endl << std::endl;
		if(cmesh->NEquations() > MaxEquations) {
			NRefs = nref+1;							// final iteration
			ErrorVecByIteration.Resize(nerros*NRefs);
			NEquations.Resize(NRefs);
			continue;
		}

		// HP REFINEMENT PROCESS - HP Case depends tabeled strategy (TABLE pre defined)
        if(sim_case.hpcase == 1)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnU_I(an.Mesh(),ErrorU,ErrorDU,Tol);
        else if(sim_case.hpcase == 2)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDUAsArticle_II(an.Mesh(),ErrorU,ErrorDU,Tol);
        else if(sim_case.hpcase == 3)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_III(an.Mesh(),ErrorU,ErrorDU,Tol);
        else if(sim_case.hpcase == 4)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_IV(an.Mesh(),ErrorU,ErrorDU,Tol);
        else if(sim_case.hpcase == 5)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_V(an.Mesh(),ErrorU,ErrorDU,Tol);
        else
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_VI(an.Mesh(),ErrorU,ErrorDU,Tol);

        out << "\n Applying Adaptive Methods... step " << nref << "\n";
		std::cout << "\n Applying Adaptive Methods... step " << nref << "\n";
		std::cout << " NElements Comp " << cmesh->NElements() << " NEquations " << cmesh->NEquations() << std::endl;
		out << " NElements Comp " << cmesh->NElements() << " NEquations " << cmesh->NEquations() << std::endl;

		// Sometimes Writing a relation between number of degree of freedom and L2 error.
        fileerrors << "H1 approximation\n";

		PrintResultsInMathematicaFormat(sim_case.eltype,sim_case.hpcase,nref,ErrorVecByIteration,NEquations,fileerrors);
        fileerrors.flush();
        fileerrors << "done\n";

		nref++;
	}while (nref < NRefs && !tolachieved);

	if (cmesh)
		delete cmesh;
	cmesh = NULL;
	if(gmesh)
		delete gmesh;
	gmesh = NULL;

	// Writing a relation between number of degree of freedom and L2 error.
    std::stringstream sout4;
//	sout4 << sim_case.dir_name.c_str() << "/ErrorsHP_Poisson.nb";
    sout4 << "ErrorsHP_Poisson.nb";
    std::ofstream finalerrors(sout4.str().c_str(),ios::app);   // To store all errors calculated by TPZAnalysis (PosProcess)
	if(!PrintResultsInMathematicaFormat(sim_case.eltype,sim_case.hpcase,nref,ErrorVecByIteration,NEquations,finalerrors))
		std::cout << "\nThe errors and nequations values in Mathematica format was not done.\n";
    finalerrors.close();
	
	fileerrors << std::endl << "Finished running for element " << sim_case.eltype << std::endl << std::endl;
	fileerrors.close();
	std::cout << std::endl << "\tFinished running for element " << sim_case.eltype << std::endl << std::endl;
	return true;
}

/**
* Get Global L2 Error for solution and the L2 error for each element.
* Return the maxime L2 error by elements. Also return in MinErrorByElement argument the minime L2 error for all elements of the mesh.
*/
bool ProcessingErrorUAndDUKnowingExactSol(TPZAnalysis &analysis, TPZVec<REAL> &ErrorVecByIteration, int nref, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU)
{
	TPZCompMesh *cmesh = analysis.Mesh();
	cmesh->LoadSolution(analysis.Solution());

	if (ModelDimension != cmesh->Dimension())
		DebugStop();

	TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
	TPZManVector<REAL, 10> errors(10);
	errors.Fill(0.0);
	long i, nel = elvec.NElements();
	long nerrors = 0L;
	ErrorU.Resize(nel, 0.0);
	// The last position will be store the maxime value of the gradient errors
	ErrorDU.Resize(nel, 0.0);

	/** Computing error for all elements with same dimension of the model */
	for (i = 0L; i<nel; i++) {
		TPZCompEl *el = (TPZCompEl *)elvec[i];
		if (!el || el->Dimension() != ModelDimension) continue;
		errors.Fill(0.0);
		el->EvaluateError(analysis.fExact, errors, 0);
		nerrors = errors.NElements();
		REAL vol = el->VolumeOfEl();
		for (int ier = 0; ier < nerrors; ier++) {
			errors[ier] *= vol;
			ErrorVecByIteration[nerrors*nref + ier] += errors[ier];
		}

		// L2 error for each element
		ErrorU[i] = errors[1];
		ErrorDU[i] = errors[2];
	}
	return true;
}

void ApplyHPRefinement(TPZCompMesh *cmesh, TPZVec<long> &PRef, TPZVec<long> &HRef) {
	long iel, nelhrefs = HRef.NElements(), nelprefs = PRef.NElements();
	long nels = cmesh->NElements();

	TPZManVector<long, 27> subels;
	TPZManVector<long, 27> subsubels;

	// Doing P Refinement
	int pelement = 0;
	TPZGeoEl *gel = 0;
	TPZInterpolationSpace *intel;
	for (iel = 0; iel<nelprefs; iel++) {
		intel = 0;
		intel = dynamic_cast<TPZInterpolationSpace* > (cmesh->Element(PRef[iel]));
		if (!intel || intel->Dimension() != cmesh->Dimension()) continue;
		pelement = intel->GetPreferredOrder(); //->PreferredSideOrder(gel->NSides() - 1);
		if (pelement < MaxPOrder)
			intel->PRefine(pelement + 1);
	}
	cmesh->ExpandSolution();

	// Doing H Refinement
	for (iel = 0; iel<nelhrefs; iel++) {
		bool twice = false;
		if (HRef[iel]<0) {
			twice = true;
			HRef[iel] *= -1;
		}
		subels.Resize(0);
		intel = 0;
		intel = dynamic_cast<TPZInterpolatedElement* > (cmesh->Element(HRef[iel]));
		if (!intel || intel->Dimension() != cmesh->Dimension()) continue;
		gel = intel->Reference();
		if (!gel) DebugStop();
		if (gel->Level() < MaxHLevel) {
			intel->Divide(intel->Index(), subels,1);
			cmesh->ElementVec().SetFree(HRef[iel]);
			//            intel = 0;
		}
		if (twice) {
			for (long isub_el = 0; isub_el<subels.NElements(); isub_el++) {
				subsubels.Resize(0);
				TPZCompEl * isub_cel = cmesh->ElementVec()[subels[isub_el]];
				if (!isub_cel || isub_cel->Dimension() != cmesh->Dimension()) continue;
				isub_cel->Divide(subels[isub_el], subsubels);
				cmesh->ElementVec().SetFree(subels[isub_el]);
			}
			twice = false;
		}
	}
	cmesh->ExpandSolution();

	// Printing information stored
	PrintNRefinementsByType(nels, cmesh->NElements(), nelhrefs, nelprefs, out);
	PrintNRefinementsByType(nels, cmesh->NElements(), nelhrefs, nelprefs);
}

bool ApplyingHPAdaptiveStrategyBasedOnU_I(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol) {
    if(!cmesh) return false;
    long iel, nelhrefs = 0, nelprefs = 0;
    long nels = cmesh->NElements();
    
    TPZVec<long> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    std::cout << " Refinando malha com " << nels  << " elementos e " << cmesh->NEquations() << " equacoes.\n";
	out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";

    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(ErrorU[iel] > Tol[2]) {
            HRef[nelhrefs++] = iel;
            HRef[nelhrefs] *= -1;
        }
        else if(ErrorU[iel] > Tol[1]) {
            PRef[nelprefs++] = iel;
            HRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] > Tol[0]) {
            PRef[nelprefs++] = iel;
        }
    }
	HRef.Resize(nelhrefs);
	PRef.Resize(nelprefs);

	// Doing h and p refinements
	ApplyHPRefinement(cmesh,PRef,HRef);
    
    // If no exists any element to refine, the tolerance was reached
    if(!nelhrefs && !nelprefs)
        return true;
    return false;
}
bool ApplyingHPAdaptiveStrategyBasedOnUAndDUAsArticle_II(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol) {
    if(!cmesh) return false;
    long iel, nelhrefs = 0, nelprefs = 0;
    long nels = cmesh->NElements();
    
    TPZVec<long> HRef(nels,0L), PRef(nels,0L);
    
    // To know laplacian values as auxiliar information to adaptive
    REAL LaplacianVal;
    REAL MaxLaplacianVal = 0., MinLaplacianVal = 1.e4;
    REAL factorLap = 0.7;
    ComputingMaxLaplacian(cmesh,MaxLaplacianVal,MinLaplacianVal);
    
    REAL LimitLaplace = factorLap*MaxLaplacianVal + (1.-factorLap)*MinLaplacianVal;

    // Applying hp refinement only for elements with dimension as model dimension
	std::cout << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
	out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";

    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(!LaplacianValue(cel,LaplacianVal)){
            DebugStop();
        }

        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[2])
                PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            else
                HRef[nelhrefs++] = iel;
        }
        else {
            HRef[nelhrefs++] = iel;
            if(ErrorDU[iel] > Tol[2] && LaplacianVal < LimitLaplace)
                HRef[nelhrefs] *= -1;
            else
                PRef[nelprefs++] = iel;
        }
    }
    
	HRef.Resize(nelhrefs);
	PRef.Resize(nelprefs);

	// Doing h and p refinements
	ApplyHPRefinement(cmesh, PRef, HRef);

	// If no exists any element to refine, the tolerance was reached
	if (!nelhrefs && !nelprefs)
		return true;
	return false;
}
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_III(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol) {
    if(!cmesh) return false;
    long iel, nelhrefs = 0, nelprefs = 0;
    long nels = cmesh->NElements();
    
    TPZVec<long> HRef(nels,0L), PRef(nels,0L);

    // Applying hp refinement only for elements with dimension as model dimension
    std::cout << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[2])
                PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            else
                HRef[nelhrefs++] = iel;
        }
        else {
            HRef[nelhrefs++] = iel;
            if(ErrorDU[iel] > Tol[2])
                PRef[nelprefs++] = iel;
        }
    }
    
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh, PRef, HRef);
    
    // If no exists any element to refine, the tolerance was reached
    if (!nelhrefs && !nelprefs)
        return true;
    return false;
}
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_IV(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol) {
    if(!cmesh) return false;
    long iel, nelhrefs = 0, nelprefs = 0;
    long nels = cmesh->NElements();
    
    TPZVec<long> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    std::cout << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[2])
                PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            else
                HRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            HRef[nelprefs++] = iel;
            if(ErrorDU[iel] > Tol[2])
                PRef[nelhrefs++] = iel;
        }
        else {
            HRef[nelhrefs++] = iel;
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            else
                HRef[nelhrefs++] *= -1;
        }
    }
    
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh, PRef, HRef);
    
    // If no exists any element to refine, the tolerance was reached
    if (!nelhrefs && !nelprefs)
        return true;
    return false;
}

bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_V(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol) {
	if(!cmesh) return false;
    long iel, nelhrefs = 0, nelprefs = 0;
	long nels = cmesh->NElements();
    
    TPZVec<long> HRef(nels,0L), PRef(nels,0L);

	// Applying hp refinement only for elements with dimension as model dimension
	std::cout << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
	out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";

    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[1])
                PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            else
                HRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            if(ErrorDU[iel] < Tol[1])
                PRef[nelprefs++] = iel;
            else {
                HRef[nelhrefs++] = iel;
                if(ErrorDU[iel] > Tol[2])
                    PRef[nelprefs++] = iel;
            }
        }
        else {
            HRef[nelhrefs++] = iel;
            if(ErrorDU[iel] > Tol[2])
                HRef[nelhrefs] *= -1;
            else if(ErrorDU[iel] > Tol[1])
                PRef[nelprefs++] = iel;
        }
    }
    
	HRef.Resize(nelhrefs);
	PRef.Resize(nelprefs);

	// Doing h and p refinements
	ApplyHPRefinement(cmesh, PRef, HRef);

	// If no exists any element to refine, the tolerance was reached
	if (!nelhrefs && !nelprefs)
		return true;
	return false;
}
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_VI(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol) {
    if(!cmesh) return false;
    long iel, nelhrefs = 0, nelprefs = 0;
    long nels = cmesh->NElements();
    
    TPZVec<long> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    std::cout << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[2])
                HRef[nelprefs++] = iel;
            else if(ErrorDU[iel] > Tol[1])
                PRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            if(ErrorDU[iel] > Tol[1]) {
                HRef[nelprefs++] = iel;
                if(ErrorDU[iel] > Tol[2])
                    PRef[nelhrefs++] = iel;
            }
            else
                PRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            HRef[nelhrefs++] = iel;
            if(ErrorDU[iel] > Tol[2])
                HRef[nelprefs++] *= -1;
            else if(ErrorDU[iel] > Tol[1])
                PRef[nelhrefs++] = iel;
        }
        else {
            HRef[nelhrefs++] = iel;
            if(ErrorDU[iel] > Tol[1])
                HRef[nelhrefs] *= -1;
            else
                PRef[nelprefs++] = iel;
        }
    }
    
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh, PRef, HRef);
    
    // If no exists any element to refine, the tolerance was reached
    if (!nelhrefs && !nelprefs)
        return true;
    return false;
}

// Writing a relation between number of degree of freedom and L2 error.
bool PrintResultsInMathematicaFormat(int typeel,int table,int ref,TPZVec<REAL> &ErrorVec,TPZVec<long> &NEquations,std::ostream &fileerrors) {
	long nref;
	long nelem = NEquations.NElements();
	if (!nelem) return false;
	int nerros = ErrorVec.NElements() / nelem;
//    STATE fact = 1.0e6;
	// setting format for ostream
	fileerrors << setprecision(20);
	fileerrors.setf(std::ios::fixed, std::ios::floatfield);
	fileerrors << "\n\n NEquations = {";

	// printing number of equations into a list
	nelem--;
	for(nref=0;nref<nelem;nref++) {
		fileerrors << NEquations[nref] << ", ";
	}
	fileerrors << NEquations[nref] << "};" << std::endl << "L2Error = {";
	// printing L2 error values into a list
	for(nref=0;nref<nelem;nref++) {
		fileerrors << ErrorVec[nref*nerros + 1] << ", ";
	}
	fileerrors << ErrorVec[nref*nerros+1] << "};" << std::endl << "SemiH1Error = {";
	// printing H1 error values into a list
	for (nref = 0; nref<nelem; nref++) {
		fileerrors << ErrorVec[nref*nerros + 2] << ", ";
	}
	fileerrors << ErrorVec[nref*nerros+2] << "};" << std::endl << "EnergyError = {";
	// printing Energy error values into a list
	for (nref = 0; nref<nelem; nref++) {
		fileerrors << ErrorVec[nref*nerros] << ", ";
	}
	fileerrors << ErrorVec[nref*nerros] << "};";
	// printing lines to create lists of logarithms
	fileerrors << std::endl << "LogNEquations = Table[Log[10,NEquations[[i]]],{i,1,Length[NEquations]}];" << std::endl;
	fileerrors << "LogL2Errors = Table[Log[10,L2Error[[i]]],{i,1,Length[L2Error]}];" << std::endl;
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "L2 = ";
	fileerrors << "ListPlot[{Table[{LogNEquations[[i]],LogL2Errors[[i]]},{i,1,Length[LogNEquations]}]";
	fileerrors << "},Joined->True,PlotRange->All,PlotStyle->";
    switch(table) {
        case 1:
            fileerrors << "Blue";
            break;
        case 2:
            fileerrors << "Gray";
            break;
        case 3:
            fileerrors << "Green";
            break;
        case 4:
            fileerrors << "Red";
            break;
        case 5:
            fileerrors << "Orange";
            break;
        default:
            fileerrors << "Yellow";
    }
    fileerrors << ",AspectRatio->1]" << std::endl;
	fileerrors << "LogSemiH1Errors = Table[Log[10,SemiH1Error[[i]]],{i,1,Length[SemiH1Error]}];";
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "SNH1 = ";
	fileerrors << "ListPlot[{Table[{LogNEquations[[i]],LogSemiH1Errors[[i]]},{i,1,Length[LogNEquations]}]";
	fileerrors << "},Joined->True,PlotRange->All,AspectRatio->1]\n" << std::endl;
	fileerrors << "LogEnergyErrors = Table[Log[10,EnergyError[[i]]],{i,1,Length[EnergyError]}];" << std::endl;
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "H1 = ";
	fileerrors << "ListPlot[{Table[{LogNEquations[[i]],LogEnergyErrors[[i]]},{i,1,Length[LogNEquations]}]";
	fileerrors << "},Joined->True,PlotRange->All,AspectRatio->1]\n" << std::endl;
    
    fileerrors << "Show[{";
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "L2,";
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "SNH1,";
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "H1}";
    fileerrors << ",PlotRange->All,AspectRatio->1]";

	return true;
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

int DefineDimensionOverElementType(int typeel) {
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
            mat->SetForcingFunction(new TPZDummyFunction<STATE>(RightTermArcTangentBad));
        }
            break;
        case 2:
            break;
        default:
            break;
    }
    cmesh->InsertMaterialObject(mat);
    // Make compatible dimension of the model and the computational mesh
    cmesh->SetDimModel(ModelDimension);
    
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

/***  Functions to manipulate and determine time (current) ***///
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

bool CreateCurrentResultDirectory(SimulationCase &sim_case) {
    // To compute processing times
    time_t sttime;
    struct tm* tmtime;
    time(&sttime);
    tmtime = localtime(&sttime);
    tmtime->tm_year += 1900;
    tmtime->tm_mon += 1;
    
    char command[512];
    memset(command,0,512);
    snprintf(command,512,"%s_%04d_%02d_%02d_%02d%02d%02d",sim_case.dir_name.c_str(),tmtime->tm_year,tmtime->tm_mon,tmtime->tm_mday,tmtime->tm_hour,tmtime->tm_min,tmtime->tm_sec);
    
	sim_case.dir_name = command;
	snprintf(command, 512, "mkdir %s", sim_case.dir_name.c_str());
    // Creating the directory
    return ((bool)system(command));
}

/* uncondense the elements unwrap the elements
void UnwrapMesh(TPZCompMesh *cmesh)
{
	long nel = cmesh->NElements();
	bool change = true;
	while (change)
	{
		change = false;
		for (long el = 0; el<nel; el++) {

			TPZCompEl *cel = cmesh->Element(el);
			TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
			if (condense) {
				condense->Unwrap();
				change = true;
			}
			cel = cmesh->Element(el);
			TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
			if (elgr) {
				elgr->Unwrap();
				change = true;
			}
		}
	}
}
*/

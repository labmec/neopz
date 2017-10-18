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
#include "pzcheckmesh.h"

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
#include "TPZSSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"

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
#include "TPZCreateHDivMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"

#include "CreateAndRefineMeshes.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;


/*** Data Types needed ***/
// Simulation Case
struct SimulationCase {
    bool  IsHdivQ;
    int   n_acc_terms;
    int   eltype;
    int   nthreads;
    std::string  dir_name;
    
    SimulationCase() : IsHdivQ(false), n_acc_terms(0), eltype(7), nthreads(0), dir_name("dump") {
        
    }
    SimulationCase(const SimulationCase &other) : IsHdivQ(other.IsHdivQ), n_acc_terms(other.n_acc_terms), eltype(other.eltype), nthreads(other.nthreads), dir_name(other.dir_name) {
        
    }
};


/**  Global variables  */
REAL GlobScale = 1.;
// Maximum number of equations allowed
long MaxEquations = 1500000;
// Input - output
ofstream out("OutPoissonArcTan.txt",ios::app);             // To store output of the console
// ABOUT H P ADAPTIVE
int MaxPOrder = 12;     // Maximum order for p refinement allowed
int MaxHLevel = 8;      // Maximum level for h refinement allowed
int MaxHUsed = 0;
int MaxPUsed = 0;

int ninitialrefs = 1;

// Poisson problem
STATE ValueK = 100000;
STATE F = sqrt(ValueK);
int ModelDimension;
// Circunference with high gradient - data
TPZManVector<REAL,3> CCircle(3,0.5);
REAL RCircle = 0.25;

//**********   Creating computational mesh with materials    *************
TPZCompMesh *CreateComputationalMesh(TPZGeoMesh *gmesh,int dim,int materialId,int hasforcingfunction,int id_bc0,int id_bc1=0,int id_bc2=0);


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

bool ApplyingHPAdaptiveStrategyBasedOnUAndDU(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol);

bool ApplyingStrategyHPAdaptiveBasedOnErrorOfSolution(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int ref,int itypeel,REAL &factorError);
bool ApplyingStrategyHPAdaptiveBasedOnErrorOfSolutionAndGradient(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int ref,int itypeel,TPZVec<REAL> &Tol);
void ApplyingStrategyPAdaptiveBasedOnExactSphereSolution(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,int ref);

// Writing a relation between number of degree of freedom and L2 error.
bool PrintResultsInMathematicaFormat(TPZVec<REAL> &ErrrVec,TPZVec<long> &NEquations,std::ostream &fileerrors);

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
	// 1 -> line		2 -> triangles		3 -> quadrilateral
    // 4 -> tetraedro	5 -> pyramid		6 -> prisma				7 -> cubo
    struct SimulationCase dummied;

	// Type of elements
	int itypeel = 3;
	// number of initial refinements over original mesh

	int count = 0;
	do {
		if(argc > 1)
			itypeel = atoi(argv[count+1]);
		if(itypeel > 7 || itypeel < 2)
			itypeel = 7;
		count++;
		dummied.eltype = itypeel;
		// Solving symmetricPoissonProblem on [0,1]^d with d=1, d=2 and d=3
		if(!SolveSymmetricPoissonProblemOnCubeMesh(dummied))
			return 1;

		itypeel++;
	} while(count < argc-1);
	
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
	TPZVec<REAL> Tol(3, 1.e-4);
//    Tol[1] = sqrt(Tol[0]); Tol[2] = sqrt(sqrt(Tol[1]));
    Tol[1] = sqrt(sqrt(sqrt(Tol[1])));
    if(Tol[1] < 1.) Tol[2] = 1.;
    else Tol[2] = 2*Tol[1];

	int materialId = 1;
	int id_bc0 = -1;
	int id_bc1 = -2;
	// Generic data for problems to solve
	int NRefs = 15;

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
	UniformRefinement(ninitialrefs, gmesh, ModelDimension);

	// Printing initial geometric mesh
    std::stringstream sout1;
    sout1 << sim_case.dir_name.c_str() << "/InitialGMesh_" << ModelDimension << "D_E" << sim_case.eltype << ".vtk";
	ofstream arg2(sout1.str().c_str());
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, arg2);

	// Printing solution on mesh with initial mesh before of adaptive process
//	TPZGeoMesh *gmeshfirst = CreateGeomMesh(typeel,materialId,id_bc0,id_bc1);
//	TPZGeoMesh gmeshfirst(*gmesh);
//	TPZCompEl::SetgOrder(1);
//	TPZCompMesh *cmeshfirst = CreateComputationalMesh(&gmeshfirst,ModelDimension,materialId,1,id_bc0,id_bc1);
//	TPZAnalysis an_sol(cmeshfirst,false);
//	LoadSolutionFirstOrder(cmeshfirst,ExactSolutionArcTangent);

    // Initializing the vectors of errors to store the errors for any iteration
    TPZMaterial *mat = new TPZMatPoisson3d();
    int nerros = mat->NEvalErrors();
    TPZVec<REAL> ErrorVecByIteration(nerros*NRefs,0.0);
    TPZVec<long> NEquations(NRefs,0L);
    TPZVec<STATE> ErrorU, ErrorDU;
    
    /** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("Pressure");

//    std::stringstream sout2;
//    sout2 << sim_case.dir_name.c_str() << "/Poisson" << ModelDimension << "D_MESHINIT_E" << sim_case.eltype << "WITHOUTREF" << ".vtk";
//	an_sol.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout2.str());

//	an_sol.PostProcess(3,ModelDimension);
//	long countels = 0;
//	for(int ii=0;ii<cmeshfirst->NElements();ii++) {
//		if(!cmeshfirst->ElementVec()[ii] || cmeshfirst->ElementVec()[ii]->Dimension()!=ModelDimension) continue;
//		countels++;
//	}
//	out << std::endl << "Number of elements 2D: " << countels << std::endl << std::endl ;
//	delete cmeshfirst;
//	gmeshfirst.CleanUp();
            
	// Creating computational mesh (approximation space and materials)
	int p = 1, pinit;
	MaxPUsed = pinit = p;
	MaxHUsed = 1;
	TPZCompEl::SetgOrder(p);
	TPZCompMesh *cmesh;
	gmesh->SetName("Malha Geometrica original");
    
    int hdivplusplus = sim_case.n_acc_terms;
    cmesh = CreateComputationalMesh(gmesh,ModelDimension,materialId,1,id_bc0,id_bc1);     // Forcing function is out 2013_07_25
    //if(HDiv) cmesh = CreateHDivMesh(gmesh, meshvec, p, ModelDimension,hdivplusplus);
    
	// To storing number of equations and errors obtained for all iterations
    AdjustFluxPolynomialOrders(cmesh, hdivplusplus);
    //if(HDiv) ReconstructHDivMesh(cmesh, meshvec, hdivplusplus);

	int nref = 0;
    bool tolachieved = false;

	// loop solving iteratively
	do {
		out << "\nSolving Poisson problem " << ModelDimension << "D. Refinement: " << nref << " TypeElement: " << sim_case.eltype << endl;
		std::cout << "\nSolving Poisson problem. Type of element: " << sim_case.eltype << std::endl;
				
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
        
        UnwrapMesh(cmesh);
		
        //if(HDiv) TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmesh);
		// Post processing
        std::stringstream sout3;
        sout3 << sim_case.dir_name.c_str() << "/" << "Poisson" << ModelDimension << "D_E" << sim_case.eltype << "H" << std::setprecision(2) << nref << "P" << pinit << ".vtk";
        an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout3.str().c_str());
		an.PostProcess(1,ModelDimension);
		std::ofstream out(sout.str().c_str());
		cmesh->LoadReferences();
		TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(),out,false);
        
		// generation mesh process finished
		time(&endtime);
		time_elapsed = endtime - sttime;
		formatTimeInSec(ptime,256,time_elapsed);
		out << "  Time elapsed " << time_elapsed << " <-> " << ptime << "\n\n";
		fileerrors << "  Time elapsed " << time_elapsed << " <-> " << ptime << "\n";
		std::cout << "  Time elapsed " << time_elapsed << " <-> " << ptime << "\n\n";

		if(!ProcessingErrorUAndDUKnowingExactSol(an,ErrorVecByIteration,nref,ErrorU,ErrorDU))
            DebugStop();

		NEquations[nref] = cmesh->NEquations();

		std::cout << "\n NRef " << nref << "\tL2 Error " << ErrorVecByIteration[3*nref] << "  NEquations: " << NEquations[nref] << " PUsed " << MaxPUsed << " HMax " << MaxHUsed << std::endl << std::endl;
		out << "\n NRef " << nref << "\tL2 Error " << ErrorVecByIteration[3*nref] << "  NEquations: " << NEquations[nref] << " PUsed " << MaxPUsed << " HMax " << MaxHUsed << std::endl << std::endl;
		if(cmesh->NEquations() > MaxEquations) {
			NRefs = nref+1;							// final iteration
			ErrorVecByIteration.Resize(nerros*NRefs);
			NEquations.Resize(NRefs);
			continue;
		}
		fileerrors.flush();
		out.flush();

		// HP Refinement Process
		if(NRefs > 1 && nref < (NRefs-1)) {
			out << "\n\nApplying Adaptive Methods... step " << nref << "\n";
			std::cout << "\n\nApplying Adaptive Methods... step " << nref << "\n";
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU(an.Mesh(),ErrorU,ErrorDU,Tol);
		}
        std::cout << "NElements " << cmesh->NElements() << " NEquations " << cmesh->NEquations() << std::endl;

        // refazer a malha multifisica
        AdjustFluxPolynomialOrders(cmesh, hdivplusplus);
        //if(HDiv) ReconstructHDivMesh(cmesh, meshvec, hdivplusplus);

		fileerrors.flush();
		out.flush();
		// Sometimes Writing a relation between number of degree of freedom and L2 error.
        fileerrors << "H1 approximation\n";
        fileerrors << "H1plusplus = " << hdivplusplus << std::endl;
        //if(HDiv) {
        //  fileerrors << "HDiv approximation\n";
        //  fileerrors << "HDivplusplus = " << hdivplusplus << std::endl;
        //}
        PrintResultsInMathematicaFormat(ErrorVecByIteration,NEquations,fileerrors);
        fileerrors.flush();
        fileerrors << "done\n";

		nref++;
	}while (nref < NRefs && !tolachieved);

	if (cmesh) {
		cmesh->SetReference(0);
		delete cmesh;
	}
	cmesh = NULL;
	if(gmesh)
		delete gmesh;
	gmesh = NULL;

	// Writing a relation between number of degree of freedom and L2 error.
	sout.clear();
	sout << sim_case.dir_name.c_str() << "/ErrorsHP_Poisson.nb";
	std::ofstream finalerrors(sout.str().c_str());   // To store all errors calculated by TPZAnalysis (PosProcess)
	if(!PrintResultsInMathematicaFormat(ErrorVecByIteration,NEquations,finalerrors))
		std::cout << "\nThe errors and nequations values in Mathematica format was not done.\n";
	
	fileerrors << std::endl << "Finished running for element " << sim_case.eltype << std::endl << std::endl;
	fileerrors.close();
	std::cout << std::endl << "\tFinished running for element " << sim_case.eltype << std::endl << std::endl;
	out.close();
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

	long neq = cmesh->NEquations();
	if (ModelDimension != cmesh->Dimension())
		DebugStop();
	//    TPZVec<REAL> ux(neq);
	//    TPZVec<REAL> sigx(neq);
	TPZManVector<REAL, 10> totalerror(10, 0.);

	TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
	TPZManVector<REAL, 10> errors(10);
	errors.Fill(0.0);
	long i, nel = elvec.NElements();
	ErrorU.Resize(nel, 0.0);
	// The last position will be store the maxime value of the gradient errors
	ErrorDU.Resize(nel, 0.0);

	/** Computing error for all elements with same dimension of the model */
	for (i = 0L; i<nel; i++) {
		TPZCompEl *el = (TPZCompEl *)elvec[i];
		if (!el || el->Dimension() != ModelDimension) continue;
		errors.Fill(0.0);
		el->EvaluateError(analysis.fExact, errors, 0);
		int nerrors = errors.NElements();
		totalerror.resize(nerrors);
		for (int ier = 0; ier < nerrors; ier++)
			totalerror[ier] += errors[ier] * errors[ier];

		// L2 error for each element
		ErrorU[i] = sqrt(errors[1] * errors[1]);
		ErrorDU[i] = sqrt(errors[2] * errors[2]);
	}

	int nerrors = errors.NElements();
	if (ErrorVecByIteration.size()<nerrors*(nref + 1))
		return false;

	// Returns the square of the calculated errors.
	for (i = 0; i<nerrors; i++)
		ErrorVecByIteration[nerrors*nref + i] = sqrt(totalerror[i]);
	return true;
}

bool ApplyingHPAdaptiveStrategyBasedOnUAndDU(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol) {
//ApplyingStrategyHPAdaptiveBasedOnErrorOfSolutionAndGradient(TPZCompMesh *cmesh,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL MaxErrorByElement,REAL &MinErrorByElement,REAL &MinGrad,int nref,int itypeel,TPZVec<REAL> &Tol ) {
	if(!cmesh) return false;
    long iel, nelhrefs = 0, nelprefs = 0;
	long nels = cmesh->NElements();
    
	TPZManVector<long,27> subels;
	TPZManVector<long,27> subsubels;
    TPZVec<long> HRef(nels,0L), PRef(nels,0L);
	TPZCompEl *cel;

    // To know laplacian values as auxiliar information to adaptive
//	REAL LaplacianVal;
//	REAL MaxLaplacianVal, MinLaplacianVal;
    
//	ComputingMaxLaplacian(cmesh,MaxLaplacianVal,MinLaplacianVal);

	// Applying hp refinement only for elements with dimension as model dimension
    std::cout << "Refinando malha com " << nels  << " elementos, e" << cmesh->NEquations() << " equacoes.\n";

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
    
    // Doing P Refinement
    int pelement = 0;
    TPZGeoEl *gel = 0;
    TPZInterpolatedElement *intel = 0;
    for(iel=0;iel<nelprefs;iel++) {
        intel = dynamic_cast<TPZInterpolatedElement* > (cmesh->Element(PRef[iel]));
        if(!intel || intel->Dimension() != cmesh->Dimension()) DebugStop();
        gel = intel->Reference();
        if(!gel) DebugStop();
        pelement = intel->PreferredSideOrder(gel->NSides() - 1);
        if(pelement < MaxPOrder)
            intel->PRefine(pelement+1);
    }
    
    // Doing H Refinement
    for(iel=0;iel<nelhrefs;iel++) {
        bool twice = false;
        if(HRef[iel]<0) {
            twice = true;
            HRef[iel] *= -1;
        }
		subels.Resize(0);
        intel = dynamic_cast<TPZInterpolatedElement* > (cmesh->Element(HRef[iel]));
		if (!intel || intel->Dimension() != cmesh->Dimension()) DebugStop();
        gel = intel->Reference();
        if(!gel) DebugStop();
        if(gel->Level() < MaxHLevel) {
			cel = cmesh->Element(intel->Index());
			intel->Divide(iel,subels);
			delete cel;
        }
        if(twice) {
            for(long isub_el=0;isub_el<subels.NElements();isub_el++) {
                subsubels.Resize(0);
                TPZCompEl * isub_cel = cmesh->ElementVec()[subels[isub_el]];
                TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement* >(isub_cel);
				if (!intel1 || intel1->Dimension() != cmesh->Dimension()) DebugStop();
				cel = cmesh->Element(intel1->Index());
				intel1->Divide(subels[isub_el],subsubels);
				//cel = cmesh->Element(intel1->Index());
				delete intel;
			}
		}
    }
    
	cmesh->ExpandSolution();
    RegularizeMesh(cmesh->Reference(),cmesh->Dimension());
//	cmesh->AutoBuild();
    cmesh->ExpandSolution();
//	cmesh->AdjustBoundaryElements();
//	cmesh->CleanUpUnconnectedNodes();
	// Printing information stored
	PrintNRefinementsByType(nels,cmesh->NElements(),nelhrefs,nelprefs,out);
	PrintNRefinementsByType(nels,cmesh->NElements(),nelhrefs,nelprefs);
    // Cleaning
    ErrorU.Resize(0);
    ErrorDU.Resize(0);
    if(!nelhrefs && !nelprefs)
        return true;
	return false;
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

	TPZVec<STATE> Laplacian(1);
	TPZFMatrix<STATE> dLap(3);
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
		RightTermArcTangentBad(center,Laplacian,dLap);
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

// Writing a relation between number of degree of freedom and L2 error.
bool PrintResultsInMathematicaFormat(TPZVec<REAL> &ErrorVec,TPZVec<long> &NEquations,std::ostream &fileerrors) {
	int nref;
    STATE fact = 1.0e6;
	// setting format for ostream
	fileerrors << setprecision(20);
	fileerrors.setf(std::ios::fixed, std::ios::floatfield);
	fileerrors << "\n\n NEquations = {";

	// printing number of equations into a list
	for(nref=0;nref<NEquations.NElements()-1;nref++) {
		fileerrors << NEquations[nref] << ", ";
	}
	fileerrors << NEquations[nref] << "};" << std::endl << "L2Error = {";
	// printing error values into a list
	for(nref=0;nref<ErrorVec.NElements()-1;nref++) {
		fileerrors << ErrorVec[nref]*fact << ", ";
	}
	fileerrors << ErrorVec[nref] << "}/1000000.0;";
	// printing lines to create lists of logarithms
	fileerrors << std::endl << "LogNEquations = Table[Log[10,NEquations[[i]]],{i,1,Length[NEquations]}];" << std::endl;
	fileerrors << "LogL2Errors = Table[Log[10,L2Error[[i]]],{i,1,Length[L2Error]}];" << std::endl;
	fileerrors << "ListPlot[{Table[{LogNEquations[[i]],LogL2Errors[[i]]},{i,1,Length[LogNEquations]}]";
	fileerrors << "},Joined->True,PlotRange->All]\n" << std::endl;
	return true;
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
    return system(command);
}


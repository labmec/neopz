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
#include <math.h>
#include <fstream>
#include <cmath>

#include "problem.h"

#include "../LibRefine/CreateAndRefineMeshes.h"
#include "../LibRefine/HPAdaptiveProcesses.h"

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
int64_t MaxEquations = 300000;
// Input - output
ofstream out("OutPoissonArcTan.txt",ios::app);             // To store output of the console
// ABOUT H P ADAPTIVE
int MaxPOrder = 12;     // Maximum order for p refinement allowed
int MaxHLevel = 9;      // Maximum level for h refinement allowed
int MaxHUsed = 0;
int MaxPUsed = 0;

int ninitialrefs = 2;

// Poisson problem
STATE ValueK = 100000;
STATE F = sqrt(ValueK);
int ModelDimension;

// Circunference with high gradient - data
TPZManVector<REAL,3> CCircle(3,0.5);
REAL RCircle = 0.25;

// To run one time
bool Once = true;

//**********   Creating computational mesh with materials    *************
TPZCompMesh *CreateComputationalMesh(TPZGeoMesh *gmesh,int dim,int materialId,int hasforcingfunction,int id_bc0,int id_bc1=0,int id_bc2=0);

int DefineDimensionOverElementType(int typeel);
void GetFilenameFromGID(MElementType typeel, std::string &name);

/** PROBLEMS */
bool SolveSymmetricPoissonProblemOnCubeMesh(SimulationCase &sim_case);

//REAL ProcessingError(TPZAnalysis &analysis,TPZVec<REAL> &ervec,TPZVec<REAL> &ervecbyel,TPZVec<REAL> &gradervecbyel,REAL &MinErrorByElement,REAL &);
void LoadSolutionFirstOrder(TPZCompMesh *cmesh, void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv,TPZVec<STATE> &ddsol));

// Writing a relation between number of degree of freedom and L2 error.
bool PrintResultsInMathematicaFormat(int typeel,int hpcase,int nref,TPZVec<REAL> &ErrrVec,TPZVec<int64_t> &NEquations,std::ostream &fileerrors);

void AdjustingOrder(TPZCompMesh *cmesh);
int MaxLevelReached(TPZCompMesh *cmesh);


/** Utilitaries Over Date And Time */
void formatTimeInSec(char *strtime,int lenstrtime,int timeinsec);
bool CreateCurrentResultDirectory(SimulationCase &sim);

#ifdef LOG4CXX
static LoggerPtr  logger(Logger::getLogger("pz.refinehp"));
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
        dummied.hpcase = 2;
        // loop over use of specific strategy hp-adaptive table
        do {
            time_t totalsttime, totalendtime;
            int time_elapsed;
            char timeformat[1024];
			memset(timeformat, 1024, 0);
			char *ptime;
            // Initial message to print computed errors
            time(&totalsttime);
            ptime = ctime(&totalsttime);
            out << "\nApproximation Error in " << ptime << std::endl << "\nType of element: " << dummied.eltype << endl;
            std::cout << "\nApproximation Error in " << ptime << std::endl << "\nType of element: " << dummied.eltype << endl;

            dummied.eltype = itypeel;
            dummied.SetDirName();
            // Solving symmetricPoissonProblem on [0,1]^d with d=1, d=2 and d=3
            if(!SolveSymmetricPoissonProblemOnCubeMesh(dummied))
                return 1;

            // generation mesh process finished
            time(&totalendtime);
            time_elapsed = totalendtime - totalsttime;
            formatTimeInSec(timeformat,1024,time_elapsed);
            out << " Time elapsed " << time_elapsed << " <-> " << timeformat << " Total solving process " << dummied.hpcase << ".\n";
            std::cout << " Time elapsed " << time_elapsed << " <-> " << timeformat << " Total solving process " << dummied.hpcase << ".\n";
            dummied.hpcase++;

        }while(dummied.hpcase < 4);

        itypeel++;
	} while(itypeel < 8 && !Once);
    out.close();
	return 0;
}

bool SolveSymmetricPoissonProblemOnCubeMesh(SimulationCase &sim_case) {
    if(CreateCurrentResultDirectory(sim_case))
        return false;
    
    // To compute processing times
    time_t sttime, endtime;
    int time_elapsed;
    char timeformat[1024];   // = time_formated;
	memset(timeformat, 1024, 0);

	// Tolerance for applying hp adaptivity
	TPZManVector<REAL,3> Tol(3, 1.e-8);
    Tol[1] = 1.e-7; Tol[2] = 1.e-5;

	int materialId = 1;
	int id_bc0 = -1;
	int id_bc1 = -2;
    
	// Generic data for problems to solve
	int NRefs = 7;

	// Output files
	std::stringstream sout;
	sout << sim_case.dir_name.c_str() << "/ErrorsHP_Poisson.txt";
	std::ofstream fileerrors(sout.str().c_str());   // To store all errors calculated by TPZAnalysis (PosProcess)
    // Initializing the generation mesh process
    time(&sttime);

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
	TPZVec<int64_t> NEquations(NRefs, 0L);
	TPZVec<STATE> ErrorU, ErrorDU;
	
    // generation mesh process finished
    time(&endtime);
    time_elapsed = endtime - sttime;
    formatTimeInSec(timeformat,1024,time_elapsed);
    out << " Time elapsed " << time_elapsed << " <-> " << timeformat << " Generating mesh.\n";
    fileerrors << " Time elapsed " << time_elapsed << " <-> " << timeformat << " Generating mesh.\n";
    std::cout << " Time elapsed " << time_elapsed << " <-> " << timeformat << " Generating mesh.\n";

	int nref = 0;
    bool tolachieved = false;

	// loop solving iteratively
	do {
        out << "\n\nSOLVING POISSON PROBLEM " << ModelDimension << "D." << " ELEMENT Type: " << sim_case.eltype << " STRATEGY: " << sim_case.hpcase << " Iteration: " << nref << std::endl;
        std::cout << "\n\nSOLVING POISSON PROBLEM " << ModelDimension << "D." << " ELEMENT Type: " << sim_case.eltype << " STRATEGY: " << sim_case.hpcase << " Iteration: " << nref << std::endl;

		// Initializing the resolution process
		time(& sttime);
				
		// Introduzing exact solution depending on the case
		// Solving adaptive process
        // To clean unwhished information
//        if(!tolachieved) {
//            TPZCompMesh *newmesh = cmesh->Clone();
//            delete cmesh;
//            cmesh = newmesh;
//            cmesh->AdjustBoundaryElements();
//            cmesh->InitializeBlock();
//            cmesh->ExpandSolution();
//        }
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
        
		// resolution process finished
		time(&endtime);
		time_elapsed = endtime - sttime;
		formatTimeInSec(timeformat,1024,time_elapsed);
		out << " Time elapsed " << time_elapsed << " <-> " << timeformat << " Resolution process.\n";
		fileerrors << " Time elapsed " << time_elapsed << " <-> " << timeformat << " Resolution process.\n";
		std::cout << " Time elapsed " << time_elapsed << " <-> " << timeformat << " Resolution process.\n";

        // Initializing the hp-adaptive process
        time(&sttime);

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
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnU_I(an.Mesh(),ErrorU,ErrorDU,Tol,MaxPOrder,MaxHLevel,out);
        else if(sim_case.hpcase == 2)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnU_II(an.Mesh(),ErrorU,ErrorDU,Tol, MaxPOrder, MaxHLevel,out);
        else if(sim_case.hpcase == 3)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnU_III(an.Mesh(),ErrorU,ErrorDU,Tol, MaxPOrder, MaxHLevel,out);
        else if(sim_case.hpcase == 4)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDUAsArticle_II(an.Mesh(),ErrorU,ErrorDU,Tol, MaxPOrder, MaxHLevel,out);
        else if(sim_case.hpcase == 5)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_III(an.Mesh(),ErrorU,ErrorDU,Tol, MaxPOrder, MaxHLevel,out);
        else if(sim_case.hpcase == 6)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_IV(an.Mesh(),ErrorU,ErrorDU,Tol, MaxPOrder, MaxHLevel,out);
        else if(sim_case.hpcase == 7)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_V(an.Mesh(),ErrorU,ErrorDU,Tol, MaxPOrder, MaxHLevel,out);
        else if(sim_case.hpcase == 8)
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_VI(an.Mesh(),ErrorU,ErrorDU,Tol, MaxPOrder, MaxHLevel,out);
        else
            tolachieved = ApplyingHPAdaptiveStrategyBasedOnUAndDU_XI(an.Mesh(),ErrorU,ErrorDU,Tol, MaxPOrder, MaxHLevel,out);

        out << "\n Applying Adaptive Methods... step " << nref << "\n";
		std::cout << "\n Applying Adaptive Methods... step " << nref << "\n";
		std::cout << " NElements Comp " << cmesh->NElements() << " NEquations " << cmesh->NEquations() << std::endl;
		out << " NElements Comp " << cmesh->NElements() << " NEquations " << cmesh->NEquations() << std::endl;

		// Sometimes Writing a relation between number of degree of freedom and L2 error.
        fileerrors << "H1 approximation\n";

        // generation mesh process finished
        time(&endtime);
        time_elapsed = endtime - sttime;
        formatTimeInSec(timeformat,1024,time_elapsed);
        out << " Time elapsed " << time_elapsed << " <-> " << timeformat << " HP-adaptive process.\n";
        fileerrors << " Time elapsed " << time_elapsed << " <-> " << timeformat << " HP-adaptive process.\n";
        std::cout << " Time elapsed " << time_elapsed << " <-> " << timeformat << " HP-adaptive process.\n";

        // Printing partial table in Mathematica format
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

	// Printing in Mathematica notekook the relation between number of degree of freedom and L2 error and others errors.
    std::stringstream sout4;
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

// Writing a relation between number of degree of freedom and L2 error.
bool PrintResultsInMathematicaFormat(int typeel,int table,int ref,TPZVec<REAL> &ErrorVec,TPZVec<int64_t> &NEquations,std::ostream &fileerrors) {
	int64_t nref;
	int64_t nelem = NEquations.NElements();
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
    fileerrors << ErrorVec[nref*nerros+1] << "};" << std::endl;
/*    fileerrors << "SemiH1Error = {";
	// printing H1 error values into a list
	for (nref = 0; nref<nelem; nref++) {
		fileerrors << ErrorVec[nref*nerros + 2] << ", ";
	}
    fileerrors << ErrorVec[nref*nerros+2] << "};" << std::endl;
    fileerrors << "EnergyError = {";
	// printing Energy error values into a list
	for (nref = 0; nref<nelem; nref++) {
		fileerrors << ErrorVec[nref*nerros] << ", ";
	}
	fileerrors << ErrorVec[nref*nerros] << "};"; */
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
    fileerrors << ",AspectRatio->2.05]" << std::endl;
    /*
    fileerrors << "LogSemiH1Errors = Table[Log[10,SemiH1Error[[i]]],{i,1,Length[SemiH1Error]}];" << std::endl;
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "SNH1 = ";
	fileerrors << "ListPlot[{Table[{LogNEquations[[i]],LogSemiH1Errors[[i]]},{i,1,Length[LogNEquations]}]";
	fileerrors << "},Joined->True,PlotRange->All,AspectRatio->2.05]\n" << std::endl;
	fileerrors << "LogEnergyErrors = Table[Log[10,EnergyError[[i]]],{i,1,Length[EnergyError]}];" << std::endl;
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "H1 = ";
	fileerrors << "ListPlot[{Table[{LogNEquations[[i]],LogEnergyErrors[[i]]},{i,1,Length[LogNEquations]}]";
	fileerrors << "},Joined->True,PlotRange->All,AspectRatio->2.05]\n" << std::endl; */
    
    fileerrors << "Show[{";
    fileerrors << "E" << typeel << "Table" << table << "Ref" << ref << "L2";
//    fileerrors << ",E" << typeel << "Table" << table << "Ref" << ref << "SNH1";
//    fileerrors << ",E" << typeel << "Table" << table << "Ref" << ref << "H1";
    fileerrors << "},PlotRange->All,AspectRatio->2.05]";

	return true;
}

void AdjustingOrder(TPZCompMesh *cmesh) {
	if(!cmesh) return;
	int64_t nels = cmesh->NElements();
	TPZInterpolatedElement *el;
	STATE Tol;
	ZeroTolerance(Tol);

	// To see where the computation is doing
	int64_t i;
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
	int64_t nels = cmesh->NElements();
	TPZVec<int64_t> subels;
	int j, k, pelement, dp = 1;
	TPZVec<int64_t> subsubels;
	TPZInterpolatedElement *el;
	REAL errorcel = 0.0;
	for(int64_t i=0L;i<nels;i++) {
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
            mat->SetForcingFunction(new TPZDummyFunction<STATE>(RightTermArcTangentBad, 5));
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
            bc1->SetForcingFunction(new TPZDummyFunction<STATE>(ExactSolLaplaceBC, 5));
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


/**
 * @file Poisson 3D on Fichera corner: First problem as Ander Solin presentation. Second problem as Rachowicz article 2006
 */

#include "pzshapelinear.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

#include "pzbstrmatrix.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzpoisson3d.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "TPZRefPatternTools.h"

#include "TPZParSkylineStructMatrix.h"

#include <stdio.h>

#include "pzlog.h"

#include "pzgeoelbc.h"
#include "TPZGeoCube.h"

#include "tpzgeoelrefpattern.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;

/**
 * @addtogroup Tutorials
 * @{
 */

// Global variable
int gLMax;
int NUniformRefs = 2;
REAL alfa = M_PI/6.;
bool anothertests = false;
int nstate = 2;

char saida[512];
ofstream out("ConsolePoisson3D.txt");   // output file from console  -> Because it has many energy faults

/** Printing level */
int gPrintLevel = 0;
bool gDebug = false;

void ExactRachowicz(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);

void ExactSolin(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
void BCSolin(const TPZVec<REAL> &x, TPZVec<STATE> &sol);
void FforcingSolin(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &dsol);

//void InitializeSolver(TPZAnalysis &an);
void InitialSolutionLinearConvection(TPZFMatrix<STATE> &InitialSol, TPZCompMesh *cmesh);
void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,REAL radius,int ntyperefs);
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &points,REAL r,REAL &distance,bool &isdefined);

TPZGeoMesh *ConstructingFicheraCorner(REAL L,int typeel=0,int problem = 1);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction,int problem=1);

void formatTimeInSec(char *strtime,int timeinsec);

// problem = 1 represents the Poisson as Solin presentation
// problem = 2 represents the Poisson as Rachowicz Ficher corner problem
int problem = 1;

// MAIN FUNCTION TO NUMERICAL SOLVE OF THE FICHERA CORNER PROBLEM
int main(int argc, char *argv[]) {
	
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
	char tempo[256];
    
	// To compute the errors
	ofstream fileerrors("ErrorsHP3D_Fichera.txt");
    
	TPZVec<REAL> ervec(100,0.0);
	// Printing computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	REAL InitialL = 1.0;
	int nref, NRefs = 5;
	int nthread, NThreads = 2;
	int dim = 3;
	// Problem
	out << "\nFICHERA CORNER:";
	if(problem==1)
		out << "\t As SOLIN Presentation.";
	else if(problem==2) 
		out << "\t As Rachowicz (2006).";
	
	for(int ntyperefs=1;ntyperefs<3;ntyperefs++) {
		fileerrors << "Type of refinement: " << ntyperefs << " Level. " << endl;
        
        // typeel = 0 -> hexahedra
        // typeel = 1 -> prisms
        // typeel = 2 -> pyramids
		for(int typeel=2;typeel>=0;typeel--) {
			REAL radius = 1.2;
			// Constructing geometric mesh as Fichera corner using hexahedra
			TPZGeoMesh *gmesh3D = ConstructingFicheraCorner(InitialL,typeel,problem);
			
			fileerrors << "Type of element: " << typeel << " (0-hexahedra, 1-three prisms, 2-four pyramids." << endl;
			for(nref=0;nref<NRefs;nref++) {
				if(nref > 4) nthread = 2*NThreads;
				else nthread = NThreads;
				
				// Initializing the generation mesh process
				time(&sttime);
				cout << "\nConstructing Fichera problem. Refinement: " << nref+1 << " Threads: " << nthread << " TypeRef: " << ntyperefs << " TypeElement: " << typeel << endl;
				
				// h_refinement
				// Refining elemenst whose center belong a circunference with radio r +- maxime distance radius
				RefiningNearCircunference(dim,gmesh3D,radius,ntyperefs);
				if(ntyperefs==2) {
					nref++;
					radius *= 0.4;
				}
				else
					radius *= 0.6;
				//			if(nref == NRefs-1) {
				//				sprintf(saida,"gmesh_3DFichera_H%dTR%dE%d.vtk",nref,ntyperefs,typeel);
				//				PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
				//			}
				
				// Creating computational mesh
				// Selecting orthogonal polynomial family to construct shape functions
				if(anothertests)
					TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;  // Setting Chebyshev polynomials as orthogonal sequence generating shape functions
				/** Set polynomial order */
				int p = 2, pinit;
				pinit = p;
				TPZCompEl::SetgOrder(1);
				TPZCompMesh *cmesh = CreateMesh(gmesh3D,dim,1,problem);
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
					p = (highlevel-level)+1;
					if(!p) p = 1;     // Fazendo os dois maiores niveis de refinamento devem ter ordem 1
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
				
				// SOLVING PROCESS
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
				
				out << "\tRefinement: " << nref+1 << " TypeRef: " << ntyperefs << " TypeElement: " << typeel << " Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << tempo << "\n\n\n";
				
				// Post processing
				char pp[300];
				std::string filename = "Poisson3DSol_";
				if(problem==1) filename += "S_";
				else if(problem==2) filename += "R_";
				sprintf(pp,"TR%1dE%1dT%02dH%02dP%02d",ntyperefs,typeel,nthread,(nref+1),pinit);
				filename += pp;
				filename += ".vtk";
				
				/** Variable names for post processing */
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
				
				// Computing error
				if(problem==1)
					an.SetExact(ExactSolin);
				else if(problem==2)
					an.SetExact(ExactRachowicz);
				fileerrors << "Refinement: " << nref+1 << "  Dimension: " << dim << "  NEquations: " << cmesh->NEquations();
                bool store_errors = false;
				an.PostProcessError(ervec,store_errors, out);
				for(int rr=0;rr<ervec.NElements();rr++)
					fileerrors << "  Error_" << rr+1 << ": " << ervec[rr]; 
				fileerrors << "  TimeElapsed: " << time_elapsed << " <-> " << tempo << std::endl;
				
				delete cmesh;
			}
			delete gmesh3D;
		}
	}
	out.close();
	fileerrors.close();
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////   FICHERA CORNER - Problem as Anders Solin Presentation   ///////////////////
////////////////////////////////////////////////////////////////////////////////////////

void ExactSolin(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r);
	sol[0] = sqrt(raiz);
	REAL den = sol[0]*sol[0]*sol[0];
	if(!IsZero(den)) {
		dsol(0,0) = .5*x[0]/den;
		dsol(1,0) = .5*x[1]/den;
		dsol(2,0) = .5*x[2]/den;
	}
	else {
		DebugStop();
	}
}

void BCSolin(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r);
	bcsol[0] = sqrt(raiz);
}

void FforcingSolin(const TPZVec<REAL> &x, TPZVec<STATE> &f,TPZFMatrix<STATE> &df) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt( sqrt(quad_r * quad_r * quad_r) );
	if(!IsZero(raiz)) {
		f[0] = 3./(4.0*raiz);
	} else {
		DebugStop();
	}
}

void ExactRachowicz(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r);
	sol[0] = sqrt(raiz);
	REAL den = sol[0]*sol[0]*sol[0];
	if(!IsZero(den)) {
		dsol(0,0) = .5*x[0]/den;
		dsol(1,0) = .5*x[1]/den;
		dsol(2,0) = .5*x[2]/den;
	}
	else {
		DebugStop();
	}
}

void BCRachowiczN(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r);
	bcsol[0] = sqrt(raiz);
}
void BCRachowiczD(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	bcsol[0] = 0.0;
}

void FforcingRachowicz(const TPZVec<REAL> &x, TPZVec<STATE> &f) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r * quad_r * quad_r);
	f[0] = 3./(4.0*sqrt(raiz));
}

TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction,int problem) {
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
	// Creating Poisson material
	TPZMaterial *mat = new TPZMatPoisson3d(1,dim);
	TPZVec<REAL> convd(3,0.);
	((TPZMatPoisson3d *)mat)->SetParameters(1.,0.,convd);
	if(hasforcingfunction) {
		switch(problem) {
			case 1:
				mat->SetForcingFunction(new TPZDummyFunction<STATE>(FforcingSolin, 5));
				break;
			case 2:
				mat->SetForcingFunction(new TPZDummyFunction<STATE>(FforcingRachowicz, 5));
				break;
		}
	}
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	
	// Boundary conditions
	TPZFMatrix<STATE> val1(dim,dim,0.),val2(dim,1,0.);
	switch(problem) {
		case 1:
		{
			// Dirichlet 
			TPZAutoPointer<TPZFunction<STATE> > FunctionBC = new TPZDummyFunction<STATE>(BCSolin, 5);
			val1.PutVal(0,0,1.);
			val1.PutVal(1,1,1.);
			val1.PutVal(2,2,1.);
			TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
			bnd->SetForcingFunction(FunctionBC);
			cmesh->InsertMaterialObject(bnd);
		}
		break;
		case 2:
		{
			// Neumann
			TPZAutoPointer<TPZFunction<STATE> > FunctionBCN = new TPZDummyFunction<STATE>(BCRachowiczN, 5);
			val1.PutVal(0,0,1.);
			val1.PutVal(1,1,1.);
			val1.PutVal(2,2,1.);
			TPZMaterial *bndN = mat->CreateBC(mat,-1,1,val1,val2);
			bndN->SetForcingFunction(FunctionBCN);
			cmesh->InsertMaterialObject(bndN);
			// Dirichlet
//			TPZAutoPointer<TPZFunction<STATE> > FunctionBCD = new TPZDummyFunction<STATE>(BCRachowiczD, 5);
			TPZMaterial *bndD = mat->CreateBC(mat,-2,0,val1,val2);
//			bndD->SetForcingFunction(FunctionBCD);
			cmesh->InsertMaterialObject(bndD);
		}
		break;
	}			
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->ExpandSolution();
	return cmesh;
}

// CONSTRUCTION OF THE FICHERA CORNER FROM A CUBE - SUBDIVIDING ITS - AND DELETING A RIGHT TOP CUBE SON
TPZGeoMesh *ConstructingFicheraCorner(REAL InitialL, int typeel,int problem) {
	
	// CREATING A CUBE WITH MASS CENTER THE ORIGIN AND VOLUME = INITIALL*INITIALL*INITIALL 
	REAL co[8][3] = {
		{-InitialL,-InitialL,-InitialL},
		{InitialL,-InitialL,-InitialL},
		{InitialL,InitialL,-InitialL},
		{-InitialL,InitialL,-InitialL},
		{-InitialL,-InitialL,InitialL},
		{InitialL,-InitialL,InitialL},
		{InitialL,InitialL,InitialL},
		{-InitialL,InitialL,InitialL}
	};
	int64_t indices[1][8] = {{0,1,2,3,4,5,6,7}};
	
	const int nelem = 1;
	int nnode = 8;
	
	TPZGeoEl *elvec[nelem];
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	int nod;
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
		TPZManVector<int64_t> nodind(8);
		for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
		int64_t index;
		elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
	}
	gmesh->BuildConnectivity();
	
	// SUBDIVIDING A CUBE
	TPZVec<TPZGeoEl*> sub;
	
	switch(typeel) {
		case 0:
			gmesh->ElementVec()[0]->Divide(sub);
			// DELETING A CUBE 6th
			delete gmesh->ElementVec()[7];
			break;
		case 1:
		{
			gmesh->ElementVec()[0]->Divide(sub);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
			// Dividing hexahedron in three prisms
			// First the hexahedra with z < 0
			std::string filename = REFPATTERNDIR;
			filename += "/3D_Hexa_Rib_Side_16_18.rpt";
			
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
			if(!gRefDBase.FindRefPattern(refpat))
			{
				gRefDBase.InsertRefPattern(refpat);
			}
			TPZGeoElRefPattern <TPZGeoCube> *gelrp;
			TPZGeoEl *gel;
			for(int i=1;i<5;i++) {
				if(i==7) i++;
				gel = gmesh->ElementVec()[i];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);
			}
			// Second the hexahedra with z > 0
			filename = REFPATTERNDIR;
			filename += "/3D_Hexa_Rib_Side_8_10.rpt";
			
			TPZAutoPointer<TPZRefPattern> refpat2 = new TPZRefPattern(filename);
			if(!gRefDBase.FindRefPattern(refpat2))
			{
				gRefDBase.InsertRefPattern(refpat2);
			}
			for(int i=5;i<9;i++) {
				if(i==7) i++;
				gel = gmesh->ElementVec()[i];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat2);
				gel->Divide(sub);
			}
			// DELETING A CUBE 6th
			delete gmesh->ElementVec()[7];
		}
			break;
		case 2:
			{
				gmesh->ElementVec()[0]->Divide(sub);
				gmesh->ResetConnectivities();
				gmesh->BuildConnectivity();
				// Dividing hexahedron in four pyramids
				TPZGeoElRefPattern <TPZGeoCube> *gelrp;
				TPZGeoEl *gel;
				TPZAutoPointer<TPZRefPattern> refpat;
				std::string filename;
				filename = REFPATTERNDIR;
				filename += "/3D_Hexa_Rib_Side_14.rpt";
				
				refpat = new TPZRefPattern(filename);
				if(!gRefDBase.FindRefPattern(refpat))
				{
					gRefDBase.InsertRefPattern(refpat);
				}
				gel = gmesh->ElementVec()[1];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);
				gel = gmesh->ElementVec()[2];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);

				filename = REFPATTERNDIR;
				filename += "/3D_Hexa_Rib_Side_12.rpt";
				
				refpat = new TPZRefPattern(filename);
				if(!gRefDBase.FindRefPattern(refpat))
				{
					gRefDBase.InsertRefPattern(refpat);
				}
				gel = gmesh->ElementVec()[3];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);
				gel = gmesh->ElementVec()[4];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);

				filename = REFPATTERNDIR;
				filename += "/3D_Hexa_Rib_Side_10.rpt";
				
				refpat = new TPZRefPattern(filename);
				if(!gRefDBase.FindRefPattern(refpat))
				{
					gRefDBase.InsertRefPattern(refpat);
				}
				gel = gmesh->ElementVec()[5];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);
				gel = gmesh->ElementVec()[6];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);

				filename = REFPATTERNDIR;
				filename += "/3D_Hexa_Rib_Side_08.rpt";
				
				refpat = new TPZRefPattern(filename);
				if(!gRefDBase.FindRefPattern(refpat))
				{
					gRefDBase.InsertRefPattern(refpat);
				}
				gel = gmesh->ElementVec()[8];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);

				gmesh->ResetConnectivities();
				gmesh->BuildConnectivity();
				// DELETING A CUBE 6th
				delete gmesh->ElementVec()[7];
			}
			break;
		case 3:
			// Refinar todo elemento cubo como quatro tetrahedros
			break;
	}
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	// INSERTING BOUNDARY ELEMENT IN THE INITIAL GEOMETRIC MESH
	switch(problem) {
		case 1:
		{
			// All boundary conditions are Dirichlet
			if(typeel==2) {    // Pyramids
				TPZGeoElBC gbc10(gmesh->ElementVec()[10],13,-1);
				TPZGeoElBC gbc11(gmesh->ElementVec()[11],13,-1);
				TPZGeoElBC gbc12(gmesh->ElementVec()[12],13,-1);
				TPZGeoElBC gbc13(gmesh->ElementVec()[13],13,-1);
				TPZGeoElBC gbc14(gmesh->ElementVec()[14],13,-1);
				TPZGeoElBC gbc15(gmesh->ElementVec()[16],13,-1);
				TPZGeoElBC gbc16(gmesh->ElementVec()[17],13,-1);
				TPZGeoElBC gbc17(gmesh->ElementVec()[17],16,-1);
				TPZGeoElBC gbc18(gmesh->ElementVec()[18],13,-1);
				TPZGeoElBC gbc19(gmesh->ElementVec()[18],16,-1);
				TPZGeoElBC gbc20(gmesh->ElementVec()[19],16,-1);
				TPZGeoElBC gbc21(gmesh->ElementVec()[20],13,-1);
				TPZGeoElBC gbc22(gmesh->ElementVec()[22],13,-1);
				TPZGeoElBC gbc23(gmesh->ElementVec()[23],13,-1);
				TPZGeoElBC gbc24(gmesh->ElementVec()[24],13,-1);
				TPZGeoElBC gbc25(gmesh->ElementVec()[26],13,-1);
				TPZGeoElBC gbc26(gmesh->ElementVec()[27],13,-1);
				TPZGeoElBC gbc27(gmesh->ElementVec()[28],13,-1);
				TPZGeoElBC gbc28(gmesh->ElementVec()[29],13,-1);
				TPZGeoElBC gbc29(gmesh->ElementVec()[29],15,-1);
				TPZGeoElBC gbc30(gmesh->ElementVec()[30],13,-1);
				TPZGeoElBC gbc32(gmesh->ElementVec()[31],15,-1);
				TPZGeoElBC gbc33(gmesh->ElementVec()[32],13,-1);
				TPZGeoElBC gbc34(gmesh->ElementVec()[32],16,-1);
				TPZGeoElBC gbc35(gmesh->ElementVec()[33],13,-1);
				TPZGeoElBC gbc36(gmesh->ElementVec()[34],13,-1);
				TPZGeoElBC gbc37(gmesh->ElementVec()[35],13,-1);
				TPZGeoElBC gbc38(gmesh->ElementVec()[36],13,-1);
			}
			else if(typeel==1) {            // Case in three Prisms
				TPZGeoElBC gbc10(gmesh->ElementVec()[9],16,-1);
				TPZGeoElBC gbc11(gmesh->ElementVec()[9],19,-1);
				TPZGeoElBC gbc12(gmesh->ElementVec()[10],15,-1);
				TPZGeoElBC gbc13(gmesh->ElementVec()[11],15,-1);
				TPZGeoElBC gbc15(gmesh->ElementVec()[11],17,-1);
				TPZGeoElBC gbc16(gmesh->ElementVec()[12],19,-1);
				TPZGeoElBC gbc17(gmesh->ElementVec()[13],15,-1);
				TPZGeoElBC gbc18(gmesh->ElementVec()[13],16,-1);
				TPZGeoElBC gbc19(gmesh->ElementVec()[14],15,-1);
				TPZGeoElBC gbc20(gmesh->ElementVec()[14],17,-1);
				TPZGeoElBC gbc21(gmesh->ElementVec()[15],15,-1);
				TPZGeoElBC gbc22(gmesh->ElementVec()[15],18,-1);
				TPZGeoElBC gbc23(gmesh->ElementVec()[16],16,-1);
				TPZGeoElBC gbc41(gmesh->ElementVec()[16],18,-1);
				TPZGeoElBC gbc42(gmesh->ElementVec()[16],19,-1);
				TPZGeoElBC gbc24(gmesh->ElementVec()[17],17,-1);
				TPZGeoElBC gbc25(gmesh->ElementVec()[17],19,-1);
				TPZGeoElBC gbc26(gmesh->ElementVec()[18],15,-1);
				TPZGeoElBC gbc27(gmesh->ElementVec()[18],16,-1);
				TPZGeoElBC gbc28(gmesh->ElementVec()[19],19,-1);
				TPZGeoElBC gbc29(gmesh->ElementVec()[20],17,-1);
				TPZGeoElBC gbc43(gmesh->ElementVec()[20],19,-1);
				TPZGeoElBC gbc30(gmesh->ElementVec()[21],16,-1);
				TPZGeoElBC gbc44(gmesh->ElementVec()[21],19,-1);
				TPZGeoElBC gbc45(gmesh->ElementVec()[22],15,-1);
				TPZGeoElBC gbc31(gmesh->ElementVec()[23],15,-1);
				TPZGeoElBC gbc32(gmesh->ElementVec()[23],17,-1);
				TPZGeoElBC gbc34(gmesh->ElementVec()[24],15,-1);
				TPZGeoElBC gbc46(gmesh->ElementVec()[24],19,-1);
				TPZGeoElBC gbc47(gmesh->ElementVec()[25],15,-1);
				TPZGeoElBC gbc35(gmesh->ElementVec()[25],16,-1);
				TPZGeoElBC gbc48(gmesh->ElementVec()[25],19,-1);
				TPZGeoElBC gbc49(gmesh->ElementVec()[26],15,-1);
				TPZGeoElBC gbc50(gmesh->ElementVec()[26],17,-1);
				TPZGeoElBC gbc36(gmesh->ElementVec()[26],19,-1);
				TPZGeoElBC gbc51(gmesh->ElementVec()[27],15,-1);
				TPZGeoElBC gbc52(gmesh->ElementVec()[27],16,-1);
				TPZGeoElBC gbc38(gmesh->ElementVec()[28],16,-1);
				TPZGeoElBC gbc54(gmesh->ElementVec()[28],19,-1);
				TPZGeoElBC gbc40(gmesh->ElementVec()[29],17,-1);
				TPZGeoElBC gbc55(gmesh->ElementVec()[29],19,-1);
			}
			else {
				// bottom condition
				TPZGeoElBC gbc10(gmesh->ElementVec()[1],20,-1);
				TPZGeoElBC gbc20(gmesh->ElementVec()[2],20,-1);
				TPZGeoElBC gbc30(gmesh->ElementVec()[3],20,-1);
				TPZGeoElBC gbc40(gmesh->ElementVec()[4],20,-1);
				// face 1 lateral left
				TPZGeoElBC gbc11(gmesh->ElementVec()[1],21,-1);
				TPZGeoElBC gbc21(gmesh->ElementVec()[2],21,-1);
				TPZGeoElBC gbc31(gmesh->ElementVec()[5],21,-1);
				TPZGeoElBC gbc41(gmesh->ElementVec()[6],21,-1);
				// face 2 lateral front
				TPZGeoElBC gbc12(gmesh->ElementVec()[2],22,-1);
				TPZGeoElBC gbc22(gmesh->ElementVec()[3],22,-1);
				TPZGeoElBC gbc32(gmesh->ElementVec()[6],22,-1);
				TPZGeoElBC gbc42(gmesh->ElementVec()[8],22,-1);
				// face 3 lateral right
				TPZGeoElBC gbc13(gmesh->ElementVec()[3],23,-1);
				TPZGeoElBC gbc23(gmesh->ElementVec()[4],23,-1);
				TPZGeoElBC gbc33(gmesh->ElementVec()[6],23,-1);
				TPZGeoElBC gbc43(gmesh->ElementVec()[8],23,-1);
				// face 4 lateral back
				TPZGeoElBC gbc14(gmesh->ElementVec()[1],24,-1);
				TPZGeoElBC gbc24(gmesh->ElementVec()[4],24,-1);
				TPZGeoElBC gbc34(gmesh->ElementVec()[5],24,-1);
				TPZGeoElBC gbc44(gmesh->ElementVec()[8],24,-1);
				// top condition
				TPZGeoElBC gbc15(gmesh->ElementVec()[3],25,-1);
				TPZGeoElBC gbc25(gmesh->ElementVec()[5],25,-1);
				TPZGeoElBC gbc35(gmesh->ElementVec()[6],25,-1);
				TPZGeoElBC gbc45(gmesh->ElementVec()[8],25,-1);
			}
		}
		break;
		case 2:
		{
			// bottom condition
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],20,-1);
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],20,-1);
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],20,-1);
			TPZGeoElBC gbc40(gmesh->ElementVec()[4],20,-1);
			// face 1 lateral left
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],21,-1);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],21,-1);
			TPZGeoElBC gbc31(gmesh->ElementVec()[5],21,-1);
			TPZGeoElBC gbc41(gmesh->ElementVec()[6],21,-1);
			// face 2 lateral front
			TPZGeoElBC gbc12(gmesh->ElementVec()[2],22,-1);
			TPZGeoElBC gbc22(gmesh->ElementVec()[3],22,-1);
			TPZGeoElBC gbc32(gmesh->ElementVec()[6],22,-1);
			TPZGeoElBC gbc42(gmesh->ElementVec()[8],22,-2);
			// face 3 lateral right
			TPZGeoElBC gbc13(gmesh->ElementVec()[3],23,-1);
			TPZGeoElBC gbc23(gmesh->ElementVec()[4],23,-1);
			TPZGeoElBC gbc33(gmesh->ElementVec()[6],23,-2);
			TPZGeoElBC gbc43(gmesh->ElementVec()[8],23,-1);
			// face 4 lateral back
			TPZGeoElBC gbc14(gmesh->ElementVec()[1],24,-1);
			TPZGeoElBC gbc24(gmesh->ElementVec()[4],24,-1);
			TPZGeoElBC gbc34(gmesh->ElementVec()[5],24,-1);
			TPZGeoElBC gbc44(gmesh->ElementVec()[8],24,-1);
			// top condition
			TPZGeoElBC gbc15(gmesh->ElementVec()[3],25,-2);
			TPZGeoElBC gbc25(gmesh->ElementVec()[5],25,-1);
			TPZGeoElBC gbc35(gmesh->ElementVec()[6],25,-1);
			TPZGeoElBC gbc45(gmesh->ElementVec()[8],25,-1);
		}
		break;
	}
	
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	return gmesh;
}

void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,REAL radius,int ntyperefs) {
	TPZManVector<REAL> point(3,0.);
	REAL r = 0.0;
	bool isdefined = true;

	if(ntyperefs==2) {
		// To refine elements with distance from center to point less than radius
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
}
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
	TPZManVector<REAL> point(3,0.);
	REAL r = 0.0, radius = 1.2;
	int i;
	bool isdefined = false;
	if(ntyperefs==2) {
		for(i=0;i<nref;i+=2) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			if(nref < 5) radius *= 0.4;
			else if(nref < 7) radius *= 0.2;
			else radius *= 0.1;
		}
		if(i==nref) {
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
		}
	}
	else {
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

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(!gel || gel->Dimension()!=dim || gel->HasSubElement()) continue;
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
	TPZStepSolver<STATE> stepstep;
	TPZBandStructMatrix matrix(an.Mesh());
	an.SetStructuralMatrix(matrix);
	stepstep.SetDirect(ELU);
	an.SetSolver(stepstep);
}



// CONSTRUCTING A GEOMETRIC MESH AS FICHERA CORNER, USING EXTRUDE AND MERGE MESHES

/// ALERT !!!!!!!!!    ALERT !!!!!!!!

// WAS DETECTED PROBLEM IN MERGE AND READ METHOD OF THE TPZGENGRID CLASS

/// ALERT !!!!!!!!!    ALERT !!!!!!!!
/// ALERT !!!!!!!!!    ALERT !!!!!!!!
/// ALERT !!!!!!!!!    ALERT !!!!!!!!

TPZGeoMesh *ConstructingFicheraCorner(REAL InitialL, REAL InitialH, bool print) {
	
	TPZManVector<REAL> x0(3,-InitialL), x1(3,-InitialL);  // Corners of the rectangular mesh.
	TPZManVector<int> nx(2,2);   // subdivisions in X and in Y. 
	
	// Mesh as a hexahedra with large edge (2 units) and was cutting a small hexahedra with edge with 1 unit:
	// 1) The bottom rectangular mesh has four corners: (1,-1,-1), (1,1,-1), (-1,1,-1) and (-1,-1,-1)
	// and was divides in four segments on X and four on Y, then hx = 0.5 and hy = 0.5
	if(print) out << "Generating (1) bottom geometric mesh bi-dimensional ...\n";
	
	
	
	TPZGeoMesh* gmesh = new TPZGeoMesh;
	x0[0] = InitialL; x0[1] = x0[2] = -InitialL;
	x1[0] = x1[2] = - InitialL; x1[1] = InitialL;
	TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
	gen.SetElementType(EQuadrilateral);       // type = 0 means rectangular elements
	gen.Read(gmesh);             // generating grid in gmesh
	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	if(print) out << "... (Extruding) first geometric mesh three-dimensional...\n";
	TPZExtendGridDimension gmeshextend(gmesh,InitialH);
	TPZGeoMesh *gmesh3D = gmeshextend.ExtendedMesh(1,1,2);
	gmesh3D->SetName("First Mesh Extruded");
	/* Dirichlet boundary condition for inferior cubes
	 TPZGeoElBC gelbcD(gmesh3D->ElementVec()[0],21,-1);   // Dirichlet condition for hexahedra//
	 TPZGeoElBC gelbcD1(gmesh3D->ElementVec()[0],24,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD2(gmesh3D->ElementVec()[1],21,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD3(gmesh3D->ElementVec()[1],22,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD4(gmesh3D->ElementVec()[2],23,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD5(gmesh3D->ElementVec()[2],24,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD6(gmesh3D->ElementVec()[2],25,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD7(gmesh3D->ElementVec()[3],22,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD8(gmesh3D->ElementVec()[3],23,-1);   // Dirichlet condition for hexahedra
	 */
	// 2) The left upper rectangular mesh has four corners: (1,-1,0), (1,0,0), (-1,0,0) and (-1,-1,0)
	// and was divides in one segment on X and two on Y, then hx = hy = 1.
	TPZGeoMesh *gmesh2 = new TPZGeoMesh;
	if(print) out << "Generating (2) left upper geometric mesh bi-dimensional ...\n";
	x0[0] = InitialL; x0[1] = -InitialL; x0[2] = 0.;
	x1[0] = -InitialL; x1[1] = x1[2] = 0.;
	nx[0] = 2; nx[1] = 1;
	gen.SetData(nx,x0,x1,EQuadrilateral);
	gen.Read(gmesh2);             // generating grid in gmesh
	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	if(print) out << "... (Extruding) first geometric mesh three-dimensional...\n";
	TPZExtendGridDimension gmeshextend2(gmesh2,InitialH);
	TPZGeoMesh *gmesh3D2 = gmeshextend2.ExtendedMesh(1,2,1);
	gmesh3D2->SetName("Second Mesh Extruded");
	/* Dirichlet boundary condition for inferior cubes
	 TPZGeoElBC gelbcD10(gmesh3D2->ElementVec()[0],21,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD11(gmesh3D2->ElementVec()[0],23,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD12(gmesh3D2->ElementVec()[0],24,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD13(gmesh3D2->ElementVec()[1],21,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD14(gmesh3D2->ElementVec()[1],22,-1);   // Dirichlet condition for hexahedra
	 */
	// 3) The right upper rectangular mesh has four corners: (0,0,0), (0,1,0), (-1,1,0) and (-1,0,0)
	// and was divides in one segment on X and two on Y, then hx = hy = 1.
	TPZGeoMesh *gmesh3 = new TPZGeoMesh;
	if(print) out << "Generating (3) left upper geometric mesh bi-dimensional ...\n";
	x0[0] = x0[1] = x0[2] = 0.;
	x1[0] = -InitialL; x1[1] = InitialL; x1[2] = 0.;
	nx[0] = nx[1] = 1;
	gen.SetData(nx,x0,x1,EQuadrilateral);       // type = 0 means rectangular elements
	gen.Read(gmesh3);             // generating grid in gmesh
	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	if(print) out << "... (Extruding) first geometric mesh three-dimensional...\n";
	TPZExtendGridDimension gmeshextend3(gmesh3,InitialH);
	TPZGeoMesh *gmesh3D3 = gmeshextend3.ExtendedMesh(1,2,1);
	gmesh3D3->SetName("Third Mesh Extruded");
	/* Dirichlet boundary condition for inferior cubes
	 TPZGeoElBC gelbcD20(gmesh3D3->ElementVec()[0],22,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD21(gmesh3D3->ElementVec()[0],23,-1);   // Dirichlet condition for hexahedra
	 TPZGeoElBC gelbcD22(gmesh3D3->ElementVec()[0],24,-1);   // Dirichlet condition for hexahedra
	 */
	// Merge two meshes
	gen.MergeGeoMesh(gmesh3D,gmesh3D2,1);
	gen.MergeGeoMesh(gmesh3D,gmesh3D3,1);
	
	// Cleaning unnecessary meshes
	delete gmesh3D2;
	delete gmesh3D3;
	
	return gmesh3D;
}

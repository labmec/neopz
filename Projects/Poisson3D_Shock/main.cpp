/**
 * @file Poisson 3D in hexahedra with shock problem
 */

#include "pzshapelinear.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

#include "pzbstrmatrix.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmaterial.h"
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

using namespace std;
using namespace pzshape;

/**
 * @addtogroup Tutorials
 * @{
 */

// Global variable
int gLMax;
int NUniformRefs = 2;
bool anothertests = false;
int nstate = 2;

// Alfa -> Coefficient of the arctang argument
REAL ALFA = 20.;

char saida[512];
ofstream out("ConsolePoisson3D.txt");   // output file from console  -> Because it has many energy faults

/** Printing level */
int gPrintLevel = 0;
bool gDebug = false;

void ExactShock(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void BCDirichletShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol);
void BCNeumannLateralFrontShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol);
void BCNeumannLateralRightShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol);
void BCNeumannTopShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol);
void FforcingShock(const TPZVec<REAL> &x, TPZVec<REAL> &f);

void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=true, const int matidtodivided=1);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &points,REAL r,REAL &distance,bool &isdefined);
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);

TPZGeoMesh *ConstructingCubePositiveOctant(REAL L,int typeel);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);

void formatTimeInSec(char *strtime,int timeinsec);

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
	ofstream fileerrors("ErrorsHPProcess3D.txt");
    
	TPZVec<REAL> ervec(100,0.0);
	// Printing computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	REAL InitialL = 1.0;
	int nref, NRefs = 10;
	int nthread, NThreads = 4;
	int dim = 3;
    for(int typeel=0;typeel<1;typeel++) {
        for(nref=0;nref<NRefs;nref++) {
            for(int ntyperefs=0;ntyperefs<2;ntyperefs++) {
				if(nref > 6) nthread = 2*NThreads;
				else nthread = 2;
				
				// Initializing the generation mesh process
				time(&sttime);
				// Constructing geometric mesh as hexahedra
				cout << "\nConstructing Shock problem in cube. Refinement: " << nref+1 << " Threads: " << nthread << " TypeRef: " << ntyperefs << " TypeElement: " << typeel << endl;
				TPZGeoMesh *gmesh3D = ConstructingCubePositiveOctant(InitialL,typeel);
				// h_refinement
				// Refining near to the origin
				RefiningNearCircunference(dim,gmesh3D,nref,ntyperefs);
				if(nref==NRefs-1) {
					sprintf(saida,"gmesh_%d_%d.vtk",nref,ntyperefs);
					PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
				}
				
				// Creating computational mesh
				/** Set polynomial order */
				int p = 8, pinit;
				pinit = p;
				TPZCompEl::SetgOrder(1);
				TPZCompMesh *cmesh = CreateMesh(gmesh3D,dim,1);
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
				
				TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
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
				sprintf(pp,"TR%dE%dP%2dH%2dP%d",ntyperefs,typeel,nthread,nref,pinit);
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
				an.SetExact(ExactShock);
				fileerrors << "Refinement: " << nref << "  Dimension: " << dim << "  NEquations: " << cmesh->NEquations();
				an.PostProcessError(ervec,out);
				for(int rr=0;rr<ervec.NElements();rr++)
					fileerrors << "  Error_" << rr+1 << ": " << ervec[rr]; 
				fileerrors << "  TimeElapsed: " << time_elapsed << " <-> " << tempo << std::endl;

				delete cmesh;
				delete gmesh3D;
			}
		}
	}
	out.close();
	fileerrors.close();
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////   FICHERA CORNER - Problem as Anders Solin Presentation   ///////////////////
////////////////////////////////////////////////////////////////////////////////////////
void ExactShock(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	sol[0] = atan(ALFA * ( R0 - sqrt(3.)) );
	REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	if(IsZero(den))
		DebugStop();
	dsol(0,0) = (ALFA*(x[0]-C0[0]))/den;
	dsol(1,0) = (ALFA*(x[1]-C0[1]))/den;
	dsol(2,0) = (ALFA*(x[2]-C0[2]))/den;
}

void BCDirichletShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	bcsol[0] = atan(ALFA * ( R0 - sqrt(3.)) );
}
void BCNeumannLateralFrontShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	if(IsZero(den))
		DebugStop();
	bcsol[0] = (ALFA*(x[0]-C0[0]))/den;
}
void BCNeumannLateralRightShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	if(IsZero(den))
		DebugStop();
	bcsol[0] = (ALFA*(x[1]-C0[1]))/den;
}
void BCNeumannTopShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	REAL den = R0 * (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	if(IsZero(den))
		DebugStop();
	bcsol[0] = (ALFA*(x[2]-C0[2]))/den;
}

/** NOTE: Forcing function in TPZMatPoisson3d is negative */
void FforcingShock(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
	TPZVec<REAL> C0(3,-0.25);
	
	REAL R0 = sqrt ((x[0]-C0[0])*(x[0]-C0[0]) + (x[1]-C0[1])*(x[1]-C0[1]) + (x[2]-C0[2])*(x[2]-C0[2]));
	REAL temp =  (1. + ALFA*ALFA*(R0-sqrt(3.))*(R0-sqrt(3.)));
	REAL den = R0*temp*temp;
	if(IsZero(den))
		DebugStop();
	f[0] = (2*ALFA*(1.+(3.*ALFA*ALFA)-(ALFA*ALFA*sqrt(3.)*R0)))/den;
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
	TPZFMatrix<REAL> val1(dim,dim,0.),val2(dim,1,0.);
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

TPZGeoMesh *ConstructingCubePositiveOctant(REAL InitialL,int typeel) {
	
	// CREATING A CUBE WITH MASS CENTER (0.5, 0.5, 0.5) AND VOLUME = INITIALL*INITIALL*INITIALL 
	REAL co[8][3] = {
		{0.,0.,0.},
		{InitialL,0.,0.},
		{InitialL,InitialL,0.},
		{0.,InitialL,0.},
		{0.,0.,InitialL},
		{InitialL,0.,InitialL},
		{InitialL,InitialL,InitialL},
		{0.,InitialL,InitialL}
	};
	int indices[1][8] = {{0,1,2,3,4,5,6,7}};
	
	const int nelem = 1;
	int nnode = 8;
	
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
		TPZManVector<int> nodind(8);
		for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
		int index;
		elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
	}
	gmesh->BuildConnectivity();
	
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

	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	return gmesh;
}

void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
	TPZManVector<REAL> point(3,-0.25);
	REAL r = sqrt(3.0), radius = .6;
	int i;
	bool isdefined = true;
	if(ntyperefs) {
		for(i=0;i<nref;i+=2) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			radius *= 0.3;
		}
		if(i==nref) {
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			radius *= 0.5;
		}
	}
	else {
		for(i=0;i<nref+1;i++) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			radius *= 0.5;
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
		sprintf(strtime,"%d a, %d m, %d d, %d:%d:%d",anos,meses,dias,horas,minutos,segundos);
	else {
		if(meses) 
			sprintf(strtime,"%d m, %d d, %d:%d:%d",meses,dias,horas,minutos,segundos);
		else {
			if(dias)
				sprintf(strtime,"%d d, %d:%d:%d",dias,horas,minutos,segundos);
			else
				sprintf(strtime,"%d:%d:%d",horas,minutos,segundos);
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
	TPZStepSolver<REAL> step;
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



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
REAL alfa = M_PI/6.;
bool anothertests = false;
int nstate = 2;

char saida[512];
ofstream out("ConsolePoisson3D.txt");   // output file from console  -> Because it has many energy faults

/** Printing level */
int gPrintLevel = 0;
bool gDebug = false;

void Exact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void ExactShock(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void BCShock(const TPZVec<REAL> &x, TPZVec<REAL> &sol);

void Ff(const TPZVec<REAL> &x, TPZVec<REAL> &f);

//void InitializeSolver(TPZAnalysis &an);
void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh *cmesh);
void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=true, const int matidtodivided=1);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &points,REAL r,REAL &distance,bool &isdefined);

TPZGeoMesh *ConstructingFicheraCorner(REAL L,bool print = false);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);

void ComputeSolutionError(REAL &error,REAL &errorL2,TPZCompMesh *cmesh,int dim);
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
	int i, nref, NRefs = 7;
	int nthread, NThreads = 4;
	int dim = 3;
    for(int typeel=0;typeel<4;typeel++) {
		for(int ntyperefs=0;ntyperefs<2;ntyperefs++) {
			for(nref=0;nref<NRefs;nref++) {
				if(nref > 5) nthread = NThreads;
				else nthread = 1;
				
				// Initializing the generation mesh process
				time(&sttime);
				// Constructing geometric mesh as Fichera corner using hexahedra
				cout << "\nConstructing Shock problem in cube. Refinement: " << nref+1 << " Threads: " << nthread << " TypeRef: " << ntyperefs << " TypeElement: " << typeel << endl;
				TPZGeoMesh *gmesh3D = ConstructingFicheraCorner(InitialL,typeel);
				// h_refinement
				TPZManVector<REAL> point(3,0.);
				REAL r = 0.0, radius = 0.9;
				bool isdefined = false;
				if(ntyperefs) {
					for(i=0;i<nref;i+=2) {
						// To refine elements with center near to points than radius
						//				RefineGeoElements(2,gmesh,Points,radius,isdefined);
						// Para refinar elementos con centro tan cerca de la circuferencia cuanto radius 
						RefineGeoElements(3,gmesh3D,point,r,radius,isdefined);
						RefineGeoElements(3,gmesh3D,point,r,radius,isdefined);
						radius *= 0.6;
					}
					if(i==nref) {
						RefineGeoElements(3,gmesh3D,point,r,radius,isdefined);
						radius *= 0.6;
					}
				}
				else {
					for(i=0;i<nref+1;i++) {
						// To refine elements with center near to points than radius
						// Para refinar elementos con centro tan cerca de la circuferencia cuanto radius 
						RefineGeoElements(3,gmesh3D,point,r,radius,isdefined);
						radius *= 0.6;
					}
				}
				// Constructing connectivities
				gmesh3D->ResetConnectivities();
				gmesh3D->BuildConnectivity();
				sprintf(saida,"meshextrudedmerged_%d_%d.vtk",nref,ntyperefs);
				PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
				
				// Creating computational mesh
				/** Set polynomial order */
				int p = 4, pinit;
				pinit = p;
				TPZCompEl::SetgOrder(1);
				TPZCompMesh *cmesh = CreateMesh(gmesh3D,dim,1);
				cmesh->SetName("Computational mesh for Fichera problem");
				dim = cmesh->Dimension();
				
				// Primeiro sera calculado o mayor nivel de refinamento. Remenber, the first level is zero level.
				// A cada nivel disminue em uma unidade o p, mas não será menor de 1.
				int level, highlevel = 0;
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
				char pp[3];
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
				
				vecnames.Push("Derivate");
				vecnames.Push("Flux");
				vecnames.Push("MinusKGradU");
				// END Determining the name of the variables
				
				an.DefineGraphMesh(dim,scalarnames,vecnames,filename);
				
				an.PostProcess(0,dim);
				
				// Computing error
				an.SetExact(ExactShock);
				fileerrors << "Refinement: " << nref << "  Dimension: " << dim << "  NEquations: " << cmesh->NEquations();
				an.PostProcess(ervec,out);
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
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	sol[0] = sqrt( sqrt (quad_r) );
	if(!IsZero(sol[0])) {
		REAL den = sol[0]*sol[0]*sol[0];
		dsol(0,0) = .5*x[0]/den;
		dsol(1,0) = .5*x[1]/den;
		dsol(1,0) = .5*x[2]/den;
	}
	else {
		dsol(0,0) = dsol(1,0) = dsol(2,0) = 0.;
	}
}

void BCShock(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	bcsol[0] = sqrt( sqrt (quad_r) );	
}

void Ff(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	REAL raiz = sqrt( sqrt(quad_r));
	f[0] = -3./(4.*(raiz*raiz*raiz));
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
		mat->SetForcingFunction(new TPZDummyFunction<STATE>(Ff));
	}
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	
	// Boundary conditions
	// Dirichlet
	TPZAutoPointer<TPZFunction<STATE> > FunctionBC = new TPZDummyFunction<STATE>(BCShock);
	TPZFMatrix<REAL> val1(dim,dim,0.),val2(dim,1,0.);
	TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
	bnd->SetForcingFunction(FunctionBC);
	cmesh->InsertMaterialObject(bnd);
	
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->ExpandSolution();
	return cmesh;
}

TPZGeoMesh *ConstructingFicheraCorner(REAL InitialL, bool print) {
	
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
	
	// SUBDIVIDING A CUBE
	TPZVec<TPZGeoEl*> sub;
	gmesh->ElementVec()[0]->Divide(sub);
	
	// DELETING A CUBE 6th
	delete gmesh->ElementVec()[7];
	
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	// INSERTING BOUNDARY ELEMENT IN THE INITIAL GEOMETRIC MESH
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
	
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	return gmesh;
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



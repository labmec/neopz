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
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

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

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "poissondesacoplados.h"


#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;

const int matId = 1;
const int dirichlet = 0;
const int neumann = 1;
const int bcDL = -1;
const int bcN = -2;
const int bcDR = -3;

const int lagrangemat = 2;
const int interfacemat = 3;

TPZGeoMesh *MalhaGeom(bool interface);

TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh*MalhaCompDois(TPZGeoMesh * gmesh, int pOrder);
TPZCompMesh *MalhaCompComInterf(TPZGeoMesh * gmesh,int pOrder);
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);
void PosProcess(TPZAnalysis &an, std::string plotfile);
void RefinamentoUniforme(TPZGeoMesh  *gMesh, int nh);
void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh, int MatId, int indexEl);
void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);
void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file);
void GeoElMultiphysicVec(TPZManVector<TPZCompMesh  *> cmeshVec,std::set <int> &geoelVec);
void AddElements(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh);
void AddConnects(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh);
void TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
void TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);

void BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, int LagrangeMat, int InterfaceMat);

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	std::string logs("log4cxx.doubleprojection1d");
	InitializePZLOG();
#endif
	
	int p =1;
	//primeira malha
	
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = MalhaGeom(false);
	ofstream arg1("gmeshZero.txt");
	gmesh->Print(arg1);
	
		
	// First computational mesh
	TPZCompMesh * cmesh1= MalhaComp(gmesh,  p);
	ofstream arg2("cmesh1.txt");
	cmesh1->Print(arg2);
			
	// Second computational mesh
	TPZCompMesh * cmesh2 = MalhaCompDois(gmesh, p+1);
	ofstream arg3("cmesh2.txt");
	cmesh2->Print(arg3);
	
	// Cleaning reference of the geometric mesh to cmesh1
	gmesh->ResetReference();
	cmesh1->LoadReferences();
	//RefinUniformElemComp(cmesh1,2);
	// Refine the 7th element of the cmesh1
	RefinElemComp(cmesh1, 7);
	// Refine the 10th element of the cmesh1
	RefinElemComp(cmesh1, 10);
	// Adjust the boundary elements after refine
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();
	
	ofstream arg4("cmesh12.txt");
	cmesh1->Print(arg4);
	ofstream arg5("gmesh2.txt");
	gmesh->Print(arg5);
	ofstream file3("malhageo1.vtk");
	PrintGMeshVTK(gmesh, file3);
	
	// Cleaning reference to cmesh2
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	//refinamento uniform
	//RefinUniformElemComp(cmesh2,3);
	// Refine 6th and 7th elements (as uniform refine)
	RefinElemComp(cmesh2, 6);
	RefinElemComp(cmesh2, 7);
	cmesh2->AdjustBoundaryElements();
	cmesh2->CleanUpUnconnectedNodes();
	
	ofstream arg6("cmesh22.txt");
	cmesh2->Print(arg6);
	ofstream arg7("gmesh3.txt");
	gmesh->Print(arg7);
	ofstream file5("malhageo2.vtk");
	PrintGMeshVTK(gmesh, file5);

	
	//--- Resolver usando a primeira malha computacional ---
	TPZAnalysis an1(cmesh1);
	SolveSist(an1, cmesh1);
	std::string plotfile("saidaSolution_cmesh1.vtk");
	PosProcess(an1, plotfile);
	//---------------------------
	
	//--- Resolver usando a segunda malha computacional ---
	TPZAnalysis an2(cmesh2);
	SolveSist(an2, cmesh2);
	std::string plotfile2("saidaSolution_cmesh2.vtk");
	PosProcess(an2, plotfile2);
	//---------------------------
	
	// List of the computational meshes
	TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh1;
	meshvec[1] = cmesh2;
	
	// Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
	mphysics->SetAllCreateFunctionsMultiphysicElem();
	
	int MatId = 1;
	TwoUncoupledPoisson *mymaterial = new TwoUncoupledPoisson(MatId, 2);
	mymaterial->SetParameters(-1., -0.1);
	mymaterial->SetInternalFlux(8.,0.);
	ofstream argm("mymaterial.txt");
	mymaterial->Print(argm);
	TPZAutoPointer<TPZMaterial> mat(mymaterial);
	mphysics->InsertMaterialObject(mat);
	
	///Inserir condicao de contorno
	TPZFMatrix val1(2,2,0.), val2(2,1,0.);
	val2(0,0)=0.;
	val2(1,0)=0.;
	TPZAutoPointer<TPZMaterial> BCondN = mymaterial->CreateBC(mat, bcN,neumann, val1, val2);
	mphysics->InsertMaterialObject(BCondN);
	
	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	val22(0,0)=0.;
	val22(1,0)=2.;
	TPZAutoPointer<TPZMaterial> BCondDL = mymaterial->CreateBC(mat, bcDL,dirichlet, val12, val22);
	mphysics->InsertMaterialObject(BCondDL);
	
	TPZFMatrix val13(2,2,0.), val23(2,1,0.);
	val23(0,0)=0.;
	val23(1,0)=1.;
	TPZAutoPointer<TPZMaterial> BCondDR = mymaterial->CreateBC(mat, bcDR,dirichlet, val13, val23);
	mphysics->InsertMaterialObject(BCondDR);
	
	
	mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
	
	// Creating multiphysic elements into mphysics computational mesh
	AddElements(meshvec, mphysics);
	AddConnects(meshvec,mphysics);
	TransferFromMeshes(meshvec, mphysics);
	
//#ifdef LOG4CXX
//    {
//        std::stringstream sout;
//        mphysics->Print(sout);
//        LOGPZ_DEBUG(logger, sout.str())
//    }
//#endif

	ofstream arg8("mphysic.txt");
	mphysics->Print(arg8);
	
	
	TPZAnalysis an(mphysics);
	SolveSist(an, mphysics);
	
	//--- pos-process -----
	TransferFromMultiPhysics(meshvec, mphysics);

	TPZManVector<std::string,10> scalnames(2), vecnames(2);
	scalnames[0] = "SolutionU";
	scalnames[1] = "SolutionP";
	vecnames[0]= "DerivateU";
	vecnames[1]= "DerivateP";
	
	std::string plotfile3("saidaSolution.vtk");
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile3);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
	//----------
	
		
	
	//------------------------------------
//	std::set<int> geoelVec;
//	std::set<int> refIndexVec;
//	TPZManVector<TPZCompMesh *,2> cmeshVec(2);
//	cmeshVec[0]=cmesh1;
//	cmeshVec[1]=cmesh2;
//	GeoElMultiphysicVec(cmeshVec, geoelVec);
//	
//	set<int>::iterator it;
//	cout << "myset contains:";
//	for (it=geoelVec.begin() ; it != geoelVec.end(); it++ )
//		cout << " " << *it;
//	cout << endl;
	
	
	//TPZMultiphysicsCompEl <pzgeom::TPZGeoQuad> *mmesh = new TPZMultiphysicsCompEl <pzgeom::TPZGeoQuad>();
//	TPZManVector<TPZTransform> tr;
//	mmesh->AffineTransform(tr);
//	
//#ifdef LOG4CXX
//    {
//        std::stringstream out;
//       tr[0].PrintInputForm(out);
//        LOGPZ_DEBUG(logger, out.str())
//    }
//#endif


	// -------- validar o metodo BuildHybridMesh() ------------
	//TPZGeoMesh * gmesh = MalhaGeom(true);
//	ofstream arg1("gmeshZero.txt");
//	gmesh->Print(arg1);
//	
//	TPZCompMesh * cmesh= MalhaCompComInterf(gmesh,  p);
//	ofstream arg2("cmesh.txt");
//	cmesh->Print(arg2);
//	
//	gmesh->ResetReference();
//	cmesh->LoadReferences();
//	RefinElemComp(cmesh, 7);
//	cmesh->AdjustBoundaryElements();
//	cmesh->CleanUpUnconnectedNodes();
//	
//	ofstream arg22("cmeshRefUniforme.txt");
//	cmesh->Print(arg22);
//	
//	ofstream arg12("gmeshRefUniforme.txt");
//	gmesh->Print(arg12);
//	ofstream file1("malhageoRefUniforme.vtk");
//	PrintGMeshVTK(gmesh, file1);
//		
//	std::set<int> MaterialIDs;
//	MaterialIDs.insert(matId);
//	MaterialIDs.insert(lagrangemat);
//	MaterialIDs.insert(interfacemat);
//		
//	set<int>::iterator it;
//	cout << "myset contains:";
//	for (it=MaterialIDs.begin() ; it != MaterialIDs.end(); it++ )
//		cout << " " << *it;
//		cout << endl;
//	
//	BuildHybridMesh(cmesh, MaterialIDs, lagrangemat, interfacemat);
//	
//	ofstream arg4("cmeshComInterface.txt");
//	cmesh->Print(arg4);
//	ofstream arg5("gmesh2.txt");
//	gmesh->Print(arg5);
//	ofstream file3("malhageo1.vtk");
//	PrintGMeshVTK(gmesh, file3);
	
	return EXIT_SUCCESS;
}


TPZGeoMesh *MalhaGeom(bool interface)
{
	
	int Qnodes = 6;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolLine(2);
	
	//indice dos nos
	int id = 0;
	REAL valx;
	REAL dx=0.5;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = 1. - xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,1. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcN,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcN,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcDR,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 4;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcN,*gmesh);
	id++;
	
	TopolLine[0] = 4;
	TopolLine[1] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcN,*gmesh);
	id++;
	
	TopolLine[0] = 5;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcDL,*gmesh);
	id++;
	
	TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 4;
	TopolQuad[3] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
	id++;
	
	TopolQuad[0] = 1;
	TopolQuad[1] = 2;
	TopolQuad[2] = 3;
	TopolQuad[3] = 4;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
	
	//criar elemento de interface
	if (interface==true) {
		gmesh->AddInterfaceMaterial(matId,lagrangemat, interfacemat); //Adicionar um material de interface associados aos elementos matId e lagrangemat do material.
		gmesh->AddInterfaceMaterial(lagrangemat,matId, interfacemat);
	}
	
	gmesh->BuildConnectivity();
			
	//ofstream arg("gmesh.txt");
//	gmesh->Print(arg);
		
	return gmesh;
	
}

TPZCompMesh*MalhaComp(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	TPZAutoPointer<TPZMaterial> mat(material);
	
	REAL diff = -1.;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
	REAL flux = 8.;
	
	material->SetParameters(diff, conv, convdir);
	material->SetInternalFlux( flux);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	

	///Inserir condicao de contorno
	TPZFMatrix val1(2,2,0.), val2(2,1,0.);
	REAL uN=0.;
	val2(0,0)=uN;
	//val2(1,0)=uN;
	TPZAutoPointer<TPZMaterial> BCondN = material->CreateBC(mat, bcN,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondN);
	
	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	REAL uD=0.;
	val22(0,0)=uD;
	//val22(1,0)=uD;
	TPZAutoPointer<TPZMaterial> BCondDL = material->CreateBC(mat, bcDL,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDL);
	
	TPZAutoPointer<TPZMaterial> BCondDR = material->CreateBC(mat, bcDR,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDR);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}

TPZCompMesh*MalhaCompDois(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	TPZAutoPointer<TPZMaterial> mat(material);
	
	REAL diff = -0.1;
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
	TPZFMatrix val1(2,2,0.), val2(2,1,0.);
	REAL uN=0.;
	val2(0,0)=uN;
	//val2(1,0)=uN;
	TPZAutoPointer<TPZMaterial> BCondN = material->CreateBC(mat, bcN,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondN);
	
	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	REAL uDL=2.;
	val22(0,0)=uDL;
	//val22(1,0)=uD;
	TPZAutoPointer<TPZMaterial> BCondDL = material->CreateBC(mat, bcDL,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDL);
	
	REAL uDR=1.;
	val22(0,0)=uDR;
	TPZAutoPointer<TPZMaterial> BCondDR = material->CreateBC(mat, bcDR,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDR);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}

TPZCompMesh*MalhaCompComInterf(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	
	TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
	TPZMatPoisson3d *matlagrange = new TPZMatPoisson3d(lagrangemat,1); 
	TPZMatPoisson3d *matinterface = new TPZMatPoisson3d(interfacemat,1); 
	
	TPZAutoPointer<TPZMaterial> mat1(material);
	TPZAutoPointer<TPZMaterial> mat2(matlagrange);
	TPZAutoPointer<TPZMaterial> mat3(matinterface);
	
	REAL diff = -1.;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
	REAL flux = 8.;
	
	material->SetParameters(diff, conv, convdir);
	material->SetInternalFlux( flux);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	
	cmesh->SetAllCreateFunctionsContinuous();
	
	cmesh->InsertMaterialObject(mat1);
	cmesh->InsertMaterialObject(mat2);
	cmesh->InsertMaterialObject(mat3);
	
	
	///Inserir condicao de contorno
	TPZFMatrix val1(2,2,0.), val2(2,1,0.);
	REAL uN=0.;
	val2(0,0)=uN;
	//val2(1,0)=uN;
	TPZAutoPointer<TPZMaterial> BCondN = material->CreateBC(mat1, bcN,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondN);
	
	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	REAL uD=0.;
	val22(0,0)=uD;
	//val22(1,0)=uD;
	TPZAutoPointer<TPZMaterial> BCondDL = material->CreateBC(mat1, bcDL,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDL);
	
	TPZAutoPointer<TPZMaterial> BCondDR = material->CreateBC(mat1, bcDR,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDR);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}

#define VTK
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
		
	//Saida de Dados: solucao e  grafico no VTK 
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
	
//#ifdef VTK
//	{
//		TPZManVector<std::string,10> scalnames(1), vecnames(1);
//		scalnames[0] = "Solution";
//		vecnames[0]= "Derivate";
//		
//		std::string plotfile("saidaSolution.vtk");
//		const int dim = 2;
//		int div = 0;
//		an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//		an.PostProcess(div,dim);
//		std::ofstream out("malha.txt");
//		an.Print("nothing",out);
//	}
//#endif	
}

void PosProcess(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
	vecnames[0]= "Derivate";
	
	const int dim = 2;
	int div = 0;
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
//			TPZGeoEl *geoel = compEl->Reference();
			//int ns = geoel->NSides();
//			TPZGeoElSide *geoside = new TPZGeoElSide(geoel,ns-1);
//			int subel = geoside->NSubElements();
			//if((geoel->Dimension()==2)/* && subel == 0*/){
				compEl->Divide(ind, subindex, 0);
			//}
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

void GeoElMultiphysicVec(TPZManVector<TPZCompMesh *> cmeshVec, std::set<int> &geoelVec){
	
	if(cmeshVec.NElements() == 0) return;
	TPZCompMesh *cmesh = cmeshVec[0];
	TPZGeoMesh *gmesh = cmesh->Reference();
	gmesh->ResetReference();
	int isub;
	int ncm = cmeshVec.NElements();
	for (isub=0; isub<ncm; isub++) {
		cmeshVec[isub]->LoadReferences();
	}
	int ncel;
	TPZStack<TPZCompElSide> sidevec;
	for(int i = 0; i< ncm; i++){
		ncel = cmeshVec[i]->NElements();
		for (int j=0; j<ncel; j++) {
			TPZCompEl * cel = cmeshVec[i]->ElementVec()[j];
			if(cel){
				TPZGeoEl *geoel = cel->Reference();
				if (!geoel) {
					std::cout << "Geoel nulo!\n";
					DebugStop();
				}
				int ns = geoel->NSides();
				TPZGeoElSide *geoside = new TPZGeoElSide(geoel,ns-1);
				sidevec.Resize(0);
				geoside->HigherLevelCompElementList2(sidevec, 1,1);
				int nel = sidevec.NElements();
				if (nel==0){
					//std::cout << "Incluindo elemento " << geoel->Index() << std::endl;
					geoelVec.insert(geoel->Index());
				}
			}
		}
	}
		//cout<<"num Elemento geom= "<<geoelVec.size()<<endl;
	//set<int>::iterator it;
//	cout << "myset contains:";
//	for (it=geoelVec.begin() ; it != geoelVec.end(); it++ )
//		cout << " " << *it;
//	cout << endl;
	
}

void AddElements(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh)
{
	TPZGeoMesh *gmesh = MFMesh->Reference();
	gmesh->ResetReference();
	int nMFEl = MFMesh->NElements();
	int nmesh = cmeshVec.size();
	int imesh;
	for(imesh = 0; imesh<nmesh; imesh++)
	{
		cmeshVec[imesh]->LoadReferences();
		int iel, is;
		for(iel=0; iel<nMFEl; iel++)
		{
			TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *> (MFMesh->ElementVec()[iel]);
			if(!mfcel)
			{
				DebugStop();
			}
			TPZGeoEl *gel = mfcel->Reference();
			TPZStack<TPZCompElSide> celstack;
			TPZGeoElSide gelside(gel,gel->NSides()-1);
			if (gel->Reference()) {
				celstack.Push(gelside.Reference());
			}
			gelside.ConnectedCompElementList(celstack, 0, 0);
			for(is=0;is<celstack.size();is++) {
				TPZGeoElSide gelside = celstack[is].Reference();
				
				if(gelside.Element()->Dimension()!=gel->Dimension()) {
					if(is!=celstack.size()-1) {
						celstack[is] = celstack.Pop();
						is--;
					}
					else {
						celstack.Pop();
					}

				}
			}
			if(celstack.size() != 1)
			{
				DebugStop();
			}
			
			mfcel->AddElement(celstack[0].Element(), imesh);
			
			//TPZManVector<TPZTransform> tr;
			//mfcel->AffineTransform(tr);
//			
//#ifdef LOG4CXX
//			{
//				int itr = tr.size();
//				std::stringstream sout;
//				for (int i = 0; i< itr; i++) {
//					sout << "Transformacao para referencia " << i << std::endl;
//					sout << tr[i] << std::endl;
//					
//				}
//				LOGPZ_DEBUG(logger, sout.str())
//			}
//#endif			
		}
		gmesh->ResetReference();
	}
		
	
}

void AddConnects(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh)
{
	int nmeshes = cmeshVec.size();
	TPZVec<int> FirstConnect(nmeshes,0);
	int nconnects = 0;
	int imesh;
	for (imesh=0; imesh<nmeshes; imesh++) 
	{
		FirstConnect[imesh] = nconnects;
		nconnects += cmeshVec[imesh]->ConnectVec().NElements();
	}
	MFMesh->ConnectVec().Resize(nconnects);
	MFMesh->Block().SetNBlocks(nconnects);
	int counter = 0;
	int seqnum = 0;
	for (imesh=0; imesh<nmeshes; imesh++) 
	{
		int ic;
		int nc = cmeshVec[imesh]->ConnectVec().NElements();
		for (ic=0; ic<nc; ic++) 
		{
			TPZConnect &refcon =  cmeshVec[imesh]->ConnectVec()[ic];
			MFMesh->ConnectVec()[counter] = refcon;
			if (refcon.SequenceNumber() >= 0) {
				MFMesh->ConnectVec()[counter].SetSequenceNumber(seqnum);
                MFMesh->ConnectVec()[counter].SetNShape(refcon.NShape());
                MFMesh->ConnectVec()[counter].SetNState(refcon.NState());
				int ndof = refcon.NDof(*cmeshVec[imesh]);
				MFMesh->Block().Set(seqnum,ndof);
				seqnum++;
			}
			counter++;
		}	
		// ajustar as dependencias
		for (ic=0; ic<nc; ic++) 
		{
			TPZConnect &cn = MFMesh->ConnectVec()[FirstConnect[imesh]+ic];
			if (cn.HasDependency()) 
			{
				TPZConnect::TPZDepend *dep = cn.FirstDepend();
				while (dep) {
					dep->fDepConnectIndex = dep->fDepConnectIndex+FirstConnect[imesh];
					dep = dep->fNext;
				}
			}
		}	
	}
	MFMesh->Block().SetNBlocks(seqnum);
	MFMesh->ExpandSolution();
	int iel;
	int nelem = MFMesh->NElements();
	for (iel = 0; iel < nelem; iel++) 
	{
		TPZMultiphysicsElement *cel = dynamic_cast<TPZMultiphysicsElement *> (MFMesh->ElementVec()[iel]);
		if (!cel) {
			DebugStop();
		}
		TPZStack<int> connectindexes;
		int imesh;
		for (imesh=0; imesh < nmeshes; imesh++) {
			TPZCompEl *celref = cel->ReferredElement(imesh);
			int ncon = celref->NConnects();
			int ic;
			for (ic=0; ic<ncon; ic++) {
				connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[imesh]);
			}
		}
		cel->SetConnectIndexes(connectindexes);
	}
}

void TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
    int imesh;
    int nmeshes = cmeshVec.size();
    TPZManVector<int> FirstConnectIndex(nmeshes+1,0);
    for (imesh = 0; imesh < nmeshes; imesh++) {
        FirstConnectIndex[imesh+1] = FirstConnectIndex[imesh]+cmeshVec[imesh]->NConnects();
    }
    TPZBlock &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
        int ncon = cmeshVec[imesh]->NConnects();
        TPZBlock &block = cmeshVec[imesh]->Block();
        int ic;
        for (ic=0; ic<ncon; ic++) {
            TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
            int seqnum = con.SequenceNumber();
			if(seqnum<0) continue;       // Whether connect was deleted by previous refined process
            int blsize = block.Size(seqnum);
            TPZConnect &conMF = MFMesh->ConnectVec()[FirstConnectIndex[imesh]+ic];
            int seqnumMF = conMF.SequenceNumber();
            int idf;
            for (idf=0; idf<blsize; idf++) {
                blockMF.Put(seqnumMF, idf, 0, block.Get(seqnum, idf, 0));
            }
        }
    }
}

void TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
    int imesh;
    int nmeshes = cmeshVec.size();
    TPZManVector<int> FirstConnectIndex(nmeshes+1,0);
    for (imesh = 0; imesh < nmeshes; imesh++) {
        FirstConnectIndex[imesh+1] = FirstConnectIndex[imesh]+cmeshVec[imesh]->NConnects();
    }
    TPZBlock &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
        int ncon = cmeshVec[imesh]->NConnects();
        TPZBlock &block = cmeshVec[imesh]->Block();
        int ic;
        for (ic=0; ic<ncon; ic++) {
            TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
            int seqnum = con.SequenceNumber();
			if(seqnum<0) continue;       // Whether connect was deleted by previous refined process
            int blsize = block.Size(seqnum);
            TPZConnect &conMF = MFMesh->ConnectVec()[FirstConnectIndex[imesh]+ic];
            int seqnumMF = conMF.SequenceNumber();
            int idf;
            for (idf=0; idf<blsize; idf++) {
                block.Put(seqnum, idf, 0, blockMF.Get(seqnumMF, idf, 0));
            }
        }
    }
    
}

void BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, int LagrangeMat, int InterfaceMat)
{
	TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh->Reference()->ElementVec();
    int meshdim = cmesh->Dimension();
    
    // cria todos os elementos sem conectar-se aos vizinhos
	int i, nelem = elvec.NElements();
	int neltocreate = 0;
	int index;
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			neltocreate++;
		}
	}
	std::set<int> matnotfound;
	int nbl = cmesh->Block().NBlocks();
	if(neltocreate > nbl) cmesh->Block().SetNBlocks(neltocreate);
	cmesh->Block().SetNBlocks(nbl);
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			int matid = gel->MaterialId();
			TPZAutoPointer<TPZMaterial> mat = cmesh->FindMaterial(matid);
			if(!mat)
			{
				matnotfound.insert(matid);
				continue;
			}
			int printing = /*0*/1;
			if (printing) {
				gel->Print(cout);
			}
			
			///checking material in MaterialIDs
            std::set<int>::const_iterator found = MaterialIDs.find(matid);
            if (found == MaterialIDs.end()) continue;
			
			if(!gel->Reference() && gel->NumInterfaces() == 0)
			{
				cmesh->CreateCompEl(gel,index);
                gel->ResetReference();
			}
		}
	}
    
    cmesh->LoadReferences();
    
    // Gera elementos geometricos para os elementos menores
    for (i=0; i<nelem; ++i) {
        TPZGeoEl *gel = elvec[i];
        if (!gel || gel->Dimension() != meshdim || !gel->Reference()) {
            continue;
        }
        int matid = gel->MaterialId();
        if(MaterialIDs.find(matid) == MaterialIDs.end())
        {
            continue;
        }
        // over the dimension-1 sides
        int nsides = gel->NSides();
        int is;
        for (is=0; is<nsides; ++is) {
            int sidedim = gel->SideDimension(is);
            if (sidedim != meshdim-1) {
                continue;
            }
			// check if there is a smaller element connected to this element
			TPZStack<TPZCompElSide> celsides;
			celsides.Resize(0);
			TPZGeoElSide gelside(gel,is);
			gelside.HigherLevelCompElementList2(celsides, 0, 0);
			int ncelsid =  celsides.NElements();
			if(celsides.NElements()) continue;
			
			//check the neighboring
			TPZCompElSide celside;
			celside = gelside.LowerLevelCompElementList2(0);
			if (celside && celside.Element()->Reference()->Dimension() != meshdim) continue;
			TPZStack<TPZGeoElSide> allneigh;
			allneigh.Resize(0);	
			gelside.AllNeighbours(allneigh);
			int nneig = allneigh.NElements();
			if(allneigh.NElements()>1) continue;
			if (nneig && allneigh[0].Element()->Dimension() != meshdim) continue;
						
			gel->CreateBCGeoEl(is, LagrangeMat);
		}
	}
   
	// now create the lagrange elements
    cmesh->Reference()->ResetReference();
    nelem = elvec.NElements();
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			int matid = gel->MaterialId();
			TPZAutoPointer<TPZMaterial> mat = cmesh->FindMaterial(matid);
			
			if(!mat)
			{
				matnotfound.insert(matid);
				continue;
			}
		
			///checking material in MaterialIDs
			if (matid != LagrangeMat) {
				continue;
			}
			int printing = /*0*/1;
			if (printing) {
				gel->Print(cout);
			}
			
			
			if(!gel->Reference())
			{
				cmesh->CreateCompEl(gel,index);
                gel->ResetReference();
			}
		}
	}
	
	//Agnaldo
	cmesh->LoadReferences();
	
    // now create the interface elements between the lagrange elements and other elements
    nelem = elvec.NElements();
    for (i=0; i<nelem; ++i) {
        TPZGeoEl *gel = elvec[i];
        if (!gel || gel->Dimension() != meshdim-1 || !gel->Reference()) {
            continue;
        }
        int matid = gel->MaterialId();
        if(matid != LagrangeMat)
        {
            continue;
        }
        // over the dimension-1 sides
        int nsides = gel->NSides();
        int is;
        for (is=0; is<nsides; ++is) {
            int sidedim = gel->SideDimension(is);
            if (sidedim != meshdim-1) {
                continue;
            }
            // check if there is a smaller element connected to this element
            TPZStack<TPZCompElSide> celsides;
            TPZGeoElSide gelside(gel,is);
            gelside.EqualLevelCompElementList(celsides, 0, 0);
            if(celsides.size() < 1)
            {
                DebugStop();
            }
            TPZCompElSide cels = gelside.LowerLevelCompElementList2(0);
            if(cels) 
            {
                celsides.Push(cels);
            }
            int nelsides = celsides.NElements();
            if(nelsides != 2) 
            {
                DebugStop();
            } 
            for (int lp=0; lp<nelsides; ++lp) {
                TPZGeoEl *interface = gel->CreateBCGeoEl(is, InterfaceMat);
                TPZCompElSide right = celsides[lp];
                TPZCompElSide left(gel->Reference(),is);
                int index;
                new TPZInterfaceElement(*cmesh,interface,index,left,right);
            }
        }
    }
	
	cmesh->InitializeBlock();

}

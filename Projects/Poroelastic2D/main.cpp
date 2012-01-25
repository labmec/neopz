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

#include "pzelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzporoelastic2d.h"

#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.poroelastic2d"));
#endif

using namespace std;

const int matId = 1;
const int dirichlet = 0;
const int neumann = 1;

const int bcNL = -1;
const int bcDR = -2;
const int bcNU = -3;
const int bcDL = -4;


TPZGeoMesh *MalhaGeom(REAL h,REAL L);
TPZCompMesh *MalhaCompPressao(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh *MalhaCompElast(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);

void PosProcess(TPZAnalysis &an, std::string plotfile);
void PosProcess2(TPZAnalysis &an, std::string plotfile);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void RefinamentoUniforme(TPZGeoMesh  *gMesh, int nh);
void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh, int MatId, int indexEl);
void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);
void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file);



int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	InitializePZLOG("../mylog.cfg");
#endif
	
	int p=1;
	//primeira malha
	
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = MalhaGeom(2.,2.);
	ofstream arg1("gmesh1.txt");
	gmesh->Print(arg1);
	ofstream file1("malhageoInicial.vtk");
	PrintGMeshVTK(gmesh, file1);
	
	// First computational mesh
	TPZCompMesh * cmesh1 = MalhaCompElast(gmesh, p);
	ofstream arg2("cmesh1.txt");
	cmesh1->Print(arg2);
	
	// Second computational mesh
	TPZCompMesh * cmesh2= MalhaCompPressao(gmesh, p);
	ofstream arg3("cmesh2.txt");
	cmesh2->Print(arg3);
	
	// Cleaning reference of the geometric mesh to cmesh1
	gmesh->ResetReference();
	cmesh1->LoadReferences();
	//RefinUniformElemComp(cmesh1,2);
	RefinElemComp(cmesh1,1);
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();
	
	ofstream arg4("cmesh12.txt");
	cmesh1->Print(arg4);
	ofstream arg5("gmesh2.txt");
	gmesh->Print(arg5);
	ofstream file3("malhageo1.vtk");
	PrintGMeshVTK(gmesh, file3);
	
	// Cleaning reference to cmesh1
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	
	//refinamento uniform
	RefinUniformElemComp(cmesh2,1);
	//RefinElemComp(cmesh2,4);
	//RefinElemComp(cmesh2,7);
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
	PosProcess2(an1, plotfile);
		
	//--- Resolver usando a segunda malha computacional ---
	TPZAnalysis an2(cmesh2);
	SolveSist(an2, cmesh2);
	std::string plotfile2("saidaSolution_cmesh2.vtk");
	PosProcess(an2, plotfile2);
		
	//--- Resolver usando a malha computacional multifisica ---
	TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh1;
	meshvec[1] = cmesh2;
	TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec);
	
	ofstream file6("mphysics.vtk");
	PrintGMeshVTK(gmesh, file6);
	
	TPZAnalysis an(mphysics);
	SolveSist(an, mphysics);
	std::string plotfile3("saidaMultphysics.vtk");
	PosProcessMultphysics(meshvec,mphysics, an, plotfile3);

	return EXIT_SUCCESS;
}


TPZGeoMesh *MalhaGeom(REAL h, REAL L)
{
	
	int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolLine(2);
	
	//indice dos nos
	int id = 0;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,xi*L);//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int yi = 0; yi < Qnodes/2; yi++)
	{
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,(1-yi)*L );//coord X
		Node[id].SetCoord(1 ,h );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcNL,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcDR,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcNU,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcDL,*gmesh);
	id++;
		
	TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 2;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
		
	gmesh->BuildConnectivity();
	
	//ofstream arg("gmesh.txt");
	//	gmesh->Print(arg);
	
	return gmesh;
	
}

TPZCompMesh*MalhaCompPressao(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	TPZAutoPointer<TPZMaterial> mat(material);
	
	REAL diff = 0.1;
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
	REAL pN=0.;
	val2(0,0)=pN;
	TPZAutoPointer<TPZMaterial> BCondNL = material->CreateBC(mat, bcNL,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondNL);
	
	TPZAutoPointer<TPZMaterial> BCondNU = material->CreateBC(mat, bcNU,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondNU);
	
	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	REAL uDL=4000.;
	val22(0,0)=uDL;
	TPZAutoPointer<TPZMaterial> BCondDL = material->CreateBC(mat, bcDL,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDL);
	
	REAL uDR=3000.;
	val22(0,0)=uDR;
	TPZAutoPointer<TPZMaterial> BCondDR = material->CreateBC(mat, bcDR,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDR);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}

TPZCompMesh*MalhaCompElast(TPZGeoMesh * gmesh,int pOrder)
{
	REAL rockrho = 2330.0; // SI system
	REAL gravity = 9.8; // SI system
	REAL overburdendepth = 2000.0; // SI system
	REAL layerthickness = 10.0;  // SI system
	
	REAL E = 100.;
	REAL poisson = 0.35;
	int dim = 2;
	TPZVec<REAL> force(2,0.);
	force[1]=gravity*rockrho;
	//force[0] = 0.;
	TPZElasticityMaterial *material;
	int planestress = 1;
	material = new TPZElasticityMaterial(matId, E, poisson, force[0], force[1], planestress); 
	
	
	TPZAutoPointer<TPZMaterial> mat(material);
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	
	
	///Inserir condicao de contorno
	
	TPZFMatrix val1(2,2,0.), val2(2,1,0.);
	REAL uNUy=rockrho*gravity*overburdendepth;
	val2(1,0)=uNUy;
	
	TPZAutoPointer<TPZMaterial> BCondNU = material->CreateBC(mat, bcNU,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondNU);
	
	REAL uNLy=rockrho*gravity*(overburdendepth+layerthickness);
	val2(1,0)=uNLy;	
	
	TPZAutoPointer<TPZMaterial> BCondNL = material->CreateBC(mat, bcNL,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondNL);
	
	
	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	REAL uDL=0.0;
	val22(0,0)=uDL;
	TPZAutoPointer<TPZMaterial> BCondDL = material->CreateBC(mat, bcDL,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDL);
	
	REAL uDR=0.0;
	val22(0,0)=uDR;
	TPZAutoPointer<TPZMaterial> BCondDR = material->CreateBC(mat, bcDR,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondDR);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
	
}

TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
{
	//Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
	mphysics->SetAllCreateFunctionsMultiphysicElem();
	
	int MatId = 1;
	int dim = 2;
	REAL Eyoung = 100.;
	REAL nu = 0.35;
	REAL alpha=1.0;
	REAL biot=1.0;
	REAL bulk=2.0;	
	REAL rockrho = 2330.0; // SI system
	REAL gravity = 9.8; // SI system
	REAL fx=0.0;
	REAL fy=gravity*rockrho;
	REAL overburdendepth = 2000.0; // SI system
	REAL layerthickness = 10.0;  // SI system
	
	REAL perm = 1.;
	REAL visc = 10.;
	int planestress = 1;
	
	TPZPoroElastic2d  *mymaterial = new TPZPoroElastic2d (MatId, dim);
	
	mymaterial->SetParameters(Eyoung, nu, alpha, fx, fy);
	mymaterial->SetParameters(perm,visc);
	mymaterial->SetfPlaneProblem(planestress);
	
	// Biot's parameter @biot
	// Bulk rock parameter @bulk
	mymaterial->SetAlpha(biot,bulk);
	
	
	ofstream argm("mymaterial.txt");
	mymaterial->Print(argm);
	TPZAutoPointer<TPZMaterial> mat(mymaterial);
	mphysics->InsertMaterialObject(mat);
	
	///Inserir condicao de contorno de Neumann
	int neum = 11;
	TPZFMatrix val1(3,2,0.), val2(3,1,0.);
	REAL uNUy=rockrho*gravity*overburdendepth;
	REAL uNLy=rockrho*gravity*(overburdendepth+layerthickness);
	REAL uNUx=0.;
	REAL uNLx=0.;
	REAL pNU=0.;
	REAL pNL=0.;
	
	val2(0,0)=uNUx;
	val2(1,0)=uNUy;
	val2(2,0)=pNU;
	TPZAutoPointer<TPZMaterial> BCondNU = mymaterial->CreateBC(mat, bcNU,neum, val1, val2);
	mphysics->InsertMaterialObject(BCondNU);
	
	val2(0,0)=uNLx;
	val2(1,0)=uNLy;
	val2(2,0)=pNL;
	TPZAutoPointer<TPZMaterial> BCondNL = mymaterial->CreateBC(mat, bcNL,neum, val1, val2);
	mphysics->InsertMaterialObject(BCondNL);
	
	///Inserir condicao de contorno de Dirichlet
	int dirich =0;
	TPZFMatrix val12(3,2,0.), val22(3,1,0.);
	REAL uDLeftx = 0.;
	REAL uDRx = 0.;
	REAL uDLefty = 0.;
	REAL uDRy = 0.;
	REAL pDLeft=4000.;
	REAL pDR=3000.;
	
	val22(0,0)=uDLeftx;
	val22(1,0)=uDLefty;
	val22(2,0)=pDLeft;
	TPZAutoPointer<TPZMaterial> BCondDL = mymaterial->CreateBC(mat, bcDL, dirich, val12, val22);
	mphysics->InsertMaterialObject(BCondDL);
	
	val22(0,0)=uDRx;
	val22(1,0)=uDRy;
	val22(2,0)=pDR;
	TPZAutoPointer<TPZMaterial> BCondDR = mymaterial->CreateBC(mat, bcDR, dirich, val12, val22);
	mphysics->InsertMaterialObject(BCondDR);
	//-----------
	
	mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
	
	// Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh * Objectdumy;
	Objectdumy->AddElements(meshvec, mphysics);
	Objectdumy->AddConnects(meshvec,mphysics);
	Objectdumy->TransferFromMeshes(meshvec, mphysics);
	
#ifdef LOG4CXX
    {
        std::stringstream sout;
        mphysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	ofstream arg("mphysic.txt");
	mphysics->Print(arg);
	
	return mphysics;
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
	
	ofstream file("Solution.out");
	an.Solution().Print("solution", file); 
}

void PosProcess(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
	vecnames[0]= "MinusKGradU";
		
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void PosProcess2(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(3), vecnames(1);
	scalnames[0] = "SigmaX";
	scalnames[1] = "SigmaY";
	scalnames[2] = "Pressure";	
	vecnames[0]= "displacement";
	
	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
{
	TPZBuildMultiphysicsMesh * Objectdumy;
	Objectdumy->TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(4), vecnames(2);
	scalnames[0] = "SigmaX";
	scalnames[1] = "SigmaY";
	scalnames[2] = "Pressure";
	scalnames[3] = "SolutionP";
	vecnames[0]= "Desplacement";
	vecnames[1]= "MinusKGradP";
			
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
			compEl->Divide(ind, subindex, 0);
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


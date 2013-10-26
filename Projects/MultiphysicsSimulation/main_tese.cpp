//
//  main_tese.cpp
//  PZ
//
//  Created by Agnaldo Farias on 2/12/13.
//
//

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
const int neumann_dirichlet = 10;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;


TPZGeoMesh *MalhaGeom();

TPZCompMesh *MalhaCompUm(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh *MalhaCompDois(TPZGeoMesh * gmesh, int pOrder);
TPZCompMesh *MalhaCompMultifisica(TPZGeoMesh * gmesh,/* TPZVec<TPZCompMesh *> meshvec,*/ TPZMatPoissonDesacoplado *&mymaterial);

void PrintGMeshVTK2(TPZGeoMesh * gmesh, std::ofstream &file);

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh);
void SaidaSolucao(TPZAnalysis &an, std::string plotfile);

void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

int main2(int argc, char *argv[])
{
    // Ordem polinomial das funções de aproximação
	int p = 2;
	gRefDBase.InitializeAllUniformRefPatterns();
    
	/*------------ Etapa 1 ------------*/
    
	// Criando a malha geométrica
	TPZGeoMesh * gmesh = MalhaGeom();
	
	// Criando a primeira malha computacional
	TPZCompMesh * cmesh1= MalhaCompUm(gmesh, p);
    
	// Criando a segunda malha computacional
	TPZCompMesh * cmesh2 = MalhaCompDois(gmesh, p+1);
	
    /*------------ Etapa 2 ------------*/
	
    // Refinando a malha da primeira equação
    gmesh->ResetReference();
	cmesh1->LoadReferences();
    // Refinando a malha com dois níveis de refinamneto uniforme
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1, 2,false);
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();
    
	// Refinando a malha da segunda equação
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	// Refinando a malha com três níveis de refinamneto uniforme
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2, 3,false);
	cmesh2->AdjustBoundaryElements();
	cmesh2->CleanUpUnconnectedNodes();
    
    /*---------- Etapa 3 ----------*/
    
    // Limpando a referência da malha geométrica
    gmesh->ResetReference();
    
	// Criando a malha computacional multifísica
    TPZMatPoissonDesacoplado * multiphysics_material;
    TPZCompMesh * mphysics = MalhaCompMultifisica(gmesh,multiphysics_material);
    
    // Criando um vetor de malhas computacionais
	TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh1;
	meshvec[1] = cmesh2;
    
    // Adicionando os elementos na malha multifísica
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
  
    /*---------- Etapa 4 ----------*/
    
	// Resolvendo o sistema linear
	TPZAnalysis an(mphysics);
	ResolverSistema(an, mphysics);
	
    /*---------- Etapa 5 ----------*/

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    // Arquivo de saida para plotar a solução
    string plotfile("Solution_mphysics.vtk");
    SaidaSolucao(an, plotfile);
    
	return EXIT_SUCCESS;
}



TPZGeoMesh *MalhaGeom()
{
	
	int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <long> TopolQuad(4);
	TPZVec <long> TopolLine(2);
	
	//indice dos nos
	long id = 0;
	REAL valx;
	REAL dx=1.;
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
    
    TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 2;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
	id++;
    
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
		
	gmesh->BuildConnectivity();
    
	return gmesh;
	
}

TPZCompMesh *MalhaCompUm(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
	TPZMaterial * mat(material);
	
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
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}

TPZCompMesh *MalhaCompDois(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
	TPZMaterial * mat(material);
	
	REAL diff = -1.;
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
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	
	return cmesh;
}

TPZCompMesh *MalhaCompMultifisica(TPZGeoMesh * gmesh,/* TPZVec<TPZCompMesh *> meshvec,*/ TPZMatPoissonDesacoplado* &mymaterial){
    
    
    // Creating computational mesh for multiphysic elements
//	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
	mphysics->SetAllCreateFunctionsMultiphysicElem();
	
    int dim = 2;
    
    mymaterial = new TPZMatPoissonDesacoplado(matId, dim);
	
    mymaterial->SetParameters(-1., -1.);
    mymaterial->SetInternalFlux(8.,0.);
    
	TPZMaterial * mat(mymaterial);
	mphysics->InsertMaterialObject(mat);
    
    
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(ForcingF);
    mymaterial->SetForcingFunction(forcef);

	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
	TPZMaterial * BCond0 = mymaterial->CreateBC(mat, bc0,neumann_dirichlet, val1, val2);
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat, bc2,neumann_dirichlet, val1, val2);
    TPZMaterial * BCond1 = mymaterial->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    
	mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);

	mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
	
	// Creating multiphysic elements into mphysics computational mesh
//    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
    return mphysics;
    
}


#define VTK
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
    
	//Saida de Dados: solucao e  grafico no VTK
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
	
}

void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}

void SaidaSolucao(TPZAnalysis &an, std::string plotfile){
    
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
	scalnames[0] = "SolutionU";
	scalnames[1] = "SolutionP";
	vecnames[0]= "DerivateU";
	vecnames[1]= "DerivateP";

	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void PrintGMeshVTK2(TPZGeoMesh * gmesh, std::ofstream &file)
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


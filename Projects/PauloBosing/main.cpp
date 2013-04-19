
#include <time.h>
#include <stdio.h>
#include <fstream>
#include <cmath>

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
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzuncoupledpoissondisc.h"


#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

using namespace std;

const int matId = 1;
const int dirichlet = 0;
const int neumann = 1;
const int neumann_dirichlet = 11;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;

TPZGeoMesh *MalhaGeom();

TPZCompMesh *MalhaCompUm(TPZGeoMesh * gmesh,int pOrder, bool isdiscontinuous);
TPZCompMesh *MalhaCompDois(TPZGeoMesh * gmesh, int pOrder, bool isdiscontinuous);
TPZCompMesh *MalhaCompMultifisica(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMatUncoupledPoissonDisc *&mymaterial);

void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix);
void SaidaSolucao(TPZAnalysis &an, std::string plotfile);
void SaidaSolucaoMultifisica(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh *mphysics, TPZAnalysis &an, std::string plotfile);
void CreatInterface(TPZCompMesh *cmesh);
void ChecarIterface(TPZCompMesh *mphysics);

void ForcingF(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

bool disc_functions = false;

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout<<"\n Elemento computacional " << std::endl;
        LOGPZ_DEBUG(logdata,sout.str())
    }
#endif
    
    
    // Ordem polinomial das funções de aproximação
	int p = 2;
    
	//---- Criando a malha geométrica ----
	TPZGeoMesh * gmesh = MalhaGeom();
    ofstream arg1("gmesh_inicial.txt");
	gmesh->Print(arg1);
	
	//---- Criando a primeira malha computacional -----
	TPZCompMesh * cmesh1= MalhaCompUm(gmesh, p,disc_functions);

	//----- Criando a segunda malha computacional ------
	TPZCompMesh * cmesh2 = MalhaCompDois(gmesh, p+1,disc_functions);

    
// -------Refinando as malhas de cada equacao-------
	
    // Refinando a malha da primeira equação
    gmesh->ResetReference();
	cmesh1->LoadReferences();
    // Refinando a malha com dois níveis de refinamneto uniforme
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,2);
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();

    ofstream arg4("cmesh_edp1_final.txt");
	cmesh1->Print(arg4);

	// Refinando a malha da segunda equação
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	// Refinando a malha com três níveis de refinamneto uniforme
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,3);
	cmesh2->AdjustBoundaryElements();
	cmesh2->CleanUpUnconnectedNodes();

    ofstream arg6("cmesh_edp2_final.txt");
	cmesh2->Print(arg6);
    
    //--- Criando a malha computacional multifísica ----
    
    // Criando um vetor de malhas computacionais
	TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh1;
	meshvec[1] = cmesh2;
    
    // Criando a malha computacional multifisica
    TPZMatUncoupledPoissonDisc * multiphysics_material;
    TPZCompMesh * mphysics = MalhaCompMultifisica(gmesh,meshvec,multiphysics_material);
    
    ofstream arg13("gmesh_multiphysics.txt");
	gmesh->Print(arg13);
        
	// Resolvendo o sistema linear
	TPZAnalysis an(mphysics);
	ResolverSistema(an, mphysics,false);
    
    ofstream arg18("mphysics_cmesh.txt");
	mphysics->Print(arg18);

    // Arquivo de saida para plotar a solução
    string plotfile3("Solution_mphysics.vtk");
    SaidaSolucaoMultifisica(meshvec, mphysics, an, plotfile3);
    
	return EXIT_SUCCESS;
}


TPZGeoMesh *MalhaGeom()
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

TPZCompMesh *MalhaCompUm(TPZGeoMesh * gmesh, int pOrder, bool isdiscontinuous)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
	TPZMaterial * mat(material);
    
    material->SetNoPenalty();
    material->SetNonSymmetric();
	
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
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
	
    //Ajuste da estrutura de dados computacional
    if (isdiscontinuous==true) {
        //Set discontinuous functions
        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->AutoBuild();
        cmesh->ExpandSolution();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
    }
    else{
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();
        cmesh->ExpandSolution();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    }
	return cmesh;
}

TPZCompMesh *MalhaCompDois(TPZGeoMesh * gmesh, int pOrder, bool isdiscontinuous)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
	TPZMaterial * mat(material);
    
    material->SetNoPenalty();
    material->SetNonSymmetric();
	
	REAL diff = -1.;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
	REAL flux = 0.;
	
	material->SetParameters(diff, conv, convdir);
	material->SetInternalFlux(flux);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(ForcingF);
    material->SetForcingFunction(forcef);
	
	///Inserir condicao de contorno
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	//Ajuste da estrutura de dados computacional
	if (isdiscontinuous==true) {
        //Set discontinuous functions
        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->AutoBuild();
        cmesh->ExpandSolution();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    }
    else{
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();
        cmesh->ExpandSolution();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    }
    
	return cmesh;
}

TPZCompMesh *MalhaCompMultifisica(TPZGeoMesh * gmesh,TPZVec<TPZCompMesh *> meshvec, TPZMatUncoupledPoissonDisc* &mymaterial){
    
    
    // Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
	mphysics->SetAllCreateFunctionsMultiphysicElem();
	
    int dim = 2;
    mphysics->SetDimModel(dim);
    mymaterial = new TPZMatUncoupledPoissonDisc(matId, mphysics->Dimension());
	
    mymaterial->SetParameters(-1., -1.);
    mymaterial->SetInternalFlux(8.,0.);
    
    mymaterial->SetNonSymmetricOne();
    mymaterial->SetNonSymmetricTwo();
    mymaterial->SetPenaltyConstant(0., 0.);
    
	TPZMaterial * mat(mymaterial);
	mphysics->InsertMaterialObject(mat);
    
    
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(ForcingF);
    mymaterial->SetForcingFunction(forcef);

	
	///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
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
	
	//Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    
    if (disc_functions==true){
        //criar elementos de interface
        int nel = mphysics->ElementVec().NElements();
        for(int el = 0; el < nel; el++)
        {
            TPZCompEl * compEl = mphysics->ElementVec()[el];
            if(!compEl) continue;
            int index = compEl ->Index();
            if(compEl->Dimension() == mphysics->Dimension())
            {
                TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(mphysics->ElementVec()[index]);
                if(!InterpEl) continue;
                InterpEl->CreateInterfaces();
                
            }
        }
    }
    
    return mphysics;
    
}

void ChecarIterface(TPZCompMesh *mphysics){
    int nel = mphysics->NElements();
    for(int i = 0; i< nel; i++){
        TPZCompEl *cel = mphysics->ElementVec()[i];

        if(!cel) continue;
        TPZMultiphysicsElement *mphel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        TPZMultiphysicsInterfaceElement *mpintel = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);

        if(mphel){
            int nsides = mphel->Reference()->NSides();
            std::cout <<"\n\n For the element of index = " << i;
            for(int is = 0; is<nsides; is++){
                bool isinterface;
                isinterface = mphel->ExistsInterface(is);
                std::cout <<"\n There is interface by side = " << is << " ? ==> " << isinterface;
            }
        }
        if(mpintel){

            TPZCompElSide leftel;
            TPZCompElSide rightel;
            mpintel->GetLeftRightElement(leftel,rightel);
            leftel.Element()->Print();
            rightel.Element()->Print();
        }
    }

}

void CreatInterface(TPZCompMesh *cmesh){
    
    for(int el = 0; el < cmesh->ElementVec().NElements(); el++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        int index = compEl ->Index();
        if(compEl->Dimension() == cmesh->Dimension())
        {
            TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces(false);
        }
    }
    
   // cmesh->AdjustBoundaryElements();
    //cmesh->CleanUpUnconnectedNodes();
}

#define VTK
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix)
{
    if(symmetric_matrix ==true){
        TPZSkylineStructMatrix skmat(fCmesh);
        an.SetStructuralMatrix(skmat);
        TPZStepSolver<REAL> direct;
        direct.SetDirect(ELDLt);
        an.SetSolver(direct);
    }
    else{
        TPZBandStructMatrix bdmat(fCmesh);
        an.SetStructuralMatrix(bdmat);
        TPZStepSolver<REAL> direct;
        direct.SetDirect(ELU);
        an.SetSolver(direct);
    }
	an.Run();
	
	//Saida de Dados: solucao txt
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);
}

void ForcingF(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}

void SaidaSolucaoMultifisica(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh *mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
	scalnames[0] = "Solution_u1";
	scalnames[1] = "Solution_u2";
	vecnames[0]= "Derivate_u1";
	vecnames[1]= "Derivate_u2";
    
	const int dim = mphysics->Dimension();
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void SaidaSolucao(TPZAnalysis &an, std::string plotfile){
    
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
	vecnames[0]= "Derivative";
    
	const int dim = an.Mesh()->Dimension();
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
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


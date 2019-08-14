#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZMatLaplacianLagrange.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"

#include "TPZDPGMeshControl.h"
#include "TPZMHMeshControl.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzbuildmultiphysicsmesh.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzgengrid.h"

#include "TPZMHMeshControl.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;


const int matfinermesh = 1;
const int matcoarsemesh = 1;
const int matskeletonmesh = 2;

const int dirichlet = 0;
const int neumann = 1;
const int mixed = 2;
const int neumanndirichlet = 10;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;


TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly,bool triang_elements);
void RefinamentoUniforme(TPZAutoPointer<TPZGeoMesh> gmesh, int nref,TPZVec<int> dims);
void InsertMaterialObjectsMHM(TPZCompMesh &cmesh, bool useDPGPhil, bool useDPG);
void InsertMaterialObjects(TPZCompMesh &cmesh);
void GetElIndexCoarseMesh(TPZAutoPointer<TPZGeoMesh>  gmesh, std::set<int64_t> &coarseindex);

void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
void ForceSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &force);
static void DirichletSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFMatrix<STATE> du(2,1);
    SolSuave(loc,result,du);
}

void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out);

///problema de Steklov
TPZGeoMesh *GMeshSteklov(bool triang_elements);
void RefinamentoSingular(TPZAutoPointer<TPZGeoMesh> gmesh,int nref);
void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
static void DirichletSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,result,du);
}


int mainDPG(int argc, char *argv[])
{
//    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    
    ofstream saidaerros("ErroMDP-mhm.txt");
    saidaerros << "\nCalculo do Erro\n";
    
    saidaerros.precision(16);
    for(int p=1; p<5; p++){
        saidaerros << "\n";
        saidaerros << "Ordens P: p_fine = "<<p+2 <<", p_coarse = " << p <<", p_interface = " << p-1 <<"\n";
        for(int h=0; h<5; h++){
            saidaerros << "\nRefinamento: h = "<< h <<"\n";
            
            TPZAutoPointer<TPZGeoMesh> gmesh = MalhaGeom(1.,1.,false);
//            ofstream arg0("gmesh0.txt");
//            gmesh->Print(arg0);
            
            //-------- construindo malha coarse ----------
            //1 refinamento uniforme
            TPZVec<int> dims(2,0);
            dims[0]=1; dims[1]=2;
            int nref = h;
            RefinamentoUniforme(gmesh, nref, dims);
        //    ofstream arg1("gmesh1.txt");
        //	gmesh->Print(arg1);
            
            //index dos elementos da malha coarse
            std::set<int64_t> coarseindex;
            GetElIndexCoarseMesh(gmesh, coarseindex);
            
//            std::set<int64_t>::iterator it;
//            for (it=coarseindex.begin(); it!=coarseindex.end(); ++it)
//                std::cout << ' ' << *it;
//            std::cout << "\n";
            
            TPZAutoPointer<TPZGeoMesh> gmesh2 = new TPZGeoMesh(gmesh);
            dims.Resize(1, 0);
            dims[0]=2;
            nref = 0;
            RefinamentoUniforme(gmesh2, nref, dims);
        //    ofstream arg2("gmesh2.txt");
        //	gmesh2->Print(arg2);
            //--------------------------------------------------
            
            TPZDPGMeshControl dpgmesh(gmesh2,coarseindex);
            int porder = p;
            dpgmesh.SetMatIds(matfinermesh, matcoarsemesh, matskeletonmesh);
            dpgmesh.SetPOrderMeshes(porder+2, porder, porder-1);
            TPZCompMesh &corsemesh = dpgmesh.PressureCoarseMesh();
            InsertMaterialObjects(corsemesh);
            
            TPZMHMeshControl &mhm = dpgmesh.MHMControl();
            bool useDPGPhil=false;
            bool useMDP=true;
            InsertMaterialObjectsMHM(mhm.CMesh(),useDPGPhil,useMDP);
            
            //refinar malha fina
        //    mhm.SetLagrangeAveragePressure(true);
        //    mhm.CreateCoarseInterfaces(matskeletonmesh);
        //    mhm.BuildComputationalMesh();
        //    TPZVec<TPZCompMesh *> meshvec;
        //    mhm.GetMeshVec(meshvec);
        //    TPZCompMesh *finemesh = meshvec[0];
        //    finemesh->Reference()->ResetReference();
        //	finemesh->LoadReferences();
        //    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(finemesh,1,false);
        //    ofstream arg3("gmesh3.txt");
        //	gmesh2->Print(arg3);
        //    ofstream arg4("cmeshfine.txt");
        //    finemesh->Print(arg4);
            
            dpgmesh.BuildComputationalMesh();
            
            TPZAnalysis an(mhm.CMesh(),false);
            TPZStepSolver<STATE> step;
            if(useMDP) {
                TPZSkylineNSymStructMatrix bst(mhm.CMesh().operator->());
                //TPZBandStructMatrix bst(mhm.CMesh().operator->());
                an.SetStructuralMatrix(bst);
                bst.SetNumThreads(8);
                step.SetDirect(ELU);
            }else{
                TPZSkylineStructMatrix skyl(mhm.CMesh());
                an.SetStructuralMatrix(skyl);
                skyl.SetNumThreads(8);
                step.SetDirect(ELDLt);
            }
            an.SetSolver(step);
            an.Assemble();
            //an.Run();
            an.Solve();

//            std::string plotfile("result.vtk");
//            TPZStack<std::string> scalnames,vecnames;
//            scalnames.Push("Solution");
//            scalnames.Push("ExactSolution");
//            scalnames.Push("ErrorEstimatorDPG");
//            an.DefineGraphMesh(mhm.CMesh()->Dimension(), scalnames, vecnames, plotfile);
//            an.PostProcess(0,2);

            
            //calculo do erro
            TPZVec<TPZCompMesh *> meshvec;
            dpgmesh.GetMeshVec(meshvec);
            int NEq1 = mhm.CMesh()->NEquations();
            int NEq2 =meshvec[3]->NEquations();
            saidaerros << "\nNumero Equacoes: Malha Multifisica = "<< NEq1 <<"  e Malha coarse = "<< NEq2 <<"\n";
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mhm.CMesh().operator->());
            
            TPZVec<REAL> erros(3);
            TPZAnalysis an_coarse(meshvec[3],false);
            an_coarse.SetExact(*SolSuave);
            an_coarse.PostProcessError(erros,saidaerros);
            //ErrorH1(meshvec[3], saidaerros);
            
            mhm.CMesh().operator->()->CleanUp();
            meshvec[0]->CleanUp();
            meshvec[1]->CleanUp();
            meshvec[2]->CleanUp();
            meshvec[3]->CleanUp();
            gmesh.operator->()->CleanUp();
            gmesh2.operator->()->CleanUp();
        }
    }
    
	return EXIT_SUCCESS;
}


//STANDARD FEM
TPZCompMesh *CompMesh(TPZGeoMesh *gmesh, int porder);
int mainfem(int argc, char *argv[])
{
    //InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    ofstream saidaerros("ErroProjecaoSemiH1.txt");
    saidaerros << "\nCalculo do Erro\n";
    int h= 0;
    int p=0;
    
    saidaerros.precision(16);
    
    for(p=1; p<5; p++)
    {
        saidaerros << "\n";
        saidaerros << "Ordens p = " << p <<"\n";
        for(h=0; h<6; h++)
        {
            saidaerros << "\nRefinamento: h = "<< h <<"\n";
            
            //TPZGeoMesh * gmesh = MalhaGeom(1.,1.,false);
            TPZAutoPointer<TPZGeoMesh> gmesh = MalhaGeom(1.,1.,false);
            TPZVec<int> dims(2,0);
            dims[0]=1; dims[1]=2;
            int nref = h;
            RefinamentoUniforme(gmesh, nref, dims);
            
            
            TPZCompMesh * cmesh = CompMesh(gmesh.operator->(),p);
            
            //analysis
            TPZAnalysis an(cmesh,false);
            TPZSkylineStructMatrix skyl(cmesh);
            an.SetStructuralMatrix(skyl);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            an.SetSolver(step);
            an.Assemble();
            an.Solve();
            
//            std::string plotfile("result.vtk");
//            TPZStack<std::string> scalnames,vecnames;
//            scalnames.Push("Solution");
//            scalnames.Push("ExactPressure");
//            an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile);
//            an.PostProcess(0,2);
            
 //           TPZVec<REAL> erros(3);
//            an.SetExact(*SolSuave);
//            an.PostProcessError(erros,saidaerros);
            ErrorH1(cmesh, saidaerros);
            
            cmesh->CleanUp();
            gmesh->CleanUp();
            delete cmesh;
        }
    }
    
    return EXIT_SUCCESS;
}

#include "TPZMDPMaterial.h"
TPZCompMesh *MalhaMDP(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh);
TPZCompMesh *CompMesh1D(TPZGeoMesh *gmesh, int porder);
TPZCompMesh *MalhaMDP1D(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh);
TPZGeoMesh *MalhaGeom1D(REAL Lx, int h);
void Sol1D(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
void Neumann1D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);

REAL fc=0.;
REAL fk=0.;
int main(int argc, char *argv[])
{
    
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    ofstream saidaerros("ErroMDP.txt");
    saidaerros << "\nCalculo do Erro\n";
    int h= 0;
    int p=0, pr;
    
    int dim = 1;
    fc = 0.1;
    fk = 0.0001;
    saidaerros.precision(8);
    
    for(p=1; p<2; p++)
    {
        pr = p+3;
        saidaerros << "\n";
        saidaerros << "Ordens p = " << p << " and Ordens pr = " << pr <<"\n";
        for(h=7; h<8; h++)
        {
            saidaerros << "\nRefinamento: h = "<< h <<"\n";
            
            //if(dim>1){
            //TPZAutoPointer<TPZGeoMesh> gmesh;
             TPZGeoMesh *gmesh;
            if(dim>1){
                gmesh = MalhaGeom(1.,1.,false);
                TPZVec<int> dims(dim,dim);
                if(dim>1) dims[0]=1;
                int nref = h;
                RefinamentoUniforme(gmesh, nref, dims);
            }else
                gmesh = MalhaGeom1D(3.,h);
            
            ofstream arg("gmesh.txt");
            gmesh->Print(arg);
            
            TPZCompMesh * cmesh1;
            TPZCompMesh * cmesh2;
            if(dim>1){
                cmesh1 = CompMesh(gmesh,pr);
                cmesh2 = CompMesh(gmesh,p);
            }else{
                cmesh1 = CompMesh1D(gmesh,pr);
                cmesh2 = CompMesh1D(gmesh,p);
            }
            
            // Criando a malha computacional multifísica
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            TPZCompMesh * mphysics;
            if(dim>1){
                mphysics = MalhaMDP(meshvec,gmesh);
            }else
                mphysics = MalhaMDP1D(meshvec,gmesh);

            //analysis
            TPZAnalysis an(mphysics,false);
            TPZStepSolver<STATE> step;
            //TPZBandStructMatrix bst(mphysics);
            //TPZFStructMatrix bst(mphysics);
            TPZSkylineNSymStructMatrix bst(mphysics);
            an.SetStructuralMatrix(bst);
            //bst.SetNumThreads(8);
            step.SetDirect(ELU);
            an.SetSolver(step);
            an.Assemble();
            an.Solve();
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            
            std::string plotfile("result.vtk");
            TPZStack<std::string> scalnames,vecnames;
            scalnames.Push("Solution");
            scalnames.Push("ExactSolution");
            scalnames.Push("OptimalTestFunction");
            an.DefineGraphMesh(mphysics->Dimension(), scalnames, vecnames, plotfile);
            an.PostProcess(0,dim);
            
            saidaerros<<"\n\nErro da simulacao MDP  para a pressao";
//            TPZVec<REAL> erros(3);
//            TPZAnalysis an2(cmesh2);
//            an2.SetExact(*SolSuave);
//            an2.PostProcessError(erros, saidaerros);
            ErrorH1(cmesh2, saidaerros);
            
            mphysics->CleanUp();
            gmesh->CleanUp();
            delete mphysics;
        }
    }
    
    return EXIT_SUCCESS;
}


TPZCompMesh * CompMesh(TPZGeoMesh *gmesh, int porder)
{
    /// criar materiais
	int dim = gmesh->Dimension();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    TPZMatLaplacian *material = new TPZMatLaplacian(1,dim);
    
//    TPZAutoPointer<TPZFunction<REAL> > forcef = new TPZDummyFunction<REAL>(ForceSuave);
//    material->SetForcingFunction(forcef);
    
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForceSuave);
    dum->SetPolynomialOrder(20);
    force = dum;
    material->SetForcingFunction(force);
    
    
    TPZAutoPointer<TPZFunction<STATE> > solExata= new TPZDummyFunction<STATE>(SolSuave);
    material->SetForcingFunctionExact(solExata);
    
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);

	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,1.), val2(2,1,0.);
	
    //BC -1
    TPZMaterial * BCondD1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD1->SetForcingFunction(bcmatDirichlet1);
    cmesh->InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet2 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    cmesh->InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZMaterial * BCondD3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet3 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD3->SetForcingFunction(bcmatDirichlet3);
    cmesh->InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet4 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    cmesh->InsertMaterialObject(BCondD4);
    
    
    //Fazendo auto build
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

//malha multifisica para o metodo da dupla projecao
TPZCompMesh *MalhaMDP(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =2;
    TPZMDPMaterial *material = new TPZMDPMaterial(1,dim);
    
    //incluindo os dados do problema
    REAL coefk = 1.;
    material->SetParameters(coefk, 0.);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolSuave);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForceSuave);
    dum->SetPolynomialOrder(20);
    force = dum;
    material->SetForcingFunction(force);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    int boundcond = dirichlet;
    //BC -1
    TPZMaterial * BCondD1 = material->CreateBC(mat, bc1,boundcond, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD1->SetForcingFunction(bcmatDirichlet1);
    mphysics->InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = material->CreateBC(mat, bc2,boundcond, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet2 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    mphysics->InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZMaterial * BCondD3 = material->CreateBC(mat, bc3,boundcond, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet3 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD3->SetForcingFunction(bcmatDirichlet3);
    mphysics->InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = material->CreateBC(mat, bc4,boundcond, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet4 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    mphysics->InsertMaterialObject(BCondD4);

    
    mphysics->InsertMaterialObject(BCondD1);
    mphysics->InsertMaterialObject(BCondD2);
    mphysics->InsertMaterialObject(BCondD3);
    mphysics->InsertMaterialObject(BCondD4);
    
    
    //set multiphysics element
    mphysics->SetDimModel(dim);
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly, bool triang_elements)
{
//    int Qnodes = 4;
//	int64_t dim = 2;
//    
//	TPZGeoMesh * gmesh = new TPZGeoMesh;
//    gmesh->SetDimension(dim);
//	gmesh->SetMaxNodeId(Qnodes-1);
//	gmesh->NodeVec().Resize(Qnodes);
//	TPZVec<TPZGeoNode> Node(Qnodes);
//
//	TPZVec <int64_t> TopolQuad(4);
//	TPZVec <int64_t> TopolLine(2);
//
//	//indice dos nos
//	int64_t id = 0;
//	REAL valx;
//	for(int xi = 0; xi < Qnodes/2; xi++)
//	{
//		valx = xi*Lx;
//		Node[id].SetNodeId(id);
//		Node[id].SetCoord(0 ,valx );//coord X
//		Node[id].SetCoord(1 ,0. );//coord Y
//		gmesh->NodeVec()[id] = Node[id];
//		id++;
//	}
//
//	for(int xi = 0; xi < Qnodes/2; xi++)
//	{
//		valx = Lx - xi*Lx;
//		Node[id].SetNodeId(id);
//		Node[id].SetCoord(0 ,valx );//coord X
//		Node[id].SetCoord(1 ,Ly);//coord Y
//		gmesh->NodeVec()[id] = Node[id];
//		id++;
//	}
//
//	//indice dos elementos
//	id = 0;
//    
//    //elementos internos
//    TopolQuad[0] = 0;
//	TopolQuad[1] = 1;
//	TopolQuad[2] = 2;
//	TopolQuad[3] = 3;
//	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matfinermesh,*gmesh);
//	id++;
//    
//    //elementos de contorno
//	TopolLine[0] = 0;
//	TopolLine[1] = 1;
//	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//	id++;
//	
//	TopolLine[0] = 1;
//	TopolLine[1] = 2;
//	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//	id++;
//	
//	TopolLine[0] = 2;
//	TopolLine[1] = 3;
//	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//	id++;
//	
//	TopolLine[0] = 3;
//	TopolLine[1] = 0;
//	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
//	id++;
//    
//    //construir a malha
//	gmesh->BuildConnectivity();
//	
//	return gmesh;

    TPZManVector<int,2> nx(2,1);
    TPZManVector<REAL,3> x0(3,0.),x1(3,0);
    x1[0] = Lx;
    x1[1] = Ly;
    TPZGenGrid gengrid(nx,x0,x1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    if(triang_elements)
    {
        gengrid.SetElementType(EOned);
    }
    gengrid.Read(gmesh);
    
    //elementos de contorno
    TPZManVector<REAL,3> firstpoint(3,0.),secondpoint(3,0.);
    
    secondpoint[0] = Lx;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc1);
    
    firstpoint = secondpoint;
    secondpoint[1] = Ly;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc2);
    
    firstpoint = secondpoint;
    secondpoint[0] = 0.;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc3);
    
    firstpoint = secondpoint;
    secondpoint[1] = 0.;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc4);
    
    return gmesh;
}

TPZGeoMesh *GMeshSteklov(bool triang_elements)
{
    TPZManVector<int,2> nx(2,2);
    REAL extent = 1;
    nx[1] =1;
    TPZManVector<REAL,3> x0(3,0.),x1(3,extent);
    x0[0] = -extent;
    TPZGenGrid gengrid(nx,x0,x1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    if(triang_elements)
    {
        gengrid.SetElementType(EOned);
    }
    gengrid.Read(gmesh);
    
    //elementos de contorno
    TPZManVector<REAL,3> firstpoint(3,0.),secondpoint(3,0.);
    firstpoint[0] = extent;
    secondpoint[0] = extent;
    secondpoint[1] = extent;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc2);
    gengrid.SetBC(gmesh,6,bc3);
    gengrid.SetBC(gmesh,7,bc4);
    firstpoint[0] = -extent;
    secondpoint[0] = 0.;
    secondpoint[1] = 0.;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc5);
    firstpoint = secondpoint;
    secondpoint[0] = extent;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc1);
    
    return gmesh;
}

void RefinamentoUniforme(TPZAutoPointer<TPZGeoMesh> gmesh, int nref,TPZVec<int> dims)
{
    
    int ir, iel, k;
    int nel=0, dim=0;
    int ndims = dims.size();
	for(ir = 0; ir < nref; ir++ )
    {
		TPZVec<TPZGeoEl *> filhos;
        nel = gmesh->NElements();
        
		for (iel = 0; iel < nel; iel++ )
        {
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            
            dim = gel->Dimension();
            
            for(k = 0; k<ndims; k++)
            {
                if(dim == dims[k])
                {
                    gel->Divide (filhos);
                    break;
                }
            }
		}
	}
    
}


void InsertMaterialObjectsMHM(TPZCompMesh &cmesh, bool useDPGPhil, bool useDPG)
{
	/// criar materiais
	int dim = cmesh.Dimension();
    TPZMatLaplacianLagrange *materialFiner = new TPZMatLaplacianLagrange(matfinermesh,dim);
    
//    TPZAutoPointer<TPZFunction<REAL> > forcef = new TPZDummyFunction<REAL>(ForceSuave);
//    materialFiner->SetForcingFunction(forcef);
    
    TPZAutoPointer<TPZFunction<STATE> > solExata= new TPZDummyFunction<STATE>(SolSuave);
    materialFiner->SetForcingFunctionExact(solExata);
    
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForceSuave);
    dum->SetPolynomialOrder(20);
    forcef = dum;
    materialFiner->SetForcingFunction(forcef);
    
    
    if(useDPGPhil==true){
        materialFiner->SetDPGPhil();
        useDPG=false;
    }
    if(useDPG==true){
        materialFiner->SetMDP();
        useDPGPhil=false;
    }
    
	TPZMaterial * mat1(materialFiner);
    
    TPZMat1dLin *materialSkeleton = new TPZMat1dLin(matskeletonmesh);
    TPZFNMatrix<1,STATE> xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
    materialSkeleton->SetMaterial(xk, xc, xb, xf);
    
	cmesh.InsertMaterialObject(mat1);
    cmesh.InsertMaterialObject(materialSkeleton);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,1.), val2(2,1,0.);
	
    int boundcond = dirichlet;
    if(useDPG) boundcond = neumanndirichlet;
    
    //BC -1
    TPZMaterial * BCondD1 = materialFiner->CreateBC(mat1, bc1,boundcond, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD1->SetForcingFunction(bcmatDirichlet1);
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = materialFiner->CreateBC(mat1, bc2,boundcond, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet2 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZMaterial * BCondD3 = materialFiner->CreateBC(mat1, bc3,boundcond, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet3 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD3->SetForcingFunction(bcmatDirichlet3);
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = materialFiner->CreateBC(mat1, bc4,boundcond, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet4 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    cmesh.InsertMaterialObject(BCondD4);
}

void InsertMaterialObjects(TPZCompMesh &cmesh)
{
	/// criar materiais
	int dim = cmesh.Dimension();
    
    TPZMatLaplacianLagrange *materialcoarse = new TPZMatLaplacianLagrange(matcoarsemesh,dim);
//    TPZAutoPointer<TPZFunction<REAL> > forcef = new TPZDummyFunction<REAL>(ForceSuave);
//    materialcoarse->SetForcingFunction(forcef);
    
	TPZMaterial * mat1(materialcoarse);
    
	cmesh.InsertMaterialObject(mat1);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,1.), val2(2,1,0.);
	
    //BC -1
    TPZMaterial * BCondD1 = materialcoarse->CreateBC(mat1, bc1,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD1->SetForcingFunction(bcmatDirichlet1);
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = materialcoarse->CreateBC(mat1, bc2,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet2 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZMaterial * BCondD3 = materialcoarse->CreateBC(mat1, bc3,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet3 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD3->SetForcingFunction(bcmatDirichlet3);
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = materialcoarse->CreateBC(mat1, bc4,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet4 = new TPZDummyFunction<REAL>(DirichletSuave);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    cmesh.InsertMaterialObject(BCondD4);
}


void GetElIndexCoarseMesh(TPZAutoPointer<TPZGeoMesh>  gmesh, std::set<int64_t> &coarseindex)
{
    int nel = gmesh->NElements();
    int iel;
    int hassubel=0;
    int dim = gmesh->Dimension();
    int eldim;
    for(iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        hassubel = gel->HasSubElement();
        eldim = gel->Dimension();
        if(!hassubel && eldim ==dim)
        {
            coarseindex.insert(gel->Index());
        }
    }
    
}

void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL sol = /*x*(1.-x)*y*(1.-y);*/sin(M_PI*x)*sin(M_PI*y);
    u[0] = sol;
    du.Resize(2, 1);
    du(0,0) = /*(1. - x)*(1. - y)*y - x*(1. - y)*y;*/M_PI*cos(M_PI*x)*sin(M_PI*y);
    du(1,0) = /*(1. - x)*x*(1. - y) - (1. - x)*x*y;*/M_PI*cos(M_PI*y)*sin(M_PI*x);
}

void ForceSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &force){
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL f = /*2.*(1. - x)*x + 2.*(1. - y)*y;*/2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);// + sin(M_PI*x)*sin(M_PI*y);
    force[0] = f;
}

void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL n = 0;
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = pow((REAL)2,0.25 + n/2.)*pow(r,0.5 + n)*cos((0.5 + n)*t);
    u[0] = sol;
    
    du(0,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    
}

void RefinamentoSingular(TPZAutoPointer<TPZGeoMesh> gmesh,int nref)
{
    int64_t nnodes = gmesh->NNodes();
    int64_t in;
    for (in=0; in<nnodes; in++) {
        TPZGeoNode *gno = &gmesh->NodeVec()[in];
        if (abs(gno->Coord(0))< 1.e-6 && abs(gno->Coord(1)) < 1.e-6) {
            break;
        }
    }
    if (in == nnodes) {
        DebugStop();
    }
    TPZGeoElSide gelside;
    int64_t nelem = gmesh->NElements();
    for (int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        int ncorner = gel->NCornerNodes();
        for (int ic=0; ic<ncorner; ic++) {
            int64_t nodeindex = gel->NodeIndex(ic);
            if (nodeindex == in) {
                gelside = TPZGeoElSide(gel, ic);
                break;
            }
        }
        if (gelside.Element()) {
            break;
        }
    }
    if (!gelside.Element()) {
        DebugStop();
    }
    for (int iref = 0; iref <nref; iref++) {
        TPZStack<TPZGeoElSide> gelstack;
        gelstack.Push(gelside);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            gelstack.Push(neighbour);
            neighbour = neighbour.Neighbour();
        }
        int64_t nstack = gelstack.size();
        for (int64_t ist=0; ist < nstack; ist++) {
            if (!gelstack[ist].Element()->HasSubElement()) {
                TPZVec<TPZGeoEl *> subel;
                gelstack[ist].Element()->Divide(subel);
            }
        }
    }
}

void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out)
{
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<STATE,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolSuave, elerror, NULL);
        
        int nerr = elerror.size();
        //globerrors.resize(nerr);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
    }
    out << "\n";
    out << "Errors associated with L2 or H1 space\n";
    out << "\nH1 Norm = "    << sqrt(globerrors[0]) << endl;
    out << "\nL2 Norm = "    << sqrt(globerrors[1]) << endl;
    out << "\nSemi H1 Norm = "    << sqrt(globerrors[2]) << endl;
    out << "\n=============================\n"<<endl;
}

TPZGeoMesh *MalhaGeom1D(REAL Lx, int h)
{
    int Qnodes = 2;
    int64_t dim = 1;

    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
    
    //indice dos nos
    int64_t id = 0;
    REAL valx;
    for(int xi = 0; xi < Qnodes; xi++)
    {
        valx = xi*Lx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    //indice dos elementos
    id = 0;

    //elementos internos
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (id,TopolLine,1,*gmesh);
    id++;
    
    //elementos de contorno
    TopolPoint[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc1,*gmesh);
    id++;

    TopolPoint[0] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bc2,*gmesh);
    
    //construir a malha
    gmesh->BuildConnectivity();
    
    //refinamento uniforme
    for(int i = 0; i < h; i++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    
    return gmesh;
}

TPZCompMesh *MalhaMDP1D(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh){
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =1;
    TPZMDPMaterial *material = new TPZMDPMaterial(1,dim);
    
    //incluindo os dados do problema
    REAL coefk = fk;
    TPZVec<REAL> Conv(3,0.);
    Conv[0] = fc;
    material->SetParameters(coefk, 0.);
    material->SetConvection(Conv);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(Sol1D);
    material->SetForcingFunctionExact(solexata);
    
//    //funcao do lado direito da equacao do problema
//    TPZAutoPointer<TPZFunction<STATE> > force;
//    TPZDummyFunction<STATE> *dum;
//    dum = new TPZDummyFunction<STATE>(ForceSuave);
//    dum->SetPolynomialOrder(20);
//    force = dum;
//    material->SetForcingFunction(force);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    int boundcond = dirichlet;
    
    //BC -2
    TPZMaterial * BCondD2 = material->CreateBC(mat, bc2,boundcond, val1, val2);
    mphysics->InsertMaterialObject(BCondD2);
    
//    val2(0,0)=0.;
//    TPZMaterial * BCondD2 = material->CreateBC(mat, bc2,neumann, val1, val2);
//    TPZAutoPointer<TPZFunction<REAL> > bcmatNeum = new TPZDummyFunction<REAL>(Neumann1D);
//    BCondD2->SetForcingFunction(bcmatNeum);
//    mphysics->InsertMaterialObject(BCondD2);
    
     //BC -1
    val2(0,0) = 10.;
    TPZMaterial * BCondD1 = material->CreateBC(mat, bc1,boundcond, val1, val2);
    mphysics->InsertMaterialObject(BCondD1);

    
    mphysics->InsertMaterialObject(BCondD1);
    mphysics->InsertMaterialObject(BCondD2);
    
    //set multiphysics element
    mphysics->SetDimModel(dim);
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

TPZCompMesh * CompMesh1D(TPZGeoMesh *gmesh, int porder)
{
    /// criar materiais
    int dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    TPZMatPoisson3d *material = new TPZMatPoisson3d(1,dim);
    TPZVec<REAL> convdir(3,0.);
    convdir[0]=1.;
    material->SetParameters(fk, fc, convdir);
    
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //BC -2
//    TPZMaterial * BCondD2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
//    cmesh->InsertMaterialObject(BCondD2);
    
    //BC -1
    val2(0,0)=10.;
    TPZMaterial * BCondD1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCondD1);
    
    //BC -2
    val2(0,0)=0.;
    TPZMaterial * BCondD2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatNeum = new TPZDummyFunction<REAL>(Neumann1D);
    BCondD2->SetForcingFunction(bcmatNeum);
    cmesh->InsertMaterialObject(BCondD2);
    
    
    //Fazendo auto build
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

void Sol1D(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    du.Resize(1,1);
    u.Resize(1);
    const REAL x = loc[0];
    const REAL Peclet = fc/fk;
    const REAL term1 = Peclet*(x-3.);
    const REAL term2 = -3.*Peclet;
    
    const REAL sol = 10.*(1.-exp(term1))/(1.-exp(term2));
    const REAL dsol = (-10.*Peclet)*(exp(term1)/(1.-exp(term2)));
    
    u[0] = sol;
    du(0,0) = dsol;
}


void Neumann1D(const TPZVec<REAL> &loc, TPZVec<STATE> &result)
{
    result.Resize(2, 0.);
    REAL dun = -((10.*fc)/(fk - exp(-((3.*fc)/fk))*fk));
    result[0] = dun;
}

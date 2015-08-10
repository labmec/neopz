//
//  main_mixed.cpp
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzinterpolationspace.h"
#include "TPZCompElDisc.h"
#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "tpzchangeel.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzgengrid.h"

#include "pzlog.h"

#include "pzl2projection.h"
#include "pzmultiphysicselement.h"
#include "pzintel.h"
#include "TPZVTKGeoMesh.h"

#include "TPZMatDualHybridPoisson.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZSkylineNSymStructMatrix.h"

#include <iostream>
#include <math.h>
using namespace std;

int const matId =1;
int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;

int const dirichlet =0;
int const neumann = 1;

REAL const Pi = 4.*atan(1.);
TPZGeoMesh *GMesh();

//with hdiv
TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson* &mymaterial);

//with hybrid method
TPZGeoMesh * MalhaGeo(const int h);
TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder,bool ismultiplierH1);
void GroupElements(TPZCompMesh *cmesh);
void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);


void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide);
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void Prefinamento(TPZCompMesh * cmesh, int ndiv,int porder);

void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numThreads=0);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out);
void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out);

void PrefinamentoRibsHybridMesh(TPZCompMesh *cmesh);
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("PiraHP.main"));
#endif


#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#define print

bool runhdiv = true;
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    
    ///Refinamento
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
    int ndiv = 1;
    int p = 2;
    
    //string outputfile("Solution_mphysics");
    std::stringstream name;
    if(runhdiv==true){
        HDivPiola = 1;
        std::ofstream myerrorfile("../Simulacao-MistaHdiv.txt");
        for(ndiv=1; ndiv<3; ndiv++){
            TPZGeoMesh *gmesh = GMesh();
            //malha inicial
            UniformRefine(gmesh,1);

            //refinamento proximo do no de id=1
            DirectionalRef(gmesh, 1, ndiv);

//            #ifdef print
//            {
//                std::ofstream malhaOut("malhageometrica.vtk");
//                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
//            }
//            #endif

            //mesh1
            TPZCompMesh * cmesh1= CMeshFlux(gmesh, p);
            Prefinamento(cmesh1, ndiv,p);

            //mesh2
            TPZCompMesh * cmesh2= CMeshPressure(gmesh, p);
            Prefinamento(cmesh2, ndiv,p);

            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            TPZMixedPoisson * mymaterial;
            TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,mymaterial);

            //resolver problema
            TPZAnalysis an(mphysics);
            ResolverSistema(an, mphysics,0);

            //pos-process
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            std::stringstream name;
            name << "Solution_mphysics" <<ndiv<< ".vtk";
            std::string paraviewfile(name.str());
            PosProcessMultphysics(meshvec,  mphysics, an, paraviewfile);
            
            //Erro global
            myerrorfile << "\nh = "<< ndiv << " p = " << p << "\n";
            myerrorfile << "neq = " << mphysics->NEquations() << "\n";
            TPZVec<REAL> erros(3);
            myerrorfile<<"\n\nErro da simulacao multifisica  para o flux";
            TPZAnalysis an1(cmesh1,false);
            //an1.SetExact(*SolExataSteklov);
            //an1.PostProcessError(erros, myerrorfile);
            ErrorHDiv(cmesh1,myerrorfile);
            
            myerrorfile<<"\n\nErro da simulacao multifisica  para a pressao";
            TPZAnalysis an2(cmesh2,false);
//            an2.SetExact(*SolExataSteklov);
//            an2.PostProcessError(erros, myerrorfile);
            ErrorL2(cmesh2,myerrorfile);
            
            {
                std::ofstream mesh1("../cmesh1.txt");
                cmesh1->Print(mesh1);
            }
            cmesh1->CleanUp();
            cmesh2->CleanUp();
            delete cmesh1;
            delete cmesh2;
            delete gmesh;
        }
    }
    else{//metodo misto: Abimael
    
        bool multiplicadorH1 = true;
        std::ofstream myerrorfile("Simulacao-Primal-LDGC.txt");
        
        for(ndiv=1; ndiv<10; ndiv++){
            TPZGeoMesh *gmesh = MalhaGeo(1);//malha geometrica
            
// ------- refinamento proximo do no de id=1 -----
            
            DirectionalRef(gmesh, 1, ndiv);
    
#ifdef print
            {
                std::ofstream malhaOut("malhageometricaHibrid.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
            }
#endif
            
            TPZCompMesh *cmesh = CreateHybridCompMesh(*gmesh, p, multiplicadorH1);//malha computacionalx
            
            {
            std::ofstream out("cmeshHib1.txt");
            cmesh->Print(out);
                
            std::ofstream out2("gmesh.txt");
            gmesh->Print(out2);
            }

            
            Prefinamento(cmesh, ndiv,p);
            PrefinamentoRibsHybridMesh(cmesh);
            
            cmesh->CleanUpUnconnectedNodes();
            cmesh->ExpandSolution();
            std::ofstream out("cmeshHib2.txt");
            cmesh->Print(out);
            GroupElements(cmesh);
            
            cmesh->LoadReferences();//mapeia para a malha geometrica lo
            
            TPZAnalysis analysis(cmesh);
            
            TPZSkylineNSymStructMatrix str(cmesh);
           // str.SetNumThreads(8);
            //TPZFStructMatrix str(cmesh);
            
            TPZAutoPointer<TPZMatrix<STATE> > mat = str.Create();
            str.EquationFilter().Reset();
            TPZAutoPointer<TPZMatrix<STATE> > mat2 = mat->Clone();
            
            analysis.SetStructuralMatrix(str);
            TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
            TPZStepSolver<STATE> *gmrs = new TPZStepSolver<STATE>(mat2);
            step->SetReferenceMatrix(mat2);
            step->SetDirect(ELU);
            gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
            TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
            TPZAutoPointer<TPZMatrixSolver<STATE> > autogmres = gmrs;
            analysis.SetSolver(autogmres);
            
            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            analysis.Assemble();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            analysis.Solve();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
            std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
            
            myerrorfile << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
#endif
            
            
            
            TPZVec<std::string> scalnames(2),vecnames(0);
            scalnames[0] = "Solution";
            scalnames[1] = "POrder";
            
            std::stringstream name;
            name << "Solution_bima" <<ndiv<< ".vtk";
            std::string paraviewfile(name.str());
            analysis.DefineGraphMesh(2,scalnames,vecnames,paraviewfile);
            analysis.PostProcess(0);
            
            myerrorfile << "\nh = "<< ndiv << " p = " << p << "\n";
            myerrorfile << "neq = " << cmesh->NEquations() << "\n";
            
            //erro global
           
            analysis.SetExact(SolExataSteklov);
            TPZVec<REAL> erros(3);
            analysis.PostProcessError(erros);
            myerrorfile << "H1 = " << erros[0];
            myerrorfile << " L2 = " << erros[1];
            myerrorfile << " semi H1 = " << erros[2] << "\n"<<std::endl;
        }
       
    }
    
    return 0;
}

TPZGeoMesh *GMesh(){
    
    int dim = 2;
    int Qnodes = 6;
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
    gmesh->SetDimension(dim);
    TPZGeoEl *gel = 0;
	
	TPZVec <long> TopolQuad(4);
	TPZVec <long> TopolLine(2);
    
    TPZStack<TPZGeoEl *> toChange;
    TPZStack<int> targetSide;
	
	//indice dos nos
	long id = 0;
	REAL valx, dx = 0.5;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = -0.5 + xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = 0.5 - xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0.5 );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 4;
    TopolQuad[3] = 5;
    gel = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
    toChange.Push(gel);
    targetSide.Push(1);
    
    id++;
    
    TopolQuad[0] = 1;
    TopolQuad[1] = 2;
    TopolQuad[2] = 3;
    TopolQuad[3] = 4;
    gel = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
    toChange.Push(gel);
    targetSide.Push(0);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    gel = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
    toChange.Push(gel);
    targetSide.Push(1);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 2;
    gel = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
    toChange.Push(gel);
    targetSide.Push(0);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
    id++;
    
    TopolLine[0] = 4;
    TopolLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
    id++;
    
    TopolLine[0] = 5;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc5,*gmesh);

	gmesh->BuildConnectivity();
    
    for (long el=0; el<toChange.size(); el++) {
        TPZGeoEl *gel = toChange[el];
        TPZChangeEl::ChangeToQuadratic(gmesh, gel->Index());
//        TPZChangeEl::ChangeToQuarterPoint(gmesh, gel->Index(), targetSide[el]);
    }
   
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
	return gmesh;
}


void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    	gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"gmesh depois de refinar uniformemente\n";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
}


TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = gmesh->Dimension();
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	TPZMaterial * mat(material);
	material->NStateVariables();
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    TPZMaterial * BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
	
    cmesh->SetAllCreateFunctionsHDiv();
    
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();//ajusta as condicoes de contorno
    cmesh->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional fluxo\n";
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif

	return cmesh;
}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder)
{
    bool triang = false;
    /// criar materiais
	int dim = gmesh->Dimension();
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
	material->NStateVariables();
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
//	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
//    TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
//    TPZMaterial * BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
//    TPZMaterial * BCond6 = material->CreateBC(mat, bc6,dirichlet, val1, val2);
//    cmesh->InsertMaterialObject(BCond1);
//    cmesh->InsertMaterialObject(BCond2);
//    cmesh->InsertMaterialObject(BCond3);
//    cmesh->InsertMaterialObject(BCond4);
//    cmesh->InsertMaterialObject(BCond5);
//    cmesh->InsertMaterialObject(BCond6);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
	newnod.SetLagrangeMultiplier(1);
    }
    
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        celdisc->SetConstC(1.);
        celdisc->SetCenterPoint(0, 0.);
        celdisc->SetCenterPoint(1, 0.);
        celdisc->SetCenterPoint(2, 0.);
        celdisc->SetTrueUseQsiEta();
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(triang==true) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
        
    }
    
    
#ifdef DEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    cmesh->AdjustBoundaryElements();//ajusta as condicoes de contorno
    cmesh->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional pressao\n";
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
	return cmesh;
}

TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson * &mymaterial){
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=1;
    int dim = gmesh->Dimension();
    
    REAL coefk=1.;
    mymaterial = new TPZMixedPoisson(MatId,dim);
    mymaterial->SetPermeability(coefk);
    
    TPZFMatrix<REAL> TensorK(dim,dim,0.);
    TPZFMatrix<REAL> InvK(dim,dim,0.);
    for(int i = 0; i<dim; i++) {
        TensorK(i,i)=coefk;
        InvK(i,i)=coefk;
    }
    mymaterial->SetPermeabilityTensor(TensorK, InvK);

    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    TPZAutoPointer<TPZFunction<STATE> > solExata;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(SolExataSteklov);
    dum->SetPolynomialOrder(20);
    solExata = dum;
    mymaterial->SetForcingFunctionExact(solExata);
    
    //Inserir condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    BCond1 = mymaterial->CreateBC(mat, bc1, dirichlet, val1, val2);//Dirichlet nulo
    BCond2 = mymaterial->CreateBC(mat, bc2, neumann, val1, val2);//Neumann nulo
    BCond3 = mymaterial->CreateBC(mat, bc3, neumann, val1, val2);
    BCond4 = mymaterial->CreateBC(mat, bc4, neumann, val1, val2);
    BCond5 = mymaterial->CreateBC(mat, bc5, neumann, val1, val2);
    
    //Set force function
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannDireito;
    bcmatNeumannDireito = new TPZDummyFunction<STATE>(NeumannDireita);
    BCond3->SetForcingFunction(bcmatNeumannDireito);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
    bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcima);
    BCond4->SetForcingFunction(bcmatNeumannAcima);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEsquerdo;
    bcmatNeumannEsquerdo = new TPZDummyFunction<STATE>(NeumannEsquerda);
    BCond5->SetForcingFunction(bcmatNeumannEsquerdo);

    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional multifisica\n";
        mphysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    return mphysics;
}


void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(3), vecnames(0);
	//vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
    
    scalnames[1] = "ExactPressure";
    //vecnames[1] = "ExactFlux";
     scalnames[2] = "POrder";
			
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
//	std::ofstream out("malha.txt");
//	an.Print("nothing",out);
    
}


//refinamento uniforme em direcao ao no
void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide){
    
    ///Refinamento
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
    for (int idivide = 0; idivide < divide; idivide++){
        const int nels = gmesh->NElements();
        TPZVec< TPZGeoEl * > allEls(nels);
        for(int iel = 0; iel < nels; iel++){
            allEls[iel] = gmesh->ElementVec()[iel];
        }
        
        for(int iel = 0; iel < nels; iel++){
            TPZGeoEl * gel = allEls[iel];
            if(!gel) continue;
            if(gel->HasSubElement()) continue;
            int nnodes = gel->NNodes();
            int found = -1;
            for(int in = 0; in < nnodes; in++){
                if(gel->NodePtr(in)->Id() == nodeAtOriginId){
                    found = in;
                    break;
                }
            }///for in
            if(found == -1) continue;
            
            MElementType gelT = gel->Type();
            TPZAutoPointer<TPZRefPattern> uniform = gRefDBase.GetUniformRefPattern(gelT);
            if(!uniform){
                DebugStop();
            }
            gel->SetRefPattern(uniform);
            TPZVec<TPZGeoEl*> filhos;
            gel->Divide(filhos);
            
        }///for iel
    }//idivide
    
    gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"gmesh depois de refinar direcionalmente\n";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif

}///void

void PrefinamentoRibsHybridMesh(TPZCompMesh *cmesh){
    const int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZCompEl * cel = cmesh->ElementVec()[iel];
        if(!cel)continue;
        TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
        if(!face)continue;
        TPZInterpolationSpace *TwoD = NULL;
        TPZInterpolationSpace *Rib = NULL;
        if(face->LeftElement()->Reference()->Dimension()==2){
            TwoD = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
            Rib = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
        }
        else{
            Rib = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
            TwoD = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
        }
        if(!Rib || !TwoD) DebugStop();
        int porder2D = TwoD->MaxOrder();
        int pOrder1D = Rib->MaxOrder();
        int neworder = pOrder1D > porder2D ? pOrder1D : porder2D;
        Rib->PRefine(neworder);
    }
    
}

void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder){
    if(ndiv<1) return;
    int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        if(!sp) continue;
        int level = sp->Reference()->Level();
        TPZGeoEl * gel = sp->Reference();
        if(gel->Dimension()==2)
            sp->PRefine(porder + (ndiv - 1) + (level-1));
   }
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional apos pRefinamento\n";
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
}

void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    u.Resize(1, 0.);
    du.Resize(2, 1);
    du(0,0)=0.;
    du(1,0)=0.;
    
    const REAL n = 0;
    REAL x = loc[0];
    REAL y = loc[1];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = pow((REAL)2,0.25 + n/2.)*pow(r,0.5 + n)*cos((0.5 + n)*t);
    u[0] = sol;
    
//    if(IsZero(x) && IsZero(y)){
//        y=y+1.e-3;
//        x=x+1.e-2;
//    }
    
    //flux = -k*grad(u), k=1 nesse problema
    if(runhdiv==true){
        du(0,0) = -pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
        du(1,0) = -pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    }else{
        du(0,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
        du(1,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    }
}

void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,result,du);
}

void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {-1,0};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {+1,0};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {0,+1};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

#include "pzbstrmatrix.h"
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numThreads)
{
    //TPZBandStructMatrix full(fCmesh);
    TPZSkylineStructMatrix full(fCmesh); //caso simetrico
    full.SetNumThreads(numThreads);
    an.SetStructuralMatrix(full);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); //caso simetrico
    //step.SetDirect(ELU);
    an.SetSolver(step);
    an.Run();
    
    //Saida de Dados: solucao e  grafico no VT
    ofstream file("../Solout.txt");
    an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out)
{
    long nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<STATE,10> elerror(10,0.);
        cel->EvaluateError(SolExataSteklov, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
}

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out)
{
    long nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    globalerrors.Fill(0.);
    for (long el=0; el<nel; el++) {
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
        cel->EvaluateError(SolExataSteklov, elerror, NULL);
        int nerr = elerror.size();
        // globalerrors.resize(nerr);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with L2 space"<< endl;
    out << "L2 Norm = "    << sqrt(globalerrors[1]) << endl;
}


//----------------------------------
TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder, bool ismultiplierH1){
    //TPZCompEl::SetgOrder(porder);
    TPZCompMesh *comp = new TPZCompMesh(&gmesh);
    comp->SetDimModel(gmesh.Dimension());
    
    comp->ApproxSpace().CreateDisconnectedElements(true);
    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    
    // Criar e inserir os materiais na malha
    REAL beta = 6;
    TPZMatDualHybridPoisson *mymaterial = new TPZMatDualHybridPoisson(1,0.,beta);
    TPZMaterial * automat(mymaterial);
    comp->InsertMaterialObject(automat);
    
    
    // Condicoes de contorno
    TPZFMatrix<STATE> val1(1,1,1.),val2(1,1,0.);
    
    TPZMaterial *bnd = automat->CreateBC (automat,-1,2,val1,val2);//misto tbem
    bnd->SetForcingFunction(Dirichlet2);
    comp->InsertMaterialObject(bnd);

    // Mixed
    bnd = automat->CreateBC (automat,-2,0,val1,val2);
    bnd->SetForcingFunction(Dirichlet);
    comp->InsertMaterialObject(bnd);

    // Mixed
    val1(0,0)=0.;
    val2(0,0)=0.;
    bnd = automat->CreateBC (automat,-3,0,val1,val2);//Dirichlet nulo
    comp->InsertMaterialObject(bnd);

    // Mixed
    bnd = automat->CreateBC (automat,-4,0,val1,val2);
    bnd->SetForcingFunction(Dirichlet);
    comp->InsertMaterialObject(bnd);

    // Mixed
    bnd = automat->CreateBC (automat,-5,0,val1,val2);
    bnd->SetForcingFunction(Dirichlet);
    comp->InsertMaterialObject(bnd);

    // Mixed
    bnd = automat->CreateBC (automat,-6,0,val1,val2);
    bnd->SetForcingFunction(Dirichlet);
    comp->InsertMaterialObject(bnd);

    // Ajuste da estrutura de dados computacional
    comp->ApproxSpace().CreateDisconnectedElements(true);
    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    std::set<int> matids;
    matids.insert(1);
    comp->AutoBuild(matids);
    
    comp->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    matids.clear();
    for (int i=2; i<=6; i++) {
        matids.insert(-i);
    }
    
    comp->SetDefaultOrder(porder);
    comp->AutoBuild(matids);
    comp->SetDimModel(2);
    
    comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
    comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
    comp->LoadReferences();
    comp->ApproxSpace().CreateInterfaceElements(comp,true);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    matids.insert(1);
    for (int i=2; i<=6; i++) {
        matids.insert(-i);
    }
    
    comp->ApproxSpace().Hybridize(*comp, matids, ismultiplierH1);
    
    comp->SetName("Malha Computacional Original");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    return comp;
}

TPZGeoMesh * MalhaGeo(const int ndiv){//malha quadrilatera
    
    ///malha geometrica
    TPZGeoMesh * gmesh = new TPZGeoMesh();
    
    ///Criando nós
    const int nnodes = 6;
    double coord[nnodes][2] = {{-0.5,0},{0,0},{0.,0.5},{-0.5,0.5},{0.5,0},{0.5,0.5}};
    for(int i = 0; i < nnodes; i++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL,3> nodeCoord(3);
        nodeCoord[0] = coord[i][0];
        nodeCoord[1] = coord[i][1];
        nodeCoord[2] = 0.;
        gmesh->NodeVec()[nodind].Initialize(i,nodeCoord,*gmesh);
    }
    
    ///Criando elementos
    const int nel = 2;
    int els[nel][4] = {{0,1,2,3},{1,4,5,2}};
    for(int iel = 0; iel < nel; iel++){
        TPZManVector<long,4> nodind(4);
        long index;
        nodind[0] = els[iel][0];
        nodind[1] = els[iel][1];
        nodind[2] = els[iel][2];
        nodind[3] = els[iel][3];
        gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    ///Criando elementos de contorno
    const int nelbc = 6;
    int bcels[nelbc][3] = {{0,1,-3},{1,4,-2},{4,5,-4},{5,2,-6},{2,3,-6},{3,0,-5}};
    for(int iel = 0; iel < nelbc; iel++){
        TPZManVector<long,4> nodind(2);
        long index;
        nodind[0] = bcels[iel][0];
        nodind[1] = bcels[iel][1];
        int matid = bcels[iel][2];
        gmesh->CreateGeoElement(EOned,nodind,matid,index);
    }
    
    ///Construindo conectividade da malha
    gmesh->BuildConnectivity();

#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"gmesh antes de refinar\n";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    ///Refinamento uniforme da malha
    
    ///Inicializando padrões de refinamento uniforme
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    for (int i = 0; i < ndiv; i++){
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++){
            TPZGeoEl * gel = gmesh->ElementVec()[iel];
            if (!gel) continue;
            if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
            TPZVec<TPZGeoEl*> filhos;
            gel->Divide(filhos);
        }///iel
    }///i

    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"gmesh depois de refinar\n";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    return gmesh;
}

void GroupElements(TPZCompMesh *cmesh)
{
    cmesh->LoadReferences();
    int nel = cmesh->NElements();
    std::set<TPZCompEl *> celset;
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();
        if (dim ==2) {
            celset.insert(cel);
        }
    }
    std::set<int> elgroupindices;
    
    for (std::set<TPZCompEl *>::iterator it = celset.begin(); it != celset.end(); it++) {
        
        std::list<TPZCompEl *> group;
        group.push_back(*it);
        TPZCompEl *cel = *it;
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZStack<TPZCompElSide> connected;
            TPZCompElSide celside(cel,is);
            celside.EqualLevelElementList(connected, false, false);
            int neq = connected.NElements();
            for (int eq=0; eq<neq; eq++) {
                TPZCompElSide eqside = connected[eq];
                TPZCompEl *celeq = eqside.Element();
                TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(celeq);
                if (!intface) {
                    continue;
                }
                TPZCompEl *left = intface->LeftElement();
                if (left == cel) {
                    //put in the group
                    group.push_back(intface);
                }
            }
        }
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
                (*it)->Print(sout);
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        long index;
        TPZElementGroup *celgroup = new TPZElementGroup(*cmesh,index);
        elgroupindices.insert(index);
        for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
            celgroup->AddElement(*it);
        }
    }
    cmesh->ComputeNodElCon();
    
    for (std::set<int>::iterator it = elgroupindices.begin(); it!=elgroupindices.end(); it++) {
        TPZCompEl *cel = cmesh->ElementVec()[*it];
        TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel);
    }
    
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}


void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFNMatrix<10,REAL> fake(2,1);
    result[0] = loc[0]*loc[0];
}


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

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzgengrid.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>
using namespace std;

int const matId =1;
int const bc0=-1;
int const bc1=-2;
int const bc2=-3;
int const bc3=-4;
int const bc4=-5;

int const dirichlet =0;
int const neumann = 1;

REAL const Pi = 4.*atan(1.);

TPZGeoMesh *GMesh(bool triang_elements);
TPZGeoMesh *GMesh2(bool triang_elements);

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson* &mymaterial);

void Forcing1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void SolExata1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &flux);
void ForcingBC1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void SolExata2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &flux);
void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
TPZCompMesh *CMeshHDivPressure(TPZGeoMesh *gmesh, int pOrder);
void PosProcessHDiv(TPZAnalysis &an, std::string plotfile);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);

//----- Boundary condition: Steklov
void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);

void Compute_dudn(TPZCompMesh *cmesh);
void Compute_dudnQuadrado(TPZCompMesh *cmesh);
TPZInterpolationSpace * FindInterpolationSpace(TPZVec<REAL> &xVec, TPZCompMesh *cmesh, TPZVec<REAL> &qsi);
void PrintToFile(std::ofstream &myfile, const std::string &title, std::map<REAL,REAL> &mymap);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
#endif

const bool triang = false;
const int teste = 3;
int mainAgnaldo(int argc, char *argv[])
{
#ifdef LOG4CXX
//	std::string logs("../logmixedproblem.cfg");
//	InitializePZLOG("../logmixedproblem.cfg");
    InitializePZLOG();
#endif
    gRefDBase.InitializeAllUniformRefPatterns();
    
    TPZVec<REAL> erros;
    ofstream arg12("Erro.txt");
    arg12<<"\n\nERRO PARA MALHA COM TRIANGULO  = " << triang;
    for (int p = 1; p< 4; p++)
    {
        arg12<<"\n\n----------- ORDEM p = " << p <<" -----------";
        for(int h = 1; h < 6;h++)
        {
            arg12<<"\nRefinamento ndiv  = " << h;
            // int p = 5;
            //primeira malha

            // geometric mesh (initial)
            TPZGeoMesh *gmesh = 0;
            if(teste == 1 || teste == 2)
            {
                gmesh = GMesh(triang);
            }
            else if(teste == 3)
            {
                gmesh = GMesh2(triang);
            }
            //UniformRefine(gmesh,h);
           // ofstream arg1("gmesh_inicial.txt");
            //gmesh->Print(arg1);

          
            // First computational mesh
            TPZCompMesh * cmesh1= CMeshFlux(gmesh, p);
            //ofstream arg2("cmesh1_inicial.txt");
            //cmesh1->Print(arg2);

#ifdef LOG4CXX
						if(logdata->isDebugEnabled())
						{
							std::stringstream sout;
							sout << "Flux mesh\n";
							cmesh1->Print(sout);
							LOGPZ_DEBUG(logdata,sout.str())
						}
#endif

            // Second computational mesh
            TPZCompMesh * cmesh2;
            if(triang==true){
               cmesh2 = CMeshPressure(gmesh, p-1);
            }
            else cmesh2 = CMeshPressure(gmesh, p);
            //ofstream arg3("cmesh2_inicial.txt");
            //cmesh2->Print(arg3);
#ifdef LOG4CXX
						if(logdata->isDebugEnabled())
						{
							std::stringstream sout;
							sout << "Pressure mesh\n";
							cmesh2->Print(sout);
							LOGPZ_DEBUG(logdata,sout.str())
						}
#endif

            // Cleaning reference of the geometric mesh to cmesh1
            gmesh->ResetReference();
            cmesh1->LoadReferences();
            TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,h,false);
            //ofstream arg4("cmesh1_final.txt");
            //cmesh1->Print(arg4);
            

            // Cleaning reference to cmesh2
            gmesh->ResetReference();
            cmesh2->LoadReferences();
            TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,h,true);
           // ofstream arg5("cmesh2_final.txt");
           // cmesh2->Print(arg5);

            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            TPZMixedPoisson * mymaterial;
            TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,mymaterial);
            //ofstream arg6("mphysic.txt");
            //mphysics->Print(arg6);
#ifdef LOG4CXX
						if(logdata->isDebugEnabled())
						{
							std::stringstream sout;
							sout << "Multiphysics mesh\n";
							mphysics->Print(sout);
							LOGPZ_DEBUG(logdata,sout.str())
						}
#endif

            int NEq = mphysics->NEquations();
            
//            ofstream arg7("gmesh_Final.txt");
//            gmesh->Print(arg7);

            TPZAnalysis an(mphysics);
            
//           ofstream arg10("mphysic_posAnalysis.txt");
//            mphysics->Print(arg10);
            
            SolveSyst(an, mphysics);
            
//            ofstream arg6("mphysic.txt");
//            mphysics->Print(arg6);
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            
//            ofstream arg9("mphysic_apos_Transferfrom.txt");
//            mphysics->Print(arg9);

            //Compute_dudnQuadrado(meshvec[0]);
            
            
            arg12<<"\nNumero de Equacoes = " << NEq;
            
            arg12<<"\n\nErro da simulacao multifisica  para o flux";
            TPZAnalysis an1(cmesh1);
            if (teste==1){
                an1.SetExact(*SolExata1);
            }
            else if(teste==2){
                an1.SetExact(*SolExata2);
            }
            else{
                an1.SetExact(*SolExataSteklov);
            }
            an1.PostProcessError(erros, arg12);
            
            
            arg12<<"\n\nErro da simulacao multifisica  para a pressao";
            TPZAnalysis an2(cmesh2);
            if (teste==1){
                an2.SetExact(*SolExata1);
            }
            else if (teste==2){
                an2.SetExact(*SolExata2);
            }
            else{
                an2.SetExact(*SolExataSteklov);
            }
            an2.PostProcessError(erros, arg12);

            //string plotfile("Solution_mphysics.vtk");
            //PosProcessMultphysics(meshvec,  mphysics, an, plotfile);

/*
            //solucao HDivPressure
            TPZCompMesh * cmesh3= CMeshHDivPressure(gmesh, p);
//            ofstream arg8("cmesh_HdivInicial.txt");
//            cmesh3->Print(arg8);

//            ofstream arg10("gmesh_final.txt");
//            gmesh->Print(arg10);

            TPZAnalysis an3(cmesh3);
            SolveSyst(an3, cmesh3);
            ofstream arg8("cmesh_HdivInicial.txt");
            cmesh3->Print(arg8);

            //TPZVec<REAL> erros;
            arg12<<" \nErro da HdivPressure  " <<endl;
            if (teste==1) an3.SetExact(*SolExata1);
            else an3.SetExact(*SolExata2);
            an3.PostProcess(erros, arg12);


            string plotile2("Solution_HDiv.vtk");
            PosProcessHDiv(an3, plotile2);
 */
            cmesh1->CleanUp();
            cmesh2->CleanUp();
            //mphysics->CleanUp();
            delete cmesh1;
            delete cmesh2;
            //delete mphysics;
            delete gmesh;
        }
    }
    return 0;
}

TPZGeoMesh *GMesh(bool triang_elements){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolTriang(3);
	TPZVec <long> TopolLine(2);
	
	//indice dos nos
	long id = 0;
	REAL valx, dx=1.;
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
    
    if(triang_elements==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    else{
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
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//        std::stringstream sout;
//        sout<<"\n\n Malha Geometrica Inicial\n ";
//        gmesh->Print(sout);
//        LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
   
	return gmesh;
}

TPZGeoMesh *GMesh2(bool triang_elements){
    TPZManVector<int,2> nx(2,2);
		nx[1] =1;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.5);
    x0[0] = -0.5;
    TPZGenGrid gengrid(nx,x0,x1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    if(triang_elements)
    {
      gengrid.SetElementType(1);
    }
    gengrid.Read(gmesh);
    TPZManVector<REAL,3> firstpoint(3,0.),secondpoint(3,0.);
    firstpoint[0] = 0.5;
    secondpoint[0] = 0.5;
    secondpoint[1] = 0.5;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc1);
    gengrid.SetBC(gmesh,6,bc2);
    gengrid.SetBC(gmesh,7,bc3);
    firstpoint[0] = -0.5;
    secondpoint[0] = 0.;
    secondpoint[1] = 0.;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc4);
    firstpoint = secondpoint;
    secondpoint[0] = 0.5;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc0);
#ifdef LOG4CXX
    if(logdata->isDebugEnabled())
    {
      std::stringstream sout;
      gmesh->Print(sout);
      LOGPZ_DEBUG(logdata,sout.str())
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
    //	gmesh->BuildConnectivity();
}


TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	TPZMaterial * mat(material);
	material->NStateVariables();
    
//    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
//	material->SetForcingFunction(force1);
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
	cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//        std::stringstream sout;
//        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
//        cmesh->Print(sout);
//        LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
	
	return cmesh;
}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
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
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
   
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
    
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
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
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_2 pressure\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
	return cmesh;
}



TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson * &mymaterial){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=1;
    int dim =2;
    
    REAL coefk=1.;
    mymaterial = new TPZMixedPoisson(MatId,dim);
    mymaterial->SetPermeability(coefk);
    
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    
    TPZAutoPointer<TPZFunction<STATE> > solExata;
    TPZAutoPointer<TPZFunction<STATE> > force;
    
    if(teste==1){
        solExata = new TPZDummyFunction<STATE>(SolExata1);
        mymaterial->SetForcingFunctionExact(solExata);
        
        force = new TPZDummyFunction<STATE>(Forcing1);
        mymaterial->SetForcingFunction(force);
        
        //Inserir condicoes de contorno
        TPZAutoPointer<TPZFunction<STATE> > fCC0;
        fCC0 = new TPZDummyFunction<STATE>(ForcingBC1);
        
        TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
        BCond0 = mymaterial->CreateBC(mat, bc0,dirichlet, val1, val2);
        BCond2 = mymaterial->CreateBC(mat, bc2,dirichlet, val1, val2);
        
        
        TPZFMatrix<STATE> val11(2,2,0.), val21(2,1,0.);
        BCond1 = mymaterial->CreateBC(mat, bc1,dirichlet, val11, val21);
        BCond3 = mymaterial->CreateBC(mat, bc3,dirichlet, val11, val21);
        BCond4 = mymaterial->CreateBC(mat, bc4,dirichlet, val11, val21);
        
        //BCond0->SetForcingFunction(fCC0);
       
    }
    else if(teste==2){
        solExata = new TPZDummyFunction<STATE>(SolExata2);
        mymaterial->SetForcingFunctionExact(solExata);
        
        TPZAutoPointer<TPZFunction<STATE> > fCC23;
        REAL fxy=8.;
        mymaterial->SetInternalFlux(fxy);
        fCC23 = new TPZDummyFunction<STATE>(ForcingBC2);
        
        //Inserir condicoes de contorno
        TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,1.);
        BCond0 = mymaterial->CreateBC(mat, bc0,dirichlet, val1, val2);
        BCond2 = mymaterial->CreateBC(mat, bc2,dirichlet, val1, val2);
        
        TPZFMatrix<STATE> val11(2,2,0.), val21(2,1,0.);
        BCond1 = mymaterial->CreateBC(mat, bc1,neumann, val11, val21);
        BCond3 = mymaterial->CreateBC(mat, bc3,neumann, val11, val21);
        BCond4 = mymaterial->CreateBC(mat, bc4,neumann, val11, val21);
        
        //BCond1->SetForcingFunction(fCC23);
        //BCond3->SetForcingFunction(fCC23);
        
    }
    else{
        solExata = new TPZDummyFunction<STATE>(SolExataSteklov);
        mymaterial->SetForcingFunctionExact(solExata);
        
        //Inserir condicoes de contorno
        TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
        BCond0 = mymaterial->CreateBC(mat, bc0, neumann, val1, val2);//bcmatNeumannZero
        BCond1 = mymaterial->CreateBC(mat, bc1, neumann, val1, val2);
        BCond2 = mymaterial->CreateBC(mat, bc2, neumann, val1, val2);
        BCond3 = mymaterial->CreateBC(mat, bc3, neumann, val1, val2);
        BCond4 = mymaterial->CreateBC(mat, bc4, dirichlet, val1, val2);//bcmaDirichletZero
        
        //Set force function
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannDireito;
        bcmatNeumannDireito = new TPZDummyFunction<STATE>(NeumannDireita);
        BCond1->SetForcingFunction(bcmatNeumannDireito);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
        bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcima);
        BCond2->SetForcingFunction(bcmatNeumannAcima);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEsquerdo;
        bcmatNeumannEsquerdo = new TPZDummyFunction<STATE>(NeumannEsquerda);
        BCond3->SetForcingFunction(bcmatNeumannEsquerdo);
        
    }//teste 3
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

void Forcing1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	//double x = pt[0];
    double y = pt[1];
    disp[0]= 4.*y*(1. - y) + 1.;
}

void SolExata2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &flux){
   // double x = pt[0];
    double y = pt[1];
    flux.Resize(3, 1);
    //disp[0] = 0.;
    disp[0]= 4.*y*(1. - y) + 1;
    flux(0,0)=0.;
    flux(1,0)=-(4.-8.*y);
    flux(2,0)=8;
}

void SolExata1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &flux){
    
    disp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    flux.Resize(3, 1);
    disp[0]= sin(Pi*x)*sin(Pi*y);
    flux(0,0)=-Pi*cos(Pi*x)*sin(Pi*y);
    flux(1,0)=-Pi*cos(Pi*y)*sin(Pi*x);
    flux(2,0)=2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}


void ForcingBC1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    //double y = 0.;
    disp[0]= Pi*sin(Pi*x);
    //disp[1]= 0.;
}

#include "pzbstrmatrix.h"
void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
	vecnames[0]  = "Flux";
    //vecnames[1]  = "GradFluxX";
    scalnames[0] = "Pressure";
    
    scalnames[1] = "ExactPressure";
    vecnames[1] = "ExactFlux";
			
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
    
}

void PosProcessHDiv(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(4), vecnames(2);
	scalnames[0] = "Pressure";
    scalnames[1] = "Divergence";
    scalnames[2] = "ExactPressure";
    scalnames[3] = "ExactDiv";
    
	vecnames[0]= "Flux";
    vecnames[1]= "ExactFlux";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}


TPZCompMesh *CMeshHDivPressure(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
    int MatId=1;
    int dim =2;
	TPZMatPoisson3d *material = new TPZMatPoisson3d(MatId,dim);
	material->NStateVariables();
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    REAL diff = 1.;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
    material->SetParameters(diff, conv, convdir);
    
    TPZAutoPointer<TPZFunction<STATE> > fCC0;
    TPZAutoPointer<TPZFunction<STATE> > force1;
    TPZAutoPointer<TPZFunction<STATE> > fCC23;
    
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    
    TPZAutoPointer<TPZFunction<STATE> > solExata;
    if(teste==1){
        solExata = new TPZDummyFunction<STATE>(SolExata1);
        force1 = new TPZDummyFunction<STATE>(Forcing1);
        material->SetForcingFunction(force1);
        material->SetForcingFunctionExact(solExata);
        // fCC0 = new TPZDummyFunction<STATE>(ForcingBC);
    }
    else{
        solExata = new TPZDummyFunction<STATE>(SolExata2);
        REAL fxy=8.;
        material->SetInternalFlux(fxy);
        material->SetForcingFunctionExact(solExata);
        fCC23 = new TPZDummyFunction<STATE>(ForcingBC2);
    }
    
    
    ///Inserir condicao de contorno
    
    if(teste==1){
        TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
        BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
        BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
        
        
        TPZFMatrix<STATE> val11(2,2,0.), val21(2,1,0.);
        BCond1 = material->CreateBC(mat, bc1,dirichlet, val11, val21);
        BCond3 = material->CreateBC(mat, bc3,dirichlet, val11, val21);
        
        //BCond0->SetForcingFunction(fCC0);
    }else{
        TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,1.);
        BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
        BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
        
        
        TPZFMatrix<STATE> val11(2,2,0.), val21(2,1,0.);
        BCond1 = material->CreateBC(mat, bc1,neumann, val11, val21);
        BCond3 = material->CreateBC(mat, bc3,neumann, val11, val21);
        
        //BCond1->SetForcingFunction(fCC23);
        //BCond3->SetForcingFunction(fCC23);
    }
    
    cmesh->SetAllCreateFunctionsHDivPressure();
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
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
    
//    if(IsZero(y) && IsZero(y)){
//        y=y+1.e-3;
//        //x=x+1.e-2;
//    }
    
    //flux = -k*grad(u), k=1 nesse problema
    du(0,0) = -pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = -pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    
}

void Compute_dudnQuadrado(TPZCompMesh *cmesh){
    std::map<REAL,REAL> dudn, dudnExato;
    const int npts = 100;
    const REAL rsize = 0.25;
    const REAL delta = rsize/(npts-1);
    std::map<REAL, std::pair<REAL,REAL> > sPath;
    std::map<REAL, std::pair<REAL,REAL> > sNormal;
    
    REAL sval;
    ///aresta y = 0 - tem tamanho 2*rsize por isso i < 2*npts
    for(int i = 1; i < 2*npts-2; i++){
        sval = delta/2.+i*delta;
        REAL x = rsize - sval;
        REAL y = 0;
        std::pair<REAL,REAL> normal(0,-1);
        sPath[sval] = std::make_pair(x,y);
        sNormal[sval] = normal;
    }
    
    ///aresta x = -rsize
    for(int i = 1; i < npts-1; i++){
        sval = 2*rsize +  i*delta;
        REAL x = -rsize;
        REAL y = 0 + i*delta;
        std::pair<REAL,REAL> normal(-1,0);
        sPath[sval] = std::make_pair(x,y);
        sNormal[sval] = normal;
    }
    
    ///aresta y = rsize - tem tamanho 2*rsize por isso i < 2*npts
    for(int i = 1; i < 2*npts-2; i++){
        sval = 3*rsize +  i*delta;
        REAL x = -rsize+i*delta ;
        REAL y = rsize;
        std::pair<REAL,REAL> normal(0,1);
        sPath[sval] = std::make_pair(x,y);
        sNormal[sval] = normal;
    }
    
    ///aresta x = +rsize
    for(int i = 1; i < npts-1; i++){
        sval = 5*rsize +  i*delta;
        REAL x = +rsize;
        REAL y = rsize - i*delta;
        std::pair<REAL,REAL> normal(1,0);
        sPath[sval] = std::make_pair(x,y);
        sNormal[sval] = normal;
    }

    std::map<REAL, std::pair<REAL,REAL> >::iterator wS, wNormal;
    int i;
    for(i = 0, wS = sPath.begin(), wNormal = sNormal.begin(); wS != sPath.end(); wS++, wNormal++, i++){
        const REAL s = wS->first;
        const REAL x = wS->second.first;
        const REAL y = wS->second.second;

        TPZManVector<REAL> xVec(3), qsi(2), normal(2);
        xVec[0] = x; xVec[1] = y; xVec[2] = 0.;
        normal[0] = wNormal->second.first;
        normal[1] = wNormal->second.second;
        
        TPZInterpolationSpace * sp = FindInterpolationSpace(xVec,cmesh,qsi);
        TPZManVector<REAL,2> flux(2,0.);
        TPZMaterialData data;
        
        sp->InitMaterialData(data);
        sp->ComputeShape(qsi,data);
        sp->ComputeSolution(qsi,data);
        
        flux[0]=data.sol[0][0];
        flux[1]=data.sol[0][1];
        
    
    
        REAL dudnval = flux[0]*normal[0] + flux[1]*normal[1];
        dudn[s] = dudnval;
        
        TPZManVector<REAL> uExato(1);
        TPZFNMatrix<100> duExato(2,1);
        SolExataSteklov(xVec, uExato, duExato);
        
        dudnExato[s] = duExato(0,0)*normal[0]+duExato(1,0)*normal[1];
    }///for i
    
    std::ofstream myfile("/Users/agnaldofarias/Documents/Fluxo-Steklov/dudnQuadradoHdiv.nb");
    PrintToFile(myfile, "dudnQuadradoHdiv", dudn);
    myfile << "\n";
    PrintToFile(myfile, "dudnQuadradoExatoHdiv", dudnExato);
    myfile << "\n\n";
    
    myfile << "ListPlot[{dudnQuadradoExatoHdiv, dudnQuadradoHdiv},PlotStyle->{Red,Black},Frame->True]\n\n";
    myfile << "R[v_] = Round[10*v]/10.;\n\n";
    myfile <<"Integrate[Interpolation[dudnQuadradoHdiv, InterpolationOrder -> 1][s], {s, 0, R[dudnQuadradoHdiv[[-1, 1]]]/12}]\n\n";
    myfile<<"Integrate[Interpolation[dudnQuadradoHdiv, InterpolationOrder -> 1][s], {s, 0,R[dudnQuadradoHdiv[[-1, 1]]]}]\n";
}///void


TPZInterpolationSpace * FindInterpolationSpace(TPZVec<REAL> &xVec, TPZCompMesh *cmesh, TPZVec<REAL> &qsi){
    const int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZInterpolationSpace * sp = dynamic_cast<TPZInterpolationSpace*>(cmesh->ElementVec()[iel]);
        if(!sp) continue;
        bool IsInDomain = sp->Reference()->ComputeXInverse(xVec, qsi, 1e-8);
        if(IsInDomain) return sp;
    }
    
    DebugStop();//nao achei ninguem
    return NULL;
    
}

void PrintToFile(std::ofstream &myfile, const std::string &title, std::map<REAL,REAL> &mymap){
    myfile << title << " = {";
    std::map<REAL,REAL>::const_iterator w;
    int n = mymap.size();
    int i;
    for(i = 0, w = mymap.begin(); w != mymap.end(); w++, i++){
        myfile << "{" << w->first << "," << w->second << "}";
        if(i != n-1) myfile << ",";
        else myfile << "};";
        if(i%4 == 0) myfile << "\n";
    }///for w
}

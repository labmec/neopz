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

#include "pzl2projection.h"

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
TPZGeoMesh *GMesh3();

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
REAL Compute_dudnQuadradoError(int ndiv, TPZCompMesh *cmesh,bool isquadradofechado);
REAL Compute_dudnQuadradoError(TPZCompMesh *cmesh);
bool EstaNoQuadrado(TPZGeoEl *gel, int ndiv);
TPZInterpolationSpace * FindInterpolationSpace(TPZVec<REAL> &xVec, TPZCompMesh *cmesh, TPZVec<REAL> &qsi);
void PrintToFile(std::ofstream &myfile, const std::string &title, std::map<REAL,REAL> &mymap);
void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out);
TPZCompMesh *L2Projection(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);

void DirichletSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannEsquerdaSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannDireitaSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannAcimaSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void SolExataSteklovSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
#endif

const bool triang = true;
const int teste = 3;
bool l2proj = false;
bool calcerroglobal = true;
int main5(int argc, char *argv[])
{
#ifdef LOG4CXX
////	std::string logs("../logmixedproblem.cfg");
////	InitializePZLOG("../logmixedproblem.cfg");
    InitializePZLOG();
#endif
    //gRefDBase.InitializeAllUniformRefPatterns();
    
    TPZVec<REAL> erros;
    ofstream arg12("ErroQuadrado.txt");
    int ndiv, p;
    arg12<<"\n\nERRO PARA MALHA COM TRIANGULO  = " << triang;
    for (p = 3; p< 4; p++)
    {
        arg12<<"\n\n----------- ORDEM p = " << p <<" -----------";
        for(ndiv = 1; ndiv < 8;ndiv++)
        {
            arg12<<"\nRefinamento ndiv  = " << ndiv;
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
                //gmesh = GMesh3();
    
            }
            UniformRefine(gmesh,ndiv);
            //ofstream arg1("gmesh_inicial.txt");
            //mesh->Print(arg1);

            if(l2proj==false)
            {
                // First computational mesh
                TPZCompMesh * cmesh1= CMeshFlux(gmesh, p);
//                ofstream arg2("cmesh1_inicial.txt");
//                cmesh1->Print(arg2);

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

//                // Cleaning reference of the geometric mesh to cmesh1
//                gmesh->ResetReference();
//                cmesh1->LoadReferences();
//                TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,ndiv,false);
//        //            ofstream arg4("cmesh1_final.txt");
//        //            cmesh1->Print(arg4);
//
//
//                // Cleaning reference to cmesh2
//                gmesh->ResetReference();
//                cmesh2->LoadReferences();
//                TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,ndiv,true);
//        //             ofstream arg5("cmesh2_final.txt");
//        //             cmesh2->Print(arg5);

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
                
                if(calcerroglobal==true){
                    
                    arg12<<"\n\nErro da simulacao multifisica  para o flux";
                    TPZAnalysis an1(cmesh1);
                    if (teste==1){
                        an1.SetExact(*SolExata1);
                    }
                    else if(teste==2){
                        an1.SetExact(*SolExata2);
                    }
                    else{
                        an1.SetExact(*SolExataSteklovSuave);
                    }
                    //an1.PostProcessError(erros, arg12);
                    ErrorHDiv(cmesh1,arg12);
                    
                    
                    arg12<<"\n\nErro da simulacao multifisica  para a pressao";
                    TPZAnalysis an2(cmesh2);
                    if (teste==1){
                        an2.SetExact(*SolExata1);
                    }
                    else if (teste==2){
                        an2.SetExact(*SolExata2);
                    }
                    else{
                        an2.SetExact(*SolExataSteklovSuave);
                    }
                    bool store_errors = false;
                    an2.PostProcessError(erros, store_errors, arg12);
                    
                }
                else
                {
                    //REAL errofluxo = Compute_dudnQuadradoError(ndiv,cmesh1,false);
                    REAL errofluxo = Compute_dudnQuadradoError(cmesh1);
                    arg12<<"\n\nErro L2 do fluxo = " << errofluxo<<"\n";
                    arg12.flush();
                }
                    
//            string plotfile("Solution_mphysics.vtk");
//            PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
                
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
            else
            {
                TPZVec<STATE> solini;
                TPZCompMesh  * cmeshL2 = L2Projection(gmesh, p, solini);
                TPZAnalysis an(cmeshL2);
                SolveSyst(an, cmeshL2);
                an.SetExact(*SolExataSteklov);
                REAL errofluxo = Compute_dudnQuadradoError(ndiv,cmeshL2,false);
                arg12<<"\n\nErro L2 da projecao L2f da solucao = " << errofluxo<<"\n";
                arg12.flush();
            }
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
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
	int64_t id = 0;
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

TPZGeoMesh *GMesh2(bool triang_elements)
{
    TPZManVector<int,2> nx(2,2);
    nx[1] =1;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.5);
    x0[0] = -0.5;
    TPZGenGrid gengrid(nx,x0,x1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    if(triang_elements)
    {
      gengrid.SetElementType(ETriangle);
    }
    gengrid.Read(gmesh);
    
    //elementos de contorno
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


TPZGeoMesh *GMesh3(){
    
    int Qnodes = 15;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx, valy, dx=0.25;
	for(int xi = 0; xi < Qnodes/3; xi++)
	{
		valx = -0.5 + xi*dx;
        valy = 0.;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx);//coord X
		Node[id].SetCoord(1 ,valy);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/3; xi++)
	{
		valx = -0.5 + xi*dx;
        valy = 0.25;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx);//coord X
		Node[id].SetCoord(1 ,valy);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
    for(int xi = 0; xi < Qnodes/3; xi++)
	{
		valx = -0.5 + xi*dx;
        valy = 0.5;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx);//coord X
		Node[id].SetCoord(1 ,valy);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    for(int i =0; i< 4; i++){
        TopolQuad[0] = i;
        TopolQuad[1] = i+1;
        TopolQuad[2] = 6+i;
        TopolQuad[3] = 5+i;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
    }

    for(int i =0; i< 4; i++){
        TopolQuad[0] = 5+i;
        TopolQuad[1] = 6+i;
        TopolQuad[2] = 11+i;
        TopolQuad[3] = 10+i;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
    }
    
    TopolLine[0] = 6;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
    id++;
    
    TopolLine[0] = 8;
    TopolLine[1] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
    id++;
    
    TopolLine[0] = 7;
    TopolLine[1] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
    id++;

//    TopolLine[0] = 1;
//    TopolLine[1] = 2;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
//    id++;
//    
//    TopolLine[0] = 2;
//    TopolLine[1] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
//    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
    id++;
    
    TopolLine[0] = 4;
    TopolLine[1] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
    id++;
    
    TopolLine[0] = 9;
    TopolLine[1] = 14;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
    id++;
    
    for(int i =0; i<4; i++){
        TopolLine[0] = 14-i;
        TopolLine[1] = 13-i;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
    }
    
    TopolLine[0] = 10;
    TopolLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    id++;
    
    TopolLine[0] = 5;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    
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
    
    cmesh->SetDimModel(dim);
    
    cmesh->InsertMaterialObject(mat);
	cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
	cmesh->SetDefaultOrder(pOrder);
    
	
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
    
    
#ifdef PZDEBUG
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
        solExata = new TPZDummyFunction<STATE>(SolExata1, 5);
        mymaterial->SetForcingFunctionExact(solExata);
        
        force = new TPZDummyFunction<STATE>(Forcing1, 5);
        mymaterial->SetForcingFunction(force);
        
        //Inserir condicoes de contorno
        TPZAutoPointer<TPZFunction<STATE> > fCC0;
        fCC0 = new TPZDummyFunction<STATE>(ForcingBC1, 5);
        
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
        solExata = new TPZDummyFunction<STATE>(SolExata2, 5);
        mymaterial->SetForcingFunctionExact(solExata);
        
        TPZAutoPointer<TPZFunction<STATE> > fCC23;
        REAL fxy=8.;
        mymaterial->SetInternalFlux(fxy);
        fCC23 = new TPZDummyFunction<STATE>(ForcingBC2, 5);
        
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
//        solExata = new TPZDummyFunction<STATE>(SolExataSteklov);
//        mymaterial->SetForcingFunctionExact(solExata);
        
        TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(SolExataSteklovSuave, 5);
        dum->SetPolynomialOrder(20);
        solExata = dum;
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
        bcmatNeumannDireito = new TPZDummyFunction<STATE>(NeumannDireitaSuave, 5);
        BCond1->SetForcingFunction(bcmatNeumannDireito);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
        bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcimaSuave, 5);
        BCond2->SetForcingFunction(bcmatNeumannAcima);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEsquerdo;
        bcmatNeumannEsquerdo = new TPZDummyFunction<STATE>(NeumannEsquerdaSuave, 5);
        BCond3->SetForcingFunction(bcmatNeumannEsquerdo);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet;
        bcmatDirichlet = new TPZDummyFunction<STATE>(DirichletSuave, 5);
        BCond4->SetForcingFunction(bcmatDirichlet);
        
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
        solExata = new TPZDummyFunction<STATE>(SolExata1, 5);
        force1 = new TPZDummyFunction<STATE>(Forcing1, 5);
        material->SetForcingFunction(force1);
        material->SetForcingFunctionExact(solExata);
        // fCC0 = new TPZDummyFunction<STATE>(ForcingBC, 5);
    }
    else{
        solExata = new TPZDummyFunction<STATE>(SolExata2, 5);
        REAL fxy=8.;
        material->SetInternalFlux(fxy);
        material->SetForcingFunctionExact(solExata);
        fCC23 = new TPZDummyFunction<STATE>(ForcingBC2, 5);
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
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {+1,0};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {0,+1};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
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
    
    if(IsZero(x) && IsZero(y)){
        y=y+1.e-3;
        x=x+1.e-2;
    }
    
    //flux = -k*grad(u), k=1 nesse problema
    du(0,0) = -pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = -pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
}

void SolExataSteklovSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    u.Resize(1, 0.);
    du.Resize(2, 1);
    du(0,0)=0.;
    du(1,0)=0.;
    
    REAL x = loc[0];
    REAL y = loc[1];
    //const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = 3.363585661014858*pow(pow(x,2) + pow(y,2),1.75)*cos(3.5*t);
    u[0] = sol;
    
    //flux = -k*grad(u), k=1 nesse problema
    du(0,0) = -pow(pow(x,2) + pow(y,2),0.75)*(11.772549813552002*x*cos(3.5*t) + 11.772549813552002*y*sin(3.5*t));
    du(1,0) = -pow(pow(x,2) + pow(y,2),0.75)*(11.772549813552002*y*cos(3.5*t) - 11.772549813552002*x*sin(3.5*t));
}

void DirichletSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklovSuave(loc,result,du);
}

void NeumannEsquerdaSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {-1,0};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
    SolExataSteklovSuave(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannDireitaSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {+1,0};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
    SolExataSteklovSuave(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannAcimaSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {0,+1};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
    SolExataSteklovSuave(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}


void Compute_dudnQuadrado(TPZCompMesh *cmesh){
    std::map<REAL,REAL> dudn, dudnExato, erro;
    const int npts = 100;
    const REAL rsize = 0.25;
    const REAL delta = rsize/(npts-1);
    std::map<REAL, std::pair<REAL,REAL> > sPath;
    std::map<REAL, std::pair<REAL,REAL> > sNormal;
    
    REAL sval;
    ///aresta y = 0 - tem tamanho 2*rsize por isso i < 2*npts
//    for(int i = 1; i < 2*npts-2; i++){
//        sval = delta/2.+i*delta;
//        REAL x = rsize - sval;
//        REAL y = 0;
//        std::pair<REAL,REAL> normal(0,-1);
//        sPath[sval] = std::make_pair(x,y);
//        sNormal[sval] = normal;
//    }
    
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
        
        TPZManVector<STATE> uExato(1);
        TPZFNMatrix<100,STATE> duExato(2,1);
        SolExataSteklov(xVec, uExato, duExato);
        
        dudnExato[s] = duExato(0,0)*normal[0]+duExato(1,0)*normal[1];
        
        erro[s] =  (dudn[s] -  dudnExato[s])*(dudn[s] -  dudnExato[s]);
    }///for i
    
    std::ofstream myfile("/Users/agnaldofarias/Documents/Fluxo-Steklov/dudnQuadradoHdiv.nb");
//    PrintToFile(myfile, "dudnQuadradoHdiv", dudn);
//    myfile << "\n";
    
    PrintToFile(myfile, "erroaoquadrado", erro);
    myfile << "\n";
    
//    PrintToFile(myfile, "dudnQuadradoExatoHdiv", dudnExato);
//    myfile << "\n\n";
    
    //myfile << "ListPlot[{dudnQuadradoExatoHdiv, dudnQuadradoHdiv},PlotStyle->{Red,Black},Frame->True]\n\n";
    myfile << "R[v_] = Round[10*v]/10.;\n\n";
    //myfile <<"Integrate[Interpolation[dudnQuadradoHdiv, InterpolationOrder -> 1][s], {s, 0, R[dudnQuadradoHdiv[[-1, 1]]]/12}]\n\n";
    //myfile<<"Integrate[Interpolation[dudnQuadradoHdiv, InterpolationOrder -> 1][s], {s, 0,R[dudnQuadradoHdiv[[-1, 1]]]}]\n";
    
    myfile <<"f = Interpolation[erroaoquadrado, InterpolationOrder -> 1]\n";
    myfile <<"Integrate[f[s], {s, f[[1,1,1]], f[[1,1,2]]}]\n\n";
}///void


///Gamma eh a uniao das regioes: 1) x=-0.25 e 0<=y<=0.25; 2) -0.25<=x<=0.25 e y=0.25; 3) x=0.25 e 0<=y<=0.25
REAL Compute_dudnQuadradoError(int ndiv, TPZCompMesh *cmesh, bool isquadradofechado){
    
    if (ndiv==0 || triang==true) {
        std::cout<<"ndiv nao pode ser menor que 1 \n";
        std::cout<<"considera-se apenas malha quadrilateral\n";
        DebugStop();
    }
    
    REAL error = 0.;
    const REAL rsize = 0.25;
    const int ninterf = pow((REAL)2,ndiv-1);

    const REAL delta = rsize/ninterf;
    
    TPZStack<TPZManVector<REAL> > coordX;
    TPZStack<TPZManVector<REAL> > solflux;
    TPZStack<REAL> soldudn;

    TPZManVector<REAL,3> xVec(3), faceNormal(2);
    xVec[0] = 0.; xVec[1] = 0.; xVec[2] = 0.;

    ///Regiao 1): x = -rsize e 0=<y<=rsize
    faceNormal[0] = -1.;
    faceNormal[1] = 0.;
    for(int i = 0; i < ninterf; i++)
    {
        xVec[0] = -rsize + 1.e-3;
        xVec[1] = (0. + delta/2.) + i*delta;
        
        TPZManVector<REAL,3> qsi(2);
        TPZInterpolationSpace * sp = FindInterpolationSpace(xVec,cmesh,qsi);
        if(!sp) DebugStop();
        
        TPZGeoEl *gel = sp->Reference();
        int dim = gel->Dimension();
        //verificando se o elemento esta no quadrado
        if(EstaNoQuadrado(gel, ndiv) == false) DebugStop();
        
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(7, 2);
        TPZManVector<int,3> prevorder(dim), maxorder(dim,intrule->GetMaxOrder());
        //intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        TPZGeoElSide geoside(gel,7);
        
        TPZMaterialData data;
        sp->InitMaterialData(data);
        
        const int npoints = intrule->NPoints();
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            qsi[0]=0.; qsi[1]=0.;
            intrule->Point(ip,qsi,weight);
            TPZTransform<> tr;
            TPZGeoElSide geosideh(gel,8);
            tr = geoside.SideToSideTransform(geosideh);
            TPZManVector<REAL,3> qsih(2);
            tr.Apply(qsi, qsih);
            geoside.Jacobian(qsi, data.jacobian, data.axes, data.detjac, data.jacinv);
            weight *= fabs(data.detjac);
            sp->ComputeShape(qsih, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            //TPZFNMatrix<660> dphi(data.dphix.Rows(), data.dphix.Cols(), 0.);
            //sp->Shape(qsih,data.phi,dphi);
            //sp->Convert2Axes(dphi, data.jacinv, data.dphix);
            sp->ComputeSolution(qsih,data);
            
            TPZManVector<REAL,2> flux(2,0.);
            flux[0]=data.sol[0][0];
            flux[1]=data.sol[0][1];
            REAL dudnval = flux[0]*faceNormal[0] + flux[1]*faceNormal[1];
            
            geoside.X(qsi,xVec);
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            SolExataSteklovSuave(xVec, uExato, duExato);
            REAL dudnExato = duExato(0,0)*faceNormal[0]+duExato(1,0)*faceNormal[1];
            
            error += weight*(dudnval - dudnExato)*(dudnval - dudnExato);
            
            coordX.push_back(xVec);
            solflux.push_back(flux);
            soldudn.push_back(dudnval);
        }
        //intrule->SetOrder(prevorder);
    }
    
    ///Regiao 2): y = rsize e -rsize <=x<=rsize
    faceNormal[0] = 0.;
    faceNormal[1] = 1.;
    for(int i = 0; i < 2*ninterf; i++)
    {
        xVec[0] = (-rsize + delta/2.) + i*delta;
        xVec[1] = rsize - 1.e-3;
        
        TPZManVector<REAL,3> qsi(2);
        TPZInterpolationSpace * sp = FindInterpolationSpace(xVec,cmesh,qsi);
        if(!sp) DebugStop();
        
        TPZGeoEl *gel = sp->Reference();
        int dim = gel->Dimension();
        //verificando se o elemento esta no quadrado
        if(EstaNoQuadrado(gel, ndiv) == false) DebugStop();

        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(6, 2);
        TPZManVector<int,3> prevorder(dim), maxorder(dim,intrule->GetMaxOrder());
        //intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        TPZGeoElSide geoside(gel,6);
  
        TPZMaterialData data;
        sp->InitMaterialData(data);
        
        const int npoints = intrule->NPoints();
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            qsi[0]=0.; qsi[1]=0.;
            intrule->Point(ip,qsi,weight);
            TPZTransform<> tr;
            TPZGeoElSide geosideh(gel,8);
            tr = geoside.SideToSideTransform(geosideh);
            TPZManVector<REAL,3> qsih(2);
            tr.Apply(qsi, qsih);
            geoside.Jacobian(qsi, data.jacobian, data.axes, data.detjac, data.jacinv);
            weight *= fabs(data.detjac);
            sp->ComputeShape(qsih, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            //TPZFNMatrix<660> dphi(data.dphix.Rows(), data.dphix.Cols(), 0.);
            //sp->Shape(qsih,data.phi,dphi);
            //sp->Convert2Axes(dphi, data.jacinv, data.dphix);
           
            sp->ComputeSolution(qsih,data);
            
            TPZManVector<REAL,2> flux(2,0.);
            flux[0]=data.sol[0][0];
            flux[1]=data.sol[0][1];
            REAL dudnval = flux[0]*faceNormal[0] + flux[1]*faceNormal[1];
            
            geoside.X(qsi,xVec);
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            SolExataSteklov(xVec, uExato, duExato);
            REAL dudnExato = duExato(0,0)*faceNormal[0]+duExato(1,0)*faceNormal[1];
            
            error += weight*(dudnval - dudnExato)*(dudnval - dudnExato);
            
            coordX.push_back(xVec);
            solflux.push_back(flux);
            soldudn.push_back(dudnval);
        }
        //intrule->SetOrder(prevorder);
    }
    
    ///Regiao 3): x = +rsize e 0<=y<=rsize
    faceNormal[0] = 1.;
    faceNormal[1] = 0.;
    for(int i = 0; i < ninterf; i++)
    {
        xVec[0] = rsize - 1.e-3;
        xVec[1] = (0. + delta/2.) + i*delta;
        
        TPZManVector<REAL,3> qsi(2);
        TPZInterpolationSpace * sp = FindInterpolationSpace(xVec,cmesh,qsi);
        if(!sp) DebugStop();
        
        TPZGeoEl *gel = sp->Reference();
        int dim = gel->Dimension();
        //verificando se o elemento esta no quadrado
        if(EstaNoQuadrado(gel, ndiv) == false) DebugStop();

        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(5, 2);
        TPZManVector<int,3> prevorder(dim), maxorder(dim,intrule->GetMaxOrder());
        //intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        TPZGeoElSide geoside(gel,5);
        
        TPZMaterialData data;
        sp->InitMaterialData(data);
        
        const int npoints = intrule->NPoints();
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            qsi[0]=0.; qsi[1]=0.;
            intrule->Point(ip,qsi,weight);
            TPZTransform<> tr;
            TPZGeoElSide geosideh(gel,8);
            tr = geoside.SideToSideTransform(geosideh);
            TPZManVector<REAL,3> qsih(2);
            tr.Apply(qsi, qsih);
            geoside.Jacobian(qsi, data.jacobian, data.axes, data.detjac, data.jacinv);
            weight *= fabs(data.detjac);
            sp->ComputeShape(qsih, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            //TPZFNMatrix<660> dphi(data.dphix.Rows(), data.dphix.Cols(), 0.);
            //sp->Shape(qsih,data.phi,dphi);
            //sp->Convert2Axes(dphi, data.jacinv, data.dphix);

            sp->ComputeSolution(qsih,data);
            
            TPZManVector<REAL,2> flux(2,0.);
            flux[0]=data.sol[0][0];
            flux[1]=data.sol[0][1];
            REAL dudnval = flux[0]*faceNormal[0] + flux[1]*faceNormal[1];
            
            geoside.X(qsi,xVec);
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            SolExataSteklov(xVec, uExato, duExato);
            REAL dudnExato = duExato(0,0)*faceNormal[0]+duExato(1,0)*faceNormal[1];
            
            error += weight*(dudnval - dudnExato)*(dudnval - dudnExato);
            
            coordX.push_back(xVec);
            solflux.push_back(flux);
            soldudn.push_back(dudnval);
        }
        //intrule->SetOrder(prevorder);
    }
    
    if(isquadradofechado==true){
        ///Regiao 4): -rsize =<x <= rsize e y=0
        faceNormal[0] = 0.;
        faceNormal[1] = -1.;
        for(int i = 0; i < 2*ninterf; i++)
        {
            xVec[0] = (-rsize + delta/2.) + i*delta;
            //xVec[0] = (0 + delta/2.) + i*delta;
            xVec[1] = 0 + 1.e-3;
            
            TPZManVector<REAL,3> qsi(2);
            TPZInterpolationSpace * sp = FindInterpolationSpace(xVec,cmesh,qsi);
            if(!sp) DebugStop();
            
            TPZGeoEl *gel = sp->Reference();
            int dim =gel->Dimension();
            //verificando se o elemento esta no quadrado
            if(EstaNoQuadrado(gel, ndiv) == false) DebugStop();
            
            TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(4, 2);
            TPZManVector<int,3> prevorder(dim), maxorder(dim,intrule->GetMaxOrder());
            intrule->GetOrder(prevorder);
            intrule->SetOrder(maxorder);
            TPZGeoElSide geoside(gel,4);
            
            TPZMaterialData data;
            sp->InitMaterialData(data);
            
            const int npoints = intrule->NPoints();
            for(int ip = 0; ip < npoints; ip++)
            {
                REAL weight;
                qsi[0]=0.; qsi[1]=0.;
                intrule->Point(ip,qsi,weight);
                TPZTransform<> tr;
                TPZGeoElSide geosideh(gel,8);
                tr = geoside.SideToSideTransform(geosideh);
                TPZManVector<REAL,3> qsih(2);
                tr.Apply(qsi, qsih);
                geoside.Jacobian(qsi, data.jacobian, data.axes, data.detjac, data.jacinv);
                weight *= fabs(data.detjac);
                sp->ComputeShape(qsih, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
                //TPZFNMatrix<660> dphi(data.dphix.Rows(), data.dphix.Cols(), 0.);
                //sp->Shape(qsih,data.phi,dphi);
                //sp->Convert2Axes(dphi, data.jacinv, data.dphix);
            
                sp->ComputeSolution(qsih,data);
                
                 TPZManVector<REAL,2> flux(2,0.);
                flux[0]=data.sol[0][0];
                flux[1]=data.sol[0][1];
                REAL dudnval = flux[0]*faceNormal[0] + flux[1]*faceNormal[1];
                
                geoside.X(qsi,xVec);
                TPZManVector<STATE> uExato(1);
                TPZFNMatrix<100,STATE> duExato(2,1);
                SolExataSteklov(xVec, uExato, duExato);
                REAL dudnExato = duExato(0,0)*faceNormal[0]+duExato(1,0)*faceNormal[1];
                
                error += weight*(dudnval - dudnExato)*(dudnval - dudnExato);
                
                coordX.push_back(xVec);
                solflux.push_back(flux);
                soldudn.push_back(dudnval);
            }
            
            intrule->SetOrder(prevorder);
        }

    }
    
#ifdef LOG4CXX
    if(logdata->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\nResultados com SideToSideTransform\n";
        int npts = coordX.size();
        for(int i=0;i<npts; i++) {
            sout << "(x,y) = (" << coordX[i][0] << ", " << coordX[i][1]<<")  e   " << "fluxo = (" << solflux[i][0] <<", " <<solflux[i][1] << ")\n";
        }
        LOGPZ_DEBUG(logdata,sout.str())
    }
#endif
    
#ifdef LOG4CXX
    if(logdata->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\nSolSideToSide = {";
        int npts = coordX.size();
        for(int i=0;i<npts-1; i++) {
            sout << "{" << coordX[i][0]<< ", " << coordX[i][1]<<", " << soldudn[i] <<"},";
        }
        sout << "{" << coordX[npts-1][0]<< ", " << coordX[npts-1][1]<<", " << soldudn[npts-1] <<"}";
        sout << "}";
        LOGPZ_DEBUG(logdata,sout.str())
    }
#endif

    
    error = sqrt(error);
    return error;
    
}///method

REAL Compute_dudnQuadradoError(TPZCompMesh *cmesh){
    
    TPZStack<TPZManVector<REAL,2> > coordX;
    TPZStack<TPZManVector<REAL,2> > solflux;
    TPZStack<REAL> soldudn;
    
    REAL error = 0.;
    int nel = cmesh->NElements();
    int iel;
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        int dim = cel->Dimension();
        if(dim != 1) continue;
        int matid = cel->Material()->Id();
        if(matid != matId) continue;
        
        TPZGeoEl *gel = cel->Reference();
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 2);
        TPZManVector<int,3> prevorder(dim), maxorder(dim,intrule->GetMaxOrder());
        intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        const int npoints = intrule->NPoints();
        TPZManVector<REAL,3> qsi(1), xVec(3), faceNormal(2);
        faceNormal[0]=0.;
        faceNormal[1]=0.;
        gel->CenterPoint(gel->NSides()-1, qsi);
        gel->X(qsi,xVec);
        if(xVec[0]==-0.25) faceNormal[0]=-1.;
        if(xVec[0]== 0.25) faceNormal[0]=1.;
        if(xVec[1]== 0.25) faceNormal[1]=1.;
        //            if(data.x[1]== 0.0) data.normal[1] = -1.;
        
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            intrule->Point(ip,qsi,weight);
            sp->ComputeShape(qsi, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            weight *= fabs(data.detjac);
            sp->ComputeSolution(qsi,data);
            
            REAL sign = faceNormal[0] + faceNormal[1];
            
            TPZManVector<REAL,2> flux(2,0.);
            flux[0]=data.sol[0][0];
            REAL dudnval = flux[0]*sign;
            
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            gel->X(qsi,xVec);
            SolExataSteklovSuave(xVec, uExato, duExato);
            REAL dudnExato = duExato(0,0)*faceNormal[0]+duExato(1,0)*faceNormal[1];
            
            error += weight*(dudnval - dudnExato)*(dudnval - dudnExato);
            
            coordX.push_back(xVec);
            solflux.push_back(flux);
            soldudn.push_back(dudnval);
        }
        intrule->SetOrder(prevorder);
    }
    
#ifdef LOG4CXX
    if(logdata->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\nResultados com elemento 1D\n";
        int npts = coordX.size();
        for(int i=0;i<npts; i++) {
            sout << "(x,y) = (" << coordX[i][0] << ", " << coordX[i][1]<<")  e   " << "fluxo = " << solflux[i][0] <<"\n";
        }
        LOGPZ_DEBUG(logdata,sout.str())
    }
#endif
    
    
#ifdef LOG4CXX
    if(logdata->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\nSol1D = {";
        int npts = coordX.size();
        for(int i=0;i<npts-1; i++) {
            sout << "{" << coordX[i][0]<< ", " << coordX[i][1]<<", " << soldudn[i] <<"},";
        }
        sout << "{" << coordX[npts-1][0]<< ", " << coordX[npts-1][1]<<", " << soldudn[npts-1] <<"}";
        sout << "}";
        LOGPZ_DEBUG(logdata,sout.str())
    }
#endif
    
    error = sqrt(error);
    return error;
    
}///method

bool EstaNoQuadrado(TPZGeoEl *gel, int ndiv){
    
    if(gel->Index()==9 || gel->Index()==12){
        return true;
    }
    
    int nfather = ndiv;
    TPZGeoEl *gelf = gel;
    while (nfather > 1) {
        gelf = gelf->Father();
        nfather--;
    }
    
    if(gelf->Index()==9 || gelf->Index()==12){
        return true;
    }
    
    return false;
    
}///bool

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

TPZCompMesh *L2Projection(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(matId, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(SolExataSteklov, 5);
    material->SetForcingFunction(forcef);
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(triang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }
    
	return cmesh;
}

void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExataSteklovSuave, elerror, 0);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<endl;
    
}



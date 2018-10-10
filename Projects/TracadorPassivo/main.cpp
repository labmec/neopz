
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <math.h>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPBSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZCompElDisc.h"
#include "pzpoisson3d.h"
#include "pzconvectionproblem.h"
#include "mixedpoisson.h"
#include "pztracerflow.h"
#include "pzl2projection.h"

#include "TPZEquationFilter.h"

#include "pzgradientreconstruction.h"

#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pztransfer.h"

#include "pzlog.h"

using namespace std;

int matId =1;
int matIdL2Proj = 2;

int dirichlet = 0;
int neumann = 1;
int inflow = 3;
int outflow = 4;
int neumanninflow = 5;
int dirichletoutflow = 6;
int dirichletinflow = 7;

int bc0=-1;
int bc1=-2;
int bc2=-3;
int bc3=-4;


TPZGeoMesh *GMesh(bool triang_elements, REAL Lx, REAL Ly);
TPZGeoMesh *GMesh2(REAL Lx, REAL Ly,bool triang_elements);

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshSaturation(TPZGeoMesh * gmesh, int pOrder);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);


TPZCompMesh *SetCondicaoInicial(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh);


void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);

void PosProcessFlux(TPZAnalysis &an, std::string plotfile);

TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZTracerFlow * mymaterial, TPZCompMesh* mphysics);

void StiffMatrixLoadVec(TPZTracerFlow *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZAutoPointer< TPZMatrix<REAL> > &matK1, TPZFMatrix<STATE> &fvec);

void ResolverComReconstGradiente(REAL deltaX,REAL maxTime,TPZManVector<TPZCompMesh *,3> meshvec, TPZCompMesh* mphysics, TPZGradientReconstruction *gradreconst, bool useRK2);

void ResolverSemReconstGradiente(REAL deltaX,REAL maxTime,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);

void FilterPressureFluxEquation(TPZTracerFlow *mymaterial, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics, TPZAnalysis &an);

void FilterSaturationEquation(TPZTracerFlow *mymaterial, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics, TPZAnalysis &an, bool currentstate);

void ForcingInicial(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Permeability(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void CondCFL(TPZFMatrix<STATE> SolutionQ, REAL deltaX, REAL maxTime, REAL &deltaT);

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
void CreateInterfaceMultiphysicsElement(TPZCompMesh* mphysics);
void CreateInterfaceElement(TPZCompMesh* cmesh);



bool ftriang = false;
bool fishomogmedium = false;
int fhetertest = 2;
bool recgrad = true;
bool useRK2 = false;
bool SaturationMeshFine = true;

REAL ftimeatual = 0.;

REAL  fLx = 400.;//m
REAL  fLy = 100.;//m
REAL fk1 =  9.86923e-13;//m2
REAL fk2 = 0.;

REAL fvazaoentrada = 60.;//m3/d
REAL fporos = 0.3;

REAL fLref = fLy;
REAL fkref = fk1;//(fk1+fk2)/2.;
REAL fpref = 1.e7;//1.96133e7;//pa
REAL fvisref = 1.e-3;

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

int main(int argc, char *argv[])
{
//#ifdef LOG4CXX
//    InitializePZLOG();
//#endif
    
    if(fishomogmedium == true)
    {
        fk2 = fk1;
    }
    else{
        fk2 = 9.86923e-16;//m2
    }
    
    REAL Lx = fLx/fLref;
    REAL Ly = fLy/fLref;
    
    int pq = 1;
    int pp;
    int ps=1;
    if(ftriang==true){
        pp = pq-1;
    }else{
        pp = pq;
    }

    int h = 4;//4
    //TPZGeoMesh *gmesh = GMesh(ftriang, Lx, Ly);
    TPZGeoMesh *gmesh = GMesh2(Lx, Ly,ftriang);
    UniformRefine(gmesh,h);
    
//    ofstream arg("gmesh1.txt");
//    gmesh->Print(arg);
//    {
//        std::ofstream filemesh("MalhaGeometricaInicial.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
//    }
    
    TPZCompMesh *cmesh1 = CMeshSaturation(gmesh, ps);
    TPZCompMesh *cmesh2 = CMeshFlux(gmesh, pq);
    TPZCompMesh *cmesh3 = CMeshPressure(gmesh, pp);
    {
        ofstream out("meshflux.txt");
        cmesh2->Print(out);
        
        ofstream out2("meshpressure.txt");
        cmesh3->Print(out2);
        
        ofstream arg1("cmeshsaturation1.txt");
        cmesh1->Print(arg1);
    }

    
    //Malha da saturacao
//    TPZCompMesh *cmesh1 = CMeshSaturation(gmesh, ps);
//    gmesh->ResetReference();
//    cmesh1->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,h,true);
//    cmesh1->AutoBuild();
//    cmesh1->AdjustBoundaryElements();
//    cmesh1->CleanUpUnconnectedNodes();
////    {
////        ofstream arg1("cmeshsaturation1.txt");
////        cmesh1->Print(arg1);
////        std::ofstream filemesh("MalhaGeometricaInicial1.vtk");
////        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
////    }
//    
//    
//    //Malha do Fluxo
//    TPZCompMesh *cmesh2 = CMeshFlux(gmesh, pq);
//    gmesh->ResetReference();
//    cmesh2->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,h,false);
//    cmesh2->AutoBuild();
//    cmesh2->AdjustBoundaryElements();
//    cmesh2->CleanUpUnconnectedNodes();
////    {
////        ofstream out("meshflux.txt");
////        cmesh2->Print(out);
////    }
////    
//    //Malha da Pressao
//    TPZCompMesh *cmesh3 = CMeshPressure(gmesh, pp);
//    gmesh->ResetReference();
//    cmesh3->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh3,h,true);
//    cmesh3->AutoBuild();
//    cmesh3->AdjustBoundaryElements();
//    cmesh3->CleanUpUnconnectedNodes();
//    {
//        ofstream out("meshpressure.txt");
//        cmesh3->Print(out);
//    }
    
    //------Malha da saturacao mais fina ----------------------------
    if(SaturationMeshFine)
    {
        gmesh->ResetReference();
        cmesh1->LoadReferences();
        TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,1,true);
        cmesh1->AdjustBoundaryElements();
        cmesh1->CleanUpUnconnectedNodes();
        
        ofstream out("cmeshsaturation2.txt");
        cmesh1->Print(out);
    }

    //inserir interface
    CreateInterfaceElement(cmesh1);
    {
        ofstream out("cmeshsaturation2ComInterface.txt");
        cmesh1->Print(out);
//        std::ofstream filemesh("MalhaGeometricaInicial2.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
    }

    
////    ofstream arg1("cmeshsaturation.txt");
////    cmesh1->Print(arg1);
    
    /*
    //-----------------------------------------------------------------------
    //Set initial conditions for saturation
    TPZAnalysis an0(cmesh1);
    int nrs = an0.Solution().Rows();
    TPZVec<STATE> solini(nrs,1);
    TPZCompMesh  * cmeshL2 = SetCondicaoInicial(gmesh, ps, solini);
    //ofstream arg6("cmeshL2.txt");
    //cmeshL2->Print(arg6);
    
    TPZAnalysis anL2(cmeshL2);
    SolveSyst(anL2, cmeshL2);
    an0.LoadSolution(anL2.Solution());
    //-----------------------------------------------------------------------
     

    
    // Cleaning reference to cmesh2
    gmesh->ResetReference();
    cmesh2->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,1,false);
    cmesh2->AdjustBoundaryElements();
    cmesh2->CleanUpUnconnectedNodes();
    //ofstream arg2("cmeshflux.txt");
    //cmesh2->Print(arg2);
    
    // Cleaning reference to cmesh3
    gmesh->ResetReference();
    cmesh3->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh3,1,true);
    cmesh3->AdjustBoundaryElements();
    cmesh3->CleanUpUnconnectedNodes();
    //ofstream arg3("cmeshpressure.txt");
    //cmesh3->Print(arg3);
     */
    

//    ofstream arg1("cmeshsaturation.txt");
//    cmesh1->Print(arg1);
//    
//    ofstream arg2("cmeshflux.txt");
//    cmesh2->Print(arg2);
//    
//    ofstream arg3("cmeshpressure.txt");
//    cmesh3->Print(arg3);
    
    //ofstream arg4("gmesh2.txt");
    //gmesh->Print(arg4);

//    //-----------------------------------------------------------------------
    //Set initial conditions for saturation
    TPZAnalysis an0(cmesh1);
    int nrs = an0.Solution().Rows();
    TPZVec<STATE> solini(nrs,1);
    TPZCompMesh  * cmeshL2 = SetCondicaoInicial(gmesh, ps, solini);
    //ofstream arg6("cmeshL2.txt");
    //cmeshL2->Print(arg6);
    
    TPZAnalysis anL2(cmeshL2);
    SolveSyst(anL2, cmeshL2);
    an0.LoadSolution(anL2.Solution());
    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    //-----------------------------------------------------------------------
    
    //malha multifisica
    TPZManVector<TPZCompMesh *,3> meshvec;
    meshvec.Resize(3);
    meshvec[0] = cmesh1;
    meshvec[1] = cmesh2;
    meshvec[2] = cmesh3;
    
    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
    {
        ofstream arg5("cmeshmultiphysics.txt");
        mphysics->Print(arg5);
    }
    
//    ofstream arg7("gmesh3.txt");
//    gmesh->Print(arg7);
//    ofstream file("malhageometrica.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file, true);
//
//    ofstream arg8("cmeshvec0.txt");
//    meshvec[0]->Print(arg8);
    
    REAL VolPoros = fLx*fLy*1.*fporos;
    REAL tempoDia = VolPoros/fvazaoentrada;
    REAL tempoSeg = tempoDia*(24.*60.*60.);
    
    REAL num = fkref*fpref;
    REAL denom = fvisref*fLref*fLref;
    
    REAL tempoAdimens = (num/denom)*tempoSeg;
    int temp = (int)tempoAdimens*fporos;
    REAL tD = (temp + 1.);

    REAL deltaX = Ly/(pow(2.,h));

    if(recgrad==true)
    {
        TPZGradientReconstruction *gradreconst = new TPZGradientReconstruction(false,1.);
        
        //if(ftriang==true){
            TPZVec<REAL> LxLyLz(2,0.);
            LxLyLz[0] = Lx; LxLyLz[1]=Ly;
            
            TPZVec<int> MatIdBC(2,0);
            MatIdBC[0] = -1; MatIdBC[1] = -3;
            
            TPZVec<REAL>  Xmin(2,0.);
            TPZVec<REAL> Xmax(2,0.);
            TPZManVector<TPZVec<REAL> > coordmin(2,0.);
            TPZManVector<TPZVec<REAL> > coordmax(2,0.);
            Xmin[0] = 0.; Xmin[1] = 0.; Xmax[0] = Lx; Xmax[1] = 0.;
            coordmin[0] = Xmin; coordmax[0] = Xmax;
            Xmin[0] = 0.; Xmin[1] = Ly; Xmax[0] = Lx; Xmax[1] = Ly;
            coordmin[1] = Xmin; coordmax[1] = Xmax;
            
            gradreconst->SetDataGhostsNeighbors(LxLyLz, MatIdBC, coordmin, coordmax);
        //}
    
        ResolverComReconstGradiente(deltaX,tD,meshvec, mphysics,gradreconst, useRK2);
    }else{
        ResolverSemReconstGradiente(deltaX,tD,meshvec, mphysics);
    }
    
	return EXIT_SUCCESS;
}

TPZGeoMesh *GMesh(bool triang_elements, REAL Lx, REAL Ly){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    if(triang_elements==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 3;
        TopolTriang[2] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
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
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
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

TPZGeoMesh *GMesh2(REAL Lx, REAL Ly, bool triang_elements){
    
    int Qnodes = 10;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx;
    int nelem = Qnodes/2 - 1;
    REAL dx = Lx/nelem;
    int hnodes = Qnodes/2;
    
	for(int xi = 0; xi < hnodes; xi++)
	{
		valx = xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < hnodes; xi++)
	{
		valx = Lx - xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    //elementos internos
    if(triang_elements==false)
    {
        for(int i = 0; i< nelem; i++)
        {
            TopolQuad[0] = i;
            TopolQuad[1] = i+1;
            TopolQuad[2] = (Qnodes-2) - i;
            TopolQuad[3] = (Qnodes-1) - i;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
            id++;
        }
        
    }
    else
    {
        for(int i = 0; i< nelem; i++)
        {
            TopolTriang[0] = i;
            TopolTriang[1] = i+1;
            TopolTriang[2] = (Qnodes-2) - i;
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
            id++;
            
            TopolTriang[0] = i;
            TopolTriang[1] = (Qnodes-2) - i;
            TopolTriang[2] = (Qnodes-1) - i;
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
            id++;
        }
    }
    
    //elemsentos do contorno
    for(int i = 0; i< nelem; i++)
    {
        TopolLine[0] = i;
        TopolLine[1] = i+1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
    }
    
    TopolLine[0] = TopolLine[1];
    TopolLine[1] = TopolLine[0]+1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
    id++;
    
    for(int i = 0; i< nelem; i++)
    {
        TopolLine[0] = hnodes + i;
        TopolLine[1] = (hnodes+1) + i ;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
    }
    
    TopolLine[0] = TopolLine[1];
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
}

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
	TPZMaterial * mat(material);
	material->NStateVariables();
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond1 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
	cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    //Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}


TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim);
	material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond1 = material->CreateBC(mat,bc0,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat,bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat,bc2,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat,bc3,dirichlet, val1, val2);
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
    ///inserir connect da pressao
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
	    newnod.SetLagrangeMultiplier(1);
    }
    
    ///set order total da shape
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
            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
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

    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}

TPZCompMesh *CMeshSaturation(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	
	TPZMatConvectionProblem *material = new TPZMatConvectionProblem(matId,dim);
	TPZMaterial * mat(material);
	
	TPZVec<REAL> convdir(dim,0.);
    convdir[0]=1.;
	REAL flux = 0.;
    REAL rho = 1.;
	
	material->SetParameters(rho,convdir);
	material->SetInternalFlux(flux);
	material->NStateVariables();
	
	//TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(mat);
    
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1 = material->CreateBC(mat,bc0,neumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat,bc1,outflow, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat,bc2,neumann, val1, val2);
    REAL uD =2.;
    val2(0,0) = uD;
	TPZMaterial * BCond4 = material->CreateBC(mat, bc3,inflow, val1, val2);
    
    //val2(0,0) = 0.;
    //TPZMaterial * BCond1 = material->CreateBC(mat, bc1,outflow, val1, val2);
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,material->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    //cmesh->CleanUpUnconnectedNodes();
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }
    
//    cmesh->CleanUpUnconnectedNodes();
//    cmesh->ExpandSolution();
//    
//    cmesh->LoadReferences();
//    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
//    cmesh->ExpandSolution();
//    cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}

TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int dim =2;
    

    REAL visc = 1.e-3;//pa.s
    REAL viscD = visc/fvisref;
    
    REAL AreaPoro = fLy*1.*fporos;
    REAL qinDia = fvazaoentrada/AreaPoro;//m/d
    REAL qin = qinDia/86400.;//m/s
    REAL num = fvisref*fLref;
    REAL denom = fkref*fpref;
    REAL qinD = (-1.)*(qin*num/denom);
    
    REAL inflow = 1.;
    
    TPZTracerFlow *material1 = new TPZTracerFlow(matId,dim);
    
    material1->SetViscosity(viscD);
    material1->SetPorosity(fporos);
    
    //funcao que define a permeabilidade
    if(fishomogmedium==false){
        TPZAutoPointer<TPZFunction<STATE> > forcefK;
        TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Permeability, 5);
        dum->SetPolynomialOrder(10);
        forcefK = dum;
        material1->SetForcingFunction(forcefK);
    }
    
    if(fishomogmedium==true){
        TPZAutoPointer<TPZFunction<STATE> > solExata;
        TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(SolExata, 5);
        dum->SetPolynomialOrder(10);
        solExata = dum;
        material1->SetForcingFunctionExact(solExata);
    }
    
    TPZMaterial *mat1(material1);
    mphysics->InsertMaterialObject(mat1);
    mphysics->SetDimModel(dim);
    
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    BCond0 = material1->CreateBC(mat1, bc0,neumann, val1, val2);
    
    val2(1,0) = 0.;//0.01019716e-3;
    BCond1 = material1->CreateBC(mat1, bc1,dirichletoutflow, val1, val2);
    
    val2(1,0) = 0.;
    BCond2 = material1->CreateBC(mat1, bc2,neumann, val1, val2);
    
    val2(0,0)= inflow;
    val2(1,0)= qinD;
    BCond3 = material1->CreateBC(mat1, bc3,neumanninflow, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);

    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
  
    
    //criar elementos de interface
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
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
//
////    fGeoMesh->ResetReference();
////    fTransportMesh->LoadReferences();
////    int64_t nel = fTransportMesh->ElementVec().NElements();
////    // Creation of interface elements
////    for(int el = 0; el < nel; el++)
////    {
////        TPZCompEl * cel = fTransportMesh->ElementVec()[el];
////        if(!cel) continue;
////        TPZGeoEl * gel = cel->Reference();
////        if(!gel) {continue;}
////        if(gel->HasSubElement()) {continue;}
////        int index = cel ->Index();
////        if(cel->Dimension() == fTransportMesh->Dimension())
////        {
////            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fTransportMesh->ElementVec()[index]);
////            if(!InterpEl) {
////                continue;
////            }
////            InterpEl->CreateInterfaces();
////        }
////    }
//    
    mphysics->CleanUpUnconnectedNodes();
    mphysics->AdjustBoundaryElements();
    mphysics->AutoBuild();
    
    return mphysics;
}

void CreateInterfaceMultiphysicsElement(TPZCompMesh* mphysics)
{
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    
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
    
    mphysics->CleanUpUnconnectedNodes();
    mphysics->AdjustBoundaryElements();
    //mphysics->AutoBuild();
}

void CreateInterfaceElement(TPZCompMesh* cmesh)
{
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();

    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
    
    cmesh->CleanUpUnconnectedNodes();
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
}

void SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
//	ofstream file("Solutout");
//	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void ResolverComReconstGradiente(REAL deltaX,REAL maxTime,TPZManVector<TPZCompMesh *,3> meshvec,TPZCompMesh* mphysics,TPZGradientReconstruction *gradreconst, bool useRK2)
{
    
    //inserir solucao inicial
    TPZAnalysis an(mphysics);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    std::string outputfile;
	outputfile = "TransientSolutionComRecGrad";
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    
    TPZMaterial *mat = mphysics->FindMaterial(matId);
    TPZTracerFlow * material = dynamic_cast<TPZTracerFlow *>(mat);
    
    //Reconstrucao do gradiente e linearizacao da solucao
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZFMatrix<REAL> datagradients;
    meshvec[0]->Reference()->ResetReference();
    meshvec[0]->LoadReferences();
    gradreconst-> ProjectionL2GradientReconstructed(meshvec[0], matIdL2Proj);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    Initialsolution = mphysics->Solution();
    
//------------- Filtrar Equacaoes da pressao-flux -----------
    FilterPressureFluxEquation(material, meshvec, mphysics, an);
    an.Solve();
    TPZAutoPointer< TPZMatrix<STATE> > matKPressureFlux;
	TPZFMatrix<STATE> fvecPressureFlux;
    matKPressureFlux = an.Solver().Matrix();
    fvecPressureFlux = an.Rhs();
    
//#ifdef LOG4CXX
//    if(logdata->isDebugEnabled())
//    {
//        std::stringstream sout;
//        matKPressureFlux->Print("matKPressureFlux = ", sout,EMathematicaInput);
//        fvecPressureFlux.Print("fvecPressureFlux = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str())
//    }
//#endif

    //posprocessar solucao
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    PosProcessMultphysics(meshvec,mphysics,an,plotfile);
    TPZFMatrix<STATE> SolutionQ = meshvec[1]->Solution();
    TPZFMatrix<STATE> SolutionP = meshvec[2]->Solution();

//--------- Calculando DeltaT maximo para a cond. CFL ------
    REAL deltaT=0.;
    CondCFL(SolutionQ, deltaX, maxTime, deltaT);
    int NDt = maxTime/deltaT;
    int nfile = 0;
    if(NDt<=100) nfile = 1;
    else nfile = int (NDt/100);
    int rest = nfile%2;
    if(rest!=0) nfile++;
    material->SetTimeStep(deltaT);
    
//------------------- Criando matriz de massa (matM) ---------------------
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(material, mphysics);
    TPZAutoPointer <TPZMatrix<STATE> > matMAux;
     if(useRK2==true){
         material->SetTrueRungeKuttaTwo();
         matMAux = MassMatrix(material, mphysics);
     }
//#ifdef LOG4CXX
//    if(logdata->isDebugEnabled())
//    {
//        std::stringstream sout;
//        matM->Print("matM = ", sout,EMathematicaInput);
//        fvecM.Print("fvecM = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str())
//    }
//#endif

    //---- Filtrar Equacaoes da saturacao ----
    //Criando matriz de rigidez (matK) e vetor de carga
    TPZAutoPointer< TPZMatrix<STATE> > matK;
	TPZFMatrix<STATE> fvecK;
    FilterSaturationEquation(material, meshvec, mphysics, an, true);
    matK = an.Solver().Matrix();
    fvecK = an.Rhs();
    
//    #ifdef LOG4CXX
//    if(logdata->isDebugEnabled())
//    {
//        std::stringstream sout;
//        matK->Print("matK = ", sout,EMathematicaInput);
//        fvecK.Print("fvecK = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str())
//    }
//    #endif
    
    
//-------- Resolver o sistema numerico no tempo -------------------
	int64_t nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp1(nrows,1,1.0);
    TPZFMatrix<STATE> TotalRhstemp2(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
	
	REAL TimeValue = 0.0;
	int cent = 1;
    TimeValue = cent*deltaT;
    
	while (TimeValue <= maxTime)
	{
        ftimeatual  = TimeValue;
		// This time solution i for Transient Analytic Solution
		material->SetTimeValue(TimeValue);

//--------- Resolver usando Runge-Kutta -------------------
        //primeiro estagio de Runge-Kutta
        matM->Multiply(Lastsolution,TotalRhstemp1);
		TotalRhs = TotalRhstemp1 + fvecK;
		an.Rhs() = TotalRhs;
		an.Solve();
        
        if(useRK2==true){
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            gradreconst-> ProjectionL2GradientReconstructed(meshvec[0], matIdL2Proj);
            TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
            an.LoadSolution(mphysics->Solution());

            //segundo estagio de Runge-Kutta
            matMAux->Multiply(Lastsolution,TotalRhstemp1);
            Lastsolution = an.Solution();
            matM->Multiply(Lastsolution,TotalRhstemp2);
            TotalRhs = TotalRhstemp1 + (TotalRhstemp2 + fvecK);
            TotalRhs = STATE(0.5)*TotalRhs;
            an.Rhs() = TotalRhs;
            an.Solution().Zero();
            an.Solve();
        }
//---------------------------------------------------------
        
        //Reconstrucao do gradiente e linearizacao da solucao
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
        gradreconst-> ProjectionL2GradientReconstructed(meshvec[0], matIdL2Proj);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        an.LoadSolution(mphysics->Solution());
        Lastsolution = an.Solution();

        
        if(cent%nfile==0){
            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            meshvec[1]->LoadSolution(SolutionQ);
            meshvec[2]->LoadSolution(SolutionP);
            PosProcessMultphysics(meshvec,mphysics,an,plotfile);
        }
		
        cent++;
		TimeValue = cent*deltaT;
        an.Solution().Zero();
    }
}

void ResolverSemReconstGradiente(REAL deltaX,REAL maxTime,TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics)
{
    
    //inserir solucao inicial
    TPZAnalysis an(mphysics);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    
    std::string outputfile;
	outputfile = "TransientSolutionSemRecGrad";
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    
    TPZMaterial *mat = mphysics->FindMaterial(matId);
    TPZTracerFlow * material = dynamic_cast<TPZTracerFlow *>(mat);
    
    //------------- Filtrar Equacaoes da pressao-flux -----------
    FilterPressureFluxEquation(material, meshvec, mphysics, an);
    an.Solve();
    TPZAutoPointer< TPZMatrix<STATE> > matKPressureFlux;
	TPZFMatrix<STATE> fvecPressureFlux;
    matKPressureFlux = an.Solver().Matrix();
    fvecPressureFlux = an.Rhs();
    
    //#ifdef LOG4CXX
    //    if(logdata->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        matKPressureFlux->Print("matKPressureFlux = ", sout,EMathematicaInput);
    //        fvecPressureFlux.Print("fvecPressureFlux = ", sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //    }
    //#endif
    
    //posprocessar solucao
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    PosProcessMultphysics(meshvec,mphysics,an,plotfile);
    TPZFMatrix<STATE> SolutionQ = meshvec[1]->Solution();
    TPZFMatrix<STATE> SolutionP = meshvec[2]->Solution();
    
    //--------- Calculando DeltaT maximo para a cond. CFL ------
    REAL deltaT=0.;
    CondCFL(SolutionQ, deltaX, maxTime, deltaT);
    int NDt = maxTime/deltaT;
    int nfile = 0;
    if(NDt<=100) nfile = 1;
    else nfile = int (NDt/100);
    int rest = nfile%2;
    if(rest!=0) nfile++;
    material->SetTimeStep(deltaT);
    //------------------- Criando matriz de massa (matM) ---------------------
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(material, mphysics);
    //#ifdef LOG4CXX
    //    if(logdata->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        matM->Print("matM = ", sout,EMathematicaInput);
    //        fvecM.Print("fvecM = ", sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //    }
    //#endif
    
    //---- Filtrar Equacaoes da saturacao ----
    //Criando matriz de rigidez (matK) e vetor de carga
    TPZAutoPointer< TPZMatrix<STATE> > matK;
	TPZFMatrix<STATE> fvecK;
    FilterSaturationEquation(material, meshvec, mphysics, an, true);
    matK = an.Solver().Matrix();
    fvecK = an.Rhs();
    
    //    #ifdef LOG4CXX
    //    if(logdata->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        matK->Print("matK = ", sout,EMathematicaInput);
    //        fvecK.Print("fvecK = ", sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //    }
    //    #endif
    
    
    //-------- Resolver o sistema numerico no tempo -------------------
	int64_t nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
	
	REAL TimeValue = 0.0;
	int cent = 1;
    
	TimeValue = cent*deltaT;
	while (TimeValue <= maxTime)
	{
        ftimeatual = TimeValue;
		// This time solution i for Transient Analytic Solution
		material->SetTimeValue(TimeValue);
		matM->Multiply(Lastsolution,TotalRhstemp);
        
        //        #ifdef LOG4CXX
        //                if(logdata->isDebugEnabled())
        //            {
        //                std::stringstream sout;
        //                sout<< " tempo = " << cent;
        //                Lastsolution.Print("\nIntial conditions = ", sout,EMathematicaInput);
        //                TotalRhstemp.Print("Mat Mass x Last solution = ", sout,EMathematicaInput);
        //                LOGPZ_DEBUG(logdata,sout.str())
        //            }
        //        #endif
        
		TotalRhs = fvecK + TotalRhstemp;
		an.Rhs() = TotalRhs;
		an.Solve();
		Lastsolution = an.Solution();
        
        if(cent%nfile==0){
            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            meshvec[1]->LoadSolution(SolutionQ);
            meshvec[2]->LoadSolution(SolutionP);
            PosProcessMultphysics(meshvec,mphysics,an,plotfile);
        }
		
        cent++;
		TimeValue = cent*deltaT;
    }
}


void CondCFL(TPZFMatrix<STATE> SolutionQ, REAL deltaX, REAL maxTime, REAL &deltaT)
{
    int nr = SolutionQ.Rows();
    int nc = SolutionQ.Cols();
    REAL maxsolq = 0.;
    REAL temp = 0.;
    for(int ir = 0; ir<nr; ir++)
    {
        for(int jc = 0; jc<nc; jc++)
        {
            temp = SolutionQ.GetVal(ir, jc);
            if(fabs(temp) > maxsolq) maxsolq = fabs(temp);
        }
    }
    maxsolq /=fporos;
    int NDt = 10;
    deltaT = 0.;
    REAL CFL = 0.;
    int count = 1;
    if(fishomogmedium==true)
    {
        deltaT = maxTime/NDt;
        CFL = maxsolq*(deltaT/deltaX);
        while (CFL > 0.2)
        {
            NDt = 2*NDt;
            deltaT = maxTime/NDt;
            CFL = maxsolq*(deltaT/deltaX);
            count++;
        }
    }
    else
    {
        deltaT = deltaX/(2.*maxsolq);
        CFL = 2.*maxsolq*(deltaT/deltaX);
        while (CFL > 0.1)
        {
            NDt = 2*NDt;
            deltaT = maxTime/NDt;
            CFL = 2.*maxsolq*(deltaT/deltaX);
            count++;
        }
    }
    std::cout<<"\nNumero de passos de tempo NDt = " << NDt << ", valor DeltaT = " << deltaT << " e CFL = " << CFL<< "\n ";
}

#include "pzfstrmatrix.h"
void FilterPressureFluxEquation(TPZTracerFlow *mymaterial, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics, TPZAnalysis &an)
{
 
    int ncon_saturation = meshvec[0]->NConnects();
    int ncon = mphysics->NConnects();
    TPZManVector<int64_t> active(0);
    for(int i = ncon_saturation; i<ncon; i++)
    {
        TPZConnect &con = mphysics->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = mphysics->Block().Position(seqnum);
        int blocksize = mphysics->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+blocksize);
        for(int ieq = 0; ieq<blocksize; ieq++){
            active[vs+ieq] = pos+ieq;
        }
    }
    
    mymaterial->SetCurrentState();
    mymaterial->SetPressureEqFilter();
    TPZSkylineStructMatrix matsk(mphysics);
    //matsk.SetNumThreads(4);
    //TPZFStructMatrix matsk(mphysics);
    //TPZBandStructMatrix matsk(mphysics);
    matsk.EquationFilter().SetActiveEquations(active);
	an.SetStructuralMatrix(matsk);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
    //step.SetDirect(ELU);
	an.SetSolver(step);
    an.Assemble();
}

void FilterSaturationEquation(TPZTracerFlow *mymaterial, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics, TPZAnalysis &an, bool currentstate)
{
    int ncon_saturation = meshvec[0]->NConnects();
    TPZManVector<int64_t> active(0);
    for(int i = 0; i<ncon_saturation; i++)
    {
        TPZConnect &con = mphysics->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        if(seqnum==-1) continue;
        int pos = mphysics->Block().Position(seqnum);
        int blocksize = mphysics->Block().Size(seqnum);
        
//        //ativar equacao
//        int vs = active.size();
//        active.Resize(vs+blocksize);
//        for(int ieq = 0; ieq<blocksize; ieq++)
//        {
//            active[vs+ieq] = pos+ieq;
//        }
        
        //Ativar apenas a ultima equacao, que corresponde a funcao de base constante
        int vs = active.size();
        active.Resize(vs+1);
        int ieq = blocksize-1;
        active[vs] = pos+ieq;
    }
    if(currentstate==true)
    {
        mymaterial->SetCurrentState();
        mymaterial->SetFalsePressureEqFilter();
        TPZSkylineStructMatrix matsk(mphysics);
        //matsk.SetNumThreads(6);
//        TPZParFrontStructMatrix<TPZFrontSym<STATE> > matsk(mphysics);
//        matsk.SetNumThreads(6);
        
        matsk.EquationFilter().SetActiveEquations(active);
        an.SetStructuralMatrix(matsk);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Assemble();
    }
    else
    {
        mymaterial->SetLastState();
        mymaterial->SetFalsePressureEqFilter();
        TPZSkylineNSymStructMatrix matst(mphysics);
        //matst.SetNumThreads(6);
        
//        TPZParFrontStructMatrix<TPZFrontNonSym<STATE> > matst(mphysics);
//        matst.SetNumThreads(6);
        
        matst.EquationFilter().SetActiveEquations(active);
        an.SetStructuralMatrix(matst);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Assemble();
    }
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    //TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(3), vecnames(2);
	
    scalnames[0] = "Saturation";
    scalnames[1] = "Pressure";
    scalnames[2] = "DivFlux";
    vecnames[0]  = "Flux";
    vecnames[1] = "SaturationFlux";
    
    if(fishomogmedium==true)
    {
        scalnames.Resize(4);
        scalnames[3] = "ExactSaturation";
    }
    
	const int dim = 2;
	int div =0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
//	std::ofstream out("malha.txt");
//	an.Print("nothing",out);
    
}

void PosProcessFlux(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(0), vecnames(1);
	vecnames[0]= "FluxL2";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malhaflux.txt");
	an.Print("nothing",out);
}

//Riemann problem
void ForcingInicial(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    //double y = pt[1];
    if(x<0)disp[0] = 1.;
    if(x>0)disp[0] = 0.;
}


void Permeability(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
	double x = pt[0];
    double y = pt[1];
    
    REAL k1D = fk1/fkref;
    REAL k2D = fk2/fkref;
    
    REAL LyD = fLy/fLref;
    REAL LxD = fLx/fLref;
    
    if(fhetertest==1)
    {
        REAL dx = LxD/(8.);
        REAL dy = LyD/(8.);
        
        //disp[0] = k1;
        if((x>= dx && x<=2.*dx) && (y>= 4.*dy && y<=LyD))
        {
            disp[0] = k2D;
        }
        else if ((x>= 6.*dx && x<=7.*dx) && (y>= 0. && y<= 4.*dy))
        {
            disp[0] = k2D;
        }
            else disp[0] = k1D;
    }
    
    if(fhetertest==2)
    {
        REAL dx = LxD/(64.);
        REAL dy = LyD/(16.);
        
        //disp[0] = k1;
        if((/*(x>=0. && x<=2.*dx) ||*/ (x>=14.*dx && x<=18.*dx) || (x>=30.*dx && x<=34.*dx) || (x>=46.*dx && x<=50.*dx) || (x>=62.*dx && x<=LxD)) && (y>= 0. && y<= 12.*dy))
        {
            disp[0] = k2D;
        }
        else if (((x>=6.*dx && x<=10.*dx) || (x>=22.*dx && x<=26.*dx) || (x>=38.*dx && x<=42.*dx) || (x>=54.*dx && x<=58.*dx)) && (y>= 4.*dy && y<=LyD))
        {
            disp[0] = k2D;
        }
        else disp[0] = k1D;
    }
    
    if(fhetertest==3)
    {
        REAL dx = LxD/(64.);
        REAL dy = LyD/(16.);
        
        if(((x>=12.*dx && x<=20.*dx)/*|| (x>=28.*dx && x<=36.*dx)*/ || (x>=44.*dx && x<=52.*dx)) && (y>= 4.*dy && y<=12.*dy))
        {
            disp[0] = k2D;
        }
        else disp[0] = k1D;
    }
    
    if(fhetertest==4)
    {
        REAL f1 = 0., f3 = 0.,f2 = 0.,f4 = 0.,f5 = 0;
        REAL temp1=0., temp2=0., temp3=0., temp4 =0., temp5;
        
        temp1 = -4.11184*x*x + (6.57895 - 65.7895*y)*y +x*(1.64474 + 29.6053*y);
        f1 = 0.0850595*exp(temp1);
        
        temp2 = -0.0868056*x*x + (50.- 138.889*y)*y + x*(-0.9375 + 5.55556*y);
        f2 = 0.00714023*exp(temp2);
        
        temp3 = -80.3472 - 0.0868056*x*x + (211.111- 138.889*y)*y + x*(-4.09722 + 5.55556*y);
        f3 = 0.663146*exp(temp3);
        
        temp4 = -62.5 - 4.11184*x*x + x*(31.25 - 29.6053*y) + (125.0 - 65.7895*y)*y;
        f4 = 2.28204*exp(temp4);
        
        temp5 = -2.864788975654116*x*x;
        f5 = 1.909859317102744*exp(temp5);
        
        REAL eps = 0.00001;
        disp[0] =  (f2 + f3 + f5/4.) + (f1 + f4)/4. + eps;
        //disp[0] = 1 - (f1 + f4)/4. + eps;
        //disp[0] = (2.*f5 + (f1 + f4))/4.5641 + eps;
    }
    
    if(fhetertest==5)
    {
        REAL f1 = 0., f2 = 0.;
        REAL temp1=0., temp2=0.;
        
        temp1 = -1.04167*x*x + (8.33333 - 16.6667*y)*y + x*(2.08333 + 4.16667*y);
        f1 =0.0178078*exp(temp1);
        
        temp2 = -12.5 - 1.04167*x*x + x*(6.25 - 4.16667*y) + (25.0 - 16.6667*y)*y;
        f2 = 1.1486*exp(temp2);
        
        disp[0] = 1.149 - (f1 + f2)/2.;
    }
    
    if(fhetertest==6)
    {
        REAL f1 = 0.;
        REAL temp1=0.;
        
        temp1 = -0.789141*x*x + x*(2.84091 + 0.631313*y) + (11.3636 - 12.6263*y)*y;
        f1 =0.00340644*exp(temp1);
        
        disp[0] = 1. - f1 ;
    }


}

TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZTracerFlow * mymaterial, TPZCompMesh* mphysics){
    
    mymaterial->SetLastState();
    mymaterial->SetFalsePressureEqFilter();
    
//    TPZParFrontStructMatrix<TPZFrontNonSym<STATE> > matsp(mphysics);
//    matsp.SetNumThreads(6);

    
    TPZSkylineNSymStructMatrix matsp(mphysics);
   // matsp.SetNumThreads(6);
    
    //TPZSpStructMatrix matsp(mphysics);
	
	std::set< int > materialid;
	int matid = mymaterial->MatId();
	materialid.insert(matid);
	matsp.SetMaterialIds (materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<STATE> Un;
    //TPZMatrix<REAL> *matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}

void StiffMatrixLoadVec(TPZTracerFlow *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZAutoPointer< TPZMatrix<STATE> > &matK1, TPZFMatrix<STATE> &fvec){
    
    mymaterial->SetCurrentState();
    
    //TPZFStructMatrix matsk(mphysics);
    TPZSkylineStructMatrix matsk(mphysics);
	an.SetStructuralMatrix(matsk);
    //matsk.SetNumThreads(30);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
	//step.SetDirect(ELU);
	an.SetSolver(step);
	
    an.Assemble();
	
	fvec = an.Rhs();
    matK1 = an.Solver().Matrix();

}

TPZCompMesh *SetCondicaoInicial(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(matId, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(ForcingInicial, 5);
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
            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }

    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
    cmesh->LoadReferences();
    TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
    cmesh->ExpandSolution();

	return cmesh;
}

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    double x = pt[0];
    //double y = pt[1];
    
    REAL tp = ftimeatual;
    u[0]=0.;
    
    REAL AreaPoro = fLy*1.*fporos;
    REAL qinDia = fvazaoentrada/AreaPoro;//m/d
    REAL qin = qinDia/86400.;//m/s
    REAL num = fvisref*fLref;
    REAL denom = fkref*fpref;
    REAL qinD = (-1.)*(qin*num/denom);
    REAL velx = -1.*qinD/fporos;
    
    REAL ptx = x-velx*tp;
    
    if(ptx < 0.0) u[0] = 1.;
    if(ptx > 0.0) u[0] = 0.;
}


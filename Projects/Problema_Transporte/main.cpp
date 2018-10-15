
#include <time.h>
#include <stdio.h>
#include <fstream>
#include <cmath>

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
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzgengrid.h"

#include "pzpoisson3d.h"

#include "pzl2projection.h"

#include "pzconvectionproblem.h"
#include "pzgradientreconstruction.h"

#include "pzbuildmultiphysicsmesh.h"

#include "TPZCompElDisc.h"
#include "pzbndcond.h"

#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

using namespace std;

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

#ifdef USING_BOOST
#include <boost/math/special_functions/erf.hpp> //Required for erfc function on windows
#endif

const int matId = 1;
const int inflow = 0;
const int outflow = 3;

int matIdL2Proj = 2;

int neumann = 1;

REAL ftimeatual = 0.;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;

TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly, bool triang_elements);
TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh,int pOrder,TPZMatConvectionProblem * &material);
TPZGeoMesh *GMesh3(REAL Lx, REAL Ly);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void CreatInterface(TPZCompMesh *cmesh);

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix);
void SolveSistTransient(REAL deltaT, TPZMatConvectionProblem * &mymaterial, TPZCompMesh* cmesh, TPZGradientReconstruction *gradreconst);
void SolveSistTransient(REAL deltaT, TPZMatConvectionProblem * &mymaterial, TPZCompMesh* cmesh);

void PosProcessSolutionTrans(TPZCompMesh* cmesh, TPZAnalysis &an, std::string plotfile);
TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZMatConvectionProblem * mymaterial, TPZCompMesh* cmesh);
void StiffMatrixLoadVec(TPZMatConvectionProblem *mymaterial, TPZCompMesh*cmesh, TPZAnalysis &an, TPZAutoPointer< TPZMatrix<STATE> > &matK1, TPZFMatrix<STATE> &fvec);

//Ativar apenas a ultima equacao, que corresponde a funcao de base constante
void FilterEquation(TPZMatConvectionProblem *mymaterial, TPZCompMesh *cmesh, TPZAnalysis &an, bool currentstate, TPZManVector<int64_t> &nonactive);
void CleanGradientSolution(TPZFMatrix<STATE> &Solution, TPZManVector<int64_t> &Gradients);

void ForcingInicial(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
TPZCompMesh *SetCondicaoInicial(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);

//problema estacionario
//variaveis iniciais:
REAL valC = 10.;
TPZGeoMesh *GMesh2(REAL Lx, REAL Ly,bool triang_elements);
TPZCompMesh *MalhaComp2(TPZGeoMesh * gmesh,int pOrder/*,TPZMatConvectionProblem * &material*/);
void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp/*,TPZFMatrix<STATE> &deriv*/);
void SolucaoExata(const TPZVec<REAL> &pt, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv);
void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
void FilterEquation(TPZCompMesh *cmesh, TPZAnalysis &an);
void DirichletCond(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void ConvGradU(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &deriv);
void ForcingInicial2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void PosProcessSolution(TPZCompMesh* cmesh, TPZAnalysis &an, std::string plotfile);

//Riemann Problem
bool triang = false;
bool userecgrad = false;
bool calcresiduo =true;
REAL teta =0.;// M_PI/6.;

int main(int argc, char *argv[])
{
//#ifdef LOG4CXX
//    InitializePZLOG();
//#endif
    
    TPZVec<REAL> erros;
    ofstream arg0("Erro.txt");
    int p = 1;
    int h = 4;
    REAL Lx=1., Ly=1.;
    

    TPZGeoMesh *gmesh = MalhaGeom(Lx,Ly,triang);
    //TPZGeoMesh *gmesh = GMesh3(Lx,Ly);
    UniformRefine(gmesh,h);
    
    TPZMatConvectionProblem *material;
    TPZCompMesh *cmesh = MalhaComp(gmesh, p, material);
    
    
    //Refinar elemnto: malha adaptada
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,3, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    //Refinar elemnto: malha adaptada
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,12, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,2, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,13, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,7, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,8, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,6, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,9, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh,2, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();

    
    //Set initial conditions
    TPZAnalysis an(cmesh,0);
    int nrs = an.Solution().Rows();
    TPZVec<STATE> solini(nrs,0.);
    TPZCompMesh  * cmeshL2 = SetCondicaoInicial(gmesh, p, solini);
    
    TPZAnalysis anL2(cmeshL2,0);
    ResolverSistema(anL2, cmeshL2, true);
    an.LoadSolution(anL2.Solution());
    //an.Solution().Print("sol_S0");
    
    gmesh->ResetReference();
    cmesh->LoadReferences();
    CreatInterface(cmesh);
    
    ofstream arg1("gmesh.txt");
    gmesh->Print(arg1);
    
    ofstream arg2("cmesh.txt");
    cmesh->Print(arg2);
    
    //--------- Calculando DeltaT maximo para a cond. CFL ------
    REAL maxTime = 2.5;
    REAL deltaX = Lx/pow(2.,h);
    REAL solV = 1.;//velocidade maxima
    int NDt = 10;
    REAL deltaT = 0.;
    REAL CFL = 0.;
    int count = 1;
    
    deltaT = maxTime/NDt;
    CFL = solV*(deltaT/deltaX);
    while (CFL > 0.1)
    {
        NDt = 2*NDt;
        deltaT = maxTime/NDt;
        CFL = solV*(deltaT/deltaX);
        count++;
    }
    
    std::cout<<"\nNdt = " << NDt <<" e "<< " CFL = " << CFL <<std::endl;
    
    material->SetTimeStep(deltaT);
    material->SetTimeValue(maxTime);

    if(userecgrad==true)
    {
        TPZGradientReconstruction *gradreconst = new TPZGradientReconstruction(false,1.);
        
        if(triang==true){
            TPZVec<REAL> LxLyLz(2,0.);
            LxLyLz[0] = 1.; LxLyLz[1]=1.;
            
            TPZVec<int> MatIdBC(2,0);
            MatIdBC[0] = -1; MatIdBC[1] = -3;
            
            TPZVec<REAL>  Xmin(2,0.);
            TPZVec<REAL> Xmax(2,0.);
            TPZManVector<TPZVec<REAL> > coordmin(2,0.);
            TPZManVector<TPZVec<REAL> > coordmax(2,0.);
            Xmin[0] = 0.; Xmin[1] = 0.; Xmax[0] = 1.; Xmax[1] = 0.;
            coordmin[0] = Xmin; coordmax[0] = Xmax;
            Xmin[0] = 0.; Xmin[1] = 1.; Xmax[0] = 1.; Xmax[1] = 1.;
            coordmin[1] = Xmin; coordmax[1] = Xmax;
            
            gradreconst->SetDataGhostsNeighbors(LxLyLz, MatIdBC, coordmin, coordmax);
        }
        SolveSistTransient(deltaT, material, cmesh, gradreconst);
    }else
    {
        SolveSistTransient(deltaT, material, cmesh);
    }
    
    return EXIT_SUCCESS;
    
}

//problema estacionario
int mainestacionario(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    TPZVec<REAL> erros;
    ofstream arg0("Erro.txt");
    int p, h;
    
     arg0<<"\n Erro problema estacionario. Malha triangular = " << triang <<endl;
    for(p=1; p<2; p++)
    {
        arg0<<"\n Ordem p = " << p <<endl;
        for(h =1; h < 8; h++)
        {
            arg0<<"\n Refinamento h  = " << h <<endl;
            
            //------------- PROBLEMA ESTACIONARIO -----------------
            TPZGeoMesh *gmesh = GMesh2(2.,2.,triang);
            UniformRefine(gmesh,h);
            //ofstream arg1("gmesh.txt");
            //gmesh->Print(arg1);
            
            TPZCompMesh *cmesh = MalhaComp2(gmesh, p);
            CreatInterface(cmesh);
            //ofstream arg2("cmesh.txt");
            //cmesh->Print(arg2);
            
            TPZAnalysis an(cmesh);
            FilterEquation(cmesh, an);
            
            if(userecgrad==false)
            {
                 an.Solve();
                 
                 //pos-process
                 std::string outputfile;
                 outputfile = "SolutionSemRecGrad";
                 std::stringstream outputfiletemp;
                 outputfiletemp << outputfile << ".vtk";
                 std::string plotfile = outputfiletemp.str();
                 PosProcessSolution(cmesh,an,plotfile);
                 
                 arg0<<" \n----ERRO SEM RECONSTRUIR GRADIENTE----" <<endl;
                 an.SetExact(*SolucaoExata);
                bool store_errors = false;
                 an.PostProcessError(erros, store_errors, arg0);
            }
            
            //Reconstruindo Gradient
            else
            {
                TPZGradientReconstruction *gradreconst = new TPZGradientReconstruction(false,1.);
                TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(ForcingInicial2, 5);
                gradreconst->fGradData->SetForcingFunctionExact(forcef);
                gradreconst->fGradData->EnableForcinFucnction();
                gradreconst-> ProjectionL2GradientReconstructed(cmesh, matIdL2Proj);
                an.Solution().Zero();
                an.LoadSolution(cmesh->Solution());
                
                //pos-process
                std::string outputfile2;
                outputfile2 = "SolutionComRecGrad";
                std::stringstream outputfiletemp2;
                outputfiletemp2 << outputfile2 << ".vtk";
                std::string plotfile2 = outputfiletemp2.str();
                PosProcessSolution(cmesh,an,plotfile2);
                
                arg0<<" \n----ERRO COM RECONSTRUCAO DO GRADIENTE----" <<endl;
                an.SetExact(*SolucaoExata);
                bool store_errors = false;
                an.PostProcessError(erros, store_errors, arg0);
            }
            
            cmesh->CleanUp();
            delete cmesh;
            delete gmesh;
        }
    }
    return EXIT_SUCCESS;
    
}


TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly, bool triang_elements)
{
	
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
	REAL valx, valx2;
    REAL valy, valy2;
    
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
        valy = 0.;
        valx2 = valx*cos(teta) - valy*sin(teta);
        valy2 = valx*sin(teta) + valy*cos(teta);
        
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0,valx2);//coord X
		Node[id].SetCoord(1,valy2);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
        valy = Ly;
        valx2 = valx*cos(teta) - valy*sin(teta);
        valy2 = valx*sin(teta) + valy*cos(teta);
        
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0,valx2);//coord X
		Node[id].SetCoord(1,valy2);//coord Y
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
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
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
    
//    TPZManVector<int,2> nx(2,1);
//    TPZManVector<REAL,3> x0(3,0),x1(3,0);
//    x0[0]=0.; x0[1]=0.;
//    x1[0]=Lx; x1[1]=Ly;
//    
//    TPZGenGrid gengrid(nx,x0,x1);
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    if(triang_elements)
//    {
//        gengrid.SetElementType(1);
//    }
//    gengrid.Read(gmesh);
//    TPZManVector<REAL,3> firstpoint(3,0.),secondpoint(3,0.);
//    
//    secondpoint[0] = 1.0;
//    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc0);
//    
//    firstpoint = secondpoint;
//    secondpoint[0] = 1.0;
//    secondpoint[1] = 1.0;
//    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc1);
//    
//    firstpoint = secondpoint;
//    secondpoint[0] = 0.0;
//    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc2);
//    
//    firstpoint = secondpoint;
//    secondpoint[1] = 0.0;
//    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc3);
    
	return gmesh;
}

TPZGeoMesh *GMesh3(REAL Lx, REAL Ly)
{
	
	int Qnodes = 5;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
	
    REAL teta = M_PI/4.;
    
	//indice dos nos
	int64_t id = 0;
	REAL valx, valx2;
    REAL valy, valy2;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
        valy = 0.;
        valx2 = valx*cos(teta) - valy*sin(teta);
        valy2 = valy*cos(teta) + valx*sin(teta);
        
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx2);//coord X
		Node[id].SetCoord(1 ,valy2);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
        valy = Ly;
        valx2 = valx*cos(teta) - valy*sin(teta);
        valy2 = valy*cos(teta) + valx*sin(teta);
        
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0,valx2);//coord X
		Node[id].SetCoord(1,valy2);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,Lx/2.);//coord X
    Node[id].SetCoord(1 ,Ly/2.);//coord Y
    gmesh->NodeVec()[id] = Node[id];
	
//----- indice dos elementos ----
	id = 0;
    TopolTriang[0] = 0;
    TopolTriang[1] = 1;
    TopolTriang[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
    id++;
    
    TopolTriang[0] = 1;
    TopolTriang[1] = 2;
    TopolTriang[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
    id++;
    
    TopolTriang[0] = 2;
    TopolTriang[1] = 3;
    TopolTriang[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
    id++;
    
    TopolTriang[0] = 3;
    TopolTriang[1] = 0;
    TopolTriang[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 0;
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


TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh, int pOrder,TPZMatConvectionProblem * &material)
{
	/// criar materiais
	int dim = 2;
	
	material = new TPZMatConvectionProblem(matId,dim);
	TPZMaterial * mat(material);
	
	TPZVec<REAL> convdir(dim,0.);
    TPZVec<REAL> convdirtemp(dim,0.);
    convdirtemp[0] = 1.;
    convdir[0] = convdirtemp[0]*cos(teta) - convdirtemp[1]*sin(teta);
    convdir[1] = convdirtemp[1]*cos(teta) + convdirtemp[0]*sin(teta);
    
	REAL flux = 0.;
    REAL rho = 1.;
	
	material->SetParameters(rho,convdir);
	material->SetInternalFlux( flux);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->InsertMaterialObject(mat);
    
	///Inserir condicao de contorno
    TPZAutoPointer<TPZFunction<STATE> > solExata;
    solExata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solExata);
    
	TPZFMatrix<STATE> val1(1,1,0.), val2(2,1,0.);
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
//	TPZMaterial * BCond1 = material->CreateBC(mat, bc1, neumann, val1, val2);
	REAL uD =2.;
    val1(0,0) = uD;
	val2(0, 0) = 1.;
	val2(1, 0) = 0.;
	TPZMaterial * BCond3 = material->CreateBC(mat, bc3, inflow, val1, val2);
	val1.Zero();
	val2(0,0) = .5;
	val2(1, 0) = 0.;
	TPZMaterial * BCond1 = material->CreateBC(mat, bc1, outflow, val1, val2);

    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond2);
    
    
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,material->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
	
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

    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;
}

#include "pzl2projection.h"
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
            if(triang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }

    
    //TPZCompElDisc::SetTotalOrderShape(cmesh);
    
	return cmesh;
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
	// Re-constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}

void CreatInterface(TPZCompMesh *cmesh){
    
    for(int64_t el = 0; el < cmesh->ElementVec().NElements(); el++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        int64_t index = compEl ->Index();
        if(compEl->Dimension() == cmesh->Dimension() || compEl->Dimension() == cmesh->Dimension()-1)
        {
            TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces(false);
        }
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
}

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix)
{
    if(symmetric_matrix ==true){
        TPZSkylineStructMatrix skmat(fCmesh);
        an.SetStructuralMatrix(skmat);
        TPZStepSolver<STATE> direct;
        direct.SetDirect(ELDLt);
        an.SetSolver(direct);
    }
    else{
        TPZBandStructMatrix bdmat(fCmesh);
        an.SetStructuralMatrix(bdmat);
        TPZStepSolver<STATE> direct;
        direct.SetDirect(ELU);
        an.SetSolver(direct);
    }
	an.Run();
	
	//Saida de Dados: solucao txt
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);
}

void SolveSistTransient(REAL deltaT, TPZMatConvectionProblem * &mymaterial, TPZCompMesh* cmesh, TPZGradientReconstruction *gradreconst){
	
    TPZAnalysis an(cmesh);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    std::string outputfile;
	outputfile = "SolutionComRecGrad";
    
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    
    //gradient reconstruction
    TPZFMatrix<REAL> datagradients;
    gradreconst-> ProjectionL2GradientReconstructed(cmesh, matIdL2Proj);
    Initialsolution = cmesh->Solution();
    
    PosProcessSolutionTrans(cmesh,an,plotfile);

    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(mymaterial, cmesh);
    mymaterial->SetTrueRungeKuttaTwo();
    TPZAutoPointer <TPZMatrix<STATE> > matMAux = MassMatrix(mymaterial, cmesh);
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
            std::stringstream sout;
        	matM->Print("matM = ", sout,EMathematicaInput);
        	LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
    
    //Criando matriz de rigidez (matK) e vetor de carga
	TPZAutoPointer< TPZMatrix<STATE> > matK;
	TPZFMatrix<STATE> fvec;
    //StiffMatrixLoadVec(mymaterial, cmesh, an, matK, fvec);
    TPZManVector<int64_t> nonactive(0);
    FilterEquation(mymaterial, cmesh, an, true,nonactive);
    matK = an.Solver().Matrix();
    fvec = an.Rhs();
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{

        std::stringstream sout;
        matK->Print("matK = ", sout,EMathematicaInput);
        fvec.Print("fvec = ", sout,EMathematicaInput);
        //Print the temporal solution
        Initialsolution.Print("Initial conditions = ", sout,EMathematicaInput);
        TPZFMatrix<STATE> Temp;
        matM->Multiply(Initialsolution,Temp);
        Temp.Print("Temp matM = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
    
    
	int64_t nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp1(nrows,1,0.0);
    TPZFMatrix<STATE> TotalRhstemp2(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
    
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT;
    REAL maxTime;
    mymaterial->GetTimeValue(maxTime);
	while (TimeValue <= maxTime)
	{
        ftimeatual  = TimeValue;
        
		// This time solution i for Transient Analytic Solution
		mymaterial->SetTimeValue(TimeValue);
        
		//--------- Resolver usando Runge-Kutta -------------------
        //primeiro estagio de Runge-Kutta
        matM->Multiply(Lastsolution,TotalRhstemp1);
		TotalRhs = TotalRhstemp1 + fvec;
		an.Rhs() = TotalRhs;
		an.Solve();
        
        if(calcresiduo == false){
            //an.Solution().Print("\nsolAntesRG = ");

            //        gradreconst-> ProjectionL2GradientReconstructed(cmesh, matIdL2Proj);
            //        an.LoadSolution(cmesh->Solution());
            //        
            //       // an.Solution().Print("\nsolDeposRG = ");
            //        
            //        //segundo estagio de Runge-Kutta
            //        matMAux->Multiply(Lastsolution,TotalRhstemp1);
            //        Lastsolution = an.Solution();
            //        matM->Multiply(Lastsolution,TotalRhstemp2);
            //        TotalRhs = TotalRhstemp1 + (TotalRhstemp2 + fvec);
            //        TotalRhs = 0.5*TotalRhs;
            //        an.Rhs() = TotalRhs;
            //        an.Solution().Zero();
            //        an.Solve();
            //---------------------------------------------------------

            gradreconst-> ProjectionL2GradientReconstructed(cmesh, matIdL2Proj);
            an.LoadSolution(cmesh->Solution());
            Lastsolution = an.Solution();
        }
        else{
            
            REAL tol = 1.e-8;
            REAL diffsol = 0.;
            TPZFMatrix<STATE> Residuo(nrows,1,0.0);
            TPZFMatrix<STATE> SolIt_k(nrows,1,0.0);
            
            gradreconst-> ProjectionL2GradientReconstructed(cmesh, matIdL2Proj);
            an.LoadSolution(cmesh->Solution());
            Lastsolution = an.Solution();
            
            matMAux->Multiply(Lastsolution,TotalRhstemp2);
            Residuo = TotalRhs - TotalRhstemp2;
//            TotalRhstemp2.Print("sol = \n");
//            TotalRhs.Print("ladodireito = \n");
            
            diffsol = Norm(Residuo);
            int k = 0;
            while(diffsol > tol)
            {
                
                Residuo *=-1.;
                an.Rhs() = Residuo;
                an.Solve();
                SolIt_k = Lastsolution + an.Solution();
                
                CleanGradientSolution(SolIt_k, nonactive);
                an.LoadSolution(SolIt_k);
                cmesh->LoadSolution(SolIt_k);
                
                gradreconst-> ProjectionL2GradientReconstructed(cmesh, matIdL2Proj);
                an.LoadSolution(cmesh->Solution());
                Lastsolution = an.Solution();
                
                matMAux->Multiply(Lastsolution,TotalRhstemp2);
                Residuo = TotalRhs - TotalRhstemp2;
                CleanGradientSolution(Residuo, nonactive);
                diffsol = Norm(Residuo);
                //Lastsolution = SolIt_k;
                
                k++;
            }
        }
        
        //an.Solution().Print("\nsolDeposRG = ");
        
        if(cent%10==0){
            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            PosProcessSolutionTrans(cmesh,an,plotfile);
        }
        
        cent++;
		TimeValue = cent*deltaT;
        //an.Solution().Zero();
	}
 
}

void SolveSistTransient(REAL deltaT, TPZMatConvectionProblem * &mymaterial, TPZCompMesh* cmesh){
	
    TPZAnalysis an(cmesh);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    std::string outputfile;
	outputfile = "SolutionSemReconstGradient";
    
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    
    PosProcessSolution(cmesh,an,plotfile);
    
    
    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(mymaterial, cmesh);
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//        std::stringstream sout;
//        matM->Print("matM = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif

    
    //Criando matriz de rigidez (matK) e vetor de carga
	TPZAutoPointer< TPZMatrix<STATE> > matK;
	TPZFMatrix<STATE> fvec;
    TPZManVector<int64_t> nonactive(0);
    FilterEquation(mymaterial, cmesh, an, true,nonactive);
    matK = an.Solver().Matrix();
    fvec = an.Rhs();
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//        
//        std::stringstream sout;
//        matK->Print("matK = ", sout,EMathematicaInput);
//        fvec.Print("fvec = ", sout,EMathematicaInput);
//        //Print the temporal solution
//        Initialsolution.Print("Initial conditions = ", sout,EMathematicaInput);
//        TPZFMatrix<STATE> Temp;
//        matM->Multiply(Initialsolution,Temp);
//        Temp.Print("Temp matM = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif

    
	int64_t nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
	
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT;
    REAL maxTime;
    mymaterial->GetTimeValue(maxTime);
    
	while (TimeValue <= maxTime)
	{
        ftimeatual  = TimeValue;
        
		// This time solution i for Transient Analytic Solution
		mymaterial->SetTimeValue(TimeValue);
		matM->Multiply(Lastsolution,TotalRhstemp);
        
//#ifdef LOG4CXX
//        if(logdata->isDebugEnabled())
//        {
//            std::stringstream sout;
//            sout<< " tempo = " << cent;
//            Lastsolution.Print("\nIntial conditions = ", sout,EMathematicaInput);
//            TotalRhstemp.Print("Mat Mass x Last solution = ", sout,EMathematicaInput);
//            LOGPZ_DEBUG(logdata,sout.str())
//        }
//#endif
        
		TotalRhs = fvec + TotalRhstemp;
		an.Rhs() = TotalRhs;
		an.Solve();
        Lastsolution = an.Solution();
        
        if(cent%10==0)
        {
                    std::stringstream outputfiletemp;
                    outputfiletemp << outputfile << ".vtk";
                    std::string plotfile = outputfiletemp.str();
                    PosProcessSolution(cmesh,an,plotfile);
        }
        
        cent++;
		TimeValue = cent*deltaT;
	}
}


void PosProcessSolutionTrans(TPZCompMesh* cmesh, TPZAnalysis &an, std::string plotfile)
{
	TPZManVector<std::string,10> scalnames(2), vecnames(0);
	scalnames[0] = "Solution";
    scalnames[1] = "ExactSolution";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}


#include "TPZSkylineNSymStructMatrix.h"
TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZMatConvectionProblem * mymaterial, TPZCompMesh* cmesh){
    
    mymaterial->SetLastState();
    TPZSkylineNSymStructMatrix matsp(cmesh);
	//TPZSpStructMatrix matsp(cmesh);
    
	std::set< int > materialid;
	int matid = mymaterial->MatId();
	materialid.insert(matid);
	matsp.SetMaterialIds (materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<STATE> Un;
    
    TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}

void StiffMatrixLoadVec(TPZMatConvectionProblem *mymaterial, TPZCompMesh*cmesh, TPZAnalysis &an, TPZAutoPointer< TPZMatrix<STATE> > &matK1, TPZFMatrix<STATE> &fvec){
    
	mymaterial->SetCurrentState();
    //TPZSkylineStructMatrix matsk(cmesh);
    TPZBandStructMatrix matsk(cmesh);
	an.SetStructuralMatrix(matsk);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//matK1 = an.StructMatrix();
    matK1 = an.Solver().Matrix();
	fvec = an.Rhs();
}

void ForcingInicial(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    //REAL theta = M_PI/4.0;
    REAL refx = x*cos(teta) + y*sin(teta);
//    REAL refy = -x*sin(teta) + y*cos(teta);
    
    disp[0]=0.;
    if(refx < 0.0) disp[0] = 1.;
    if(refx > 0.0) disp[0] = 0.;
}

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){

    double x = pt[0];
    double y = pt[1];
    
    //REAL theta = M_PI/4.0;
    REAL refx = x*cos(teta) + y*sin(teta);
    
    REAL tp = ftimeatual;
    REAL velx = 1.;
    u[0]=0.;
    
    REAL ptx = refx-velx*tp;
    

    
    if(ptx <= 0.0) u[0] = 1.;
    if(ptx > 0.0) u[0] = 0.;
}

//Ativar apenas a ultima equacao, que corresponde a funcao de base constante
void FilterEquation(TPZMatConvectionProblem *mymaterial, TPZCompMesh *cmesh, TPZAnalysis &an, bool currentstate, TPZManVector<int64_t> &nonactive)
{
    int ncon_saturation = cmesh->NConnects();
    TPZManVector<int64_t> active(0);
    
    for(int i = 0; i<ncon_saturation; i++)
    {
        TPZConnect &con = cmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        if(seqnum==-1) continue;
        int pos = cmesh->Block().Position(seqnum);
        int blocksize = cmesh->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+1);
        int ieq = blocksize-1;
        active[vs] = pos+ieq;
        
    }
    if(currentstate==true)
    {
        mymaterial->SetCurrentState();
        TPZSkylineStructMatrix matsk(cmesh);
        matsk.SetNumThreads(4);
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
        TPZSkylineNSymStructMatrix matst(cmesh);
        //TPZSpStructMatrix matst(mphysics);
        //matsp.SetNumThreads(30);
        matst.EquationFilter().SetActiveEquations(active);
        an.SetStructuralMatrix(matst);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Assemble();
    }
    
    
    for(int i = 0; i<ncon_saturation; i++)
    {
        TPZConnect &con = cmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = cmesh->Block().Position(seqnum);
        int blocksize = cmesh->Block().Size(seqnum);
        int vs = nonactive.size();
        nonactive.Resize(vs+blocksize-1);
        for(int ieq = 0; ieq<blocksize-1; ieq++)
		{
            nonactive[vs+ieq] = pos+ieq;
        }
    }
}

void CleanGradientSolution(TPZFMatrix<STATE> &Solution, TPZManVector<int64_t> &Gradients)
{
    for(int i=0; i < Gradients.size(); i++ )
    {
        Solution(Gradients[i],0)=0.0;
    }
}

//----------------------------------------------------------------------------
#include "pzgengrid.h"
TPZGeoMesh *GMesh2(REAL Lx, REAL Ly, bool triang_elements){
//    TPZManVector<int,2> nx(2,1);
//
//    TPZManVector<REAL,3> x0(3,0.),x1(3,2.0);
//    TPZGenGrid gengrid(nx,x0,x1);
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    if(triang_elements)
//    {
//        gengrid.SetElementType(1);
//    }
//    gengrid.Read(gmesh);
//    
//    gengrid.SetBC(gmesh,4,bc0);
//    gengrid.SetBC(gmesh,5,bc1);
//    gengrid.SetBC(gmesh,6,bc2);
//    gengrid.SetBC(gmesh,7,bc3);
//    
//    gmesh->BuildConnectivity();
    
    
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
	REAL valx, dx=Lx;
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
		valx = Ly - xi*dx;
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
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
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

TPZCompMesh *MalhaComp2(TPZGeoMesh * gmesh, int pOrder/*,TPZMatConvectionProblem * &material*/)
{
	/// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    REAL diff= 0.;
    REAL conv =1.;
    TPZVec<REAL> convdir;
    convdir.Resize(3,0.);
    convdir[0]=1.;
    material-> SetParameters(diff, conv, convdir);
    material->SetSD(0.);
    material->SetNoPenalty();
    material->SetNonSymmetric();
    //material->SetSolutionPenalty();
    //material->SetSymmetric();
    TPZMaterial * mat(material);
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
//    int dim = 2;
//	material = new TPZMatConvectionProblem(matId,dim);
//	TPZMaterial * mat(material);
//	
//	TPZVec<REAL> convdir(dim,0.);
//    convdir[0] = 1.;
//	REAL flux = 0.;
//    REAL rho = 1.;
//	
//	material->SetParameters(rho,convdir);
//	material->SetInternalFlux( flux);
//	material->NStateVariables();
//	
//	TPZCompEl::SetgOrder(pOrder);
//	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//	cmesh->SetDimModel(dim);
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZAutoPointer<TPZFunction<STATE> > myforce = new TPZDummyFunction<STATE>(ConvGradU, 5);
    material->SetForcingFunction(myforce);
    
    TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(SolucaoExata, 5);
	material->SetForcingFunctionExact(solExata);
    
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,1, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,1, val1, val2);
    
    TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(DirichletCond, 5);
    TPZAutoPointer<TPZFunction<STATE> > fCC3 = new TPZDummyFunction<STATE>(DirichletCond, 5);
    
    //val2(0,0) =  0.0886226925452758;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,3, val1, val2);
    
    //val2(0,0) = -0.0886226925452758;
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,0, val1, val2);
    
    BCond1->SetForcingFunction(fCC1);
    BCond3->SetForcingFunction(fCC3);
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,material->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
    
    cmesh->SetDefaultOrder(pOrder);
	
	//Ajuste da estrutura de dados computacional
    cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
//    int nel = cmesh->NElements();
//    for(int i=0; i<nel; i++){
//        TPZCompEl *cel = cmesh->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        celdisc->SetConstC(1.);
//        celdisc->SetCenterPoint(0, 0.);
//        celdisc->SetCenterPoint(1, 0.);
//        celdisc->SetCenterPoint(2, 0.);
//        celdisc->SetTrueUseQsiEta();
//        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//        {
//            if(triang==true) celdisc->SetTotalOrderShape();
//            else celdisc->SetTensorialShape();
//        }
//    }
    
	return cmesh;
}


void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp/*, TPZFMatrix<STATE> &deriv*/){
    
    double x = pt[0];
    double aux1 = -valC*(1 - 2.*x + x*x);
    double aux2 = exp(aux1);
    disp[0]=-(valC)*(x - 1.)*aux2;
}

void ForcingInicial2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	REAL x = pt[0];
    
    REAL temp1 = 0.;
    REAL temp2 = sqrt(valC);
    
#ifdef USING_BOOST
    temp1 = 0.5*sqrt(M_PI)*(boost::math::erf(temp2*(x-1.)));
    temp1 /=temp2;
#else
    temp1 = 0.5*sqrt(M_PI)*erf(temp2*(x-1.));
    temp1 /=temp2;
#endif
    
    disp[0] = temp1;
    //if(x>=1.0) disp[0] = temp1;
}


void ConvGradU(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &deriv){
    
    double x = pt[0];
    disp[0]=0.;
    REAL val = -valC*(x-1.)*(x-1.);
    disp[0]=-exp(val);
}

void SolucaoExata(const TPZVec<REAL> &pt, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv){
    
    deriv(0,0)=0.;
    deriv(1,0)=0.;
    sol[0]=0.;
    
    
    double x = pt[0];
    //double y = pt[1];
    REAL tp = ftimeatual;
    REAL vel = 1.;
    REAL temp = sqrt(valC);
    
    REAL px = x-vel*tp;
#ifdef USING_BOOST
    sol[0] = 0.5*sqrt(M_PI)*(boost::math::erf(temp*(px-1.)));
    sol[0] /= temp;
#else
    sol[0] = 0.5*sqrt(M_PI)*erf(temp*(px-1.));
    sol[0] /= temp;
#endif
    REAL val = -valC*(px-1.)*(px-1.);
    deriv(0,0)=exp(val);
}

#include "TPZSkylineNSymStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "TPZFrontNonSym.h"
void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
    TPZSkylineNSymStructMatrix full(Cmesh);//caso nao-simetrico
    //TPZBandStructMatrix full(Cmesh);
    //TPZSkylineStructMatrix full(Cmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELU);//caso nao simetrico
	an.SetSolver(step);
	an.Run();
	
//	//Saida de Dados: solucao e  grafico no VT
//	ofstream file("Solutout");
//	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

//Ativar apenas a ultima equacao, que corresponde a funcao de base constante
void FilterEquation(TPZCompMesh *cmesh, TPZAnalysis &an)
{
    int ncon_saturation = cmesh->NConnects();
    TPZManVector<int64_t> active(0);
    for(int i = 0; i<ncon_saturation; i++)
    {
        TPZConnect &con = cmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        if(seqnum==-1) continue;
        int pos = cmesh->Block().Position(seqnum);
        int blocksize = cmesh->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+1);
        int ieq = blocksize-1;
        active[vs] = pos+ieq;
    }
   
    TPZSkylineNSymStructMatrix matst(cmesh);
    //TPZBandStructMatrix matst(cmesh);
    //TPZSkylineStructMatrix matst(cmesh);
    matst.SetNumThreads(5);
    matst.EquationFilter().SetActiveEquations(active);
    an.SetStructuralMatrix(matst);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    an.Assemble();
}

void DirichletCond(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
   
    //TPZManVector<REAL> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
    SolucaoExata(loc,result,du);
}

void PosProcessSolution(TPZCompMesh* cmesh, TPZAnalysis &an, std::string plotfile)
{
	TPZManVector<std::string,10> scalnames(2), vecnames(0);
	scalnames[0] = "Solution";
    scalnames[1] = "ExactSolution";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}



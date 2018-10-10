
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

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

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

#include "pzhdivfull.h"

#include <iostream>
#include <math.h>
using namespace std;

int matId =1;

int dirichlet = 0;
int neumann = 1;

int bc0=-1;
int bc1=-2;
int bc2=-3;
int bc3=-4;


TPZGeoMesh *GMesh(bool triang_elements, REAL Lx, REAL Ly);

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcessFlux(TPZAnalysis &an, std::string plotfile);

//solucao exata
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);

//lado direito da equacao
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno de Neumann
void ForcingBC0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

REAL const Pi = 4.*atan(1.);
// nao esta rodando com estas configuracoes..aguardar Agnaldo
bool ftriang = true;
bool isStab = false;
bool iscontinuou = false;
bool useh2 = false;
REAL delta1 = 0.5;
REAL delta2 = 0.5;
bool isFullHdiv=false;
 /*
bool ftriang = true;//seta polinomios completos ou totais
bool isStab = false;//ativa ou nao estabilizacao
bool iscontinuou = false;//ativa h1 ou l2 para pressao
bool useh2 = false;//ativa o termo h2 para penalizacao
REAL delta1 = 0.5;//seta os valores de delta 1 e 2 para a penalizacao
REAL delta2 = 0.5;
bool isFullHdiv=false;//seta espaco completo ou nao para o fluxo
  */

void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out);
void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out);

STATE EPSILON = 1000.;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.main"));
#endif

#include "pztransfer.h"
int main2(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    REAL Lx = 1.;
    REAL Ly = 1.;
    
    {// Checando funcoes
        TPZVec<REAL> pto(2,0.);
        TPZVec<STATE> disp(1,0.),solp(1,0.);
        pto[0]=0.200935;
        pto[1]=0.598129;
        TPZFMatrix<STATE>flux(3,1,0.);
        SolExata(pto, solp,flux);
        Forcing(pto, disp);
        
        solp.Print(std::cout);
        flux.Print(std::cout);
        disp.Print(std::cout);
    }
    
   // ofstream saidaerro("../ErroPoissonHdivMalhaQuad.txt",ios::app);
    ofstream saidaerro("../ErroPoissonHdivMalhaTriang.txt",ios::app);

    for(int p = 1; p<2; p++)
    {
        int pq = p;
        int pp;
        if(ftriang==true){
            pp = pq-1;
        }else{
            pp = pq;
        }
        
        int ndiv;
        saidaerro<<"\n CALCULO DO ERRO, COM ORDEM POLINOMIAL pq = " << pq << " e pp = "<< pp <<endl;
        for (ndiv = 1; ndiv< 2; ndiv++)
        {
            
            //std::cout << "p order " << p << " number of divisions " << ndiv << std::endl;
            
            saidaerro<<"\n<<<<<< Numero de divisoes uniforme ndiv = " << ndiv <<" >>>>>>>>>>> "<<endl;
            
            TPZGeoMesh *gmesh = GMesh(ftriang, Lx, Ly);
            ofstream arg("gmesh1.txt");
            gmesh->Print(arg);
            
            UniformRefine(gmesh, ndiv);
            
            TPZCompMesh *cmesh1 = CMeshFlux(gmesh, pq);
            TPZCompMesh *cmesh2 = CMeshPressure(gmesh, pp);
            
//            ofstream arg1("cmeshflux.txt");
//            cmesh1->Print(arg1);
//            
//            ofstream arg2("cmeshpressure.txt");
//            cmesh2->Print(arg2);
//            
//            ofstream arg4("gmesh2.txt");
//            gmesh->Print(arg4);
            
            
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            
            TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
            
            std::cout << "Number of equations " << mphysics->NEquations() << std::endl;
            int numthreads = 8;
            std::cout << "Number of threads " << numthreads << std::endl;

            TPZAnalysis an(mphysics);
            SolveSyst(an, mphysics,numthreads);
            
            
//            ofstream arg5("cmeshmultiphysics.txt");
//            mphysics->Print(arg5);

            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            TPZVec<REAL> erros;
    
            saidaerro << "Valor de epsilone " << EPSILON << std::endl;
            saidaerro << "Numero de threads " << numthreads << std::endl;
            saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
            ErrorHDiv(cmesh1, saidaerro);
            
            saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
            ErrorL2(cmesh2, saidaerro);
            
            std::cout << "Postprocessed\n";
            
            //Plot da solucao aproximada
//            string plotfile("Solution_mphysics.vtk");
//            PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
        }
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
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
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
    	if(logger->isDebugEnabled())
    	{
            std::stringstream sout;
            sout<<"\n\n Malha Geometrica Inicial\n ";
            gmesh->Print(sout);
          LOGPZ_DEBUG(logger,sout.str())
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
    cmesh->SetDimModel(dim);
	
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
	

    if(isFullHdiv){
        cmesh->SetAllCreateFunctionsHDivFull();
    }
    else{
		cmesh->SetAllCreateFunctionsHDiv();
    }
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
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
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    if(iscontinuou==false)
    {
        cmesh->SetAllCreateFunctionsDiscontinuous();

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
                if(ftriang==true) celdisc->SetTotalOrderShape();
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
    }
    else{
        cmesh->SetAllCreateFunctionsContinuous();
        
        //Ajuste da estrutura de dados computacional
        cmesh->AutoBuild();
    }
    
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

TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =2;
    TPZMixedPoisson *material = new TPZMixedPoisson(matId,dim);
    
    //incluindo os dados do problema
    REAL coefk = 1.;
    REAL coefvisc = 1.;
    material->SetPermeability(coefk);
    material->SetViscosity(coefvisc);
    
    if(isStab==true){
        material->SetStabilizedMethod();
        material->SetStabilizationCoeficients(delta1,delta2);
	}
    if(isStab==true && useh2==true) material->SetHdois();
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing, 5);
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    
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

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads)
{
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
    full.SetNumThreads(numthreads);
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

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(2), vecnames(3);
	vecnames[0]  = "Flux";
    vecnames[1]  = "GradFluxX";
    vecnames[2]  = "GradFluxY";
    scalnames[0] = "Pressure";
    scalnames[1] = "DivFlux";
    
	const int dim = 2;
	int div =1;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
//	std::ofstream out("malha.txt");
//	an.Print("nothing",out);
    
}

void SolExata2(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    solp[0] = sin(Pi*x)*sin(Pi*y);
    flux(0,0)=-Pi*cos(Pi*x)*sin(Pi*y);
    flux(1,0)=-Pi*cos(Pi*y)*sin(Pi*x);
    flux(2,0)=2*Pi*Pi*sin(Pi*y)*sin(Pi*x);
}

#define Power pow
#define ArcTan atan
#define Sqrt sqrt

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1);
    REAL eps = EPSILON;
    REAL x = pt[0];
    REAL y = pt[1];
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    solp[0] = 8*(1. - x)*x*(1. - y)*y*ArcTan(0.0625 +
                                             2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)));
    STATE a = (-64*Sqrt(eps)*(1. - x)*(-0.5 + x)*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) - (32*Sqrt(eps)*(1. - x)*x*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) + (64*Sqrt(eps)*(-0.5 + x)*x*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) - (256*eps*(1. - x)*Power(-0.5 + x,2)*x*
                      (0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))*(1. - y)*
                      y)/Power(1 + Power(0.0625 + 2*Sqrt(eps)*
                                         (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),2),2) -
    16*(1. - y)*y*ArcTan(0.0625 + 2*Sqrt(eps)*
                         (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)));
    STATE b = (-64*Sqrt(eps)*(1. - x)*x*(1. - y)*(-0.5 + y))/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) - (32*Sqrt(eps)*(1. - x)*x*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) + (64*Sqrt(eps)*(1. - x)*x*(-0.5 + y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) - (256*eps*(1. - x)*x*(0.0625 +
                                          2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))*(1. - y)*
                      Power(-0.5 + y,2)*y)/
    Power(1 + Power(0.0625 + 2*Sqrt(eps)*
                    (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),2),2) -
    16*(1. - x)*x*ArcTan(0.0625 + 2*Sqrt(eps)*
                         (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)));
    flux(2,0) = -a-b;
    flux(0,0) = (-32*Sqrt(eps)*(1. - x)*(-0.5 + x)*x*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) + 8*(1. - x)*(1. - y)*y*
    ArcTan(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2))) -
    8*x*(1. - y)*y*ArcTan(0.0625 + 2*Sqrt(eps)*
                          (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)));
    flux(1,0) = (-32*Sqrt(eps)*(1. - x)*x*(1. - y)*(-0.5 + y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) + 8*(1. - x)*x*(1. - y)*
    ArcTan(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2))) -
    8*(1. - x)*x*y*ArcTan(0.0625 + 2*Sqrt(eps)*
                          (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)));
    flux(0,0) *= -1.;
    flux(1,0) *= -1.;
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    REAL eps = EPSILON;
    STATE a = (-64*Sqrt(eps)*(1. - x)*(-0.5 + x)*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) - (32*Sqrt(eps)*(1. - x)*x*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) + (64*Sqrt(eps)*(-0.5 + x)*x*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) - (256*eps*(1. - x)*Power(-0.5 + x,2)*x*
                      (0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))*(1. - y)*
                      y)/Power(1 + Power(0.0625 + 2*Sqrt(eps)*
                                         (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),2),2) -
    16*(1. - y)*y*ArcTan(0.0625 + 2*Sqrt(eps)*
                         (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)));
    STATE b = (-64*Sqrt(eps)*(1. - x)*x*(1. - y)*(-0.5 + y))/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) - (32*Sqrt(eps)*(1. - x)*x*(1. - y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) + (64*Sqrt(eps)*(1. - x)*x*(-0.5 + y)*y)/
    (1 + Power(0.0625 + 2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),
               2)) - (256*eps*(1. - x)*x*(0.0625 +
                                          2*Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))*(1. - y)*
                      Power(-0.5 + y,2)*y)/
    Power(1 + Power(0.0625 + 2*Sqrt(eps)*
                    (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)),2),2) -
    16*(1. - x)*x*ArcTan(0.0625 + 2*Sqrt(eps)*
                         (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)));

    disp[0]= -a-b;
}

void Forcing2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

void ForcingBC0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    disp[0]= Pi*cos(Pi*y)*sin(Pi*x);
}

void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    disp[0]= -Pi*cos(Pi*y)*sin(Pi*x);
}

void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
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
        cel->EvaluateError(SolExata, elerror, false);
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

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out)
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
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, false);
        int nerr = elerror.size();
        globerrors.resize(nerr);
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
    out << "Errors associated with L2 space\n";
    out << "L2 Norm = "    << sqrt(globerrors[1]) << endl;
}

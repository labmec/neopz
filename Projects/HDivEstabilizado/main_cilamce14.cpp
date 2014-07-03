
#ifdef HAVE_CONFIG_H
#include <config.h>
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

int MatId =1;

int bcdirichlet = 0;
int bcneumann = 1;

int BC0=-1;
int BC1=-2;
int BC2=-3;
int BC3=-4;


TPZGeoMesh *GMesh2(REAL Lx, REAL Ly,bool triang_elements);

TPZCompMesh *CMeshFlux2(int pOrder, TPZGeoMesh *gmesh);
TPZCompMesh *CMeshPressure2(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *CMeshMixed2(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh);

void RefinamentoUnif(TPZGeoMesh* gmesh, int nDiv);
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads);
void PosProcessMultph(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcessFluxo(TPZAnalysis &an, std::string plotfile);

//solucao exata
void SolExataMista(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void SolExataPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);

//lado direito da equacao
void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno de Neumann
void NeumannAcima(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAbaixo(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

REAL const pi = 4.*atan(1.);
// nao esta rodando com estas configuracoes..aguardar Agnaldo
bool fTriang = true;
bool IsStab = false;
bool IsContinuou = false;
bool Useh2 = false;
REAL Delta1 = 0.5;
REAL Delta2 = 0.5;
bool IsFullHdiv=false;

//erros
void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out);
void ErrorL22(TPZCompMesh *l2mesh, std::ostream &out);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.main"));
#endif

#include "pztransfer.h"
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    REAL Lx = 1.;
    REAL Ly = 0.5;
    
    ofstream saidaerro( "../erros-hdiv-estab.txt",ios::app);

    for(int p = 2; p<3; p++)
    {
        int pq = p;
        int pp;
        if(fTriang==true){
            pp = pq-1;
        }else{
            pp = pq;
        }
        
        int ndiv;
        saidaerro<<"\n CALCULO DO ERRO, COM ORDEM POLINOMIAL pq = " << pq << " e pp = "<< pp <<endl;
        for (ndiv = 4; ndiv< 5; ndiv++)
        {
            
            //std::cout << "p order " << p << " number of divisions " << ndiv << std::endl;
            
            saidaerro<<"\n<<<<<< Numero de divisoes uniforme ndiv = " << ndiv <<" >>>>>>>>>>> "<<endl;
            
            TPZGeoMesh *gmesh = GMesh2(Lx, Ly,fTriang);
           // ofstream arg("gmesh1.txt");
            //gmesh->Print(arg);
            
            RefinamentoUnif(gmesh, ndiv);
            
            TPZCompMesh *cmesh1 = CMeshFlux2(pq,gmesh);
            TPZCompMesh *cmesh2 = CMeshPressure2(pp,gmesh);
            
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
            
            TPZCompMesh * mphysics = CMeshMixed2(meshvec,gmesh);
            
            std::cout << "Number of equations " << mphysics->NEquations() << std::endl;
            int numthreads = 8;
            std::cout << "Number of threads " << numthreads << std::endl;

            TPZAnalysis an(mphysics);
            ResolverSistema(an, mphysics,numthreads);
            
            
//            ofstream arg5("cmeshmultiphysics.txt");
//            mphysics->Print(arg5);

            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            TPZVec<REAL> erros;
    
            saidaerro<<"\nErro da simulacao multifisica do fluxo (q)" <<endl;
            ErrorHDiv2(cmesh1, saidaerro);
            
            saidaerro<<"Erro da simulacao multifisica da pressao (p)" <<endl;
            ErrorL22(cmesh2, saidaerro);
            
            std::cout << "Postprocessed\n";
            
            //Plot da solucao aproximada
            //string plotfile("Solution_mphysics.vtk");
           //PosProcessMultph(meshvec,  mphysics, an, plotfile);
        }
    }
    
	return EXIT_SUCCESS;
}

TPZGeoMesh *GMesh2(REAL Lx, REAL Ly,bool triang_elements){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolTriang(3);
	TPZVec <long> TopolLine(2);
    TPZVec <long> TopolPoint(1);
	
	//indice dos nos
	long id = 0;
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
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,MatId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,MatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
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

void RefinamentoUnif(TPZGeoMesh* gmesh, int nDiv)
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


TPZCompMesh *CMeshFlux2(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(MatId,dim);
	TPZMaterial * mat(material);
	material->NStateVariables();
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
	
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
	

    if(IsFullHdiv){
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

TPZCompMesh *CMeshPressure2(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(MatId,dim);
	material->NStateVariables();
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    if(IsContinuou==false)
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
                if(fTriang==true) celdisc->SetTotalOrderShape();
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

TPZCompMesh *CMeshMixed2(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =2;
    TPZMixedPoisson *material = new TPZMixedPoisson(MatId,dim);
    
    //incluindo os dados do problema
    REAL coefk = 1.;
    REAL coefvisc = 1.;
    material->SetPermeability(coefk);
    material->SetViscosity(coefvisc);
    
    if(IsStab==true){
        material->SetStabilizedMethod();
        material->SetStabilizationCoeficients(Delta1,Delta2);
	}
    if(IsStab==true && Useh2==true) material->SetHdois();
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExataMista);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingMista);
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
    BCond0 = material->CreateBC(mat, BC0,bcneumann, val1, val2);
    BCond1 = material->CreateBC(mat, BC1,bcneumann, val1, val2);
    BCond2 = material->CreateBC(mat, BC2,bcneumann, val1, val2);
    BCond3 = material->CreateBC(mat, BC3,bcneumann, val1, val2);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
    bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcima);
    BCond2->SetForcingFunction(bcmatNeumannAcima);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAbaixo;
    bcmatNeumannAbaixo = new TPZDummyFunction<STATE>(NeumannAbaixo);
    BCond0->SetForcingFunction(bcmatNeumannAbaixo);
    
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

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads)
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

void PosProcessMultph(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(3), vecnames(4);
	vecnames[0]  = "Flux";
    vecnames[1]  = "GradFluxX";
    vecnames[2]  = "GradFluxY";
    vecnames[3]  = "ExactFlux";
    scalnames[0] = "Pressure";
    scalnames[1] = "DivFlux";
    scalnames[2] = "ExactPressure";
    
	const int dim = 2;
	int div =1;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
//	std::ofstream out("malha.txt");
//	an.Print("nothing",out);
    
}

void SolExataMista(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = (-2./pi)*cos(pi*x)*exp(y/2.);
    flux(0,0)= -2*sin(pi*x)*exp(y/2.);
    flux(1,0)= (1./pi)*cos(pi*x)*exp(y/2.);
    flux(2,0)= (1.- 4.*pi*pi)/(2.*pi)*(cos(pi*x)*exp(y/2.));
}

void SolExataPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = (-2./pi)*cos(pi*x)*exp(y/2.);
    flux(0,0)= 2*sin(pi*x)*exp(y/2.);
    flux(1,0)= -(1./pi)*cos(pi*x)*exp(y/2.);
}

void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= (1.- 4.*pi*pi)/(2.*pi)*(cos(pi*x)*exp(y/2.));
}

void NeumannAbaixo(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    //double y = pt[1];
    disp[0] = -cos(pi*x)/pi;
}

void NeumannAcima(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    //double y = pt[1];
    disp[0] = cos(pi*x)/pi;
}

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out)
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
        cel->EvaluateError(SolExataMista, elerror, NULL);
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

void ErrorL22(TPZCompMesh *l2mesh, std::ostream &out)
{
    long nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
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
        cel->EvaluateError(SolExataPressao, elerror, NULL);
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
    out << "\n";
    out << "Errors associated with L2 or H1 space\n";
    out << "H1 Norm = "    << sqrt(globerrors[0]) << endl;
    out << "L2 Norm = "    << sqrt(globerrors[1]) << endl;
    out << "Semi H1 Norm = "    << sqrt(globerrors[2]) << endl;
}

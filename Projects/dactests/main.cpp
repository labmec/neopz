
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
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
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

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "TPZReadGIDGrid.h"
#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

#include "pzhdivfull.h"

#include "pzgeopyramid.h"

#include "PZMatPoissonD3.h"
#include "pzcompelhdivcurved.h"


#include <iostream>
#include <math.h>
using namespace std;

int matId = 1;

int dirichlet = 0;
int neumann = 1;

int bc0 = -1;
int bc1 = -2;
int bc2 = -3;
int bc3 = -4;


TPZGeoMesh *GMesh(int dimensao, int tipo, REAL Lx, REAL Ly, REAL Lz);

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcessFlux(TPZAnalysis &an, std::string plotfile);

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv);
void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);

//solucao exata
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);

//lado direito da equacao
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno de Neumann
void ForcingBC0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

REAL const Pi = 4.*atan(1.);

bool ftriang = false;
bool iscontinuous = false;
int Maps();

/*
 bool ftriang = true;//seta polinomios completos ou totais
 bool isStab = false;//ativa ou nao estabilizacao
 bool iscontinuou = false;//ativa h1 ou l2 para pressao
 bool useh2 = false;//ativa o termo h2 para penalizacao
 REAL delta1 = 0.5;//seta os valores de delta 1 e 2 para a penalizacao
 REAL delta2 = 0.5;
 bool isFullHdiv=false;//seta espaco completo ou nao para o fluxo
 */

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif


#include "pztransfer.h"
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    //int outcome = Maps();
    
    
    REAL Lx = 1.; // limite inferior e superior do intervalo em x
    REAL Ly = 1.; // limite inferior e superior do intervalo em y
    REAL Lz = 1.; // limite inferior e superior do intervalo em z
    int p = 1;
    int ndiv = 0;
    
    
    // Para dimensao 2
    // tipo 1 triangulo
    // tipo 2 quadrilatero
    
    int dimensao = 2;
    int tipo = 2;
    ofstream saidaerro("../ErroPoissonHdivMalhaQuad.txt",ios::app);
    //int tipo = 1;
    //ofstream saidaerro("../ErroPoissonHdivMalhaTriang.txt",ios::app);
    
    for(p=1;p<2;p++)
    {
        int pq = p;
        int pp;
        if(ftriang==true){
            pp = pq-1;
        }else{
            pp = pq;
        }

        for (ndiv=0; ndiv<1; ndiv++)
        {
            
            TPZGeoMesh *gmesh = GMesh(dimensao, tipo, Lx, Ly, Lz);
            ofstream arg("gmesh1.txt");
            gmesh->Print(arg);
            // gmesh->
            
            UniformRefine(gmesh, ndiv);
            
            {
                //	Print Geometrical Base Mesh
                std::ofstream Dummyfile("GeometricMesh.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
            }
            
            TPZCompMesh *cmesh1 = CMeshFlux(gmesh, pq);
            TPZCompMesh *cmesh2 = CMeshPressure(gmesh, pp);
            
            ofstream arg1("cmeshflux.txt");
            cmesh1->Print(arg1);
            
            ofstream arg2("cmeshpressure.txt");
            cmesh2->Print(arg2);
            
            ofstream arg4("gmesh2.txt");
            gmesh->Print(arg4);
            
            
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            
            TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
            ofstream arg5("cmeshmultiphysics.txt");
            mphysics->Print(arg5);
            
            
            TPZAnalysis an(mphysics);
            string plotfile("Solution_mphysics.vtk");
            //    PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
            
            SolveSyst(an, mphysics);
            
            
            //
            //    TPZAutoPointer< TPZMatrix<REAL> > matKdoAna;
            //    matKdoAna = an.Solver().Matrix();
            //
            //#ifdef LOG4CXX
            //    if(logdata->isDebugEnabled())
            //    {
            //        std::stringstream sout;
            //        matKdoAna->Print("Eglobal = ", sout,EMathematicaInput);
            //        an.Rhs().Print("Fglobal = ", sout,EMathematicaInput);
            //        an.Solution().Print("Sglobal = ", sout,EMathematicaInput);
            //        LOGPZ_DEBUG(logdata,sout.str())
            //    }
            //#endif
            
            
            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            TPZVec<REAL> erros;
            
//            TPZManVector<REAL,3> myerrors(3,0.);
//            an.SetExact(SolExata);
//            an.PostProcessError(myerrors);
            
//            saidaerro << "Valor de epsilone " << EPSILON << std::endl;
//            saidaerro << "Numero de threads " << numthreads << std::endl;
            saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
            ErrorHDiv(cmesh1, saidaerro, p, ndiv);
            
            saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
            ErrorL2(cmesh2, saidaerro, p, ndiv);
            
            std::cout << "Postprocessed\n";
            
            //Plot da solucao aproximada
            // para estudo do ero nao precisa gerar dados de saida
            //PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
            
            std::cout<< " grau  polinomio " << p << " numero de divisoes " << ndiv << std::endl;
         }
        
    }
    
    std::cout<< " fim " << std::endl;
    
	return EXIT_SUCCESS;
}

TPZGeoMesh *GMesh(int d, int tipo, REAL Lx, REAL Ly, REAL Lz){
    
    int dim = d;
    if(d<2 || d>3)
    {
        std::cout << "dimensao errada" << std::endl;
        dim = 2;
        DebugStop();
    }
    
    int Qnodes = dim == 2 ? 4 : 8;
    
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
    
    if(tipo==1) // triangulo
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
    //	gmesh->BuildConnectivity();
}


TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
	TPZMaterial * mat(material);
	material->NStateVariables();
    
    //  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
	
    cmesh->SetDimModel(dim);
    
//    cmesh->SetAllCreateFunctionsHDivFull();
    cmesh->SetAllCreateFunctionsHDiv();
    
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
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
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    if(iscontinuous == true){
        cmesh->SetAllCreateFunctionsContinuous(); 
    }
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    if (iscontinuous==false) {
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
            //celdisc->SetFalseUseQsiEta();
            
//            TPZVec<REAL> qsi(3,0.);
//            qsi[0] = 0.5;
//            qsi[1] = 0.5;
//            TPZFMatrix<REAL> phi;
//            TPZFMatrix<REAL> dphi;
//            celdisc->Shape(qsi, phi,dphi);
//            phi.Print("phi = ");
            
            
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                if(ftriang==true) celdisc->SetTotalOrderShape();
                else celdisc->SetTensorialShape();
            }
            
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

TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =2;
    TPZMatPoissonD3 *material = new TPZMatPoissonD3(matId,dim);
    
    //incluindo os dados do problema
    //incluindo os dados do problema
    REAL coefk = 1.;
    REAL coefvisc = 1.;
//    material->SetPermeability(coefk);
//    material->SetViscosity(coefvisc);
//    
//    if(isStab==true){
//        material->SetStabilizedMethod();
//        material->SetStabilizationCoeficients(delta1,delta2);
//	}
//    if(isStab==true && useh2==true) material->SetHdois();
//    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing);
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
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
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
    
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    
	// Creation of interface elements
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
    
    return mphysics;
}

void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
//    TPZSkylineNSymStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt); //caso simetrico
//	step.SetDirect(ELU);
	an.SetSolver(step);
//    an.Assemble();
	an.Run();
    
      // Nao entendi como isso funciona, nao consegui fazer o erro do material ser usado. Falta algo?
//    an.SetExact(SolExata);
//    TPZVec<STATE> error(2,0.);
//    an.PostProcessError(error);
    
    
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(3), vecnames(3);
	vecnames[0]  = "Flux";
    vecnames[1]  = "GradFluxX";
    vecnames[2]  = "GradFluxY";
    scalnames[0] = "Pressure";
    scalnames[1] = "DivFlux";
    scalnames[2] = "ExactPressure";
    
	const int dim = 2;
	int div =2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
    
//    mphysics->Solution().Print("Solucao");
    
}

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
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

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0] = 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
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

int Maps()
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(2);
    
    // Create nodes
    int MatID = 1;
    int nodeId=0;
    
    long index = 0;
    TPZVec<REAL> coord(3,0.);
    TPZVec<long> TopologyPoint(1);
    TPZVec<long> TopologyLinear(2);
    
    coord[0]=0.0; // xcoord
    coord[1]=0.0; // ycoord
    coord[2]=0.0; // zcoord
    gmesh->NodeVec()[nodeId].SetCoord(coord);
    gmesh->NodeVec()[nodeId].SetNodeId(nodeId);
    TopologyPoint[0] = nodeId;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > *Mypoint1=new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (nodeId++, TopologyPoint, MatID,*gmesh);
    
    coord[0]=1.0; // xcoord
    coord[1]=1.0; // ycoord
    coord[2]=1.0; // zcoord
    gmesh->NodeVec()[nodeId].SetCoord(coord);
    gmesh->NodeVec()[nodeId].SetNodeId(nodeId);
    TopologyPoint[0] = nodeId;
    TPZGeoElRefPattern < pzgeom::TPZGeoPoint > *Mypoint2=new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (nodeId++, TopologyPoint, MatID,*gmesh);
    
    // Creating line
    TopologyLinear[0] = 0;
    TopologyLinear[1] = 1;
    TPZGeoElRefPattern < pzgeom::TPZGeoLinear > *Myline=new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (nodeId++, TopologyLinear, 2,*gmesh);
    
    
    
    
    TPZVec<REAL> qsi(1,1.0);
    TPZVec<REAL> coordX(3,1.0);
    
    Myline->X(qsi, coordX);
    std::cout << "Xcoord = " << coordX[0] << std::endl;
    std::cout << "Ycoord = " << coordX[1] << std::endl;
    std::cout << "Zcoord = " << coordX[2] << std::endl;
    TPZFMatrix<REAL> jac,axes,jacinverse;
    REAL detjac;
    
    Myline->Jacobian(qsi, jac, axes, detjac, jacinverse);
    
    std::cout << "jacx = " << jac(0,0) << std::endl;
    std::cout << "jacinversex = " << jacinverse(0,0) << std::endl;    
    
//	std::string GeoGridFile;
//    //	GeoGridFile = "SQDomain.dump";
//	GeoGridFile = "MeshTest.dump";
//    TPZReadGIDGrid GIDMesh;
//    TPZGeoMesh *gmesh = GIDMesh.GeometricGIDMesh(GeoGridFile);
//    
//	{
//		//	Print Geometrical Base Mesh
//		std::ofstream Dummyfile("GeometricGID.vtk");
//		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
//	}
    
    gmesh->BuildConnectivity();
    gmesh->Print();
    
    return 0;
    
}


void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv)
{
    long nel = hdivmesh->NElements();
    //int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        TPZManVector<STATE,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    out << "L2 Norm for flux - L2 Norm for divergence - Hdiv Norm for flux " << endl;
    out <<  setw(16) << sqrt(globalerrors[1]) << setw(25)  << sqrt(globalerrors[2]) << setw(21)  << sqrt(globalerrors[3]) << endl;

//    
//    out << "L2 Norm for flux = "    << sqrt(globalerrors[1]) << endl;
//    out << "L2 Norm for divergence = "    << sqrt(globalerrors[2])  <<endl;
//    out << "Hdiv Norm for flux = "    << sqrt(globalerrors[3])  <<endl;
//    
}

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv)
{
    long nel = l2mesh->NElements();
    //int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<STATE,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, NULL);
        int nerr = elerror.size();
        globalerrors.resize(nerr);
#ifdef LOG4CXX
        if (logdata->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logdata, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with L2 space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    out << "L2 Norm = "    << sqrt(globalerrors[1]) << endl;
}

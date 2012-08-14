//
//  main.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/20/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>
using namespace std;

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzinterpolationspace.h"

#include "TPZCompElDisc.h"
#include "pzpoisson3d.h"
#include "pzelasmat.h"
#include "pzporoelasticmf2d.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"

#include "pzfunction.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>
using namespace std;

int const matId =1;
const int bcBottom = -1;
const int bcRight = -2;
const int bcTop = -3;
const int bcLeft = -4;

int const dirichlet =0;
int const neumann = 1;
int const neumdir=10;
int const dirfreey_neum=300;
int const dirneum = 1;
int const mixedneum = 21;

REAL const Pi = 4.*atan(1.);

TPZGeoMesh *GMesh(bool triang_elements, REAL w, REAL L);
TPZCompMesh *MalhaCompElast(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZPoroElasticMF2d * &mymaterial);

void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);

TPZAutoPointer <TPZMatrix<REAL> > MassMatrix(TPZPoroElasticMF2d *mymateria, TPZCompMesh *mphysics);
void StiffMatrixLoadVec(TPZPoroElasticMF2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec);

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);

void SolveSistTransient(REAL deltaT,REAL maxTime, TPZPoroElasticMF2d * &mymaterial,
                        TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);

void SolucaoExata1D(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);

bool const triang=false;
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.porolasticmf2d.data"));
#endif

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
	std::string logs("../logporoelastc2d.cfg");
	InitializePZLOG("../logporoelastc2d.cfg");
#endif
    
    int pu = 3;
    int pq = 3;
    int pp = 2;
	//primeira malha
	
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = GMesh(triang,1.,1.);
    ofstream arg1("gmesh_inicial.txt");
    gmesh->Print(arg1);
    
    // First computational mesh
    TPZCompMesh * cmesh1 = MalhaCompElast(gmesh,pu);
    ofstream arg2("cmesh1_inicial.txt");
    cmesh1->Print(arg2);
    
    // second computational mesh
	TPZCompMesh * cmesh2= CMeshFlux(gmesh, pq);
    ofstream arg3("cmesh2_inicial.txt");
    cmesh2->Print(arg3);
    
    
	// Third computational mesh
	TPZCompMesh * cmesh3 = CMeshPressure(gmesh, pp);
    ofstream arg4("cmesh3_inicial.txt");
    cmesh3->Print(arg4);
    
    
    // Cleaning reference of the geometric mesh to cmesh1
	gmesh->ResetReference();
	cmesh1->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,3);
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();
    ofstream arg5("cmesh1_final.txt");
    cmesh1->Print(arg5);
	
	
	// Cleaning reference to cmesh2
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,3);
	cmesh2->AdjustBoundaryElements();
	cmesh2->CleanUpUnconnectedNodes();
    ofstream arg6("cmesh2_final.txt");
    cmesh2->Print(arg6);
    
    // Cleaning reference to cmesh3
	gmesh->ResetReference();
	cmesh3->LoadReferences();
	TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh3,3);
	cmesh3->AdjustBoundaryElements();
	cmesh3->CleanUpUnconnectedNodes();
    ofstream arg7("cmesh3_final.txt");
    cmesh3->Print(arg7);
    
    
    //	Set initial conditions for pressure
    TPZAnalysis an3(cmesh3);
	SolveSist(an3, cmesh3);
	int nrs = an3.Solution().Rows();
	TPZFMatrix<REAL> solucao1(nrs,1,1000.);
	cmesh3->Solution() = solucao1;
    
    //malha multifisica
    TPZVec<TPZCompMesh *> meshvec(3);
	meshvec[0] = cmesh1;
	meshvec[1] = cmesh2;
    meshvec[2] = cmesh3;
    TPZPoroElasticMF2d * mymaterial;
    TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,mymaterial);
    ofstream arg8("mphysic.txt");
	mphysics->Print(arg8);
    
    REAL deltaT=1.; //second
    mymaterial->SetTimeStep(deltaT);
    REAL maxTime = 10.;
    
    //condicao inicial
   // mymaterial->SetLastState();
//    TPZAnalysis an(mphysics);
//	TPZFMatrix<REAL> Initialsolution = an.Solution();
    
    SolveSistTransient(deltaT, maxTime, mymaterial, meshvec, mphysics);
    
    return 0;
}

TPZGeoMesh *GMesh(bool triang_elements, REAL w, REAL L){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
    TPZVec <int> TopolTriang(3);
	TPZVec <int> TopolLine(2);
	
	//indice dos nos
	int id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*L;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = L - xi*L;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,w);//coord Y
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
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcBottom,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcRight,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcTop,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcLeft,*gmesh);
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
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcBottom,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcRight,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcTop,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcLeft,*gmesh);
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

TPZCompMesh*MalhaCompElast(TPZGeoMesh * gmesh,int pOrder)
{
    /// criar material
	int dim = 2;
	TPZVec<REAL> force(dim,0.);
    REAL E = 100.;
	REAL poisson = 0.35;
    int planestress = -1;
	TPZElasticityMaterial *material;
	material = new TPZElasticityMaterial(matId, E, poisson, force[0], force[1], planestress); 
	TPZMaterial * mat(material);
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
	///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1 = material->CreateBC(mat, bcBottom,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bcLeft,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bcTop,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bcRight,dirichlet, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);
	cmesh->InsertMaterialObject(BCond4);
	
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();	
	
	return cmesh;
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
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond1 = material->CreateBC(mat, bcBottom,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bcRight,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bcTop,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bcLeft,dirichlet, val1, val2);
    
	cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
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
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond1 = material->CreateBC(mat, bcBottom,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bcRight,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bcTop,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bcLeft,dirichlet, val1, val2);
    
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
        newnod.SetPressure(true);
    }
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
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
    
	return cmesh;
}

TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZPoroElasticMF2d * &mymaterial){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=1;
    int dim =2;
    
    //	criar material
    REAL Eyoung = 3.e4;
    REAL poisson = 0.2;
    REAL alpha = 1.0;
    REAL Se = 0.0;
    REAL rockrho = 2330.0; // SI system
    REAL gravity = 0.0;//-9.8; // SI system
    REAL fx=0.0;
    REAL fy=gravity*rockrho;
    REAL perm = 1.e-10;
    REAL visc = 1.e-3;
    int planestress = 0; // This is a Plain strain problem
    
    
    mymaterial = new TPZPoroElasticMF2d(MatId,dim);
    
    mymaterial->SetfPlaneProblem(planestress);
    mymaterial->SetParameters(perm, visc);
    mymaterial->SetParameters(Eyoung, poisson, fx, fy);
    mymaterial->SetBiotParameters(alpha, Se);  
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(SolucaoExata1D);
	mymaterial->SetForcingFunctionExact(solExata);
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
    
    REAL sig0 = -1000.;
    REAL ptop = 0.;
    val2(0,0)= 0.;
    val2(1,0)= sig0;
    val2(2,0)= ptop;
	TPZMaterial * BCond1 = mymaterial->CreateBC(mat, bcTop,neumdir, val1, val2);
    
    val2.Redim(3,1);
    REAL big = mymaterial->gBigNumber;
    val1(0,0) = big;
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat,bcRight, mixedneum, val1, val2);
    TPZMaterial * BCond4 = mymaterial->CreateBC(mat,bcLeft, mixedneum, val1, val2);
    
    val1.Redim(3,2);
    val2(2,0)= 1000.;
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat,bcBottom,dirichlet, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
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

void SolucaoExata1D(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
	//REAL x = ptx[0];
	REAL x = ptx[1];
	
	REAL pini = 1000.;
	REAL lamb = 8333.33;
	REAL mi = 12500.;
	REAL visc =0.001; 
	REAL perm =  1.e-10;
	REAL H=1.;
	REAL tp = 10.;
	int in;
	REAL pD, uD, sigD;
    REAL VDy;
	REAL PI = atan(1.)*4.;
	
	sol[0]=0.;
	sol[1]=0.;
	sol[2]=0.;
    flux(0,0)=0.;
    flux(1,0)=0.;
	
	REAL tD = (lamb+2.*mi)*perm*tp/(visc*H);
	REAL xD = fabs(1.-x)/H;
	for (in =0; in<1000; in++) {
		
		REAL M = PI*(2.*in+1.)/2.;
		pD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
		uD += (2./(M*M))*cos(M*xD)*exp(-1.*M*M*tD);
		sigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
        
        VDy+= 2.*cos(M*xD)*exp(-1.*M*M*tD);
	}
	
	sol[0] = pD*pini;
	sol[1] = (1.- xD - uD)*(-pini*H)/(lamb+2.*mi);
	sol[2] = (-1.+ sigD)*pini;
    
    flux(1,0) = (-1.)*pini*(perm/visc)*VDy;
}

void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
}

TPZAutoPointer <TPZMatrix<REAL> > MassMatrix(TPZPoroElasticMF2d * mymaterial, TPZCompMesh* mphysics){
    
    mymaterial->SetLastState();
    //TPZSkylineStructMatrix matsp(mphysics);
	TPZSpStructMatrix matsp(mphysics);
		
	std::set< int > materialid;
	int matid = mymaterial->MatId();
	materialid.insert(matid);
	matsp.SetMaterialIds (materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<REAL> Un;
    //TPZMatrix<REAL> *matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}

void StiffMatrixLoadVec(TPZPoroElasticMF2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec){
    
	mymaterial->SetCurrentState();
    //TPZFStructMatrix matsk(mphysics);
    TPZSkylineStructMatrix matsk(mphysics);
	an.SetStructuralMatrix(matsk); 
	TPZStepSolver<REAL> step; 
	step.SetDirect(ELDLt); 
	//step.SetDirect(ELU);
	an.SetSolver(step); 
	an.Run(); 
	
	matK1 = an.StructMatrix();
	fvec = an.Rhs();
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
{
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(9), vecnames(4);
	scalnames[0] = "DisplacementX";
	scalnames[1] = "DisplacementY";
    scalnames[2] = "SigmaX";
	scalnames[3] = "SigmaY";
	scalnames[4] = "PressureStress";
	scalnames[5] = "PorePressure";
    
	vecnames[0]  = "Displacement";
	vecnames[1]  = "DarcyVelocity";
	vecnames[2]  = "MinusKMuGradP";	
    
    scalnames[6] = "ExactPressure";
    scalnames[7] = "ExactDisplacementY";
    scalnames[8] = "ExactSigmaY";
    vecnames[3]  = "ExactDarcyVelocity";
	
	const int dim = 2;
	int div =0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void SolveSistTransient(REAL deltaT,REAL maxTime, TPZPoroElasticMF2d * &mymaterial ,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics){
	
    
    TPZAnalysis an(mphysics);
	TPZFMatrix<REAL> Initialsolution = an.Solution();
    
    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<REAL> > matM = MassMatrix(mymaterial, mphysics);
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
            std::stringstream sout;
        	matM->Print("matM = ", sout,EMathematicaInput);
        	LOGPZ_DEBUG(logdata,sout.str())
	}
#endif   
    
    //Criando matriz de rigidez (matK) e vetor de carga
	TPZFMatrix<REAL> matK;	
	TPZFMatrix<REAL> fvec; 
    StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fvec);
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		
        std::stringstream sout;
        matK.Print("matK = ", sout,EMathematicaInput);
        fvec.Print("fvec = ", sout,EMathematicaInput);		
        //Print the temporal solution
        Initialsolution.Print("Intial conditions = ", sout,EMathematicaInput);
        TPZFMatrix<REAL> Temp;
        TPZFMatrix<REAL> Temp2;
        matM->Multiply(Initialsolution,Temp);
        Temp.Print("Temp matM = ", sout,EMathematicaInput);	
        LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
    
    
	int nrows;
	nrows = matM->Rows();
	TPZFMatrix<REAL> TotalRhs(nrows,1,0.0);
	TPZFMatrix<REAL> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<REAL> Lastsolution = Initialsolution;
	
	std::string outputfile;
	outputfile = "TransientSolution";
	
	REAL TimeValue = 0.0;
	int cent = 0;
	TimeValue = cent*deltaT; 
	while (TimeValue < maxTime)
	{	
		// This time solution i for Transient Analytic Solution
		mymaterial->SetTimeValue(TimeValue);
		matM->Multiply(Lastsolution,TotalRhstemp);
        
#ifdef LOG4CXX
        if(logdata->isDebugEnabled())
        {
            std::stringstream sout;
            sout<< " tempo = " << cent;
            Lastsolution.Print("\nIntial conditions = ", sout,EMathematicaInput);
            TotalRhstemp.Print("Mat Mass x Last solution = ", sout,EMathematicaInput);	
            LOGPZ_DEBUG(logdata,sout.str())
        }
#endif

		TotalRhs = fvec + TotalRhstemp;
		an.Rhs() = TotalRhs;
		an.Solve(); 
		Lastsolution = an.Solution();
		
		std::stringstream outputfiletemp;
		outputfiletemp << outputfile << ".vtk";
		std::string plotfile = outputfiletemp.str();
		PosProcessMultphysics(meshvec,mphysics,an,plotfile);		
		
        cent++;
		TimeValue = cent*deltaT;
	}
}

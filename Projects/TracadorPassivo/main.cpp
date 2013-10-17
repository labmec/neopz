
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzinterpolationspace.h"
#include "TPZCompElDisc.h"
#include "pzpoisson3d.h"
#include "pzconvectionproblem.h"
#include "mixedpoisson.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pztracerflow.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>
using namespace std;

int matId =1;
int inflow = 0;
int outflow = 3;

int bc0=-1;
int bc1=-2;
int bc2=-3;
int bc3=-4;

int const dirichlet =0;
int const neumann = 1;


TPZGeoMesh *GMesh(bool triang_elements, REAL Lx, REAL Ly);

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, bool twomaterial);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder,bool twomaterial);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, bool twomaterial);
TPZCompMesh *CMeshSaturation(TPZGeoMesh * gmesh, int pOrder,TPZTracerFlow * &material);


void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcessFlux(TPZAnalysis &an, std::string plotfile);

void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref,TPZVec<REAL> pt, bool changeMatId, int newmatId, REAL &Area);
void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref);

void SolExata(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);

bool ftriang = false;

#include "pztransfer.h"
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    REAL Lx = 200.;
    REAL Ly = 200.;
    
    int pq = 1;
    int pp;
    if(ftriang==true){
        pp = pq-1;
    }else{
        pp = pq;
    }

    TPZVec<REAL> pt(3);
    pt[0] = Lx/2.;
    pt[1] = Ly/2.;
    REAL Area;

    
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    TPZGeoMesh *gmesh = GMesh(ftriang, Lx, Ly);
    ofstream arg("gmesh1.txt");
    gmesh->Print(arg);
    
    RefinamentoPadrao3x3(gmesh,1,pt, true, matId+1, Area);
    //UniformRefine(gmesh, 3);
    
    TPZCompMesh *cmesh1 = CMeshFlux(gmesh, pq,true);
    TPZCompMesh *cmesh2 = CMeshPressure(gmesh, pp,true);
    //TPZCompMesh *cmeshfluxcoarse = CMeshFlux(gmesh, pq,true);
    
    //---------------------------------------------------------------------
    // Cleaning reference of the geometric mesh to cmesh1
    gmesh->ResetReference();
    cmesh1->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,3,false);
    cmesh1->AdjustBoundaryElements();
    cmesh1->CleanUpUnconnectedNodes();
    ofstream arg1("cmeshflux.txt");
    cmesh1->Print(arg1);
    
    
    // Cleaning reference to cmesh2
    gmesh->ResetReference();
    cmesh2->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,3,true);
    cmesh2->AdjustBoundaryElements();
    cmesh2->CleanUpUnconnectedNodes();
    ofstream arg2("cmeshpressure.txt");
    cmesh2->Print(arg2);
    
    // Cleaning reference to cmesh3
//    gmesh->ResetReference();
//    cmeshfluxcoarse->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmeshfluxcoarse,2,false);
//    
//    cmeshfluxcoarse->AdjustBoundaryElements();
//    cmeshfluxcoarse->CleanUpUnconnectedNodes();
    //        ofstream arg7("cmesh3_final.txt");
    //        cmesh3->Print(arg7);

    //-----------------------------------------------------------------------
    
    
    ofstream arg4("gmesh2.txt");
    gmesh->Print(arg4);
    
    
    //malha multifisica
    TPZVec<TPZCompMesh *> meshvec(2);
    meshvec[0] = cmesh1;
    meshvec[1] = cmesh2;
    
    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec,true);
    ofstream arg3("cmeshmultiphysics.txt");
    mphysics->Print(arg3);
    
    TPZAnalysis an(mphysics);
    SolveSyst(an, mphysics);
    
    string plotfile("Solution_mphysics.vtk");
    PosProcessMultphysics(meshvec,  mphysics, an, plotfile);
    
    
    //cmesh1->LoadReferences();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZAnalysis anflux(cmesh1);
    cmesh1->Solution().Print("");
    anflux.Solution().Print("");
    string plotfile2("sol_flux.vtk");
    //cmesh1->LoadSolution(meshvec[0]->Solution());
    PosProcessFlux(anflux, plotfile2);
//    
//    TPZVec<REAL> erros;
//    ofstream saidaerro("Erro.txt");
//    anflux.SetExact(*SolExata);
//     anflux.PostProcessError(erros, saidaerro);
//
//    TPZTransfer<REAL> transfer;
//    cmesh1->BuildTransferMatrix(*cmeshfluxcoarse,transfer);
//    TPZFMatrix<STATE> coarsesol = cmeshfluxcoarse->Solution();
//    TPZAnalysis anfluxcoarse(cmeshfluxcoarse);
//    string plotfile3("sol_fluxcoarse.vtk");
//    PosProcessFlux(anfluxcoarse, plotfile3);
    
    
    
	return EXIT_SUCCESS;
}


void SolExata(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
    
    REAL x = ptx[0];
	REAL y = ptx[1];
    
        
    
    sol.Resize(1, 0.);
	
    flux(0,0)=0.;
    flux(1,0)=0.;
    
}


TPZGeoMesh *GMesh(bool triang_elements, REAL Lx, REAL Ly){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
    TPZVec <int> TopolTriang(3);
	TPZVec <int> TopolLine(2);
    TPZVec <int> TopolPoint(1);
	
	//indice dos nos
	int id = 0;
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
    
    cmesh->InsertMaterialObject(mat);
	cmesh->SetAllCreateFunctionsHDiv();
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

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, bool twomaterial)
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
    
    if(twomaterial==true){
        TPZMatPoisson3d *material2 = new TPZMatPoisson3d(matId+1,dim);
        TPZMaterial * mat2(material2);
        cmesh->InsertMaterialObject(mat2);
    }
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond1 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
	cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    
    if(twomaterial==true){
        
        //Inserir materiais
        std::set<int> MaterialIDs;
        MaterialIDs.insert(matId);
        MaterialIDs.insert(matId+1);
        MaterialIDs.insert(bc0);
        MaterialIDs.insert(bc1);
        MaterialIDs.insert(bc2);
        MaterialIDs.insert(bc3);
        
        cmesh->AutoBuild(MaterialIDs);
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
        return cmesh;
    }
    
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
    
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
    
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
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
            if(ftriang==true) celdisc->SetTotalOrderShape();
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


TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, bool twomaterial)
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
    
    if(twomaterial==true){
        TPZMatPoisson3d *material2 = new TPZMatPoisson3d(matId+1,dim);
        TPZMaterial * mat2(material2);
        cmesh->InsertMaterialObject(mat2);
    }
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
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
    
    
#ifdef DEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    if(twomaterial==true){
        
        //Inserir materiais
        std::set<int> MaterialIDs;
        MaterialIDs.insert(matId);
        MaterialIDs.insert(matId+1);
        MaterialIDs.insert(bc0);
        MaterialIDs.insert(bc1);
        MaterialIDs.insert(bc2);
        MaterialIDs.insert(bc3);
        
        cmesh->AutoBuild(MaterialIDs);
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
        return cmesh;
    }
    
    
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}



TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int dim =2;
    
    REAL coefk = 9.869233e-07;//[mˆ2]==[1000 mD], 1 mD = 9.869233e-13 mˆ2
    REAL coefvisc = 0.000894;//[Pa.s]==[0.894 cP],  1 cp = 10ˆ3 Pa.s
    REAL vazao = -1.;//[m/s]
    
    TPZMixedPoisson *material = new TPZMixedPoisson(matId,dim);
    
    material->SetPermeability(coefk);
    material->SetViscosity(coefvisc);
    
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    
    val2(0,0)=vazao;
    BCond3 = material->CreateBC(mat, bc3,neumann, val1, val2);
    
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
    
    return mphysics;
}

TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, bool twomaterial){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int dim =2;
    
    REAL coefk = 9.869233e-08;//[mˆ2]==[100 mD], 1 mD = 9.869233e-13 mˆ2
    REAL coefvisc = 0.000894;//[Pa.s]==[0.894 cP],  1 cp = 10ˆ3 Pa.s
    REAL vazao = -1.;//[m/s]
    
    TPZMixedPoisson *material1 = new TPZMixedPoisson(matId,dim);
    
    material1->SetPermeability(coefk);
    material1->SetViscosity(coefvisc);
    
    TPZMaterial *mat1(material1);
    mphysics->InsertMaterialObject(mat1);
    
    //criar segundo material
    if(twomaterial==true){
        TPZMixedPoisson *material2 = new TPZMixedPoisson(matId+1,dim);
        material2->SetPermeability(coefk/10.);
        material2->SetViscosity(coefvisc);
        TPZMaterial * mat2(material2);
        mphysics->InsertMaterialObject(mat2);
    }

    
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    BCond0 = material1->CreateBC(mat1, bc0,neumann, val1, val2);
    BCond1 = material1->CreateBC(mat1, bc1,dirichlet, val1, val2);
    BCond2 = material1->CreateBC(mat1, bc2,neumann, val1, val2);
    
    val2(0,0)=vazao;
    BCond3 = material1->CreateBC(mat1, bc3,neumann, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    
    
    //Inserir materiais
     if(twomaterial==true){
        std::set<int> MaterialIDs;
        MaterialIDs.insert(matId);
        MaterialIDs.insert(matId+1);
        MaterialIDs.insert(bc0);
        MaterialIDs.insert(bc1);
        MaterialIDs.insert(bc2);
        MaterialIDs.insert(bc3);
        
        mphysics->AutoBuild(MaterialIDs);
        mphysics->AdjustBoundaryElements();
        mphysics->CleanUpUnconnectedNodes();
         
         // Creating multiphysic elements into mphysics computational mesh
         TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
         TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
         TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
         
         return mphysics;
     }
    
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

TPZCompMesh *CMeshSaturation(TPZGeoMesh * gmesh, int pOrder,TPZTracerFlow * &material)
{
	/// criar materiais
	int dim = 2;
	
	material = new TPZMatConvectionProblem(matId,dim);
	TPZMaterial * mat(material);
	
	TPZVec<REAL> convdir(dim,0.);
    convdir[0]=1.;
	REAL flux = 0.;
    REAL rho = 1.;
	
	material->SetParameters(rho,convdir);
	material->SetInternalFlux(flux);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->InsertMaterialObject(mat);
    
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    REAL uD =2.;
    val2(0,0) = uD;
	TPZMaterial * BCond3 = material->CreateBC(mat, bc3,inflow, val1, val2);
    val2(0,0) = 0.;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,outflow, val1, val2);
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond3);
    
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
    //TPZCompElDisc::SetTotalOrderShape(cmesh);
	
	return cmesh;
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
	TPZManVector<std::string,10> scalnames(1), vecnames(2);
	vecnames[0]  = "Flux";
    vecnames[1]  = "GradFluxX";
    scalnames[0] = "Pressure";
    
	const int dim = 2;
	int div =0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
    
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


void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref){
    
    {//NAO APAGAR!!!
        char buf[] =
        "16 10 "
        "-50 Qua000022224 "
        "-1 -1 0 "
        "1 -1 0 "
        "1 1 0  "
        "-1 1 0 "
        "-0.33333 -1 0  "
        "+0.33333 -1 0  "
        "+1 -0.33333 0  "
        "+1 +0.33333 0  "
        "+0.33333 +1 0  "
        "-0.33333 +1 0  "
        "-1 +0.33333 0  "
        "-1 -0.33333 0  "
        "-0.33333 -0.33333 0    "
        "+0.33333 -0.33333 0    "
        "+0.33333 +0.33333 0    "
        "-0.33333 +0.33333 0    "
        "3 4 0  1  2  3 "
        "3 4 0 4 12 11  "
        "3 4 4 5 13 12  "
        "3 4 5 1 6 13   "
        "3 4 11 12 15 10    "
        "3 4 12 13 14 15    "
        "3 4 13 6 7 14  "
        "3 4 10 15 9 3  "
        "3 4 15 14 8 9  "
        "3 4 14 7 2 8   ";
        std::istringstream str(buf);
        TPZAutoPointer<TPZRefPattern> refp = new TPZRefPattern(str);
        refp->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp);
        if(!refp)
        {
            DebugStop();
        }
    }
    
    TPZAutoPointer<TPZRefPattern> refp2D = gRefDBase.FindRefPattern("Qua000022224");
    //TPZAutoPointer<TPZRefPattern> refp1D = gRefDBase.FindRefPattern("Lin000022224");
    
    if(!refp2D) DebugStop();
    
    TPZGeoEl * gel = NULL;
    for(int r = 0; r < nref; r++)
    {
        int nels = gmesh->NElements();
        for(int iel = 0; iel < nels; iel++)
        {
            gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            if(gel->Dimension()==2)
            {
                gel->SetRefPattern(refp2D);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
            if(gel->Dimension()==1)
            {
                TPZAutoPointer<TPZRefPattern> refp1D = TPZRefPatternTools::PerfectMatchRefPattern(gel);
                if(refp1D)
                {
                    gel->SetRefPattern(refp1D);
                    TPZVec<TPZGeoEl*> sons;
                    gel->Divide(sons);
                }
                else{
                    DebugStop();//nao conseguiu refinar elementos 1D baseados no refp do pai
                }
            }
            
        }
    }
    
    std::ofstream malhaOut("malhaOut.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
}

void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref,TPZVec<REAL> pt, bool changeMatId, int newmatId, REAL &Area){
    
    {//NAO APAGAR!!!
        char buf[] =
        "16 10 "
        "-50 Qua000022224 "
        "-1 -1 0 "
        "1 -1 0 "
        "1 1 0  "
        "-1 1 0 "
        "-0.33333 -1 0  "
        "+0.33333 -1 0  "
        "+1 -0.33333 0  "
        "+1 +0.33333 0  "
        "+0.33333 +1 0  "
        "-0.33333 +1 0  "
        "-1 +0.33333 0  "
        "-1 -0.33333 0  "
        "-0.33333 -0.33333 0    "
        "+0.33333 -0.33333 0    "
        "+0.33333 +0.33333 0    "
        "-0.33333 +0.33333 0    "
        "3 4 0  1  2  3 "
        "3 4 0 4 12 11  "
        "3 4 4 5 13 12  "
        "3 4 5 1 6 13   "
        "3 4 11 12 15 10    "
        "3 4 12 13 14 15    "
        "3 4 13 6 7 14  "
        "3 4 10 15 9 3  "
        "3 4 15 14 8 9  "
        "3 4 14 7 2 8   ";
        std::istringstream str(buf);
        TPZAutoPointer<TPZRefPattern> refp = new TPZRefPattern(str);
        refp->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp);
        if(!refp)
        {
            DebugStop();
        }
    }

    
    TPZAutoPointer<TPZRefPattern> refpOutroLugar = gRefDBase.FindRefPattern("Qua000022224");
    if(!refpOutroLugar) DebugStop();
    
    int iniEl = 0;
    TPZVec<REAL> qsi(2,0.);
    
    TPZGeoEl * gel = NULL;
    for(int r = 0; r < nref; r++)
    {
        gel = gmesh->FindElement(pt, qsi, iniEl,2);
        if(!gel) DebugStop();
        if(gel->Dimension()==2)
        {
            gel->SetRefPattern(refpOutroLugar);
            TPZVec<TPZGeoEl*> sons;
            gel->Divide(sons);
            for(int edg = gel->NNodes(); edg < gel->NSides(); edg++)
            {
                TPZGeoElSide gelS(gel,edg);
                if(gelS.Dimension() > 1)
                {
                    break;
                }
                TPZGeoElSide neighS(gelS.Neighbour());
                while(neighS != gelS)
                {
                    if(neighS.Element()->Dimension() == 1)
                    {
                        TPZAutoPointer<TPZRefPattern> refp = gel->GetRefPattern();
                        TPZAutoPointer<TPZRefPattern> refp2 = refp->SideRefPattern(gelS.Side());
                        if(refp2)
                        {
                            TPZGeoEl * gel1D = neighS.Element();
                            gel1D->SetRefPattern(refp2);
                            TPZVec<TPZGeoEl*> sons2;
                            gel1D->Divide(sons2);
                        }
                    }
                    neighS = neighS.Neighbour();
                }
            }
        }
    }
    
    if(changeMatId==true)
    {
        gel = gmesh->FindElement(pt, qsi, iniEl,2);
        if(!gel) DebugStop();
        gel->SetMaterialId(newmatId);
        Area = gel->Volume();
    }
    
    std::ofstream malhaOut("malhaOut2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
}


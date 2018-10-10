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
#include "pzgeopoint.h"
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
#include "pzaxestools.h"
#include "TPZCopySolve.h"
#include "pzstrmatrix.h"

#include <iostream>
#include <math.h>
using namespace std;

int MatId = 1;

#define KConstant

int bcdirichlet = 0;
int bcneumann = 1;

int BC0 = -1;
int BC1 = -2;
int BC2 = -3;
int BC3 = -4;
int BC4 = -5;
int BC5 = -6;

//Meio altamente heterogeneo
TPZGeoMesh *GMesh(REAL Lx,REAL Ly);

TPZCompMesh *CMeshFluxo(int pOrder, TPZGeoMesh * gmesh);
TPZCompMesh *CMeshPressure(int pOrder, TPZGeoMesh * gmesh);
TPZCompMesh *CMeshMixed(TPZVec<TPZCompMesh *> meshvec, TPZGeoMesh * gmesh);
void PermeabilityTensor(const TPZVec<REAL> &pt,TPZVec<STATE> &kabs,TPZFMatrix<STATE> &tensorK);
void ForcingHighHeter(const TPZVec<REAL> &pt,TPZVec<STATE> &disp);
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads, bool direct);
TPZFMatrix<STATE> * ComputeInverse(TPZCompMesh * mphysics);

void RefinamentoUnif(TPZGeoMesh * gmesh,int nDiv);
void ResolverSistema(TPZAnalysis &an, TPZCompMesh * fCmesh, int numthreads,bool direct);
void PosProcessMultph(TPZVec<TPZCompMesh * >meshvec,TPZCompMesh* mphysics,TPZAnalysis &an, std::string plotfile);
void PosProcessFluxo(TPZAnalysis &an,std::string plotfile);

// Lado direito da equacao - meio altamente heterogeneo.
void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &dips);

// Condicao de Contorno - meio altamente heterogeneo.
void NeummanEntPres(const TPZVec<REAL> &pt,TPZVec<STATE> &disp);
void NeumannSaiPres(const TPZVec<REAL> &pt,TPZVec<STATE> &disp);
void NeummanEntFlux(const TPZVec<REAL> &pt,TPZVec<STATE> &disp);
void NeumannSaiFlux(const TPZVec<REAL> &pt,TPZVec<STATE> &disp);

// Analisis Solucao Exata
TPZCompMesh *CmeshFluxTeo = NULL;
TPZCompMesh *CmeshPresTeo = NULL;


//erros
void SolExataFluxo(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void SolExataPress(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void ComputeFluxError(TPZCompMesh *CMixedCoarse,int64_t NivelRef, std::ostream &out);
void ComputePressureError(TPZCompMesh *CmeshPres,int64_t NivelRef, std::ostream &out);
void ComputeErrorNomrs(TPZCompMesh *CMixedCoarse,int64_t NivelRef, std::ostream &out, bool IsFlux);

void ComputeFatherqsi(TPZGeoEl *Son, TPZGeoEl *Father,TPZManVector<REAL,3> &Sonqsi, TPZManVector<REAL,3> &Fatherqsi, TPZTransform<> &Transformation);


//void NeumannBound1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
//void DirichletXIgualDeis(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//REAL const pi = 4.*atan(1.);

bool fTriang = false;
bool IsStab = true;
bool IsContinuou = false;
bool Useh2 = false;
REAL Delta1 = 0.5;
REAL Delta2 = 0.5;
REAL FontValue = 1.0e5;
bool IsFullHdiv=false;
bool IsHomogeneo=false;
bool IsHHeterogeneo = true;
int64_t NivelMaxi = 2;

REAL Lx = 16.0;
REAL Ly = 16.0;

int main(int argc, char *argv[]){

    InitializePZLOG();
    std::cout.precision(5);
    
    ofstream saidaerro( "erros-hdiv-estab.txt");

    for(int p =1; p<2;p++){
        int ndiv;
        int pq = p+1;
        int pp = p;
        TPZGeoMesh * gmesh = NULL;
        
        for(ndiv = NivelMaxi;ndiv >= NivelMaxi-1;ndiv--){
            
            gmesh = GMesh(Lx, Ly);
            
            ofstream ArgGeo("GeoMesh.txt");
            gmesh->Print(ArgGeo);
            RefinamentoUnif(gmesh, ndiv);
            
            TPZCompMesh * cmesh1 = CMeshFluxo(pq, gmesh);
            TPZCompMesh * cmesh2 = CMeshPressure(pp, gmesh);
            
            //malha multifisica
            TPZVec<TPZCompMesh*> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            
            TPZCompMesh * mphysics = CMeshMixed(meshvec, gmesh);
            {
                ofstream arg4("gmeshMulti.txt");
                mphysics->Print(arg4);
            }
            std::cout << "Number of equations " << mphysics->NEquations() << std::endl;
            int numthreads = 1;
            std::cout << "Number of threads " << numthreads << std::endl;
            
            int64_t neq = mphysics->NEquations();
            TPZFMatrix<STATE> rhs(neq,1);
            TPZBandStructMatrix full(mphysics);
            
//            TPZMatrix<STATE> *Stiff = full.CreateAssemble(rhs, 0);
//            ofstream ArgStiff("Stiff.nb");
//            Stiff->Print("Stiff=",ArgStiff,EMathematicaInput);
            
            bool opti = true;
            TPZAnalysis  * analysis = new TPZAnalysis(mphysics,opti);
            ResolverSistema(*analysis, mphysics,numthreads, true);

            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            
            if( ndiv == NivelMaxi ){
                CmeshFluxTeo = mphysics;
                CmeshPresTeo = mphysics;
            }
            else{
//            ofstream SolPressFile("SolPress.nb");
//            SolPress.Print("Press=",SolPressFile,EMathematicaInput);
            
                //saidaerro<<"\nErro da simulacao multifisica do fluxo (q)" <<endl;
                ComputeErrorNomrs(mphysics,ndiv,saidaerro,true);
                
                //saidaerro<<"\nErro da simulacao multifisica da pressao (p)" <<endl;
                ComputeErrorNomrs(mphysics,ndiv,saidaerro,false);
                
            //Plot da solucao aproximada
                string plotfile("Solution_mphysics.vtk");
                char buf[256] ;
                sprintf(buf,"ProblemaJuanGC_porder%d_h%d.vtk",p,ndiv);
                PosProcessMultph(meshvec,  mphysics, *analysis, buf);
            }
        }
    }
    return 0;
}

TPZGeoMesh * GMesh(REAL Lx,REAL Ly){
    
    int Qnodes = 6;
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec<int64_t> TopolQuad(4);
    TPZVec<int64_t> TopolTriangle(4);
    TPZVec<int64_t> TopolLine(2);
    TPZVec<int64_t> TopolPoint(1);
    
    //indice dos nos
    int64_t id = 0 ;
    REAL valx, valy;
    REAL dx = Lx;
    
    //Node 1
    for(int ix=0;ix<2;ix++){
        valx = ix*dx;
        valy = 0.;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0,valx);
        Node[id].SetCoord(1,valy);
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    for(int ix=0;ix<2;ix++){
        valx = Lx-ix*dx;
        valy = Ly;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0,valx);
        Node[id].SetCoord(1,valy);
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    valx = Lx;
    valy = (1.0/8.0)*Ly;
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0,valx);
    Node[id].SetCoord(1,valy);
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    valx = (1.0/8.0)*Lx;
    valy = Ly;
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0,valx);
    Node[id].SetCoord(1,valy);
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    // Nodos de elementos
    id = 0;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 4;
    TopolQuad[2] = 5;
    TopolQuad[3] = 3;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(id,TopolQuad,MatId,*gmesh);
    id++;

    TopolTriangle[0] = 0;
    TopolTriangle[1] = 1;
    TopolTriangle[2] = 4;
    new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(id,TopolTriangle,MatId,*gmesh);
    id++;

    TopolTriangle[0] = 4;
    TopolTriangle[1] = 2;
    TopolTriangle[2] = 5;
    new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(id,TopolTriangle,MatId,*gmesh);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(id,TopolLine,BC0,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 4;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(id,TopolLine,BC1,*gmesh);
    id++;
    
    TopolLine[0] = 4;
    TopolLine[1] = 2;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(id,TopolLine,BC2,*gmesh);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 5;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(id,TopolLine,BC3,*gmesh);
    id++;
    
    TopolLine[0] = 5;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(id,TopolLine,BC4,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(id,TopolLine,BC5,*gmesh);
    
    gmesh->BuildConnectivity();
    return gmesh;
}

void RefinamentoUnif(TPZGeoMesh *gmesh,int nDiv){
    for(int D = 0; D<nDiv; D++){
        int nels = gmesh->NElements();
        for(int elem = 0;elem<nels;elem++){
            TPZVec<TPZGeoEl *> filhos;
            TPZGeoEl * gel = gmesh ->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
}

TPZCompMesh *CMeshFluxo(int pOrder, TPZGeoMesh * gmesh){

    // criar materiais
    int dim=2;
    TPZMatPoisson3d * material;
    material = new TPZMatPoisson3d(MatId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    if(IsFullHdiv){
        cmesh->SetAllCreateFunctionsHDivFull();
    }
    else{
        cmesh->SetAllCreateFunctionsHDiv();
    }
    
    cmesh->InsertMaterialObject(mat);
    
    // Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);

    TPZMaterial * BCond0 = material->CreateBC(mat, BC0, bcdirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, BC1, bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, BC2, bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, BC3, bcdirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, BC4, bcdirichlet, val1, val2);
    TPZMaterial * BCond5 = material->CreateBC(mat, BC5, bcdirichlet, val1, val2);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEnt;
    bcmatNeumannEnt = new TPZDummyFunction<STATE>(NeummanEntFlux, 5);
    BCond1->SetForcingFunction(bcmatNeumannEnt);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeummanSai;
    bcmatNeummanSai = new TPZDummyFunction<STATE>(NeumannSaiFlux, 5);
    BCond4->SetForcingFunction(bcmatNeummanSai);
 
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);

    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    {
        std::ofstream out("MeshFlux.txt");
        cmesh->Print(out);
    }
#endif
    return cmesh;
}

TPZCompMesh * CMeshPressure(int pOrder,TPZGeoMesh * gmesh){

    //Criar materiais
    int dim=2;
    TPZMatPoisson3d * material;
    material = new TPZMatPoisson3d(MatId,dim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    
    // Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);
    
    TPZMaterial * BCond0 = material->CreateBC(mat, BC0, bcneumann, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, BC1, bcneumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, BC2, bcneumann, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, BC3, bcneumann, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, BC4, bcneumann, val1, val2);
    TPZMaterial * BCond5 = material->CreateBC(mat, BC5, bcneumann, val1, val2);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEnt;
    bcmatNeumannEnt = new TPZDummyFunction<STATE>(NeummanEntPres, 5);
    BCond1->SetForcingFunction(bcmatNeumannEnt);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeummanSai;
    bcmatNeummanSai = new TPZDummyFunction<STATE>(NeumannSaiPres, 5);
    BCond4->SetForcingFunction(bcmatNeummanSai);

    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);

    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    // Ajuste de Estrutura de dados computacional
    if(IsContinuou==false)
    {
        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->AutoBuild();
        
        int ncon = cmesh->NConnects();
        for(int i = 0; i<ncon; i++)
        {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
        
        int nel = cmesh->NElements();
        
        for(int i=0;i<nel;i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc*>(cel);
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
    
    return cmesh;
}

TPZCompMesh *CMeshMixed(TPZVec<TPZCompMesh*> meshvec,TPZGeoMesh * gmesh){

    //Creating computational mesh for multiphysics elements
    gmesh->ResetReference();
    TPZCompMesh * mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim=2;
    TPZMixedPoisson * material = new TPZMixedPoisson(MatId,dim);
    //Incluindo os dados do problema;
    REAL coefk = 1.;
    material->SetPermeability(coefk);
    REAL coefvisc = 1.;
    material->SetViscosity(coefvisc);
    
    //permeabilidade
    TPZAutoPointer<TPZFunction<STATE> > tensorK;
    tensorK = new TPZDummyFunction<STATE>(PermeabilityTensor, 5);
    material->SetPermeabilityFunction(tensorK);
    
    if(IsStab == true){
        material->SetStabilizedMethod();
        material->SetStabilizationCoeficients(Delta1, Delta2);
    }
    
    if(IsStab== true && Useh2 == true) material->SetHdois();
    ///// Solucao Exata ???
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForcingHighHeter, 5);
    
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    //Criadno condicoes de contorno;
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * BCond0 = material->CreateBC(mat, BC0, bcneumann, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, BC1, bcneumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, BC2, bcneumann, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, BC3, bcneumann, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, BC4, bcneumann, val1, val2);
    TPZMaterial * BCond5 = material->CreateBC(mat, BC5, bcneumann, val1, val2);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEnt;
    bcmatNeumannEnt = new TPZDummyFunction<STATE>(NeummanEntFlux, 5);
    BCond1->SetForcingFunction(bcmatNeumannEnt);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeummanSai;
    bcmatNeummanSai = new TPZDummyFunction<STATE>(NeumannSaiFlux, 5);
    BCond4->SetForcingFunction(bcmatNeummanSai);

    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->SetDimModel(dim);
    
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

void PermeabilityTensor(const TPZVec<REAL> &pt,TPZVec<STATE> &kabs,TPZFMatrix<STATE> &tensorK){
    
#ifndef KConstant
    
    REAL x = pt[0];
    REAL y = pt[1];
    tensorK.Resize(4,2);
    kabs.Resize(1, 0);
    bool Omega1,Omega2,Omega3,Omega4;
 //   Omega1 = (3<x) && (x<7);
    
    Omega1 = ((2. <= x) && (x < 4.) && (0. <= y) && (y< 2. )) || ((2. <= x) && (x < 4.) && (3. <= y ) && (y< 7.)) || ((4. <= x ) && (x < 9.) && (4. <= y) && (y< 7.)) || ((6. <= x) && (x < 9.) && (7. <= y) && (y < 11.)) || ((9. <= x )&& (x < 11.) && (9. <= y )&& (y < 11.))||  ((11. <= x) && (x <= Lx) && (11. <= y )&& (y < 14.)) ;
    Omega2 =  (2. <= x) && (x < 4.)  && (2. <= y ) && (y< 3.);
    Omega3 = ((0. <= x) && (x < 2.)  && (y >= 0.)) || ((2. <= x )&& (x < 6.)  && (y >= 7.)) || ((6. <= x )&& (x < 11.)  && (y >= 11.)) || ((11. <= x)&&(x <= Lx)  && (y >= 14.)) ;
    Omega4 = ((4. <= x)&& (x < 9.)  && (y <= 4.)) || ((9. <= x ) && (x < 11.)  && (y <= 9.)) || ((11. <= x ) && (x <= Lx)  && (y <= 11.)) ;
    
    if(Omega1){
        //K
        tensorK(0,0) = 10.;     tensorK(0,1) = 0.0;
        tensorK(1,0) = 0.0;     tensorK(1,1) = 10.;
        
        //Kinv
        tensorK(2,0) = 1.0/0.1;     tensorK(2,1) = 0.0;
        tensorK(3,0) = 0.0;     tensorK(3,1) = 1.0/0.1;
    }
    else if(Omega2){
        //K
        tensorK(0,0) = 100.0;     tensorK(0,1) = 0.0;
        tensorK(1,0) = 0.0;     tensorK(1,1) = 0.1;
        
        //Kinv
        tensorK(2,0) = 1.0/100.0;     tensorK(2,1) = 0.0;
        tensorK(3,0) = 0.0;     tensorK(3,1) = 1.0/0.1;
    }
    else if (Omega3){
        //K
        tensorK(0,0) = 0.1;     tensorK(0,1) = 0.0;
        tensorK(1,0) = 0.0;     tensorK(1,1) = 0.1;
        
        //Kinv
        tensorK(2,0) = 1.0/.1;     tensorK(2,1) = 0.0;
        tensorK(3,0) = 0.0;     tensorK(3,1) = 1.0/.1;
    }else if (Omega4){
        //K
        tensorK(0,0) = 0.1;     tensorK(0,1) = 0.0;
        tensorK(1,0) = 0.0;     tensorK(1,1) = 0.1;
        
        //Kinv
        tensorK(2,0) = 1.0/0.1;     tensorK(2,1) = 0.0;
        tensorK(3,0) = 0.0;     tensorK(3,1) = 1.0/0.1;
    }
    else{
        DebugStop();
    }
    
#else
    
    tensorK.Resize(4,2);
    tensorK.Zero();
    
    //K
    tensorK(0,0) = 1.0;
    tensorK(1,1) = 1.0;
    
    //Kinv
    tensorK(2,0) = 1.0;
    tensorK(3,1) = 1.0;
    
#endif
    
}

void ForcingHighHeter(const TPZVec<REAL> &pt,TPZVec<STATE> &disp){
    disp[0]=1.0;
}

void NeummanEntPres(const TPZVec<REAL> &pt,TPZVec<STATE> &disp){
    TPZVec<STATE> kabs(2,0.);
    TPZFMatrix<STATE> tensorK(2,2,0.);
    PermeabilityTensor(pt,kabs,tensorK);
    disp[0]=+FontValue/tensorK(1,1);
}

void NeumannSaiPres(const TPZVec<REAL> &pt,TPZVec<STATE> &disp){
    TPZVec<STATE> kabs;
    TPZFMatrix<STATE> tensorK(2,2,0.);
    PermeabilityTensor(pt,kabs,tensorK);
    disp[0]=-FontValue/tensorK(0,0);
}

void NeummanEntFlux(const TPZVec<REAL> &pt,TPZVec<STATE> &disp){
    disp[0]=-FontValue;
}

void NeumannSaiFlux(const TPZVec<REAL> &pt,TPZVec<STATE> &disp){
    disp[0]=+FontValue;
}


void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads, bool direct)
{
    
    if (direct)
    {
        TPZSkylineStructMatrix full(fCmesh); //caso simetrico
        full.SetNumThreads(numthreads);
        an.SetStructuralMatrix(full);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt); //caso simetrico
        //  step.SetDirect(ELU);
        an.SetSolver(step);
        an.Run();
    }
    else{
        
        TPZSkylineStructMatrix full(fCmesh); //caso simetrico
        full.SetNumThreads(numthreads);
        an.SetStructuralMatrix(full);
        TPZStepSolver<STATE> step;
        TPZAutoPointer<TPZMatrix<STATE> > Inverse = ComputeInverse(fCmesh);
        TPZStepSolver<STATE> precond(Inverse);
        precond.SetMultiply();
        //        step.SetGMRES(10, 40, precond, 1.0e-10, 0);
        step.SetCG(10, precond, 1.0e-10, 0);
        an.SetSolver(step);
        an.Run();
    }
    
    //Saida de Dados: solucao e  grafico no VT
    //	ofstream file("Solutout");
    //	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

TPZFMatrix<STATE> * ComputeInverse(TPZCompMesh * mphysics)
{
    int neq = mphysics->NEquations();
    TPZFMatrix<STATE> * PreInverse =  new TPZFMatrix<STATE> (neq,neq,0.0);
    TPZSkylineStructMatrix skyl(mphysics);
    std::set<int> matids; // to be computed
    matids.insert(MatId);
    matids.insert(BC0);
    matids.insert(BC1);
    matids.insert(BC2);
    matids.insert(BC3);
    skyl.SetMaterialIds(matids);
    TPZFMatrix<STATE> rhsfrac;
    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
    TPZAutoPointer<TPZMatrix<STATE> > matfrac = skyl.CreateAssemble(rhsfrac, gui);
    TPZFMatrix<STATE> oldmat = *matfrac.operator->();
    matfrac->Inverse( * PreInverse,ELDLt);
    return PreInverse;
    
}

void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &dips){
    dips[0]=0.;
}

void SolExataFluxo(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
}

void SolExataPress(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
}


void ComputeErrorNomrs(TPZCompMesh *CMixedCoarse,int64_t NivelRef, std::ostream &out, bool IsFlux){
    
    
    TPZCompMesh * CMixedFine = CmeshFluxTeo;
    CMixedFine->Reference()->ResetReference();
    CMixedFine->LoadReferences();
//    TPZGeoMesh * GmeshFine = CMixedFine->Reference();
    
    bool noisymode = false;
    
    TPZMaterial *material = CMixedFine->FindMaterial(1); // Here Juan fixed for MatId  == 1
    TPZMixedPoisson *matp = dynamic_cast<TPZMixedPoisson *>(material);
    
    if (!matp) {
        DebugStop();
    }

    CMixedCoarse->Reference()->ResetReference();
    CMixedCoarse->LoadReferences();
    TPZGeoMesh * GmeshCoarse = CMixedCoarse->Reference();
    
    int nFcel = CMixedFine->NElements();
//    int nCcel = CMixedCoarse->NElements();
//    int nFgel = GmeshFine->NElements();
//    int nCgel = GmeshCoarse->NElements();

    TPZManVector<REAL,2> errors(3,0.);
    int dimmesh = CMixedFine->Dimension();
//    int nel = CMixedFine->NElements();
//    int iel;
//    int64_t nivel = NivelMaxi-1;
    
    for (int  icel = 0; icel < nFcel; icel++) {
        TPZCompEl *celF = CMixedFine->Element(icel);
//        int celFdim= celF->Dimension();
        
        // Conditions to be avoided
        if (!celF) continue;
        if (dimmesh != celF->Dimension()) continue;
        
        TPZGeoEl * gelF = celF->Reference();
        if(!gelF) continue;
        
        TPZGeoEl * gelC = GmeshCoarse->Element(gelF->Father()->Index());
        if(!gelC) continue;
        
        TPZCompEl *celC = gelC->Reference();
        if (!celC) continue;
        
        if (CMixedCoarse->Element(celC->Index()) != celC) continue;
        if (dimmesh != celC->Dimension()) continue;
        
//        celC->Print();
//        CMixedCoarse->Element(celC->Index())->Print();
        
        
        TPZAutoPointer<TPZIntPoints> intrule = gelF->CreateSideIntegrationRule(gelF->NSides()-1, 2);
        TPZManVector<int,3> prevorder(dimmesh), maxorder(dimmesh,intrule->GetMaxOrder());
        intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);

        const int npoints = intrule->NPoints();
        TPZManVector<REAL,3> qsiC(dimmesh), qsiF(dimmesh), xVecC(3), xVecF(3);
        
        for (int ipoint = 0; ipoint < npoints; ipoint++) {
            REAL weight;
            intrule->Point(ipoint, qsiF, weight);
            
            TPZTransform<> Transformation;
            ComputeFatherqsi(gelF, gelF->Father(), qsiF, qsiC, Transformation);
            
            
            if (IsFlux) {
                
                int qfluxindex = matp->VariableIndex("Flux");
                int Gqxfluxindex = matp->VariableIndex("GradFluxX");
                int Gqyfluxindex = matp->VariableIndex("GradFluxY");
                int Divqfluxindex = matp->VariableIndex("DivFlux");
                
                TPZVec<STATE> qfluxF, dqfluxdxF, dqfluxdyF, divqfluxF;
                TPZVec<STATE> qfluxC, dqfluxdxC, dqfluxdyC, divqfluxC;
                
                celF->Solution(qsiF, qfluxindex, qfluxF);
                celC->Solution(qsiC, qfluxindex, qfluxC);
                
                celF->Solution(qsiF, Divqfluxindex, divqfluxF);
                celC->Solution(qsiC, Divqfluxindex, divqfluxC);
                
                celF->Solution(qsiF, Gqxfluxindex, dqfluxdxF);
                celC->Solution(qsiC, Gqxfluxindex, dqfluxdxC);
                
                celF->Solution(qsiF, Gqyfluxindex, dqfluxdyF);
                celC->Solution(qsiC, Gqyfluxindex, dqfluxdyC);
                

                if (noisymode) {
                    
                    gelF->X(qsiF, xVecF);
                    gelC->X(qsiC, xVecC);
                    
                    std::cout << "qsiF " << qsiF[0] << "  " << qsiF[1] << std::endl;
                    std::cout << "qsiC " << qsiC[0] << "  " << qsiC[1] << std::endl;
                    std::cout << "xVecF " << xVecF[0] << "  " << xVecF[1] << "  " << xVecF[2] << std::endl;
                    std::cout << "xVecC " << xVecC[0] << "  " << xVecC[1] << "  " << xVecC[2] << std::endl;
                    
                    
                    
                    std::cout << "qfluxF " << qfluxF << std::endl;
                    std::cout << "qfluxC " << qfluxC << std::endl;
                    
                    std::cout << "divqfluxF " << divqfluxF << std::endl;
                    std::cout << "divqfluxC " << divqfluxC << std::endl;
                    
                    std::cout << "dqfluxdxF " << dqfluxdxF << std::endl;
                    std::cout << "dqfluxdxC " << dqfluxdxC << std::endl;
                    
                    std::cout << "dqfluxdyF " << dqfluxdyF << std::endl;
                    std::cout << "dqfluxdyC " << dqfluxdyC << std::endl;
                }
                
                TPZManVector<REAL,2> diff(2,0.);
                diff[0] = qfluxC[0] - qfluxF[0] ;
                diff[1] = qfluxC[1] - qfluxF[1] ;
                
                //erro L2 do fluxo
                errors[0] += weight*(diff[0]*diff[0] + diff[1]*diff[1]);
                
                //erro L2 do divergente do fluxo
                REAL diffDiv = abs(divqfluxC[0] - divqfluxF[0]);
                errors[1] += weight*diffDiv*diffDiv;
                
            }else{

                int Pressureindex = matp->VariableIndex("Pressure");
                int GradPressureindex = matp->VariableIndex("GradPressure");
                TPZVec<STATE> PressureF, GradPressureF;
                TPZVec<STATE> PressureC, GradPressureC;
                
                celF->Solution(qsiF, Pressureindex, PressureF);
                celC->Solution(qsiC, Pressureindex, PressureC);
                
                celF->Solution(qsiF, GradPressureindex, GradPressureF);
                celC->Solution(qsiC, GradPressureindex, GradPressureC);
                
                
                if (noisymode) {
                    
                    gelF->X(qsiF, xVecF);
                    gelC->X(qsiC, xVecC);
                    
                    std::cout << "qsiF " << qsiF[0] << "  " << qsiF[1] << std::endl;
                    std::cout << "qsiC " << qsiC[0] << "  " << qsiC[1] << std::endl;
                    std::cout << "xVecF " << xVecF[0] << "  " << xVecF[1] << "  " << xVecF[2] << std::endl;
                    std::cout << "xVecC " << xVecC[0] << "  " << xVecC[1] << "  " << xVecC[2] << std::endl;
                    
                    std::cout << "PressureF " << PressureF << std::endl;
                    std::cout << "PressureC " << PressureC << std::endl;
                    
                    std::cout << "GradPressureF " << GradPressureF << std::endl;
                    std::cout << "GradPressureC " << GradPressureC << std::endl;
                }
                
                //erro L2 da pressao
                REAL diffP = 0;
                diffP = PressureC[0] - PressureF[0];
                errors[1] += weight*(diffP*diffP);

                //erro semi H1 da pressao
                TPZManVector<REAL,2> diffGrad(2,0.);
                diffGrad[0] = GradPressureC[0]-GradPressureF[0];
                diffGrad[1] = GradPressureC[1]-GradPressureF[1];
                errors[2] += weight*(diffGrad[0]*diffGrad[0] + diffGrad[1]*diffGrad[1]);

                //erro H1 para a pressao
                //errors[0] += errors[1] + errors[2];
                
            }
            

            
        }
        
    }
    
    if (IsFlux) {
        errors[2] = errors[0] + errors[1];
        errors[0] = sqrt(errors[0]);
        errors[1] = sqrt(errors[1]);
        errors[2] = sqrt(errors[2]);
        out << "\n";
        out << "Erros associados ao fluxo na norma L2\n";
        out << "Norma L2  para fluxo = " << errors[0] << endl;
        out << "Norma L2 para divergente = " << errors[1] << endl;
        out << "Norma Hdiv para o fluxo = " << errors[2] << endl;
    }
    else
    {
        //erro H1 para a pressao
        errors[0] = errors[1] + errors[2];
        errors[0] = sqrt(errors[0]);
        errors[1] = sqrt(errors[1]);
        errors[2] = sqrt(errors[2]);
        out << "\n";
        out << "Erros associados a pressao nas normas L2 e H1\n";
        out << "Norma H1 para a pressao = " << errors[0] << endl;
        out << "Norma L2 para a pressao = " << errors[1] << endl;
        out << "Norma semi-H1 para a pressao = " << errors[2] << endl;
    }
    
}

void ComputeFatherqsi(TPZGeoEl *Son, TPZGeoEl *Father,TPZManVector<REAL,3> &Sonqsi, TPZManVector<REAL,3> &Fatherqsi, TPZTransform<> &Transformation)
{
    
    if ((!Son && !Father) && (Sonqsi.size() != Fatherqsi.size() && Father != Son->Father())) {
        DebugStop();
    }
    
    int side, dimSon, dimFather;
    dimSon      = Son->Dimension();
    dimFather   = Father->Dimension();
    side = Father->NSides()-1;
    
    if (dimSon == dimFather) {
        TPZTransform<> tr(dimSon);
        Transformation = Son->BuildTransform2(side, Father, tr);
        //        int son = Father->WhichSubel();
        //        Transformation = Father->GetTransform(side, son);
    }
    else
    {
        std::cout<< "This method works just for volumetri-volumetric transformation ";
        DebugStop();
        
    }
    
    Transformation.Apply(Sonqsi, Fatherqsi);
    
    
}

//void ComputeFluxError(TPZCompMesh *CmeshFlux,int64_t NivelRef, std::ostream &out){
//    
//    TPZCompMesh * cmesh = CmeshFluxTeo;
//
//
//    TPZManVector<REAL,2> errors(3,0.);
//    int dimmesh = cmesh->Dimension();
//    int nel = cmesh->NElements();
//    int iel;
//    int64_t nivel = NivelMaxi-1;
//    
//    for(iel=0; iel<nel; iel++)
//    {
//        TPZCompEl *cel = cmesh->ElementVec()[iel];
//        if(!cel) continue;
//        
//        int dimcel = cel->Dimension();
//        if(dimcel != dimmesh) continue;
//        
//        TPZGeoEl *gel = cel->Reference();
//        TPZGeoEl *FatherGel=gel->Father();
//        
//        while (nivel > NivelRef) {
//            TPZGeoEl *FatherGel=FatherGel->Father();
//            nivel--;
//        }
//        
//        int64_t ielt = FatherGel->Id();
//        std::cout<<iel<<std::endl;
//        std::cout<<ielt<<std::endl;
//        
//        
//        TPZCompEl *celn = CmeshFlux->Element(FatherGel->Id());
//        
//        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 2);
//        TPZManVector<int,3> prevorder(dimcel), maxorder(dimcel,intrule->GetMaxOrder());
//        intrule->GetOrder(prevorder);
//        intrule->SetOrder(maxorder);
//        
//        TPZInterpolationSpace *spt = dynamic_cast<TPZInterpolationSpace *>(cel);
//        TPZInterpolationSpace *spn = dynamic_cast<TPZInterpolationSpace *>(celn);
//        
//        TPZMaterialData datat,datan;
//        spt->InitMaterialData(datat);
//        spn->InitMaterialData(datan);
//        
//        const int npoints = intrule->NPoints();
//        TPZManVector<REAL,3> qsi(dimcel), qsiFather(dimcel), xVec(3), XVecFather(3);
//        
//        
//        int Matid = cel->Reference()->MaterialId();
//        TPZMaterial  *material = cmesh->FindMaterial(MatId);
//        TPZMatPoisson3d *matp = dynamic_cast<TPZMatPoisson3d *>(material);
//        int qindex = matp->VariableIndex("Flux");
//        
//        TPZVec<STATE> QfluxSon,QfluxFather;
//        
//        for(int ip = 0; ip < npoints; ip++)
//        {
//            REAL weight;
//            intrule->Point(ip,qsi,weight);
//            
//            TPZTransform<> Transformation;
//            ComputeFatherqsi(gel, FatherGel, qsi, qsiFather, Transformation);
//            
//            gel->X(qsi, xVec);
//            FatherGel->X(qsiFather, XVecFather);
//            
//            std::cout << "qsi " << qsi[0] << "  " << qsi[1] << std::endl;
//            std::cout << "qsiFather " << qsiFather[0] << "  " << qsiFather[1] << std::endl;
//            std::cout << "xVec " << xVec[0] << "  " << xVec[1] << "  " << xVec[2] << std::endl;
//            std::cout << "XVecFather " << XVecFather[0] << "  " << XVecFather[1] << "  " << XVecFather[2] << std::endl;
//            
////            spt->ComputeShape(qsi, datat.x, datat.jacobian, datat.axes, datat.detjac, datat.jacinv, datat.phi, datat.dphix);
////            spn->ComputeShape(qsiFather, datan.x, datan.jacobian, datan.axes, datan.detjac, datan.jacinv, datan.phi, datan.dphix);
////            
////            weight *= fabs(datat.detjac);
////            spt->ComputeSolution(qsi,datat);
//            
//            
////            cel->Solution(qsi, qindex, QfluxSon);
////            celn->Solution(qsiFather, qindex, QfluxFather);
//            
//            std::cout << "QfluxSon " << QfluxSon[0] << "  " << QfluxSon[1] << std::endl;
//            std::cout << "QfluxFather " << QfluxFather[0] << "  " << QfluxFather[1] << std::endl;
//            
//            
////            TPZManVector<REAL,2> fluxn(2,0.);
////            REAL divfluxon;
////            fluxn[0]=datan.sol[0][0];
////            fluxn[1]=datan.sol[0][1];
////            divfluxon =  datan.dsol[0](0,0)+datan.dsol[0](1,1);
////
////            TPZManVector<REAL,2> fluxt(2,0.);
////            REAL divfluxot;
////            fluxt[0]=datat.sol[0][0];
////            fluxt[1]=datat.sol[0][1];
////            divfluxot =  datat.dsol[0](0,0)+datat.dsol[0](1,1);
////            
////            TPZManVector<REAL,2> diff(2,0.);
////            diff[0] = fluxn[0] - fluxt[0] ;
////            diff[1] = fluxn[1] - fluxt[1] ;
////            
////            //erro L2 do fluxo
////            errors[0] += weight*(diff[0]*diff[0] + diff[1]*diff[1]);
////            
////            //erro L2 do divergente do fluxo
////            REAL diffDiv = abs(divfluxon - divfluxot);
////            errors[1] += weight*diffDiv*diffDiv;
//            
//            //erro Hdiv para o fluxo
//            //            errors[2] += errors[0] + errors[1];
//        }
//        intrule->SetOrder(prevorder);
//    }
//    
//    //erro Hdiv para o fluxo
//    errors[2] = errors[0] + errors[1];
//    
//    errors[0] = sqrt(errors[0]);
//    errors[1] = sqrt(errors[1]);
//    errors[2] = sqrt(errors[2]);
//    
//    out << "\n";
//    out << "Erros associados ao fluxo na norma L2\n";
//    out << "Norma L2  para fluxo = " << errors[0] << endl;
//    out << "Norma L2 para divergente = " << errors[1] << endl;
//    out << "Norma Hdiv para o fluxo = " << errors[2] << endl;
//    
//}///method
//

void ComputePressureError(TPZCompMesh *CmeshPres,int64_t NivelRef, std::ostream &out){
 
    TPZCompMesh * cmesh = CmeshPresTeo;
    TPZManVector<REAL,2> errors(3,0.);
    errors.Fill(0.);
    
    int dimmesh = cmesh->Dimension();
    int nel = cmesh->NElements();
    int iel;
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        int dimcel = cel->Dimension();
        if(dimcel != dimmesh) continue;
        
        TPZGeoEl *gel = cel->Reference();
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 2);
        TPZManVector<int,3> prevorder(dimcel), maxorder(dimcel,intrule->GetMaxOrder());
        intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        const int npoints = intrule->NPoints();
        TPZManVector<REAL,3> qsi(dimcel), xVec(3);
        
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            intrule->Point(ip,qsi,weight);
            sp->ComputeShape(qsi, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            weight *= fabs(data.detjac);
            sp->ComputeSolution(qsi,data);
            
            TPZManVector<STATE,2> gradP(2,0.);
            REAL diffP;
            
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            gel->X(qsi,xVec);
            SolExataPress(xVec, uExato, duExato);

            //erro L2 da pressao
            STATE solP = data.sol[0][0];
            diffP = solP - uExato[0];
            errors[1] += weight*(diffP*diffP);
            
            //erro semi H1 da pressao
            TPZFNMatrix<3,REAL> dsoldx;
            TPZFMatrix<REAL> dsoldaxes(2,1);
            dsoldaxes(0,0) = data.dsol[0][0];
            dsoldaxes(1,0) = data.dsol[0][1];
            
            TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, data.axes);
            TPZManVector<REAL,2> diffGrad(2,0.);
            diffGrad[0] = dsoldx(0,0)-duExato(0,0);
            diffGrad[1] = dsoldx(1,0)-duExato(1,0);
            errors[2] += weight*(diffGrad[0]*diffGrad[0] + diffGrad[1]*diffGrad[1]);
            
            //erro H1 para a pressao
            //errors[0] += errors[1] + errors[2];
        }
        intrule->SetOrder(prevorder);
    }
    
    //erro H1 para a pressao
    errors[0] = errors[1] + errors[2];
    
    errors[0] = sqrt(errors[0]);
    errors[1] = sqrt(errors[1]);
    errors[2] = sqrt(errors[2]);
    
    out << "\n";
    out << "Erros associados a pressao nas normas L2 e H1\n";
    out << "Norma H1 para a pressao = " << errors[0] << endl;
    out << "Norma L2 para a pressao = " << errors[1] << endl;
    out << "Norma semi-H1 para a pressao = " << errors[2] << endl;
    
    
}///method


void PosProcessMultph(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZManVector<std::string,10> scalnames(2), vecnames(3);
    vecnames[0]  = "Flux";
    vecnames[1]  = "GradFluxX";
    vecnames[2]  = "GradFluxY";
//    vecnames[3]  = "ExactFlux";
    scalnames[0] = "Pressure";
    scalnames[1] = "DivFlux";
 //   scalnames[2] = "ExactPressure";
    

    const int dim = 2;
    int div =1;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    //	std::ofstream out("malha.txt");
    //	an.Print("nothing",out);
    
}



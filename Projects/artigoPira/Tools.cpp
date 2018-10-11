//
//  Tools.cpp
//  pz
//
//  Created by Denise De Siqueira on 04/05/2018.
//

#include "Tools.h"

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
#include "tpzchangeel.h"
#include "tpzarc3d.h"
#include "tpzquadraticquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzgengrid.h"

#include "pzlog.h"

#include "pzl2projection.h"
#include "pzmultiphysicselement.h"
#include "pzintel.h"
#include "TPZVTKGeoMesh.h"

#include "TPZMatDualHybridPoisson.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "tpzgeoelmapped.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZIntQuadQuarterPoint.h"
#include "tpzintpoints.h"
#include "pzquad.h"
#include "tpzhierarquicalgrid.h"
#include "TPZVecL2.h"

#include <iostream>
#include <math.h>


using namespace std;

int const matId =1;
int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;

int const bc7=-7;
int const bc8=-8;

int const dirichlet =0;
int const neumann = 1;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("Tools.main"));
#endif

TPZGeoMesh *GMesh(){
    
    int dim = 2;
    int Qnodes = 6;
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    gmesh->SetDimension(dim);
    TPZGeoEl *gel = 0;
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolLine(2);
    
    TPZStack<TPZGeoEl *> toChange;
    TPZStack<int> targetSide;
    
    //indice dos nos
    long id = 0;
    REAL valx, dx = 0.5;
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = -0.5 + xi*dx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,0. );//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = 0.5 - xi*dx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,0.5 );//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    //indice dos elementos
    id = 0;
    
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 4;
    TopolQuad[3] = 5;
    gel = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
    toChange.Push(gel);
    targetSide.Push(1);
    
    id++;
    
    TopolQuad[0] = 1;
    TopolQuad[1] = 2;
    TopolQuad[2] = 3;
    TopolQuad[3] = 4;
    gel = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
    toChange.Push(gel);
    targetSide.Push(0);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    gel = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
    toChange.Push(gel);
    targetSide.Push(1);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 2;
    gel = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
    toChange.Push(gel);
    targetSide.Push(0);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
    id++;
    
    TopolLine[0] = 4;
    TopolLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
    id++;
    
    TopolLine[0] = 5;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc5,*gmesh);
    
    gmesh->BuildConnectivity();
    
    //    if(QuarterPoint){
    //        for (long el=0; el<toChange.size(); el++) {
    //            TPZGeoEl *gel = toChange[el];
    //    //        TPZChangeEl::ChangeToQuadratic(gmesh, gel->Index());
    //            TPZChangeEl::ChangeToQuarterPoint(gmesh, gel->Index(), targetSide[el]);
    //        }
    //    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    return gmesh;
}


TPZGeoMesh *GMeshTeste(){
    
    REAL Lx=0.5;
    REAL Ly=0.5;
    int dim = 2;
    
    int Qnodes = 4;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    
    
    gmesh->SetDimension(dim);
    
    
    TPZGeoEl *gel = 0;
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
    TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
    
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
    
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 2;
    TopolQuad[3] = 3;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 2;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
    
    
    gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
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
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"gmesh depois de refinar uniformemente\n";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    
}


TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
    int dim = gmesh->Dimension();
    
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    
    //    TPZVecL2 *material = new TPZVecL2(matId);
    //    material->SetDimension(dim);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    TPZAutoPointer<TPZFunction<STATE> > solExata;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(SolExataSteklov, 5);
    dum->SetPolynomialOrder(20);
    solExata = dum;
    material->SetForcingFunctionExact(solExata);
    
    //Inserir condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    BCond1 = material->CreateBC(mat, bc1, dirichlet, val1, val2);//Dirichlet nulo
    BCond2 = material->CreateBC(mat, bc2, dirichlet, val1, val2);//Neumann nulo
    BCond3 = material->CreateBC(mat, bc3, dirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, bc4, dirichlet, val1, val2);
    BCond5 = material->CreateBC(mat, bc5, dirichlet, val1, val2);
    
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    
    TPZMaterial * BCond7;
    TPZMaterial * BCond8;
    if (dim ==3) {
        BCond7 = material->CreateBC(mat, bc7, dirichlet, val1, val2);//Dirichlet nulo
        BCond8 = material->CreateBC(mat, bc8, dirichlet, val1, val2);//Dirichlet nulo
        cmesh->InsertMaterialObject(BCond7);
        cmesh->InsertMaterialObject(BCond8);
    }
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();//ajusta as condicoes de contorno
    cmesh->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional fluxo\n";
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    return cmesh;
}

TPZCompMesh *CMeshFluxL2(TPZGeoMesh *gmesh, int pOrder, int nodeAtOriginId)
{
    /// criar materiais
    int dim = gmesh->Dimension();
    
    
    TPZVecL2 *material = new TPZVecL2(matId);
    material->SetDimension(dim);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(SolExataSteklov, 5);
    material->SetForcingFunctionExact(solExata);
    
    //Inserir condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    BCond1 = material->CreateBC(mat, bc1, dirichlet, val1, val2);//Dirichlet nulo
    BCond2 = material->CreateBC(mat, bc2, neumann, val1, val2);//Neumann nulo
    BCond3 = material->CreateBC(mat, bc3, neumann, val1, val2);
    BCond4 = material->CreateBC(mat, bc4, neumann, val1, val2);
    BCond5 = material->CreateBC(mat, bc5, neumann, val1, val2);
    
    //Set force function
    int intorder = 10;
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannDireito;
    TPZDummyFunction<STATE> *dumBC3 = new TPZDummyFunction<STATE>(NeumannDireita, 5);
    dumBC3->SetPolynomialOrder(intorder);
    bcmatNeumannDireito = dumBC3;
    BCond3->SetForcingFunction(bcmatNeumannDireito);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
    TPZDummyFunction<STATE> *dumBC4 = new TPZDummyFunction<STATE>(NeumannAcima, 5);
    dumBC4->SetPolynomialOrder(intorder);
    bcmatNeumannAcima = dumBC4;
    BCond4->SetForcingFunction(bcmatNeumannAcima);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEsquerdo;
    TPZDummyFunction<STATE> *dumBC5 = new TPZDummyFunction<STATE>(NeumannEsquerda, 5);
    dumBC5->SetPolynomialOrder(intorder);
    bcmatNeumannEsquerdo = dumBC5;
    BCond5->SetForcingFunction(bcmatNeumannEsquerdo);
    
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    
    
    cmesh->SetAllCreateFunctionsHDiv();
    //cmesh->SetAllCreateFunctionsHDivFull();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();//ajusta as condicoes de contorno
    cmesh->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional fluxo\n";
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    return cmesh;
}


TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder)
{
    bool triang = false;
    /// criar materiais
    int dim = gmesh->Dimension();
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(matId,dim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    //cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    
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
        if(!celdisc) continue;
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
    
    cmesh->AdjustBoundaryElements();//ajusta as condicoes de contorno
    cmesh->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional pressao\n";
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    return cmesh;
}

TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson * &mymaterial,bool QuarterPointRule){
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=1;
    int dim = gmesh->Dimension();
    
    REAL coefk=1.;
    mymaterial = new TPZMixedPoisson(MatId,dim);
    mymaterial->SetPermeability(coefk);
    
    TPZFMatrix<REAL> TensorK(dim,dim,0.);
    TPZFMatrix<REAL> InvK(dim,dim,0.);
    for(int i = 0; i<dim; i++) {
        TensorK(i,i)=coefk;
        InvK(i,i)=coefk;
    }
    mymaterial->SetPermeabilityTensor(TensorK, InvK);
    
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    TPZAutoPointer<TPZFunction<STATE> > solExata = NULL;
    
    if (dim ==2) {
        solExata = new TPZDummyFunction<STATE>(SolExataSteklov, 5);
        mymaterial->SetForcingFunctionExact(solExata);
    }else{
        solExata = new TPZDummyFunction<STATE>(SolExata3D, 5);
        mymaterial->SetForcingFunctionExact(solExata);
        
        TPZAutoPointer<TPZFunction<STATE>> source = new TPZDummyFunction<STATE> (f_3D, 5);
        mymaterial->SetForcingFunction(source);
        
    }
    
    //Inserir condicoes de contorno
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
#ifdef SmoothSol
    
    TPZMaterial * BCond7;
    TPZMaterial * BCond8;
    
    BCond1 = mymaterial->CreateBC(mat, bc1, dirichlet, val1, val2);//Dirichlet nulo
    BCond2 = mymaterial->CreateBC(mat, bc2, dirichlet, val1, val2);//Neumann nulo
    BCond3 = mymaterial->CreateBC(mat, bc3, dirichlet, val1, val2);
    BCond4 = mymaterial->CreateBC(mat, bc4, dirichlet, val1, val2);
    BCond5 = mymaterial->CreateBC(mat, bc5, dirichlet, val1, val2);
    BCond7 = mymaterial->CreateBC(mat, bc7, dirichlet, val1, val2);
    BCond8 = mymaterial->CreateBC(mat, bc8, dirichlet, val1, val2);
    
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);
    mphysics->InsertMaterialObject(BCond7);
    mphysics->InsertMaterialObject(BCond8);
    
#else
    
    BCond1 = mymaterial->CreateBC(mat, bc1, dirichlet, val1, val2);//Dirichlet nulo
    BCond2 = mymaterial->CreateBC(mat, bc2, neumann, val1, val2);//Neumann nulo
    BCond3 = mymaterial->CreateBC(mat, bc3, neumann, val1, val2);
    BCond4 = mymaterial->CreateBC(mat, bc4, neumann, val1, val2);
    BCond5 = mymaterial->CreateBC(mat, bc5, neumann, val1, val2);
    
    //Set force function
    int intorder = 10;
    
    if (dim == 2 ) {
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannDireito;
        TPZDummyFunction<STATE> *dumBC3 = new TPZDummyFunction<STATE>(NeumannDireita, 5);
        dumBC3->SetPolynomialOrder(intorder);
        bcmatNeumannDireito = dumBC3;
        BCond3->SetForcingFunction(bcmatNeumannDireito);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
        TPZDummyFunction<STATE> *dumBC4 = new TPZDummyFunction<STATE>(NeumannAcima, 5);
        dumBC4->SetPolynomialOrder(intorder);
        bcmatNeumannAcima = dumBC4;
        BCond4->SetForcingFunction(bcmatNeumannAcima);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEsquerdo;
        TPZDummyFunction<STATE> *dumBC5 = new TPZDummyFunction<STATE>(NeumannEsquerda, 5);
        dumBC5->SetPolynomialOrder(intorder);
        bcmatNeumannEsquerdo = dumBC5;
        BCond5->SetForcingFunction(bcmatNeumannEsquerdo);
    }else{
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannDireito;
        TPZDummyFunction<STATE> *dumBC3 = new TPZDummyFunction<STATE>(NeumannDireita_3D, 5);
        dumBC3->SetPolynomialOrder(intorder);
        bcmatNeumannDireito = dumBC3;
        BCond3->SetForcingFunction(bcmatNeumannDireito);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
        TPZDummyFunction<STATE> *dumBC4 = new TPZDummyFunction<STATE>(NeumannAcima_3D, 5);
        dumBC4->SetPolynomialOrder(intorder);
        bcmatNeumannAcima = dumBC4;
        BCond4->SetForcingFunction(bcmatNeumannAcima);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEsquerdo;
        TPZDummyFunction<STATE> *dumBC5 = new TPZDummyFunction<STATE>(NeumannEsquerda_3D, 5);
        dumBC5->SetPolynomialOrder(intorder);
        bcmatNeumannEsquerdo = dumBC5;
        BCond5->SetForcingFunction(bcmatNeumannEsquerdo);
    }
    
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);
    
    TPZMaterial * BCond7;
    TPZMaterial * BCond8;
    
    if (dim == 3) {
        BCond7 = mymaterial->CreateBC(mat, bc7, dirichlet, val1, val2);//Dirichlet nulo
        BCond8 = mymaterial->CreateBC(mat, bc8, dirichlet, val1, val2);//Dirichlet nulo
        mphysics->InsertMaterialObject(BCond7);
        mphysics->InsertMaterialObject(BCond8);
    }
    
#endif
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    
    //Condensacao Estatica
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    mphysics->SetDimModel(dim);
    
    // create condensed elements
    // increase the NumElConnected of one pressure connects in order to prevent condensation
    if(!QuarterPointRule){
        mphysics->ComputeNodElCon();
        for (long icel=0; icel < mphysics->NElements(); icel++) {
            TPZCompEl  * cel = mphysics->Element(icel);
            if(!cel) continue;
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    break;
                }
            }
            new TPZCondensedCompEl(cel);
        }
    }
    mphysics->CleanUpUnconnectedNodes();
    mphysics->ExpandSolution();
    
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional multifisica\n";
        mphysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    return mphysics;
}


void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZManVector<std::string,10> scalnames(3), vecnames(0);
    //vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
    
    scalnames[1] = "ExactPressure";
    // vecnames[0] = "ExactFlux";
    scalnames[2] = "POrder";
    
    const int dim = mphysics->Reference()->Dimension();
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    //    std::ofstream out("malha.txt");
    //    an.Print("nothing",out);
    
}


//refinamento uniforme em direcao ao no
void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide){
    
    for (int idivide = 0; idivide < divide; idivide++){
        const int nels = gmesh->NElements();
        TPZManVector< TPZGeoEl * > allEls(nels);
        for(int iel = 0; iel < nels; iel++){
            allEls[iel] = gmesh->ElementVec()[iel];
        }
        
        // std::cout << "idivide = " << idivide << std::endl;
        
        for(int iel = 0; iel < nels; iel++){
            TPZGeoEl * gel = allEls[iel];
            if(!gel) continue;
            if(gel->HasSubElement()) continue;
            int nnodes = gel->NNodes();
            int found = -1;
            for(int in = 0; in < nnodes; in++){
                if(gel->NodePtr(in)->Id() == nodeAtOriginId){
                    found = in;
                    break;
                }
            }///for in
            if(found == -1) continue;
            
            MElementType gelT = gel->Type();
            TPZAutoPointer<TPZRefPattern> uniform = gRefDBase.GetUniformRefPattern(gelT);
            if(!uniform){
                DebugStop();
            }
            
            gel->SetRefPattern(uniform);
            //std::cout << "Dividing element " << gel->Index() << " idivide = "<< idivide << std::endl;
            TPZManVector<TPZGeoEl*,4> filhos;
            gel->Divide(filhos);
            
        }///for iel
    }//idivide
    
    gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"gmesh depois de refinar direcionalmente\n";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
}///void
////

void DirectUniRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide){
    
    ///Refinamento
    //    gRefDBase.InitializeUniformRefPattern(EOned);
    //    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
    for (int idivide = 0; idivide < divide; idivide++){
        const int nels = gmesh->NElements();
        TPZManVector< TPZGeoEl * > allEls(nels);
        for(int iel = 0; iel < nels; iel++){
            allEls[iel] = gmesh->ElementVec()[iel];
        }
        
        // std::cout << "idivide = " << idivide << std::endl;
        
        for(int iel = 0; iel < nels; iel++){
            TPZGeoEl * gel = allEls[iel];
            if(!gel) continue;
            if(gel->HasSubElement()) continue;
            int nnodes = gel->NNodes();
            int found = -1;
            for(int in = 0; in < nnodes; in++){
                if(gel->NodePtr(in)->Id() != nodeAtOriginId){
                    found = in;
                    break;
                }
            }///for in
            if(found == -1) continue;
            
            MElementType gelT = gel->Type();
            TPZAutoPointer<TPZRefPattern> uniform = gRefDBase.GetUniformRefPattern(gelT);
            if(!uniform){
                DebugStop();
            }
            
            gel->SetRefPattern(uniform);
            //std::cout << "Dividing element " << gel->Index() << " idivide = "<< idivide << std::endl;
            TPZManVector<TPZGeoEl*,4> filhos;
            gel->Divide(filhos);
            
        }///for iel
    }//idivide
    
    gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"gmesh depois de refinar direcionalmente\n";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
}///void


/////

void QuarterPointRef(TPZGeoMesh *gmesh, int nodeAtOriginId)
{
    TPZStack<TPZGeoEl *> toChange;
    TPZStack<int> targetSide;
    int found = 0;
    
    int nels = gmesh->NElements();
    
    for(int iel = 0; iel < nels; iel++){
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        
        if(!gel) continue;
        if(gel->HasSubElement()) continue;
        int nnodes = gel->NNodes();
        
        for(int in = 0; in < nnodes; in++){
            if(gel->NodePtr(in)->Id() == nodeAtOriginId){
                
                toChange.Push(gel);
                targetSide.Push(in);
                found++;
                break;
            }
        }
        if(found == 4) break;
    }
    
    
    for (int el=0; el<toChange.size(); el++) {
        TPZGeoEl *gel = toChange[el];
        long index = gel->Index();
        TPZChangeEl::ChangeToQuarterPoint(gmesh, gel->Index(), targetSide[el]);
        gel = gmesh->Element(index);
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"gmesh com quarterpoint\n";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
}



void PrefinamentoRibsHybridMesh(TPZCompMesh *cmesh){
    const int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZCompEl * cel = cmesh->ElementVec()[iel];
        if(!cel)continue;
        TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
        if(!face)continue;
        TPZInterpolationSpace *TwoD = NULL;
        TPZInterpolationSpace *Rib = NULL;
        if(face->LeftElement()->Reference()->Dimension()==2){
            TwoD = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
            Rib = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
        }
        else{
            Rib = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
            TwoD = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
        }
        if(!Rib || !TwoD) DebugStop();
        int porder2D = TwoD->MaxOrder();
        int pOrder1D = Rib->MaxOrder();
        int neworder = pOrder1D > porder2D ? pOrder1D : porder2D;
        Rib->PRefine(neworder);
    }
    
}

void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder){
    if(ndiv<1) return;
    int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        if(!sp) continue;
        int level = sp->Reference()->Level();
        TPZGeoEl * gel = sp->Reference();
        //int fatherlevel=gel->Father()->Level();
        if(gel->Dimension()==2){
            int ordem= 0;
            ordem=porder + (ndiv - 1) + (level);
            std::cout<<"level "<< level<<" ordem "<<ordem<<std::endl;
            sp->PRefine(porder + (ndiv - 1) + (level));
        }
    }
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional apos pRefinamento\n";
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
}

void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    u.Resize(1, 0.);
    du.Resize(2, 1);
    du(0,0)=0.;
    du(1,0)=0.;
    
    const REAL n = 0;//n=3 para solucao suave e n=0 para solucao singular
    REAL x = loc[0];
    REAL y = loc[1];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = pow((REAL)2,0.25 + n/2.)*pow(r,0.5 + n)*cos((0.5 + n)*t);
    u[0] = sol;

        du(0,0) = -pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
        du(1,0) = -pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));

}

void SolExata3D(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    u.Resize(1, 0.);
    du.Resize(4, 1);
    du(0,0)=0.;
    du(1,0)=0.;
    du(2,0)=0.;
    
    
    REAL x = loc[0];
    REAL y = loc[1];
    REAL z = loc[2];
    
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    
#ifdef SmoothSol
    
    u[0] = sqrt(r)*(0.5 - z)*z*cos(t/2.);
    
#else
    
    u[0] = pow(pow(x,2) + pow(y,2),0.25)*(0.5 - z)*z*cos(atan2(y,x)/2.);
    
    
    du(0,0) = (x*z*(-1 + 2*z)*cos(t/2.))/(4.*sqrt(r)*sqrt(pow(x,2) + pow(y,2))) +
    (sqrt(r)*y*z*(-1 + 2*z)*sin(t/2.))/(4.*sqrt(pow(x,2) + pow(y,2)));
    
    du(1,0) = (y*z*(-1 + 2*z)*cos(t/2.))/(4.*sqrt(r)*sqrt(pow(x,2) + pow(y,2))) -
    (sqrt(r)*x*z*(-1 + 2*z)*sqrt(t/2.))/(4.*sqrt(pow(x,2) + pow(y,2)));
    
    du(2,0) = (sqrt(r)*(-1 + 4*z)*cos(t/2.))/2.;
    
    //termo fonte
    du(3,0)= ((z - 2*pow(z,2) + pow(r,2)*(16 + z - 2*pow(z,2)))*cos(t/2.))/(8.*pow(r,1.5));
#endif
    
}

void f_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &f){
    
    f.Resize(1, 0.);
    REAL x = loc[0];
    REAL y = loc[1];
    REAL z = loc[2];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    
#ifdef SmoothSol
    f[0] = -((-1 + 4*pow(x,2))*z*(-1 + 2*z))/4. +
    pow(y,2)*(0.5 - 2*pow(x,2) + z - 2*pow(z,2)) +
    y*(-0.25 + pow(x,2) - z/2. + pow(z,2));
#else
    //    f[0] = ((pow(pow(x,2) + pow(y,2),0.25)*(-1 + 4*z)*cos(atan2(y,x)/2.))/2.);
    f[0] = ((z - 2*pow(z,2) + pow(r,2)*(16 + z - 2*pow(z,2)))*cos(t/2.))/(8.*pow(r,1.5));
#endif
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
    
    //result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {+1,0};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    // result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {0,+1};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    //result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannEsquerda_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[3] = {-1,0,0};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(4,1);
    SolExata3D(loc, u, du);
    
    //    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
    result[0] *= -1.0;
}

void NeumannDireita_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[3] = {+1,0,0};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(4,1);
    SolExata3D(loc, u, du);
    
    //    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
    result[0] *= -1.0;
}

void NeumannAcima_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[3] = {0,+1,0};
    
    TPZManVector<STATE> u(1);
    TPZFNMatrix<10,STATE> du(4,1);
    SolExata3D(loc,u,du);
    
    //    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1]+du(2,0)*normal[2];
    result[0] *= -1.0;
}

#include "pzbstrmatrix.h"
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numThreads)
{
    if(numThreads==0){
        
        TPZSkylineStructMatrix strmat(fCmesh); //caso simetrico
        an.SetStructuralMatrix(strmat);
        
    }else
    {
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
        strmat.SetDecomposeType(ELDLt);
        strmat.SetNumThreads(numThreads);
        an.SetStructuralMatrix(strmat);
    }
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Run();
    
    //Saida de Dados: solucao e  grafico no VT
    //ofstream file("../Solout.txt");
    //an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out)
{
    long nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        if (dim == 2) {
            cel->EvaluateError(SolExataSteklov, elerror, false);
        }else{
            cel->EvaluateError(SolExata3D, elerror, false);
        }
        
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "\nErrors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
}
void ErroL2NoElemento(TPZCompMesh *Pmesh, std::ostream &out,  int nodeAtOriginId)
{
    long nel = Pmesh->NElements();
    int dim = Pmesh->Dimension();
    TPZManVector<REAL,10> globerrors(10,0.);
    int found = 0;
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = Pmesh->ElementVec()[el];
        if (!cel) continue;
        
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) continue;
        
        
        int nnodes = gel->NNodes();
        for(int in = 0; in < nnodes; in++)
        {
            if(gel->NodePtr(in)->Id() != nodeAtOriginId){
                continue;
            }else
            {
                TPZManVector<REAL,10> elerror(10,0.);
                cel->EvaluateError(SolExataSteklov, elerror, false);
                int nerr = elerror.size();
                for (int i=0; i<nerr; i++)
                {
                    globerrors[i] += elerror[i]*elerror[i];
                }
                found++;
            }
            if(found==2) break;
        }
        if(found==2) break;
    }
    out << "\nErrors associated with L2 space\n";
    out << "L2 Norm for pressure = "    << sqrt(globerrors[1]) << endl;
    
    if(found!=2)
    {
        DebugStop();
    }
}

void ErroHDivNoElemento(TPZCompMesh *hdivmesh, std::ostream &out,  int nodeAtOriginId)
{
    long nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globerrors(10,0.);
    int found = 0;
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) continue;
        
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) continue;
        
        
        int nnodes = gel->NNodes();
        for(int in = 0; in < nnodes; in++)
        {
            if(gel->NodePtr(in)->Id() != nodeAtOriginId){
                continue;
            }else
            {
                TPZManVector<REAL,10> elerror(10,0.);
                cel->EvaluateError(SolExataSteklov, elerror, false);
                int nerr = elerror.size();
                for (int i=0; i<nerr; i++)
                {
                    globerrors[i] += elerror[i]*elerror[i];
                }
                found++;
            }
            if(found==2) break;
        }
        if(found==2) break;
    }
    out << "\nErrors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
    
    if(found!=2)
    {
        DebugStop();
    }
}


void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out)
{
    long nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    globalerrors.Fill(0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        if (dim == 2) {
            cel->EvaluateError(SolExataSteklov, elerror, false);
        }
        else{
            cel->EvaluateError(SolExata3D, elerror, false);
        }
        int nerr = elerror.size();
        // globalerrors.resize(nerr);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "\nErrors associated with L2 space"<< endl;
    out << "L2 Norm for Pressure = "    << sqrt(globalerrors[1]) << endl;
}


//----------------------------------
TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder, bool ismultiplierH1){
    //TPZCompEl::SetgOrder(porder);
    TPZCompMesh *comp = new TPZCompMesh(&gmesh);
    int dim = gmesh.Dimension();
    comp->SetDimModel(dim);
    comp->SetDefaultOrder(porder);
    
    comp->ApproxSpace().CreateDisconnectedElements(true);
    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    
    // Criar e inserir os materiais na malha
    REAL beta = 6;
    TPZMatDualHybridPoisson *mymaterial = new TPZMatDualHybridPoisson(1,0.,beta);
    TPZMaterial * automat(mymaterial);
    comp->InsertMaterialObject(automat);
    mymaterial->SetDimension(dim);
    
    
    // Condicoes de contorno
    TPZFMatrix<STATE> val1(1,1,1.),val2(1,1,0.);
    
    TPZMaterial *bnd = automat->CreateBC (automat,-1,2,val1,val2);//misto tbem
    bnd->SetForcingFunction(Dirichlet2,porder);
    comp->InsertMaterialObject(bnd);
    
    // Mixed
    bnd = automat->CreateBC (automat,-2,0,val1,val2);
    bnd->SetForcingFunction(Dirichlet,porder);
    comp->InsertMaterialObject(bnd);
    
    // Mixed
    val1(0,0)=0.;
    val2(0,0)=0.;
    bnd = automat->CreateBC (automat,-3,0,val1,val2);//Dirichlet nulo
    comp->InsertMaterialObject(bnd);
    
    // Mixed
    bnd = automat->CreateBC (automat,-4,0,val1,val2);
    bnd->SetForcingFunction(Dirichlet,porder);
    comp->InsertMaterialObject(bnd);
    
    // Mixed
    bnd = automat->CreateBC (automat,-5,0,val1,val2);
    bnd->SetForcingFunction(Dirichlet,porder);
    comp->InsertMaterialObject(bnd);
    
    // Mixed
    bnd = automat->CreateBC (automat,-6,0,val1,val2);
    bnd->SetForcingFunction(Dirichlet,porder);
    comp->InsertMaterialObject(bnd);
    
    // Ajuste da estrutura de dados computacional
    comp->ApproxSpace().CreateDisconnectedElements(true);
    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    std::set<int> matids;
    matids.insert(1);
    comp->AutoBuild(matids);
    
    comp->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    matids.clear();
    for (int i=2; i<=6; i++) {
        matids.insert(-i);
    }
    
    comp->SetDefaultOrder(porder);
    comp->AutoBuild(matids);
    comp->SetDimModel(2);
    
    comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
    comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
    comp->LoadReferences();
    comp->ApproxSpace().CreateInterfaceElements(comp,true);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    matids.insert(1);
    for (int i=2; i<=6; i++) {
        matids.insert(-i);
    }
    
    comp->ApproxSpace().Hybridize(*comp, matids, ismultiplierH1);
    comp->AdjustBoundaryElements();
    comp->CleanUpUnconnectedNodes();
    comp->ExpandSolution();
    
    comp->SetName("Malha Computacional Original");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    return comp;
}


void GroupElements(TPZCompMesh *cmesh)
{
    long neq = cmesh->NEquations();
    
    cmesh->LoadReferences();
    
    int nel = cmesh->NElements();
    std::set<TPZCompEl *> celset;
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();
        if (dim ==2) {
            celset.insert(cel);
        }
    }
    std::set<int> elgroupindices;
    
    for (std::set<TPZCompEl *>::iterator it = celset.begin(); it != celset.end(); it++) {
        
        std::list<TPZCompEl *> group;
        group.push_back(*it);
        TPZCompEl *cel = *it;
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZStack<TPZCompElSide> connected;
            TPZCompElSide celside(cel,is);
            celside.EqualLevelElementList(connected, false, false);
            int neq = connected.NElements();
            for (int eq=0; eq<neq; eq++) {
                TPZCompElSide eqside = connected[eq];
                TPZCompEl *celeq = eqside.Element();
                TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(celeq);
                if (!intface) {
                    continue;
                }
                TPZCompEl *left = intface->LeftElement();
                if (left == cel) {
                    //put in the group
                    group.push_back(intface);
                }
            }
        }
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
                (*it)->Print(sout);
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        int64_t index;
        TPZElementGroup *celgroup = new TPZElementGroup(*cmesh,index);
        elgroupindices.insert(index);
        for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
            celgroup->AddElement(*it);
        }
    }
    
    cmesh->ComputeNodElCon();
    
    neq = cmesh->NEquations();
    
    for (std::set<int>::iterator it = elgroupindices.begin(); it!=elgroupindices.end(); it++) {
        TPZCompEl *cel = cmesh->ElementVec()[*it];
        TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel);
    }
    
    neq = cmesh->NEquations();
    
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}


void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFNMatrix<10,REAL> fake(2,1);
    result[0] = loc[0]*loc[0];
}

TPZCompMesh *MalhaCompH1QP(TPZGeoMesh * gmesh,int ordem){
    /// criar materiais
    int dim = gmesh->Dimension();
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(matId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(ordem);
    cmesh->SetDimModel(dim);
    
    TPZAutoPointer<TPZFunction<STATE> > solExata;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(SolExataSteklov, 5);
    dum->SetPolynomialOrder(20);
    solExata = dum;
    material->SetForcingFunctionExact(solExata);
    
    //Inserir condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    BCond1 = material->CreateBC(mat, bc1, dirichlet, val1, val2);//Dirichlet nulo
    BCond2 = material->CreateBC(mat, bc2, neumann, val1, val2);//Neumann nulo
    BCond3 = material->CreateBC(mat, bc3, neumann, val1, val2);
    BCond4 = material->CreateBC(mat, bc4, neumann, val1, val2);
    BCond5 = material->CreateBC(mat, bc5, neumann, val1, val2);
    
    //Set force function
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannDireito;
    bcmatNeumannDireito = new TPZDummyFunction<STATE>(NeumannDireita, 5);
    BCond3->SetForcingFunction(bcmatNeumannDireito);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
    bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcima, 5);
    BCond4->SetForcingFunction(bcmatNeumannAcima);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannEsquerdo;
    bcmatNeumannEsquerdo = new TPZDummyFunction<STATE>(NeumannEsquerda, 5);
    BCond5->SetForcingFunction(bcmatNeumannEsquerdo);
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Fazendo auto build
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"malha computacional H1\n";
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif
    
    
    return cmesh;
}


void ChangeInternalConnectOrder(TPZCompMesh *mesh){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        int nshape2 = 0;
        
        if(cel->Dimension()== dim)
        {
            TPZConnect &conel = cel->Connect(ncon-1);
            corder = conel.Order();
            nshape = conel.NShape();
            int neworder = corder + 1;
            conel.SetOrder(neworder,cel->ConnectIndex(ncon-1));
            nshape = 2*(corder + 1)*(corder + 2);
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            intel->SetPreferredOrder(neworder);
            nshape2 = intel->NConnectShapeF(ncon-1,neworder);
            
            conel.SetNShape(nshape);
            mesh->Block().Set(conel.SequenceNumber(),nshape);
        }
    }
    mesh->CleanUpUnconnectedNodes();
    mesh->ExpandSolution();
}

void ChangeSideConnectOrderConnects(TPZCompMesh *mesh, int reduction_value){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    int neworder =0;
    
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        
        if(cel->Dimension()== dim)
        {
            for (int icon=0; icon<ncon-1; icon++)
            {
                TPZConnect &conel = cel->Connect(icon);
                corder = conel.Order();
                nshape = conel.NShape();
                
                if(corder!=neworder)
                {
                    neworder = corder - reduction_value;
                    conel.SetOrder(neworder, cel->ConnectIndex(icon));
                    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                    nshape = intel->NConnectShapeF(icon,neworder);
                    conel.SetNShape(nshape);
                    mesh->Block().Set(conel.SequenceNumber(),nshape);
                }
                
                if(conel.HasDependency())
                {
                    int order = corder - reduction_value;
                    conel.SetOrder(order,cel->ConnectIndex(icon));
                    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                    nshape = intel->NConnectShapeF(icon,order);
                    conel.SetNShape(nshape);
                    mesh->Block().Set(conel.SequenceNumber(),nshape);
                }
                
            }
        }
    }
    
    mesh->CleanUpUnconnectedNodes();
    mesh->ExpandSolution();
    
}



void IntegrationRuleConvergence(bool intQuarterPoint){
    
    std::ofstream saidaMat("../ConvIntRuleP2.nb");
    
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    HDivPiola = 1;
    int p = 0;
    
    for(p=2; p<3; p++)
    {
        TPZGeoMesh *gmesh = GMesh();
        UniformRefine(gmesh,1);
        for(long el=0; el < gmesh->NElements(); el++)
        {
            TPZGeoEl *gel = gmesh->Element(el);
            gel->SetFather(-1);
        }
        
        //refinamento quarter point proximo do no de id=1
        int nodeAtOriginId = 1;
        QuarterPointRef(gmesh, nodeAtOriginId);
        
        //refinamento proximo do no de id=1
        DirectionalRef(gmesh, nodeAtOriginId, 0);
        
        
        //cmeshL2 flux
        TPZCompMesh *cmesh = CMeshFluxL2(gmesh, p, nodeAtOriginId);
        //        {
        //            std::ofstream malhaOut("cmeshL2.txt");
        //            cmesh->Print(malhaOut);
        //        }
        
        
        /////--------------------------------
        int nelem = cmesh->NElements();
        int dim  = cmesh->Dimension();
        
        TPZCompEl *cel = NULL;
        TPZGeoEl * gel = NULL;
        
        int cornerpoint = -1;
        for(int i=0; i<nelem; i++)
        {
            cel = cmesh->ElementVec()[i];
            if(!cel || cel->Dimension()!=dim) continue;
            
            //            if(cel->Index()==0){
            //                cornerpoint=1;
            //                break;
            //            }
            
            gel = cel->Reference();
            int nnodes = gel->NNodes();
            
            int found = 0;
            for(int in = 0; in < nnodes; in++)
            {
                if(gel->NodePtr(in)->Id() != nodeAtOriginId){
                    continue;
                }else{
                    cornerpoint = in;
                    found++;
                }
                if(found==1) break;
            }
            if(found==1) break;
        }
        
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        
        TPZElementMatrix ek(cel->Mesh(), TPZElementMatrix::EK);
        TPZElementMatrix ef(cel->Mesh(), TPZElementMatrix::EF);
        
        TPZFNMatrix<1000, STATE> RefMat1, RefMat2, RefMat;
        TPZFNMatrix<1000, STATE> DifMat;
        RefMat1.Zero();
        RefMat2.Zero();
        RefMat.Zero();
        DifMat.Zero();
        //
        //        //Matriz de Referencia
        //        TPZManVector<int,3> order(2,100);
        //        intel->GetIntegrationRule().SetOrder(order);
        //        int npoints1 = intel->GetIntegrationRule().NPoints();
        //        cel->CalcStiff(ek,ef);
        //        RefMat1 = ek.fMat;
        //        RefMat = RefMat1;
        //        //ek.fMat.Print("\n\nEK1 = ", saidaMat, EMathematicaInput);
        
        
        //        TPZIntQuad *QIntRule = new TPZIntQuad(100);
        //        intel->SetIntegrationRule(QIntRule);
        //        intel->GetIntegrationRule().SetType(2, 100);
        //        npoints = intel->GetIntegrationRule().NPoints();
        //        cel->CalcStiff(ek,ef);
        //        RefMat2 = ek.fMat;
        //        ek.fMat.Print("\n\nEK2 = ", saidaMat, EMathematicaInput);
        //
        //        REAL erro = Norm(RefMat1-RefMat2);
        //        saidaMat <<"\nErro =  "<< erro <<";"<<std::endl;
        
        
        //matriz calculada com a regra QP com ordem 100
        TPZIntQuadQuarterPoint *QPIntRule = new TPZIntQuadQuarterPoint(100);
        QPIntRule->SetCorner(cornerpoint);
        {
            TPZCompEl *cel = intel;
            cel->SetIntegrationRule((TPZIntPoints *)QPIntRule);
        }
        int npoints1 = intel->GetIntegrationRule().NPoints();
        cel->CalcStiff(ek,ef);
        RefMat2 = ek.fMat;
        RefMat = RefMat2;
        //                ek.fMat.Print("\n\nEKref = ", saidaMat, EMathematicaInput);
        //                erro = Norm(RefMat1-RefMat2);
        //                saidaMat <<"\nErro =  "<< erro <<";"<<std::endl;
        
        if(intQuarterPoint)
        {
            int order =1;
            TPZIntQuadQuarterPoint *QPIntRule = new TPZIntQuadQuarterPoint(order);
            QPIntRule->SetCorner(cornerpoint);
            {
                TPZCompEl *cel = intel;
                cel->SetIntegrationRule(QPIntRule);
            }
        }else{
            TPZIntQuad *QIntRule = new TPZIntQuad(1);
            {
                TPZCompEl *cel = intel;
                cel->SetIntegrationRule(QIntRule);
            }
            intel->GetIntegrationRule().SetType(0, 1);
        }
        
        if(intQuarterPoint){
            saidaMat <<"\n(*Usando QP Rule*) "<<std::endl;
            
        }
        else{
            saidaMat <<"\n(*Usando QL* Rule*) "<<std::endl;
        }
        
        saidaMat <<"\n(*Para p = "<< p <<"*)"<<std::endl;
        //RefMat.Print("\n\nEKRef = ", saidaMat, EMathematicaInput);
        saidaMat <<"\n(*Numero de Pontos da Matriz de referencia: "<<npoints1<<"*)"<<std::endl;
        saidaMat << "VecErros = {};"<<std::endl;
        saidaMat << "VecErrosRelat = {};"<<std::endl;
        saidaMat << "NumPontosInt = {};"<<std::endl;
        
        int pin = 10;
        int pmax = 200;
        if(intQuarterPoint) pmax = 48;
        for(int ip=pin; ip<pmax; ip+=2)
        {
            TPZManVector<int,3> order(2,ip);
            intel->GetIntegrationRule().SetOrder(order);
            
            int npoints = intel->GetIntegrationRule().NPoints();
            
            //load the matrix ek and vector ef of the element
            cel->CalcStiff(ek,ef);
            
            //            ek.fMat.Print("\n\nEKcalc = ", saidaMat, EMathematicaInput);
            
            //if(ip==pin) TempMat.Resize(ek.fMat.Rows(), ek.fMat.Cols()); b
            DifMat = ek.fMat - RefMat;
            
            REAL Erro;
            int indexi, indexj;
            NormMax(DifMat, Erro, indexi, indexj);
            //REAL norma = NormMax(RefMat);
            //REAL ErroRelat = Erro/norma;
            
            saidaMat<< "\nNumPontos = "<<npoints<<";"<<std::endl;
            saidaMat<< "OrdemNecessaria = "<<ip<<";"<<std::endl;
            //            ek.fMat.Print("\n\nEK = ", saidaMat, EMathematicaInput);
            DifMat.Print("\n\nDifMat = ", saidaMat, EMathematicaInput);
            saidaMat << "Erro = Norm[DifMat,Infinity];"<<std::endl;
            //           saidaMat << "Erro = " << Erro<<"; (*Ocorrido na Posicao (i,j) = ("<<indexi<<","<<indexj<<") da Matriz*)"<<std::endl;
            //saidaMat << "\nErroRelativo = " << ErroRelat<<std::endl;
            saidaMat << "AppendTo[VecErros, Erro];" << std::endl;
            //saidaMat << "\nAppendTo[VecErrosRelat, ErroRelat];" << std::endl;
            saidaMat << "AppendTo[NumPontosInt, NumPontos];" << std::endl;
        }
        saidaMat <<"\n\n(*------------------------------------------------ *)"<<std::endl;
        cmesh->CleanUp();
        delete cmesh;
        delete gmesh;
    }
    
}

void NormMax(TPZFMatrix<STATE> A, REAL &val, int &indexi, int &indexj){
    
    val = 0.;
    indexi = 1;
    indexj = 1;
    
    REAL maxAij = 0.;
    REAL temp = 0.;
    
    int i, j;
    int nr = A.Rows();
    int nc = A.Cols();
    
    for(i=0; i<nr; i++){
        for(j=0; j<nc; j++){
            
            maxAij = val;
            temp = fabs(A(i,j));
            val = Max(temp, maxAij);
            
            REAL dif = fabs(maxAij - val);
            
            if(dif > 1.e-15){
                indexi = i;
                indexj = j;
            }
        }
    }
}

void ChangeIntegrationRule(TPZCompMesh *cmesh, int porder, bool IsQPrule){
    
    
    //quarter point proximo do no de id=1
    int nodeAtOriginId = 1;
    
    int nelem = cmesh->NElements();
    int dim  = cmesh->Dimension();
    
    TPZCompEl *cel = NULL;
    TPZGeoEl * gel = NULL;
    int found = 0;
    int cornerpoint = -1;
    for(int i=0; i<nelem; i++)
    {
        cel = cmesh->ElementVec()[i];
        if(!cel || cel->Dimension()!=dim) continue;
        
        gel = cel->Reference();
        int nnodes = gel->NNodes();
        
        
        for(int in = 0; in < nnodes; in++)
        {
            if(gel->NodePtr(in)->Id() != nodeAtOriginId) continue;
            
            if(IsQPrule){
                cornerpoint = in;
                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                TPZIntQuadQuarterPoint *QPIntRule = new TPZIntQuadQuarterPoint(porder);
                QPIntRule->SetCorner(cornerpoint);
                {
                    TPZCompEl *cel = intel;
                    cel->SetIntegrationRule(QPIntRule);
                }
                
                TPZManVector<int,3> order(2, porder);
                intel->GetIntegrationRule().SetOrder(order);
                found++;
            }
            else{
                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                TPZIntQuad *GIntRule = new TPZIntQuad(1);
                //intel->F(GIntRule);
                intel->GetIntegrationRule().SetType(0, 1);
                
                TPZManVector<int,3> order(2, porder);
                intel->GetIntegrationRule().SetOrder(order);
                found++;
            }
            if(found!=0) break;
        }
        if(found==2) break;
    }
    if(found!=2) DebugStop();
}

void TransferMatrixFromMeshes(TPZCompMesh *cmesh, TPZCompMesh *MPMesh, TPZAutoPointer< TPZMatrix<STATE> > matF,TPZAutoPointer< TPZMatrix<STATE> > matMP, int nodeAtOriginId)
{
    int dim = cmesh->Reference()->Dimension();
    
    TPZBlock<STATE> blockMP = MPMesh->Block();
    TPZBlock<STATE> blockF = cmesh->Block();
    
    
    blockF.SetMatrix(matF.operator->());
    blockMP.SetMatrix(matMP.operator->());
    
    long nel = cmesh->NElements();
    long nelMP = MPMesh->NElements();
    if(nel!=nelMP) DebugStop();
    int found = 0;
    
    long connectindex = 0, connectindexMP = 0;
    long seqnumi=0, seqnumMPi=0, seqnumj=0, seqnumMPj=0;
    
    for(int iel = 0; iel<nel; iel++){
        
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        TPZCompEl * celMP = MPMesh->ElementVec()[iel];
        if(cel->Index()!=celMP->Index()) DebugStop();
        
        TPZGeoEl * gel = cel->Reference();
        if(!gel) continue;
        if(gel->HasSubElement()) continue;
        if(gel->Dimension()!=dim) continue;
        int nnodes = gel->NNodes();
        
        for(int in = 0; in < nnodes; in++){
            
            int nodeid = gel->NodePtr(in)->Id();
            if(nodeid == nodeAtOriginId){
                
                found++;
                int ncon = cel->NConnects();
                int ic, jc;
                
                for(ic=0; ic<ncon; ic++){
                    
                    connectindex = cel->ConnectIndex(ic);
                    connectindexMP = celMP->ConnectIndex(ic);
                    // if(connectindex!=connectindexMP) DebugStop();
                    
                    TPZConnect &coni = cmesh->ConnectVec()[connectindex];
                    seqnumi = coni.SequenceNumber();
                    TPZConnect &conMPi = MPMesh->ConnectVec()[connectindexMP];
                    seqnumMPi = conMPi.SequenceNumber();
                    
                    if (coni.NDof() != conMPi.NDof()) {
                        DebugStop();
                    }
                    
                    for(jc=0; jc<ncon; jc++){
                        
                        connectindex = cel->ConnectIndex(jc);
                        connectindexMP = celMP->ConnectIndex(jc);
                        
                        TPZConnect &conj = cmesh->ConnectVec()[connectindex];
                        seqnumj = conj.SequenceNumber();
                        TPZConnect &conMPj =  MPMesh->ConnectVec()[connectindexMP];
                        seqnumMPj = conMPj.SequenceNumber();
                        if (conj.NDof() != conMPj.NDof()) {
                            DebugStop();
                        }
                        
                        TPZFMatrix<STATE> blocktemp;
                        blockF.GetBlock(seqnumi, seqnumj, blocktemp);
                        blockMP.PutBlock(seqnumMPi, seqnumMPj, blocktemp);
                    }
                }
                break;
            }//end if
        }//end for nnnode
        if(found==2) break;
    }
    if(found!=2) DebugStop();
}

void GlobalSubMatrix(TPZCompMesh *cmesh, TPZAutoPointer< TPZMatrix<STATE> > mat, int nodeAtOriginId, bool matInicial, std::ofstream &subMat)
{
    int dim = cmesh->Reference()->Dimension();
    TPZBlock<STATE> blockMat = cmesh->Block();
    blockMat.SetMatrix(mat.operator->());
    
    long nel = cmesh->NElements();
    
    long connectindex = 0;
    long seqnumi=0, seqnumj=0;
    
    TPZFNMatrix<500,STATE> substiff;
    
    for(int iel = 0; iel<nel; iel++){
        
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        
        TPZGeoEl * gel = cel->Reference();
        if(!gel) continue;
        if(gel->HasSubElement()) continue;
        if(gel->Dimension()!=dim) continue;
        if(gel->NodePtr(1)->Id() != nodeAtOriginId) continue;
        
        
        //int ncon = cel->NConnects();
        int ic, jc;
        
        // calcular o tamanho
        // redimensionar
        
        int icount = 0;
        long row=0,col=0;
        for(ic=0; ic<5; ic++){
            
            connectindex = cel->ConnectIndex(ic);
            
            TPZConnect &coni = cmesh->ConnectVec()[connectindex];
            seqnumi = coni.SequenceNumber();
            long iblsize = blockMat.Size(seqnumi);
            // int iblpos = blockMat.Position(seqnumi);
            if (iblsize != coni.NDof()) {
                DebugStop();
            }
            
            row +=iblsize;
            int jcount = 0;
            
            for(jc=0; jc<5; jc++){
                
                connectindex = cel->ConnectIndex(jc);
                
                TPZConnect &conj = cmesh->ConnectVec()[connectindex];
                seqnumj = conj.SequenceNumber();
                
                long jblsize = blockMat.Size(seqnumj);
                if (conj.NDof() != jblsize) {
                    DebugStop();
                }
                
                if(ic==0) col +=jblsize;
                substiff.Resize(row,col);
                
                for (int r = 0; r < iblsize; r++) {
                    for (int c = 0; c < jblsize; c++) {
                        substiff.PutVal(icount+r, jcount+c, blockMat.GetVal(seqnumi, seqnumj,r,c));
                    }
                }
                jcount += conj.NDof();
            }
            icount += coni.NDof();
        }
        break;
    }
    
    if(matInicial){
        subMat<<"NRows = " <<substiff.Rows()<<";"<<std::endl;
        subMat<<"NCols = " <<substiff.Cols()<<";"<<std::endl;
        substiff.Print("MatAntesTransfer = ", subMat,EMathematicaInput);
        subMat<<"(* ---fim  ---*)"<<std::endl;
    }else{
        subMat<<"\n\nNRows = " <<substiff.Rows()<<";"<<std::endl;
        subMat<<"NCols = " <<substiff.Cols()<<";"<<std::endl;
        substiff.Print("MatAposTransfer = ", subMat,EMathematicaInput);
        subMat<<"\n\nMatrixForm[MatAntesTransfer-MatAposTransfer]"<<std::endl;
        //subMat<<"(* ---fim  ---*)"<<std::endl;
    }
}

void ComputeCharacteristicHElSize(TPZGeoMesh * geometry, REAL & aspectratio/*h_min*/, REAL & h_max){
    
    
    h_max=0;
    
    REAL h;
    REAL rho;
    int nel = geometry->NElements();
    TPZVec<REAL> ratio(nel,0);
    TPZVec<REAL> vecH(nel,0);
    TPZVec<REAL> vecRho(nel,0);
    
    int count=0;
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if(!gel){
            DebugStop();
        }
#endif
        if (gel->Dimension() != geometry->Dimension() || gel->HasSubElement() == 1) {
            continue;
        }
        
        h = gel->CharacteristicSize();
        
        //     std::cout << "h= " <<h<<std::endl;
        rho = 2.0*(gel->ElementRadius());
        
        REAL auxratio=h/rho;
        vecH[count]=h;
        ratio[count]=auxratio;
        count ++;
        
    }
    
    int nvec= ratio.NElements();
    
    for(int j=0;j< nvec;j++){
        if(ratio[j]> aspectratio){
            aspectratio=ratio[j];
        }
    }
    
    for(int j=0;j< nvec;j++){
        if( vecH[j]>h_max){
            h_max=vecH[j];
            //  std::cout << "h= " <<h_max<< std::endl;
            
        }
    }
    
    std::cout << "h= " <<h_max<<" apectratio = "    << aspectratio << std::endl;
    //
    //    std::cout << " rhomin = "    << rho_min << std::endl;
    
}





TPZGeoMesh *CurvedMesh()
{
    int dim=2;
    
    REAL coord[23][2] =
    {
        {0,0},
        {0,0},
        {0,0},
        {0,0},
        {0,0},
        {0,0},
        {0,0},
        {0,0},
        {0,0},
        {-1./8,0},
        {-1./8,1./8},
        {0,1./8},
        {1./8,1./8},
        {1./8,0},
        {-0.5,0},
        {-0.5,0.25},
        {-0.5,0.5},
        {-0.25,0.5},
        {0,0.5},
        {0.25,0.5},
        {0.5,0.5},
        {0.5,0.25},
        {0.5,0}
        
    };
    
    int elno[4][8]=
    {
        {0,2,16,14,1,10,15,9},
        {2,4,18,16,3,11,17,10},
        {4,6,20,18,5,12,19,11},
        {6,8,22,20,7,13,21,12}
    };
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nodes =  23;
    REAL radius = 0.5;
    //REAL innerradius = radius/2.0;
    geomesh->NodeVec().Resize(nodes);
    geomesh->SetDimension(dim);
    
    for(int no=0; no<nodes; no++)
    {
        TPZManVector<REAL,3> xco(3,0.);
        xco[0] = coord[no][0];
        xco[1] = coord[no][1];
        geomesh->NodeVec()[no].Initialize(xco, *geomesh);
    }
    
    for(int el=0; el<4; el++)
    {
        //TPZManVector<long,8> elnodes(8);
        TPZVec<int64_t>elnodes(8);
        for (int i=0; i<8; i++) {
            elnodes[i] = elno[el][i];
        }
        TPZGeoEl *gel = new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad > (elnodes, matId,*geomesh);
    }
    geomesh->BuildConnectivity();
    for(int el=0; el<4; el++)
    {
        TPZGeoEl *gel= geomesh->Element(el);
        TPZGeoElBC(gel,4,-6);
        switch (el) {
            case 0:
                TPZGeoElBC(gel,6,-5);
                TPZGeoElBC(gel,7,-1);
                break;
            case 1:
                TPZGeoElBC(gel,6,-4);
                break;
            case 2:
                TPZGeoElBC(gel,6,-4);
                break;
            case 3:
                TPZGeoElBC(gel,5,-2);
                TPZGeoElBC(gel,6,-3);
                break;
            default:
                break;
        }
    }
    
    geomesh->BuildConnectivity();
    
    
    return geomesh;
}

TPZGeoMesh *VerticalExtrusion(TPZGeoMesh *malha2D ){
    
    int n = 1;
    REAL t0 = 0.0;
    REAL tf = 0.5;
    REAL dt = (tf-t0)/REAL(n);
    int bottom_id = -7;
    int top_id = -8;
    TPZHierarquicalGrid grid;
    grid.SetGeometricMesh(malha2D);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(ParametricfunctionZ, 5);
    grid.SetParametricFunction(ParFunc);
    grid.SetFrontBackMatId(bottom_id, top_id);
    TPZGeoMesh * geometry3D = grid.ComputeExtrusion(t0, dt, n);
    geometry3D->SetName("L_shape deformed");
    geometry3D->SetDimension(3);
    return geometry3D;
    
}

void ParametricfunctionZ(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}

long ComputeAverageBandWidth(TPZCompMesh *cmesh){
    
    TPZVec<int64_t> skyline;
    cmesh->Skyline(skyline);
    int nelemSkyline = 0;
    long nel = skyline.NElements();
    for(long i=0; i<nel; i++){
        nelemSkyline += i-skyline[i]+1;
        if(nelemSkyline < 0) std::cout << "negativo\n";
    }
    long averageband = nelemSkyline/skyline.NElements();
    return averageband;
}


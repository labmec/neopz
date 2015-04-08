//
//  main_mixed.cpp
//  PZ
//
//  Created by Agnaldo Farias on 3/05/15.
//  Copyright (c) 2015 LabMec-Unicamp. All rights reserved.
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
#include "pzmultiphysicselement.h"
#include "pzintel.h"
#include "TPZVTKGeoMesh.h"

#include "TPZMatDualHybridPoisson.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"

#include <iostream>
#include <math.h>
using namespace std;

int const matId =1;
//int const bc1=-1;
//int const bc2=-2;
//int const bc3=-3;
//int const bc4=-4;
//int const bc5=-5;

REAL const Pi = 4.*atan(1.);

//with hybrid method
TPZGeoMesh * MalhaGeo();
TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder,bool ismultiplierH1);
TPZGeoMesh *GMesh();

void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide);
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void Prefinamento(TPZCompMesh * cmesh, int ndiv,int porder);
void PrefinamentoRibsHybridMesh(TPZCompMesh *cmesh);
void SetPOrderRibsHybridMesh(TPZCompMesh *cmesh, int porder);

void GroupElements(TPZCompMesh *cmesh);

//problema Steklov
void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);

//problema Suave
void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &df);
void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

bool fsolsuave=false;
int main(int argc, char *argv[])
{
    //#ifdef LOG4CXX
    //    InitializePZLOG();
    //#endif
    
    ///Refinamento
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
    bool multiplicadorH1 = true;
    bool MudarOrdemPdoMultiplicador = false;
    std::ofstream myerrorfile("Simulacao-Primal-LDGD.txt");
    
    if(multiplicadorH1){
        myerrorfile<<"\nDADOS PARA O REFINAMENTO hp E MULTIPLICADOR CONTINUO"<<std::endl;
        myerrorfile<<std::endl;
        
    }else {
        myerrorfile<<"\nDADOS PARA O REFINAMENTO hp E MULTIPLICADOR DESCONTINUO"<<std::endl;
        myerrorfile<<std::endl;
    }
    
    int pini = 2;
    for(int p = pini; p<3; p++)
    {
        
        myerrorfile<<"\nORDEM p = "<<p <<"\n\n";
        for(int ndiv=9; ndiv<10; ndiv++){
            TPZGeoMesh *gmesh = MalhaGeo();//malha geometrica
            
            if(fsolsuave){
                UniformRefine(gmesh, ndiv);
            }
            
            if(!fsolsuave){
                //refinamento proximo do no de id=1
                DirectionalRef(gmesh, 1, ndiv);
            }
            
            TPZCompMesh *cmesh = CreateHybridCompMesh(*gmesh, p, multiplicadorH1);//malha computacional
            
            if(fsolsuave==false){
                Prefinamento(cmesh, ndiv,p);
                //PrefinamentoRibsHybridMesh(cmesh);
            }
            //            {
            //                                std::ofstream malhaOut("malhageometricaHibrid.vtk");
            //                                 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
            //                                std::ofstream out("cmeshHib1.txt");
            //                                cmesh->Print(out);
            //
            //                                std::ofstream out2("gmesh.txt");
            //                                gmesh->Print(out2);
            //            }
            
            //------- Criar elementos de Lagrange (Ribs)--------
            //materiais do problema
            std::set<int>matids;
            matids.insert(matId);
            if(!fsolsuave){
                for (int i=2; i<=6; i++) {
                    matids.insert(-i);
                }
            }else{
                matids.insert(-1);
            }
            cmesh->ApproxSpace().Hybridize(*cmesh, matids, multiplicadorH1);
            
            if(MudarOrdemPdoMultiplicador){
                SetPOrderRibsHybridMesh(cmesh, 2/*p-(pini-1)*/);
                cmesh->CleanUpUnconnectedNodes();
                cmesh->ExpandSolution();
            }
            
            //            std::ofstream out2("cmeshHib2.txt");
            //            cmesh->Print(out2);
            
            
            myerrorfile << "\nRefinamento h = "<< ndiv <<"\n";
            myerrorfile << "\nDOF Total = "<< cmesh->NEquations() << "\n";
            
            //condesacao estatica
            GroupElements(cmesh);
            cmesh->LoadReferences();//mapeia para a malha geometrica lo
            
            myerrorfile << "DOF Condensado = " << cmesh->NEquations() << "\n";
            
            
            //Resolver problema
            TPZAnalysis analysis(cmesh);
            
            //            TPZSkylineNSymStructMatrix str(cmesh);
            //            str.SetNumThreads(8);
            //            //TPZFStructMatrix str(cmesh);
            //            TPZAutoPointer<TPZMatrix<STATE> > mat = str.Create();
            //            str.EquationFilter().Reset();
            //            TPZAutoPointer<TPZMatrix<STATE> > mat2 = mat->Clone();
            //            analysis.SetStructuralMatrix(str);
            //            TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
            //            TPZStepSolver<STATE> *gmrs = new TPZStepSolver<STATE>(mat2);
            //            step->SetReferenceMatrix(mat2);
            //            step->SetDirect(ELU);
            //            gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
            //            TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
            //            TPZAutoPointer<TPZMatrixSolver<STATE> > autogmres = gmrs;
            //            analysis.SetSolver(autogmres);
            
            
            TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh);
            strmat.SetDecomposeType(ELDLt);
            strmat.SetNumThreads(16);
            analysis.SetStructuralMatrix(strmat);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt); //caso simetrico
            //	step.SetDirect(ELU);
            analysis.SetSolver(step);
            
            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            analysis.Assemble();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            analysis.Solve();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
#endif
            
            //            std::ofstream out("cmeshHib22.txt");
            //            cmesh->Print(out);
#ifdef USING_BOOST
            myerrorfile <<"Tempo para Assemblagem:  " << t2-t1 << std::endl;
            myerrorfile <<"Tempo para Resolução:  " << t3-t2 << std::endl;
#endif
            
            TPZVec<std::string> scalnames(2),vecnames(0);
            scalnames[0] = "Solution";
            scalnames[1] = "POrder";

            std::stringstream name;
            name << "Solution_bima" <<ndiv<< ".vtk";
            std::string paraviewfile(name.str());
            analysis.DefineGraphMesh(2,scalnames,vecnames,paraviewfile);
            analysis.PostProcess(0);
            
            
            if(fsolsuave){
                analysis.SetExact(SolSuave);
            }else{
                analysis.SetExact(SolExataSteklov);
            }
            TPZVec<REAL> erros(3);
            analysis.PostProcessError(erros);
            myerrorfile << "Errro H1 = " << erros[0];
            myerrorfile << "   Errro L2 = " << erros[1];
            myerrorfile << "   Errro semi H1 = " << erros[2] << "\n"<<std::endl;
            
            
            cmesh->CleanUp();
            delete cmesh;
            delete gmesh;
        }
        myerrorfile <<"======================================"<<std::endl;
    }
    
    
    return 0;
}

TPZGeoMesh *GMesh(){
    
    int dim = 2;
    int Qnodes = 4;
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    gmesh->SetDimension(dim);
    
    TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolLine(2);
    
    //indice dos nos
    long id = 0;
    REAL valx, dx = 1.;
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = xi*dx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,0.);//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = 1. - xi*dx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx);//coord X
        Node[id].SetCoord(1 ,dx);//coord Y
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
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
    id++;
    
    TopolLine[0] = 3;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
    
    gmesh->BuildConnectivity();
    
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


//refinamento uniforme em direcao ao no
void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide){
    
    ///Refinamento
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
    for (int idivide = 0; idivide < divide; idivide++){
        const int nels = gmesh->NElements();
        TPZVec< TPZGeoEl * > allEls(nels);
        for(int iel = 0; iel < nels; iel++){
            allEls[iel] = gmesh->ElementVec()[iel];
        }
        
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
            TPZVec<TPZGeoEl*> filhos;
            gel->Divide(filhos);
            
        }///for iel
    }//idivide
    
    gmesh->BuildConnectivity();
}///void

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
        if(pOrder1D<0) continue;//elemento de contorno
        int neworder = pOrder1D > porder2D ? pOrder1D : porder2D;
        Rib->PRefine(neworder);
    }
    
}

void SetPOrderRibsHybridMesh(TPZCompMesh *cmesh, int porder){
    const int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZCompEl * cel = cmesh->ElementVec()[iel];
        if(!cel)continue;
        TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
        if(!face)continue;
        //TPZInterpolationSpace *TwoD = NULL;
        TPZInterpolationSpace *Rib = NULL;
        if(face->LeftElement()->Reference()->Dimension()==2){
            //TwoD = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
            Rib = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
        }
        else{
            Rib = dynamic_cast<TPZInterpolationSpace*>(face->LeftElement());
            //TwoD = dynamic_cast<TPZInterpolationSpace*>(face->RightElement());
        }
        if(!Rib) DebugStop();
        Rib->PRefine(porder);
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
        if(gel->Dimension()==2)
            sp->PRefine(porder + (ndiv - 1) + (level-1));
    }
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
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
    
    du(0,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    
}

void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,result,du);
}

void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {-1,0};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {+1,0};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {0,+1};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFNMatrix<10,REAL> fake(2,1);
    result[0] = loc[0]*loc[0];
}

//----------------------------------
TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder, bool ismultiplierH1){
    //TPZCompEl::SetgOrder(porder);
    TPZCompMesh *comp = new TPZCompMesh(&gmesh);
    comp->SetDimModel(gmesh.Dimension());
    comp->SetDefaultOrder(porder);
    
    comp->ApproxSpace().CreateDisconnectedElements(true);
    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    
    // Criar e inserir os materiais na malha
    REAL beta = 6;
    TPZMatDualHybridPoisson *mymaterial = new TPZMatDualHybridPoisson(matId,0.,beta);
    TPZMaterial * automat(mymaterial);
    comp->InsertMaterialObject(automat);
    
    //Condicoes de contorno
    TPZFMatrix<STATE> val1(1,1,1.),val2(1,1,0.);
    
    if(!fsolsuave){
        TPZMaterial *bnd = automat->CreateBC (automat,-1,2,val1,val2);//misto tbem
        bnd->SetForcingFunction(Dirichlet2);
        comp->InsertMaterialObject(bnd);
        
        // Mixed
        bnd = automat->CreateBC (automat,-2,0,val1,val2);
        bnd->SetForcingFunction(Dirichlet);
        comp->InsertMaterialObject(bnd);
        
        // Mixed
        val1(0,0)=0.;
        val2(0,0)=0.;
        bnd = automat->CreateBC (automat,-3,0,val1,val2);//Dirichlet nulo
        comp->InsertMaterialObject(bnd);
        
        // Mixed
        bnd = automat->CreateBC (automat,-4,0,val1,val2);
        bnd->SetForcingFunction(Dirichlet);
        comp->InsertMaterialObject(bnd);
        
        // Mixed
        bnd = automat->CreateBC (automat,-5,0,val1,val2);
        bnd->SetForcingFunction(Dirichlet);
        comp->InsertMaterialObject(bnd);
        
        // Mixed
        bnd = automat->CreateBC (automat,-6,0,val1,val2);
        bnd->SetForcingFunction(Dirichlet);
        comp->InsertMaterialObject(bnd);
    }
    else{
        
        TPZAutoPointer<TPZFunction<STATE> > forcefunction;
        TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingF);
        dum->SetPolynomialOrder(10);
        forcefunction = dum;
        mymaterial->SetForcingFunction(forcefunction);
        
        TPZMaterial *bnd = automat->CreateBC (automat,-1,0,val1,val2);//Dirichlet Nulo
        comp->InsertMaterialObject(bnd);
    }
    
    // Ajuste da estrutura de dados computacional
    comp->ApproxSpace().CreateDisconnectedElements(true);
    comp->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    std::set<int> matids;
    matids.insert(1);
    comp->AutoBuild(matids);
    
    comp->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    matids.clear();
    if(!fsolsuave){
        for (int i=2; i<=6; i++) {
            matids.insert(-i);
        }
    }else{
        matids.insert(-1);
    }
    
    comp->SetDefaultOrder(porder);
    comp->AutoBuild(matids);
    comp->SetDimModel(2);
    
    comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
    comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
    comp->LoadReferences();
    comp->ApproxSpace().CreateInterfaceElements(comp,true);
    
    //
    //    matids.insert(1);
    //    if(!fsolsuave){
    //        for (int i=2; i<=6; i++) {
    //            matids.insert(-i);
    //        }
    //    }else{
    //        matids.insert(-1);
    //    }
    //
    //    comp->ApproxSpace().Hybridize(*comp, matids, ismultiplierH1);
    
    return comp;
}

TPZGeoMesh * MalhaGeo(){//malha quadrilatera
    
    if(fsolsuave){
        TPZGeoMesh * gmesh = GMesh();
        return gmesh;
    }
    
    ///malha geometrica
    TPZGeoMesh * gmesh = new TPZGeoMesh();
    
    ///Criando nós
    const int nnodes = 6;
    double coord[nnodes][2] = {{-0.5,0},{0,0},{0.,0.5},{-0.5,0.5},{0.5,0},{0.5,0.5}};
    for(int i = 0; i < nnodes; i++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL,3> nodeCoord(3);
        nodeCoord[0] = coord[i][0];
        nodeCoord[1] = coord[i][1];
        nodeCoord[2] = 0.;
        gmesh->NodeVec()[nodind].Initialize(i,nodeCoord,*gmesh);
    }
    
    ///Criando elementos
    const int nel = 2;
    int els[nel][4] = {{0,1,2,3},{1,4,5,2}};
    for(int iel = 0; iel < nel; iel++){
        TPZManVector<long,4> nodind(4);
        long index;
        nodind[0] = els[iel][0];
        nodind[1] = els[iel][1];
        nodind[2] = els[iel][2];
        nodind[3] = els[iel][3];
        gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    ///Criando elementos de contorno
    const int nelbc = 6;
    int bcels[nelbc][3] = {{0,1,-3},{1,4,-2},{4,5,-4},{5,2,-6},{2,3,-6},{3,0,-5}};
    for(int iel = 0; iel < nelbc; iel++){
        TPZManVector<long,4> nodind(2);
        long index;
        nodind[0] = bcels[iel][0];
        nodind[1] = bcels[iel][1];
        int matid = bcels[iel][2];
        gmesh->CreateGeoElement(EOned,nodind,matid,index);
    }
    
    ///Construindo conectividade da malha
    gmesh->BuildConnectivity();
    
    ///Refinamento uniforme da malha
    
    ///Inicializando padrões de refinamento uniforme
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    int numel = gmesh->NElements();
    for (int iel = 0; iel < numel; iel++){
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if (!gel) continue;
        if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
        TPZVec<TPZGeoEl*> filhos;
        gel->Divide(filhos);
    }
    
    return gmesh;
}

void GroupElements(TPZCompMesh *cmesh)
{
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
        
        //#ifdef LOG4CXX
        //        if (logger->isDebugEnabled())
        //        {
        //            std::stringstream sout;
        //            for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
        //                (*it)->Print(sout);
        //            }
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //#endif
        
        long index;
        TPZElementGroup *celgroup = new TPZElementGroup(*cmesh,index);
        elgroupindices.insert(index);
        for (std::list<TPZCompEl *>::iterator it = group.begin(); it != group.end(); it++) {
            celgroup->AddElement(*it);
        }
    }
    cmesh->ComputeNodElCon();
    
    for (std::set<int>::iterator it = elgroupindices.begin(); it!=elgroupindices.end(); it++) {
        TPZCompEl *cel = cmesh->ElementVec()[*it];
        TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel);
    }
    
    cmesh->CleanUpUnconnectedNodes();
}


//solucao suave:sin*sin
void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp,  TPZFMatrix<STATE> &df){
    double x = pt[0];
    double y = pt[1];
    disp.Resize(1);
    disp[0]=0.;
    df.Redim(1, 1);
    disp[0]= 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}

void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    u.Resize(1, 0.);
    du.Resize(3, 1.);
    du(0,0)=du(1,0)=du(2,0)=0.;
    
    const REAL sol = sin(M_PI*x)*sin(M_PI*y);
    u[0] = sol;
    
    du(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y);
    du(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x);
}


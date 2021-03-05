////
////  mainHpJan2015.cpp
////  PZ
////
////  Created by Denise de Siqueira on 1/16/15.
////
////
//
//#ifdef HAVE_CONFIG_H
//#include <pz_config.h>
//#endif
//
//#include "pzvec.h"
//#include "pzstack.h"
//#include "pzfmatrix.h"
//#include "pzfstrmatrix.h"
//#include "pzlog.h"
//
//#include "pzgmesh.h"
//#include "pzcmesh.h"
//#include "pzcompel.h"
//#include "TPZInterfaceEl.h"
//#include "pzgeoelside.h"
//#include "TPZGeoLinear.h"
//#include "pzgeopoint.h"
//
//#include "TPZRefPattern.h"
//#include "tpzgeoelrefpattern.h"
//#include "tpzcompmeshreferred.h"
//#include "tpzautopointer.h"
//#include "pzbndcond.h"
//#include "pzanalysis.h"
//
//#include "pzskylstrmatrix.h"
//#include "pzstepsolver.h"
//#include "pzstrmatrix.h"
//#include "TPZFrontNonSym.h"
//#include "TPZFrontSym.h"
//#include "TPBSpStructMatrix.h"
//#include "TPZSpStructMatrix.h"
//#include "pzbstrmatrix.h"
//
//#include "pzpoisson3d.h"
//#include "pzpoisson3dreferred.h"
//#include "TPZVTKGeoMesh.h"
//
//#include "pzmultiphysicselement.h"
//#include "pzmultiphysicscompel.h"
//#include "pzbuildmultiphysicsmesh.h"
//#include "pzcondensedcompel.h"
//#include "pzelementgroup.h"
//#include "TPZLagrangeMultiplier.h"
//#include "pzmatmixedpoisson3d.h"
//
//#include "mixedpoisson.h"
//
//#include "pzlog.h"
//
//#include <iostream>
//#include <string>
//
//#include <math.h>
//#include <cmath>
//#include <set>
//
//#include "TPZParFrontStructMatrix.h"
//
//# include "Tools.h"
//
//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
//#endif
//
//using namespace std;
//
//const int matId = 1;
//const int matskeleton = -5;
//
//const int bcdirichlet = 0;
//const int bcneumann = 1;
//
//const int bc0 = -1;
//const int bc1 = -2;
//const int bc2 = -3;
//const int bc3 = -4;
//bool fTriang = false;
//const int eps=100000;
//const REAL Pi=M_PI;
//TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly,bool ftriang);
//TPZGeoMesh *MalhaGeomTeste();
//
//
//TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh);
//TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh);
//TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh, bool hdivskeleton);
//TPZGeoMesh * MalhaGeo4Element(const int h);
//TPZCompMesh *MalhaCompH1(TPZGeoMesh * gmesh,int ordem);
//
//void TesteMesh(TPZCompMesh *mesh,std::ostream &out);
//
//void UniformRefine2(TPZGeoMesh* gmesh, int nDiv);
//
//void PrintGMeshVTK2(TPZGeoMesh * gmesh, std::ofstream &file);
//
//void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh,int numthreads);
//void SaidaSolucao(TPZAnalysis &an, std::string plotfile);
//
//void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out);
//void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out);
//
//void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
//void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
//
//void NeumannBC1(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannBC3(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//void NeumannBC4(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
//
/////Metodos para o ArcTangente
//void RefineGeoElements2(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined);
//void RegularizeMesh2(TPZGeoMesh *gmesh);
//void GetPointsOnCircunference2(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points);
//void RefiningNearCircunference2(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);
//void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder);
//void SetPDifferentOrder(TPZCompMesh *comp,int porder);
//
//void SolArcTan2(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux);
//void ForcingTang2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
//void ForcingTang3(const TPZVec<REAL> &pt, TPZVec<REAL> &res,TPZFMatrix<STATE> &disp);
//
//int EquationsElementosConectados(TPZCompMesh *cmesh);
//void PrefinamentoTeste(TPZCompMesh * cmesh, int porder);
//void ChangeOrderExternalConnects(TPZCompMesh *mesh);
//
//bool fmetodomisto;
//
//int flevel=0;
//
//int main()
//{
//    //InitializePZLOG();
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    
//    TPZVec<REAL> erros;
//    
//    fmetodomisto = true;//false --> formulacao H1
//    
//    ofstream saidaerrosHdiv("../ErroHP-Misto.txt",ios::app);
//    ofstream saidaerrosH1("../ErroHP-H1.txt",ios::app);
//    
//    
//    int p = 5;
//    int pq = p;
//    int pp = p;
//    
//    flevel = 2;
//    //Apenas quando estou usando o balanceamento antigo: PkPk-1
//    if(fTriang==true && fmetodomisto==true){
//        pq = p+1;
//    }
//    
//    if(fmetodomisto==true){
//        saidaerrosHdiv<< "\nSaida do erro hp-adptativo para formulacão hdiv, com ordem p inicial = " << p << "\n\n";
//    }
//    else{
//        saidaerrosH1<< "\nSaida do erro hp-adptativo para formulacão H1, com ordem p inicial = " << p << "\n\n";
//    }
//    //refinamentos h adptativos
//    for(int ref_adpth = 2; ref_adpth<6; ref_adpth++){
//        
//        TPZGeoMesh * gmesh = MalhaGeom(1,1,fTriang);
//        //TPZGeoMesh * gmesh = MalhaGeomTeste();
//        UniformRefine(gmesh,/*flevel+*/ref_adpth);
//        
////        std::ofstream filemesh("MalhaGeometricaInicial.vtk");
////        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
//        
//        //RefiningNearCircunference2(2,gmesh,ref_adpth,1);
//        
////         std::ofstream filemesh51("MalhaGeoHConsPRef.vtk");
////         TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh51, true);
//        
////        ofstream filemesh1("malhageometrica.txt");
////        gmesh->Print(filemesh1);
//        
//        //rodar formulacao mista
//        if(fmetodomisto==true){
//            
//            // Criando a primeira malha computacional
//            TPZCompMesh * cmesh1= CMeshFlux(pq, gmesh);
//            //Prefinamento(cmesh1, ref_adpth, pq);
//            //PrefinamentoTeste(cmesh1, pq);
//            ChangeOrderExternalConnects(cmesh1);
//            
////            std::ofstream filemesh2("MalhaFluxAposPRef.txt");
////            cmesh1->Print(filemesh2);
//            
//            // Criando a segunda malha computacional
//            TPZCompMesh * cmesh2 = CMeshPressure(pp, gmesh);
//            //Prefinamento(cmesh2, ref_adpth, pp);
//            //PrefinamentoTeste(cmesh2, pp);
//            
////            ofstream filemesh3("MalhaPressaoPosP_Ref.txt");
////            cmesh2->Print(filemesh3);
//            
//            
//            // Criando a malha computacional multifísica
//            //malha multifisica
//            TPZManVector<TPZCompMesh *,2> meshvec(2);
//            meshvec[0] = cmesh1;
//            meshvec[1] = cmesh2;
//            TPZCompMesh * mphysics = MalhaCompMultifisica(meshvec, gmesh,true);
//            
////            ofstream filemesh4("MalhaMultifisicaPosP_Ref.txt");
////            mphysics->Print(filemesh4);
//            
//            //TesteMesh(mphysics,testeOrdem);
//            
//            int nDofTotal;
//            nDofTotal = meshvec[0]->NEquations() + meshvec[1]->NEquations();
//            
//            int nDofCondensed;
//            nDofCondensed = mphysics->NEquations();
//            
//            saidaerrosHdiv<< "\nNRefinamento h adptativo = "<< ref_adpth <<std::endl;
//            saidaerrosHdiv<< "\nGrau de Liberdade Total = " << nDofTotal<<std::endl;
//            saidaerrosHdiv<< "\nGrau de Liberdade Condensado = " << nDofCondensed<<std::endl;
//            
//            // Resolvendo o sistema linear
//            TPZAnalysis an(mphysics);
//            ResolverSistema(an, mphysics,8);
//            
////            ofstream filemesh4("MalhaMultifisicaPosP_Ref.txt");
////            mphysics->Print(filemesh4);
//            
//            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
//            
//            //Arquivo de saida para plotar a solução
////            string plotfile("Solution_ArcTanHP_Hdiv.vtk");
////            SaidaSolucao(an, plotfile);
//            
//            saidaerrosHdiv<<"\n\nErro da simulacao multifisica  para o Fluxo";
//            //TPZAnalysis an1(cmesh1);
//            //an1.SetExact(*SolSuave);
//            //an1.PostProcessError(erros, saidaerros);
//            ErrorHDiv2(cmesh1,saidaerrosHdiv);
//            
//            saidaerrosHdiv<<"\n\nErro da simulacao multifisica  para a Pressao";
//            //            TPZAnalysis an2(cmesh2);
//            //            an2.SetExact(*SolSuave);
//            //            an2.PostProcessError(erros, saidaerros);
//            ErrorH1(cmesh2, saidaerrosHdiv);
//        }
//        else{
//            
//            TPZCompMesh * cmesh= MalhaCompH1(gmesh,p);
//            //Prefinamento(cmesh, ref_adpth, p);
//            
//            ofstream arg("MalhaCompH1.txt");
//            cmesh->Print(arg);
//            
//            int nDofTotal;
//            nDofTotal=cmesh->NEquations();
//            //nDofCondensed = EquationsElementosConectados(cmesh) - nDofTotal;
//            
//            saidaerrosH1<< "\nNRefinamento h adptativo = "<< ref_adpth <<std::endl;
//            saidaerrosH1<< "\nGrau de Liberdade Total = " << nDofTotal<<std::endl;
//            //saidaerrosH1<< "\nGrau de Liberdade Condensado = " << nDofCondensed<<std::endl;
//            
//            //ELementos por gau polinomial
////            int NElk2 = 0;
////            int NElk3 = 0;
////            int NElk4 = 0;
////            int NElk5 = 0;
////            int NElk6 = 0;
////            int NElk7 = 0;
////            
////            int nel = cmesh->NElements();
////            for(int iel = 0; iel < nel; iel++){
////                TPZCompEl *cel = cmesh->ElementVec()[iel];
////                if(!cel) continue;
////                
////                TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
////                if(!sp) continue;
////                TPZGeoEl * gel = sp->Reference();
////                
////                if(gel->Dimension()==2)
////                {
////                    int level = sp->Reference()->Level();
////                    
////                    if(level==flevel) NElk2 +=1;
////                    if(level==flevel+1) NElk3 +=1;
////                    if(level==flevel+2) NElk4 +=1;
////                    if(level==flevel+3) NElk5 +=1;
////                    if(level==flevel+4) NElk6 +=1;
////                    if(level==flevel+5) NElk7 +=1;
////                }
////            }
//
////            saidaerrosH1<< "\nQuantidade de Elementos de Ordem k = 2: "<< NElk2 <<std::endl;
////            saidaerrosH1<< "\nQuantidade de Elementos de Ordem k = 3: "<< NElk3 <<std::endl;
////            saidaerrosH1<< "\nQuantidade de Elementos de Ordem k = 4: "<< NElk4 <<std::endl;
////            saidaerrosH1<< "\nQuantidade de Elementos de Ordem k = 5: "<< NElk5 <<std::endl;
////            saidaerrosH1<< "\nQuantidade de Elementos de Ordem k = 6: "<< NElk6 <<std::endl;
////            saidaerrosH1<< "\nQuantidade de Elementos de Ordem k = 7: "<< NElk7 <<std::endl;
//            
//            
//            // Resolvendo o sistema linear
//            TPZAnalysis an(cmesh);
//            ResolverSistema(an, cmesh,8);
//            
//            //Arquivo de saida para plotar a solução
//            string plotfile("Solution_ArcTanHP_H1.vtk");
//            SaidaSolucao(an, plotfile);
//            
//            //Pos-processamento calculo do erro
//            an.SetExact(*SolArcTan2);
//            an.PostProcessError(erros, saidaerrosH1);
//            saidaerrosH1<<"==========================\n\n";
//            
//        }
//    }
//    return EXIT_SUCCESS;
//}
//
//
//
//TPZGeoMesh *MalhaGeom( REAL Lx, REAL Ly,bool ftriang){
//    
//    int Qnodes = 4;
//    
//    TPZGeoMesh * gmesh = new TPZGeoMesh;
//    gmesh->SetMaxNodeId(Qnodes-1);
//    gmesh->NodeVec().Resize(Qnodes);
//    TPZVec<TPZGeoNode> Node(Qnodes);
//    
//    TPZVec <int64_t> TopolQuad(4);
//    TPZVec <int64_t> TopolTriang(3);
//    TPZVec <int64_t> TopolLine(2);
//    TPZVec <int64_t> TopolPoint(1);
//    
//    //indice dos nos
//    int64_t id = 0;
//    REAL valx;
//    for(int xi = 0; xi < Qnodes/2; xi++)
//    {
//        valx = xi*Lx;
//        Node[id].SetNodeId(id);
//        Node[id].SetCoord(0 ,valx );//coord X
//        Node[id].SetCoord(1 ,0. );//coord Y
//        gmesh->NodeVec()[id] = Node[id];
//        id++;
//    }
//    
//    for(int xi = 0; xi < Qnodes/2; xi++)
//    {
//        valx = Lx - xi*Lx;
//        Node[id].SetNodeId(id);
//        Node[id].SetCoord(0 ,valx );//coord X
//        Node[id].SetCoord(1 ,Ly);//coord Y
//        gmesh->NodeVec()[id] = Node[id];
//        id++;
//    }
//    
//    //indice dos elementos
//    
//    
//    if(ftriang==true)
//    {
//        gmesh->SetMaxNodeId(Qnodes);
//        gmesh->NodeVec().Resize(Qnodes+1);
//        Node.Resize(Qnodes+1);
//        
//        //No central
//        Node[id].SetNodeId(id);
//        Node[id].SetCoord(0, Lx/2.);//coord X
//        Node[id].SetCoord(1, Ly/2.);//coord Y
//        gmesh->NodeVec()[id] = Node[id];
//        
////        id = 0;
////        TopolTriang[0] = 0;
////        TopolTriang[1] = 1;
////        TopolTriang[2] = 3;
////        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
////        id++;
////        
////        TopolTriang[0] = 2;
////        TopolTriang[1] = 1;
////        TopolTriang[2] = 3;
////        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
////        id++;
//
//        id = 0;
//        TopolTriang[0] = 0;
//        TopolTriang[1] = 1;
//        TopolTriang[2] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
//        id++;
//        
//        TopolTriang[0] = 1;
//        TopolTriang[1] = 2;
//        TopolTriang[2] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
//        id++;
//
//        
//        TopolTriang[0] = 2;
//        TopolTriang[1] = 3;
//        TopolTriang[2] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
//        id++;
//        
//        
//        TopolTriang[0] = 3;
//        TopolTriang[1] = 0;
//        TopolTriang[2] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
//        id++;
//        
//        
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
//        id++;
//        
////        TopolLine[0] = 2;
////        TopolLine[1] = 1;
////        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
////        id++;
//
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//        id++;
//        
////        TopolLine[0] = 3;
////        TopolLine[1] = 2;
////        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
////        id++;
//        
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//        id++;
//        
//        TopolLine[0] = 3;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//    }
//    else{
//        TopolQuad[0] = 0;
//        TopolQuad[1] = 1;
//        TopolQuad[2] = 2;
//        TopolQuad[3] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//        id++;
//        
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
//        id++;
//        
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//        id++;
//        
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//        id++;
//        
//        TopolLine[0] = 3;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//    }
//    
//    gmesh->BuildConnectivity();
//    
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout<<"\n\n Malha Geometrica Inicial\n ";
//        gmesh->Print(sout);
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
//    
//    return gmesh;
//}
//
//TPZGeoMesh *MalhaGeomTeste(){
//    
//    int Qnodes = 6;
//    
//    TPZGeoMesh * gmesh = new TPZGeoMesh;
//    gmesh->SetMaxNodeId(Qnodes-1);
//    gmesh->NodeVec().Resize(Qnodes);
//    TPZVec<TPZGeoNode> Node(Qnodes);
//    
//    TPZVec <int64_t> TopolQuad(4);
//    TPZVec <int64_t> TopolLine(2);
//    
//    //indice dos nos
//    int64_t id = 0;
//    REAL valx;
//    REAL dx=0.5;
//    for(int xi = 0; xi < 3; xi++)
//    {
//        valx = xi*dx;
//        Node[id].SetNodeId(id);
//        Node[id].SetCoord(0 ,valx );//coord X
//        Node[id].SetCoord(1 ,0. );//coord Y
//        gmesh->NodeVec()[id] = Node[id];
//        id++;
//    }
//    
//    for(int xi = 0; xi < 3; xi++)
//    {
//        valx = 1. - xi*dx;
//        Node[id].SetNodeId(id);
//        Node[id].SetCoord(0 ,valx );//coord X
//        Node[id].SetCoord(1 ,1.);//coord Y
//        gmesh->NodeVec()[id] = Node[id];
//        id++;
//    }
//    
//    //indice dos elementos
//    id = 0;
//    
//    TopolQuad[0] = 0;
//    TopolQuad[1] = 1;
//    TopolQuad[2] = 4;
//    TopolQuad[3] = 5;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//    id++;
//    
//    TopolQuad[0] = 1;
//    TopolQuad[1] = 2;
//    TopolQuad[2] = 3;
//    TopolQuad[3] = 4;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//    id++;
//    
//    
//    TopolLine[0] = 0;
//    TopolLine[1] = 1;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
//    id++;
//    
//    TopolLine[0] = 1;
//    TopolLine[1] = 2;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
//    id++;
//    
//    TopolLine[0] = 2;
//    TopolLine[1] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//    id++;
//    
//    
//    TopolLine[0] = 3;
//    TopolLine[1] = 4;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//    id++;
//    
//    TopolLine[0] = 4;
//    TopolLine[1] = 5;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//    id++;
//    
//    TopolLine[0] = 5;
//    TopolLine[1] = 0;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//    
//    
//    gmesh->BuildConnectivity();
//    
//    
//    ///Refinamento uniforme
//    TPZVec<TPZGeoEl *> filhos;
//    int64_t n = gmesh->NElements();
//    for ( int64_t i = 0; i < n; i++ )
//    {
//        TPZGeoEl * gel = gmesh->ElementVec() [i];
//        int64_t index = gel->Index();
//        if((index == 0) || (index == 2) || (index == 6) || (index == 7)) continue;
//        if(!gel->HasSubElement())
//        {
//            gel->Divide(filhos);
//        }
//    }//for i
//    
//    gmesh->ResetConnectivities();
//    gmesh->BuildConnectivity();
//
//    
//    
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout<<"\n\n Malha Geometrica Inicial\n ";
//        gmesh->Print(sout);
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
//    
//    return gmesh;
//}
//
//
//
//TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh)
//{
//    /// criar materiais
//    int dim = 2;
//    TPZMatPoisson3d *material;
//    material = new TPZMatPoisson3d(matId,dim);
//    TPZMaterial * mat(material);
//    material->NStateVariables();
//    
//    
//    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//    cmesh->SetDimModel(dim);
//    
//    cmesh->SetAllCreateFunctionsHDiv();
//    
//    cmesh->InsertMaterialObject(mat);
//    
//    ///Criar condicoes de contorno
//    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//    TPZMaterial * BCond0 = material->CreateBC(mat, bc0, bcdirichlet, val1, val2);
//    TPZMaterial * BCond1 = material->CreateBC(mat, bc1, bcdirichlet, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(mat, bc2, bcdirichlet, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(mat, bc3, bcdirichlet, val1, val2);
//    
//    ///Inserir condicoes de contorno
//    cmesh->InsertMaterialObject(BCond0);
//    cmesh->InsertMaterialObject(BCond1);
//    cmesh->InsertMaterialObject(BCond2);
//    cmesh->InsertMaterialObject(BCond3);
//    
//    cmesh->SetDefaultOrder(pOrder);
//    cmesh->SetDimModel(dim);
//    
////    TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(matskeleton, dim-1, 1);
////    TPZMaterial * mat2(matskelet);
////    cmesh->InsertMaterialObject(mat2);
//    
//    cmesh->AutoBuild();//Ajuste da estrutura de dados computacional
//    
//    //#ifdef LOG4CXX
//    //	if(logdata->isDebugEnabled())
//    //	{
//    //        std::stringstream sout;
//    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
//    //        cmesh->Print(sout);
//    //        LOGPZ_DEBUG(logdata,sout.str())
//    //	}
//    //#endif
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    return cmesh;
//}
//
//TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh)
//{
//    /// criar materiais
//    int dim = 2;
//    TPZMatPoisson3d *material;
//    material = new TPZMatPoisson3d(matId,dim);
//    material->NStateVariables();
//    
//    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//    cmesh->SetDimModel(dim);
//    TPZMaterial * mat(material);
//    cmesh->InsertMaterialObject(mat);
//    
//        ///Inserir condicao de contorno
//        TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//    
//    
////        TPZMaterial * BCond0 = material->CreateBC(mat, bc0,bcdirichlet, val1, val2);
////        TPZMaterial * BCond1 = material->CreateBC(mat, bc1,bcdirichlet, val1, val2);
////        TPZMaterial * BCond2 = material->CreateBC(mat, bc2,bcdirichlet, val1, val2);
////        TPZMaterial * BCond3 = material->CreateBC(mat, bc3,bcdirichlet, val1, val2);
////        //para malha com 4 elementos usar mais estes contornos
////        //    TPZMaterial *BCond4 = mat->CreateBC (mat,-5,0,val1,val2);
////        //    TPZMaterial *BCond5 = mat->CreateBC (mat,-6,0,val1,val2);
////        //    TPZMaterial *BCond6 = mat->CreateBC (mat,-7,0,val1,val2);
////        //    TPZMaterial *BCond7 = mat->CreateBC (mat,-8,0,val1,val2);
////    
////    
////    
////    
////        ///Inserir condicoes de contorno
////        cmesh->InsertMaterialObject(BCond0);
////        cmesh->InsertMaterialObject(BCond1);
////        cmesh->InsertMaterialObject(BCond2);
////        cmesh->InsertMaterialObject(BCond3);
////        //    cmesh->InsertMaterialObject(BCond4);
////        //    cmesh->InsertMaterialObject(BCond5);
////        //    cmesh->InsertMaterialObject(BCond6);
////        //    cmesh->InsertMaterialObject(BCond7);
//    
//    
//    cmesh->SetDefaultOrder(pOrder);
//    cmesh->SetDimModel(dim);
//    
//    
//    cmesh->SetAllCreateFunctionsContinuous();
//    cmesh->ApproxSpace().CreateDisconnectedElements(true);
//   //cmesh->SetAllCreateFunctionsDiscontinuous();
//    
//    //Ajuste da estrutura de dados computacional
//    cmesh->AutoBuild();
//    
//    int ncon = cmesh->NConnects();
//    for(int i=0; i<ncon; i++)
//    {
//        TPZConnect &newnod = cmesh->ConnectVec()[i];
//        //newnod.SetPressure(true);
//        newnod.SetLagrangeMultiplier(1);
//    }
//    
//    int nel = cmesh->NElements();
//    for(int i=0; i<nel; i++){
//        TPZCompEl *cel = cmesh->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        if(!celdisc) continue;
//        celdisc->SetConstC(1.);
//        celdisc->SetCenterPoint(0, 0.);
//        celdisc->SetCenterPoint(1, 0.);
//        celdisc->SetCenterPoint(2, 0.);
//        celdisc->SetTrueUseQsiEta();
//        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//        {
//            if(fTriang==true) celdisc->SetTotalOrderShape();
//            else celdisc->SetTensorialShape();
//        }
//        
//    }
//    
//    cmesh->CleanUpUnconnectedNodes();
//    return cmesh;
//}
//
//TPZCompMesh *MalhaCompH1(TPZGeoMesh * gmesh,int ordem){
//    
//    //Creating computational mesh for multiphysic elements
//    gmesh->ResetReference();
//    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
//    
//    //criando material
//    int dim =2;
//    cmesh->SetDefaultOrder(ordem);
//    cmesh->SetDimModel(dim);
//    
//    TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
//    
//    
//    //funcao do lado direito da equacao do problema
//    TPZAutoPointer<TPZFunction<STATE> > force;
//    TPZDummyFunction<STATE> *dum;
//    dum = new TPZDummyFunction<STATE>(ForcingTang3);
//    dum->SetPolynomialOrder(20);
//    force = dum;
//    material->SetForcingFunction(force);
//    
//    //solucao exata
//    TPZAutoPointer<TPZFunction<STATE> > solexata;
//    solexata = new TPZDummyFunction<STATE>(SolArcTan2);
//    material->SetForcingFunctionExact(solexata);
//    
//    //inserindo o material na malha computacional
//    TPZMaterial *mat(material);
//    cmesh->InsertMaterialObject(mat);
//    
//    //Criando condicoes de contorno
//    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//    
//    
//    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,bcdirichlet, val1, val2);
//    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,bcdirichlet, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,bcdirichlet, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,bcdirichlet, val1, val2);
//    
//    
//    ///Inserir condicoes de contorno
//    cmesh->InsertMaterialObject(BCond0);
//    cmesh->InsertMaterialObject(BCond1);
//    cmesh->InsertMaterialObject(BCond2);
//    cmesh->InsertMaterialObject(BCond3);
//    
//    
//    cmesh->SetAllCreateFunctionsContinuous();
//    
//    //Fazendo auto build
//    cmesh->AutoBuild();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    
//    return cmesh;
//}
//
//#include "pzintel.h"
//TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh, bool hdivskeleton){
//    
//    //Creating computational mesh for multiphysic elements
//    gmesh->ResetReference();
//    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
//    
//    //criando material
//    int dim =2;
//    
//    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(matId,dim);
//    //TPZMixedPoisson *material = new TPZMixedPoisson(matId,dim);
//    
//    //incluindo os dados do problema
//    //    REAL coefk = 1.;
//    //    material->SetPermeability(coefk);
//    REAL coefvisc = 1.;
//    material->SetViscosity(coefvisc);
//    
//    //permeabilidade
//    TPZFMatrix<REAL> Ktensor(dim,dim,0.);
//    TPZFMatrix<REAL> InvK(dim,dim,0.);
//    Ktensor(0,0)=1.; Ktensor(1,1)=1.;
//    InvK=Ktensor;
//    material->SetPermeabilityTensor(Ktensor,InvK);
//    
//    //solucao exata
//    TPZAutoPointer<TPZFunction<STATE> > solexata;
//    solexata = new TPZDummyFunction<STATE>(SolArcTan2);
//    material->SetForcingFunctionExact(solexata);
//    
//    //funcao do lado direito da equacao do problema
//    TPZAutoPointer<TPZFunction<STATE> > force;
//    TPZDummyFunction<STATE> *dum;
//    dum = new TPZDummyFunction<STATE>(ForcingTang2);
//    dum->SetPolynomialOrder(20);
//    force = dum;
//    material->SetForcingFunction(force);
//    
//    //inserindo o material na malha computacional
//    TPZMaterial *mat(material);
//    mphysics->InsertMaterialObject(mat);
//    mphysics->SetDimModel(dim);
//    
//    //Criando condicoes de contorno
//    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//    
//    
//    TPZMaterial * BCond0 = material->CreateBC(mat, bc0, bcdirichlet, val1, val2);
//    TPZMaterial * BCond1 = material->CreateBC(mat, bc1, bcdirichlet, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(mat, bc2, bcdirichlet, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(mat, bc3, bcdirichlet, val1, val2);
//    
//    ///Inserir condicoes de contorno
//    mphysics->InsertMaterialObject(BCond0);
//    mphysics->InsertMaterialObject(BCond1);
//    mphysics->InsertMaterialObject(BCond2);
//    mphysics->InsertMaterialObject(BCond3);
//    
//    mphysics->SetAllCreateFunctionsMultiphysicElem();
//    
//    //Fazendo auto build
//    mphysics->AutoBuild();
//    mphysics->AdjustBoundaryElements();
//    mphysics->CleanUpUnconnectedNodes();
//    
//    if(hdivskeleton==false){
//        // Creating multiphysic elements into mphysics computational mesh
//        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
//    }
//    
//    //Creating multiphysic elements containing skeletal elements.
//    else
//    {
//        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//        mphysics->Reference()->ResetReference();
//        mphysics->LoadReferences();
//        
//        int64_t nel = mphysics->ElementVec().NElements();
//        
//        std::map<int64_t, int64_t> bctoel, eltowrap;
//        for (int64_t el=0; el<nel; el++) {
//            TPZCompEl *cel = mphysics->Element(el);
//            TPZGeoEl *gel = cel->Reference();
//            int matid = gel->MaterialId();
//            if (matid < 0) {
//                TPZGeoElSide gelside(gel,gel->NSides()-1);
//                TPZGeoElSide neighbour = gelside.Neighbour();
//                while (neighbour != gelside) {
//                    if (neighbour.Element()->Dimension() == dim && neighbour.Element()->Reference()) {
//                        // got you!!
//                        bctoel[el] = neighbour.Element()->Reference()->Index();
//                        break;
//                    }
//                    neighbour = neighbour.Neighbour();
//                }
//                if (neighbour == gelside) {
//                    DebugStop();
//                }
//            }
//        }
//        
//        TPZStack< TPZStack< TPZMultiphysicsElement *,7> > wrapEl;
//        for(int64_t el = 0; el < nel; el++)
//        {
//            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(mphysics->Element(el));
//            if(mfcel->Dimension()==dim) TPZBuildMultiphysicsMesh::AddWrap(mfcel, matId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
//        }
//        
//        for (int64_t el =0; el < wrapEl.size(); el++) {
//            TPZCompEl *cel = wrapEl[el][0];
//            int64_t index = cel->Index();
//            eltowrap[index] = el;
//        }
//        
//        meshvec[0]->CleanUpUnconnectedNodes();
//        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
//        
//        std::map<int64_t, int64_t>::iterator it;
//        for (it = bctoel.begin(); it != bctoel.end(); it++) {
//            int64_t bcindex = it->first;
//            int64_t elindex = it->second;
//            if (eltowrap.find(elindex) == eltowrap.end()) {
//                DebugStop();
//            }
//            int64_t wrapindex = eltowrap[elindex];
//            TPZCompEl *bcel = mphysics->Element(bcindex);
//            TPZMultiphysicsElement *bcmf = dynamic_cast<TPZMultiphysicsElement *>(bcel);
//            if (!bcmf) {
//                DebugStop();
//            }
//            wrapEl[wrapindex].Push(bcmf);
//
//        }
//        
//        //------- Create and add group elements -------
//        int64_t index, nenvel;
//        nenvel = wrapEl.NElements();
//        TPZStack<TPZElementGroup *> elgroups;
//        for(int64_t ienv=0; ienv<nenvel; ienv++){
//            TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
//            elgroups.Push(elgr);
//            nel = wrapEl[ienv].NElements();
//            for(int jel=0; jel<nel; jel++){
//                elgr->AddElement(wrapEl[ienv][jel]);
//            }
//        }
//        
//        mphysics->ComputeNodElCon();
//        // create condensed elements
//        // increase the NumElConnected of one pressure connects in order to prevent condensation
//        for (int64_t ienv=0; ienv<nenvel; ienv++) {
//            TPZElementGroup *elgr = elgroups[ienv];
//            int nc = elgr->NConnects();
//            for (int ic=0; ic<nc; ic++) {
//                TPZConnect &c = elgr->Connect(ic);
//                if (c.LagrangeMultiplier() > 0) {
//                    c.IncrementElConnected();
//                    break;
//                }
//            }
//            TPZCondensedCompEl *condense = new TPZCondensedCompEl(elgr);
//        }
//    }
//    
//    mphysics->CleanUpUnconnectedNodes();
//    mphysics->ExpandSolution();
//    
//    return mphysics;
//}
//
//
//#define VTK
//void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads)
//{
//    
//    //TPZSkylineStructMatrix strmat(fCmesh);
//    
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
//    strmat.SetDecomposeType(ELDLt);
//    strmat.SetNumThreads(numthreads);
//    
//    an.SetStructuralMatrix(strmat);
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELDLt); //caso simetrico
//    //	step.SetDirect(ELU);
//    an.SetSolver(step);
//    //    an.Assemble();
//    an.Run();
//    
//    //Saida de Dados: solucao e  grafico no VTK
//    //ofstream file("Solution.out");
//    //an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
//}
//
//void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
//    double x = pt[0];
//    double y = pt[1];
//    disp[0]= 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
//}
//
//void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
//    
//    const REAL x = loc[0];
//    const REAL y = loc[1];
//    u.Resize(1, 0.);
//    du.Resize(3, 1.);
//    du(0,0)=du(1,0)=du(2,0)=0.;
//    
//    const REAL sol = sin(M_PI*x)*sin(M_PI*y);
//    u[0] = sol;
//    
//    du(0,0) = -M_PI*cos(M_PI*x)*sin(M_PI*y);//-k*Grad(u)[[0]]
//    du(1,0) = -M_PI*cos(M_PI*y)*sin(M_PI*x);////-k*Grad(u)[[1]]
//    du(2,0) = 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
//}
//
//void NeumannBC1(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    REAL normal[2] = {0,-1.};
//    
//    TPZManVector<REAL> u(1);
//    TPZFNMatrix<10> du(2,1);
//    SolSuave(loc,u,du);
//    
//    result.Resize(1);
//    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
//}
//
//void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    REAL normal[2] = {1.,0};
//    
//    TPZManVector<REAL> u(1);
//    TPZFNMatrix<10> du(2,1);
//    SolSuave(loc,u,du);
//    
//    result.Resize(1);
//    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
//}
//
//void NeumannBC3(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    REAL normal[2] = {0,1.};
//    
//    TPZManVector<REAL> u(1);
//    TPZFNMatrix<10> du(2,1);
//    SolSuave(loc,u,du);
//    
//    result.Resize(1);
//    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
//}
//
//void NeumannBC4(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
//    REAL normal[2] = {-1.,0};
//    
//    TPZManVector<REAL> u(1);
//    TPZFNMatrix<10> du(2,1);
//    SolSuave(loc,u,du);
//    
//    result.Resize(1);
//    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
//}
//
//
//void SaidaSolucao(TPZAnalysis &an, std::string plotfile){
//    
//    if(fmetodomisto==true){
//        TPZManVector<std::string,10> scalnames(3), vecnames(2);
//        scalnames[0] = "Pressure";
//        scalnames[1] = "ExactPressure";
//        scalnames[2]="POrder";
//        vecnames[0]= "Flux";
//        vecnames[1]= "ExactFlux";
//        
//        const int dim = 2;
//        int div = 0;
//        an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//        an.PostProcess(div,dim);
//    }
//    
//    //Rodar H1
//    else{
//        TPZManVector<std::string,10> scalnames(3), vecnames(1);
//        scalnames[0] = "Solution";
//        scalnames[1]="POrder";
//        scalnames[2] = "ExactPressure";
//        vecnames[0]= "Derivative";
//        
//        const int dim = an.Mesh()->Dimension();
//        int div = 0;
//        an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//        an.PostProcess(div,dim);
//    }
//}
//
//
//void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out)
//{
//    int64_t nel = hdivmesh->NElements();
//    int dim = hdivmesh->Dimension();
//    TPZManVector<STATE,10> globerrors(10,0.);
//    TPZStack<REAL> vech;
//    
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = hdivmesh->ElementVec()[el];
//        if (!cel) {
//            continue;
//        }
//        //--metodo para pegar tamanho da malha
//        //        TPZMaterialData data;
//        //        data.fNeedsHSize=true;
//        //        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace * >(cel);
//        //        TPZVec<REAL> qsi(2,0.);
//        //        sp->ComputeRequiredData(data, qsi);
//        //        REAL &hsize = data.HSize;
//        //        vech.push_back(hsize);
//        
//        //----
//        
//        TPZGeoEl *gel = cel->Reference();
//        if (!gel || gel->Dimension() != dim) {
//            continue;
//        }
//        TPZManVector<STATE,10> elerror(10,0.);
//        cel->EvaluateError(SolArcTan2, elerror, NULL);
//        int nerr = elerror.size();
//        for (int i=0; i<nerr; i++) {
//            globerrors[i] += elerror[i]*elerror[i];
//        }
//        
//    }
//    
//    
//    //    int nh = vech.size();
//    //    REAL hmax=0;
//    //    for(int i=0; i<nh; i++){
//    //        if(vech[i]>hmax){
//    //            hmax=vech[i];
//    //        }
//    //    }
//    // out << "Errors associated with HDiv space\n";
//    //out << " Hmax = "    << hmax << std::endl;
//    
//    
//    //out << "Errors associated with HDiv space\n";
//    out << "\nL2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
//    //out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<endl;
//    //out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<endl;
//    
//}
//
//void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out)
//{
//    int64_t nel = l2mesh->NElements();
//    int dim = l2mesh->Dimension();
//    TPZManVector<STATE,10> globerrors(10,0.);
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = l2mesh->ElementVec()[el];
//        if (!cel) {
//            continue;
//        }
//        TPZGeoEl *gel = cel->Reference();
//        if (!gel || gel->Dimension() != dim) {
//            continue;
//        }
//        TPZManVector<STATE,10> elerror(10,0.);
//        elerror.Fill(0.);
//        cel->EvaluateError(SolArcTan2, elerror, NULL);
//        
//        int nerr = elerror.size();
//        globerrors.resize(nerr);
//#ifdef LOG4CXX
//        if (logger->isDebugEnabled()) {
//            std::stringstream sout;
//            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
//            LOGPZ_DEBUG(logger, sout.str())
//        }
//#endif
//        for (int i=0; i<nerr; i++) {
//            globerrors[i] += elerror[i]*elerror[i];
//        }
//    }
//    out << "\n";
//    //out << "Errors associated with L2 or H1 space\n";
//    out << "\nH1 Norm = "    << sqrt(globerrors[0]) << endl;
//    out << "\nL2 Norm = "    << sqrt(globerrors[1]) << endl;
//    out << "\nSemi H1 Norm = "    << sqrt(globerrors[2]) << endl;
//    out << "\n=============================\n"<<endl;
//}
//
//
//void UniformRefine2(TPZGeoMesh* gmesh, int nDiv)
//{
//    for(int D = 0; D < nDiv; D++)
//    {
//        int nels = gmesh->NElements();
//        for(int elem = 0; elem < nels; elem++)
//        {
//            TPZVec< TPZGeoEl * > filhos;
//            TPZGeoEl * gel = gmesh->ElementVec()[elem];
//            gel->Divide(filhos);
//        }
//    }
//    //	gmesh->BuildConnectivity();
//}
/////Metodos para o problema ArcTangente
//void RefineGeoElements2(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
//    
//    TPZManVector<REAL> centerpsi(3), center(3);
//    // Refinamento de elementos selecionados
//    TPZGeoEl *gel;
//    TPZVec<TPZGeoEl *> sub;
//    
//    int nelem = 0;
//    int ngelem=gmesh->NElements();
//    // na esquina inferior esquerda Nó = (0,-1,0)
//    while(nelem<ngelem) {
//        gel = gmesh->ElementVec()[nelem++];
//        if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
//        gel->CenterPoint(gel->NSides()-1,centerpsi);
//        gel->X(centerpsi,center);
//        if(!isdefined) {
//            TPZVec<REAL> FirstNode(3,0.);
//            gel->CenterPoint(0,centerpsi);
//            gel->X(centerpsi,FirstNode);
//            REAL distancia = TPZGeoEl::Distance(center,FirstNode);
//            if(distancia > distance) distance = distancia;
//            isdefined = true;
//        }
//        REAL centerdist = TPZGeoEl::Distance(center,point);
//        if(fabs(r-centerdist) < distance) {
//            gel->Divide(sub);
//        }
//    }
//    
//}
//void GetPointsOnCircunference2(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points) {
//    Points.Resize(npoints);
//    TPZManVector<REAL> point(3,0.);
//    REAL angle = (2*M_PI)/npoints;
//    for(int i=0;i<npoints;i++) {
//        point[0] = center[0]+radius*cos(i*angle);
//        point[1] = center[1]+radius*sin(i*angle);
//        Points[i] = point;
//    }
//}
//void RefiningNearCircunference2(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
//    /*int i;
//     bool isdefined = false;
//     
//     // Refinando no local desejado
//     int npoints = 1000;
//     TPZVec<REAL> point(3);
//     point[0] = point[1] = 0.5; point[2] = 0.0;
//     REAL r = 0.25;
//     TPZVec<TPZManVector<REAL> > Points(npoints);
//     GetPointsOnCircunference2(npoints,point,r,Points);
//     
//     if(ntyperefs==2) {
//     REAL radius = 0.19;
//     for(i=0;i<nref;i+=2) {
//     // To refine elements with center near to points than radius
//     RefineGeoElements2(dim,gmesh,point,r,radius,isdefined);
//     RefineGeoElements2(dim,gmesh,point,r,radius,isdefined);
//     if(i < 5) radius *= 0.35;
//     else if(i < 7) radius *= 0.2;
//     else radius *= 0.1;
//     }
//     if(i==nref) {
//     RefineGeoElements2(dim,gmesh,point,r,radius,isdefined);
//     }
//     }
//     else {
//     REAL radius = 0.2;
//     for(i=0;i<nref+1;i++) {
//     // To refine elements with center near to points than radius
//     RefineGeoElements2(dim,gmesh,point,r,radius,isdefined);
//     if(i < 8) radius *= 0.6;
//     else if(i < 7) radius *= 0.3;
//     else radius *= 0.15;
//     }
//     }
//     // Constructing connectivities
//     //gmesh->ResetConnectivities();
//     RegularizeMesh2(gmesh);
//     gmesh->BuildConnectivity();
//     */
//    int i;
//    bool isdefined = false;
//    
//    // Refinando no local desejado
//    int npoints = 1000;
//    TPZVec<REAL> point(3);
//    point[0] = point[1] = 0.5; point[2] = 0.0;
//    REAL r = 0.25;
//    TPZVec<TPZManVector<REAL> > Points(npoints);
//    GetPointsOnCircunference(npoints,point,r,Points);
//    
//    if(ntyperefs==2) {
//        REAL radius = 0.19;
//        for(i=0;i<nref;i+=2) {
//            // To refine elements with center near to points than radius
//            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
//            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
//            if(nref < 5) radius *= 0.35;
//            else if(nref < 7) radius *= 0.2;
//            else radius *= 0.1;
//        }
//        if(i==nref) {
//            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
//        }
//    }
//    else {
//        REAL radius = 0.03;//0.2;
//        for(i=0;i<nref;i++) {
//            // To refine elements with center near to points than radius
//            RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
//            if(nref < 6) radius *= 0.6;
//            else if(nref < 7) radius *= 0.3;
//            else radius *= 0.15;
//        }
//    }
//    // Constructing connectivities
//    //  gmesh->ResetConnectivities();
//    RegularizeMesh2(gmesh);
//    gmesh->BuildConnectivity();
//    
//    
//}
//
//void RegularizeMesh2(TPZGeoMesh *gmesh)
//{
//    bool changed = true;
//    while (changed)
//    {
//        changed = false;
//        int nel = gmesh->NElements();
//        for (int64_t el=0; el<nel; el++) {
//            TPZGeoEl *gel = gmesh->ElementVec()[el];
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            int dim = gel->Dimension();
//            if (dim != 2) {
//                continue;
//            }
//            int nsides = gel->NSides();
//            int nrefined = 0;
//            int nsidedim = 0;
//            for (int is=0; is<nsides; is++) {
//                int sidedim = gel->SideDimension(is);
//                if (sidedim != dim-1) {
//                    continue;
//                }
//                nsidedim++;
//            }
//            for (int is=0; is<nsides; is++) {
//                int sidedim = gel->SideDimension(is);
//                if (sidedim != dim-1) {
//                    continue;
//                }
//                TPZGeoElSide thisside(gel,is);
//                TPZGeoElSide neighbour = thisside.Neighbour();
//                if (neighbour != thisside) {
//                    TPZStack<TPZGeoElSide> subelements;
//                    neighbour.GetSubElements2(subelements);
//                    int nsub = subelements.size();
//                    if (nsub > 0) {
//                        nrefined++;
//                    }
//                    for (int isub=0; isub<nsub; isub++) {
//                        TPZGeoElSide sub = subelements[isub];
//                        if (sub.Dimension() != dim-1) {
//                            continue;
//                        }
//                        if (sub.HasSubElement()) {
//                            TPZManVector<TPZGeoEl *> newsub;
//                            gel->Divide(newsub);
//                            changed = true;
//                            break;
//                        }
//                    }
//                }
//                if (gel->HasSubElement()) {
//                    break;
//                }
//            }
//            if (nrefined >= nsidedim-1) {
//                TPZManVector<TPZGeoEl *> newsub;
//                gel->Divide(newsub);
//                changed = true;
//            }
//        }
//    }
//}
//
//#define Power pow
//#define ArcTan atan
//#define Sqrt sqrt
//
//void SolArcTan2(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux){
//    REAL x = pt[0];
//    REAL y = pt[1];
//    p[0]=0;
//    flux(0,0)=0;
//    flux(1,0)=0;
//    //flux(2,0)=0;
//    
//    p[0]= 5.*(-1. + x)*x*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))));
//    
//    
//    //px
//    flux(0,0)=(12.732395447351628*Sqrt(eps)*(-1. + x)*(-0.5 + x)*x*(-1. + y)*y)/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) -
//    5.*(-1. + x)*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))) -
//    5.*x*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))));
//    
//    
//    //py
//    flux(1,0)= (12.732395447351628*Sqrt(eps)*(-1. + x)*x*(-1. + y)*(-0.5 + y)*y)/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) -
//    5.*(-1. + x)*x*(-1. + y)*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))) -
//    5.*(-1. + x)*x*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))));
//    
//    
//    if(fmetodomisto==false){
//        flux(0,0) *=-1.;
//        flux(1,0) *=-1.;
//    }
//}
//
//
//void ForcingTang2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
//    double x = pt[0];
//    double y = pt[1];
//    
//    disp[0] = 0.;
//    
//    disp[0]= (-1.)*((3.183098861837907*Sqrt(eps)*(-1. + x)*x*(-4.*(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) -
//                                                              64.*eps*Power(-0.5 + x,2)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))*(-1. + y)*y)/
//                    Power(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2),2) -
//                    (5.092958178940651*Sqrt(eps)*(-0.5 + x)*(-5. + 10.*x)*(-1. + y)*y)/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) +
//                    10.*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))) +
//                    5.*(-1. + x)*x*(2. + (0.6366197723675814*Sqrt(eps)*(-4.*(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) -
//                                                                        64.*eps*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))*Power(-0.5 + y,2))*(-1. + y)*y)/
//                                    Power(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2),2) -
//                                    (5.092958178940651*Sqrt(eps)*(-0.5 + y)*(-1. + 2*y))/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) +
//                                    1.2732395447351628*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))));
//    
//}
//
//void ForcingTang3(const TPZVec<REAL> &pt, TPZVec<REAL> &res,TPZFMatrix<STATE> &disp){
//    //   disp.Redim(2,1);
//    double x = pt[0];
//    double y = pt[1];
//    disp.Resize(2, 2);
//    
//    res[0] = 0.;
//    res[0]= ((3.183098861837907*Sqrt(eps)*(-1. + x)*x*(-4.*(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) -
//                                                       64.*eps*Power(-0.5 + x,2)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))*(-1. + y)*y)/
//             Power(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2),2) -
//             (5.092958178940651*Sqrt(eps)*(-0.5 + x)*(-5. + 10.*x)*(-1. + y)*y)/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) +
//             10.*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))) +
//             5.*(-1. + x)*x*(2. + (0.6366197723675814*Sqrt(eps)*(-4.*(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) -
//                                                                 64.*eps*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))*Power(-0.5 + y,2))*(-1. + y)*y)/
//                             Power(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2),2) -
//                             (5.092958178940651*Sqrt(eps)*(-0.5 + y)*(-1. + 2*y))/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) +
//                             1.2732395447351628*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))));
//}
//
////void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder){
////    if(ndiv<1) return;
////    int nel = cmesh->NElements();
////    for(int iel = 0; iel < nel; iel++){
////        TPZCompEl *cel = cmesh->ElementVec()[iel];
////        if(!cel) continue;
////
////        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
////        if(!sp) continue;
////        int level = sp->Reference()->Level();
////        TPZGeoEl * gel = sp->Reference();
////        if(gel->Dimension()==2)
////            sp->PRefine(porder + (ndiv - 1) + (level-1));
////    }
////    cmesh->AdjustBoundaryElements();
////    cmesh->CleanUpUnconnectedNodes();
////    cmesh->ExpandSolution();
////
////#ifdef LOG4CXX
////    if (logger->isDebugEnabled())
////    {
////        std::stringstream sout;
////        sout<<"malha computacional apos pRefinamento\n";
////        cmesh->Print(sout);
////        LOGPZ_DEBUG(logger, sout.str());
////    }
////#endif
////}
//
//void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder){
//    if(ndiv<1) return;
//    int nel = cmesh->NElements();
//    
//    for(int iel = 0; iel < nel; iel++){
//        TPZCompEl *cel = cmesh->ElementVec()[iel];
//        if(!cel) continue;
//        
//        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
//        if(!sp) continue;
//        int level = sp->Reference()->Level();
//        
//        TPZGeoEl * gel = sp->Reference();
//        
//        
//        //        int ordemaux=2;
//        //
//        //        if (level>4) {
//        //
//        //            ordemaux=3;
//        //        }
//        //
//        //
//        //
//        //
//        ////        if (level > 4 && level <7) {
//        ////            ordemaux=3;
//        ////        }
//        ////
//        //////        else
//        //////        if (level>6) {
//        //////            ordemaux=4;
//        //////        }
//        ////
//        ////        else ordemaux=2;
//        //
//        //        if(gel->Dimension()==2)
//        //            sp->PRefine(ordemaux/*porder + (level-3)*/);
//        
//        if(gel->Dimension()==2)
//            sp->PRefine(porder + (level-flevel));
//    }
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    cmesh->ExpandSolution();
//    
//    
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout<<"malha computacional apos pRefinamento\n";
//        cmesh->Print(sout);
//        LOGPZ_DEBUG(logger, sout.str());
//    }
//#endif
//}
//
//void PrefinamentoTeste(TPZCompMesh * cmesh, int porder){
//    int nel = cmesh->NElements();
//    
//    for(int iel = 0; iel < nel; iel++){
//        TPZCompEl *cel = cmesh->ElementVec()[iel];
//        if(!cel) continue;
//        
//        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
//        if(!sp) continue;
//        
//        TPZGeoEl * gel = sp->Reference();
//        
//        if(gel->Dimension()==2 && gel->Index()!= 0){
//            sp->PRefine(porder + 1);
//        }
//    }
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    cmesh->ExpandSolution();
//    
//}
//
//
//void SetPDifferentOrder(TPZCompMesh *comp,int porder){
//    int nel = comp->NElements();
//    int iel;
//    for(iel=0; iel<nel; iel++){
//        
//        TPZInterpolationSpace *intel;
//        TPZCompEl *cel = comp->ElementVec()[iel];
//        
//        
//        
//        intel = dynamic_cast<TPZInterpolationSpace *>(cel);
//        if(intel){
//            
//            int fator=iel%2;
//            
//            if (cel->Dimension()==2 && fator==0) {
//                
//                intel->PRefine(porder+1);
//                
//            }
//            
//            
//            if(cel->Dimension()==2 && fator!=0) {
//                
//                intel->PRefine(porder);
//                
//                
//            }
//            
//            
//        }
//    }
//    
//    
//    
//    comp->LoadReferences();
//    comp->ExpandSolution();
//    comp->AdjustBoundaryElements();
//    
//    comp->SetName("Malha Computacional cOm diferentes Ordens");
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        comp->Print(sout);
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
//    
//}
//
//TPZGeoMesh * MalhaGeo4Element(const int h){//malha quadrilatero com 4 elementos
//    TPZGeoMesh *gmesh = new TPZGeoMesh();
//    const int nelem=12;
//    TPZGeoEl *elvec[nelem];
//    //Criar ns
//    const int nnode = 9;//AQUI
//    const int dim = 2;//AQUI
//    
//    REAL co[nnode][dim] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.},{0.5,0},{1.,0.5},{0.5,1.},{0.,0.5},{0.5,0.5}};//{{-1.,-1},{1.,-1},{1.,1.},{-1.,1.},{0.,-1.},{0.,1.}};
//    
//    
//    int nodindAll[4][4]={{0,4,8,7},{4,1,5,8},{8,5,2,6},{7,8,6,3}};//como serao enumerados os nos
//    
//    
//    int nod;
//    TPZVec<REAL> coord(dim);
//    for(nod=0; nod<nnode; nod++) {
//        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
//        
//        for(int d = 0; d < dim; d++)
//        {
//            coord[d] = co[nod][d];
//        }
//        gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
//    }
//    
//    
//    
//    int matId=1;
//    int64_t id=0;
//    TPZVec <int64_t> TopolQuad(4);
//    TPZVec <int64_t> TopolLine(2);
//    //-----
//    
//    //indice dos nos
//    
//    //Criacao de elementos
//    id=0.;
//    TopolQuad[0] = 0;
//    TopolQuad[1] = 4;
//    TopolQuad[2] = 8;
//    TopolQuad[3] = 7;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//    id++;
//    // matId++;
//    
//    TopolQuad[0] = 4;
//    TopolQuad[1] = 1;
//    TopolQuad[2] = 5;
//    TopolQuad[3] = 8;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//    id++;
//    // matId++;
//    
//    TopolQuad[0] = 8;
//    TopolQuad[1] = 5;
//    TopolQuad[2] = 2;
//    TopolQuad[3] = 6;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//    id++;
//    // matId++;
//    
//    TopolQuad[0] = 7;
//    TopolQuad[1] = 8;
//    TopolQuad[2] = 6;
//    TopolQuad[3] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//    id++;
//    //  matId++;
//    
//    
//    TopolLine[0] = 0;
//    TopolLine[1] = 4;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
//    id++;
//    
//    TopolLine[0] = 4;
//    TopolLine[1] = 1;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-2,*gmesh);
//    id++;
//    
//    TopolLine[0] = 1;
//    TopolLine[1] = 5;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-3,*gmesh);
//    id++;
//    
//    TopolLine[0] = 5;
//    TopolLine[1] = 2;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-4,*gmesh);
//    id++;
//    
//    TopolLine[0] = 2;
//    TopolLine[1] = 6;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-5,*gmesh);
//    id++;
//    
//    TopolLine[0] = 6;
//    TopolLine[1] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-6,*gmesh);
//    id++;
//    
//    TopolLine[0] = 3;
//    TopolLine[1] = 7;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-7,*gmesh);
//    id++;
//    
//    TopolLine[0] = 7;
//    TopolLine[1] = 0;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-8,*gmesh);
//    
//    
//    gmesh->BuildConnectivity();
//    
//    
//#ifdef LOG4CXX
//    {
//        std::stringstream sout;
//        gmesh->Print(sout);
//        LOGPZ_DEBUG(logger, sout.str().c_str());
//    }
//#endif
//    return gmesh;
//    
//}
//
//void TesteMesh(TPZCompMesh *mesh,std::ostream &out){
//    int nEl= mesh-> NElements();
//    
//    out<< "Num. Elementos "<<nEl<<std::endl;
//    
//    for (int iel=0; iel<nEl; iel++) {
//        TPZCompEl *cel = mesh->ElementVec()[iel];
//        if (!cel) continue;
//        int ncon=cel->NConnects();
//        
//        out<< "El "<< iel<<" Num. de Connects "<<ncon<<std::endl;
//        int ordemMax=0;
//        for (int icon=0; icon<ncon; icon++) {
//            
//            int cOrder=cel->Connect(icon).Order();
//            
//            if (cel->Connect(icon).LagrangeMultiplier()) {
//                out<< "connect "<<icon<<"------connect Pressao ------ "<<" Ordem "<< cOrder<<" NShape "<< cel->Connect(icon).NShape()<<std::endl;
//            }
//            out<< "connect "<<icon<<" Ordem "<< cOrder<<" NShape "<< cel->Connect(icon).NShape()<< std::endl;
//            
//        }
//        
//        
//        if (cel->Dimension()!=1) {
//            
//            for(int icon=0;icon<ncon-2;icon++){
//                
//                int cOrder2=cel->Connect(icon).Order();
//                out<< "------connect  "<<icon<<" Ordem "<< cOrder2<<std::endl;
//                if (cOrder2>ordemMax) {
//                    cOrder2=ordemMax;
//                }
//            }
//            
//            int interCon=cel->Connect(ncon-2).Order();
//            out<< "------connect Interno "<<ncon-2<<" Ordem "<< interCon<<std::endl;
//            if (interCon < ordemMax)  DebugStop();
//        }
//        
//    }
//}
//
//void ChangeOrderExternalConnects(TPZCompMesh *mesh){
//    
//    int nEl= mesh-> NElements();
//    int dim = mesh->Dimension();
//    
//    int cordermin = -1;
//    for (int iel=0; iel<nEl; iel++) {
//        TPZCompEl *cel = mesh->ElementVec()[iel];
//        if (!cel) continue;
//        int ncon = cel->NConnects();
//        int corder = 0;
//        int nshape = 0;
//        
//        if(cel->Dimension()== dim){
//            for (int icon=0; icon<ncon-1; icon++){
//                TPZConnect &co  = cel->Connect(icon);
//                corder = co.Order();
//                nshape = co.NShape();
//                if(corder!=cordermin){
//                    cordermin = corder-1;
//                    co.SetOrder(cordermin);
//                    co.SetNShape(nshape-1);
//                    mesh->Block().Set(co.SequenceNumber(),nshape-1);
//                }
//            }
//        }
//    }
//    mesh->ExpandSolution();
//    mesh->CleanUpUnconnectedNodes();
//}
//
//
//
//int EquationsElementosConectados(TPZCompMesh *cmesh){
//    
//    int ncon = cmesh->NConnects();
//    int neqglob = 0;
//    int64_t neqcond = 0;
//    for(int i = 0; i< ncon; i++){
//        TPZConnect &co  = cmesh->ConnectVec()[i];
//        if(co.HasDependency()) continue;
//        
//        int nelc = co.NElConnected();
//        
//        if(nelc>0){
//            int dofsize = co.NShape()*co.NState();
//            if (nelc == 1) {
//                neqcond += dofsize;
//            }
//            neqglob += dofsize;
//        }
//    }
//    return neqcond;
//}
//

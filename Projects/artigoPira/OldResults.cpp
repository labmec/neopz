//
//  OldResults.cpp
//  pz
//
//  Created by Denise De Siqueira on 04/05/2018.
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
#include "tpzchangeel.h"
#include "tpzarc3d.h"
#include "tpzquadraticquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZGenGrid2D.h"

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

#include "Tools.h"

#include <iostream>
#include <math.h>
using namespace std;

//#define SmoothSol

#ifdef PZ_LOG
static TPZLogger logger("PiraHP.main");
#endif


#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#define print





int mainOld(int argc, char *argv[])
{
    
    
    
    //       IntegrationRuleConvergence(true);
    //       DebugStop();
    bool QuarterPoint = false; //geometria qp
    bool QuarterPointRule = false; //para regra de integracao qp
    
    bool CondenseQuarterPointRule = false;//para condensar os elementos usando a regra de integracao QP
    
    bool RefAfim=true;//refinamento afim
    
    bool HpRefine=false;//para p refinamento
    
    
    bool HDivMaisMais = false;//para RT+
    int order_reduce = 0;
    
    int p = 1;
    int pq = p;
    int pp = p;
    // int order=0;
    if(HDivMaisMais){
        //pq +=1;
        pp +=1;
        order_reduce = 1;
    }
    
    ///Refinamento
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
    
    
    //string outputfile("Solution_mphysics");
    std::stringstream name;
    
    
    //std::ofstream myerrorfile("MalhasHpQP.txt");
    
    std::ofstream myerrorfile("MalhaQuadUniP2.txt");
    
    for(int ndiv=1; ndiv<6; ndiv++)
    {
        
        
        TPZGeoMesh *gmesh=GMesh();//CurvedMesh2(ndiv);
        
        //UniformRefine(gmesh, 1);
        
        for(long el=0; el < gmesh->NElements(); el++)
        {
            TPZGeoEl *gel = gmesh->Element(el);
            gel->SetFather(nullptr);
            
        }
        
        
        //se refinamento  QP
        int nodeAtOriginId=1;
        if (QuarterPoint){
            
            QuarterPointRef(gmesh, nodeAtOriginId);
            
            if(RefAfim){
                UniformRefine(gmesh, ndiv);
                //DirectUniRef(gmesh, nodeAtOriginId,ndiv);
                
            }
            else{DirectionalRef(gmesh, nodeAtOriginId,ndiv);}
            
            
            
            
            
        }
        //se refinamento uniforme
        else{
            
            if(RefAfim){
                UniformRefine(gmesh,ndiv);}
            else{
                DirectionalRef(gmesh, nodeAtOriginId,ndiv); //para malha direcional
            }
            
            
            
            // #ifdef print
            //            {
            //                std::ofstream out("gmesh.txt");
            //                gmesh->Print(out);
            //            }
            //#endif
            
            
            
            //            REAL rho_min=0;
            //            REAL aspectratio=0;
            //
            //            ComputeCharacteristicHElSize(gmesh,aspectratio , rho_min);
            //            myerrorfile << "ndiv= "<< ndiv<< " aspectratio = " << aspectratio<< "\n";
            
        }
        
#ifdef print
        {
            std::ofstream malhaOut("GMesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
        }
#endif
        
        //        #ifdef print
        //                    {
        //                        std::ofstream malhaOut("malhageometrica.vtk");
        //                        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
        //                    }
        //        #endif
        
        //refinamento quarter point proximo do no de id=1
        //        int nodeAtOriginId = 1;
        //        if(QuarterPoint){
        //            int nodeAtOriginId = 1;
        //            QuarterPointRef(gmesh, nodeAtOriginId);
        //            //DirectionalRef(gmesh, nodeAtOriginId,iddivide);
        //            DirectUniRef(gmesh, nodeAtOriginId,ndiv);
        //
        //        }
        
        //#ifdef print
        //        {
        //            std::ofstream malhaOut("malhageometrica2.vtk");
        //            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
        //
        ////            std::ofstream out("gmesh.txt");
        ////            gmesh->Print(out);
        //        }
        //#endif
        
        //Testando isotropia da mallha
        
        //        REAL h_min=0;
        //        REAL rho_min=0;
        //
        //        ComputeCharacteristicHElSize(gmesh,h_min , rho_min);
        //        myerrorfile << "h_min = " << h_min << " rho_min= "<<rho_min<<"\n";
        
        
        
        //cmeshL2 flux
        TPZCompMesh * cmeshL2 = NULL;
        if(QuarterPointRule)
        {
            cmeshL2= CMeshFluxL2(gmesh, pq, nodeAtOriginId);
            if(HpRefine)  Prefinamento(cmeshL2, ndiv,pq);
            if(HDivMaisMais)
            {
                ChangeInternalConnectOrder(cmeshL2);
                //ChangeSideConnectOrderConnects(cmeshL2, order_reduce);
            }
            
            //order = pp+1;
            int order = (pp+2*ndiv - 1)+1;
            int max_order = 2*order;
            int Gauss_order = max_order+1;
            int QP_order =2*Gauss_order;
            
            ChangeIntegrationRule(cmeshL2, QP_order,true);
            //            {
            //                std::ofstream out("cmeshL2.txt");
            //                cmeshL2->Print(out);
            //
            //                std::ofstream out2("gmeshL2.txt");
            //                cmeshL2->Reference()->Print(out2);
            //            }
        }
        
        //mesh1
        TPZCompMesh * cmesh1= CMeshFlux(gmesh, pq);
        //        {
        //            std::ofstream out("cmeshAntesPRef.txt");
        //            cmesh1->Print(out);
        //        }
        
        if(HpRefine)   Prefinamento(cmesh1, ndiv,pq);
        
        //        {
        //            std::ofstream out("cmeshPosPRef.txt");
        //            cmesh1->Print(out);
        //        }
        if(HDivMaisMais)
        {
            ChangeInternalConnectOrder(cmesh1);
            //ChangeSideConnectOrderConnects(cmesh1, order_reduce);
        }
        
        //        {
        //            std::ofstream out("cmesh1.txt");
        //            cmesh1->Print(out);
        //        }
        
        //mesh2
        TPZCompMesh * cmesh2= CMeshPressure(gmesh, pp);
        if(HpRefine) Prefinamento(cmesh2, ndiv,pp);
        //        {
        //            std::ofstream out("cmesh2.txt");
        //            cmesh2->Print(out);
        //        }
        
        //malha multifisica
        TPZVec<TPZCompMesh *> meshvec(2);
        meshvec[0] = cmesh1;
        meshvec[1] = cmesh2;
        TPZMixedPoisson * mymaterial;
        TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,mymaterial,QuarterPointRule);
        
        
        //std::cout << "NEquations " << mphysics->NEquations() << std::endl;
        //        {
        //            std::ofstream out("cmeshMPhysics.txt");
        //            mphysics->Print(out);
        //
        ////            std::ofstream out2("gmeshMultifisic.txt");
        ////            mphysics->Reference()->Print(out2);
        //        }
        
        //resolver problema
        TPZAnalysis anMP(mphysics);
        if(!QuarterPointRule){
            ResolverSistema(anMP, mphysics,3);
            //long averageband_Depois = ComputeAverageBandWidth(mphysics);
            //cout<<"averageband Depois de Reenumerar = " <<averageband_Depois<<std::endl;
        }
        else
        {//Transferir solucao de uma malha para outra
            
            TPZSkylineStructMatrix strMP(mphysics); //caso simetrico
            //TPZParFrontStructMatrix<TPZFrontSym<STATE> > strMP(mphysics);
            //strMP.SetDecomposeType(ELDLt);
            strMP.SetNumThreads(6);
            anMP.SetStructuralMatrix(strMP);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt); //caso simetrico
            anMP.SetSolver(step);
            anMP.Assemble();
            
            
            //Matriz projecao L2
            TPZAnalysis anFlux(cmeshL2);
            TPZSkylineStructMatrix strFlux(cmeshL2);
            //TPZParFrontStructMatrix<TPZFrontSym<STATE> > strFlux(cmeshL2);
            //strFlux.SetDecomposeType(ELDLt);
            strFlux.SetNumThreads(6);
            anFlux.SetStructuralMatrix(strFlux);
            TPZStepSolver<STATE> step2;
            anFlux.SetSolver(step2);
            anFlux.Assemble();
            
            // Transferir EkL2 para EkMultiPhisica
            TPZAutoPointer< TPZMatrix<STATE> > matF =  anFlux.Solver().Matrix();
            TPZAutoPointer< TPZMatrix<STATE> > matMP =  anMP.Solver().Matrix();
            
            // std::ofstream saidamatriz("../saidamatriz.nb");
            //            matF->Print("MatFlux = ", saidamatriz,EMathematicaInput);
            //            matMP->Print("\n\nMatAntesTransfer = ", saidamatriz,EMathematicaInput);
            //            anMP.Rhs().Print("\n\nRhsAntesTransfer = ", saidamatriz,EMathematicaInput);
            
            
            
            //#ifdef PZ_LOG
            //            if(logger.isDebugEnabled())
            //            {
            //                std::stringstream sout;
            //                matF->Print("MatFlux = ", sout,EMathematicaInput);
            //                matMP->Print("MatAntesTransfer = ", sout,EMathematicaInput);
            //                LOGPZ_DEBUG(logger,sout.str())
            //            }
            //#endif
            //             {
            //             std::ofstream out("cmeshMPhysics.txt");
            //             mphysics->Print(out);
            //
            //             std::ofstream out2("cmeshL2.txt");
            //             cmeshL2->Print(out2);
            //             }
            
            //            std::ofstream subMat("../SubMat.nb");
            //            GlobalSubMatrix(mphysics, matMP, nodeAtOriginId, true, subMat);
            
            
            TransferMatrixFromMeshes(cmeshL2, mphysics, matF, matMP,nodeAtOriginId);
            
            //GlobalSubMatrix(mphysics, anMP.Solver().Matrix(), nodeAtOriginId, false, subMat);
            
            if(CondenseQuarterPointRule)//resolve com condensacao estatica
            {
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
                
                mphysics->CleanUpUnconnectedNodes();
                mphysics->ExpandSolution();
                
                TPZAnalysis newMP(mphysics);
                ResolverSistema(newMP, mphysics,3);
            }
            else//resolve sem condensar
            {
                anMP.Solve();
            }
            
            
            
            //   anMP.Solver().Matrix()->Print("\n\nMatAposTransfer = ", saidamatriz,EMathematicaInput);
            //   anMP.Rhs().Print("\n\nRhsAposTransfer = ", saidamatriz,EMathematicaInput);
            
            {
                //             std::ofstream out("cmeshMPhysicsfinal.txt");
                //             mphysics->Print(out);
                //
                //             std::ofstream out2("cmeshL2final.txt");
                //             cmeshL2->Print(out2);
            }
            //            saidamatriz << "\n\nMatrixForm[MatFlux]"<<std::endl;
            //            saidamatriz << "\nMatrixForm[MatAntesTransfer]"<<std::endl;
            //            saidamatriz << "\nMatrixForm[MatAposTransfer]"<<std::endl;
            //            saidamatriz << "\nMatrixForm[MatAntesTransfer - MatAposTransfer]"<<std::endl;
            //            saidamatriz << "\nMatrixForm[RhsAntesTransfer - RhsAposTransfer]"<<std::endl;
            
            
            
            //#ifdef PZ_LOG
            //        if(logger.isDebugEnabled())
            //        {
            //            std::stringstream sout;
            //            anMP.Solver().Matrix()->Print("MatAposTransfer = ", sout,EMathematicaInput);
            //            LOGPZ_DEBUG(logger,sout.str())
            //        }
            //#endif
        }
        
        //pos-process
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
        std::stringstream name;
        name << "TesteUni" <<ndiv<< ".vtk";
        std::string paraviewfile(name.str());
        PosProcessMultphysics(meshvec,  mphysics, anMP, paraviewfile);
        //Erro global
        myerrorfile << "\n-------------------------------------------------- \n";
        myerrorfile << "h = "<< ndiv << " p = " << p << "\n";
        myerrorfile << "DOF Total = " << cmesh1->NEquations() + cmesh2->NEquations()<< "\n";
        myerrorfile << "DOF Condensado = " << mphysics->NEquations() << "\n";
        TPZVec<REAL> erros(3);
        //myerrorfile<<"\n\nErro da simulacao multifisica  para o flux";
        //saidamatriz<<"\nP = " <<pq<<";"<<std::endl;
        //cmesh1->Solution().Print("SolfluxQP = ", saidamatriz,EMathematicaInput);
        //saidamatriz<<"(*---*)"<<std::endl;
        
        
        
        
        if(QuarterPoint)
        {
            int order = (pp+2*ndiv - 1)+1;
            ChangeIntegrationRule(cmesh1,2*order,true);
            ChangeIntegrationRule(cmesh2, 2*order,true);
        }
        
        myerrorfile<<"\n\nErro da simulacao multifisica  para o fluxo";
        ErrorHDiv(cmesh1,myerrorfile);
        // ErroHDivNoElemento(cmesh1,myerrorfile, nodeAtOriginId);
        
        myerrorfile<<"\n\nErro da simulacao multifisica  para a pressao";
        ErrorL2(cmesh2,myerrorfile);
        //erro L2 no elemento
        //ErroL2NoElemento(cmesh2,myerrorfile, nodeAtOriginId);
        
        cmesh1->CleanUp();
        cmesh2->CleanUp();
        delete cmesh1;
        delete cmesh2;
        delete gmesh;
        //
        
        
    }
    
    
    return 0;
}

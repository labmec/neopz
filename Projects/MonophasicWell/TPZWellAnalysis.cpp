///*
// *  pznondarcyanalysis.cpp
// *  PZ
// *
// *  Created by Omar Duran Triana on 5/21/14.
// *  Copyright 2014 __MyCompanyName__. All rights reserved.
// *
// */
//
//#include "TPZWellAnalysis.h"
//#include "pzcheckgeom.h"
//#include "pzlog.h"
//
//#include "TPZReadGIDGrid.h"
//#include "tpzgeoelrefpattern.h"
//#include "TPZRefPatternTools.h"
//#include "TPZGeoLinear.h"
//#include "tpztriangle.h"
//#include "pzgeoquad.h"
//#include "pzgeopoint.h"
//
//#include "tpzhierarquicalgrid.h"
//
//#include "TPZVTKGeoMesh.h"
//#include "TPZMonoPhaseWell.h"
//#include "pzpoisson3d.h"
//#include "pzbuildmultiphysicsmesh.h"
//#include "TPZSkylineNSymStructMatrix.h"
//#include "pzfstrmatrix.h"
//#include "math.h"
//
//#ifdef PZDEBUG
//    #ifdef LOG4CXX
//    static PZLogger logger("pz.DarcyFlow");
//    #endif
//#endif
//
//TPZWellAnalysis::TPZWellAnalysis()
//{
//    
//    fmeshvec.Resize(4);
//    
//    fgmesh=NULL;
//    fcmeshdarcy=NULL;
//    
//    // Vector which will store tha residuum in the last state (n)
//    fResidualAtn.Resize(0, 0);
//    
//    // Vector which will store tha residuum in the last state (n+1)
//    fResidualAtnplusOne.Resize(0, 0);
//    
//    /** @brief unknowns for n time step */
//    falphaAtn.Resize(0, 0);
//    
//    /** @brief unknowns for n+1 time step */
//    falphaAtnplusOne.Resize(0, 0);
//    
//}
//
//
//TPZWellAnalysis::~TPZWellAnalysis()
//{
//    
//}
//
//void TPZWellAnalysis::SetLastState()
//{
//
//}
//
//void TPZWellAnalysis::SetNextState()
//{
//
//}
//
//void TPZWellAnalysis::Assemble()
//{
//    
//}
//
//void TPZWellAnalysis::AssembleLastStep(TPZAnalysis *an)
//{
//    fcmeshdarcy->LoadSolution(falphaAtn);
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
//    SetLastState();
//    an->AssembleResidual();
//    fResidualAtn = an->Rhs();
//    
//    // #ifdef PZDEBUG
//    //     #ifdef LOG4CXX
//    //         if(logger.isDebugEnabled())
//    //         {
//    //             std::stringstream sout;
//    //             fResidualAtn.Print("fResidualAtn = ", sout,EMathematicaInput);
//    //             LOGPZ_DEBUG(logger,sout.str())
//    //         }
//    //     #endif
//    // #endif
//    
//}
//
//void TPZWellAnalysis::AssembleNextStep(TPZAnalysis *an)
//{
//    fcmeshdarcy->LoadSolution(falphaAtnplusOne);
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
//    SetNextState();
//    an->Assemble();
//    fResidualAtnplusOne = an->Rhs();
//    
//    // #ifdef PZDEBUG
//    //     #ifdef LOG4CXX
//    //         if(logger.isDebugEnabled())
//    //         {
//    //             std::stringstream sout;
//    //             fResidualAtnplusOne.Print("fResidualAtnplusOne = ", sout,EMathematicaInput);
//    //             LOGPZ_DEBUG(logger,sout.str())
//    //         }
//    //     #endif
//    // #endif
//    
//}
//
//void TPZWellAnalysis::UpDateAlphaVec(TPZFMatrix<REAL> &alpha)
//{
//    falphaAtn = alpha;
//    falphaAtnplusOne = alpha;
//    
//    // #ifdef PZDEBUG
//    //     #ifdef LOG4CXX
//    //         if(logger.isDebugEnabled())
//    //         {
//    //             std::stringstream sout;
//    //             falphaAtn.Print("falphaAtn = ", sout,EMathematicaInput);
//    //             falphaAtnplusOne.Print("falphaAtnplusOne = ", sout,EMathematicaInput);
//    //             LOGPZ_DEBUG(logger,sout.str())
//    //         }
//    //     #endif
//    // #endif
//    
//}
//
//
//
//void TPZWellAnalysis::InitializeSolution(TPZAnalysis *an)
//{
//    
//    // Compute the intial saturation distribution
//    int nalpha = fcmeshdarcy->Solution().Rows();
//    falphaAtn.Resize(nalpha, 1);
//    falphaAtnplusOne.Resize(nalpha, 1);
//    falphaAtn.Zero();
//    falphaAtnplusOne.Zero();
//    an->LoadSolution(falphaAtn);
//    
//    int nsoil = fmeshvec[3]->Solution().Rows();
//    TPZFMatrix<REAL> SOil(nsoil,1,1.0);
//    fmeshvec[3]->LoadSolution(SOil);
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshdarcy);
//    
//    falphaAtn = fcmeshdarcy->Solution();
//    falphaAtnplusOne = fcmeshdarcy->Solution();
//    
//    // #ifdef PZDEBUG
//    //     #ifdef LOG4CXX
//    //         if(logger.isDebugEnabled())
//    //         {
//    //             std::stringstream sout;
//    //             falphaAtn.Print("falphaAtn = ", sout,EMathematicaInput);
//    //             falphaAtnplusOne.Print("falphaAtnplusOne = ", sout,EMathematicaInput);
//    //             LOGPZ_DEBUG(logger,sout.str())
//    //         }
//    //     #endif
//    // #endif
//    
//}
//
//void TPZWellAnalysis::Run()
//{
//    
//    std::string dirname = PZSOURCEDIR;
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    
//    
//    //  Reading mesh
//    std::string GridFileName;
//    GridFileName = dirname + "/Projects/DarcyflowAxisymmetricHdiv/";
//    GridFileName += "SingleLayer.dump";
//    //GridFileName += "MixLayer.dump";
//    //GridFileName += "BatatacoarseQ.dump";
//    //GridFileName += "QUAD4.dump";
//    
//    int nx = 100;
//    int ny = 1;
//    
//    
//    REAL deg = 0.0;
//    int hcont = 0;
//    RotateGeomesh(deg * M_PI/180.0);
//    this->UniformRefinement(fSimulationData->GetHrefinement());
//    
//    std::set<int> matidstoRef;
//    //    matidstoRef.insert(2);
//    //    matidstoRef.insert(3);
//    //    matidstoRef.insert(4);
//    matidstoRef.insert(3);
//    matidstoRef.insert(5);
//    
//    this->UniformRefinement(hcont, matidstoRef);
//    this->PrintGeoMesh();
//    
//    
//    
//    int q = fSimulationData->Getqorder();
//    int p = fSimulationData->Getporder();
//    int s = 0;
//    
//    //    if (fSimulationData->GetIsH1approx())
//    if (false)
//    {
//        //        CmeshH1(p);
//    }
//    else
//    {
//        CreateMultiphysicsMesh(q,p,s);
//        CreateInterfaces();
//    }
//    
//    
//    // Analysis
//    bool mustOptimizeBandwidth = false;
//    TPZAnalysis *an = new TPZAnalysis(fcmeshdarcy,mustOptimizeBandwidth);
//    int numofThreads = 0;
//    
//    bool IsDirecSolver = fSimulationData->GetIsDirect();
//    
//    if (IsDirecSolver) {
//        
//        if (fSimulationData->GetIsBroyden()) {
//            TPZFStructMatrix fullMatrix(fcmeshdarcy);
//            TPZStepSolver<STATE> step;
//            fullMatrix.SetNumThreads(numofThreads);
//            step.SetDirect(ELU);
//            an->SetStructuralMatrix(fullMatrix);
//            an->SetSolver(step);
//        }
//        else{
//            
//            TPZSkylineNSymStructMatrix skylnsym(fcmeshdarcy);
//            TPZStepSolver<STATE> step;
//            skylnsym.SetNumThreads(numofThreads);
//            step.SetDirect(ELU);
//            an->SetStructuralMatrix(skylnsym);
//            an->SetSolver(step);
//        }
//        
//    }
//    else
//    {
//        if (fSimulationData->GetIsBroyden()) {
//            TPZFStructMatrix fullMatrix(fcmeshdarcy);
//            fullMatrix.SetNumThreads(numofThreads);
//            
//            TPZAutoPointer<TPZMatrix<STATE> > fullMatrixa = fullMatrix.Create();
//            TPZAutoPointer<TPZMatrix<STATE> > fullMatrixaClone = fullMatrixa->Clone();
//            
//            TPZStepSolver<STATE> *stepre = new TPZStepSolver<STATE>(fullMatrixaClone);
//            TPZStepSolver<STATE> *stepGMRES = new TPZStepSolver<STATE>(fullMatrixa);
//            TPZStepSolver<STATE> *stepGC = new TPZStepSolver<STATE>(fullMatrixa);
//            stepre->SetDirect(ELU);
//            stepre->SetReferenceMatrix(fullMatrixa);
//            stepGMRES->SetGMRES(10, 20, *stepre, 1.0e-10, 0);
//            stepGC->SetCG(10, *stepre, 1.0e-10, 0);
//            if (fSimulationData->GetIsCG()) {
//                an->SetSolver(*stepGC);
//            }
//            else{
//                an->SetSolver(*stepGMRES);
//            }
//            
//        }
//        else{
//            
//            TPZSkylineNSymStructMatrix skylnsym(fcmeshdarcy);
//            skylnsym.SetNumThreads(numofThreads);
//            
//            TPZAutoPointer<TPZMatrix<STATE> > skylnsyma = skylnsym.Create();
//            TPZAutoPointer<TPZMatrix<STATE> > skylnsymaClone = skylnsyma->Clone();
//            
//            TPZStepSolver<STATE> *stepre = new TPZStepSolver<STATE>(skylnsymaClone);
//            TPZStepSolver<STATE> *stepGMRES = new TPZStepSolver<STATE>(skylnsyma);
//            TPZStepSolver<STATE> *stepGC = new TPZStepSolver<STATE>(skylnsyma);
//            
//            stepre->SetDirect(ELU);
//            stepre->SetReferenceMatrix(skylnsyma);
//            stepGMRES->SetGMRES(10, 20, *stepre, 1.0e-10, 0);
//            stepGC->SetCG(10, *stepre, 1.0e-10, 0);
//            if (fSimulationData->GetIsCG()) {
//                an->SetSolver(*stepGC);
//            }
//            else{
//                an->SetSolver(*stepGMRES);
//            }
//        }
//        
//    }
//    
//    this->InitializeSolution(an);
//    this->TimeForward(an);
//    
//}
//
//
//void TPZWellAnalysis::CreateInterfaces()
//{
//    fgmesh->ResetReference();
//    fcmeshdarcy->LoadReferences();
//    
//    // Creation of interface elements
//    int nel = fcmeshdarcy->ElementVec().NElements();
//    for(int el = 0; el < nel; el++)
//    {
//        TPZCompEl * compEl = fcmeshdarcy->ElementVec()[el];
//        if(!compEl) continue;
//        TPZGeoEl * gel = compEl->Reference();
//        if(!gel) {continue;}
//        if(gel->HasSubElement()) {continue;}
//        int index = compEl ->Index();
//        if(compEl->Dimension() == fcmeshdarcy->Dimension())
//        {
//            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fcmeshdarcy->ElementVec()[index]);
//            if(!InterpEl) continue;
//            InterpEl->CreateInterfaces();
//        }
//    }
//}
//
//
//void TPZWellAnalysis::PrintLS(TPZAnalysis *an)
//{
//    //    an->Assemble();
//    TPZAutoPointer< TPZMatrix<REAL> > KGlobal;
//    TPZFMatrix<STATE> FGlobal;
//    KGlobal =   an->Solver().Matrix();
//    FGlobal =   an->Rhs();
//    
//#ifdef PZDEBUG
//    #ifdef LOG4CXX
//        if(logger.isDebugEnabled())
//        {
//            std::stringstream sout;
//            KGlobal->Print("KGlobal = ", sout,EMathematicaInput);
//            FGlobal.Print("FGlobal = ", sout,EMathematicaInput);
//            LOGPZ_DEBUG(logger,sout.str())
//        }
//    #endif
//#endif
//    
//}
//
//void TPZWellAnalysis::CreateMultiphysicsMesh(int q, int p, int s)
//{
//    fmeshvec[0] = CmeshFlux(q);
//    fmeshvec[1] = CmeshPressure(p);
//    fmeshvec[2] = CmeshSw(s);
//    fmeshvec[3] = CmeshSo(s);
//    
//    fcmeshdarcy = CmeshMixed();
//    
//    // Transferindo para a multifisica
//    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshdarcy);
//    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshdarcy);
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshdarcy);
//    
//#ifdef PZDEBUG
//    std::ofstream dumpfile("ComputationaMeshMultiphysics.txt");
//    fcmeshdarcy->Print(dumpfile);
//#endif
//    
//}
//
//void TPZWellAnalysis::TimeForward(TPZAnalysis *an)
//{
//    int nsteps = fSimulationData->GetMaxTime() / fSimulationData->GetDeltaT();
//    
//    REAL tk = 0;
//    
//    this->fSimulationData->SetTime(tk);
//    this->PostProcessVTK(an);
//    
//    for (int istep = 1 ; istep <=  nsteps; istep++) {
//        tk = istep*this->fSimulationData->GetDeltaT();
//        this->fSimulationData->SetTime(tk);
//        
//        std::cout << "Begin of time (days): " << tk/86400.0 << std::endl;
//        std::cout<<  "Time step: " << istep << std::endl;
//        
//        
//        this->AssembleLastStep(an);
//        this->AssembleNextStep(an);
//        
//        const clock_t tinia = clock();
//        NewtonIterations(an);
//        const clock_t tenda = clock();
//        const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
//        std::cout << "Time for Newton: " << timea << std::endl;
//        std::cout << "Number of DOF = " << fcmeshdarcy->Solution().Rows() << std::endl;
//        this->PostProcessVTK(an);
//        
//        
//        
//    }
//    
//    
//    
//}
//
//
//void TPZWellAnalysis::NewtonIterations(TPZAnalysis *an)
//{
//    
//    TPZFMatrix<STATE> Residual(an->Rhs().Rows(),1,0.0);
//
//    //    TPZManVector<int64_t> Actives(0),NoActives(0);
//    //
//    //    this->FilterSaturations(Actives, NoActives);
//    //    an->StructMatrix()->EquationFilter().Reset();
//    //    an->StructMatrix()->EquationFilter().SetActiveEquations(Actives);
//    
//    Residual = fResidualAtn + fResidualAtnplusOne;
//    
//    TPZFMatrix<STATE> X = falphaAtn;
//    TPZFMatrix<STATE> DeltaX = falphaAtn;
//    
//    STATE error     =   1;
//    STATE normdx    =   1;
//    int iterations  =   0;
//    int centinel    =   0;
//    int fixed       =   fSimulationData->GetFixediterations();
//    
//    while (error >= fSimulationData->GetToleranceRes() && iterations <= fSimulationData->GetMaxiterations()) {
//        
//        an->Rhs() = Residual;
//        an->Rhs() *= -1.0;
//        
//        an->Solve();
//
//        
//        DeltaX = an->Solution();
//        normdx = Norm(DeltaX);
//        X += DeltaX;
//        
//        fcmeshdarcy->LoadSolution(X);
//        if (true)
//        {
//            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
//        }
//        
//        if (((fixed+1) * (centinel) == iterations)) {
//            an->Assemble();
//            centinel++;
//        }
//        else{
//            an->AssembleResidual();
//            
//        }
//        fResidualAtnplusOne = an->Rhs();
//        
//        Residual = fResidualAtn + fResidualAtnplusOne;
//        error = Norm(Residual);
//        iterations++;
//        
////#ifdef PZDEBUG
////    #ifdef LOG4CXX
////            if(logger.isDebugEnabled())
////            {
////                std::stringstream sout;
////                fResidualAtn.Print("fResidualAtn = ", sout,EMathematicaInput);
////                fResidualAtnplusOne.Print("fResidualAtnplusOne = ", sout,EMathematicaInput);
////                DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
////                X.Print("X = ", sout,EMathematicaInput);
////                LOGPZ_DEBUG(logger,sout.str())
////            }
////    #endif
////#endif
//        
//        if(error < fSimulationData->GetToleranceRes() || normdx < fSimulationData->GetToleranceDX())
//        {
//            std::cout << "Converged with iterations:  " << iterations << std::endl;
//            std::cout << "error norm: " << error << std::endl;
//            std::cout << "error of dx: " << normdx << std::endl;
//            this->UpDateAlphaVec(X);
//            break;
//        }
//        
//        if (iterations == fSimulationData->GetMaxiterations()) {
//            std::cout << "Out max iterations " << iterations << std::endl;
//            std::cout << "error norm " << error << std::endl;
//            this->UpDateAlphaVec(X);
//            break;
//        }
//        
//    }
//    
//    
//}
//
//void TPZWellAnalysis::BroydenIterations(TPZAnalysis *an)
//{
//    
//    
//    int m = an->Solution().Rows();
//    
//    TPZFMatrix<STATE> Residual(m,1,0.0);
//    TPZFMatrix<STATE> Rank(m,m,0.0);
//    TPZFMatrix<STATE> X(m,1,0.0);
//    TPZFMatrix<STATE> DeltaX(m,1,0.0);
//    
//    TPZAutoPointer<TPZMatrix<STATE> >  D;
//    TPZAutoPointer<TPZMatrix<STATE> >  Dk;
//    
//    TPZFMatrix<STATE> * DInverse =  new TPZFMatrix<STATE> (m,m,0.0);
//    
//    bool IsShermanMorrison = false;
//    STATE ck = 0.0;
//    //    TPZFMatrix<STATE> DInverse(m,m,0.0);
//    TPZFMatrix<STATE> Identity(m,m,0.0);
//    TPZFMatrix<STATE> DInverseT(m,m,0.0);
//    TPZFMatrix<STATE> u(m,m,0.0);
//    TPZFMatrix<STATE> du(m,m,0.0);
//    TPZFMatrix<STATE> duT(m,m,0.0);
//    TPZFMatrix<STATE> Ck(1,1,0.0);
//    
//    STATE error=1;
//    STATE normdx=1;
//    int iterations=0;
//    int centinel    =0;
//    int fixed       =fSimulationData->GetFixediterations();
//    
//    // Computing the first Newton iteration
//    
//    
//    an->Rhs() = fResidualAtn + fResidualAtnplusOne;       // g(X0)
//    
//    
//    an->Rhs() *= -1.0;
//    
//    if (IsShermanMorrison) {
//        DInverse = ComputeInverse();
//        DInverse->Multiply(Residual, DeltaX);
//        D.operator->()->Multiply(*DInverse, Identity);
//        
//#ifdef PZDEBUG
//#ifdef LOG4CXX
//        if(logger.isDebugEnabled())
//        {
//            std::stringstream sout;
//            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
//            Residual.Print("Residual = ", sout,EMathematicaInput);
//            Identity.Print("Identity = ", sout,EMathematicaInput);
//            DInverse->Print("DInverse = ", sout,EMathematicaInput);
//            LOGPZ_DEBUG(logger,sout.str())
//        }
//#endif
//#endif
//        
//    }
//    else{
//        D = an->Solver().Matrix();  // J(X0)
//        an->Solve();
//        DeltaX = an->Solution();    // d(X0)
//    }
//    normdx = Norm(DeltaX);      // d(X0)*d(X0)
//    X += DeltaX;                // X1
//    
//    // End of newton iteration
//    
//    // Procedure without Inverse computation
//    
//    iterations++;
//    
//    while (error >= fSimulationData->GetToleranceRes() && iterations <= fSimulationData->GetMaxiterations()) {
//        
//        
//        fcmeshdarcy->LoadSolution(X);
//        if (!fSimulationData->GetIsH1approx())
//        {
//            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
//        }
//        
//        this->AssembleResidual();
//        Residual = fResidualAtn + fResidualAtnplusOne;       // g(Xk)
//        error = Norm(Residual);     // g(Xk)*g(Xk)
//        
//        if(error <= fSimulationData->GetToleranceRes() || normdx <= fSimulationData->GetToleranceDX())
//        {
//            std::cout << "Converged with iterations:  " << iterations << std::endl;
//            std::cout << "error norm: " << error << std::endl;
//            std::cout << "error of dx: " << normdx << std::endl;
//            break;
//        }
//#ifdef PZDEBUG
//#ifdef LOG4CXX
//        if(logger.isDebugEnabled())
//        {
//            std::stringstream sout;
//            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
//            Residual.Print("Residual = ", sout,EMathematicaInput);
//            LOGPZ_DEBUG(logger,sout.str())
//        }
//#endif
//#endif
//        
//        if (((fixed+1) * (centinel) == iterations)) {
//            
//            if (IsShermanMorrison) {
//                DInverse->Multiply(Residual, u);
//                u.Add(DeltaX, du);
//                u.Transpose(&duT);
//                DeltaX.Multiply(duT,Ck);
//                ck = Ck(0,0);
//                Rank = TensorProduct(u, DeltaX);
//                Rank *= -(1/ck);
//                Rank.Multiply(*DInverse, DInverseT);
//                DInverse->Add(DInverseT, *DInverse);
//                
//            }else
//            {
//                // Application of the Secant condition
//                Rank = TensorProduct(Residual, DeltaX);
//                Rank *= 1.0/(normdx*normdx);
//                D.operator->()->Add(Rank, *D.operator->());
//                an->Solver().Matrix() = D;
//            }
//            
//            
//            
//            centinel++;
//            
//        }
//        
//        
//        
//        //#ifdef LOG4CXX
//        //        if(logger.isDebugEnabled())
//        //        {
//        //            std::stringstream sout;
//        //            an->Solver().Matrix().operator->()->Print("*an->Solver().Matrix().operator->() = ", sout,EMathematicaInput);
//        //            D->Print("D = ", sout,EMathematicaInput);
//        //            Rank.Print("Rank = ", sout,EMathematicaInput);
//        //            LOGPZ_DEBUG(logger,sout.str())
//        //        }
//        //#endif
//        an->Rhs() = fResidualAtn + fResidualAtnplusOne;
//        an->Rhs() *= -1.0;
//        
//        if (IsShermanMorrison) {
//            DInverse->Multiply(Residual, DeltaX);
//        }
//        else{
//            an->Solve();
//            DeltaX = an->Solution();    // d(Xk)
//        }
//        
//        
//        normdx = Norm(DeltaX);      // d(Xk)*d(Xk)
//        X += DeltaX;                // Xk+1
//        
//        
//        //#ifdef LOG4CXX
//        //        if(logger.isDebugEnabled())
//        //        {
//        //            std::stringstream sout;
//        //            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
//        //            X.Print("X = ", sout,EMathematicaInput);
//        //            LOGPZ_DEBUG(logger,sout.str())
//        //        }
//        //#endif
//        iterations++;
//        
//        if (iterations == fSimulationData->GetMaxiterations()) {
//            std::cout << "Out max iterations " << iterations << std::endl;
//            std::cout << "error norm " << error << std::endl;
//            break;
//        }
//        
//    }
//    
//    
//}
//
//
//TPZFMatrix<STATE>  TPZWellAnalysis::TensorProduct(TPZFMatrix<STATE> &g, TPZFMatrix<STATE> &d)
//{
//    TPZFMatrix<STATE> dT=d;
//    d.Transpose(&dT);
//    TPZFMatrix<STATE> RankOne;
//    g.Multiply(dT, RankOne);
//    
//    //#ifdef LOG4CXX
//    //    if(logger.isDebugEnabled())
//    //    {
//    //        std::stringstream sout;
//    //        g.Print("g = ", sout,EMathematicaInput);
//    //        dT.Print("dT = ", sout,EMathematicaInput);
//    //        d.Print("d = ", sout,EMathematicaInput);
//    //        RankOne.Print("RankOne = ", sout,EMathematicaInput);
//    //        LOGPZ_DEBUG(logger,sout.str())
//    //    }
//    //#endif
//    
//    return RankOne;
//    
//}
//
//TPZCompMesh * TPZWellAnalysis::CmeshMixed()
//{
//    int dim = 2;
//    int ilayer = 0;
//    int RockId = fLayers[ilayer]->GetMatIDs()[0];
//    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
//    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
//    int topId = fLayers[ilayer]->GetMatIDs()[3];
//    int leftId = fLayers[ilayer]->GetMatIDs()[4];
//    
//    int typeFluxin = 1, typePressurein = 0;
//    int typeFluxout = 3, typePressureout = 2;
//    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
//    
//    // Malha computacional
//    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
//    
//    // Material medio poroso
//    TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
//    mat->SetSimulationData(fSimulationData);
//    mat->SetReservoirData(fLayers[ilayer]);
//    mat->SetPetroPhysicsData(fRockPetroPhysic[ilayer]);
//    mat->SetFluidModelData(fFluidData);
//    cmesh->InsertMaterialObject(mat);
//    
//    
//    // Rigth hand side function
//    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Ffunction);
//    TPZAutoPointer<TPZFunction<STATE> > forcef;
//    dum->SetPolynomialOrder(0);
//    forcef = dum;
//    mat->SetForcingFunction(forcef);
//    
//    // Setting up linear tracer solution
//    TPZDummyFunction<STATE> *Ltracer = new TPZDummyFunction<STATE>(LinearTracer);
//    TPZAutoPointer<TPZFunction<STATE> > fLTracer = Ltracer;
//    mat->SetTimeDependentFunctionExact(fLTracer);
//    
//    
//    // Bc Bottom
//    val2(0,0) = 0.0;
//    val2(1,0) = 0.0;
//    val2(2,0) = 0.0;
//    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFluxin, val1, val2);
//    
//    // Bc Right
//    val2(0,0) = 10.0*1e6;
//    val2(1,0) = 0.0;
//    val2(2,0) = 0.0;
//    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressureout, val1, val2);
//    
//    // Bc Top
//    val2(0,0) = 0.0;
//    val2(1,0) = 0.0;
//    val2(2,0) = 0.0;
//    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFluxin, val1, val2);
//    
//    // Bc Left
//    val2(0,0) = -0.00001;
//    val2(1,0) = 1.0;
//    val2(2,0) = 0.0;
//    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFluxin, val1, val2);
//    
//    cmesh->InsertMaterialObject(bcBottom);
//    cmesh->InsertMaterialObject(bcRight);
//    cmesh->InsertMaterialObject(bcTop);
//    cmesh->InsertMaterialObject(bcLeft);
//    
//    
//    cmesh->SetDimModel(dim);
//    cmesh->SetAllCreateFunctionsMultiphysicElem();
//    cmesh->AutoBuild();
//    
//    
//    return cmesh;
//}
//
//void TPZWellAnalysis::CmeshH1(int porder)
//{
//    int dim = 2;
//    int ilayer = 0;
//    int RockId = fLayers[ilayer]->GetMatIDs()[0];
//    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
//    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
//    int topId = fLayers[ilayer]->GetMatIDs()[3];
//    int leftId = fLayers[ilayer]->GetMatIDs()[4];
//    
//    const int typeFlux = 1, typePressure = 0;
//    TPZFMatrix<STATE> val1(1,2,0.), val2(1,1,0.);
//    
//    // Malha computacional
//    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
//    
//    // Material medio poroso
//    TPZAxiSymmetricDarcyFlowH1 * mat = new TPZAxiSymmetricDarcyFlowH1(RockId);
//    mat->SetReservoirData(fLayers[ilayer]);
//    cmesh->InsertMaterialObject(mat);
//    
//    // Rigth hand side function
//    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Ffunction);
//    TPZAutoPointer<TPZFunction<STATE> > forcef;
//    dum->SetPolynomialOrder(0);
//    forcef = dum;
//    mat->SetForcingFunction(forcef);
//    
//    // Bc Bottom
//    val2(0,0) = 0.0;
//    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcBottom);
//    
//    // Bc Right
//    val2(0,0) = 1.0;
//    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
//    cmesh->InsertMaterialObject(bcRight);
//    
//    // Bc Top
//    val2(0,0) = 0.0;
//    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcTop);
//    
//    // Bc Left
//    val2(0,0) = -1.0;
//    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcLeft);
//    
//    
//    cmesh->SetDimModel(dim);
//    cmesh->SetDefaultOrder(porder);
//    cmesh->SetAllCreateFunctionsContinuous();
//    cmesh->AutoBuild();
//    
//#ifdef PZDEBUG
//    std::ofstream out("cmeshPressureH1.txt");
//    cmesh->Print(out);
//#endif
//    
//    fcmeshdarcy = cmesh;
//    
//}
//
//
//
//TPZCompMesh * TPZWellAnalysis::CmeshFlux(int qorder)
//{
//    
//    int dim = 2;
//    int ilayer = 0;
//    int RockId = fLayers[ilayer]->GetMatIDs()[0];
//    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
//    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
//    int topId = fLayers[ilayer]->GetMatIDs()[3];
//    int leftId = fLayers[ilayer]->GetMatIDs()[4];
//    
//    const int typeFlux = 0, typePressure = 1;
//    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
//    
//    // Malha computacional
//    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
//    
//    TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
//    cmesh->InsertMaterialObject(mat);
//    
//    // Bc Bottom
//    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcBottom);
//    
//    // Bc Right
//    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
//    cmesh->InsertMaterialObject(bcRight);
//    
//    // Bc Top
//    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcTop);
//    
//    // Bc Left
//    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcLeft);
//    
//    // Setando Hdiv
//    cmesh->SetDimModel(2);
//    cmesh->SetDefaultOrder(qorder);
//    cmesh->SetAllCreateFunctionsHDiv();
//    
//    
//    cmesh->AutoBuild();
//    
//    
//#ifdef PZDEBUG
//    std::ofstream out("cmeshFlux.txt");
//    cmesh->Print(out);
//#endif
//    
//    return cmesh;
//}
//
//TPZCompMesh * TPZWellAnalysis::CmeshPressure(int porder)
//{
//    
//    int dim = 2;
//    int ilayer = 0;
//    int RockId = fLayers[ilayer]->GetMatIDs()[0];
//    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
//    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
//    int topId = fLayers[ilayer]->GetMatIDs()[3];
//    int leftId = fLayers[ilayer]->GetMatIDs()[4];
//    
//    const int typeFlux = 0, typePressure = 1;
//    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
//    
//    // Malha computacional
//    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
//    
//    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
//    cmesh->InsertMaterialObject(mat);
//    
//    // Bc Bottom
//    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcBottom);
//    
//    // Bc Right
//    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
//    cmesh->InsertMaterialObject(bcRight);
//    
//    // Bc Top
//    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcTop);
//    
//    // Bc Left
//    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcLeft);
//    
//    // Setando L2
//    cmesh->SetDimModel(dim);
//    cmesh->SetDefaultOrder(porder);
//    
//    cmesh->SetAllCreateFunctionsContinuous();
//    cmesh->ApproxSpace().CreateDisconnectedElements(false);
//    //    cmesh->SetAllCreateFunctionsDiscontinuous();
//    cmesh->AutoBuild();
//    
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    int ncon = cmesh->NConnects();
//    for(int i=0; i<ncon; i++)
//    {
//        TPZConnect &newnod = cmesh->ConnectVec()[i];
//        newnod.SetLagrangeMultiplier(1);
//    }
//    
//    //    int nel = cmesh->NElements();
//    //    for(int i=0; i<nel; i++){
//    //        TPZCompEl *cel = cmesh->ElementVec()[i];
//    //        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//    //        celdisc->SetConstC(1.);
//    //        celdisc->SetCenterPoint(0, 0.);
//    //        celdisc->SetCenterPoint(1, 0.);
//    //        celdisc->SetCenterPoint(2, 0.);
//    //        celdisc->SetFalseUseQsiEta();
//    //    }
//    
//    
//    
//    
//#ifdef PZDEBUG
//    std::ofstream out("cmeshPress.txt");
//    cmesh->Print(out);
//#endif
//    
//    return cmesh;
//    
//}
//
//TPZCompMesh * TPZWellAnalysis::CmeshSw(int Sworder)
//{
//    
//    int dim = 2;
//    int ilayer = 0;
//    int RockId = fLayers[ilayer]->GetMatIDs()[0];
//    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
//    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
//    int topId = fLayers[ilayer]->GetMatIDs()[3];
//    int leftId = fLayers[ilayer]->GetMatIDs()[4];
//    
//    const int typeFlux = 0, typePressure = 1;
//    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
//    
//    // Malha computacional
//    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
//    
//    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
//    cmesh->InsertMaterialObject(mat);
//    
//    // Bc Bottom
//    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcBottom);
//    
//    // Bc Right
//    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
//    cmesh->InsertMaterialObject(bcRight);
//    
//    // Bc Top
//    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcTop);
//    
//    // Bc Left
//    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcLeft);
//    
//    // Setando L2
//    cmesh->SetDimModel(dim);
//    cmesh->SetDefaultOrder(Sworder);
//    
//    cmesh->SetAllCreateFunctionsDiscontinuous();
//    cmesh->AutoBuild();
//    
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    
//#ifdef PZDEBUG
//    std::ofstream out("cmeshSw.txt");
//    cmesh->Print(out);
//#endif
//    
//    return cmesh;
//    
//}
//
//TPZCompMesh * TPZWellAnalysis::CmeshSo(int Soorder)
//{
//    
//    int dim = 2;
//    int ilayer = 0;
//    int RockId = fLayers[ilayer]->GetMatIDs()[0];
//    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
//    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
//    int topId = fLayers[ilayer]->GetMatIDs()[3];
//    int leftId = fLayers[ilayer]->GetMatIDs()[4];
//    
//    const int typeFlux = 0, typePressure = 1;
//    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
//    
//    // Malha computacional
//    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
//    
//    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
//    cmesh->InsertMaterialObject(mat);
//    
//    // Bc Bottom
//    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcBottom);
//    
//    // Bc Right
//    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
//    cmesh->InsertMaterialObject(bcRight);
//    
//    // Bc Top
//    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcTop);
//    
//    // Bc Left
//    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
//    cmesh->InsertMaterialObject(bcLeft);
//    
//    // Setando L2
//    cmesh->SetDimModel(dim);
//    cmesh->SetDefaultOrder(Soorder);
//    
//    cmesh->SetAllCreateFunctionsDiscontinuous();
//    cmesh->AutoBuild();
//    
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    
//#ifdef PZDEBUG
//    std::ofstream out("cmeshSo.txt");
//    cmesh->Print(out);
//#endif
//    
//    return cmesh;
//}
//
//
//void TPZWellAnalysis::ReadGeoMesh(std::string GridFileName)
//{
//    TPZReadGIDGrid GeometryInfo;
//    GeometryInfo.SetfDimensionlessL(1.0);
//    fgmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
//    fgmesh->SetDimension(3);
//}
//
//void TPZWellAnalysis::CreatedGeoMesh()
//{
//    
//    
//    int64_t Qnodes = 4;
//    int ilayer = 0;
//    
//    TPZGeoMesh *gmesh= new TPZGeoMesh;
//    
//    gmesh->SetMaxNodeId(Qnodes-1);
//    gmesh->NodeVec().Resize(Qnodes);
//    TPZVec<TPZGeoNode> Node(Qnodes);
//    
//    TPZVec <int64_t> TopolQuad(4);
//    TPZVec <int64_t> TopolLine(2);
//    REAL r     = fLayers[ilayer]->Layerr();
//    REAL rw    = fLayers[ilayer]->Layerrw();
//    REAL h     = fLayers[ilayer]->Layerh();
//    REAL top   = fLayers[ilayer]->LayerTop();
//    
//    int RockId = fLayers[ilayer]->GetMatIDs()[0];
//    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
//    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
//    int topId = fLayers[ilayer]->GetMatIDs()[3];
//    int leftId = fLayers[ilayer]->GetMatIDs()[4];
//    
//    // Nodes
//    int64_t id = 0;
//    
//    Node[id].SetNodeId(id);
//    Node[id].SetCoord(0 ,  rw);         //coord r
//    Node[id].SetCoord(1 , top - h);     //coord z
//    gmesh->NodeVec()[id] = Node[id];
//    id++;
//    
//    Node[id].SetNodeId(id);
//    Node[id].SetCoord(0 ,  rw + r);         //coord r
//    Node[id].SetCoord(1 , top - h);     //coord z
//    gmesh->NodeVec()[id] = Node[id];
//    id++;
//    
//    Node[id].SetNodeId(id);
//    Node[id].SetCoord(0 ,  rw + r);         //coord r
//    Node[id].SetCoord(1 ,  top);     //coord z
//    gmesh->NodeVec()[id] = Node[id];
//    id++;
//    
//    Node[id].SetNodeId(id);
//    Node[id].SetCoord(0 , rw);         //coord r
//    Node[id].SetCoord(1 , top);     //coord z
//    gmesh->NodeVec()[id] = Node[id];
//    id++;
//    
//    
//    //  Geometric Elements
//    int elid = 0;
//    
//    TopolLine[0] = 0;
//    TopolLine[1] = 1;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,bottomId,*gmesh);
//    id++;
//    
//    TopolLine[0] = 1;
//    TopolLine[1] = 2;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,rigthId,*gmesh);
//    id++;
//    
//    TopolLine[0] = 2;
//    TopolLine[1] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,topId,*gmesh);
//    id++;
//    
//    TopolLine[0] = 0;
//    TopolLine[1] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,leftId,*gmesh);
//    id++;
//    
//    
//    TopolQuad[0] = 0;
//    TopolQuad[1] = 1;
//    TopolQuad[2] = 2;
//    TopolQuad[3] = 3;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elid,TopolQuad,RockId,*gmesh);
//    
//    
//    
//    
//    gmesh->BuildConnectivity();
//    fgmesh = gmesh;
//    
//}
//
//void TPZWellAnalysis::Parametricfunction(const TPZVec<REAL> &par, TPZVec<REAL> &X)
//{
//    X[0] = par[0];//cos(par[0]);
//    X[1] = 0.0;//sin(par[0]);
//    X[2] = 0.0;
//}
//
//void TPZWellAnalysis::Parametricfunction2(const TPZVec<REAL> &par, TPZVec<REAL> &X)
//{
//    X[0] = 0.0;//par[0];
//    X[1] = par[0];
//    X[2] = 0.0;
//}
//
//void TPZWellAnalysis::GeometryLine(int nx, int ny)
//{
//    REAL t=0.0;
//    REAL dt = 50.0;
//    int n = nx;
//    
//    // Creating a 0D element to be extruded
//    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
//    GeoMesh1->NodeVec().Resize(1);
//    TPZGeoNode Node;
//    TPZVec<REAL> coors(3,0.0);
//    Node.SetCoord(coors);
//    Node.SetNodeId(0);
//    GeoMesh1->NodeVec()[0]=Node;
//    
//    TPZVec<int64_t> Topology(1,0);
//    int elid=0;
//    int matid=1;
//    
//    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
//    GeoMesh1->BuildConnectivity();
//    GeoMesh1->SetDimension(0);
//    
//    TPZHierarquicalGrid CreateGridFrom(GeoMesh1);
//    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(Parametricfunction);
//    CreateGridFrom.SetParametricFunction(ParFunc);
//    CreateGridFrom.SetFrontBackMatId(5,3);
//    
//    
//    // Computing Mesh extruded along the parametric curve Parametricfunction
//    TPZGeoMesh * GeoMesh2 = CreateGridFrom.ComputeExtrusion(t, dt, n);
//    
//    TPZHierarquicalGrid CreateGridFrom2(GeoMesh2);
//    TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(Parametricfunction2);
//    CreateGridFrom2.SetParametricFunction(ParFunc2);
//    CreateGridFrom2.SetFrontBackMatId(2,4);
//    
//    dt = 100.0;
//    n = ny;
//    
//    
//    // Computing Mesh extruded along the parametric curve Parametricfunction2
//    fgmesh = CreateGridFrom2.ComputeExtrusion(t, dt, n);
//    
//    TPZCheckGeom * GeometryTest = new TPZCheckGeom;
//    int isBadMesh = 0;
//    bool CheckGeometry = false;
//    
//    
//}
//
//void TPZWellAnalysis::PrintGeoMesh()
//{
//    
//    
//#ifdef PZDEBUG
//    //  Print Geometrical Base Mesh
//    std::ofstream argument("GeometicMesh.txt");
//    fgmesh->Print(argument);
//    std::ofstream Dummyfile("GeometricMesh.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(fgmesh,Dummyfile, true);
//    
//#endif
//}
//
//void TPZWellAnalysis::RotateGeomesh(REAL CounterClockwiseAngle)
//{
//    REAL theta = CounterClockwiseAngle;
//    // It represents a 3D rotation around the z axis.
//    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
//    RotationMatrix(0,0) =   +cos(theta);
//    RotationMatrix(0,1) =   -sin(theta);
//    RotationMatrix(1,0) =   +sin(theta);
//    RotationMatrix(1,1) =   +cos(theta);
//    RotationMatrix(2,2) = 1.0;
//    TPZVec<STATE> iCoords(3,0.0);
//    TPZVec<STATE> iCoordsRotated(3,0.0);
//    
//    RotationMatrix.Print("Rotation = ");
//    
//    int NumberofGeoNodes = fgmesh->NNodes();
//    for (int inode = 0; inode < NumberofGeoNodes; inode++)
//    {
//        TPZGeoNode GeoNode = fgmesh->NodeVec()[inode];
//        GeoNode.GetCoordinates(iCoords);
//        // Apply rotation
//        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
//        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
//        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
//        GeoNode.SetCoord(iCoordsRotated);
//        fgmesh->NodeVec()[inode] = GeoNode;
//    }
//}
//
//void TPZWellAnalysis::UniformRefinement(int nh)
//{
//    for ( int ref = 0; ref < nh; ref++ ){
//        TPZVec<TPZGeoEl *> filhos;
//        int64_t n = fgmesh->NElements();
//        for ( int64_t i = 0; i < n; i++ ){
//            TPZGeoEl * gel = fgmesh->ElementVec() [i];
//            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
//        }//for i
//    }//ref
//    fgmesh->BuildConnectivity();
//}
//
//void TPZWellAnalysis::UniformRefinement(int nh, std::set<int> &MatToRef)
//{
//    for ( int ref = 0; ref < nh; ref++ ){
//        TPZVec<TPZGeoEl *> filhos;
//        int64_t n = fgmesh->NElements();
//        for ( int64_t i = 0; i < n; i++ ){
//            TPZGeoEl * gel = fgmesh->ElementVec() [i];
//            if(!gel){continue;}
//            //            int reflevel = gel->Level();
//            //            if (reflevel == ref + 1) {
//            //                continue;
//            //            }
//            TPZRefPatternTools::RefineDirectional(gel,MatToRef);
//        }//for i
//    }//ref
//    fgmesh->BuildConnectivity();
//}
//
//void TPZWellAnalysis::UniformRefinement(int nh, int MatId)
//{
//    //    for ( int ref = 0; ref < nh; ref++ ){
//    //        TPZVec<TPZGeoEl *> filhos;
//    //        int64_t n = fgmesh->NElements();
//    //        for ( int64_t i = 0; i < n; i++ ){
//    //            TPZGeoEl * gel = fgmesh->ElementVec() [i];
//    //            if(!gel){continue;}
//    //            if (gel->Dimension() == 1){
//    //                if (gel->MaterialId() == MatId) {
//    //                    gel->Divide(filhos);
//    //                }
//    //
//    //            }
//    //        }//for i
//    //    }//ref
//    
//    ///Refinamento
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    
//    
//    for (int idivide = 0; idivide < nh; idivide++){
//        const int nels = fgmesh->NElements();
//        TPZVec< TPZGeoEl * > allEls(nels);
//        for(int iel = 0; iel < nels; iel++){
//            allEls[iel] = fgmesh->ElementVec()[iel];
//        }
//        
//        for(int iel = 0; iel < nels; iel++){
//            TPZGeoEl * gel = allEls[iel];
//            if(!gel) continue;
//            if(gel->HasSubElement()) continue;
//            int nnodes = gel->NNodes();
//            int found = -1;
//            for(int in = 0; in < nnodes; in++){
//                if(gel->NodePtr(in)->Id() == MatId){
//                    found = in;
//                    break;
//                }
//            }///for in
//            if(found == -1) continue;
//            
//            MElementType gelT = gel->Type();
//            TPZAutoPointer<TPZRefPattern> uniform = gRefDBase.GetUniformRefPattern(gelT);
//            if(!uniform){
//                DebugStop();
//            }
//            gel->SetRefPattern(uniform);
//            TPZVec<TPZGeoEl*> filhos;
//            gel->Divide(filhos);
//            
//        }///for iel
//    }//idivide
//    
//    fgmesh->BuildConnectivity();
//    
//}
//
//////refinamento uniforme em direcao ao no
////void DirectionalRef(int nh, int MatId){
////
////    ///Refinamento
////    gRefDBase.InitializeUniformRefPattern(EOned);
////    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
////
////
////    for (int idivide = 0; idivide < nh; idivide++){
////        const int nels = fgmesh->NElements();
////        TPZVec< TPZGeoEl * > allEls(nels);
////        for(int iel = 0; iel < nels; iel++){
////            allEls[iel] = gmesh->ElementVec()[iel];
////        }
////
////        for(int iel = 0; iel < nels; iel++){
////            TPZGeoEl * gel = allEls[iel];
////            if(!gel) continue;
////            if(gel->HasSubElement()) continue;
////            int nnodes = gel->NNodes();
////            int found = -1;
////            for(int in = 0; in < nnodes; in++){
////                if(gel->NodePtr(in)->Id() == nodeAtOriginId){
////                    found = in;
////                    break;
////                }
////            }///for in
////            if(found == -1) continue;
////
////            MElementType gelT = gel->Type();
////            TPZAutoPointer<TPZRefPattern> uniform = gRefDBase.GetUniformRefPattern(gelT);
////            if(!uniform){
////                DebugStop();
////            }
////            gel->SetRefPattern(uniform);
////            TPZVec<TPZGeoEl*> filhos;
////            gel->Divide(filhos);
////
////        }///for iel
////    }//idivide
////
////    gmesh->BuildConnectivity();
////
////#ifdef LOG4CXX
////    if (logger.isDebugEnabled())
////    {
////        std::stringstream sout;
////        sout<<"gmesh depois de refinar direcionalmente\n";
////        gmesh->Print(sout);
////        LOGPZ_DEBUG(logger, sout.str());
////    }
////#endif
////
////}///void
//
//
//
//void TPZWellAnalysis::PostProcessVTK(TPZAnalysis *an)
//{
//    const int dim = 2;
//    int div = fSimulationData->GetHPostrefinement();
//    TPZStack<std::string> scalnames, vecnames;
//    std::string plotfile;
//    if (fSimulationData->GetIsH1approx()) {
//        plotfile = "2DH1Darcy.vtk";
//    }
//    else{
//        plotfile = "2DMixedDarcy.vtk";
//    }
//    
//    scalnames.Push("WeightedPressure");
//    scalnames.Push("WaterSaturation");
//    scalnames.Push("OilSaturation");
//    scalnames.Push("WaterDensity");
//    scalnames.Push("OilDensity");
//    scalnames.Push("Porosity");
//    scalnames.Push("DivOfBulkVeclocity");
//    scalnames.Push("ExactSaturation");
//    vecnames.Push("BulkVelocity");
//    an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
//    an->PostProcess(div);
//}
//
//void TPZWellAnalysis::Ffunction(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
//{
//    
//    ff[0] = 0.0*0.0000001;
//}
//
//TPZFMatrix<STATE> * TPZWellAnalysis::ComputeInverse()
//{
//    int neq = fcmeshdarcy->NEquations();
//    TPZFMatrix<STATE> * PreInverse =  new TPZFMatrix<STATE> (neq,neq,0.0);
//    TPZFStructMatrix skyl(fcmeshdarcy);
//    std::set<int> matids; // to be computed
//    matids.insert(1);
//    matids.insert(2);
//    matids.insert(3);
//    matids.insert(4);
//    matids.insert(5);
//    skyl.SetMaterialIds(matids);
//    TPZFMatrix<STATE> rhsfrac;
//    TPZFMatrix<STATE> Identity;
//    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
//    TPZAutoPointer<TPZMatrix<STATE> > MatG = skyl.CreateAssemble(rhsfrac, gui);
//    TPZFMatrix<STATE> oldmat = *MatG.operator->();
//    oldmat.Inverse( * PreInverse);
//    oldmat.Multiply(*PreInverse, Identity);
//    
//#ifdef PZDEBUG
//#ifdef LOG4CXX
//    if(logger.isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout << "Is decomposed=  " << MatG->IsDecomposed() << std::endl;
//        oldmat.Print("oldmat = ", sout,EMathematicaInput);
//        PreInverse->Print("PreInverse = ", sout,EMathematicaInput);
//        Identity.Print("Identity = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
//#endif
//    
//    return PreInverse;
//    
//}
//
//
//
///**
// * Computes the saturation at shock using the Welge method
// */
//REAL TPZDarcyAnalysis::SwatShock(REAL epsilon, REAL ds){
//    REAL muw= 0.001;
//    REAL muo= 0.001;
//    REAL value = 0.0;
//    int npoints  = int(1.0/ds);
//    int pos = 0;
//    TPZFMatrix<REAL> points(npoints,2,0.0);
//    TPZManVector<REAL> secants(npoints,0.0);
//    REAL S = 0.0;
//    for (int ip = 0; ip < npoints; ip++) {
//        S = ds*ip;
//        points(ip,0) = S;
//        points(ip,1) = ((S*S*muo)/(muw - 2.0 * S * muw + S * S * (muw + muo)));
//        
//    }
//    
//    for (int ip = 0; ip < npoints - 1; ip++) {
//        secants[ip] = (points(ip + 1,1)-points(0,1))/(points(ip+1,0)-points(0,0));
//        
//    }
//    
//    // computing the max value
//    
//    
//    
//    pos = Extract(epsilon, secants, value);
//    return points(pos,1);
//    
//}
//
///**
// * Extract a value from a given list
// */
//int TPZDarcyAnalysis::Extract(REAL epsilon, TPZManVector<REAL> &list, REAL value){
//    
//    for (int i = 0; i<list.size(); i++) {
//        if (fabs(list[i]-value) <= epsilon) {
//            return i;
//        }
//    }
//    return 0;
//}
//

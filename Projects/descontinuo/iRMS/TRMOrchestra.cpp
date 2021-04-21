//
//  TRMOrchestra.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMOrchestra.h"
#include "TPZMaterial.h"
#include "TPZVecL2.h"
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TRMFlowConstants.h"
#include "pzquad.h"
#include "pzaxestools.h"

#include "tpzintpoints.h"
#include "TPZMatWithMem.h"
#include "TRMMemory.h"
#include "TRMMixedDarcy.h"
#include "pzinterpolationspace.h"
#include "pzysmp.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

TRMOrchestra::TRMOrchestra(){
    
    fgmesh                          = NULL;
    fSpaceGenerator                 = new TRMSpaceOdissey;
    fSimulationData                 = NULL;
    fPrimalMultiphaseAnalysis       = new TRMPrimalMultiphaseAnalysis;
    fMonolithicMultiphaseAnalysis_I = new TRMMonolithicMultiphaseAnalysis;
    fMonolithicMultiphaseAnalysis   = new TRMMonolithicMultiphaseAnalysis;
    fSegregatedAnalysis_I           = new TRMSegregatedAnalysis;
    fSegregatedAnalysis             = new TRMSegregatedAnalysis;
//    fparabolic                      = new TRMFluxPressureAnalysis;
//    fhyperbolic                     = new TRMTransportAnalysis;
    
    fIsMonolithicQ         =  false;
    fIsSegregatedQ         =  false;
    fIsSegregatedwithCGQ   =  false;
    
    
}


TRMOrchestra::~TRMOrchestra(){
    
}

/** @brief Create geometric mesh being used by space odissey */
void TRMOrchestra::BuildGeometry(bool Is3DGeometryQ){
    
    bool IsReservoirBoxQ = false;
    
    if (Is3DGeometryQ) {
        
        int nel_x = 5;
        int nel_y = 1;
        int nel_z = 1;
        
        TPZManVector<REAL,2> dx(2,nel_x), dy(2,nel_y), dz(2,nel_z);
        dx[0] = 1000.0/REAL(nel_x);
        dy[0] = 100.0/REAL(nel_y);
        dz[0] = 20.0/REAL(nel_z);
        
        if (IsReservoirBoxQ) {
            fSpaceGenerator->CreateGeometricBoxMesh(dx, dy, dz);
        }
        else
        {
            std::string dirname = PZSOURCEDIR;
            std::string file;
//            file = dirname + "/Projects/iRMS/Meshes/BarriesGeo.dump";
//            file = dirname + "/Projects/iRMS/Meshes/Ciruclar_ReservoirC.dump";
            file = dirname + "/Projects/iRMS/Meshes/CircularMeshVerticalWellQII.dump";
//            file = dirname + "/Projects/iRMS/Meshes/FiveSpotQ.dump";
//            file = dirname + "/Projects/iRMS/Meshes/FiveSpotBarriesQ.dump";
            fSpaceGenerator->CreateGeometricExtrudedGIDMesh(file, dz);
        }
        
        fSpaceGenerator->Gmesh()->SetDimension(3);
    }
    else{
        
        int nel_x = 8;
        int nel_y = 4;
        
        TPZManVector<REAL,2> dx(2,nel_x), dy(2,nel_y);
        dx[0] = 1000.0/REAL(nel_x);
        dy[0] = 100.0/REAL(nel_y);
        
        if (IsReservoirBoxQ) {
            fSpaceGenerator->CreateGeometricBoxMesh2D(dx, dy);
        }
        else
        {
            std::string dirname = PZSOURCEDIR;
            std::string file;
//            file = dirname + "/Projects/iRMS/Meshes/BarriesGeo.dump";
//            file = dirname + "/Projects/iRMS/Meshes/Ciruclar_ReservoirC.dump";
            file = dirname + "/Projects/iRMS/Meshes/CircularMeshVerticalWellQII.dump";
//            file = dirname + "/Projects/iRMS/Meshes/FiveSpotQ.dump";
//            file = dirname + "/Projects/iRMS/Meshes/FiveSpotBarriesQ.dump";
//            file = dirname + "/Projects/iRMS/Meshes/TwoWellQ.dump";
            fSpaceGenerator->CreateGeometricGIDMesh(file);
        }
        
        fSpaceGenerator->Gmesh()->SetDimension(2);

    }
    

    int ref = 1;
    fSpaceGenerator->UniformRefinement(ref);
//    fSpaceGenerator->UniformRefinement_Around_MaterialId(ref, 11);
//    fSpaceGenerator->UniformRefinement_Around_MaterialId(ref, 12);
    fSpaceGenerator->PrintGeometry();
//    int father_index = 9;
//    fSpaceGenerator->UniformRefinement_at_Father(1, father_index);
//    fSpaceGenerator->PrintGeometry();
}

/** @brief Create geometric mesh being used by space odissey */
void TRMOrchestra::BuildGeometry2(){
    
    std::string dirname = PZSOURCEDIR;
    std::string file;
    file = dirname + "/Projects/iRMS/Meshes/Gmsh/reservoir.msh";
    fSpaceGenerator->CreateGeometricGmshMesh(file);
    
    int ref = 0;
    fSpaceGenerator->UniformRefinement(ref);
    fSpaceGenerator->PrintGeometry();
}

/** @brief Create a primal analysis using space odissey */
void TRMOrchestra::CreateAnalysisPrimal()
{
    
    TPZManVector<REAL,2> dx(2,1), dy(2,1), dz(2,1);
    dx[0] = 100;
    dy[0] = 100;
    dz[0] = 100;
    
    fSpaceGenerator->CreateGeometricBoxMesh(dx, dy, dz);
//    spacegenerator.CreateGeometricReservoirMesh();
    fSpaceGenerator->PrintGeometry();
    fgmesh = fSpaceGenerator->Gmesh();
    fSpaceGenerator->CreateH1Cmesh();
    
    TPZCompMesh * Cmesh = fSpaceGenerator->H1CMesh();
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis * AnalysisPrimal = new TPZAnalysis(Cmesh,mustOptimizeBandwidth);
    int numofThreads = 8;
    
    TPZSkylineStructMatrix skylnsym(Cmesh);
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ECholesky);
    AnalysisPrimal->SetStructuralMatrix(skylnsym);
    AnalysisPrimal->SetSolver(step);;
    AnalysisPrimal->Run();
    std::cout << "Primal dof: " << AnalysisPrimal->Rhs().Rows() << std::endl;
    
    const int dim = 3;
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "PrimalDarcy.vtk";
    scalnames.Push("Pressure");
    vecnames.Push("MinusKGradU");
    AnalysisPrimal->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    AnalysisPrimal->PostProcess(div);
    
}

/** @brief Create a dual analysis using space odissey */
void TRMOrchestra::CreateAnalysisDualonBox(bool IsInitialQ)
{

//    this->BuildGeometry(false);
    this->BuildGeometry2();
    
    fSimulationData->SetInitialStateQ(IsInitialQ);
    TRMFluxPressureAnalysis * parabolic = new TRMFluxPressureAnalysis;
    TRMTransportAnalysis * hyperbolic = new TRMTransportAnalysis;
    
#ifdef PZDEBUG
    if (!fSpaceGenerator->Gmesh()) {
        std::cout << "iRMS:: Call BuildGeometry " << std::endl;
        DebugStop();
    }
    
    TPZCheckGeom check(fSpaceGenerator->Gmesh());
    if (check.PerformCheck() != 0){
        DebugStop();
    }
    
#endif
    
    fSpaceGenerator->SetDefaultPOrder(1);
    fSpaceGenerator->SetDefaultSOrder(0);

    bool UseMHMQ = false;
    
    if(UseMHMQ){
        int skeleton_id = 0;
        fSpaceGenerator->InsertSkeletonInterfaces(skeleton_id); // @omar:: Primitive use of the mhm capabilities
        fSpaceGenerator->BuildMHM_Mesh();
    }
    else{
        fSpaceGenerator->BuildMixed_Mesh();
    }
    
    if(fSimulationData->IsTwoPhaseQ()){
        fSpaceGenerator->CreateAlphaTransportMesh();
        hyperbolic->Meshvec().Resize(1);
        hyperbolic->Meshvec()[0] = fSpaceGenerator->AlphaSaturationMesh();
        fSpaceGenerator->CreateTransportMesh();        
    }
    
    if(fSimulationData->IsThreePhaseQ()){
        fSpaceGenerator->CreateAlphaTransportMesh();
        fSpaceGenerator->CreateBetaTransportMesh();
        hyperbolic->Meshvec().Resize(2);
        hyperbolic->Meshvec()[0] = fSpaceGenerator->AlphaSaturationMesh();
        hyperbolic->Meshvec()[1] = fSpaceGenerator->BetaSaturationMesh();
        fSpaceGenerator->CreateTransportMesh();
    }
        
    int numofThreads_p = 16;
    bool mustOptimizeBandwidth_parabolic = true;
    
    /////////////////////////////////////////// No subtructures ///////////////////////////////////////////
    // Analysis for parabolic part
    parabolic->Meshvec()[0] = fSpaceGenerator->FluxCmesh();
    parabolic->Meshvec()[1] = fSpaceGenerator->PressureCmesh();
    parabolic->SetCompMesh(fSpaceGenerator->MixedFluxPressureCmesh(), mustOptimizeBandwidth_parabolic);
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat_p(fSpaceGenerator->MixedFluxPressureCmesh());
//    strmat_p.SetDecomposeType(ELDLt);

    TPZSSpStructMatrix strmat_p(fSpaceGenerator->MixedFluxPressureCmesh());
    
    TPZStepSolver<STATE> step_p;
    step_p.SetDirect(ELDLt);
    strmat_p.SetNumThreads(numofThreads_p);
    parabolic->SetStructuralMatrix(strmat_p);
    parabolic->SetSolver(step_p);
    parabolic->AdjustVectors();
    parabolic->SetSimulationData(fSimulationData);
    
    std::cout << "ndof parabolic = " << parabolic->Solution().Rows() << std::endl;
    
    if (fSimulationData->IsTwoPhaseQ() || fSimulationData->IsThreePhaseQ()) {
    
        // Analysis for hyperbolic part
        int numofThreads_t = 16;
        bool mustOptimizeBandwidth_hyperbolic = true;
        hyperbolic->SetCompMesh(fSpaceGenerator->TransportMesh(), mustOptimizeBandwidth_hyperbolic);
//        TPZSkylineNSymStructMatrix strmat_t(fSpaceGenerator->TransportMesh());
//        strmat_t.SetNumThreads(numofThreads_t);
//        step_t.SetDirect(ELU);
//        hyperbolic->SetStructuralMatrix(strmat_t);
//        hyperbolic->SetSolver(step_t);
//        hyperbolic->AdjustVectors();
//        hyperbolic->SetSimulationData(fSimulationData);
//        hyperbolic->FilterEquations();
        
        TPZSpStructMatrix strmat_t(fSpaceGenerator->TransportMesh());
        TPZStepSolver<STATE> step_t;

        const int64_t numiterations = 20;
        const REAL tol = 1.0e-8;
        step_t.SetJacobi(numiterations, tol, 1);

        strmat_t.SetNumThreads(numofThreads_t);
        hyperbolic->SetStructuralMatrix(strmat_t);
        hyperbolic->SetSolver(step_t);
        hyperbolic->AdjustVectors();
        hyperbolic->SetSimulationData(fSimulationData);
        hyperbolic->FilterEquations();
        
        
//        TPZAutoPointer<TPZMatrix<STATE> > nsyma = strmat_t.Create();
//        TPZAutoPointer<TPZMatrix<STATE> > nsymaClone = nsyma->Clone();
//        
//        TPZStepSolver<STATE> *stepre = new TPZStepSolver<STATE>(nsymaClone);
//        TPZStepSolver<STATE> *stepGMRES = new TPZStepSolver<STATE>(nsyma);
//        
//        const int64_t numiterations = 20;
//        const REAL tol = 1.0e-6;
//        step_t.SetJacobi(numiterations, tol, 1);
//        stepre->SetReferenceMatrix(nsyma);
//        stepGMRES->SetGMRES(10, 20, *stepre, 1.0e-10, 0);
//        strmat_t.SetNumThreads(numofThreads_t);
//        hyperbolic->SetStructuralMatrix(strmat_t);
//        hyperbolic->SetSolver(*stepGMRES);
//        hyperbolic->AdjustVectors();
//        hyperbolic->SetSimulationData(fSimulationData);
//        hyperbolic->FilterEquations();

        
        
        std::cout << "ndof hyperbolic = " << hyperbolic->Solution().Rows() << std::endl;
    }

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Transfer object
    TRMBuildTransfers * Transfer = new TRMBuildTransfers;
    Transfer->SetSimulationData(fSimulationData);
    Transfer->Fill_u_To_Mixed(parabolic->Mesh(), 0);
    Transfer->Fill_p_To_Mixed(parabolic->Mesh(), 1);
    
    if(fSimulationData->IsOnePhaseQ()){
        Transfer->FillComputationalElPairs(parabolic->Mesh(),parabolic->Mesh());
    }
    
    if(fSimulationData->IsTwoPhaseQ()){
        Transfer->FillComputationalElPairs(parabolic->Mesh(),hyperbolic->Mesh());
        Transfer->Fill_s_To_Transport(hyperbolic->Mesh(), 0);
        Transfer->ComputeLeftRight(hyperbolic->Mesh());
        Transfer->Fill_un_To_Transport(parabolic->Mesh(),hyperbolic->Mesh(),true);
        Transfer->Fill_un_To_Transport(parabolic->Mesh(),hyperbolic->Mesh(),false);
    }
    
    if(fSimulationData->IsThreePhaseQ()){
        Transfer->FillComputationalElPairs(parabolic->Mesh(),hyperbolic->Mesh());
        Transfer->Fill_s_To_Transport(hyperbolic->Mesh(), 0);
        Transfer->Fill_s_To_Transport(hyperbolic->Mesh(), 1);
        Transfer->ComputeLeftRight(hyperbolic->Mesh());
        Transfer->Fill_un_To_Transport(parabolic->Mesh(),hyperbolic->Mesh(),true);
        Transfer->Fill_un_To_Transport(parabolic->Mesh(),hyperbolic->Mesh(),false);
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    parabolic->SetTransfer(Transfer);
    if (fSimulationData->IsTwoPhaseQ() || fSimulationData->IsThreePhaseQ()) {
        hyperbolic->SetTransfer(Transfer);
    }
    
    TRMSegregatedAnalysis * segregated = new TRMSegregatedAnalysis;
    segregated->SetTransfer(Transfer);
    segregated->SetSimulationData(fSimulationData);
    segregated->SetParabolic(parabolic);
    segregated->SetHyperbolic(hyperbolic);
    
    if (IsInitialQ) {
        fSegregatedAnalysis_I = segregated;
    }
    else{
        fSegregatedAnalysis   =  segregated;
    }
    
    
    
}

/** @brief Create a monolithic dual analysis on box geometry using space odissey */
void TRMOrchestra::CreateMonolithicAnalysis(bool IsInitialQ){
    
    this->BuildGeometry(false);
    
    fSimulationData->SetInitialStateQ(IsInitialQ);
    
    TRMMonolithicMultiphaseAnalysis * mono_analysis = new TRMMonolithicMultiphaseAnalysis;
    
#ifdef PZDEBUG
    if (!fSpaceGenerator->Gmesh()) {
        std::cout << "iRMS:: Call BuildGeometry " << std::endl;
        DebugStop();
    }
    fSpaceGenerator->PrintGeometry();
#endif
    
    fSpaceGenerator->SetDefaultUOrder(2);
    fSpaceGenerator->SetDefaultPOrder(2);
    fSpaceGenerator->SetDefaultSOrder(0);

    // Structure for one-phase flow
    if(fSimulationData->IsOnePhaseQ()){
        
        fSpaceGenerator->CreateBiotCmesh();
        fSpaceGenerator->CreateFluxCmesh();
        fSpaceGenerator->CreatePressureCmesh();
        fSpaceGenerator->CreateMultiphaseCmesh();

        mono_analysis->Meshvec()[0] = fSpaceGenerator->BiotCMesh();
        mono_analysis->Meshvec()[1] = fSpaceGenerator->FluxCmesh();
        mono_analysis->Meshvec()[2] = fSpaceGenerator->PressureCmesh();
        
    }
    
    if(fSimulationData->IsTwoPhaseQ()){
        mono_analysis->Meshvec().Resize(4);
        fSpaceGenerator->CreateBiotCmesh();
        fSpaceGenerator->CreateFluxCmesh();
        fSpaceGenerator->CreatePressureCmesh();
        fSpaceGenerator->CreateAlphaTransportMesh();
        fSpaceGenerator->CreateMultiphaseCmesh();
        fSpaceGenerator->CreateInterfacesInside(fSpaceGenerator->MonolithicMultiphaseCmesh());

        mono_analysis->Meshvec()[0] = fSpaceGenerator->BiotCMesh();
        mono_analysis->Meshvec()[1] = fSpaceGenerator->FluxCmesh();
        mono_analysis->Meshvec()[2] = fSpaceGenerator->PressureCmesh();
        mono_analysis->Meshvec()[3] = fSpaceGenerator->AlphaSaturationMesh();
        
    }
    
    
    if(fSimulationData->IsThreePhaseQ()){
        
        mono_analysis->Meshvec().Resize(5);
        fSpaceGenerator->CreateBiotCmesh();        
        fSpaceGenerator->CreateFluxCmesh();
        fSpaceGenerator->CreatePressureCmesh();
        fSpaceGenerator->CreateAlphaTransportMesh();
        fSpaceGenerator->CreateBetaTransportMesh();
        fSpaceGenerator->CreateMultiphaseCmesh();
        fSpaceGenerator->CreateInterfacesInside(fSpaceGenerator->MonolithicMultiphaseCmesh());
        
        mono_analysis->Meshvec()[0] = fSpaceGenerator->BiotCMesh();
        mono_analysis->Meshvec()[1] = fSpaceGenerator->FluxCmesh();
        mono_analysis->Meshvec()[2] = fSpaceGenerator->PressureCmesh();
        mono_analysis->Meshvec()[3] = fSpaceGenerator->AlphaSaturationMesh();
        mono_analysis->Meshvec()[4] = fSpaceGenerator->BetaSaturationMesh();
        
    }
    
    bool mustOptimizeBandwidth = true;
    mono_analysis->SetCompMesh(fSpaceGenerator->MonolithicMultiphaseCmesh(), mustOptimizeBandwidth);
    std::cout << "Total dof: " << mono_analysis->Solution().Rows() << std::endl;
    
    // Use this matrix for a linear tracer
    TPZSkylineNSymStructMatrix skyns_mat(fSpaceGenerator->MonolithicMultiphaseCmesh());
    TPZStepSolver<STATE> step;
    int numofThreads = 0;
    skyns_mat.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    mono_analysis->SetStructuralMatrix(skyns_mat);
    mono_analysis->SetSolver(step);
    mono_analysis->SetSimulationData(fSimulationData);
    mono_analysis->AdjustVectors();
    
    if (IsInitialQ) {
        fMonolithicMultiphaseAnalysis_I = mono_analysis;
    }
    else{
        fMonolithicMultiphaseAnalysis   =  mono_analysis;
    }
    
}

/** @brief Run the static problem over a single large time step */
void TRMOrchestra::RunStaticProblem(){
    
    std::cout<< "iRMS:: Finding Initial State" << std::endl;
    
    int n = 2;
    bool draw_mixed_mapQ = false;
    REAL dt = fSimulationData->dt();
    fSimulationData->Setdt(1.0e10);
    
    for (int i = 0; i < n; i++) {
        if (IsMonolithicQ()) {
            fMonolithicMultiphaseAnalysis_I->ExcecuteOneStep();
            fMonolithicMultiphaseAnalysis_I->PostProcessStep();
        }
        
        if (IsSegregatedQ()) {

            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            fSegregatedAnalysis_I->ExcecuteOneStep();
            fSegregatedAnalysis_I->PostProcessStep(draw_mixed_mapQ);
            
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            
#ifdef USING_BOOST
            std::cout  << "iRMS:: OneStep execution time = " << (t2-t1) << std::endl;
#endif
            
        }
        
    }
    
    fSimulationData->Setdt(dt);
    
    
    // @omar:: initial conditions
    
    if(fSimulationData->IsOnePhaseQ()){
        return;
    }
    
//    int neq_sa = fSegregatedAnalysis_I->Hyperbolic()->Meshvec()[0]->Solution().Rows();
//    int neq_sb = fSegregatedAnalysis_I->Hyperbolic()->Meshvec()[1]->Solution().Rows();
//    for (int i = 0; i < neq_sb; i++) {
//        fSegregatedAnalysis_I->Hyperbolic()->Meshvec()[0]->Solution()(i,0) = 0.0;
//        fSegregatedAnalysis_I->Hyperbolic()->Meshvec()[1]->Solution()(i,0) = 0.8;
//    }

    TPZBuildMultiphysicsMesh::TransferFromMeshes(fSegregatedAnalysis_I->Hyperbolic()->Meshvec(), fSegregatedAnalysis_I->Hyperbolic()->Mesh());
    fSegregatedAnalysis_I->Hyperbolic()->SetX_n(fSegregatedAnalysis_I->Hyperbolic()->Mesh()->Solution());
    
//    fSegregatedAnalysis_I->Hyperbolic()->Mesh()->Solution().Print("ah = ");
//    fSegregatedAnalysis_I->Hyperbolic()->X_n().Print("ah = ");
    
}

/** @brief Run the evolutionary problem for all steps set in the simulation data */
void TRMOrchestra::RunEvolutionaryProblem(){
    
    std::cout<< "iRMS:: Running Evolutionary problem" << std::endl;
    
    if (IsMonolithicQ()) {
        fMonolithicMultiphaseAnalysis->SetX(fMonolithicMultiphaseAnalysis_I->X_n());
        fMonolithicMultiphaseAnalysis->SetX_n(fMonolithicMultiphaseAnalysis_I->X_n());
        fMonolithicMultiphaseAnalysis->LoadSolution(fMonolithicMultiphaseAnalysis_I->X_n());     
    }
    
    if (IsSegregatedQ()) {
        fSegregatedAnalysis->Parabolic()->SetX(fSegregatedAnalysis_I->Parabolic()->X_n());
        fSegregatedAnalysis->Parabolic()->SetX_n(fSegregatedAnalysis_I->Parabolic()->X_n());
        fSegregatedAnalysis->Parabolic()->LoadSolution(fSegregatedAnalysis_I->Parabolic()->X_n());
        
        fSegregatedAnalysis->Hyperbolic()->SetX(fSegregatedAnalysis_I->Hyperbolic()->X_n());
        fSegregatedAnalysis->Hyperbolic()->SetX_n(fSegregatedAnalysis_I->Hyperbolic()->X_n());
        fSegregatedAnalysis->Hyperbolic()->LoadSolution(fSegregatedAnalysis_I->Hyperbolic()->X_n());
        
    }
    
    // Evolutionary problem
    bool draw_mixed_mapQ = false;
    int n = fSimulationData->n_steps();
    REAL time = 0.0;
    fSimulationData->SetTime(time);
    bool MustReportQ = false;
    
    if (fSimulationData->ReportingTimes().size() == 0) {
        return;
    }
    
    time = fSimulationData->t();
    MustReportQ = MustResporTimeQ(time,draw_mixed_mapQ);
    
    if (MustReportQ) {
#ifdef USING_BOOST
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
        std::cout << "iRMS:: Reporting at: " << fSimulationData->t()/86400.0 << "; (day): " << std::endl;
        fSegregatedAnalysis->PostProcessStep(draw_mixed_mapQ);
        
#ifdef USING_BOOST
        boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
        
#ifdef USING_BOOST
        std::cout  << "iRMS:: PostProcess execution time = " << (t2-t1) << std::endl;
#endif
        
    }
    
    for (int i = 0; i < n; i++) {
        
        if (IsMonolithicQ()) {
            
            if (MustReportQ) {
#ifdef USING_BOOST
                boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
                std::cout << "iRMS:: Reporting at: " << fSimulationData->t()/86400.0 << "; (day): " << std::endl;
                fMonolithicMultiphaseAnalysis->PostProcessStep();
#ifdef USING_BOOST
                boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
                
#ifdef USING_BOOST
                std::cout  << "iRMS:: PostProcess execution time = " << (t2-t1) << std::endl;
#endif
            }
            
            if (fSimulationData->ReportingTimes().size() == 0) {
                return;
            }
            
            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            fMonolithicMultiphaseAnalysis->ExcecuteOneStep();
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            
#ifdef USING_BOOST
            std::cout  << "iRMS:: OneStep execution time = " << (t2-t1) << std::endl;
#endif

        }
        
        if (IsSegregatedQ()) {
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            fSegregatedAnalysis->ExcecuteOneStep();
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            
#ifdef USING_BOOST
            std::cout  << "iRMS:: OneStep execution time = " << (t2-t1) << std::endl;
#endif

            time = fSimulationData->t();
            MustReportQ = MustResporTimeQ(time,draw_mixed_mapQ);
            
            if (MustReportQ) {

                
#ifdef USING_BOOST
                boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
                
                std::cout << "iRMS:: Reporting at: " << fSimulationData->t()/86400.0 << "; (day): " << std::endl;
                fSegregatedAnalysis->PostProcessStep(draw_mixed_mapQ);
                
#ifdef USING_BOOST
                boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
                
#ifdef USING_BOOST
                std::cout  << "iRMS:: PostProcess execution time = " << (t2-t1) << std::endl;
#endif
            }
            
            if (fSimulationData->ReportingTimes().size() == 0) {
                return;
            }
        }

    }
}

/** @brief Must report time */
bool TRMOrchestra::MustResporTimeQ(REAL time, bool & draw_mixed_mapQ){
    
    int index = fSimulationData->ReportingTimes().size();
    REAL time_r = fSimulationData->ReportingTimes()[index-1];
    draw_mixed_mapQ = fSimulationData->ReportingTimesMixedQ()[index-1];
    REAL dt = fSimulationData->dt();
    REAL t_range = dt;
    REAL deltat = time-time_r;
    REAL tolerance = 1.0e-3;
    
    if(index > 1){
        REAL c_dt_max = fSimulationData->ReportingTimes()[index-2] - time_r;
        if(dt > c_dt_max){
            fSimulationData->Setdt(c_dt_max);
            std::cout << "Segregated:: Reporting time:: Set time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
        }
    }
    
    if(fabs(deltat) <= tolerance){
        fSimulationData->ReportingTimes().Pop();
        return true;
    }
    
    if (fabs(deltat) < t_range) {
        fSimulationData->Setdt(fabs(deltat));
            std::cout << "Segregated:: Reporting time:: Set time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;        
    }
    
    return false;
}

/** @brief Computes the post processed results */
void TRMOrchestra::PostProcess(){
    
    fMonolithicMultiphaseAnalysis->PostProcessStep();
    
}

void TRMOrchestra::ComputationalMeshUniformRefinement(TPZCompMesh  *cmesh, int ndiv){
    TPZVec<int64_t > subindex;
    for (int64_t iref = 0; iref < ndiv; iref++) {
        TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
        int64_t nel = elvec.NElements();
        for(int64_t el=0; el < nel; el++){
            TPZCompEl * cel = elvec[el];
            if(!cel) continue;
            int64_t ind = cel->Index();
            cel->Divide(ind, subindex, 0);
        }
    }
}


///** @brief Compute the system of equations using transfer matrixces */
//TPZFMatrix<STATE> TRMOrchestra::IntegrateResidue(TPZAutoPointer<TPZCompMesh> cmesh_multiphysics, TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer< TPZCompMesh> cmesh_pressure, TPZAutoPointer<TRMBuildTransfers> transfer){
//    
//    TPZFMatrix<STATE> rhs_flux(cmesh_flux->NEquations(),1,0.0);
//    TPZFMatrix<STATE> rhs_pressure(cmesh_pressure->NEquations(),1,0.0);
//    
//    // once we have
//    
//    TPZFMatrix<STATE> u_x_at_intpoints, u_y_at_intpoints, u_z_at_intpoints, divu_at_intpoints;
//    transfer->GetTransfer_X_Flux_To_Mixed_V().Multiply(cmesh_flux->Solution(), u_x_at_intpoints);
//    transfer->GetTransfer_Y_Flux_To_Mixed_V().Multiply(cmesh_flux->Solution(), u_y_at_intpoints);
//    transfer->GetTransfer_Z_Flux_To_Mixed_V().Multiply(cmesh_flux->Solution(), u_z_at_intpoints);
//    transfer->GetTransferDivergenceTo_Mixed_V().Multiply(cmesh_flux->Solution(),divu_at_intpoints);
//    
//    TPZFMatrix<STATE> p_at_intpoints;
//    transfer->GetTransferPressure_To_Mixed_V().Multiply(cmesh_pressure->Solution(), p_at_intpoints);
//
//    
//    // Integrate the volumetric forms
////    transfer->GetTransfer_X_Flux_To_Mixed().Print("Flux_x = ");
////    transfer->GetTransfer_Y_Flux_To_Mixed().Print("Flux_y = ");
////    transfer->GetTransfer_Z_Flux_To_Mixed().Print("Flux_z = ");
////    transfer->GetTransferPressure_To_Mixed().Print("Pressure = ");
////    transfer->GetTransferDivergenceTo_Mixed().Print("DivFlux = ");
////    transfer->GetJacobianDet_To_Mixed().Print("det = ");
////    transfer->GetWeightsTo_Mixed().Print("w = ");
////    transfer->GetRhs_To_Mixed().Print("Rhs = ");
////    divu_at_intpoints.Print("Div u = ");
//    
//
//    
//    int64_t nel = cmesh_multiphysics->NElements();
//    
//    // Getting the total integration point of the destination cmesh
//    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
//    TPZMatWithMem<TRMMemory,TPZMaterial> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZMaterial> *>(material);
//    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
//
//    
//    int origin      = 1;
////    int destination = 1;
//    int volumetric_elements = 0;
//    TPZCompMesh * cmesh_o;
////    TPZCompMesh * cmesh_d;
//    int mesh_o_nequ = 0;
//    
//    for (int64_t icel = 0; icel < nel; icel++) {
//        
//        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
//        if (!cel) {
//            DebugStop();
//        }
//        
//        //         TPZCondensedCompEl * mf_cel_condensed = dynamic_cast<TPZCondensedCompEl *> (cel);
//        //         if(!mf_cel_condensed){
//        //             DebugStop();
//        //         }
//        //        TPZCompEl * mf_cel_cond_ref = mf_cel_condensed->ReferenceCompEl();
//        
//        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
//        if(!mf_cel)
//        {
//            DebugStop();
//        }
//        
//        // Avoiding all materials that don't correspond to volumetric ones
//        if(mf_cel->Material()->Id() != _ReservMatId)
//        {
//            continue;
//        }
//        
//        TPZInterpolationSpace * intel_o = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(origin));
////        TPZInterpolationSpace * intel_d = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(destination));
//        if (!intel_o /*|| !intel_d*/) {
//            DebugStop();
//        }
//        
//        volumetric_elements++;
//        cmesh_o = intel_o->Mesh();
////        cmesh_d = intel_d->Mesh();
//        if (!cmesh_o /*|| !cmesh_d*/) {
//            DebugStop();
//        }
//        
//        mesh_o_nequ = cmesh_o->NEquations();
//        
//        // Computing the global integration points indexes
//        TPZManVector<int64_t> globindexes;
//        mf_cel->GetMemoryIndices(globindexes);
//        
//        int nshapes = intel_o->NShapeF();
//        int npoints = globindexes.size();
//        
//
//        
//        // Computing the global dof indexes for pressure mesh
//        TPZManVector<int64_t> pressure_dofs(nshapes);
//        
//        int nconnect = intel_o->NConnects();
//        int iphicount = 0;
//        // Compute all the phi values and push inside corresponding j-destination;
//        for (int icon = 0; icon < nconnect; icon++) {
//            TPZConnect  & con = intel_o->Connect(icon);
//            int64_t seqnumber = con.SequenceNumber();
//            int64_t position = cmesh_o->Block().Position(seqnumber);
//            int nconnectshape = con.NShape();
//            
//            for (int ish=0; ish < nconnectshape; ish++) {
//                int64_t i_equ = position + ish;
//                // Computing the global dof for pressure mesh
//                pressure_dofs[iphicount] = i_equ;
//                iphicount++;
//            }
//            
//        }
//
//
//        TPZFMatrix<double> point_value_phi(npoints,1);
//        TPZFMatrix<double> point_value_det(npoints,1);
//        TPZFMatrix<double> point_value_w(npoints,1);
//        TPZFMatrix<double> point_value_divu(npoints,1);
//        TPZFMatrix<double> point_value_rhs(npoints,1);
//
//        TPZManVector<int64_t> phi_index(1);
//        TPZManVector<int64_t> el_index(1);
//        TPZManVector<int64_t> unique_index(1,0);
//        
//        divu_at_intpoints.GetSub(globindexes, unique_index, point_value_divu);
//        
//        for (int iphi = 0; iphi < nshapes; iphi++) {
//            // Testing the volumetric forms of the weak statement
//            
//            int64_t equ = pressure_dofs[iphi];
//            phi_index[0] = equ;
//            el_index[0] = volumetric_elements-1;
//            
////            transfer->GetWeightsTo_Mixed().GetSub(globindexes, el_index, point_value_w);
////            transfer->GetJacobianDet_To_Mixed().GetSub(globindexes, el_index, point_value_det);
////            transfer->GetRhs_To_Mixed().GetSub(globindexes, el_index, point_value_rhs);
////            
////            transfer->GetTransferPressure_To_Mixed().GetSub(globindexes,phi_index,point_value_phi);
//            
//            
//            STATE i_integral = 0.0;
//            for (int ip = 0; ip < npoints; ip++) {
//                i_integral += point_value_w(ip,0)*point_value_det(ip,0)*(- point_value_divu(ip,0)+ point_value_rhs(ip,0)) * point_value_phi(ip,0);
//            }
//            rhs_pressure(equ,0) += i_integral;
//        }
//        
//    }
//    
//    rhs_pressure.Print("rhs_pressure = ");
//    
//    return rhs_pressure;
//
//}

/** @brief Compute gradient of the system of equations using transfer matrixces */
void TRMOrchestra::IntegrateGradientOfResidue(TPZAutoPointer< TPZCompMesh> cmesh_pressure, TPZAutoPointer<TRMBuildTransfers> transfer){
   
    DebugStop();
}

/** @brief Create a dual analysis using space odissey */
void TRMOrchestra::CreateAnalysisDual(){
    DebugStop();
//    int nel = 2;
//    TPZManVector<REAL,2> dx(2,nel), dy(2,nel), dz(2,nel);
//    dx[0] = 1;
//    dy[0] = 1;
//    dz[0] = 1;
//    
//    fSpaceGenerator->CreateGeometricBoxMesh(dx, dy, dz);
////    fSpaceGenerator.CreateGeometricReservoirMesh();
//#ifdef PZDEBUG
//    fSpaceGenerator->PrintGeometry();
//#endif
//    
//    fSpaceGenerator->SetDefaultPOrder(2);
//    
//    fSpaceGenerator->CreateFluxCmesh();
//    fSpaceGenerator->CreatePressureCmesh();
//
//    
//    fSpaceGenerator->CreateMixedCmesh();
//    
////    ProjectExactSolution();
//    
//    
//    fSpaceGenerator->IncreaseOrderAroundWell(2);
//
//    fSpaceGenerator->ConfigureWellConstantPressure(0., 1000.);
//
//    fSpaceGenerator->StaticallyCondenseEquations();
//    
//    // transfer the solution from the meshes to the multiphysics mesh
//    TPZManVector<TPZAutoPointer<TPZCompMesh>,3 > meshvec(2);
//    meshvec[0] = fSpaceGenerator->FluxCmesh();
//    meshvec[1] = fSpaceGenerator->PressureCmesh();
//    TPZAutoPointer<TPZCompMesh > Cmesh = fSpaceGenerator->MixedFluxPressureCmesh();
//    
//    
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, Cmesh);
//    
//    
//    // Analysis
//    bool mustOptimizeBandwidth = true;
//    fFluxPressureAnalysis->SetCompMesh(Cmesh.operator->(), mustOptimizeBandwidth);
////    TPZAnalysis * AnalysisDual = new TPZAnalysis(Cmesh.operator->(),mustOptimizeBandwidth);
//    int numofThreads = 8;
//#ifdef PZDEBUG
//    {
//        std::ofstream out("../MFCompMesh.txt");
//        Cmesh->Print(out);
//    }
//#endif
//    fFluxPressureAnalysis->Solution().Zero();
//    fFluxPressureAnalysis->LoadSolution();
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, Cmesh);
//    TPZFMatrix<STATE> prevsol = fFluxPressureAnalysis->Solution();
//    std::cout << "Total dof: " << prevsol.Rows() << std::endl;
//    
//    std::map<REAL,STATE> RatebyPosition;
//    STATE TotalFlux = 0.;
//     ComputeProductionRate(RatebyPosition, TotalFlux);
//
//    TPZSkylineStructMatrix strmat(Cmesh.operator->());
////    TPZSkylineNSymStructMatrix strmat(Cmesh.operator->());
//    TPZStepSolver<STATE> step;
//    strmat.SetNumThreads(numofThreads);
//    step.SetDirect(ELDLt);
//    fFluxPressureAnalysis->SetStructuralMatrix(strmat);
//    fFluxPressureAnalysis->SetSolver(step);
//    fFluxPressureAnalysis->Run();
//    std::cout << "Rhs norm " << Norm(fFluxPressureAnalysis->Rhs()) << std::endl;
//    prevsol -= fFluxPressureAnalysis->Solution();
//    Cmesh->LoadSolution(prevsol);
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, Cmesh);
//
//    RatebyPosition.clear();
//    TotalFlux = 0.;
//    ComputeProductionRate(RatebyPosition, TotalFlux);    
//    std::cout << "Total flux " << TotalFlux << std::endl;
//    for (std::map<REAL,STATE>::iterator it = RatebyPosition.begin(); it != RatebyPosition.end(); it++) {
//        std::cout << "Y_position " << it->first << " rate " << it->second << std::endl;
//    }
//    
//    if (0)
//    {
//        fFluxPressureAnalysis->Run();
//        std::cout << "Rhs norm " << Norm(fFluxPressureAnalysis->Rhs()) << std::endl;
//        prevsol -= fFluxPressureAnalysis->Solution();
//        Cmesh->LoadSolution(prevsol);
//        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, Cmesh);
//        fFluxPressureAnalysis->AssembleResidual();
//        std::cout << "Rhs norm " << Norm(fFluxPressureAnalysis->Rhs()) << std::endl;
//    }
//    const int dim = 3;
//    int div = 0;
//    TPZStack<std::string> scalnames, vecnames;
//    std::string plotfile =  "../DualDarcy.vtk";
//    scalnames.Push("WeightedPressure");
//    scalnames.Push("DivOfBulkVeclocity");
//    scalnames.Push("POrder");
//    vecnames.Push("BulkVelocity");
//    fFluxPressureAnalysis->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
//    fFluxPressureAnalysis->PostProcess(div);
    
}

/** @brief Create computational meshes using space odissey */
void TRMOrchestra::CreateCompMeshes(){
    DebugStop();
}


/** @brief Project an exact solution */
void TRMOrchestra::ProjectExactSolution()
{
    TPZAutoPointer<TPZCompMesh> mesh = fSpaceGenerator->FluxCmesh();
    if (!mesh) {
        DebugStop();
    }
    // copiar os materiais
    std::map<int, TPZMaterial *> matmap;// = mesh->MaterialVec();
    std::map<int, TPZMaterial *>::iterator it;
    for (it= mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        it->second->Clone(matmap);
    }
    
    // put L2 projection material
    TPZVecL2 *vecmat = new TPZVecL2(_ReservMatId);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(TRMOrchestra::ExactFlux);
    vecmat->SetForcingFunction(force);
    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        delete it->second;
    }
    mesh->MaterialVec().clear();
    mesh->InsertMaterialObject(vecmat);
    {
        TPZAnalysis an(mesh,false);
        TPZSkylineStructMatrix strmat(mesh);
        std::set<int> matids;
        matids.insert(_ReservMatId);
        strmat.SetMaterialIds(matids);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Run();
        TPZStack<std::string> scalnames,vecnames;
        vecnames.Push("state");
        an.DefineGraphMesh(3, scalnames, vecnames, "../ProjectFlux.vtk");
        an.PostProcess(0);
    }
    delete vecmat;
    mesh->MaterialVec() = matmap;
    
    mesh = fSpaceGenerator->PressureCmesh();
    matmap.clear();// = mesh->MaterialVec();
    for (it= mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        it->second->Clone(matmap);
    }
    int nstate = 1;
    TPZManVector<STATE,1> sol(1,1.);
    int dim = 3;
    TPZL2Projection *l2proj = new TPZL2Projection(_ReservMatId,dim,nstate,sol);
    TPZAutoPointer<TPZFunction<STATE> > force2 = new TPZDummyFunction<STATE>(TRMOrchestra::ExactPressure);
    l2proj->SetForcingFunction(force2);

    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        delete it->second;
    }
    mesh->MaterialVec().clear();
    mesh->InsertMaterialObject(l2proj);
    {
        TPZAnalysis an(mesh,false);
        TPZSkylineStructMatrix strmat(mesh);
        std::set<int> matids;
        matids.insert(_ReservMatId);
        strmat.SetMaterialIds(matids);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Run();
//        mesh->Print(std::cout);
        TPZStack<std::string> scalnames,vecnames;
        scalnames.Push("state");
        an.DefineGraphMesh(3, scalnames, vecnames, "../ProjectPressure.vtk");
        an.PostProcess(0);
    }
    delete l2proj;
    mesh->MaterialVec() = matmap;
    
    mesh = fSpaceGenerator->MixedFluxPressureCmesh();
    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        TPZBndCond * bnd = dynamic_cast<TPZBndCond *>(it->second);
        if (bnd)
        {
            bnd->SetForcingFunction(0,force2);
        }
    }

}

/** @brief exact pressure */
void TRMOrchestra::ExactPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &pressure)
{
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    pressure[0] = (1. - x)*x + (1. - y)*y + (1. - z)*z;
}

/** @brief exact flux */
void TRMOrchestra::ExactFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &flux)
{
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    flux[0] = -1. + 2*x;//2.*(-0.5 + x)*(-1. + y)*y*(-1. + z)*z;
    flux[1] = -1. + 2*y;//2.*(-1. + x)*x*(-0.5 + y)*(-1. + z)*z;
    flux[2] = -1. + 2*z;//2.*(-1. + x)*x*(-1. + y)*y*(-0.5 + z);
    
}

/** @brief exact laplacian */
void TRMOrchestra::ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &f)
{
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    f[0] = 6;//2.*(-1. + x)*x*(-1. + y)*y + 2.*(-1. + x)*x*(-1. + z)*z + 2.*(-1. + y)*y*(-1. + z)*z;
}

/** @brief Compute the production rate of the reservoir */
void TRMOrchestra::ComputeProductionRate(std::map<REAL,STATE> &RatebyPosition, STATE &TotalIntegral)
{
    fSpaceGenerator->Gmesh()->ResetReference();
    fSpaceGenerator->FluxCmesh()->LoadReferences();
    TPZAutoPointer<TPZGeoMesh> gmesh = fSpaceGenerator->Gmesh();
    TPZInt1d intrule(14);
    int np = intrule.NPoints();
    TotalIntegral = 0.;
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != _WellMatId3D) {
            continue;
        }
//        std::cout << "Accumulating for well3d = " << gel->Index() << std::endl;
        TPZStack<int64_t> faceElements;
        for (int side=21; side < 25 ; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour.Element()->MaterialId() != _Well3DReservoirFaces && neighbour != gelside) {
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                DebugStop();
            }
            faceElements.Push(neighbour.Element()->Index());
        }
        for (int i=0; i<np; i++) {
            REAL yIntegral = 0.;
            REAL lineIntegral = 0.;
            STATE flux1DIntegral = 0.;
            TPZManVector<REAL,2> posi(1);
            TPZStack<REAL> yvals;
            REAL wi;
            intrule.Point(i,posi,wi);
            for (int iface=0; iface< faceElements.size(); iface++) {
                for (int j=0; j<np; j++) {
                    TPZGeoEl *gelface = gmesh->Element(faceElements[iface]);
                    TPZCompEl *celface = gelface->Reference();
                    TPZFNMatrix<9,REAL> jac(2,2),jacinv(2,2),axes(2,3),gradx(3,2);
                    REAL detjac,wj;
                    TPZManVector<REAL,3> posj(1),intpoint(2),xco(3);
                    intrule.Point(j, posj, wj);
                    intpoint[1] = posi[0];
                    intpoint[0] = posj[0];
                    gelface->Jacobian(intpoint, jac, axes, detjac, jacinv);
                    gelface->X(intpoint, xco);
//                    std::cout << "el index " << gelface->Index() << " intpoint " << intpoint << " xco " << xco << std::endl;
                    TPZAxesTools<REAL>::ComputeGradX(jac, axes, gradx);
                    int xdir = 0;
                    REAL detjac1d = sqrt(gradx(0,xdir)*gradx(0,xdir)+gradx(1,xdir)*gradx(1,xdir)+gradx(2,xdir)*gradx(2,xdir));
                    lineIntegral += detjac1d*wj;
                    yIntegral += detjac1d*wj*xco[1];
                    yvals.Push(xco[1]);
                    TPZManVector<STATE,1> sol(1);
                    celface->Solution(intpoint, 0, sol);
                    flux1DIntegral += detjac1d*wj*sol[0];
                    TotalIntegral += fabs(detjac)*wi*wj*sol[0];
                }
            }
//            std::cout << "Y values " << yvals << std::endl;
            REAL averageY = yIntegral/lineIntegral;
            RatebyPosition[averageY] = flux1DIntegral;
        }
    }
}

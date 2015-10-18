/*
 *  pznondarcyanalysis.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZDarcyAnalysis.h"
#include "pzcheckgeom.h"
#include "pzlog.h"

#include "TPZReadGIDGrid.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "tpzhierarquicalgrid.h"
#include "pzelementgroup.h"

#include "TPZVTKGeoMesh.h"
#include "TPZAxiSymmetricDarcyFlow.h"
#include "TPZAxiSymmetricDarcyFlowH1.h"
#include "pzpoisson3d.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzfstrmatrix.h"
#include "math.h"

#include "TPZMultiphysicsInterfaceEl.h"


#ifdef DEBUG
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.DarcyFlow"));
#endif
#endif

//#define SolutionI
//#define SolutionII

TPZDarcyAnalysis::TPZDarcyAnalysis(TPZAutoPointer<SimulationData> DataSimulation, TPZVec<TPZAutoPointer<ReservoirData> > Layers, TPZVec<TPZAutoPointer<PetroPhysicData> > PetroPhysic)
{
    
    fmeshvecini.Resize(2);
    fmeshvec.Resize(2);
    
    fgmesh=NULL;
    fcmeshinitialdarcy=NULL;
    fcmeshdarcy=NULL;
    fcmesh = NULL;
    
    falpha_fluid=NULL;
    fbeta_fluid=NULL;
    fgamma_fluid=NULL;
    
    // Vector which will store tha residuum in the last state (n)
    fResidualAtn.Resize(0, 0);
    
    // Vector which will store tha residuum in the last state (n+1)
    fResidualAtnplusOne.Resize(0, 0);
    
    fSimulationData     = DataSimulation;
    fLayers             = Layers;
    fRockPetroPhysic    = PetroPhysic;
    
    /** @brief unknowns for n time step */
    falphaAtn.Resize(0, 0);
    
    /** @brief unknowns for n+1 time step */
    falphaAtnplusOne.Resize(0, 0);
    
    /** @brief Store DOF associated with Constant Saturations */
    fConstantSaturations.Resize(0);
    
    /** @brief Store DOF associated with  Saturation gradients */
    fGradientSaturations.Resize(0);
    
}


TPZDarcyAnalysis::~TPZDarcyAnalysis()
{
    
}


void TPZDarcyAnalysis::SetFluidData(TPZVec< TPZAutoPointer<ReducedPVT> > PVTData){
    
    falpha_fluid = PVTData[0];
    fbeta_fluid = PVTData[1];
    fgamma_fluid = PVTData[2];
    
    
    
    //    if(fSimulationData->IsOnePhaseQ()){
    //
    //        falpha_fluid = PVTData[0];
    //
    //    }
    
    if(fSimulationData->IsTwoPhaseQ()){
        
        fmeshvecini.Resize(3);
        fmeshvec.Resize(3);
        
    }
    
    if(fSimulationData->IsThreePhaseQ()){
        
        fmeshvecini.Resize(4);
        fmeshvec.Resize(4);
        
    }
    
}

void TPZDarcyAnalysis::SetLastState()
{
    fSimulationData->SetnStep(true);
}

void TPZDarcyAnalysis::SetNextState()
{
    fSimulationData->SetnStep(false);
}

void TPZDarcyAnalysis::Assemble()
{
    
}

/**
 * Push the initial cmesh
 */
void TPZDarcyAnalysis::PushInitialCmesh(){
    fcmesh = NULL;
    fcmesh = fcmeshinitialdarcy;
    
}

/**
 * Push the current cmesh
 */
void TPZDarcyAnalysis::PushCmesh(){
    fcmesh = NULL;
    fcmesh = fcmeshdarcy;
    
}

void TPZDarcyAnalysis::AssembleLastStep(TPZAnalysis *an)
{
    fcmesh->LoadSolution(falphaAtn);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmesh);
    SetLastState();
    an->AssembleResidual();
    fResidualAtn = an->Rhs();
    
    //#ifdef DEBUG
    //   #ifdef LOG4CXX
    //       if(logger->isDebugEnabled())
    //       {
    //           std::stringstream sout;
    //           falphaAtn.Print("falphaAtn = ", sout,EMathematicaInput);
    //           fResidualAtn.Print("fResidualAtn = ", sout,EMathematicaInput);
    //           LOGPZ_DEBUG(logger,sout.str())
    //       }
    //   #endif
    //#endif
    
}

void TPZDarcyAnalysis::AssembleResNextStep(TPZAnalysis *an)
{
    fcmesh->LoadSolution(falphaAtnplusOne);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmesh);
    SetNextState();
    an->AssembleResidual();
    fResidualAtnplusOne = an->Rhs();
    
    //#ifdef DEBUG
    //  #ifdef LOG4CXX
    //      if(logger->isDebugEnabled())
    //      {
    //          std::stringstream sout;
    //          falphaAtnplusOne.Print("falphaAtnplusOne = ", sout,EMathematicaInput);
    //          fResidualAtnplusOne.Print("fResidualAtnplusOne = ", sout,EMathematicaInput);
    //          LOGPZ_DEBUG(logger,sout.str())
    //      }
    //  #endif
    //#endif
    
}

void TPZDarcyAnalysis::AssembleNextStep(TPZAnalysis *an)
{
    fcmesh->LoadSolution(falphaAtnplusOne);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmesh);
    SetNextState();
    an->Assemble();
    fResidualAtnplusOne = an->Rhs();
    
    // #ifdef DEBUG
    //     #ifdef LOG4CXX
    //         if(logger->isDebugEnabled())
    //         {
    //             std::stringstream sout;
    //             falphaAtnplusOne.Print("falphaAtnplusOne = ", sout,EMathematicaInput);
    //             fResidualAtnplusOne.Print("fResidualAtnplusOne = ", sout,EMathematicaInput);
    //             LOGPZ_DEBUG(logger,sout.str())
    //         }
    //     #endif
    // #endif
    
}

void TPZDarcyAnalysis::UpDateAlphaVec(TPZFMatrix<REAL> &alpha)
{
    falphaAtn = alpha;
    falphaAtnplusOne = alpha;
    
}

void TPZDarcyAnalysis::Residual(TPZFMatrix<STATE> &residual, int icase)
{
    //    TPZNonLinearAnalysis::Residual(residual, icase);
    //    residual = fResidualLastState + residual;
}

void TPZDarcyAnalysis::ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase)
{
    this->SetNextState();
    TPZDarcyAnalysis::ComputeTangent(tangent, coefs, icase);
}

void TPZDarcyAnalysis::InitializeSolution(TPZAnalysis *an)
{
    
    if (fSimulationData->IsTwoPhaseQ()) {
        TPZVec<STATE> Swlini(fmeshvec[2]->Solution().Rows());
        TPZCompMesh * L2Sw = L2ProjectionCmesh(Swlini);
        TPZAnalysis * L2Analysis = new TPZAnalysis(L2Sw,false);
        SolveProjection(L2Analysis,L2Sw);
        fmeshvec[2]->LoadSolution(L2Analysis->Solution());
    }
    
    
    if (fSimulationData->IsThreePhaseQ()) {
        TPZFMatrix<REAL> Soini = fmeshvec[2]->Solution();
        int SoDOF = Soini.Rows();
        for (int iso = 0; iso < SoDOF; iso++)
        {
            Soini(iso,0) = 1.0 - Soini(iso,0);
        }
        fmeshvec[3]->LoadSolution(Soini);
    }
    
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshinitialdarcy);
    
    falphaAtn = fcmeshinitialdarcy->Solution();
    falphaAtnplusOne = fcmeshinitialdarcy->Solution();
    
    REAL TimeStep = fSimulationData->GetDeltaT();
    REAL BigTimeStep = 100000.0;
    fSimulationData->SetDeltaT(BigTimeStep);
    this->AssembleLastStep(an);
    this->AssembleNextStep(an);
    const clock_t tinia = clock();
    NewtonIterations(an);
    const clock_t tenda = clock();
    const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
    std::cout << "Time for Newton: " << timea << std::endl;
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshinitialdarcy);
    
    
    
    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshdarcy);
    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshdarcy);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshdarcy);
    
    fcmeshdarcy->LoadSolution(fcmeshinitialdarcy->Solution());
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    
    fSimulationData->SetDeltaT(TimeStep);
    falphaAtn = fcmeshdarcy->Solution();
    falphaAtnplusOne = fcmeshdarcy->Solution();
    
}

void TPZDarcyAnalysis::Run()
{
    
    std::string dirname = PZSOURCEDIR;
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
#ifdef DEBUG
    #ifdef LOG4CXX
        
        std::string FileName = dirname;
        FileName = dirname + "/Projects/DarcyflowAxisymmetricHdiv/";
        FileName += "DarcyFlowLog.cfg";
        InitializePZLOG(FileName);
        
    #endif
#endif
    
    //  Reading mesh
    std::string GridFileName;
    GridFileName = dirname + "/Projects/DarcyflowAxisymmetricHdiv/";
    GridFileName += "SingleLayer.dump";
    //GridFileName += "MixLayer.dump";
    //GridFileName += "BatatacoarseQ.dump";
    //GridFileName += "QUAD4.dump";
    
    if(fLayers[0]->GetIsGIDGeometry())
    {
        ReadGeoMesh(GridFileName);
    }
    else
    {
        
        int nx = fSimulationData->GetnElementsx();
        int ny = fSimulationData->GetnElementsy();
        GeometryLine(nx,ny);
        //        CreatedGeoMesh();
    }
    
    
    REAL deg = 0.0;
    int hcont = 0;
    RotateGeomesh(deg * M_PI/180.0);
    this->UniformRefinement(fSimulationData->GetHrefinement());
    
    std::set<int> matidstoRef;
    //    matidstoRef.insert(2);
    //    matidstoRef.insert(3);
    //    matidstoRef.insert(4);
    matidstoRef.insert(3);
    matidstoRef.insert(5);
    
    this->UniformRefinement(hcont, matidstoRef);
    this->PrintGeoMesh();
    
    int q = fSimulationData->Getqorder();
    int p = fSimulationData->Getporder();
    int s = fSimulationData->Getsorder();
    
    CreateMultiphysicsMesh(q,p,s);
    CreateInterfaces();
    
    this->PushInitialCmesh();
    this->PrintCmesh();
    
    if(fSimulationData->GetSC())
    {
        this->ApplyStaticCondensation();
    }
    
    // Analysis
    
    TPZAnalysis * Initialan = CreateAnalysis(fcmeshinitialdarcy);
    this->PushInitialCmesh();
    this->PrintCmesh();
    this->InitializeSolution(Initialan);
    delete Initialan;
    
    TPZAnalysis * an = CreateAnalysis(fcmeshdarcy);
    this->PushCmesh();
    this->PrintCmesh();
    //    this->CheckGlobalJacobian(an);
    //    this->CheckGlobalConvergence(an);
    this->TimeForward(an);
    
    return;
}

TPZAnalysis * TPZDarcyAnalysis::CreateAnalysis(TPZCompMesh * cmesh){
    
    if (!cmesh) {
        std::cout << "NUll computational mesh " << std::endl;
        DebugStop();
    }
    
    bool mustOptimizeBandwidth = fSimulationData->GetOptband();
    TPZAnalysis *an = new TPZAnalysis(cmesh,mustOptimizeBandwidth);
    int numofThreads = 0;
    
    bool IsDirecSolver = fSimulationData->GetIsDirect();
    
    if (IsDirecSolver) {
        
        if (fSimulationData->GetIsBroyden()) {
            TPZFStructMatrix fullMatrix(cmesh);
            TPZStepSolver<STATE> step;
            fullMatrix.SetNumThreads(numofThreads);
            step.SetDirect(ELU);
            an->SetStructuralMatrix(fullMatrix);
            an->SetSolver(step);
        }
        else{
            
            TPZSkylineNSymStructMatrix skylnsym(cmesh);
            TPZStepSolver<STATE> step;
            skylnsym.SetNumThreads(numofThreads);
            step.SetDirect(ELU);
            an->SetStructuralMatrix(skylnsym);
            an->SetSolver(step);
        }
        
    }
    else
    {
        if (fSimulationData->GetIsBroyden()) {
            TPZFStructMatrix fullMatrix(cmesh);
            fullMatrix.SetNumThreads(numofThreads);
            
            TPZAutoPointer<TPZMatrix<STATE> > fullMatrixa = fullMatrix.Create();
            TPZAutoPointer<TPZMatrix<STATE> > fullMatrixaClone = fullMatrixa->Clone();
            
            TPZStepSolver<STATE> *stepre = new TPZStepSolver<STATE>(fullMatrixaClone);
            TPZStepSolver<STATE> *stepGMRES = new TPZStepSolver<STATE>(fullMatrixa);
            TPZStepSolver<STATE> *stepGC = new TPZStepSolver<STATE>(fullMatrixa);
            stepre->SetDirect(ELU);
            stepre->SetReferenceMatrix(fullMatrixa);
            stepGMRES->SetGMRES(10, 20, *stepre, 1.0e-10, 0);
            stepGC->SetCG(10, *stepre, 1.0e-10, 0);
            if (fSimulationData->GetIsCG()) {
                an->SetSolver(*stepGC);
            }
            else{
                an->SetSolver(*stepGMRES);
            }
            
        }
        else{
            
            TPZSkylineNSymStructMatrix skylnsym(cmesh);
            skylnsym.SetNumThreads(numofThreads);
            
            TPZAutoPointer<TPZMatrix<STATE> > skylnsyma = skylnsym.Create();
            TPZAutoPointer<TPZMatrix<STATE> > skylnsymaClone = skylnsyma->Clone();
            
            TPZStepSolver<STATE> *stepre = new TPZStepSolver<STATE>(skylnsymaClone);
            TPZStepSolver<STATE> *stepGMRES = new TPZStepSolver<STATE>(skylnsyma);
            TPZStepSolver<STATE> *stepGC = new TPZStepSolver<STATE>(skylnsyma);
            
            stepre->SetDirect(ELU);
            stepre->SetReferenceMatrix(skylnsyma);
            stepGMRES->SetGMRES(10, 20, *stepre, 1.0e-10, 0);
            stepGC->SetCG(10, *stepre, 1.0e-10, 0);
            if (fSimulationData->GetIsCG()) {
                an->SetSolver(*stepGC);
            }
            else{
                an->SetSolver(*stepGMRES);
            }
        }
        
    }
    
    return an;
    
}

void TPZDarcyAnalysis::CreateInterfaces()
{
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    fgmesh->ResetReference();
    fcmeshinitialdarcy->LoadReferences();
    
    // Creation of interface elements
    int nel = fcmeshdarcy->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = fcmeshinitialdarcy->ElementVec()[el];
        if(!compEl) continue;
        TPZGeoEl * gel = compEl->Reference();
        if(!gel) {continue;}
        if(gel->HasSubElement()) {continue;}
        int index = compEl ->Index();
        if(compEl->Dimension() == fcmeshinitialdarcy->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fcmeshinitialdarcy->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }
    }
    
    if (fSimulationData->IsTwoPhaseQ()) {
        return;
    }
    
    fgmesh->ResetReference();
    fcmeshdarcy->LoadReferences();
    
    // Creation of interface elements
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = fcmeshdarcy->ElementVec()[el];
        if(!compEl) continue;
        TPZGeoEl * gel = compEl->Reference();
        if(!gel) {continue;}
        if(gel->HasSubElement()) {continue;}
        int index = compEl ->Index();
        if(compEl->Dimension() == fcmeshdarcy->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fcmeshdarcy->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }
    }
    
    return;
}


void TPZDarcyAnalysis::PrintLS(TPZAnalysis *an)
{
    TPZAutoPointer< TPZMatrix<REAL> > KGlobal;
    TPZFMatrix<STATE> FGlobal;
    KGlobal =   an->Solver().Matrix();
    FGlobal =   an->Rhs();
#ifdef DEBUG
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        KGlobal->Print("KGlobal = ", sout,EMathematicaInput);
        FGlobal.Print("FGlobal = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
#endif
    
}

void TPZDarcyAnalysis::PrintCmesh(){
    
    
#ifdef DEBUG
#ifdef LOG4CXX
    
    std::ofstream out("ComputationalMesh.txt");
    fcmesh->Print(out);
    
#endif
#endif
    
}

void TPZDarcyAnalysis::CreateMultiphysicsMesh(int q, int p, int s)
{
    fmeshvec[0] = CmeshFlux(q);
    fmeshvec[1] = CmeshPressure(p);
    
    if (fSimulationData->IsTwoPhaseQ()) {
        fmeshvec[2] = CmeshSw(s);
    }
    
    if (fSimulationData->IsThreePhaseQ()) {
        fmeshvec[2] = CmeshSw(s);
        fmeshvec[3] = CmeshSo(s);
    }
    
    fcmeshinitialdarcy = CmeshMixedInitial();
    fcmeshdarcy = CmeshMixed();
    
    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshinitialdarcy);
    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshinitialdarcy);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshinitialdarcy);
    
    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshdarcy);
    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshdarcy);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshdarcy);
    
    int DOF = fcmeshdarcy->Solution().Rows();
    std::cout << "Degrees of freedom: " << DOF << std::endl;
    
}

void TPZDarcyAnalysis::TimeForward(TPZAnalysis *an)
{
    int nsteps = fSimulationData->GetMaxTime() / fSimulationData->GetDeltaT();
    
    REAL tk = 0;
    
    
    //  Read the binary file
    //    {
    //      TPZBFileStream load;
    //      load.OpenRead("Solution.bin");
    //      TPZFMatrix<REAL> SolutiontoLoad = fcmeshdarcy->Solution();
    //      SolutiontoLoad.Read(load,0);
    //      fcmeshdarcy->LoadSolution(SolutiontoLoad);
    //      TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    //    }
    
    //    TPZFMatrix<REAL> TrueSolution = fcmeshdarcy->Solution();
    //  TrueSolution *=0.0;
    ////    TrueSolution(0,0) = -0.5 ;
    ////    TrueSolution(1,0) = -0.5 ;
    ////    TrueSolution(2,0) = 0.5;
    ////    TrueSolution(3,0) = 0.5;
    //
    ////    TrueSolution(12,0) = 0.1;//10.29367992;
    ////    TrueSolution(13,0) = 0.1;
    ////
    ////    TrueSolution(14,0) = 0.1;
    ////    TrueSolution(15,0) = 0.1;//10.29367992;
    ////
    ////    TrueSolution(16,0) = 1.0;
    ////    TrueSolution(17,0) = 0.0;
    //
    //    falphaAtn= TrueSolution;
    //    falphaAtnplusOne= TrueSolution;
    //
    //
    
    if (fSimulationData->GetGR())
    {
        this->FilterSaturationGradients(fConstantSaturations,fGradientSaturations);
        an->StructMatrix()->EquationFilter().Reset();
        an->StructMatrix()->EquationFilter().SetActiveEquations(fConstantSaturations);
        
        this->SaturationReconstruction(an);
    }
    
    this->fSimulationData->SetTime(tk);
    this->PostProcessVTK(an);
    
    for (int istep = 1 ; istep <=  nsteps; istep++) {
        tk = istep*this->fSimulationData->GetDeltaT();
        this->fSimulationData->SetTime(tk);
        
        std::cout << "Time value (days): " << tk << std::endl;
        std::cout<<  "Time step: " << istep << std::endl;
        
        
        this->AssembleLastStep(an);
        this->AssembleNextStep(an);
        
        if (fSimulationData->GetIsBroyden())
        {
            
            const clock_t tinia = clock();
            BroydenIterations(an);
            const clock_t tenda = clock();
            const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
            std::cout << "Time for Broyden: " << timea << std::endl;
            std::cout << "Number of DOF = " << fcmeshdarcy->Solution().Rows() << std::endl;
            this->PostProcessVTK(an);
            
        }
        else
        {
            const clock_t tinia = clock();
            NewtonIterations(an);
            const clock_t tenda = clock();
            const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
            std::cout << "Time for Newton: " << timea << std::endl;
            std::cout << "Number of DOF = " << fcmeshdarcy->Solution().Rows() << std::endl;
            this->PostProcessVTK(an);
            
            //          TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
            //          TPZFMatrix<REAL> SolutiontoSave = fcmeshdarcy->Solution();
            //          TPZBFileStream save;
            //          save.OpenWrite("Solution.bin");
            //          SolutiontoSave.Write(save,0);
            
        }
        
    }
    
    //    this->PostProcessVTK(an);
    
    
}


void TPZDarcyAnalysis::CleanUpGradients(TPZAnalysis *an){
    
    long numofdof = fGradientSaturations.size();
    TPZFMatrix<REAL> SolToLoad = an->Solution();
    for(long i=0; i < numofdof; i++)
    {
        SolToLoad(fGradientSaturations[i],0) = 0.0;
    }
    an->LoadSolution(SolToLoad);
    //   TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec,fcmeshdarcy);
    
}

void TPZDarcyAnalysis::SaturationReconstruction(TPZAnalysis *an)
{
    TPZGradientReconstruction *gradreconst = new TPZGradientReconstruction(false,1.);
    
    fmeshvec[2]->Reference()->ResetReference();
    fmeshvec[2]->LoadReferences();
    gradreconst->ProjectionL2GradientReconstructed(fmeshvec[2], fSimulationData->fMatL2);
    
    fmeshvec[3]->Reference()->ResetReference();
    fmeshvec[3]->LoadReferences();
    gradreconst->ProjectionL2GradientReconstructed(fmeshvec[3], fSimulationData->fMatL2);
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshdarcy);
    an->LoadSolution(fcmeshdarcy->Solution());
}

void TPZDarcyAnalysis::NewtonIterations(TPZAnalysis *an)
{
    
    TPZFMatrix<STATE> Residual(an->Rhs().Rows(),1,0.0);
    Residual = fResidualAtn + fResidualAtnplusOne;
    
    TPZFMatrix<STATE> X = falphaAtn;
    TPZFMatrix<STATE> DeltaX = falphaAtn;
    
    STATE error     =   1.0;
    STATE normdx    =   1.0;
    int iterations  =   0;
    int centinel    =   0;
    int fixed       =   fSimulationData->GetFixediterations();
    
    
    while (error >= fSimulationData->GetToleranceRes() && iterations <= fSimulationData->GetMaxiterations()) {
        
        an->Rhs() = Residual;
        an->Rhs() *= -1.0;
        
        this->PrintLS(an);
        
        an->Solve();
        
        
        DeltaX = an->Solution();
        normdx = Norm(DeltaX);
        X += DeltaX;
        
        falphaAtnplusOne=X;
        
        if (fSimulationData->GetGR())
        {
            this->SaturationReconstruction(an);
            CleanUpGradients(an);
        }
        
        if (((fixed+1) * (centinel) == iterations)) {
            
            this->AssembleNextStep(an);
            centinel++;
        }
        else{
            this->AssembleResNextStep(an);
        }
        
        fResidualAtnplusOne = an->Rhs();
        
        Residual = fResidualAtn + fResidualAtnplusOne;
        error = Norm(Residual);
        iterations++;
        
        
#ifdef DEBUG
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
            X.Print("X = ", sout,EMathematicaInput);
            Residual.Print("Residual = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
#endif
        
        if(error < fSimulationData->GetToleranceRes() || normdx < fSimulationData->GetToleranceDX())
        {
            std::cout << "Converged with iterations:  " << iterations << std::endl;
            std::cout << "error norm: " << error << std::endl;
            std::cout << "error of dx: " << normdx << std::endl;
            this->UpDateAlphaVec(X);
            break;
        }
        
        
        if (iterations == fSimulationData->GetMaxiterations()) {
            std::cout << "Out max iterations " << iterations << std::endl;
            std::cout << "error norm " << error << std::endl;
            this->UpDateAlphaVec(X);
            break;
        }
        
    }
    
    
}

void TPZDarcyAnalysis::BroydenIterations(TPZAnalysis *an)
{
    
    
    int m = an->Solution().Rows();
    
    TPZFMatrix<STATE> Residual(m,1,0.0);
    TPZFMatrix<STATE> Rank(m,m,0.0);
    TPZFMatrix<STATE> X(m,1,0.0);
    TPZFMatrix<STATE> DeltaX(m,1,0.0);
    
    TPZAutoPointer<TPZMatrix<STATE> >  D;
    TPZAutoPointer<TPZMatrix<STATE> >  Dk;
    
    TPZFMatrix<STATE> * DInverse =  new TPZFMatrix<STATE> (m,m,0.0);
    
    bool IsShermanMorrison = false;
    STATE ck = 0.0;
    //    TPZFMatrix<STATE> DInverse(m,m,0.0);
    TPZFMatrix<STATE> Identity(m,m,0.0);
    TPZFMatrix<STATE> DInverseT(m,m,0.0);
    TPZFMatrix<STATE> u(m,m,0.0);
    TPZFMatrix<STATE> du(m,m,0.0);
    TPZFMatrix<STATE> duT(m,m,0.0);
    TPZFMatrix<STATE> Ck(1,1,0.0);
    
    STATE error=1;
    STATE normdx=1;
    int iterations=0;
    int centinel    =0;
    int fixed       =fSimulationData->GetFixediterations();
    
    // Computing the first Newton iteration
    
    
    an->Rhs() = fResidualAtn + fResidualAtnplusOne;       // g(X0)
    
    
    an->Rhs() *= -1.0;
    
    if (IsShermanMorrison) {
        DInverse = ComputeInverse();
        DInverse->Multiply(Residual, DeltaX);
        D.operator->()->Multiply(*DInverse, Identity);
        
#ifdef DEBUG
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
            Residual.Print("Residual = ", sout,EMathematicaInput);
            Identity.Print("Identity = ", sout,EMathematicaInput);
            DInverse->Print("DInverse = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
#endif
        
    }
    else{
        D = an->Solver().Matrix();  // J(X0)
        an->Solve();
        DeltaX = an->Solution();    // d(X0)
    }
    normdx = Norm(DeltaX);      // d(X0)*d(X0)
    X += DeltaX;                // X1
    
    // End of newton iteration
    
    // Procedure without Inverse computation
    
    iterations++;
    
    while (error >= fSimulationData->GetToleranceRes() && iterations <= fSimulationData->GetMaxiterations()) {
        
        
        fcmeshdarcy->LoadSolution(X);
        //        if (!fSimulationData->GetIsH1approx())
        //        {
        //            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
        //        }
        
        this->AssembleResidual();
        Residual = fResidualAtn + fResidualAtnplusOne;       // g(Xk)
        error = Norm(Residual);     // g(Xk)*g(Xk)
        
        if(error <= fSimulationData->GetToleranceRes() || normdx <= fSimulationData->GetToleranceDX())
        {
            std::cout << "Converged with iterations:  " << iterations << std::endl;
            std::cout << "error norm: " << error << std::endl;
            std::cout << "error of dx: " << normdx << std::endl;
            break;
        }
#ifdef DEBUG
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
            Residual.Print("Residual = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
#endif
        
        if (((fixed+1) * (centinel) == iterations)) {
            
            if (IsShermanMorrison) {
                DInverse->Multiply(Residual, u);
                u.Add(DeltaX, du);
                u.Transpose(&duT);
                DeltaX.Multiply(duT,Ck);
                ck = Ck(0,0);
                Rank = TensorProduct(u, DeltaX);
                Rank *= -(1/ck);
                Rank.Multiply(*DInverse, DInverseT);
                DInverse->Add(DInverseT, *DInverse);
                
            }else
            {
                // Application of the Secant condition
                Rank = TensorProduct(Residual, DeltaX);
                Rank *= 1.0/(normdx*normdx);
                D.operator->()->Add(Rank, *D.operator->());
                an->Solver().Matrix() = D;
            }
            
            
            
            centinel++;
            
        }
        
        
        
        //#ifdef LOG4CXX
        //        if(logger->isDebugEnabled())
        //        {
        //            std::stringstream sout;
        //            an->Solver().Matrix().operator->()->Print("*an->Solver().Matrix().operator->() = ", sout,EMathematicaInput);
        //            D->Print("D = ", sout,EMathematicaInput);
        //            Rank.Print("Rank = ", sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger,sout.str())
        //        }
        //#endif
        an->Rhs() = fResidualAtn + fResidualAtnplusOne;
        an->Rhs() *= -1.0;
        
        if (IsShermanMorrison) {
            DInverse->Multiply(Residual, DeltaX);
        }
        else{
            an->Solve();
            DeltaX = an->Solution();    // d(Xk)
        }
        
        
        normdx = Norm(DeltaX);      // d(Xk)*d(Xk)
        X += DeltaX;                // Xk+1
        
        
        //#ifdef LOG4CXX
        //        if(logger->isDebugEnabled())
        //        {
        //            std::stringstream sout;
        //            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
        //            X.Print("X = ", sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger,sout.str())
        //        }
        //#endif
        iterations++;
        
        if (iterations == fSimulationData->GetMaxiterations()) {
            std::cout << "Out max iterations " << iterations << std::endl;
            std::cout << "error norm " << error << std::endl;
            break;
        }
        
    }
    
    
}


TPZFMatrix<STATE>  TPZDarcyAnalysis::TensorProduct(TPZFMatrix<STATE> &g, TPZFMatrix<STATE> &d)
{
    TPZFMatrix<STATE> dT=d;
    d.Transpose(&dT);
    TPZFMatrix<STATE> RankOne;
    g.Multiply(dT, RankOne);
    
    //#ifdef LOG4CXX
    //    if(logger->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        g.Print("g = ", sout,EMathematicaInput);
    //        dT.Print("dT = ", sout,EMathematicaInput);
    //        d.Print("d = ", sout,EMathematicaInput);
    //        RankOne.Print("RankOne = ", sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logger,sout.str())
    //    }
    //#endif
    
    return RankOne;
    
}

TPZCompMesh * TPZDarcyAnalysis::CmeshMixedInitial()
{
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    int dim = 2;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    
    
    // Material medio poroso
    TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
    mat->SetSimulationData(fSimulationData);
    mat->SetReservoirData(fLayers[ilayer]);
    mat->SetPetroPhysicsData(fRockPetroPhysic[ilayer]);
    mat->SetFluidAlpha(falpha_fluid);
    mat->SetFluidBeta(fbeta_fluid);
    mat->SetFluidGamma(fgamma_fluid);
    int nvars = 2 + fSimulationData->GetsystemType().size() - 1;
    mat->SetNvars(nvars);
    cmesh->InsertMaterialObject(mat);
    
    
    // Rigth hand side function
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Ffunction);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
//     mat->SetTimeDependentForcingFunction(forcef);
    
    // Setting up linear tracer solution
    //  TPZDummyFunction<STATE> *Ltracer = new TPZDummyFunction<STATE>(LinearTracer);
    TPZDummyFunction<STATE> *Ltracer = new TPZDummyFunction<STATE>(BluckleyAndLeverett);
    TPZAutoPointer<TPZFunction<STATE> > fLTracer = Ltracer;
    mat->SetTimeDependentFunctionExact(fLTracer);
    
    TPZManVector<REAL,4> Bottom     = fSimulationData->GetBottomBCini();
    TPZManVector<REAL,4> Right      = fSimulationData->GetRightBCini();
    TPZManVector<REAL,4> Top        = fSimulationData->GetTopBCini();
    TPZManVector<REAL,4> Left       = fSimulationData->GetLeftBCini();
    
    
    // Bc Bottom
    val2(0,0) = Bottom[1];
    val2(1,0) = Bottom[2];
    val2(2,0) = Bottom[3];
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, int(Bottom[0]), val1, val2);
    
    // Bc Right
    val2(0,0) = Right[1];
    val2(1,0) = Right[2];
    val2(2,0) = Right[3];
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, int(Right[0]), val1, val2);
    
    // Bc Top
    val2(0,0) = Top[1];
    val2(1,0) = Top[2];
    val2(2,0) = Top[3];
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, int(Top[0]), val1, val2);
    
    // Bc Left
    val2(0,0) = Left[1];
    val2(1,0) = Left[2];
    val2(2,0) = Left[3];
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, int(Left[0]), val1, val2);
    
    cmesh->InsertMaterialObject(bcBottom);
    cmesh->InsertMaterialObject(bcRight);
    cmesh->InsertMaterialObject(bcTop);
    cmesh->InsertMaterialObject(bcLeft);
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CmeshMixed()
{
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    int dim = 2;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    
    
    // Material medio poroso
    TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
    mat->SetSimulationData(fSimulationData);
    mat->SetReservoirData(fLayers[ilayer]);
    mat->SetPetroPhysicsData(fRockPetroPhysic[ilayer]);
    mat->SetFluidAlpha(falpha_fluid);
    mat->SetFluidBeta(fbeta_fluid);
    mat->SetFluidGamma(fgamma_fluid);
    int nvars = 2 + fSimulationData->GetsystemType().size() - 1;
    mat->SetNvars(nvars);
    cmesh->InsertMaterialObject(mat);
    
    
    // Rigth hand side function
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Ffunction);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
    mat->SetTimeDependentForcingFunction(forcef);
    
    // Setting up linear tracer solution
    //  TPZDummyFunction<STATE> *Ltracer = new TPZDummyFunction<STATE>(LinearTracer);
    TPZDummyFunction<STATE> *Ltracer = new TPZDummyFunction<STATE>(BluckleyAndLeverett);
    TPZAutoPointer<TPZFunction<STATE> > fLTracer = Ltracer;
    mat->SetTimeDependentFunctionExact(fLTracer);

/*    TPZDummyFunction<STATE> *bcnfunction = new TPZDummyFunction<STATE>(BCNfunction);
    TPZAutoPointer<TPZFunction<STATE> > bcN;
    bcnfunction->SetPolynomialOrder(20);
    bcN = bcnfunction;

    TPZDummyFunction<STATE> *bcdfunction = new TPZDummyFunction<STATE>(BCDfunction);
    TPZAutoPointer<TPZFunction<STATE> > bcD;
    bcdfunction->SetPolynomialOrder(20);
    bcD = bcdfunction;   */ 
    
    TPZManVector<REAL,4> Bottom     = fSimulationData->GetBottomBC();
    TPZManVector<REAL,4> Right      = fSimulationData->GetRightBC();
    TPZManVector<REAL,4> Top        = fSimulationData->GetTopBC();
    TPZManVector<REAL,4> Left       = fSimulationData->GetLeftBC();
    
    
    // Bc Bottom
    val2(0,0) = Bottom[1];
    val2(1,0) = Bottom[2];
    val2(2,0) = Bottom[3];
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, int(Bottom[0]), val1, val2);
    
    // Bc Right
    val2(0,0) = Right[1];
    val2(1,0) = Right[2];
    val2(2,0) = Right[3];
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, int(Right[0]), val1, val2);
//     mat->SetTimeDependentFunctionExact(bcD);
    
    // Bc Top
    val2(0,0) = Top[1];
    val2(1,0) = Top[2];
    val2(2,0) = Top[3];
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, int(Top[0]), val1, val2);
    
    // Bc Left
    val2(0,0) = Left[1];
    val2(1,0) = Left[2];
    val2(2,0) = Left[3];
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, int(Left[0]), val1, val2);
//     mat->SetTimedependentBCForcingFunction(bcN);
    
    cmesh->InsertMaterialObject(bcBottom);
    cmesh->InsertMaterialObject(bcRight);
    cmesh->InsertMaterialObject(bcTop);
    cmesh->InsertMaterialObject(bcLeft);
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    return cmesh;
}

void TPZDarcyAnalysis::ApplyStaticCondensation(){
    
    if (!fcmeshdarcy) {
        std::cout<< "No multiphysic computational mesh " << std::endl;
        DebugStop();
    }
    
    fcmeshinitialdarcy->Reference()->ResetReference();
    fcmeshinitialdarcy->LoadReferences();
    
    fcmeshinitialdarcy->ComputeNodElCon();
    // create condensed elements
    // increase the NumElConnected of one pressure connects in order to prevent condensation
    for (long icel=0; icel < fcmeshinitialdarcy->NElements(); icel++) {
        TPZCompEl  * cel = fcmeshinitialdarcy->Element(icel);
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
    
    int DOF = fmeshvec[0]->NEquations() + fmeshvec[1]->NEquations();// + fmeshvec[2]->NEquations() + fmeshvec[3]->NEquations();
    REAL PercentCondensedDOF = 100.0*(1.0 - REAL(fcmeshinitialdarcy->NEquations())/REAL(DOF));
    std::cout << "Percent of condensed Degrees of freedom: " << PercentCondensedDOF << std::endl;
    std::cout << "Condensed degrees of freedom: " << fcmeshinitialdarcy->NEquations() << std::endl;
    
    fcmeshinitialdarcy->CleanUpUnconnectedNodes();
    fcmeshinitialdarcy->ExpandSolution();
    
}

void TPZDarcyAnalysis::CmeshH1(int porder)
{
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,2,0.), val2(1,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    int dim = fgmesh->Dimension();
    int ilayer = 0;
    
    
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    // Material medio poroso
    TPZAxiSymmetricDarcyFlowH1 * mat = new TPZAxiSymmetricDarcyFlowH1(RockId);
    mat->SetReservoirData(fLayers[ilayer]);
    cmesh->InsertMaterialObject(mat);
    
    // Rigth hand side function
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Ffunction);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(0);
    forcef = dum;
    mat->SetForcingFunction(forcef);
    
    // Bc Bottom
    val2(0,0) = 0.0;
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    val2(0,0) = 1.0;
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    val2(0,0) = 0.0;
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    val2(0,0) = -1.0;
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef DEBUG
    std::ofstream out("cmeshPressureH1.txt");
    cmesh->Print(out);
#endif
    
    fcmeshdarcy = cmesh;
    
}



TPZCompMesh * TPZDarcyAnalysis::CmeshFlux(int qorder)
{
    
    int dim = fgmesh->Dimension();
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    int ilayer = 0;
    
    
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    // Setando Hdiv
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(qorder);
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->AutoBuild();
    
    
#ifdef DEBUG
    std::ofstream out("cmeshFlux.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CmeshPressure(int porder)
{
    
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    int dim = fgmesh->Dimension();
    
    
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    
    
    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
    cmesh->InsertMaterialObject(mat);
    
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
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    //    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef DEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CmeshSw(int Sworder)
{
    
    int dim = fgmesh->Dimension();
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    // Void material
    int matIdL2Proj = fSimulationData->fMatL2;
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,mat->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(Sworder);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
#ifdef DEBUG
    std::ofstream out("cmeshSw.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CmeshSo(int Soorder)
{
    
    int dim = fgmesh->Dimension();
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    // Void material
    int matIdL2Proj = fSimulationData->fMatL2;
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,mat->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(Soorder);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
#ifdef DEBUG
    std::ofstream out("cmeshSo.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}


void TPZDarcyAnalysis::ReadGeoMesh(std::string GridFileName)
{
    TPZReadGIDGrid GeometryInfo;
    GeometryInfo.SetfDimensionlessL(1.0);
    fgmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    fgmesh->SetDimension(3);
}

void TPZDarcyAnalysis::CreatedGeoMesh()
{
    
    
    long Qnodes = 4;
    int ilayer = 0;
    
    TPZGeoMesh *gmesh= new TPZGeoMesh;
    
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolLine(2);
    REAL r     = fLayers[ilayer]->Layerr();
    REAL rw    = fLayers[ilayer]->Layerrw();
    REAL h     = fLayers[ilayer]->Layerh();
    REAL top   = fLayers[ilayer]->LayerTop();
    
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    // Nodes
    long id = 0;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,  rw);         //coord r
    Node[id].SetCoord(1 , top - h);     //coord z
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,  rw + r);         //coord r
    Node[id].SetCoord(1 , top - h);     //coord z
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,  rw + r);         //coord r
    Node[id].SetCoord(1 ,  top);     //coord z
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 , rw);         //coord r
    Node[id].SetCoord(1 , top);     //coord z
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    
    //  Geometric Elements
    int elid = 0;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,bottomId,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,rigthId,*gmesh);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,topId,*gmesh);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,leftId,*gmesh);
    id++;
    
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 2;
    TopolQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elid,TopolQuad,RockId,*gmesh);
    
    
    
    
    gmesh->BuildConnectivity();
    fgmesh = gmesh;
    
}

void TPZDarcyAnalysis::Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0];//cos(par[0]);
    X[1] = 0.0 * par[0];//sin(par[0]);
    X[2] = 0.0;
}

void TPZDarcyAnalysis::Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;//par[0];
    X[1] = par[0];
    X[2] = 0.0;
}

void TPZDarcyAnalysis::GeometryLine(int nx, int ny)
{
    REAL t=0.0;
    REAL dt = fSimulationData->GetLengthElementx();
    int n = nx;
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
    GeoMesh1->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh1->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
    GeoMesh1->BuildConnectivity();
    GeoMesh1->SetDimension(0);
    
    TPZHierarquicalGrid *CreateGridFrom = new TPZHierarquicalGrid(GeoMesh1);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc = new TPZDummyFunction<STATE>(Parametricfunction);
    CreateGridFrom->SetParametricFunction(ParFunc);
    CreateGridFrom->SetFrontBackMatId(5,3);
    
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh2 = CreateGridFrom->ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid * CreateGridFrom2 = new TPZHierarquicalGrid(GeoMesh2);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc2 = new TPZDummyFunction<STATE>(Parametricfunction2);
    CreateGridFrom2->SetParametricFunction(ParFunc2);
    CreateGridFrom2->SetFrontBackMatId(2,4);
    //    CreateGridFrom2->SetTriangleExtrusion();
    
    dt = fSimulationData->GetLengthElementy();
    n = ny;
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fgmesh = CreateGridFrom2->ComputeExtrusion(t, dt, n);
    
    
}

void TPZDarcyAnalysis::PrintGeoMesh()
{
    
    
#ifdef DEBUG
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    fgmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fgmesh,Dummyfile, true);
    
#endif
}

void TPZDarcyAnalysis::RotateGeomesh(REAL CounterClockwiseAngle)
{
    REAL theta = CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    RotationMatrix(0,0) =   +cos(theta);
    RotationMatrix(0,1) =   -sin(theta);
    RotationMatrix(1,0) =   +sin(theta);
    RotationMatrix(1,1) =   +cos(theta);
    RotationMatrix(2,2) = 1.0;
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = fgmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = fgmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        fgmesh->NodeVec()[inode] = GeoNode;
    }
}

void TPZDarcyAnalysis::UniformRefinement(int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = fgmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fgmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
    fgmesh->BuildConnectivity();
}

void TPZDarcyAnalysis::UniformRefinement(int nh, std::set<int> &MatToRef)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = fgmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fgmesh->ElementVec() [i];
            if(!gel){continue;}
            //            int reflevel = gel->Level();
            //            if (reflevel == ref + 1) {
            //                continue;
            //            }
            TPZRefPatternTools::RefineDirectional(gel,MatToRef);
        }//for i
    }//ref
    fgmesh->BuildConnectivity();
}

void TPZDarcyAnalysis::UniformRefinement(int nh, int MatId)
{
    //    for ( int ref = 0; ref < nh; ref++ ){
    //        TPZVec<TPZGeoEl *> filhos;
    //        long n = fgmesh->NElements();
    //        for ( long i = 0; i < n; i++ ){
    //            TPZGeoEl * gel = fgmesh->ElementVec() [i];
    //            if(!gel){continue;}
    //            if (gel->Dimension() == 1){
    //                if (gel->MaterialId() == MatId) {
    //                    gel->Divide(filhos);
    //                }
    //
    //            }
    //        }//for i
    //    }//ref
    
    ///Refinamento
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
    for (int idivide = 0; idivide < nh; idivide++){
        const int nels = fgmesh->NElements();
        TPZVec< TPZGeoEl * > allEls(nels);
        for(int iel = 0; iel < nels; iel++){
            allEls[iel] = fgmesh->ElementVec()[iel];
        }
        
        for(int iel = 0; iel < nels; iel++){
            TPZGeoEl * gel = allEls[iel];
            if(!gel) continue;
            if(gel->HasSubElement()) continue;
            int nnodes = gel->NNodes();
            int found = -1;
            for(int in = 0; in < nnodes; in++){
                if(gel->NodePtr(in)->Id() == MatId){
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
    
    fgmesh->BuildConnectivity();
    
}

////refinamento uniforme em direcao ao no
//void DirectionalRef(int nh, int MatId){
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
//            allEls[iel] = gmesh->ElementVec()[iel];
//        }
//
//        for(int iel = 0; iel < nels; iel++){
//            TPZGeoEl * gel = allEls[iel];
//            if(!gel) continue;
//            if(gel->HasSubElement()) continue;
//            int nnodes = gel->NNodes();
//            int found = -1;
//            for(int in = 0; in < nnodes; in++){
//                if(gel->NodePtr(in)->Id() == nodeAtOriginId){
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
//    gmesh->BuildConnectivity();
//
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout<<"gmesh depois de refinar direcionalmente\n";
//        gmesh->Print(sout);
//        LOGPZ_DEBUG(logger, sout.str());
//    }
//#endif
//
//}///void



void TPZDarcyAnalysis::PostProcessVTK(TPZAnalysis *an)
{
    const int dim = 2;
    int div = fSimulationData->GetHPostrefinement();
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile;
    if (fSimulationData->GetGR()) {
        plotfile = "2DMixedDarcyGR.vtk";
    }
    else{
        plotfile = "2DMixedDarcy.vtk";
    }

    scalnames.Push("WeightedPressure");
    vecnames.Push("BulkVelocity");
    scalnames.Push("Porosity");
    scalnames.Push("Rhs");

    
    if (fSimulationData->IsOnePhaseQ()) {
        scalnames.Push("Rhoalpha");
        an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        an->PostProcess(div);
        return;
    }
    
    
    if (fSimulationData->IsTwoPhaseQ()) {
        scalnames.Push("Rhoalpha");
        scalnames.Push("Rhobeta");
        scalnames.Push("Salpha");
        an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        an->PostProcess(div);
        return;
    }
    
    
    if (fSimulationData->IsThreePhaseQ()) {
        scalnames.Push("Salpha");
        scalnames.Push("Sbeta");
        scalnames.Push("Rhoalpha");
        scalnames.Push("Rhobeta");
        scalnames.Push("Rhogamma");
        an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        an->PostProcess(div);
        return;
    }


}

void TPZDarcyAnalysis::BCDfunction(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &ff, TPZFMatrix<REAL> &Grad){

    REAL t = time;
    REAL pd;
    REAL x = pt[0];
    REAL y = pt[1];

    if(time <= 0.0){
        t = 0.0001;
    }
    
    pd = exp(-x/t) + x;
    ff[0] = pd;
}

void TPZDarcyAnalysis::BCNfunction(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &ff, TPZFMatrix<REAL> &Grad){
    
    REAL t = time;
    REAL qn;
    REAL x = pt[0];
    REAL y = pt[1];
    
    if(time <= 0.0){
        t = 0.0001;
    }
    
    qn = -1.*(1 - 1/(exp(x/t)*t))*(1 + (-0.1 + exp(-x/t) + x)/10.);
    ff[0] = -qn;
}

void TPZDarcyAnalysis::Ffunction(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &ff, TPZFMatrix<REAL> &Grad)
{
    
    REAL epsilon = 0.01;
    REAL xc = 0.5;
    REAL yc = 0.5;
    REAL x = pt[0];
    REAL y = pt[1];
    REAL rho = 1000.0/((1.0e7)/(1000.0*9.81)); //check this value
//    REAL c = 1.0;
    REAL Pref = 0.1;

    REAL f  =0.0;/* -((9.*exp(10.0*(-(1.0/1000.0) + log(1.0 + x))))/((1 + x)*(1 + x)))*/;
    REAL t = time;

    if(time <= 0.0){
        t = 0.0001;
    }
    REAL pD = 0.1 + log(1.0 + x);
    REAL dpDdx = 1.0/(1.0 + x);
    REAL dpDdx2 = -pow(1.0 + x,-2.0);
    REAL b = (1.0e7);
    REAL c = 1000.0;
    REAL e = 1.0e-6;
    REAL rhoD = (6.4290583420332425e-6*(101325.353 + b*pD))/
    (c*(0.5017121086995703 + 0.4982878913004297*
        exp(-5.1483395031556254e-8*(101325.353 + b*pD) -
              2.1184388752031904e-15*pow(101325.353 + b*pD,2.0) -
              3.2120945504363473e-47*pow(101325.353 + b*pD,6.0)) +
        1.2777444873399675e-8*pow(101325.353 + b*pD,1.0023796430661953)));
    REAL rhoDpe = (6.4290583420332425e-6*(101325.353 + b*pD+e))/
    (c*(0.5017121086995703 + 0.4982878913004297*
        exp(-5.1483395031556254e-8*(101325.353 + b*pD+e) -
            2.1184388752031904e-15*pow(101325.353 + b*pD+e,2.0) -
            3.2120945504363473e-47*pow(101325.353 + b*pD+e,6.0)) +
        1.2777444873399675e-8*pow(101325.353 + b*pD+e,1.0023796430661953)));
    REAL rhoDme = (6.4290583420332425e-6*(101325.353 + b*pD-e))/
    (c*(0.5017121086995703 + 0.4982878913004297*
        exp(-5.1483395031556254e-8*(101325.353 + b*pD-e) -
            2.1184388752031904e-15*pow(101325.353 + b*pD-e,2.0) -
            3.2120945504363473e-47*pow(101325.353 + b*pD-e,6.0)) +
        1.2777444873399675e-8*pow(101325.353 + b*pD-e,1.0023796430661953)));
    REAL drhoD = b * ((rhoDpe) - (rhoDme)) / (2.0 * e);
    
    f = -rhoD*dpDdx2 - dpDdx*dpDdx*drhoD;

//     f=-0.1*pow(1.0 - 1.0/(exp(x/t)*t),2.0) + (0.01*x)/(exp(x/t)*pow(t,2.0)) - (1.*(1.0 + (-0.1 + exp(-x/t) + x)/10.))/(exp(x/t)*pow(t,2.0));
    
#ifdef SolutionI
    
    REAL Pressure = exp(-((x-xc)*(x-xc)+(y-yc)*(y-yc))/epsilon);
    REAL Coefglob = exp(-2.0*(x*x-x*xc+xc*xc+(y-yc)*(y-yc))/epsilon);
    REAL Coef1 = (c*Pref - 1.0)*exp((x*x+xc*xc + (y-yc)*(y-yc))/epsilon);
    REAL Coef2 = (c)*exp((2.0*x*xc)/epsilon);
    REAL numerator1      = (x*x - 2.0*x*xc + xc*xc + y*y - 2.0*y*yc + yc*yc - epsilon) ;
    REAL numerator2      = (2.0*x*x - 4.0*x*xc + 2.0*xc*xc + 2.0*y*y - 4.0*y*yc + 2.0*yc*yc - epsilon) ;
    REAL denominator    =epsilon*epsilon;
    
    ff[0] = 4.0*rho*Coefglob*((Coef1*numerator1-Coef2*numerator2)/denominator);
    return;
#endif
    
#ifdef SolutionII
    
    REAL Pressure = x*y*(1-x)*(1-y)*(sin(M_PI*x))*(cos(M_PI*y));
    REAL cospix = cos(M_PI*x);
    REAL sinpiy = sin(M_PI*y);
    REAL cospiy = cos(M_PI*y);
    REAL sinpix = sin(M_PI*x);
    REAL x_1 = x-1;
    REAL y_1 = y-1;
    REAL x2_1 = 2.0*x-1.0;
    REAL y2_1 = 2.0*y-1.0;
    
    REAL cPlusSomething = (c*(x_1*x*y_1*y*sinpix*cospiy-Pref));
    REAL term1_1 = (2.0-M_PI*M_PI*x_1*x)*sinpix + 2.0*M_PI*x2_1*cospix;
    REAL term2_1 = 2.0*M_PI*y2_1*sinpiy + (M_PI*M_PI*y_1*y - 2.0)* cospiy;
    REAL term3_1 = M_PI*y_1*y*sinpiy + (1.0-2.0*y)*cospiy;
    REAL term3_2 = x2_1*sinpix+M_PI*x_1*x*cospix;
    
    REAL Term1 = -y_1*y*cospiy*(term1_1)*(cPlusSomething + 1.0) ;
    REAL Term2 = x_1*x*sinpix*(term2_1)*(cPlusSomething + 1.0);
    REAL Term3 = -c*x_1*x_1*x*x*sinpix*sinpix*(term3_1*term3_1) - c*y_1*y_1*y*y*cospiy*cospiy * (term3_2*term3_2);
    
    
    ff[0] = rho*(Term1+Term2+Term3);
    return;
#endif
    
    ff[0] = f;
    return;
    
}

TPZFMatrix<STATE> * TPZDarcyAnalysis::ComputeInverse()
{
    int neq = fcmeshdarcy->NEquations();
    TPZFMatrix<STATE> * PreInverse =  new TPZFMatrix<STATE> (neq,neq,0.0);
    TPZFStructMatrix skyl(fcmeshdarcy);
    std::set<int> matids; // to be computed
    matids.insert(1);
    matids.insert(2);
    matids.insert(3);
    matids.insert(4);
    matids.insert(5);
    skyl.SetMaterialIds(matids);
    TPZFMatrix<STATE> rhsfrac;
    TPZFMatrix<STATE> Identity;
    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
    TPZAutoPointer<TPZMatrix<STATE> > MatG = skyl.CreateAssemble(rhsfrac, gui);
    TPZFMatrix<STATE> oldmat = *MatG.operator->();
    oldmat.Inverse( * PreInverse);
    oldmat.Multiply(*PreInverse, Identity);
    
#ifdef DEBUG
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Is decomposed=  " << MatG->IsDecomposed() << std::endl;
        oldmat.Print("oldmat = ", sout,EMathematicaInput);
        PreInverse->Print("PreInverse = ", sout,EMathematicaInput);
        Identity.Print("Identity = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
#endif
    
    return PreInverse;
    
}

void TPZDarcyAnalysis::FilterSaturationGradients(TPZManVector<long> &active, TPZManVector<long> &nonactive)
{
    int ncon_flux       = fmeshvec[0]->NConnects();
    int ncon_pressure   = fmeshvec[1]->NConnects();
    int ncon_sw = fmeshvec[2]->NConnects();
    int ncon_so    = fmeshvec[3]->NConnects();
    int ncon = fcmeshdarcy->NConnects();
    
    
    // DOF related with the Q-P system
    for(int i = 0; i < ncon-ncon_sw-ncon_so; i++)
    {
        TPZConnect &con = fcmeshdarcy->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmeshdarcy->Block().Position(seqnum);
        int blocksize = fcmeshdarcy->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+blocksize);
        for(int ieq = 0; ieq<blocksize; ieq++)
        {
            active[vs+ieq] = pos+ieq;
        }
    }
    
    // DOF Related with S constant
    for(int i = ncon-ncon_sw-ncon_so; i< ncon; i++)
    {
        TPZConnect &con = fcmeshdarcy->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmeshdarcy->Block().Position(seqnum);
        int blocksize = fcmeshdarcy->Block().Size(seqnum);
        int vs = active.size();
        active.Resize(vs+1);
        
        int ieq = blocksize-1;
        active[vs] = pos+ieq;
    }
    
    // DOF Related with S grandients
    for(int i = ncon-ncon_sw-ncon_so; i< ncon; i++)
    {
        TPZConnect &con = fcmeshdarcy->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmeshdarcy->Block().Position(seqnum);
        int blocksize = fcmeshdarcy->Block().Size(seqnum);
        int vs = nonactive.size();
        nonactive.Resize(vs+blocksize-1);
        for(int ieq = 0; ieq<blocksize-1; ieq++)
        {
            nonactive[vs+ieq] = pos+ieq;
        }
        
    }
    
}

void TPZDarcyAnalysis::LinearTracer(const TPZVec< REAL >& pt, REAL time, TPZVec< STATE >& Saturation, TPZFMatrix< STATE >& Grad)
{
    REAL x = pt[0];
    REAL v = 0.00001;
    REAL Porosity = 0.1;
    REAL xshock = v*time/Porosity;
    
    if(x <= xshock)
    {
        Saturation[0] = 1.0;
    }
    else
    {
        Saturation[0] = 0.0;
    }
    
    return;
    
}


void TPZDarcyAnalysis::BluckleyAndLeverett(const TPZVec< REAL >& pt, REAL time, TPZVec< STATE >& Saturation, TPZFMatrix< STATE >& Grad)
{
    
    //definir globalmente:
    REAL Porosity = 0.1;
    
    REAL x = pt[0];
    
    REAL xshock = 0.0;
    REAL Swshock = 0.0;
    REAL epsilon = 1.0e-3;
    REAL ds = 1.0e-1;
    REAL qt = 0.00001;
    REAL muw=0.001;
    REAL muo= 0.001;//FluidData->GetMuOil();
    
    
    REAL lambdaw, lambdao;
    REAL lambdat, fw, dfwdSw;
    
    Swshock = SwatShock(epsilon,ds);
    dfwdSw = -((2.0*(Swshock-1.0) * Swshock * muw * muo)/((muw - 2.0 * Swshock * muw + Swshock * Swshock * (muw + muo))*(muw - 2.0 * Swshock * muw + Swshock * Swshock *(muw + muo))));//CUadratico
    // dfwdSw = -(( muw * muo)/(((Swshock *(muo-muw)+muw)*(Swshock *(muo-muw)+muw));//LInear
    //     dfwdSw = 1.0; //???
    
    // REAL j;
    // REAL xp;
    // int indexsat=0;
    // int nsx=(1.0-Swshock)/0.000001;
    //TPZFMatrix<REAL> xs(nsx,2,0.0);
    REAL Area=1;
    xshock      = (qt*time)/(Porosity*Area)*dfwdSw;
    //j=Ssaturation(0.1);
    //xp     = (qt*time)/(Porosity*Area)*j;
    //std::cout<<"RESPUE "<<j<<"RESPUE "<<xp<<std::endl;
    
    //for (int sx = 1; sx < nsx ; sx++) {
    //xs(sx,0)=Swshock+0.000001*sx;
    //xs(sx,1)=(qt*time)/(Porosity*Area)*Ssaturation(xs(sx,0));
    //}
    
    
    
    if(x <= xshock)
    {
        // for (int index = 0; index < nsx; index++) {
        //if (x< xs(index,1)){
        //indexsat=index;
        //}
        //}
        
        Saturation[0] =SaturationNewton(x,time,muo,muw,1,qt);
        //indexsat=0;
    }
    else
    {
        Saturation[0] = 0.0;
    }
    // std::cout<<x<<std::endl;
    return;
    
}


/**
 * Computes the saturation at shock using the Welge method
 */
REAL TPZDarcyAnalysis::SwatShock(REAL epsilon, REAL &ds){
    REAL muw= 0.001;
    REAL muo= 0.001;
    REAL value = 0.0;
    int npoints  = int(1.0/ds);
    int pos = 0;
    TPZFMatrix<REAL> points(npoints,2,0.0);
    TPZManVector<REAL> secants(npoints,0.0);
    REAL S = 0.0;
    for (int ip = 0; ip < npoints; ip++) {
        S = ds*ip;
        points(ip,0) = S;
        points(ip,1) = ((S*S*muo)/(muw - 2.0 * S * muw + S * S * (muw + muo))); //Cuadratico
        //points(ip,1) = ((S*muo)/(S*(muo-muw) +muw )); //Linear
    }
    
    
    for (int ip = 0; ip < npoints - 1; ip++) {
        secants[ip] = (points(ip + 1,1)-points(0,1))/(points(ip+1,0)-points(0,0));
        
    }
    
    // computing the max
    int maxpos =0;
    int value2=0;
    REAL maxval=0;
    for (int i=0;i<=npoints-1;i++){
        if (secants[maxpos]< secants[i]){
            maxval=secants[i];
            maxpos=i;
        }
        
    }
    
    
    pos = Extract(epsilon, secants, maxval);
    return points(pos,1);
    
}

/**
 * Extract a value from a given list
 */
int TPZDarcyAnalysis::Extract(REAL epsilon, TPZManVector<REAL> &list, REAL value){
    
    for (int i = 0; i<list.size(); i++) {
        if (fabs(list[i]-value) <= epsilon) {
            return i;
        }
    }
    return 0;
}

//REAL TPZDarcyAnalysis::Ssaturation( REAL S){

//REAL muw= 0.001;
// REAL dfwdSwS;
// REAL muo= 0.001;
// REAL deriv;
//    deriv = -((2.0*(S-1.0) * S * muw * muo)/((muw - 2.0 * S * muw + S * S * (muw + muo))*(muw - 2.0 * S * muw + S * S *(muw + muo))));
// deriv=-(( muw * muo)/(((Swshock *(muo-muw)+muw)*(Swshock *(muo-muw)+muw)); //Linear
//  return (deriv);
//}
REAL TPZDarcyAnalysis::SaturationNewton(REAL x,REAL t,REAL muo, REAL muw, REAL Area, REAL q)
{
    REAL Sin=0.8;
    REAL Swplusone=0;
    REAL fx;
    REAL derivatS;
    REAL deriv2atS;
    REAL sin2=0.0;
    
    
    
    while ( (Sin-sin2) > 0.0001){
        sin2=Sin;
        derivatS = dfdsw(Sin,muo,muw);
        deriv2atS = df2dsw(Sin,muo,muw);
        fx= dfdsw(Sin,muo,muw) - (Area * 0.1 * x)/(q * t);
        Swplusone= Sin - (fx)/(df2dsw(Sin,muo,muw));
        Sin = Swplusone;
        
    }
    
    return (Swplusone);
}

REAL TPZDarcyAnalysis::dfdsw(REAL Sw,REAL muo,REAL muw)
{
    REAL dfwdSwS;
    
    dfwdSwS = -((2.0*(Sw-1.0) * Sw * muw * muo)/((muw - 2.0 * Sw * muw + Sw * Sw * (muw + muo))*(muw - 2.0 * Sw * muw + Sw * Sw *(muw + muo))));
    
    return (dfwdSwS);
}
REAL TPZDarcyAnalysis::df2dsw(REAL Sw, REAL muo,REAL muw)
{
    REAL dfw2dSwS2;
    dfw2dSwS2 = (2* muo * muw *(muw + (Sw * Sw)*((2 * Sw) - 3)*(muo + muw)))/((muw - (2 * Sw * muw)+(Sw*Sw)*(muo + muw))*(muw - (2*Sw*muw)+(Sw*Sw)*(muo + muw))*(muw - (2*Sw*muw)+(Sw*Sw)*(muo + muw)));
    return (dfw2dSwS2);
}


// L2 projection

void TPZDarcyAnalysis::SolveProjection(TPZAnalysis *an, TPZCompMesh *Cmesh)
{
    TPZSkylineStructMatrix full(Cmesh);
    an->SetStructuralMatrix(full);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    //step.SetDirect(ELU);
    an->SetSolver(step);
    an->Run();
    
}


TPZCompMesh * TPZDarcyAnalysis::L2ProjectionCmesh(TPZVec<STATE> &solini)
{
    /// criar materiais
    int dim = 2;
    TPZL2Projection *material;
    material = new TPZL2Projection(1, dim, 1, solini, fSimulationData->Getsorder());
    
    TPZCompMesh * cmesh = new TPZCompMesh(this->fgmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialWaterSaturation);
    material->SetForcingFunction(forcef);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(fSimulationData->Getsorder());
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
    return cmesh;
    
}

void TPZDarcyAnalysis::InitialWaterSaturation(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    //    REAL x = pt[0];
    REAL y = pt[1];
    
    disp.resize(1);
    disp[0] = 0.0;
    
    //     if (y >=  50.0 ) {
    //         disp[0] = 0.0;
    //     }
    //   disp[0]= (rand() / (double)RAND_MAX);
    
}


/**
 * Computes convergence rate for an element
 */
void TPZDarcyAnalysis::CheckElementConvergence(int wichelement)
{
    
    TPZFMatrix<STATE> UAtn          = falphaAtn;
    TPZFMatrix<STATE> UAtnPlusOne   = fcmeshdarcy->Solution();
    
    std::ofstream outsol("Solutions.txt");
    UAtn.Print("UAtn =",outsol,EMathematicaInput);
    UAtnPlusOne.Print("UAtnPlusOne =",outsol,EMathematicaInput);
    outsol.flush();
    
    
    long neq = fcmeshdarcy->NEquations();
    TPZElementMatrix Jacobian(fcmeshdarcy, TPZElementMatrix::EK),Residual(fcmeshdarcy, TPZElementMatrix::EF);
    
    int nsteps = 4;
    REAL du=0.0001;
    
    std::ofstream outfile("CheckElementConvergence.txt");
    
    TPZCompEl * icel = fcmeshdarcy->ElementVec()[wichelement];
    if (!icel) {
        std::cout << "There's no element with index = " << wichelement << std::endl;
        DebugStop();
    }
    icel->Print(outfile);
    
    this->SetLastState();
    fcmeshdarcy->LoadSolution(UAtn);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    
    icel->CalcStiff(Jacobian,Residual);
    TPZFNMatrix<9,REAL> ResidualAtUn = Residual.fMat;
    
    
    this->SetNextState();
    fcmeshdarcy->LoadSolution(UAtnPlusOne);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    
    icel->CalcStiff(Jacobian,Residual);
    TPZFNMatrix<9,REAL> JacobianAtUnPlusOne = Jacobian.fMat;
    TPZFNMatrix<9,REAL> ResidualAtUnPlusOne = Residual.fMat;
    
    int SizeOfElMat = JacobianAtUnPlusOne.Rows();
    TPZFNMatrix<9,STATE> DeltaU(SizeOfElMat,1,du),JacobianDeltaU(SizeOfElMat,1,0.0);
    TPZFNMatrix<9,STATE> ResidualU(SizeOfElMat,1,0.0),Alpha_ResidualAtUnPlusOne_DeltaU(SizeOfElMat,1,0.0);
    TPZFNMatrix<9,STATE> ResidualUPerturbed(SizeOfElMat,1,0.0);
    
    ResidualU = ResidualAtUnPlusOne + ResidualAtUn;
    
    std::ofstream outRoot("Root.txt");
    std::ofstream outPerturbedRoot("PerturbedRoot.txt");
    
    JacobianAtUnPlusOne.Print("JacobianAtUnPlusOne = ",outRoot,EMathematicaInput);
    ResidualAtUn.Print("ResidualAtUn = ",outRoot,EMathematicaInput);
    ResidualAtUnPlusOne.Print("ResidualAtUnPlusOne = ",outRoot,EMathematicaInput);
    outRoot.flush();
    
    REAL alpha = 0;
    TPZFNMatrix<4,REAL> alphas(nsteps,1,0.0),ElConvergenceOrder(nsteps-1,1,0.0);
    TPZFNMatrix<9,REAL> res(nsteps,1,0.0);
    TPZFMatrix<REAL> DeltaUGlobal(neq,1,du);
    
    for(int j = 0; j < nsteps; j++)
    {
        
        alpha = (1.0*j+1.0)/100.0;
        JacobianAtUnPlusOne.Multiply(alpha*DeltaU,Alpha_ResidualAtUnPlusOne_DeltaU);
        alphas(j,0) = log10(alpha);
        
        fcmeshdarcy->LoadSolution(UAtnPlusOne+alpha*DeltaUGlobal);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
        
        this->SetNextState();
        icel->CalcStiff(Jacobian,Residual);
        
        TPZFNMatrix<9,REAL> ResidualAtUnPlusOne_AlphaDeltaU = Residual.fMat;
        ResidualAtUnPlusOne_AlphaDeltaU.Print("ResidualAtUnPlusOne_AlphaDeltaU =",outPerturbedRoot,EMathematicaInput);
        ResidualUPerturbed = ResidualAtUnPlusOne_AlphaDeltaU + ResidualAtUn;
        STATE NormValue = Norm(ResidualUPerturbed-(ResidualU+Alpha_ResidualAtUnPlusOne_DeltaU));
        res(j) = log10(NormValue);
        
    }
    
    for(int j = 1; j < nsteps ; j++){
        ElConvergenceOrder(j-1,0)=(res(j,0)-res(j-1,0))/(alphas(j,0)-alphas(j-1,0));
    }
    
    ElConvergenceOrder.Print("CheckConv = ",outfile,EMathematicaInput);
    outfile.flush();
    
}




void TPZDarcyAnalysis::CheckGlobalConvergence(TPZAnalysis * an)
{
    
    long neq = fcmeshdarcy->NEquations();
    int nsteps = 5;
    REAL du=1.0;
    REAL alpha = 0;
    
    TPZFMatrix<STATE> U0(neq,1,0.5);//(rand()/double(RAND_MAX))
    TPZFMatrix<STATE> UPlusAlphaDU   = fcmeshdarcy->Solution();
    
    std::ofstream outfile("CheckConvergence.txt");
    
    fcmeshdarcy->LoadSolution(U0);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    
    this->SetLastState();
    an->AssembleResidual();
    TPZFNMatrix<9,REAL> ResidualAtU0n = an->Rhs();
    
    fcmeshdarcy->LoadSolution(U0);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    
    this->SetNextState();
    an->Assemble();
    TPZAutoPointer< TPZMatrix<REAL> > Jacobian;
    Jacobian =   an->Solver().Matrix();
    TPZFMatrix<REAL> JacobianAtU0 = *Jacobian.operator->();
    TPZFNMatrix<9,REAL> ResidualAtU0 = an->Rhs();
    
    
    int SizeOfGlobalMat = JacobianAtU0.Rows();
    TPZFNMatrix<9,STATE> DeltaU(SizeOfGlobalMat,1,du),JacobianDeltaU(SizeOfGlobalMat,1,0.0);
    TPZFNMatrix<9,STATE> ResidualU0(SizeOfGlobalMat,1,0.0),Alpha_JacAtU0_DeltaU(SizeOfGlobalMat,1,0.0);
    TPZFNMatrix<9,STATE> ResidualUPerturbed(SizeOfGlobalMat,1,0.0);
    
    ResidualU0 = ResidualAtU0 + ResidualAtU0n;
    
    //    JacobianAtU0.Print("JacobianAtU0 = ",outfile,EMathematicaInput);
    //    ResidualAtU0n.Print("ResidualAtU0n = ",outfile,EMathematicaInput);
    //    ResidualAtU0.Print("ResidualAtU0 = ",outfile,EMathematicaInput);
    //    ResidualU0.Print("ResidualU0 = ",outfile,EMathematicaInput);
    //    U0.Print("U0 = ",outfile,EMathematicaInput);
    
    TPZFNMatrix<4,REAL> alphas(nsteps,1,0.0),ElConvergenceOrder(nsteps-1,1,0.0);
    TPZFNMatrix<9,REAL> res(nsteps,1,0.0);
    this->SetNextState();
    for(int j = 0; j < nsteps; j++)
    {
        
        alpha = (j+1.0)/1000.0;
        
        JacobianAtU0.Multiply(alpha*DeltaU,Alpha_JacAtU0_DeltaU);
        
        fcmeshdarcy->LoadSolution(U0+alpha*DeltaU);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
        
        
        an->AssembleResidual();
        TPZFNMatrix<9,REAL> ResidualAtU0PlusAlphaDU = an->Rhs();
        ResidualUPerturbed = ResidualAtU0PlusAlphaDU + ResidualAtU0n;
        TPZFNMatrix<9,REAL> error = ResidualUPerturbed-(ResidualU0+Alpha_JacAtU0_DeltaU);
        
        //        outfile << alpha <<std::endl;
        //        error.Print("error = ",outfile,EMathematicaInput);
        
        REAL NormValue = Norm(error);
        res(j) = log(NormValue);
        alphas(j,0) = log(alpha);
        
    }
    
    for(int j = 1; j < nsteps ; j++){
        ElConvergenceOrder(j-1,0)=(res(j,0)-res(j-1,0))/(alphas(j,0)-alphas(j-1,0));
    }
    
    ElConvergenceOrder.Print("CheckConv = ",outfile,EMathematicaInput);
    outfile.flush();
    
}

void TPZDarcyAnalysis::CheckGlobalJacobian(TPZAnalysis * an)
{
    
    long neq = fcmeshdarcy->NEquations();
    REAL du=1.0;
    
    TPZFMatrix<STATE> U(neq,1,0.0);//(rand()/double(RAND_MAX))
    
    U(0,0) = -0.5 ;
    U(1,0) = -0.5 ;
    U(2,0) = 0.5;
    U(3,0) = 0.5;
    
    U(12,0) = 10.29367992;
    U(13,0) = 0.1;
    
    U(14,0) = 0.1;
    U(15,0) = 10.29367992;
    
    U(16,0) = 1.0;
    U(17,0) = 0.0;
    
    TPZFMatrix<STATE> DeltaU(neq,1,du);
    DeltaU = du*U;
    TPZFMatrix<STATE> Uk = U + DeltaU;
    
    
    Uk(16,0) = 1.0 -  0.0 * du;
    Uk(17,0) = 0.0 +  0.0 * du;
    Uk.Print("Uk = ");
    std::ofstream outfile("CheckJacobian.txt");
    
    fcmeshdarcy->LoadSolution(U);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    
    this->SetLastState();
    an->AssembleResidual();
    TPZFNMatrix<9,REAL> ResidualAtUn = an->Rhs();
    
    fcmeshdarcy->LoadSolution(U);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    
    this->SetNextState();
    an->Assemble();
    TPZAutoPointer< TPZMatrix<REAL> > Jacobian;
    Jacobian =   an->Solver().Matrix();
    TPZFMatrix<REAL> JacobianAtU = *Jacobian.operator->();
    TPZFNMatrix<9,REAL> ResidualAtU = an->Rhs();
    
    
    int SizeOfGlobalMat = JacobianAtU.Rows();
    TPZFNMatrix<9,STATE> JacobianDeltaU(SizeOfGlobalMat,1,0.0);
    TPZFNMatrix<9,STATE> ResidualU(SizeOfGlobalMat,1,0.0),Jac_DeltaU(SizeOfGlobalMat,1,0.0);
    TPZFNMatrix<9,STATE> Residual(SizeOfGlobalMat,1,0.0);
    
    ResidualU = ResidualAtU + ResidualAtUn;
    
    fcmeshdarcy->LoadSolution(Uk);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    
    an->AssembleResidual();
    TPZFNMatrix<9,REAL> ResidualAtUk = an->Rhs() + ResidualAtUn;
    JacobianAtU.Multiply(Uk-U,Jac_DeltaU);
    Residual = ResidualAtUk - (ResidualU + JacobianDeltaU);
    
    JacobianAtU.Print("JacobianAtU = ",outfile,EMathematicaInput);
    JacobianDeltaU.Print("JacobianDeltaU = ",outfile,EMathematicaInput);
    ResidualU.Print("ResidualU = ",outfile,EMathematicaInput);
    ResidualAtUk.Print("ResidualAtUk = ",outfile,EMathematicaInput);
    Residual.Print("Residual = ",outfile,EMathematicaInput);
    //    U.Print("U = ",outfile,EMathematicaInput);
    //    Uk.Print("Uk = ",outfile,EMathematicaInput);
    outfile << std::endl;
    outfile << "Norm DeltaU " << Norm(DeltaU) << std::endl;
    outfile << "Norm Residual " << Norm(Residual) << std::endl;
    outfile.flush();
    
}
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
#include "pzpoisson3d.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzfstrmatrix.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"

#include "math.h"
//#include <boost/cstdfloat.hpp> // For float_64_t. Must be first include!
#include <cmath>  // for pow function.
//#include <boost/math/special_functions.hpp> // For gamma function.

#include "TPZMultiphysicsInterfaceEl.h"


#ifdef PZDEBUG
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
    
    /** @brief unknowns saturations for n time step */
    fSAtn.Resize(0, 0);
    
    /** @brief unknowns saturations for n+1 time step */
    fSAtnplusOne.Resize(0, 0);
    
    /** @brief Store DOF associated with active */
    fActiveEquations.Resize(0);
    
    /** @brief Store DOF associated with  non active */
    fNonactiveEquations.Resize(0);
    
    /** @brief L2 norm */
    fL2_norm.Resize(1,0.0);

    /** @brief L2 norm */
    fL2_norm_s.Resize(1,0.0);
    
    /** @brief Hdiv norm */
    fl2_norm_flux.Resize(1,0.0);
    
}


TPZDarcyAnalysis::~TPZDarcyAnalysis()
{
    
}


void TPZDarcyAnalysis::SetFluidData(TPZVec< TPZAutoPointer<Phase> > PVTData){
    
    TPZStack<std::string> System =  fSimulationData->GetsystemType();
    int nphases = System.size();
    
    if (fSimulationData->IsOnePhaseQ()) {
        
        for(int iphase = 0; iphase < nphases; iphase++){
        
            if (!strcmp("Water", System[iphase].c_str())){
                falpha_fluid = PVTData[0];
                fbeta_fluid = PVTData[1];
                fgamma_fluid = PVTData[2];
            }
            
            if (!strcmp("Oil", System[iphase].c_str())){
                falpha_fluid = PVTData[1];
                fbeta_fluid  = PVTData[0];
                fgamma_fluid = PVTData[2];                
            }
            
            if (!strcmp("Gas", System[iphase].c_str())){
                falpha_fluid = PVTData[2];
                fbeta_fluid = PVTData[0];
                fgamma_fluid = PVTData[1];
            }

        }

    }
    
    if(fSimulationData->IsTwoPhaseQ()){
        
        fmeshvecini.Resize(3);
        fmeshvec.Resize(3);
        
        for(int iphase = 0; iphase < nphases; iphase++){
            
            switch (iphase) {
                case 0:
                {
                    if (!strcmp("Water", System[iphase].c_str())){
                        falpha_fluid = PVTData[0];
                    }
                    
                    if (!strcmp("Oil", System[iphase].c_str())){
                        falpha_fluid = PVTData[1];
                    }
                    
                    if (!strcmp("Gas", System[iphase].c_str())){
                        falpha_fluid = PVTData[2];
                    }
                    
                    fgamma_fluid = PVTData[2];
                    
                }
                    break;
                    
                case 1:
                {
                    if (!strcmp("Water", System[iphase].c_str())){
                        fbeta_fluid = PVTData[0];
                    }
                    
                    if (!strcmp("Oil", System[iphase].c_str())){
                        fbeta_fluid = PVTData[1];
                    }
                    
                    if (!strcmp("Gas", System[iphase].c_str())){
                        fbeta_fluid = PVTData[2];
                    }
                    
                }
                    fgamma_fluid = PVTData[2];                    
                    break;
                default:
                {
                    DebugStop();
                }
                    break;
            }
            
        }
        
    }
    
    if(fSimulationData->IsThreePhaseQ()){
        
        fmeshvecini.Resize(4);
        fmeshvec.Resize(4);
        
        std::cout << "System not impelmented " << System << std::endl;
        DebugStop();
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
    
//    fResidualAtn.Print("fResidualAtn = ", cout,EMathematicaInput);
    
//#ifdef PZDEBUG
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
    
    //#ifdef PZDEBUG
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
    
//    fResidualAtnplusOne.Print("fResidualAtnplusOne = ", cout,EMathematicaInput);
    
    // #ifdef PZDEBUG
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
    
    REAL TimeStep = fSimulationData->GetDeltaT();
    fcmeshinitialdarcy->Solution().Zero();
    
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
    an->LoadSolution(fcmeshinitialdarcy->Solution());
    
    int n_dt = 10;
    int n_sub_dt = 10;
    int i_time = 0;
    REAL dt = (fSimulationData->GetMaxTime())/REAL(n_sub_dt);
    fSimulationData->SetDeltaT(dt);
    
    if (fSimulationData->GetGR())
    {
        if (fSimulationData->IsImpesQ()) {
            FilterSaturations(fActiveEquations,fNonactiveEquations);
            an->StructMatrix()->EquationFilter().Reset();
            an->StructMatrix()->EquationFilter().SetActiveEquations(fActiveEquations);
        }
        else
        {
            
            FilterSaturationGradients(fActiveEquations,fNonactiveEquations);
            an->StructMatrix()->EquationFilter().Reset();
            an->StructMatrix()->EquationFilter().SetActiveEquations(fActiveEquations);
            CleanUpGradients(an);
            SaturationReconstruction(an);
            
        }
        
    }
    
    falphaAtn = an->Solution();
    falphaAtnplusOne = an->Solution();
    
    while (i_time < n_dt) {
        
        
        
        this->AssembleLastStep(an);
        this->AssembleNextStep(an);
        
        
        if (fSimulationData->IsImpesQ())
        {
            
            const clock_t tinia = clock();
            PicardIterations(an);
            const clock_t tenda = clock();
            const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
            std::cout << "Time for Picard: " << timea << std::endl;
            
        }
        else
        {
            const clock_t tinia = clock();
            NewtonIterations(an);
            const clock_t tenda = clock();
            const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
            std::cout << "Time for Newton: " << timea << std::endl;
            
        }
        i_time++;
    }

    fcmeshinitialdarcy->LoadSolution(an->Solution());
    falphaAtn = fcmeshinitialdarcy->Solution();
    falphaAtnplusOne = fcmeshinitialdarcy->Solution();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshinitialdarcy);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshdarcy);    
    fSimulationData->SetDeltaT(TimeStep);

    
}

void TPZDarcyAnalysis::RunAnalysis()
{
    
//  Creating ouput directory
    std::string ouputname("mkdir ");
    ouputname += fSimulationData->GetDirectory();
    system(ouputname.c_str());
    
    std::string dirname = PZSOURCEDIR;
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
#ifdef PZDEBUG
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

    GridFileName += fSimulationData->GetGIDFile();
    
    if(fLayers[0]->GetIsGIDGeometry())
    {
        ReadGeoMesh(GridFileName);
    }
    else
    {
        
        int nx = fSimulationData->GetnElementsx();
        int ny = fSimulationData->GetnElementsy();
        Geometry2D(nx,ny);
        //        CreatedGeoMesh();
    }
    

    REAL deg = fSimulationData->GetRotationAngle();
    RotateGeomesh( deg * M_PI/180.0);
//    ApplyShear(2.0);
    this->UniformRefinement(fSimulationData->GetHrefinement());
    
    // Not refine  for read properties map
//    this->ReadRasterized();

//    int hcont = 1;
//    std::set<int> matidstoRef;
//    //    matidstoRef.insert(2);
//    //    matidstoRef.insert(3);
//    //    matidstoRef.insert(4);
//    //    matidstoRef.insert(3);
//    matidstoRef.insert(5);
//    this->UniformRefinement(hcont, matidstoRef);
    
//    TPZCheckGeom check(fgmesh);
//    int isbad = check.PerformCheck();
//    if (isbad) {
//        DebugStop();
//    }
    
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
    
    // Analyses
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Computing the Initial conditions " << std::endl;
    
    TPZAnalysis * Initialan = CreateAnalysis(fcmeshinitialdarcy);
    this->PushInitialCmesh();
    const clock_t t_beg_steady = clock();
    this->InitializeSolution(Initialan);
    const clock_t t_end_steady = clock();
    const REAL time_steady = REAL(REAL(t_end_steady - t_beg_steady)/CLOCKS_PER_SEC);
    std::cout << "Time for computing initial condition (minutes): " << time_steady/60.0 << std::endl;
    delete Initialan;

    std::cout << std::endl;
    std::cout << std::endl;    
    std::cout << "Computing the transient solution. " << std::endl;
    
    TPZAnalysis * an = CreateAnalysis(fcmeshdarcy);
    this->PushCmesh();
    //    this->CheckGlobalJacobian(an);
    //    this->CheckGlobalConvergence(an);
    const clock_t t_beg_transient = clock();
    this->TimeForward(an);
    const clock_t t_end_transient = clock();
    const REAL time_transient = REAL(REAL(t_end_transient - t_beg_transient)/CLOCKS_PER_SEC);
    std::cout << "Time for computing transient solution (minutes): " << time_transient/60.0 << std::endl;
    std::cout << std::endl;    

    
    return;
}

/**
 * read the rasterized list of k values
 */
void TPZDarcyAnalysis::ReadRasterized(){
    
    if (fLayers.size() != 1) {
        std::cout << " It works for a single hydraulic unit and a given rasterized list of k values. " << std::endl;
        DebugStop();
    }
    
    int ilayer = 0;
    REAL kappa[465] = {0.2627313748963993,0.5458605764803388,0.715765724284004,0.8650888663781195,\
         0.8638456579795561,0.8496178285293305,0.8896307210608712,0.9189151855603647,\
         0.9534947969426283,0.9475550234828253,0.7788470393222213,0.6593608987936276,\
         0.5546090800257851,0.8696472971728519,0.626530988120453,0.599917119440096,\
         0.6421862049912516,0.8315222396169077,0.7389262363016852,0.6237222580348099,\
         0.8163274703011327,0.9067593701077447,0.8272400773551893,0.8005341191638274,\
         0.8159591122571139,0.8942812413666084,0.9032139239340637,0.8899990791048901,\
         0.8460263376001472,0.7646192098719957,0.9264665254627499,0.16870798416060412,\
         0.4725112809650981,0.680863799613224,0.7968965834791417,0.7779261442121742,\
         0.7952389722810572,0.8523805138594714,0.893176167234552,0.9983884335574178,\
         0.9160143659637168,0.8486969334192835,0.7631918224514227,0.48595634957178374,\
         0.6906713325352242,0.4744451606961968,0.44456211437517273,0.6453172483654112,\
         0.8032047149829634,0.6417718021917304,0.7926144212174233,0.8562022285661663,\
         0.8624643153144858,0.7568376461920987,0.7821622617183903,0.7428400405193848,\
         0.8542683488350677,0.9001749700709091,0.8839672161340822,0.8506308131503821,\
         0.7890689750437426,0.9354452527857078,0.09070816833962612,0.2776959204346625,\
         0.5281333456119348,0.5486232618104797,0.5929643613592412,0.7533382447739203,\
         0.8365411179666636,0.8922092273690027,0.964821806796206,0.8754489363661478,\
         0.9513306934340178,0.8436320103140252,0.6340823280228383,0.501197163643061,\
         0.37752095036375355,0.5577401233999447,0.8783958007182983,0.8693249838843355,\
         0.7326181047978636,0.9603554655124782,0.8179390367437149,0.7734598029284465,\
         0.8148079933695552,0.8098812045308039,0.5876231697209688,0.7101943088682199,\
         0.865273045400129,0.8078091905331982,0.8357123123676212,0.7931669582834515,\
         0.9278018233723179,0.0860116032783866,0.3082235933327194,0.4033060134450686,\
         0.46924210332443134,0.5721981766276821,0.7565613776590848,0.8369555207661847,\
         0.8993001197163644,0.9406022654019707,0.8649046873561103,0.9603554655124782,\
         0.8956625840316789,0.8173404549221844,0.5460908002578506,0.37784326365227006,\
         0.6088037572520489,0.9065751910857354,0.7808730085643246,0.48130582926604665,\
         0.9099364582374069,0.7017681186112901,0.7064646836725298,0.743392577585413,\
         0.774058384749977,0.5586149737544894,0.594806151579335,0.7137857997974031,\
         0.6970255087945484,0.8098812045308039,0.7619486140528594,0.8938668385670873,\
         0.09199742149369187,0.24615526291555392,0.20628050465052034,0.482410903398103,\
         0.714430426374436,0.8114006814623814,0.8593793166958283,0.8834607238235564,\
         0.8958928078091905,0.9061607882862143,0.8344691039690579,0.8819872916474814,\
         0.9134819044110877,0.5465512478128741,0.3217147066949075,0.45874389906989593,\
         0.7495165300672254,0.6702735058476839,0.1976701353715812,0.5513859471406207,\
         0.6442121742333549,0.7787549498112165,0.43166958283451523,0.5525370660281794,\
         0.691684317156276,0.6185192006630444,0.4561193480062621,0.5600423611750621,\
         0.7492863062897136,0.706234459895018,0.8634772999355373,0.0768486969334193,\
         0.07399392209227368,0.12441292936734508,0.5762961598673911,0.8544064831015746,\
         0.8720876692144764,0.8842895294225989,0.6512109770697118,0.7480891426466525,\
         0.9017404917579889,0.604982042545354,0.6707799981582098,0.8031126254719587,\
         0.493783958007183,0.27764987567916016,0.2636983147619486,0.48494336495073215,\
         0.5673174325444332,0.1565061239524818,0.21152960677778798,0.7881020351781931,\
         0.7225803480983516,0.13643061055345795,0.19117782484575008,0.645823740675937,\
         0.6812321576572429,0.3904595266599134,0.4968689566258403,0.6368910581084813,\
         0.5849986186573349,0.8171562759001749,0.08992540749608621,0.22575743622801364,\
         0.4279399576388249,0.7850630813150382,0.8752187125886361,0.8970439266967493,\
         0.9081867575283176,0.4692421033244314,0.6169536789759646,0.8710286398379224,\
         0.6330693434017866,0.5856432452343678,0.5731651164932314,0.45814531724836544,\
         0.26986831199926326,0.22115296067777881,0.3326273137489639,0.4797863523344691,\
         0.432176075145041,0.19532185284096143,0.8422046228934524,0.45865180955889123,\
         0.07887466617552262,0.1020351781932038,0.4372870430058016,0.5017036559535869,\
         0.4372409982502993,0.5249102127267704,0.582235933327194,0.4097983239708997,\
         0.7643889860944838,0.09959480615157934,0.5378487890229302,0.7589096601897044,\
         0.8600239432728612,0.8910120637259417,0.8437701445805322,0.6566903029744913,\
         0.7363016852380513,0.7970347177456488,0.892117137857998,0.6299843447831293,\
         0.6693526107376371,0.4570862878718114,0.37148908739294595,0.3189520213647666,\
         0.2584492126346809,0.33819872916474814,0.5869785431439359,0.7615802560088406,\
         0.5669030297449121,0.909568100193388,0.5586149737544894,0.22216594529883046,\
         0.10779077263099734,0.3205635878073488,0.2215213187217976,0.18514596187494245,\
         0.33921171378579984,0.4955797034717746,0.23008564324523437,0.7513122755318171,\
         0.082097799060687,0.6189336034625655,0.8290818675752832,0.8780274426742793,\
         0.882862142002026,0.8275623906437056,0.4292292107928907,0.8746201307671057,\
         0.9450225619301961,0.9868772446818308,0.9050557141541579,0.7226724376093563,\
         0.44727875494981123,0.47729993553734223,0.4333271940325997,0.4191454093378764,\
         0.5172667833133806,0.5686527304540013,0.8879270651072843,0.9285385394603554,\
         0.9286306289713603,0.873883414679068,0.6810940233907358,0.5403352058200571,\
         0.33465328299106734,0.08214384381618933,0.07836817386499678,0.20001841790220096,\
         0.35758357123123674,0.12155815452619947,0.7542130951284649,0.08651809558891242,\
         0.4165208582742426,0.6370752371304909,0.829588359885809,0.8181692605212266,\
         0.5974767473984713,0.7650336126715168,0.8776590846302605,0.803803296804494,\
         0.8178929919882124,0.6113362188046781,0.613960769868312,0.6143751726678331,\
         0.7266783313380606,0.5901556312735978,0.6387788930840779,0.7557325720600423,\
         0.4581453172483654,0.7453264573165117,0.9471406206833042,0.8600699880283634,\
         0.8371857445436965,0.8355741781011142,0.795745464591583,0.4058845197532001,\
         0.07468459342480892,0.07463854866930657,0.16677410442950547,0.35569573625564055,\
         0.12993830002762685,0.7528317524633945,0.19053319826871723,0.26761211897964826,\
         0.5002302237775117,0.6437056819228288,0.8142554563035269,0.5811308591951377,\
         0.8431715627590016,0.8521502900819596,0.6975780458605764,0.4515609172115296,\
         0.5262915553918409,0.6320563587807349,0.6653467170089327,0.7866746477576205,\
         0.8051385947140621,0.7199097522792155,0.8300488074408324,0.6135924118242932,\
         0.7312367621327931,0.9429965926880929,0.7898517358872824,0.821668661939405,\
         0.8304632102403535,0.6545261994658808,0.5057555944377936,0.17160880375725207,\
         0.13629247628695093,0.10963256285109127,0.32618104797863523,0.16672805967400317,\
         0.7518187678423429,0.5945298830463209,0.709549682291187,0.7091813242471683,\
         0.7248365411179667,0.8190441108757711,0.7951008380145501,0.7903582281978083,\
         0.7848789022930289,0.8799152776498755,0.47310986278662864,0.4613224053780275,\
         0.6246431531448569,0.6725757436228013,0.9852196334837461,1.,0.4774380698038493,\
         0.521548945575099,0.8850722902661387,0.7410442950547932,0.826549406022654,\
         0.8291739570862879,0.8924394511465145,0.8749424440556222,0.5626208674831936,\
         0.6351413573993923,0.6016207753936826,0.47831292015839394,0.1678331338060595,\
         0.19265125702182523,0.20167602910028548,0.7513122755318171,0.6518095588912423,\
         0.8308315682843723,0.7568836909476012,0.7225803480983516,0.7691315959112258,\
         0.7022285661663136,0.7044387144304264,0.7980477023667004,0.6929275255548394,\
         0.7323878810203518,0.8146238143475457,0.8349295515240814,0.6065936089879364,\
         0.8486969334192835,0.8253522423795929,0.5403812505755595,0.25158854406483105,\
         0.7916935261073763,0.834607238235565,0.835389999079105,0.8824477392025049,\
         0.9184547380053412,0.8700156552168707,0.8097430702642969,0.727921539736624,\
         0.7017681186112903,0.6497375448936367,0.47877336771341744,0.3129662031494613,\
         0.27834054701169536,0.7513122755318171,0.6926052122663229,0.8442766368910581,\
         0.7862602449580994,0.7252969886729902,0.6840408877428861,0.5662584031678791,\
         0.6214660650151947,0.7023206556773184,0.22603370476102774,0.7782945022561929,\
         0.9798323970899715,0.9222764527120361,0.3350676857905885,0.4847591859287227,\
         0.4377474905608251,0.7240537802744267,0.32567455566810943,0.6999723731466987,\
         0.7995211345427756,0.8883875126623076,0.9471406206833042,0.6204070356386407,\
         0.8906437056819229,0.900313104337416,0.8284372409982503,0.7828989778064279,\
         0.7434846670964177,0.6981766276821071,0.5492678883875127,0.40639101206372596,\
         0.815959112257114,0.7438990698959389,0.8479141725757436,0.7939957638824937,\
         0.7285201215581546,0.6370291923749885,0.4145409337876416,0.38221751542499305,\
         0.6014365963716732,0.3702919237498849,0.7680265217791693,0.947508978727323,\
         0.7764527120360991,0.09057003407311907,0.6038769684132977,0.6674647757620407,\
         0.9336495073211161,0.8903213923934064,0.7662307763145779,0.5184179022009393,\
         0.8639837922460633,0.8152223961690764,0.7673818952021366,0.9274795100838014,\
         0.9151855603646746,0.8853025140436505,0.8914264665254628,0.8589649138963072,\
        0.8241550787365319,0.6718850722902662,0.47117598305553,0.8522423795929643};
    
    TPZFMatrix<STATE> Kabsolute(2,2);
    Kabsolute.Zero();
    int i_k = 0;
    int64_t pos_id = 0;
    
    int nel = fgmesh->NElements();
    TPZManVector< TPZFMatrix<REAL> > kvector(nel);
    
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl * gel  = fgmesh->Element(iel);
        if (!gel) {
            std::cout<< "Gel not found." << std::endl;
            DebugStop();
        }
        
        if (gel->Dimension() != 2) {
            continue;
        }

    
        Kabsolute(0,0) = kappa[i_k];
        Kabsolute(1,1) = kappa[i_k];
        
        pos_id = gel->Id();
        
        kvector[pos_id] = Kabsolute;
        i_k++;
        
    }
    
    fLayers[ilayer]->SetKvector(kvector);
    std::cout<< "Assigned values of permeabilities =  " << i_k << std::endl;
}

TPZAnalysis * TPZDarcyAnalysis::CreateAnalysis(TPZCompMesh * cmesh){
    
    if (!cmesh) {
        std::cout << "NUll computational mesh " << std::endl;
        DebugStop();
    }
    
    
    bool mustOptimizeBandwidth = fSimulationData->GetOptband();
    TPZAnalysis *an = new TPZAnalysis(cmesh,mustOptimizeBandwidth);
    int numofThreads = fSimulationData->GetNthreads();
    
    bool IsDirecSolver = fSimulationData->GetIsDirect();
    
    if (IsDirecSolver) {
        
        TPZSkylineNSymStructMatrix skylnsym(cmesh);
//            TPZParFrontStructMatrix<TPZFrontNonSym<STATE> > skylnsym(cmesh);
//            skylnsym.SetDecomposeType(ELU);
//            skylnsym.SetQuiet(1);
        TPZStepSolver<STATE> step;
        skylnsym.SetNumThreads(numofThreads);
        step.SetDirect(ELU);
        an->SetStructuralMatrix(skylnsym);
        an->SetSolver(step);
        
    }
    else
    {
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
    
    return an;
    
}

void TPZDarcyAnalysis::CreateInterfaces()
{
    if (fSimulationData->IsOnePhaseQ()) {
        return;
    }
    
    if(fLayers.size() > 1){
      fgmesh->AddInterfaceMaterial(1,2, 1);
      fgmesh->AddInterfaceMaterial(2,1, 1);
        
      fgmesh->AddInterfaceMaterial(2,3, 1);
      fgmesh->AddInterfaceMaterial(3,2, 1);
        
      fgmesh->AddInterfaceMaterial(2,4, 1);
      fgmesh->AddInterfaceMaterial(4,2, 1);
        
      fgmesh->AddInterfaceMaterial(3,4, 1);
      fgmesh->AddInterfaceMaterial(4,3, 1);
      fgmesh->BuildConnectivity();
    }
    
    fgmesh->ResetReference();
    fcmeshinitialdarcy->LoadReferences();
    
    // Creation of interface elements
    int nel = fcmeshinitialdarcy->ElementVec().NElements();
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
    
    
    
    fgmesh->ResetReference();
    fcmeshdarcy->LoadReferences();
    nel = fcmeshdarcy->ElementVec().NElements();
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
    FGlobal =   fResidualAtn+fResidualAtnplusOne;
    
#ifdef PZDEBUG
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
    
    
#ifdef PZDEBUG
    std::string out_name = fSimulationData->GetDirectory();
    out_name += "/ComputationalMesh.txt";
    std::ofstream out(out_name.c_str());
    fcmesh->Print(out);
    
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
    int n_sub_dt = fSimulationData->GetNSubSteps();
    TPZManVector<REAL> reporting_times = fSimulationData->GetTimes();
    REAL tk = fSimulationData->GetTime();
    REAL current_dt = (reporting_times[0] - tk)/n_sub_dt;
    
    if (fSimulationData->GetGR())
    {
        if (fSimulationData->IsImpesQ()) {
            FilterSaturations(fActiveEquations,fNonactiveEquations);
            an->StructMatrix()->EquationFilter().Reset();
            an->StructMatrix()->EquationFilter().SetActiveEquations(fActiveEquations);
        }
        else
        {
            
            FilterSaturationGradients(fActiveEquations,fNonactiveEquations);
            an->StructMatrix()->EquationFilter().Reset();
            an->StructMatrix()->EquationFilter().SetActiveEquations(fActiveEquations);
            CleanUpGradients(an);
            SaturationReconstruction(an);
            
        }

    }
    
    if (fSimulationData->IsImpesQ()) {
        FilterSaturations(fActiveEquations,fNonactiveEquations);
        an->StructMatrix()->EquationFilter().Reset();
        an->StructMatrix()->EquationFilter().SetActiveEquations(fActiveEquations);
    }

    // Out file for initial condition
    this->PostProcessVTK(an);
    std::cout << "Reported time " << tk/(86400.0) << "; dt = " << current_dt/(86400.0) << std::endl;
    std::cout << std::endl;
    
    falphaAtn = an->Solution();
    falphaAtnplusOne = an->Solution();
    
    int i_time = 0;
    int i_sub_dt = 0;
    int idt = 0;
    
    TPZManVector<REAL> velocities(3,0.0);

    TPZFNMatrix<100,REAL> current_v(n_sub_dt*reporting_times.size(),4,0.0);
    TPZFNMatrix<100,REAL> accumul_v(n_sub_dt*reporting_times.size(),4,0.0);
    
    while (i_time < reporting_times.size()) {

        fSimulationData->SetDeltaT(current_dt);
        this->AssembleLastStep(an);
        this->AssembleNextStep(an);
        
        
        if (fSimulationData->IsImpesQ())
        {
            
            const clock_t tinia = clock();
            PicardIterations(an);
            const clock_t tenda = clock();
            const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
            std::cout << "Time for Broyden: " << timea << std::endl;
            
        }
        else
        {
            const clock_t tinia = clock();
            NewtonIterations(an);
            const clock_t tenda = clock();
            const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
            std::cout << "Time for Newton: " << timea << std::endl;
            
        }
        REAL scale = 1.0/fSimulationData->Time_Scale();
        REAL v_scale = fSimulationData->Velocity_Scale();
        
        // Computing the rates at reporting times
        IntegrateVelocities(velocities);
        current_v(idt,0) = idt*current_dt*scale;
        current_v(idt,1) += velocities[0]*v_scale;
        current_v(idt,2) += velocities[1]*v_scale;
        current_v(idt,3) += velocities[2]*v_scale;
        idt++;
        
        i_sub_dt++;

        if (i_sub_dt == n_sub_dt) {
            tk = reporting_times[i_time];
            this->fSimulationData->SetTime(tk);
            std::cout << "Reported time " << tk/(86400.0) << "; dt = " << current_dt/(86400.0) << "; File number = " << i_time + 1 << std::endl;

            if (fSimulationData->GetGR())
            {
                std::cout << "Writing Saturations with reconstruction." << std::endl;
                fcmesh->LoadSolution(an->Solution());
                TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmesh);
            }
            std::cout << std::endl;
            
            this->PostProcessVTK(an);
            
            if (i_time == reporting_times.size()-1) {
                break;
            }
            current_dt = (reporting_times[i_time+1] - reporting_times[i_time])/n_sub_dt;
            i_time++;
            i_sub_dt=0;


        }
        
    }
    
    
    current_v(0,2) = current_v(1,2);
    

    accumul_v = TrapezoidalRule(current_v);
    
    std::string out_name_current = fSimulationData->GetDirectory();
    std::string out_name_accumul = fSimulationData->GetDirectory();
    out_name_current += "/current_production.txt";
    out_name_accumul += "/accumulated_production.txt";
    
    std::ofstream q_out(out_name_current.c_str());
    std::ofstream qt_out(out_name_accumul.c_str());
    current_v.Print("q = ", q_out,EMathematicaInput);
    accumul_v.Print("qt = ", qt_out,EMathematicaInput);
    
    return;
    
}

TPZFMatrix<REAL> TPZDarcyAnalysis::RiemmanRule(TPZFMatrix<REAL>  current_production){
    
    int length  = current_production.Rows();
    TPZFMatrix<REAL> accumulated(length,4,0.0);
    REAL dt = 0.0;
    REAL intw = 0.0;
    REAL into = 0.0;
    REAL intg = 0.0;
    
    for (int i = 1; i < length ; i++) {
        dt = current_production(i,0) - current_production(i-1,0);
        
        intw = current_production(i,1)*dt;
        into = current_production(i,2)*dt;
        intg = current_production(i,3)*dt;
        
        accumulated(i,0) = current_production(i,0);
        accumulated(i,1) += intw + accumulated(i-1,1);
        accumulated(i,2) += into + accumulated(i-1,2);
        accumulated(i,3) += intg + accumulated(i-1,3);
        
    }
    
    
    return accumulated;
    
}

TPZFMatrix<REAL> TPZDarcyAnalysis::TrapezoidalRule(TPZFMatrix<REAL>  current_production){
    
    int length  = current_production.Rows();
    TPZFMatrix<REAL> accumulated(length,4,0.0);
    REAL dt = 0.0;
    REAL intw = 0.0;
    REAL into = 0.0;
    REAL intg = 0.0;
    
    for (int i = 1; i < length ; i++) {
        dt = current_production(i,0) - current_production(i-1,0);
        
        intw = current_production(i,1)*dt - 0.5*(current_production(i,1) - current_production(i-1,1))*dt;
        into = current_production(i,2)*dt - 0.5*(current_production(i,2) - current_production(i-1,2))*dt;
        intg = current_production(i,3)*dt - 0.5*(current_production(i,3) - current_production(i-1,3))*dt;;
        
        accumulated(i,0) = current_production(i,0);
        accumulated(i,1) += intw + accumulated(i-1,1);
        accumulated(i,2) += into + accumulated(i-1,2);
        accumulated(i,3) += intg + accumulated(i-1,3);
        
    }
    
    
    return accumulated;
    
}


void TPZDarcyAnalysis::IntegrateL2SError(TPZManVector<REAL> & l2_norm){
    
    int mat_id = 1;
    int int_order = 36;
    int int_typ = 0;
    l2_norm[0] = 0.0;
    TPZManVector<STATE> saturation;
    TPZManVector<STATE> s_exact;
    
    fcmesh->Reference()->ResetReference();
    fcmesh->LoadReferences();
    
    int64_t n_elements = fcmesh->NElements();
    
    for (int64_t iel = 0 ; iel < n_elements; iel++) {
        TPZCompEl * cel = fcmesh->Element(iel);
        if (!cel) {
            std::cout << "Computational element non-exist " << std::endl;
            DebugStop();
        }
        
        TPZGeoEl * gel = cel->Reference();
        if (!gel) {
            std::cout << "Geomtric element non-exist " << std::endl;
            DebugStop();
        }
        
        if (gel->Dimension() != 2) {
            continue;
        }
        
        
        if (gel->MaterialId() == mat_id && gel->NumInterfaces() == 0) {
            
            // Creating the integration rule
            int gel_side = gel->NSides() - 1;
            TPZIntPoints * NumericIntegral = gel->CreateSideIntegrationRule(gel_side, int_order);
            NumericIntegral->SetType(int_typ, int_order);
            
            int dimension   = NumericIntegral->Dimension();
            int npoints     = NumericIntegral->NPoints();
            
            if (dimension != gel->Dimension()) {
                std::cout << "Incompatible dimensions." << std::endl;
                DebugStop();
            }
            
            // compute the integrals
            TPZManVector<REAL,2> xi_eta_duplet(2,0.0);
            TPZManVector<REAL,3> x(3,0.0);
            
            REAL weight = 0.0;
            REAL s = 1.0;
            for (int it = 0 ; it < npoints; it++) {
                
                NumericIntegral->Point(it, xi_eta_duplet, weight);
                TPZFMatrix<REAL> jac;
                TPZFMatrix<REAL> axes;
                REAL detjac;
                TPZFMatrix<REAL> jacinv;
                gel->Jacobian(xi_eta_duplet, jac, axes, detjac, jacinv);
                gel->X(xi_eta_duplet, x);
                
                if (fSimulationData->IsAxisymmetricQ()) {
                    s = 2.0*M_PI*x[0];
                }
                
                cel->Solution(xi_eta_duplet, 2, saturation);
                cel->Solution(xi_eta_duplet, 10, s_exact);
//                l2_norm[0] += s*weight * detjac * (saturation[0]);
                l2_norm[0] += s*weight * detjac * (s_exact[0] - saturation[0])*(s_exact[0] - saturation[0]);
                
            }
        }
        
    }
    
    
}

void TPZDarcyAnalysis::IntegrateFluxPError(TPZManVector<REAL> & l2_norm_flux,TPZManVector<REAL> & l2_norm){
    
    int mat_id = 1;
    int int_order = 36;
    int int_typ = 0;
    l2_norm[0] = 0.0;
    l2_norm_flux[0] = 0.0;
    REAL s = 1.0;    
    TPZManVector<STATE> p;
    TPZManVector<STATE> p_exact;
    TPZManVector<STATE> u;
    TPZManVector<STATE> u_exact;
    TPZManVector<STATE> divu;
    TPZManVector<STATE> divu_exact;
    
    TPZFMatrix<REAL> normals;
    TPZManVector<REAL,3> n(3,0.0);
    
    TPZManVector<int> v_sides;
    
    fcmesh->Reference()->ResetReference();
    fcmesh->LoadReferences();
    
    int64_t n_elements = fcmesh->NElements();
    
    for (int64_t iel = 0 ; iel < n_elements; iel++) {
        TPZCompEl * cel = fcmesh->Element(iel);
        if (!cel) {
            std::cout << "Computational element non-exist " << std::endl;
            DebugStop();
        }
        
        TPZGeoEl * gel = cel->Reference();
        if (!gel) {
            std::cout << "Geomtric element non-exist " << std::endl;
            DebugStop();
        }
        
        if (gel->Dimension() != 2) {
            continue;
        }
        
        
        if (gel->MaterialId() == mat_id && gel->NumInterfaces() == 0) {
            
            // Creating the integration rule
            int gel_side = gel->NSides() - 1;
            TPZIntPoints * NumericIntegral = gel->CreateSideIntegrationRule(gel_side, int_order);
            NumericIntegral->SetType(int_typ, int_order);
            
            int dimension   = NumericIntegral->Dimension();
            int npoints     = NumericIntegral->NPoints();
            
            if (dimension != gel->Dimension()) {
                std::cout << "Incompatible dimensions." << std::endl;
                DebugStop();
            }
            
            // compute the integrals
            TPZManVector<REAL,2> xi_eta_duplet(2,0.0);
            TPZManVector<REAL,3> x(3,0.0);
            
            REAL weight = 0.0;
            REAL dot = 0.0;
            for (int it = 0 ; it < npoints; it++) {
                
                NumericIntegral->Point(it, xi_eta_duplet, weight);
                TPZFMatrix<REAL> jac;
                TPZFMatrix<REAL> axes;
                REAL detjac;
                TPZFMatrix<REAL> jacinv;
                gel->Jacobian(xi_eta_duplet, jac, axes, detjac, jacinv);
                gel->X(xi_eta_duplet, x);
                
                if (fSimulationData->IsAxisymmetricQ()) {
                    // Computing the radius
                    TPZFMatrix<REAL> x_spatial(3,1,0.0);
                    x_spatial(0,0) = x[0];
                    REAL r = Norm(x_spatial);
                    s = 2.0*M_PI*r;
                }
                
                cel->Solution(xi_eta_duplet, 0, p);
                cel->Solution(xi_eta_duplet, 10, p_exact);
                
                cel->Solution(xi_eta_duplet, 1, u);
                cel->Solution(xi_eta_duplet, 11, u_exact);
                
                cel->Solution(xi_eta_duplet, 9, divu);
                cel->Solution(xi_eta_duplet, 12, divu_exact);

                dot = (u[0]-u_exact[0])*(u[0]-u_exact[0]) + (u[1]-u_exact[1])*(u[1]-u_exact[1]);
                l2_norm[0] += s*weight * detjac * (p_exact[0] - p[0])*(p_exact[0] - p[0]);
                l2_norm_flux[0] += s*weight * detjac * (dot);
                
            }
            
        }
        
    }
    
    
}

void TPZDarcyAnalysis::IntegrateVelocities(TPZManVector<REAL> & velocities){
    
    TPZStack<int> MaterialsToIntegrate = fSimulationData->MaterialsToIntegrate();
    int n_materials = MaterialsToIntegrate.size();
    
    if (n_materials == 0) {
        std::cout << "There is not material ids to identify." << std::endl;
        DebugStop();
    }
    
    int int_order = 10;
    int int_typ = 0;
    velocities[0] = 0.0;
    velocities[1] = 0.0;
    velocities[2] = 0.0;
    TPZManVector<STATE> sol;
    TPZFMatrix<REAL> normals;
    TPZManVector<REAL,3> n(3,0.0);
    
    TPZManVector<int> v_sides;
    
    fcmesh->Reference()->ResetReference();
    fcmesh->LoadReferences();
    
    int64_t n_elements = fcmesh->NElements();
    
    for (int64_t iel = 0 ; iel < n_elements; iel++) {
        TPZCompEl * cel = fcmesh->Element(iel);
        if (!cel) {
            std::cout << "Computational element non-exist " << std::endl;
            DebugStop();
        }
        
        TPZGeoEl * gel = cel->Reference();
        if (!gel) {
            std::cout << "Geomtric element non-exist " << std::endl;
            DebugStop();
        }
        
        if (gel->Dimension() != 1) {
            continue;
        }
        
        bool IntegrateQ = false;
        for (int imat = 0; imat < n_materials; imat++) {
            int mat_id = MaterialsToIntegrate[imat];
            if (gel->MaterialId() == mat_id) {
                IntegrateQ = true;
                break;
            }
        }
        
        
        
        if (IntegrateQ && gel->NumInterfaces() == 0) {
            
            
            TPZGeoEl * gel_2D = GetVolElement(gel);
            TPZGeoElSide intermediate_side;
            TPZTransform<> afine_transformation = Transform_1D_To_2D(gel,gel_2D,intermediate_side);
            
            int itself_2d = gel_2D->NSides()-1;
            TPZGeoElSide gel_side_2D(gel_2D,itself_2d);
            TPZCompEl * cel_2D = gel_2D->Reference();
            if (!cel_2D) {
                DebugStop();
            }
            
            
            int gel_side = gel->NSides() - 1;
            TPZIntPoints * NumericIntegral = gel->CreateSideIntegrationRule(gel_side, int_order);
            NumericIntegral->SetType(int_typ, int_order);
            
            // Creating the integration rule
            int dimension   = NumericIntegral->Dimension();
            int npoints     = NumericIntegral->NPoints();
            
            if (dimension != gel->Dimension()) {
                std::cout << "Incompatible dimensions." << std::endl;
                DebugStop();
            }
            
            
            gel_2D->ComputeNormals(normals, v_sides);
            for (int i = 0 ; i < v_sides.size(); i++) {
                if(v_sides[i] ==  intermediate_side.Side()){
                    n[0] = normals(0,i);
                    n[1] = normals(1,i);
                    n[2] = normals(2,i);
                    break;
                }
            }
            
            // compute the integrals
            TPZManVector<REAL,1> xi_singlet(1,0.0);
            TPZManVector<REAL,2> xi_eta_duplet(2,0.0);
            TPZManVector<REAL,2> x(3,0.0);
            REAL weight = 0.0;
            for (int it = 0 ; it < npoints; it++) {

                TPZFMatrix<REAL> jac;
                TPZFMatrix<REAL> axes;
                REAL detjac;
                TPZFMatrix<REAL> jacinv;
                gel->Jacobian(xi_singlet, jac, axes, detjac, jacinv);
                
                NumericIntegral->Point(it, xi_singlet, weight);
                afine_transformation.Apply(xi_singlet, xi_eta_duplet);
                gel_2D->X(xi_eta_duplet, x);
                
                REAL cross_area = (detjac*2.0)*(1.0);
                if (fSimulationData->IsAxisymmetricQ()) {
                    REAL rw = fabs(x[0]);
                    cross_area *= 2.0*M_PI*rw;
                    weight *= 2.0*M_PI*rw;
                }
                

                
                cel_2D->Solution(xi_eta_duplet, 19, sol);
                velocities[0] += weight * detjac * (sol[0] * n[0] + sol[1] * n[1]);

                
                cel_2D->Solution(xi_eta_duplet, 20, sol);
                velocities[1] += weight * detjac * (sol[0] * n[0] + sol[1] * n[1]);
                
            }
        }
        
    }
    
    
}

TPZGeoEl * TPZDarcyAnalysis::GetVolElement(TPZGeoEl * gel){
    
    if (gel->Dimension() != 1) {
        DebugStop();
    }
    TPZGeoEl * gel_2D;
    TPZGeoElSide gel_side(gel,2);
    while (gel_side != gel_side.Neighbour()) {

        TPZGeoElSide gel_2d_side = gel_side.Neighbour();
        gel_2D =  gel_2d_side.Element();
        if(gel_2D->Dimension() ==2)
        {
            break;
        }
        gel_side = gel_2d_side;
    }
    

    
    return gel_2D;
}

TPZTransform<>  TPZDarcyAnalysis::Transform_1D_To_2D(TPZGeoEl * gel_o, TPZGeoEl * gel_d, TPZGeoElSide & intermediate_side){
    
    int itself_o = gel_o->NSides()-1;
    int itself_d = gel_d->NSides()-1;
    
    TPZGeoElSide gel_side_o(gel_o,itself_o);
    TPZGeoElSide gel_side_d(gel_d,itself_d);
    TPZGeoElSide neigh = gel_side_o.Neighbour();
    TPZGeoEl * gel_2D;
    while(neigh != neigh.Neighbour()){
        gel_2D =  neigh.Element();
        if(gel_2D->Dimension() ==2)
        {
            break;
        }
        neigh = neigh.Neighbour();
    }
    
    intermediate_side = neigh;
    TPZTransform<> t1 = gel_side_o.NeighbourSideTransform(neigh);
    TPZTransform<> t2 = neigh.SideToSideTransform(gel_side_d);
    TPZTransform<> t3 = t2.Multiply(t1);
    
    return t3;
}

void TPZDarcyAnalysis::CleanUpGradients(TPZAnalysis *an){
    
    int64_t numofdof = fNonactiveEquations.size();
    TPZFMatrix<REAL> SolToLoad = an->Solution();
    for(int64_t i=0; i < numofdof; i++)
    {
        SolToLoad(fNonactiveEquations[i],0) = 0.0;
    }
    an->LoadSolution(SolToLoad);
    
}

void TPZDarcyAnalysis::PrintSaturations(TPZAnalysis *an){
    
    int64_t numofdof = fNonactiveEquations.size();
    TPZFMatrix<REAL> SolToLoad = an->Solution();
    TPZFMatrix<REAL> gradS(numofdof,1,0.0);
    for(int64_t i=0; i < numofdof; i++)
    {
       gradS(i,0) = SolToLoad(fNonactiveEquations[i],0);
    }
    
    TPZManVector<int64_t> AverageS;
    
    int ncon_sw = fmeshvec[2]->NConnects();
    int ncon = fcmesh->NConnects();
    
    // DOF Related with S constant
    for(int i = ncon-ncon_sw; i< ncon; i++)
    {
        TPZConnect &con = fcmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmesh->Block().Position(seqnum);
        int blocksize = fcmesh->Block().Size(seqnum);
        int vs = AverageS.size();
        AverageS.Resize(vs+1);
        
        int ieq = blocksize-1;
        AverageS[vs] = pos+ieq;
    }
    
    int64_t numofSdof = AverageS.size();
    TPZFMatrix<REAL> S(numofSdof,1,0.0);
    for(int64_t i=0; i < numofSdof; i++)
    {
        S(i,0) = SolToLoad(AverageS[i],0);
    }

    std::cout << "S = " << S << std::endl;
    std::cout << "Gradients = " << gradS << std::endl;
}

void TPZDarcyAnalysis::SaturationReconstruction(TPZAnalysis *an)
{
    fcmesh->LoadSolution(an->Solution());
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmesh);
    TPZGradientReconstruction *gradreconst = new TPZGradientReconstruction(false,1.);

    
    fmeshvec[2]->Reference()->ResetReference();
    fmeshvec[2]->LoadReferences();
    gradreconst->ProjectionL2GradientReconstructed(fmeshvec[2], fSimulationData->fMatL2);
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmesh);
    an->LoadSolution(fcmesh->Solution());
    

}

void TPZDarcyAnalysis::NewtonIterations(TPZAnalysis *an)
{
    
    TPZFMatrix<STATE> Residual(an->Rhs().Rows(),1,0.0);
    Residual = fResidualAtn + fResidualAtnplusOne;
    
    TPZFMatrix<STATE> X = falphaAtn;
    TPZFMatrix<STATE> DeltaX = falphaAtn;
    
    STATE error     =   1.0;
    STATE normdx    =   1.0;//Norm(Residual);
    int iterations  =   0;
    int centinel    =   0;
    int fixed       =   fSimulationData->GetFixediterations();
    
    while (iterations <= fSimulationData->GetMaxiterations()) {
        
        an->Rhs() = Residual;
        an->Rhs() *= -1.0;
        
//        this->PrintLS(an);
        
        an->Solve();
        
        DeltaX = an->Solution();
        normdx = Norm(DeltaX);
        X += DeltaX;
        
//        X.Print("X = ");
        
        falphaAtnplusOne=X;
        an->LoadSolution(X);
        
        if (fSimulationData->GetGR())
        {
            CleanUpGradients(an);
            this->SaturationReconstruction(an);
            falphaAtnplusOne = an->Solution();
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
//        Residual.Print("Residual = ");
        
#ifdef PZDEBUG
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
            std::cout << "Converged; iterations:  " << iterations << "; error: " << error <<  "; dx: " << normdx << std::endl;
            break;
        }
        
        
        if (iterations == fSimulationData->GetMaxiterations()) {
            std::cout << "Out; iterations:  " << iterations << "; error: " << error <<  "; dx: " << normdx << std::endl;
            break;
        }
        
    }
    
    this->UpDateAlphaVec(X);
    
}

void TPZDarcyAnalysis::UpdateSaturations(TPZAnalysis *an)
{
    
    
}

void TPZDarcyAnalysis::ComputeSaturations(TPZAnalysis *an){
    
    REAL scale = 0.2 * 1.0 /(fSimulationData->GetDeltaT());
    SetLastState();
    an->AssembleResidual();
    TPZFMatrix<REAL> residual = an->Rhs();

    
    SetNextState();
}


void TPZDarcyAnalysis::PicardIterations(TPZAnalysis *an)
{
    
    TPZFMatrix<STATE> Residual(an->Rhs().Rows(),1,0.0);
    Residual = fResidualAtn + fResidualAtnplusOne;
    
    TPZFMatrix<STATE> X = falphaAtn;
    TPZFMatrix<STATE> DeltaX = falphaAtn;
    
    STATE error     =   1.0;
    STATE normdx    =   1.0;//Norm(Residual);
    int iterations  =   0;
    int centinel    =   0;
    int fixed       =   fSimulationData->GetFixediterations();
    
    while (iterations <= fSimulationData->GetMaxiterations()) {
        
        an->Rhs() = Residual;
        an->Rhs() *= -1.0;
        
        this->PrintLS(an);
        
        an->Solve();
        
        
        DeltaX = an->Solution();
        normdx = Norm(DeltaX);
        X += DeltaX;
        
//        DeltaX.Print("DeltaX = ", std::cout,EMathematicaInput);
        X.Print("X = ", std::cout,EMathematicaInput);
//        Residual.Print("Residual = ", std::cout,EMathematicaInput);
        
        // compute Saturations
        ComputeSaturations(an);
        
        // update Saturations
        UpdateSaturations(an);
        
        falphaAtnplusOne=X;
        
        if (fSimulationData->GetGR())
        {
            an->Solution() = falphaAtnplusOne;
            CleanUpGradients(an);
            this->SaturationReconstruction(an);
            fcmesh->LoadSolution(an->Solution());
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
        
        
#ifdef PZDEBUG
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
            std::cout << "Converged; iterations:  " << iterations << "; error: " << error <<  "; dx: " << normdx << std::endl;
            break;
        }
        
        
        if (iterations == fSimulationData->GetMaxiterations()) {
            std::cout << "Out; iterations:  " << iterations << "; error: " << error <<  "; dx: " << normdx << std::endl;
            break;
        }
        
    }
    
    this->UpDateAlphaVec(X);
    
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
    int n_hydraulic_units = fLayers.size();
    std::set<int> set;
    
    for (int ilayer = 0; ilayer < n_hydraulic_units; ilayer++) {
        
        int RockId = fLayers[ilayer]->GetMatIDs()[0]; // Volumetric material
        set.insert(RockId);
        
        // Material medio poroso
        TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
        mat->SetSimulationData(fSimulationData);
        mat->SetReservoirData(fLayers[ilayer]);
        mat->SetPetroPhysicsData(fRockPetroPhysic[0]);
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
        //    mat->SetTimeDependentForcingFunction(forcef);
        
        // Setting up linear tracer solution
        TPZDummyFunction<STATE> *Ltracer = new TPZDummyFunction<STATE>(LinearTracer);
        //    TPZDummyFunction<STATE> *Ltracer = new TPZDummyFunction<STATE>(BluckleyAndLeverett);
        TPZAutoPointer<TPZFunction<STATE> > fLTracer = Ltracer;
        mat->SetTimeDependentFunctionExact(fLTracer);
        
        TPZDummyFunction<STATE> * P_hydrostatic = new TPZDummyFunction<STATE>(P_Hydrostatic);
        TPZAutoPointer<TPZFunction<STATE> > P_hydrostatic_ptr;
        P_hydrostatic_ptr = P_hydrostatic;
        mat->SetTimedependentBCForcingFunction(P_hydrostatic_ptr);
        
        int n_boundaries = fLayers[ilayer]->GetInitialBC().size();
        TPZManVector<REAL,4> ibctype;
        
        for (int ibc = 0; ibc < n_boundaries; ibc++) {
            ibctype  = fLayers[ilayer]->GetInitialBC()[ibc];
            
            int BcId = fLayers[ilayer]->GetMatIDs()[ibc+1];
            set.insert(BcId);
            
            // Bc
            val2(0,0) = ibctype[1];
            val2(1,0) = ibctype[2];
            val2(2,0) = ibctype[3];
            TPZBndCond * bc_condition = mat->CreateBC(mat, BcId, int(ibctype[0]), val1, val2);
            cmesh->InsertMaterialObject(bc_condition);
        }
        

    }
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    cmesh->AutoBuild(set);
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
    int n_hydraulic_units = fLayers.size();
    std::set<int> set;
    
    for (int ilayer = 0; ilayer < n_hydraulic_units; ilayer++) {
        
        int RockId = fLayers[ilayer]->GetMatIDs()[0]; // Volumetric material
        set.insert(RockId);
        
        // Material medio poroso
        TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
        mat->SetSimulationData(fSimulationData);
        mat->SetReservoirData(fLayers[ilayer]);
        mat->SetPetroPhysicsData(fRockPetroPhysic[0]);
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
    //    TPZDummyFunction<STATE> *Load_function= new TPZDummyFunction<STATE>(Cylindrical_Elliptic);
    //    TPZDummyFunction<STATE> *Load_function = new TPZDummyFunction<STATE>(Dupuit_Thiem);
        TPZDummyFunction<STATE> *Load_function = new TPZDummyFunction<STATE>(LinearTracer);
//        TPZDummyFunction<STATE> *Load_function = new TPZDummyFunction<STATE>(BluckleyAndLeverett);
        TPZAutoPointer<TPZFunction<STATE> > fLoad_function = Load_function;
        mat->SetTimeDependentFunctionExact(fLoad_function);
        
        TPZDummyFunction<STATE> * P_hydrostatic = new TPZDummyFunction<STATE>(P_Hydrostatic);
        TPZAutoPointer<TPZFunction<STATE> > P_hydrostatic_ptr;
        P_hydrostatic_ptr = P_hydrostatic;
        mat->SetTimedependentBCForcingFunction(P_hydrostatic_ptr);
        
        int n_boundaries = fLayers[ilayer]->GetBC().size();
        TPZManVector<REAL,4> ibctype;
        
        for (int ibc = 0; ibc < n_boundaries; ibc++) {
            ibctype  = fLayers[ilayer]->GetBC()[ibc];
            
            int BcId = fLayers[ilayer]->GetMatIDs()[ibc+1];
            set.insert(BcId);
            
            // Bc
            val2(0,0) = ibctype[1];
            val2(1,0) = ibctype[2];
            val2(2,0) = ibctype[3];
            TPZBndCond * bc_condition = mat->CreateBC(mat, BcId, int(ibctype[0]), val1, val2);
            cmesh->InsertMaterialObject(bc_condition);
        }
    }
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    cmesh->AutoBuild(set);
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
    for (int64_t icel=0; icel < fcmeshinitialdarcy->NElements(); icel++) {
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


TPZCompMesh * TPZDarcyAnalysis::CmeshFlux(int qorder)
{
    
    int dim = fgmesh->Dimension();
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    int n_hydraulic_units = fLayers.size();
    std::set<int> set;
    
    for (int ilayer = 0; ilayer < n_hydraulic_units; ilayer++) {
        
        int RockId = fLayers[ilayer]->GetMatIDs()[0]; // Volumetric material
        set.insert(RockId);
        
        // Material medio poroso
        TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
        cmesh->InsertMaterialObject(mat);
        
        int n_boundaries = fLayers[ilayer]->GetBC().size();
        TPZManVector<REAL,4> ibctype;
        
        for (int ibc = 0; ibc < n_boundaries; ibc++) {
            ibctype  = fLayers[ilayer]->GetBC()[ibc];
            
            int BcId = fLayers[ilayer]->GetMatIDs()[ibc+1];
            set.insert(BcId);
            
            // Bc
            val2(0,0) = ibctype[1];
            val2(1,0) = ibctype[2];
            val2(2,0) = ibctype[3];
            TPZBndCond * bc_condition = mat->CreateBC(mat, BcId, int(ibctype[0]), val1, val2);
            cmesh->InsertMaterialObject(bc_condition);
        }
    }
    
    
    // Setando Hdiv
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(qorder);
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->AutoBuild(set);
    
    
#ifdef PZDEBUG
    std::string out_name = fSimulationData->GetDirectory();
    out_name += "/cmeshFlux.txt";
    std::ofstream out(out_name.c_str());
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
    
    int n_hydraulic_units = fLayers.size();
    std::set<int> set;
    
    for (int ilayer = 0; ilayer < n_hydraulic_units; ilayer++) {
        
        int RockId = fLayers[ilayer]->GetMatIDs()[0]; // Volumetric material
        set.insert(RockId);
        
        // Material medio poroso
        TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
        cmesh->InsertMaterialObject(mat);
    
    }
    
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild(set);
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::string out_name = fSimulationData->GetDirectory();
    out_name += "/cmeshPress.txt";
    std::ofstream out(out_name.c_str());
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CmeshSw(int Sworder)
{
    
    int dim = fgmesh->Dimension();
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    int n_hydraulic_units = fLayers.size();
    std::set<int> set;
    
    for (int ilayer = 0; ilayer < n_hydraulic_units; ilayer++) {
        
        int RockId = fLayers[ilayer]->GetMatIDs()[0]; // Volumetric material
        set.insert(RockId);
        
        // Material medio poroso
        TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
        cmesh->InsertMaterialObject(mat);
        
        // Void material
        int matIdL2Proj = fSimulationData->fMatL2;
        TPZVec<STATE> sol(1,0.);
        TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,mat->NStateVariables(),sol);
        cmesh->InsertMaterialObject(matl2proj);
        
        int n_boundaries = fLayers[ilayer]->GetBC().size();
        TPZManVector<REAL,4> ibctype;
        
        for (int ibc = 0; ibc < n_boundaries; ibc++) {
            ibctype  = fLayers[ilayer]->GetBC()[ibc];
            
            int BcId = fLayers[ilayer]->GetMatIDs()[ibc+1];
            set.insert(BcId);
            
            // Bc
            val2(0,0) = ibctype[1];
            val2(1,0) = ibctype[2];
            val2(2,0) = ibctype[3];
            TPZBndCond * bc_condition = mat->CreateBC(mat, BcId, int(ibctype[0]), val1, val2);
            cmesh->InsertMaterialObject(bc_condition);
        }
    }
    

    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(Sworder);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild(set);
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
#ifdef PZDEBUG
    std::string out_name = fSimulationData->GetDirectory();
    out_name += "/cmeshSw.txt";
    std::ofstream out(out_name.c_str());
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CmeshSo(int Soorder)
{
    
    int dim = fgmesh->Dimension();
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    int n_hydraulic_units = fLayers.size();
    std::set<int> set;
    
    for (int ilayer = 0; ilayer < n_hydraulic_units; ilayer++) {
        
        int RockId = fLayers[ilayer]->GetMatIDs()[0]; // Volumetric material
        set.insert(RockId);
        
        // Material medio poroso
        TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
        cmesh->InsertMaterialObject(mat);
        
        // Void material
        int matIdL2Proj = fSimulationData->fMatL2;
        TPZVec<STATE> sol(1,0.);
        TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,mat->NStateVariables(),sol);
        cmesh->InsertMaterialObject(matl2proj);
        
        int n_boundaries = fLayers[ilayer]->GetBC().size();
        TPZManVector<REAL,4> ibctype;
        
        for (int ibc = 0; ibc < n_boundaries; ibc++) {
            ibctype  = fLayers[ilayer]->GetBC()[ibc];
            
            int BcId = fLayers[ilayer]->GetMatIDs()[ibc+1];
            set.insert(BcId);
            
            // Bc
            val2(0,0) = ibctype[1];
            val2(1,0) = ibctype[2];
            val2(2,0) = ibctype[3];
            TPZBndCond * bc_condition = mat->CreateBC(mat, BcId, int(ibctype[0]), val1, val2);
            cmesh->InsertMaterialObject(bc_condition);
        }
    }
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(Soorder);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild(set);
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
#ifdef PZDEBUG
    std::string out_name = fSimulationData->GetDirectory();
    out_name += "/cmeshSo.txt";
    std::ofstream out(out_name.c_str());
    cmesh->Print(out);
#endif
    
    return cmesh;
}


void TPZDarcyAnalysis::ReadGeoMesh(std::string GridFileName)
{
    TPZReadGIDGrid GeometryInfo;
    REAL s = fSimulationData->Length_Scale();
    GeometryInfo.SetfDimensionlessL(s);
    fgmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    fgmesh->SetDimension(2);
}

void TPZDarcyAnalysis::CreatedGeoMesh()
{
    
    if (fLayers.size() != 1) {
        std::cout << " It works for a single hydraulic unit." << std::endl;
        DebugStop();
    }
    
    int64_t Qnodes = 4;
    int ilayer = 0;
    
    TPZGeoMesh *gmesh= new TPZGeoMesh;
    
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolLine(2);
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
    int64_t id = 0;
    
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

void TPZDarcyAnalysis::ParametricfunctionX(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = par[0];//cos(par[0]);
    X[1] = 0.0 * par[0];//sin(par[0]);
    X[2] = 0.0;
}

void TPZDarcyAnalysis::ParametricfunctionY(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;//par[0];
    X[1] = par[0];
    X[2] = 0.0;
}

void TPZDarcyAnalysis::ApplyPG(TPZGeoMesh * geomesh){

    std::cout << "Applyging geometric progression. " << std::endl;
    if (!geomesh) {
        DebugStop();
    }
    
    int n = 0;
    REAL dx = fSimulationData->GetLengthElementx();
    REAL ratio = fSimulationData->GetPGRatio();
    REAL l = REAL(fSimulationData->GetnElementsx())*dx;
    int n_nodes = geomesh->NNodes();
    REAL dr_0,dr_n;
    REAL r_powers = 0.0;
    for (int i=1; i <= n_nodes-1; i++) {
        r_powers += std::pow(ratio, REAL(i-1));
    }
    dr_0 = l/r_powers;
    TPZManVector<REAL,3> coor_n(3,0.0);
    TPZManVector<REAL,3> coor(3,0.0);
    for (int inode = 1; inode < n_nodes-1; inode++) {
        geomesh->NodeVec()[inode-1].GetCoordinates(coor);
        geomesh->NodeVec()[inode].GetCoordinates(coor_n);
        
        n = (n_nodes-2) - inode + 1;
        dr_n = dr_0*std::pow(ratio, REAL(n));
        coor_n[0] = dr_n + coor[0];
        geomesh->NodeVec()[inode].SetCoord(coor_n);
        
    }
    
}

void TPZDarcyAnalysis::Geometry2D(int nx, int ny)
{
    REAL t=0.0;
    REAL dt = fSimulationData->GetLengthElementx();
    int n = nx;
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
    GeoMesh1->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    coors[0] = fLayers[0]->Layerrw();
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    
    
    GeoMesh1->NodeVec()[0]=Node;
    
    TPZVec<int64_t> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
    GeoMesh1->BuildConnectivity();
    GeoMesh1->SetDimension(0);
    GeoMesh1->SetMaxNodeId(0);
    GeoMesh1->SetMaxElementId(0);
    
    TPZHierarquicalGrid *CreateGridFrom = new TPZHierarquicalGrid(GeoMesh1);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(ParametricfunctionX);
    CreateGridFrom->SetParametricFunction(ParFunc);
    CreateGridFrom->SetFrontBackMatId(5,3);
    
    // Computing Mesh extruded along the parametric curve ParametricfunctionX
    TPZGeoMesh * GeoMesh2 = CreateGridFrom->ComputeExtrusion(t, dt, n);
    if (fSimulationData->IsMeshwithPGQ()) {
        ApplyPG(GeoMesh2);
    }

    TPZHierarquicalGrid * CreateGridFrom2 = new TPZHierarquicalGrid(GeoMesh2);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(ParametricfunctionY);
    CreateGridFrom2->SetParametricFunction(ParFunc2);
    CreateGridFrom2->SetFrontBackMatId(2,4);
    
    if (fSimulationData->IsTriangularMeshQ())
    {
        CreateGridFrom2->SetTriangleExtrusion();
    }
    
    dt = fSimulationData->GetLengthElementy();
    n = ny;
    
    // Computing Mesh extruded along the parametric curve ParametricfunctionY
    fgmesh = CreateGridFrom2->ComputeExtrusion(t, dt, n);
    
    
}

void TPZDarcyAnalysis::PrintGeoMesh()
{
    
//#ifdef PZDEBUG
    
    //  Print Geometrical Base Mesh
    std::string out_name_text = fSimulationData->GetDirectory();
    std::string out_name_vtk = fSimulationData->GetDirectory();
    out_name_text += "/GeometicMesh.txt";
    out_name_vtk += "/GeometricMesh.vtk";
    std::ofstream argument(out_name_text.c_str());
    std::ofstream Dummyfile(out_name_vtk.c_str());
    
    fgmesh->Print(argument);
    TPZVTKGeoMesh::PrintGMeshVTK(fgmesh,Dummyfile, true);
    
//#endif
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

void TPZDarcyAnalysis::ApplyShear(REAL CounterClockwiseAngle)
{
    REAL theta = CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<STATE> ShearMatrix(3,3,0.0);
    ShearMatrix(0,0) =   1.0;
    ShearMatrix(0,1) =   sin(theta);
    ShearMatrix(1,0) =   0.0;
    ShearMatrix(1,1) =   1.0;
    ShearMatrix(2,2) = 1.0;
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    ShearMatrix.Print("ShearMatrix = ");
    
    int NumberofGeoNodes = fgmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = fgmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = ShearMatrix(0,0)*iCoords[0]+ShearMatrix(0,1)*iCoords[1]+ShearMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = ShearMatrix(1,0)*iCoords[0]+ShearMatrix(1,1)*iCoords[1]+ShearMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = ShearMatrix(2,0)*iCoords[0]+ShearMatrix(2,1)*iCoords[1]+ShearMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        fgmesh->NodeVec()[inode] = GeoNode;
    }
}

void TPZDarcyAnalysis::UniformRefinement(int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = fgmesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = fgmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
    fgmesh->BuildConnectivity();
}

void TPZDarcyAnalysis::UniformRefinement(int nh, std::set<int> &MatToRef)
{
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = fgmesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = fgmesh->ElementVec() [i];
            if(!gel || gel->HasSubElement())
            {
                continue;
            }
//                int reflevel = gel->Level();
//                if (reflevel == ref + 1) {
//                    continue;
//                }
            TPZRefPatternTools::RefineDirectional(gel,MatToRef);
        }//for i
    }//ref
    fgmesh->BuildConnectivity();
}

void TPZDarcyAnalysis::UniformRefinement(int nh, int MatId)
{
    //    for ( int ref = 0; ref < nh; ref++ ){
    //        TPZVec<TPZGeoEl *> filhos;
    //        int64_t n = fgmesh->NElements();
    //        for ( int64_t i = 0; i < n; i++ ){
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
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
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
    
    std::string plotfile = fSimulationData->GetDirectory();
    plotfile += "/2DMixed.vtk";
    
    scalnames.Push("P");
    vecnames.Push("u");
    vecnames.Push("kappa");    
    scalnames.Push("Porosity");
    scalnames.Push("Rhs");
    scalnames.Push("div_u");
    
//    scalnames.Push("Exact_S");
//    vecnames.Push("Exact_GradS");
    
    
    if (fSimulationData->IsOnePhaseQ()) {
        scalnames.Push("Rho_alpha");
        an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        an->PostProcess(div);
        return;
    }
    
    
    if (fSimulationData->IsTwoPhaseQ()) {
        scalnames.Push("Rho_alpha");
        scalnames.Push("Rho_beta");
        scalnames.Push("S_alpha");
        scalnames.Push("S_beta");
        scalnames.Push("f_alpha");
        scalnames.Push("f_beta");
        scalnames.Push("P_alpha");
        scalnames.Push("P_beta");
        scalnames.Push("Pc_beta_alpha");
        vecnames.Push("u_alpha");
        vecnames.Push("u_beta");
        vecnames.Push("u_alpha_sc");
        vecnames.Push("u_beta_sc");
        an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        an->PostProcess(div);
        return;
    }
    
    
    if (fSimulationData->IsThreePhaseQ()) {
        scalnames.Push("S_alpha");
        scalnames.Push("S_beta");
        scalnames.Push("S_gamma");
        scalnames.Push("Rho_alpha");
        scalnames.Push("Rho_beta");
        scalnames.Push("Rho_gamma");
        an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        an->PostProcess(div);
        return;
    }
    
    
}

void TPZDarcyAnalysis::BCDfunction(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &ff, TPZFMatrix<REAL> &Grad){
    
    REAL t = time;
    REAL pd;
    REAL x = pt[0];
//    REAL y = pt[1];
    
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
//    REAL y = pt[1];
    
    if(time <= 0.0){
        t = 0.0001;
    }
    
    qn = -1.*(1 - 1/(exp(x/t)*t))*(1 + (-0.1 + exp(-x/t) + x)/10.);
    ff[0] = -qn;
}

void TPZDarcyAnalysis::Ffunction(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &ff, TPZFMatrix<REAL> &Grad)
{
    REAL rwD = 1.0*0.127/10.0;
    REAL rD  = pt[0];
    REAL zD  = pt[1];
    REAL a = 1.0;
    REAL b = 0.1;
    
//    REAL f  = ((std::pow(a,2.0) + std::pow(b,2.0))*std::pow(M_PI,2.0)*sin((M_PI*rD)/a)*sin((M_PI*zD)/b))/(std::pow(a,2.0)*std::pow(b,2.));
//    REAL f = (M_PI*(-(a*std::pow(b,2.0)*cos((M_PI*(rD - rwD))/a)) + (std::pow(a,2.0) + std::pow(b,2.0))*M_PI*rD*sin((M_PI*(rD - rwD))/a))*sin((M_PI*zD)/b))/(std::pow(a,2.0)*std::pow(b,2.0)*rD);
//    REAL f = 2.0+ 2.0*(rwD+rD)/rD;
//    REAL f = (5.0*cos(5.0*(rwD+rD))/rD)-25.0*sin(5.0*(rwD+rD));
//    REAL f = -1.0/((rwD+rD)*(rwD+rD))+1.0/((rwD+rD)*rD);
//    ff[0] = f;
    ff[0] = 0.0;
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
//    oldmat.Inverse( * PreInverse);
    DebugStop();
    oldmat.Multiply(*PreInverse, Identity);
    
#ifdef PZDEBUG
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

void TPZDarcyAnalysis::FilterSaturations(TPZManVector<int64_t> &active, TPZManVector<int64_t> &nonactive){
    
    int ncon_sw = fmeshvec[2]->NConnects();
    int ncon = fcmesh->NConnects();
    
    
    // DOF related with the Q-P system
    for(int i = 0; i < ncon-ncon_sw; i++)
    {
        TPZConnect &con = fcmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmesh->Block().Position(seqnum);
        int blocksize = fcmesh->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+blocksize);
        for(int ieq = 0; ieq<blocksize; ieq++)
        {
            active[vs+ieq] = pos+ieq;
        }
    }
    
   
    // DOF Related with S
    for(int i = ncon-ncon_sw; i< ncon; i++)
    {
        TPZConnect &con = fcmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmesh->Block().Position(seqnum);
        int blocksize = fcmesh->Block().Size(seqnum);
        int vs = nonactive.size();
        nonactive.Resize(vs+blocksize);
        for(int ieq = 0; ieq<blocksize; ieq++)
        {
            nonactive[vs+ieq] = pos+ieq;
        }
        
    }
    
    
}

void TPZDarcyAnalysis::FilterSaturationGradients(TPZManVector<int64_t> &active, TPZManVector<int64_t> &nonactive)
{

    int ncon_sw = fmeshvec[2]->NConnects();
    int ncon = fcmesh->NConnects();
    
    
    // DOF related with the Q-P system
    for(int i = 0; i < ncon-ncon_sw; i++)
    {
        TPZConnect &con = fcmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmesh->Block().Position(seqnum);
        int blocksize = fcmesh->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+blocksize);
        for(int ieq = 0; ieq<blocksize; ieq++)
        {
            active[vs+ieq] = pos+ieq;
        }
    }
    
    // DOF Related with S constant
    for(int i = ncon-ncon_sw; i< ncon; i++)
    {
        TPZConnect &con = fcmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmesh->Block().Position(seqnum);
        int blocksize = fcmesh->Block().Size(seqnum);
        int vs = active.size();
        active.Resize(vs+1);
        
        int ieq = blocksize-1;
        active[vs] = pos+ieq;
    }
    
    // DOF Related with S grandients
    for(int i = ncon-ncon_sw; i< ncon; i++)
    {
        TPZConnect &con = fcmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmesh->Block().Position(seqnum);
        int blocksize = fcmesh->Block().Size(seqnum);
        int vs = nonactive.size();
        nonactive.Resize(vs+blocksize-1);
        for(int ieq = 0; ieq<blocksize-1; ieq++)
        {
            nonactive[vs+ieq] = pos+ieq;
        }
        
    }
    
}

void TPZDarcyAnalysis::Cylindrical_Elliptic(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &GradSol){
    
    REAL rwD = 1.0*0.127/100.0;
    REAL rD  = pt[0];
    REAL zD  = pt[1];
    REAL a = 1.0;
    REAL b = 0.1;

    Sol[0] = sin(M_PI* (rD-rwD)/a)*sin(M_PI*zD/b);
    GradSol(0,0) = -((M_PI*cos((M_PI*(rD-rwD))/a)*sin((M_PI*zD)/b))/a);
    GradSol(1,0) = -((M_PI*cos((M_PI*zD)/b)*sin((M_PI*(rD-rwD))/a))/b);
    
}

void TPZDarcyAnalysis::Dupuit_Thiem(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &GradSol){
    
    REAL rstar = 1000.0;
    REAL r = pt[0]*rstar;
    REAL rw = 0.127;
    REAL h = 10.0;
    REAL re = 1000.0 + rw;

    REAL Pstar = 20.0*1.0e6;
    REAL day = 86400.0;
    
    REAL Q = 158.99/day;//158.99/day;
    REAL u = Q/(2.0*M_PI*rw*h);
    REAL rho = 1000.0;
    REAL mu = 0.001;
    REAL k = 1.0e-13;
    REAL muD = mu/mu;
    REAL rhoD = rho/rho;
    REAL kD = k/k;

    REAL m = rho * u ;
    REAL reD = re/rstar;
    REAL rwD = rw/rstar;
    REAL rD  = r/rstar;

    REAL mD = (m*mu*rstar)/(k*rho*Pstar);
    Sol[0] = 1 + mD * rwD * log(rD/reD);
    GradSol(0,0) = -((kD*rhoD*mD*rwD)/muD)/(rD);

}

void TPZDarcyAnalysis::Morris_Muskat(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &GradSol){
    
    REAL rw = 0.127;
    REAL re = 1000.0 + rw;
    REAL reD = re/rw;
    REAL rD = pt[0];
    REAL rhoD = 1.0;
    REAL muD = 1.0;
    REAL kD = 1.0;
    Sol[0] = - log(rD/reD)/log(re);
    GradSol(0,0) =  (kD*rhoD/muD)/(rD*log(re));
    
    double x = 0.2;
    int p = 0;
    //double J = boost::math::cyl_bessel_j(p,x);
    //double Y = boost::math::cyl_neumann(p,x);
    //std::cout << "x = " << x << ", J = " << J << ", Y = " << Y << std::endl;
    int p2 = 0;
}

void TPZDarcyAnalysis::LinearTracer(const TPZVec< REAL >& pt, REAL time, TPZVec< STATE >& Saturation, TPZFMatrix< STATE >& Grad)
{
    REAL x = pt[0];
    REAL v = 1.0;
    REAL Porosity = 0.2;
    REAL Sor            = 0.20;
    REAL Swr            = 0.20;
    REAL rho = 1.0;
    REAL xshock = (v*time/((1.0-Swr-Sor)*Porosity*rho));
    
    if(x < xshock)
    {
        Saturation[0] = 0.8;
    }
    else
    {
        Saturation[0] = 0.2;
    }
    
    return;
    
}


void TPZDarcyAnalysis::BluckleyAndLeverett(const TPZVec< REAL >& pt, REAL time, TPZVec< STATE >& Saturation, TPZFMatrix< STATE >& Grad)
{
    
    REAL Porosity = 0.2;
    REAL x = pt[0];
    
    REAL x_shock        = 0.0;
    REAL S_shock        = 0.0;
    REAL epsilon        = 1.0e-8;
    REAL u              = 0.1;
    REAL mu_alpha       = 1.0;
    REAL mu_beta        = 1.0;
    REAL rho_alpha      = 1.0;
    REAL rho_beta       = 1.0;
    REAL Sor            = 0.2;
    REAL Swr            = 0.2;
    S_shock = Swr + sqrt(mu_alpha*rho_beta*(mu_beta*rho_alpha + mu_alpha*rho_beta)*std::pow(-1.0 + Sor + Swr,2.0))/(mu_beta*rho_alpha + mu_alpha*rho_beta);
    x_shock  = (u*time)/(Porosity*rho_alpha)*dfdsw(S_shock, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta);
  
    if(x < x_shock)
    {
        REAL Sw = S_Newton(x, time, u, Swr, Sor, Porosity, S_shock, mu_alpha, mu_beta, rho_alpha, rho_beta, epsilon);
        Saturation[0] = Sw;
        
    }
    else
    {
        Saturation[0] = Swr;
    }

    return;
    
}


REAL TPZDarcyAnalysis::S_Newton(REAL x, REAL t, REAL u, REAL Swr, REAL Sor, REAL phi, REAL s_shok, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta, REAL epsilon)
{
    REAL S_trial = ((1.0-Sor) + s_shok)/2.0;
    REAL jac    = 0.0;
    REAL r      = 1.0;
    REAL delta_S;
    REAL ds;
    REAL S_k = S_trial;
    int max = 20;
    int it = 0;
    
    
    
    while ( fabs(r) > epsilon && it < max){

        ds = dfdsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta);
        r = x - (u * t * ds)/(phi * rho_alpha);

        jac = - (u * t * df2dsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta))/(phi * rho_alpha);
        
        delta_S = -r/jac;
        S_k += delta_S;
        
        ds = dfdsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta);
        r = x - (u * t * ds)/(phi * rho_alpha);
        it++;
        
    }
    
    return S_k;
}

REAL TPZDarcyAnalysis::dfdsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta)
{
    REAL dfwdSwS;
    
    
    dfwdSwS = (2.0*mu_alpha*mu_beta*rho_alpha*rho_beta*(-1.0 + Sor + Sw)*(Sw - Swr)*(-1.0 + Sor + Swr))/std::pow(mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0) + mu_beta*rho_alpha*std::pow(Sw - Swr,2.0),2.0);
    
    return (dfwdSwS);
}
REAL TPZDarcyAnalysis::df2dsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta)
{
    REAL dfw2dSwS2;

  
    dfw2dSwS2 = (2.0*mu_alpha*mu_beta*rho_alpha*rho_beta*(-1.0 + Sor + Swr)*(-(mu_beta*rho_alpha*std::pow(Sw - Swr,2.0)*(-3.0 + 3.0*Sor + 2.0*Sw + Swr)) + mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0)*(-1.0 + Sor - 2.0*Sw + 3.0*Swr)))/
    std::pow(mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0) + mu_beta*rho_alpha*std::pow(Sw - Swr,2.0),3.0);

    return dfw2dSwS2;
}


// L2 projection

void TPZDarcyAnalysis::SolveProjection(TPZAnalysis *an, TPZCompMesh *Cmesh)
{
    TPZSkylineStructMatrix full(Cmesh);
    an->SetStructuralMatrix(full);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    an->Run();
    
}


TPZCompMesh * TPZDarcyAnalysis::L2ProjectionCmesh(TPZVec<STATE> &solini)
{
    /// criar materiais
    int dim = 2;

    TPZCompMesh * cmesh = new TPZCompMesh(this->fgmesh);
    
    int n_hydraulic_units = fLayers.size();
    
    for (int ilayer = 0; ilayer < n_hydraulic_units; ilayer++) {
        
        int RockId = fLayers[ilayer]->GetMatIDs()[0]; // Volumetric material
        
        // Material medio poroso
        TPZL2Projection *material = new TPZL2Projection(RockId, dim, 1, solini, fSimulationData->Getsorder());
        cmesh->SetDimModel(dim);
        TPZMaterial * mat(material);
        TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialS_alpha);
        material->SetForcingFunction(forcef);
        cmesh->InsertMaterialObject(mat);
    }
    

    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(fSimulationData->Getsorder());
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
    return cmesh;
    
}

void TPZDarcyAnalysis::InitialS_alpha(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
//    REAL x = pt[0];
    REAL y = pt[1];
    REAL S_wett_nc = 0.2;
    REAL S_nwett_ir = 0.0;
    disp[0] = S_wett_nc;

//    bool inside_x = ( x >= 0.3 ) && (x <= 0.7);
//    bool inside_y = ( y >= 0.3 ) && (y <= 0.7);
    
//    if (y<=0.5) {
//         disp[0] = 0.0;
//    }
//
//    if ( y >= 50.0 ) {
//         disp[0] = 1.0;
//    }
    
}

/**
 * Computes computational mesh for L2 projection
 */
void TPZDarcyAnalysis::P_Hydrostatic(const TPZVec< REAL >& pt, REAL time, TPZVec< STATE >& P_Hydro, TPZFMatrix< STATE >& GradP_Hydro){
    
//    REAL x = pt[0];
    REAL y = pt[1];
    
//    REAL Kstr           = 1.0e-13;
    REAL Pstr           = 2.0e7;
//    REAL Tstr           = 355.37;
//    REAL Tres           = 355.37;
    REAL Lstr           = 500.0;
//    REAL Mustr          = 0.001;
    REAL Rhostr         = 800.0;

    REAL rho_beta = 900.0/Rhostr;
    REAL P_at_datum = 0.5240;//2.0*1.0e7;
    REAL g = -10.0*((Lstr*Rhostr)/Pstr);
    P_Hydro[0] = (rho_beta * g * y)+P_at_datum;
    
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
    
    
    int64_t neq = fcmeshdarcy->NEquations();
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
    
    int64_t neq = fcmeshdarcy->NEquations();
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
    
    int64_t neq = fcmeshdarcy->NEquations();
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
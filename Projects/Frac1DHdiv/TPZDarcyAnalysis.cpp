 
//  TPZDarcyAnalysis.cpp
//  PZ
//
//  Created by Nathan Shauer and Omar Duran on 9/8/14.
//
//

#include "pzlog.h"
#include "TPZDarcyAnalysis.h"
#include "TPZMatDarcy2dhdiv.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZReadGIDGrid.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZCompElDisc.h"
#include "pzl2projection.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzcompelwithmem.h"
#include "pzelchdiv.h"

#ifdef USING_BOOST
#include <boost/math/special_functions/erf.hpp>
#endif

#include "TPZFracAnalysis.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.frac"));
#endif

TPZDarcyAnalysis::TPZDarcyAnalysis(TPZAutoPointer<TPZFracData> Data)
{
    fData = Data;
    fmeshvec.Resize(2);
    fgmesh = NULL;
    fcmeshMixed = NULL;
    for (int i = 0; i < 2; i++) {
      fmeshvec[i] = NULL;
    }
    fLastStepRhs.Redim(0, 0);
    fmustStop = false;
}


TPZDarcyAnalysis::~TPZDarcyAnalysis()
{
    fData = NULL;
    fmeshvec.Resize(2);
    fgmesh = NULL;
    fcmeshMixed = NULL;
    for (int i = 0; i < 2; i++) {
        fmeshvec[i] = NULL;
    }
    fLastStepRhs.Redim(0, 0);
    fmustStop = false;
    
}

/** @brief Initial pressure field */
void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    //    REAL x = pt[0];
    //    REAL y = pt[1];
    disp[0] = 15.;
}

/** @brief Analytic pressure field */
void PressureAnalytic(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux)
{
#ifdef USING_BOOST
    REAL x = pt[0], t=time;
    if (time <= 1.0e-8){t=1.0e-8;}
    sol[0]      =   (sqrt((4.0*t)/(M_PI))*exp(-1.0*(x*x)/(4.0*t))) - x*(1.0-boost::math::erf(x/sqrt(4.0*t)));
    flux(0,0)   =   (1.0-boost::math::erf(x/sqrt(4.0*t)));
#endif
}

void TPZDarcyAnalysis::Run()
{
    fData->DebugMap()[fData->Time()] = 0.;
    // Computing vl as the first memory value on each integration point
    REAL vl = this->RunUntilOpen();
    TPZFMatrix<REAL> vlMatrix(1,1,vl);
    
    const REAL lFrac = fData->ElSize();
    fData->SetLfrac(lFrac);
    fData->DebugMap()[fData->Time()] = lFrac;
    
    // Solving initial darcy
  
    // Malha geometrica
//    fgmesh = CreateGMesh(nel);
  
    const int npropag = fData->NPropagations();
    this->CreateGeoMeshQuad(npropag,fData->Hy(),fData->Ly());
    this->PrintGeometricMesh(fgmesh);
    

    //Indcluding the fist ghost element of frac
//    DarcyGmesh(fgmesh);
//    this->PrintGeometricMesh(fgmesh);
//    CreateMultiphysicsMesh(vlMatrix);
  

    this->InsertFracGeoMesh();
    this->PrintGeometricMesh(fgmesh);
    CreateMultiphysicsMesh(vlMatrix);
  
    // Analysis
    bool mustOptimizeBandwidth = false;
    TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
    TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
    skyl.SetNumThreads(fData->NThreadsForAssemble());
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an->SetSolver(step);
    an->SetStructuralMatrix(skyl);
//    SolveSistTransient(an,true);
    if(fData->IsPlotVTK())
    {
        this->PostProcessVTK(an);
    }
    
    while (fmustStop == false)
    {
        bool propagate = SolveSistTransientWithFracture(an);
        const int nfracel = this->HowManyFracElement();
        if (nfracel == npropag - 1){ // So it can not propagate further than domain in x
          propagate = false;
          fmustStop = true;
        }
        if (propagate)
        {
            
          // Novo comprimento de fratura
          REAL newLfrac = fData->Lfrac() + fData->ElSize();
          fData->SetLfrac(newLfrac);
          std::cout << " ------> New Lfrac = " << newLfrac << std::endl;
          fData->DebugMap()[fData->Time()] = newLfrac;
          
          this->InsertFracGeoMesh();
          this->InsertFracCompMesh();
          this->PrintGeometricMesh(fgmesh);
          
          fcmeshMixed->ConnectVec().Resize(0);
          TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshMixed);
          
          // Ajustando a estrutura de dados
          fmeshvec[0]->ComputeNodElCon();
          fmeshvec[1]->ComputeNodElCon();
          fcmeshMixed->ComputeNodElCon();
          fmeshvec[0]->CleanUpUnconnectedNodes();
          fmeshvec[1]->CleanUpUnconnectedNodes();
          fcmeshMixed->CleanUpUnconnectedNodes();
          fmeshvec[0]->ExpandSolution();
          fmeshvec[1]->ExpandSolution();
          fcmeshMixed->ExpandSolution();
          
          this->SetInterfaceConnects(); // seta o vetor de connects dos elementos de interface multifisicos
          
          // Transferindo para a multifisica
          TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);

                
          TPZEquationFilter newEquationFilter(fcmeshMixed->NEquations());
          an->StructMatrix()->EquationFilter() = newEquationFilter; //AQUINATHAN
          
          an->LoadSolution(fcmeshMixed->Solution());
            if (fData->IsPlotVTK()) {
              this->PostProcessVTK(an);
            }
        }
    }
    

    
    std::string filename = "lfracNoCoupled.nb";
    if (fData->IsCoupled()) { filename = "lfracCoupled.nb";}
    fData->PrintDebugMapForMathematica(filename);
    
    delete an;
    
}


bool TPZDarcyAnalysis::SolveSistTransientWithFracture(TPZAnalysis *an)
{
    
    bool propagate = false;
    int nfracel = this->HowManyFracElement();
    int it = 0;
    this->SetPressureOnNewElement(an);
    while (fmustStop == false && propagate == false) {
        AssembleLastStep(an);
        TPZFMatrix<STATE> lastSol = an->Solution();
    
        if (it == 0)
        { // aqui inicializo chutes iniciais para newton depois da propagacao
        
            if(nfracel == 1){ // esse caso eh para o primeiro elemento da simulacao
                this->ComputeFirstSolForOneELement(an);
            }
            else{ // aqui eh quando ha mais de 1 elemento
                this->SetPressureOnLastElement(an);
            }
            
            fData->SetAccumVl(0.); // Zerando accumvl para os proximos
            
        }
      
        IterativeProcess(an, std::cout, 50);
      
        const REAL qtip = this->Qtip();
        std::cout << "\nqtip = " << qtip << std::endl;
        
        if (qtip < 0.) {
            DebugStop();
            propagate = false;
        }
        else{
          propagate = VerifyIfPropagate(qtip);
        }
        
        
        if (propagate) {
            fData->SetLastQtip(qtip);
            an->Solution() = lastSol;
            an->LoadSolution();
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
            std::cout << "\n******************************** FRACTURE PROPAGATED ********************************" << std::endl;
        }
        else{
            fData->SetNextTime();
            AcceptSolution(an); // updates leak off
            if (fData->IsPlotVTK()){
                this->PostProcessVTK(an);
            }
        }
        
        REAL peteleco = 1.e-8;
        if( fData->Time() > (fData->TotalTime() - peteleco) )
        {
            fmustStop = true;
        }
        it++;
    }
    return propagate;
}


TPZGeoMesh * TPZDarcyAnalysis::CreateGMesh(const int nel)
{
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    
    //  Reading mesh
    std::string GridFileName;
    GridFileName = dirname + "/Projects/Frac1DHdiv/";
    GridFileName += "OilWaterSystemUnit.dump";
//    GridFileName += "BaseGeometryDakeThin.dump";//"FiveSpot.dump";
//    GridFileName += "FiveSpot.dump";//"FiveSpot.dump";
 //   REAL angle = 0.0*M_PI/4.0;
    
    TPZReadGIDGrid GeometryInfo;
    GeometryInfo.SetfDimensionlessL(1.0);
    TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    
//    // Inserting
//    TPZGeoEl *  Quad = gmesh->ElementVec()[0];
//    TPZGeoElSide gelside(Quad,0);
//    TPZGeoEl *  InletPoint = Quad->CreateBCGeoEl(0,7);
//    int64_t index;
//    TPZVec<int64_t> TopologyPoint(1,InletPoint->Node(0).Id());
//    gmesh->CreateGeoElement(EPoint, TopologyPoint, 7, index);
//    RotateGeomesh(gmesh, angle);
    
    UniformRefinement(gmesh, nel);

    //  Print Geometrical Base Mesh
    this->PrintGeometricMesh(gmesh);
    
    return gmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshFluxHdiv()
{
    const int matId2d = TPZFracData::EMatDarcy, bcBottomId = TPZFracData::EBCBottom, bcRightId = TPZFracData::EBCRight, bcTopId = TPZFracData::EBCTop, bcLeftId = TPZFracData::EBCLeft, bcBottomIdAux = TPZFracData::EBCAuxBottom;
    const int matId1d = TPZFracData::EMatFrac, MatINterface = TPZFracData::EMatInterFrac;
    const int bcinlet = TPZFracData::EBCInlet, bcoutlet = TPZFracData::EBCOutlet;
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId2d);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
  
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    // Bc Bottom being used for insert Qtip
    TPZBndCond * Auxbcbottom = mat->CreateBC(mat, bcBottomIdAux, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(Auxbcbottom);

    // Interface material
    TPZMatfrac1dhdiv *mat1d = new TPZMatfrac1dhdiv(matId1d);
    cmesh->InsertMaterialObject(mat1d);

    // Material da fratura
    TPZMatfrac1dhdiv *matInterface1d = new TPZMatfrac1dhdiv(MatINterface);
    cmesh->InsertMaterialObject(matInterface1d);
    
    // Condicao de contorno na esquerda
    val2(0,0) = fData->Q();
    TPZBndCond * bcin = mat1d->CreateBC(mat1d, bcinlet, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcin);
    
    // Condicao de contorno na direita
    val2(0,0) = fData->SigmaConf();
    TPZBndCond * bcout = mat1d->CreateBC(mat1d, bcoutlet, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcout);
    
    
    // Setando Hdiv
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(fData->PorderDarcyFlow());
    cmesh->SetAllCreateFunctionsHDiv();
    
    std::set<int> BuildGroup2D,BuildGroup1D;
    BuildGroup2D.insert(matId2d);
    BuildGroup2D.insert(bcLeftId);
    BuildGroup2D.insert(bcRightId);
    BuildGroup2D.insert(bcTopId);
    BuildGroup2D.insert(bcBottomId);
    BuildGroup2D.insert(bcBottomIdAux);
    cmesh->AutoBuild(BuildGroup2D);

    // Setando H1
    fgmesh->ResetReference();
    cmesh->SetDimModel(1);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(fData->PorderFlow());
    BuildGroup1D.insert(matId1d);
    BuildGroup1D.insert(bcinlet);
    BuildGroup1D.insert(bcoutlet);
    cmesh->AutoBuild(BuildGroup1D);
 
    cmesh->SetDimModel(2);
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshPressureL2()
{

    const int matId2d = TPZFracData::EMatDarcy, bcBottomId = TPZFracData::EBCBottom, bcRightId = TPZFracData::EBCRight, bcTopId = TPZFracData::EBCTop, bcLeftId = TPZFracData::EBCLeft, bcBottomIdAux = TPZFracData::EBCAuxBottom;
    const int matId1d = TPZFracData::EMatFrac, MatINterface = TPZFracData::EMatInterFrac;
    const int bcinlet = TPZFracData::EBCInlet, bcoutlet = TPZFracData::EBCOutlet;
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material da fratura
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId2d);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
    
    // Bc Bottom being used for insert Qtip
    TPZBndCond * Auxbcbottom = mat->CreateBC(mat, bcBottomIdAux, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(Auxbcbottom);
    
    // Material da fratura
    TPZMatfrac1dhdiv *matInterface1d = new TPZMatfrac1dhdiv(MatINterface);
    cmesh->InsertMaterialObject(matInterface1d);
    
    // Material da fratura
    TPZMatfrac1dhdiv *mat1d = new TPZMatfrac1dhdiv(matId1d);
    cmesh->InsertMaterialObject(mat1d);
    
    // Condicao de contorno na esquerda
    val2(0,0) = fData->Q();
    TPZBndCond * bcin = mat1d->CreateBC(mat1d, bcinlet, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcin);
    
    // Condicao de contorno na direita
    val2(0,0) = fData->SigmaConf();
    TPZBndCond * bcout = mat1d->CreateBC(mat1d, bcoutlet, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcout);
    
    
    // Setando L2
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(fData->PorderDarcyPressure());
    cmesh->SetAllCreateFunctionsDiscontinuous();

    std::set<int> BuildGroup2D,BuildGroup1D;
    BuildGroup2D.insert(matId2d);
    BuildGroup2D.insert(bcLeftId);
    BuildGroup2D.insert(bcRightId);
    BuildGroup2D.insert(bcTopId);
    BuildGroup2D.insert(bcBottomId);
    BuildGroup2D.insert(bcBottomIdAux);
    cmesh->AutoBuild(BuildGroup2D);
    
    fgmesh->ResetReference();
    cmesh->SetDimModel(1);
    cmesh->SetDefaultOrder(fData->PorderPressure());
    BuildGroup1D.insert(matId1d);
    BuildGroup1D.insert(bcinlet);
    BuildGroup1D.insert(bcoutlet);
    cmesh->AutoBuild(BuildGroup1D);
    
    cmesh->SetDimModel(2);
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
#ifdef PZDEBUG
  std::ofstream out("cmeshPress.txt");
  cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshMixed(TPZFMatrix<REAL> Vl)
{
    // Definicao de ids e tipos
    const int matId2d = TPZFracData::EMatDarcy, bcBottomId = TPZFracData::EBCBottom, bcRightId = TPZFracData::EBCRight, bcTopId = TPZFracData::EBCTop, bcLeftId = TPZFracData::EBCLeft, bcBottomIdAux = TPZFracData::EBCAuxBottom;
    const int matId1d = TPZFracData::EMatFrac, MatINterface = TPZFracData::EMatInterFrac;
    const int bcinlet = TPZFracData::EBCInlet, bcoutlet = TPZFracData::EBCOutlet;
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material medio poroso
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId2d);
    mat->SetSimulationData(fData);
    cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > TimeDepFExact = new TPZDummyFunction<STATE>(PressureAnalytic, 5);
    mat->SetTimeDependentFunctionExact(TimeDepFExact);
    
    // Bc Bottom
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = fData->Pe();
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    val2(0,0) = 0.0;// Massic flux 5.0 kg/s over 100000 m2
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
    
    // Bc Bottom being used for insert Qtip
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * Auxbcbottom = mat->CreateBC(mat, bcBottomIdAux, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(Auxbcbottom);
    
    // Material da fratura
    TPZMatfrac1dhdiv *matInterface1d = new TPZMatfrac1dhdiv(MatINterface);
    matInterface1d->SetSimulationData(fData);
    matInterface1d->SetDefaultMem(Vl);
    cmesh->InsertMaterialObject(matInterface1d);
    
    // Material da fratura
    TPZMatfrac1dhdiv *mat1d = new TPZMatfrac1dhdiv(matId1d);
    mat1d->SetSimulationData(fData);
    cmesh->InsertMaterialObject(mat1d);
    mat1d->SetDefaultMem(Vl);
    
    // Condicao de contorno na esquerda
    val2(0,0) = fData->Q();
    TPZBndCond * bcin = mat1d->CreateBC(mat1d, bcinlet, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcin);
    
    // Condicao de contorno na direita
    val2(0,0) = fData->SigmaConf();
    TPZBndCond * bcout = mat1d->CreateBC(mat1d, bcoutlet, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcout);
    
    // Setando Multifisico
    cmesh->SetDimModel(2);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    std::set<int> BuildGroup2D,BuildGroup1D;
    BuildGroup2D.insert(matId2d);
    BuildGroup2D.insert(bcLeftId);
    BuildGroup2D.insert(bcRightId);
    BuildGroup2D.insert(bcTopId);
    BuildGroup2D.insert(bcBottomId);
    BuildGroup2D.insert(bcBottomIdAux);
    cmesh->AutoBuild(BuildGroup2D);
    
    fgmesh->ResetReference();
    cmesh->SetDimModel(1);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    BuildGroup1D.insert(matId1d);
    BuildGroup1D.insert(bcinlet);
    BuildGroup1D.insert(bcoutlet);
    cmesh->AutoBuild(BuildGroup1D);
    
    cmesh->SetDimModel(2);
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
    int dim = 2;
    TPZL2Projection *material;
    material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialPressure, 5);
  
    material->SetForcingFunction(forcef);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
//    ///set order total da shape HERE when Triangles are used
//    int nel = cmesh->NElements();
//    for(int i=0; i<nel; i++){
//        TPZCompEl *cel = cmesh->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//        {
//            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
//            else celdisc->SetTensorialShape();
//        }
//    }
    
    return cmesh;
    
}


void TPZDarcyAnalysis::CreateInterfaces()
{
    fgmesh->ResetReference();
   
    // Creation of interface elements
    int nel = fcmeshMixed->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = fcmeshMixed->ElementVec()[el];
        if(!compEl) continue;
        int index = compEl ->Index();
        if(compEl->Dimension() == fcmeshMixed->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fcmeshMixed->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }

    }
}

void TPZDarcyAnalysis::IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter)
{
  int iter = 0;
  REAL error = 1.e10, NormResLambdaLast = 1.e10;;
  const REAL tol = 5.e-10; // because the SI unit system
  
  fData->SetCurrentState();
  int numeq = an->Mesh()->NEquations();
  
  TPZFMatrix<STATE> Uatk0(an->Solution());
  TPZFMatrix<STATE> Uatk(Uatk0),DeltaU(Uatk0);
  if(Uatk0.Rows() != numeq) Uatk0.Redim(numeq,1);
  
  an->Assemble();
  an->Rhs() += fLastStepRhs;
  an->Rhs() *= -1.0; //- [R(U0)];
  
  TPZAutoPointer< TPZMatrix<STATE> > matK; // getting X(Uatn)
  
  bool notconverged = true;
  while(notconverged && iter < numiter) {
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
      std::stringstream sout;
      matK=an->Solver().Matrix();
      matK->Print("matK = ", sout,EMathematicaInput);
      an->Rhs().Print("Rhs = ", sout, EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    // Computing Uatk = Uatn + DeltaU;
    an->Solve();
    DeltaU= an->Solution();
    Uatk = Uatk0 + DeltaU;
    
    //Computing ||DeltaU||
    REAL NormOfDeltaU = Norm(DeltaU);
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
      std::stringstream sout;
      DeltaU.Print("DeltaU = ", sout,EMathematicaInput);
      Uatk.Print("Uatk = ", sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    // Loading Sol
    an->LoadSolution(Uatk); // Loading Uatk
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
    
    // Putting last Qtip as bc to darcy mesh
    TPZMaterial *mat = fcmeshMixed->FindMaterial(TPZFracData::EBCAuxBottom);
    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
    bnd->Val2()(1,0) = this->Qtip()/fData->ElSize()/2.0;
    //std::cout << "qtip = " << this->Qtip() << std::endl;
    
    // Assembling
    an->Assemble();
    an->Rhs() += fLastStepRhs;
    an->Rhs() *= -1.0; //- [R(U0)];
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
      std::stringstream sout;
      an->Rhs().Print("Res = ", sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    // Computing ||[R(Uatk)]||
    double ResidualNorm = Norm(an->Rhs());
    double norm = NormOfDeltaU; //ResidualNorm;
    out << "Iteration n : " << (iter+1) << " : norms ||DeltaU|| e ||[R(Uatk)]|| : " << NormOfDeltaU << " / " << ResidualNorm << std::endl;
    
    if(norm < tol) {
      out << "\nNewton Converged! Tolerance Of Norm(DeltaU) at n : " << (iter+1) << std::endl;
      out << "Norm ||DeltaU|| - USED : " << NormOfDeltaU << std::endl;
      out << "Norm ||[R(Uatk)]||  : " << ResidualNorm << std::endl;
      out << "--------------------  --------------------" << std::endl;
      notconverged = false;
    }

    else if( (ResidualNorm - NormResLambdaLast) > 1.e-4 ) {
      out << "\nDivergent Method\n" << "You can try implementing Line Search" << std::endl;
      out << " ***** BE AWARE ***** - It may be also a problem due to delayed Qtip. It will most likely fix itself on next iterations" << std::endl;
    }
    
    NormResLambdaLast = ResidualNorm;
    error = norm;
    iter++;
    Uatk0 = Uatk;
    out.flush();
  }
  
  if (error > tol) {
    DebugStop(); // Something is very wrong (spooky :O)
  }
  
}

void TPZDarcyAnalysis::PrintGeometricMesh(TPZGeoMesh * gmesh)
{
    std::ofstream argument("GeometicMesh.txt");
    gmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
}

void TPZDarcyAnalysis::PrintComputatiolaMeshInfo(TPZCompMesh * cmesh)
{
  std::ofstream argument("ComputationalMesh.txt");
  cmesh->Print(argument);
}

void TPZDarcyAnalysis::PostProcessVTK(TPZAnalysis *an)
{
    const int dim = 2;
    int div = fData->HrefPostPro();
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile = "2DMixedDarcy.vtk";
    scalnames.Push("Pressure");
    scalnames.Push("PressureExact");
    vecnames.Push("MassVelocity");
    vecnames.Push("MassVelocityExact");
    an->DefineGraphMesh(dim, scalnames, vecnames, fData->PostProcessFileName());
    an->PostProcess(div);
}


void TPZDarcyAnalysis::AssembleLastStep(TPZAnalysis *an)
{
    fData->SetLastState();
    an->Assemble();
    fLastStepRhs = an->Rhs();
}

void TPZDarcyAnalysis::SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
    
    TPZSkylineStructMatrix skymat(Cmesh);
    an.SetStructuralMatrix(skymat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Run();

}

bool TPZDarcyAnalysis::SolveSistTransient(TPZAnalysis *an, bool initial)
{

    this->PostProcessVTK(an);

    while (fmustStop == false) {
        
        AssembleLastStep(an);
        IterativeProcess(an, std::cout);
        fData->SetNextTime();
        this->PostProcessVTK(an);
        
        if (initial) {
            fmustStop=true;
        }
        
        REAL tricky = 1.E-8;
        if( fData->Time() > (fData->TotalTime() - tricky) )
        {
            fmustStop = true;
        }
    }
  
    return true; // Old SolveSist just for darcy mesh
}



void TPZDarcyAnalysis::UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

void TPZDarcyAnalysis::RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle)
{
    REAL theta = CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    RotationMatrix(0,0) =   +cos(theta);
    RotationMatrix(0,1) =   -sin(theta);
    RotationMatrix(1,0) =   +sin(theta);
    RotationMatrix(1,1) =   +cos(theta);
    RotationMatrix(2,2) = 1.0;
    TPZVec<REAL> iCoords(3,0.0);
    TPZVec<REAL> iCoordsRotated(3,0.0);
    
    RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

void TPZDarcyAnalysis::InsertFracGeoMesh()
{
    int nelements = fgmesh->NElements();
    fgmesh->ResetReference();
    for (int iel = 0 ; iel < nelements; iel++)
    {
        TPZGeoEl *igel = fgmesh->ElementVec()[iel];
        if (igel->HasSubElement()) {
            continue;
        }
      if(igel->MaterialId()== TPZFracData::EBCBottom || igel->MaterialId()== TPZFracData::EBCAuxBottom)
        {
            TPZGeoElSide igelside(igel,0);
            TPZGeoElSide neigh = igelside.Neighbour();
            
            while (igelside != neigh) {
                if (neigh.Element()->Dimension() == 0 || (neigh.Element()->MaterialId()== TPZFracData::EBCLeft && neigh.Element()->Dimension() == 1)) {
                    break;
                }
                neigh = neigh.Neighbour();
            }
            if (igelside == neigh) {
                continue;
            }
            
            TPZGeoEl *Neigel = neigh.Element();
            if (Neigel->MaterialId() == TPZFracData::EBCLeft) {
                igel->SetMaterialId(TPZFracData::EMatFrac);
                int64_t rightnode=igel->NodeIndex(1);
                TPZVec<int64_t> topopoint(1,rightnode);
                int64_t OutLetindex;
                fgmesh->CreateGeoElement(EPoint, topopoint, TPZFracData::EBCOutlet, OutLetindex);
                
                int64_t leftnode=igel->NodeIndex(0);
                topopoint=leftnode;
                fgmesh->CreateGeoElement(EPoint, topopoint, TPZFracData::EBCInlet, OutLetindex);
                igel->CreateBCGeoEl(2, TPZFracData::EMatInterFrac);
                fgmesh->AddInterfaceMaterial(TPZFracData::EMatDarcy, TPZFracData::EMatFrac, TPZFracData::EMatInterFrac);
                fgmesh->AddInterfaceMaterial(TPZFracData::EMatFrac, TPZFracData::EMatDarcy, TPZFracData::EMatInterFrac);
                this->SwitchBCInFrontOfFrac(igel);
                break;
            }else if (Neigel->MaterialId()==TPZFracData::EBCOutlet)
            {
                igel->SetMaterialId(TPZFracData::EMatFrac);
                igel->CreateBCGeoEl(2, TPZFracData::EMatInterFrac);
//                int64_t rightnode=igel->NodeIndex(1);
              
                fcmeshMixed->LoadReferences();
                TPZCompEl *celneighel = Neigel->Reference();
                fgmesh->DeleteElement(Neigel);
                //Neigel->SetMaterialId(1000);
                TPZGeoEl *newbcgel = igel->CreateBCGeoEl(1, TPZFracData::EBCOutlet);
              
                TPZCompEl * cel = igel->Reference();
                SwitchTipElement(cel, celneighel, newbcgel);
                this->SwitchBCInFrontOfFrac(igel);

                break;
            }
            
        }
    }
  
    fgmesh->BuildConnectivity();
}

void TPZDarcyAnalysis::InsertFracCompMesh()
{
  TPZCompMesh *cmesh = this->fcmeshMixed;
  cmesh->LoadReferences();
  const int64_t nel = fgmesh->NElements();
  cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
  cmesh->SetDefaultOrder(fData->PorderPressure());
  // Creation of interface elements for frac coupling by searching interfaces
  
  for(int el = 0; el < nel; el++)
  {
    TPZGeoEl *gel = fgmesh->ElementVec()[el];
    if(!gel) continue;
    
    if(gel->MaterialId() == TPZFracData::EMatInterFrac)
    {
      TPZGeoElSide rightside(gel,1);
      TPZGeoElSide rightneigh = rightside.Neighbour();
      while (rightside != rightneigh) {
        if (rightneigh.Element()->MaterialId() == TPZFracData::EBCOutlet) {
          break;
        }
        rightneigh = rightneigh.Neighbour();
      }
      if (rightside == rightneigh) {
        continue;
      }
      
      int side = gel->NSides()-1;
      TPZStack < TPZCompElSide > neigh;
      TPZGeoElSide gelside(gel,side);
      
      const int onlyinterpolated = 0, removeduplicates = 0;
      gelside.EqualLevelCompElementList(neigh, onlyinterpolated, removeduplicates);
      TPZCompElSide darcyside, fracside;
      for (int i = 0; i < neigh.size(); i++) {
        if (neigh[i].Element()->Reference()->MaterialId() == TPZFracData::EMatDarcy) {
          darcyside = neigh[i];
        }
        if (neigh[i].Element()->Reference()->MaterialId() == TPZFracData::EMatFrac) {
          fracside = neigh[i];
        }
      }
      if (!darcyside || !fracside) {
        DebugStop();
      }
      
      if (neigh.size() == 0) {
        DebugStop();
      }
      int64_t gelindex;
      
      // Preparando os index dos pontos de integracao.
      TPZFMatrix<> Vl(1,1,fData->AccumVl());
      TPZMatfrac1dhdiv *matfrac = dynamic_cast<TPZMatfrac1dhdiv *> (fcmeshMixed->FindMaterial(TPZFracData::EMatInterFrac));
      matfrac->SetDefaultMem(Vl);
      new TPZCompElWithMem <TPZMultiphysicsInterfaceElement >(*cmesh, gel, gelindex, darcyside, fracside);
      break;
    }
  }
  
// this->PrintComputatiolaMeshInfo(fcmeshMixed);
  
    // Colocando a dependencia entre a ponta da fratura e o elemento de darcy a frente
  if(0){
    fgmesh->ResetReference();
    fcmeshMixed->LoadReferences();
    for(int el = 0; el < nel; el++)
    {
        TPZGeoEl *gel = fgmesh->ElementVec()[el];
        if(!gel) continue;
        if(gel->MaterialId() == TPZFracData::EMatDarcy)
        {
            TPZGeoElSide rightside(gel,1);
            TPZGeoElSide rightneigh = rightside.Neighbour();
            while (rightside != rightneigh) {
                if (rightneigh.Element()->MaterialId() == TPZFracData::EBCOutlet) {
                    break;
                }
                rightneigh = rightneigh.Neighbour();
            }
            if (rightside == rightneigh) {
                continue;
            }
            TPZCompEl *cel = gel->Reference();
            TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mcel) {
                DebugStop();
            }
            if (mcel->Element(0)->NConnects() != 3) {
                DebugStop(); // Should be H1 flux element
            }
            int64_t ctipindex = mcel->Element(0)->ConnectIndex(1);
            TPZConnect &ctip = mcel->Element(0)->Connect(1);
            
            rightneigh = rightside.Neighbour();
            TPZMultiphysicsElement *mcelhdiv = NULL;
            while (rightside != rightneigh) {
                if (rightneigh.Element()->Reference()) {
                    rightneigh.Element()->Reference()->Print();
                }
                if (rightneigh.Element()->MaterialId() == TPZFracData::EBCBottom && rightneigh.Element()->Dimension() == 1) {
                    mcelhdiv = dynamic_cast<TPZMultiphysicsElement *>(rightneigh.Element()->Reference());
                    if (mcelhdiv) {
                        break;
                    }
                }
                rightneigh = rightneigh.Neighbour();
            }
            if (rightside == rightneigh) {
                DebugStop();
            }
            
            if (mcelhdiv->NConnects() != 1) {
                DebugStop(); // should be hdiv bc
            }
            
            int64_t cdarcyindex = mcelhdiv->ConnectIndex(0);
            TPZConnect &cdarcy = mcelhdiv->Connect(0);
            
            if (cdarcy.NShape() != 2 || ctip.NShape() != 1) {
                DebugStop();
            }
            REAL area = mcelhdiv->Reference()->SideArea(2);
            TPZFNMatrix<2,REAL> depend(2,1,0.5/area);
            cdarcy.AddDependency(cdarcyindex, ctipindex, depend, 0, 0, 2, 1);
            
            
        }
    }
  }

}

void TPZDarcyAnalysis::DarcyGmesh(TPZGeoMesh * gmesh)
{
    int nelements = fgmesh->NElements();
    for (int iel = 0 ; iel < nelements; iel++)
    {
        TPZGeoEl *igel = fgmesh->ElementVec()[iel];
        if (igel->HasSubElement()) {
            continue;
        }
        if(igel->MaterialId()== TPZFracData::EBCBottom)
        {
            TPZGeoElSide igelside(igel,0);
            TPZGeoElSide neigh = igelside.Neighbour();
            
            while (igelside != neigh) {
                if (neigh.Element()->MaterialId()==TPZFracData::EBCLeft && neigh.Element()->Dimension() == 1) {
                    break;
                }
                neigh = neigh.Neighbour();
            }
            if (igelside == neigh) {
                continue;
            }
          igel->SetMaterialId(TPZFracData::EBCAuxBottom); // Not being used
        }
    }
    
}

// Computing the flux on the right tip of the fracture
REAL TPZDarcyAnalysis::Qtip()
{
    fgmesh->ResetReference();
    fcmeshMixed->LoadReferences();
    const int bcOfTip = TPZFracData::EBCOutlet;
    const int64_t nel = fcmeshMixed->NElements();
    TPZCompEl *cel = NULL;
    for (int64_t iel = 0; iel < nel; iel++) {
        cel = fcmeshMixed->Element(iel);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != 1) {
            continue;
        }
        TPZGeoElSide gelside(gel,1);
        TPZGeoElSide neigh = gelside.Neighbour();
        while (neigh != gelside) {
            if (neigh.Element()->MaterialId() == bcOfTip) {
                break;
            }
            neigh = neigh.Neighbour();
        }
        if (neigh == gelside) {
            continue;
        }
        break;
    }
    
    //Here mat 6 means fracture elements
    TPZMaterial *Material = fcmeshMixed->FindMaterial(TPZFracData::EMatFrac);
    TPZMatfrac1dhdiv * MaterialOfFract = dynamic_cast<TPZMatfrac1dhdiv *>(Material);

    
    TPZVec<REAL> qsi(3,1.);
    TPZManVector<STATE,2> sol(1,0.);
    const int varQ = MaterialOfFract->VariableIndex("Flow");
    cel->Solution(qsi, varQ, sol);
    const REAL qTip = sol[0];
    
    return qTip;
}

REAL TPZDarcyAnalysis::PropagationFlowCriteria(REAL qFreshNewEl, REAL ql){
  return MAX(qFreshNewEl/3., 4.*ql);
}

bool TPZDarcyAnalysis::VerifyIfPropagate(REAL qtip)
{
  const REAL dt = fData->TimeStep();
  //const REAL AccumVolThroughTip = fData->AccumVl() * fData->ElSize() * 2.;
  //const REAL volThroughTip = qtip * dt + AccumVolThroughTip;
  const REAL pfrac = fData->SigmaConf();
  const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
  REAL vl = fData->VlFtau(pfrac, tstar+dt);
  const REAL totalLeakOff = 2. * fData->ElSize() * vl;
  const REAL totalLeakOffprev = 2. * fData->ElSize() * fData->AccumVl();
  const REAL ql = (totalLeakOff - totalLeakOffprev)/dt;
  
  const REAL qFreshNewEl = this->QOfAFreshNewElement();
  const REAL crit = this->PropagationFlowCriteria(qFreshNewEl,ql);
  if (qtip > crit) { // AQUINATHAN
    return true;
  }
  else{
    vl = fData->AccumVl() + qtip*dt/fData->ElSize()/2.;
    fData->SetAccumVl(vl);
    return false;
  }
}

REAL  TPZDarcyAnalysis::RunUntilOpen()
{
  const int maxinitialit = 10000;
  const REAL qtip = fData->Q();
  int it = 0;
  for (it = 0; it < maxinitialit; it++) {
    const REAL dt = fData->TimeStep();
    //const REAL AccumVolThroughTip = fData->AccumVl() * fData->ElSize() * 2.;
    //const REAL volThroughTip = qtip * dt + AccumVolThroughTip;
    const REAL pfrac = fData->SigmaConf();
    const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
    REAL vlnext = fData->VlFtau(pfrac, tstar+dt);
    const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
    const REAL totalLeakOffPrev = 2. * fData->ElSize() * fData->AccumVl();
    const REAL ql = (totalLeakOff - totalLeakOffPrev)/dt;
    
    const REAL qFreshNewEl = this->QOfAFreshNewElement();
    const REAL crit = this->PropagationFlowCriteria(qFreshNewEl,ql);
    if (qtip > crit) { //AQUINATHAN
      break;
    }
    
    vlnext = fData->AccumVl() + qtip*dt/fData->ElSize()/2.;
    fData->SetAccumVl(vlnext);
    fData->SetNextTime();
  }
  if (it == maxinitialit) {
    DebugStop();
  }
  
  std::cout << "#################### Opening of the fracture occured at time t = " << fData->Time() << " s ####################" << std::endl;
  std::cout << "Total vol injected = " << qtip *fData->Time() << std::endl;
  std::cout << "\nStarting First Simulation" << std::endl;
  
  return fData->AccumVl();
}

REAL TPZDarcyAnalysis::QOfAFreshNewElement()
{
  const REAL dt = fData->TimeStep();
  const REAL pfrac = fData->SigmaConf();
  REAL vlnext = fData->VlFtau(pfrac, dt);
  const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
  const REAL qFresh = (totalLeakOff)/dt;
  
  return qFresh;
}


void TPZDarcyAnalysis::SetPressureOnLastElement(TPZAnalysis *an)
{
    fgmesh->ResetReference();
    fmeshvec[1]->LoadReferences();
    
    int64_t nfracel = this->HowManyFracElement();
    int nel = fmeshvec[1]->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl * cel = fmeshvec[1]->Element(iel);
        if (!cel) continue;
        TPZGeoEl * gel = cel->Reference();
        
      if (gel->Dimension() != 1 || gel->MaterialId() != TPZFracData::EMatFrac) {
            continue;
        }
        
        int side = 1;
        TPZGeoElSide gelside(gel,side);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        // Seeking for condition 1
        while (gelside!= neigh) {
            if (neigh.Element()->Dimension()==0 && neigh.Element()->MaterialId() == TPZFracData::EBCOutlet) {
                break;
            }
            neigh = neigh.Neighbour();
        }
        if (gelside == neigh) {
            continue;
        }
        
        TPZBlock<STATE> & block = fmeshvec[1]->Block();
        TPZGeoElSide gelsideleft(gel,0);
        TPZGeoElSide neighTip = gelsideleft.Neighbour();
        
        // Seeking for condition 2
        while (gelsideleft != neighTip) {
            if (neighTip.Element()->Dimension() == 1 && neighTip.Element()->MaterialId() == TPZFracData::EMatFrac ) {
                break;
            }
            neighTip = neighTip.Neighbour();
        }
        if (gelsideleft == neighTip) {
            DebugStop();
        }
        TPZGeoEl * gelleft = neighTip.Element();
        TPZCompEl * celleft = gelleft->Reference();
        
        
        // Chanching value
#ifdef PZDEBUG
        if (celleft->NConnects() != 1 && cel->NConnects() != 1) {
            DebugStop();
        }
#endif
        TPZConnect &connectleft =  celleft->Connect(0);
        TPZConnect &connect =  cel->Connect(0);
        
        int seqleft = connectleft.SequenceNumber();
        int seq = connect.SequenceNumber();
#ifdef PZDEBUG
        if (block.Size(seqleft) != 1 && block.Size(seq) != 1) {
            DebugStop();
        }
#endif
        int posleft = block.Position(seqleft);
        int pos = block.Position(seq);
        fmeshvec[1]->Solution()(pos,0) = fmeshvec[1]->Solution()(posleft,0);
        
        if (nfracel < 5) {
            TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
            an->LoadSolution(fcmeshMixed->Solution());
            return;
        }
        
        // Seeking for condition 3
        TPZGeoElSide gelsidesecondleft(gelleft,0);
        TPZGeoElSide neighsecondeleft = gelsidesecondleft.Neighbour();
        while (gelsidesecondleft != neighsecondeleft) {
            if (neighsecondeleft.Element()->Dimension() == 1 && neighsecondeleft.Element()->MaterialId() == TPZFracData::EMatFrac) {
                break;
            }
            neighsecondeleft = neighsecondeleft.Neighbour();
        }
        if (gelsidesecondleft == neighsecondeleft) {
            DebugStop();
        }
        TPZGeoEl * gelsecondleft = neighsecondeleft.Element();
        TPZCompEl * celsecondleft = gelsecondleft->Reference();
#ifdef PZDEBUG
        if (celsecondleft->NConnects() != 1) {
            DebugStop();
        }
#endif
        TPZConnect &connectsecondleft = celsecondleft->Connect(0);
        int seqsecondleft = connectsecondleft.SequenceNumber();
        
#ifdef PZDEBUG
        if (block.Size(seqsecondleft) != 1) {
            DebugStop();
        }
#endif
        int possecondleft = block.Position(seqsecondleft);
        
        fmeshvec[1]->Solution()(posleft,0) = fmeshvec[1]->Solution()(possecondleft,0);
        
        TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
        an->LoadSolution(fcmeshMixed->Solution());
        return;
    }
}

int TPZDarcyAnalysis::HowManyFracElement()
{
    int64_t nel = fgmesh->NElements();
    int64_t nfracel = 0;
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = fgmesh->Element(iel);
        if (gel->MaterialId() == TPZFracData::EMatFrac) {
            nfracel++;
        }
    }
    return nfracel;
}

void TPZDarcyAnalysis::CreateMultiphysicsMesh(TPZFMatrix<REAL> Vl)
{
    fmeshvec[0] = CreateCMeshFluxHdiv();
    fmeshvec[1] = CreateCMeshPressureL2();
    
    // Initial Pressure
    TPZVec<STATE> solini(1,0.0);
    TPZCompMesh  * cmeshL2 = L2ProjectionP(fgmesh, fData->PorderDarcyPressure(), solini);
    TPZAnalysis anL2(cmeshL2,0);
    SolveSyst(anL2, cmeshL2);
    fmeshvec[1]->LoadSolution(anL2.Solution());
  
    fcmeshMixed = CreateCMeshMixed(Vl);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshMixed);
    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshMixed);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
    
    // Create Interfaces
    this->CreateInterfaces();
    this->InsertFracCompMesh();

    this->SetIntPointMemory();
  
    std::ofstream dumpfile("ComputationaMeshMultiphysics.txt");
    fcmeshMixed->Print(dumpfile);
    
}

void TPZDarcyAnalysis::SetIntPointMemory()
{
  // Preparando os index dos pontos de integracao.
  int64_t nel = fcmeshMixed->NElements();
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZCompEl *cel = fcmeshMixed->ElementVec()[iel];
    cel->PrepareIntPtIndices();
  }
  
}

void TPZDarcyAnalysis::AcceptSolution(TPZAnalysis *an)
{
    //Here mat 6 means fracture elements
    TPZMaterial *material = fcmeshMixed->FindMaterial(TPZFracData::EMatInterFrac);
    TPZMatfrac1dhdiv * materialOfFract = dynamic_cast<TPZMatfrac1dhdiv *>(material);
    material = fcmeshMixed->FindMaterial(TPZFracData::EMatDarcy);
    TPZMatDarcy2dhdiv * materialOfDarcy = dynamic_cast<TPZMatDarcy2dhdiv *>(material);
    materialOfFract->SetUpdateMem();
    materialOfDarcy->SetNotContribute();
    an->AssembleResidual();
    materialOfFract->SetUpdateMem(false);
    materialOfDarcy->SetNotContribute(false);
}

void TPZDarcyAnalysis::ComputeFirstSolForOneELement(TPZAnalysis * an)
{
    int64_t nfracel = this->HowManyFracElement();
    
    if (nfracel != 1) {
        PZError << "This method sould only be called when the mesh has a single frac element " << std::endl;
        DebugStop();
    }
    
    int nel = fgmesh->NElements();
    TPZGeoEl *gel = NULL;
    for (int iel = 0; iel < nel; iel++) {
        gel = fgmesh->Element(iel);
        if (gel->MaterialId() == TPZFracData::EMatFrac) {
            break;
        }
    }
    
    TPZBlock<STATE> &blockQ = fmeshvec[0]->Block(), &blockP = fmeshvec[1]->Block();
    
    
    // Setando os valores dos fluxos na fratura
    fgmesh->ResetReference();
    fmeshvec[0]->LoadReferences();
    
    const REAL pfrac = fData->SigmaConf();
    const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
    const REAL vlnext = fData->VlFtau(pfrac, tstar+fData->TimeStep());
    const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
    const REAL totalLeakOffPrev = 2. * fData->ElSize() * fData->AccumVl();
    const REAL ql = (totalLeakOff - totalLeakOffPrev)/fData->TimeStep();
    REAL qout = fData->Q() - ql;
    
    TPZCompEl *celQ = gel->Reference();
    if (celQ->NConnects() != 3) {
        DebugStop(); // Mesh H1 1D p = 1
    }
    TPZConnect &c1Q = celQ->Connect(0), &c2Q = celQ->Connect(1);
    int seq1Q = c1Q.SequenceNumber(), seq2Q = c2Q.SequenceNumber();
    int pos1Q = blockQ.Position(seq1Q), pos2Q = blockQ.Position(seq2Q);
    fmeshvec[0]->Solution()(pos1Q,0) = fData->Q();
    fmeshvec[0]->Solution()(pos2Q,0) = qout;
    
    
    // Setando a pressao
    fgmesh->ResetReference();
    fmeshvec[1]->LoadReferences();
    const REAL dwdp = fData->GetDwDp();
    const REAL pini = fData->SigmaConf() + pow(12. * fData->Viscosity() * qout * fData->ElSize() / (dwdp*dwdp*dwdp),1./4.);
    
    TPZCompEl *celP = gel->Reference();
    if (celP->NConnects() != 1) {
        DebugStop(); // Mesh L2 1D p = 0
    }
    TPZConnect &cP = celP->Connect(0);
    int seqP = cP.SequenceNumber();
    int posP = blockP.Position(seqP);
    
    fmeshvec[1]->Solution()(posP,0) = pini;
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
    an->LoadSolution(fcmeshMixed->Solution());
    
}

void TPZDarcyAnalysis::SetPressureOnNewElement(TPZAnalysis *an)
{
    fgmesh->ResetReference();
    fmeshvec[1]->LoadReferences();
    
    int nel = fmeshvec[1]->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl * cel = fmeshvec[1]->Element(iel);
        if (!cel) continue;
        TPZGeoEl * gel = cel->Reference();
        
        if (gel->Dimension() != 1) {
            continue;
        }
        
        int side = 1;
        TPZGeoElSide gelside(gel,side);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        // Seeking for condition 1
        while (gelside!= neigh) {
            if (neigh.Element()->Dimension()==0 && neigh.Element()->MaterialId() == TPZFracData::EBCOutlet) {
                break;
            }
            neigh = neigh.Neighbour();
        }
        if (gelside == neigh) {
            continue;
        }
        
        TPZBlock<STATE> & block = fmeshvec[1]->Block();
        TPZGeoElSide gelsideleft(gel,0);

#ifdef PZDEBUG
      {
        const int ncon = cel->NConnects();
        if (ncon != 1) {
          DebugStop();
        }
      }
#endif
        TPZConnect &connect =  cel->Connect(0);
        int seq = connect.SequenceNumber();
#ifdef PZDEBUG
        if (block.Size(seq) != 1) {
            DebugStop();
        }
#endif
        int pos = block.Position(seq);
        fmeshvec[1]->Solution()(pos,0) = fData->SigmaConf();
        
        TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
        an->LoadSolution(fcmeshMixed->Solution());
        return;
    }
}

void TPZDarcyAnalysis::CreateGeoMeshQuad(int npropag, int nrefy, REAL ly)
{
  fgmesh = new TPZGeoMesh;
  const int matiddarcy = TPZFracData::EMatDarcy, bcbot = TPZFracData::EBCBottom, bcright = TPZFracData::EBCRight, bctop = TPZFracData::EBCTop, bcleft = TPZFracData::EBCLeft;
  
  
  int ndivV = npropag;
  int ndivH = nrefy;
  
  int64_t ncols = ndivV + 1;
  int64_t nrows = ndivH + 1;
  int64_t nnodes = nrows*ncols;
  
  fgmesh->NodeVec().Resize(nnodes);
  
  REAL deltadivV = fData->ElSize();
  REAL deltandivH = ly/nrefy;
  
  int64_t nid = 0;
  for(int64_t r = 0; r < nrows; r++)
  {
    for(int64_t c = 0; c < ncols; c++)
    {
      REAL x = c*deltadivV;
      REAL y = r*deltandivH;
      
      TPZVec<REAL> coord(3,0.);
      coord[0] = x;
      coord[1] = y;
      fgmesh->NodeVec()[r*ncols + c].SetCoord(coord);
      fgmesh->NodeVec()[r*ncols + c].SetNodeId(nid);
      nid++;
    }
  }
  
  TPZGeoEl * gel = NULL;
  TPZVec<int64_t> topol(4);
  int64_t indx = 0;
  for(int64_t r = 0; r < nrows-1; r++)
  {
    for(int64_t c = 0; c < ncols-1; c++)
    {
      topol[0] = r*(ncols) + c;
      topol[1] = r*(ncols) + c + 1;
      topol[2] = r*(ncols) + c + 1 + ncols;
      topol[3] = r*(ncols) + c + ncols;
      
      gel = fgmesh->CreateGeoElement(EQuadrilateral, topol, matiddarcy, indx);
      gel->SetId(indx);
      indx++;
    }
  }
  
  fgmesh->BuildConnectivity();
  
  int64_t nelem = fgmesh->NElements();
  for(int64_t el = 0; el < nelem; el++)
  {
    TPZGeoEl * gel = fgmesh->ElementVec()[el];
    
    //south BC
    TPZGeoElSide sideS(gel,4);
    TPZGeoElSide neighS(sideS.Neighbour());
    if(sideS == neighS)
    {
      if(gel->MaterialId() == matiddarcy)
      {
        gel->CreateBCGeoEl(4, bcbot);
      }
      else
      {
        DebugStop();
      }
    }
    
    //east BC
    TPZGeoElSide sideE(gel,5);
    TPZGeoElSide neighE(sideE.Neighbour());
    if(sideE == neighE)
    {
      gel->CreateBCGeoEl(5, bcright);
    }
    
    //north BC
    TPZGeoElSide sideN(gel,6);
    TPZGeoElSide neighN(sideN.Neighbour());
    if(sideN == neighN)
    {
      if(gel->MaterialId() == matiddarcy)
      {
        gel->CreateBCGeoEl(6, bctop);
      }
      else
      {
        DebugStop();
      }
    }
    
    //west BC
    TPZGeoElSide sideW(gel,7);
    TPZGeoElSide neighW(sideW.Neighbour());
    if(sideW == neighW)
    {
      gel->CreateBCGeoEl(7, bcleft);
    }
  }
  
  fgmesh->BuildConnectivity();
}

void TPZDarcyAnalysis::SwitchTipElement(TPZCompEl * cel, TPZCompEl *celpoint, TPZGeoEl *gelpoint)
{
  
  TPZMultiphysicsElement * mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
  TPZMultiphysicsElement * mcelpoint = dynamic_cast<TPZMultiphysicsElement *>(celpoint);
  if (!mcel || !mcelpoint){
    DebugStop();
  }
  
  TPZGeoEl *gel = mcel->Reference();
  TPZGeoElSide gelside(gel, gel->NSides()-1);
  TPZGeoElSide neigh = gelside.Neighbour();
  
  // Cleanning interfaces elements
  fgmesh->ResetReference();
  fcmeshMixed->LoadReferences();
  while (gelside != neigh) {
    TPZMultiphysicsInterfaceElement * intel = dynamic_cast<TPZMultiphysicsInterfaceElement *> (neigh.Element()->Reference());
    if (intel) {
      delete intel;
      break;
    }
    neigh = neigh.Neighbour();
  }
  
  
  int64_t nmesh = mcel->NMeshes();
  int64_t nmeshpt = mcelpoint->NMeshes();
  if (nmesh != 2 || nmeshpt != 2) {
    DebugStop();
  }
  
  if (!mcel->Element(0) || mcel->Element(0)->NConnects() != 1) {
    DebugStop();
  }
  
  delete mcel->Element(0);
  delete mcelpoint->Element(0);
  
  fgmesh->ResetReference();
  int nel = fmeshvec[0]->NElements();
  for (int iel=0; iel < nel; iel++) {
    TPZCompEl *cel = fmeshvec[0]->Element(iel);
    if (!cel) continue;
    if (cel->Reference()->MaterialId() == TPZFracData::EMatFrac) {
      cel->LoadElementReference();
    }
  }
  
  fmeshvec[0]->SetAllCreateFunctionsContinuous();
  fmeshvec[0]->SetDefaultOrder(fData->PorderFlow());

  int64_t index;
  TPZCompEl *celh1 = fmeshvec[0]->CreateCompEl(gel, index);
  TPZCompEl *celpth1 = fmeshvec[0]->CreateCompEl(gelpoint, index);

  
  fgmesh->ResetReference();
  if (mcel->Element(1))
  {
    delete mcel->Element(1);
  }
  if (mcelpoint->Element(1))
  {
    delete mcelpoint->Element(1);
  }
  
  int dim = 1;

  fmeshvec[1]->SetDimModel(dim);
  fmeshvec[1]->SetAllCreateFunctionsDiscontinuous();
  fmeshvec[1]->SetDefaultOrder(fData->PorderPressure());
  TPZCompEl *celDisc = fmeshvec[1]->CreateCompEl(gel, index);
  TPZCompEl *celDiscpt = fmeshvec[1]->CreateCompEl(gelpoint, index);
  fmeshvec[1]->SetDimModel(2);
  
  mcel->AddElement(celh1, 0);
  mcel->AddElement(celDisc, 1);
  
  fcmeshMixed->SetAllCreateFunctionsMultiphysicElem();
  fgmesh->ResetReference();
  fcmeshMixed->LoadReferences();
  TPZCompEl *celmpbc = fcmeshMixed->CreateCompEl(gelpoint, index);
  TPZMultiphysicsElement *mcelptnew = dynamic_cast<TPZMultiphysicsElement *>(celmpbc);
  mcelptnew->AddElement(celpth1,0);
  mcelptnew->AddElement(celDiscpt,1);
  delete mcelpoint;
  
}

void TPZDarcyAnalysis::SetInterfaceConnects()
{
  for (int64_t iel =0 ; iel < fcmeshMixed->NElements(); iel++) {
    TPZCompEl *cel = fcmeshMixed->Element(iel);
    if (!cel) {
      continue;
    }
    
    TPZMultiphysicsInterfaceElement *interfaceEl = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
    if (!interfaceEl) {
      continue;
    }
    
    TPZCompElSide celleft, celright;
    interfaceEl->GetLeftRightElement(celleft, celright);
    interfaceEl->SetLeftRightElement(celleft, celright);
    
  }
}

void TPZDarcyAnalysis::SwitchBCInFrontOfFrac(TPZGeoEl * gel)
{
  if (!gel) {
    DebugStop();
  }
  
  TPZGeoElSide gelside(gel,1);
  TPZGeoElSide neigh = gelside.Neighbour();
  fgmesh->ResetReference();
  int count = 0;
  if (fcmeshMixed){
    fcmeshMixed->LoadReferences();
  }
  while (gelside != neigh)
  {
    if (neigh.Element()->MaterialId() == TPZFracData::EBCBottom && neigh.Element()->Dimension() == 1 ) {
      if (fcmeshMixed) {
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(neigh.Element()->Reference());
        TPZMultiphysicsInterfaceElement * intmcel = dynamic_cast<TPZMultiphysicsInterfaceElement *>(neigh.Element()->Reference());
        if (mcel) {
          neigh.Element()->SetMaterialId(TPZFracData::EBCAuxBottom);
          count++;
        }
        if (intmcel) {
          neigh.Element()->SetMaterialId(TPZFracData::EBCAuxBottom);
          count++;
        }
        
        if (count == 2) {
          break;
        }
        
      }
      else{
        neigh.Element()->SetMaterialId(TPZFracData::EBCAuxBottom);
        break;
      }
    }
    neigh = neigh.Neighbour();
  }
  if (gelside == neigh){
    DebugStop();
  }
  
}


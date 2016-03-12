//
//  TRMSpaceOdissey.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpaceOdissey.h"
#include "TRMFlowConstants.h"

#include "TRMMixedDarcy.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"
#include "TRMPhaseTransport.h"
#include "TRMPhaseInterfaceTransport.h"

#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "tpzhierarquicalgrid.h"
#include "pzgeopoint.h"
#include "TRMSimworxMeshGenerator.h"
#include "TPZCompMeshTools.h"
#include "pzelchdivbound2.h"
#include "pzshapequad.h"


void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
    f.Resize(1,0.);
    f[0] = pt[2];
    return;
}

static void CreateExampleRawData(TRMRawData &data)
{
    data.fLw = 500.;
    data.fHasLiner = false; //AQUINATHAN esta false para gerar uma malha sem os refinamentos do meio que geram hangnodes
    data.fHasCasing = false; //AQUINATHAN esta false para gerar uma malha sem os refinamentos do meio que geram hangnodes
    
    data.fReservoirWidth = 500.;
    data.fReservoirLength = 1000.;
    data.fReservoirHeight = 50.;
    data.fProdVertPosition = 25;
    data.fWellDiam = 0.2159;
}


/** @brief Default constructor */
TRMSpaceOdissey::TRMSpaceOdissey() : fMeshType(TRMSpaceOdissey::EBox), fPOrder(4)
{
    fGeoMesh                = NULL;
    fH1Cmesh                = NULL;
    fSimulationData         = NULL;
    fFluxCmesh              = NULL;
    fPressureCmesh          = NULL;
    fMixedFluxPressureCmesh = NULL;
    fSaturationMesh         = NULL;
    fGeoMechanicsCmesh      = NULL;
    fTransferGenerator      = new TRMBuildTransfers;
    
}

/** @brief Default desconstructor */
TRMSpaceOdissey::~TRMSpaceOdissey(){
    
    /** @brief H1 computational mesh for validation */
    if(fH1Cmesh)
    {
        fH1Cmesh->CleanUp();
    }
    
    /** @brief Mixed computational mesh for a dual analysis */
    if(fMixedFluxPressureCmesh)
    {
        fMixedFluxPressureCmesh->CleanUp();
    }
    
    /** @brief Autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief Hdiv computational mesh conservative vector field */
    if(fFluxCmesh)
    {
        fFluxCmesh->CleanUp();
    }
    
    /** @brief L2 computational mesh the restriction equation */
    if(fPressureCmesh)
    {
        fPressureCmesh->CleanUp();
    }
    
    /** @brief H1 computational mesh for Maurice Biot Linear Poroelasticity */
    if(fGeoMechanicsCmesh)
    {
        fGeoMechanicsCmesh->CleanUp();
    }

}


/** @brief Create a Hdiv computational mesh Hdiv */
void TRMSpaceOdissey::CreateFluxCmesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int qorder = fPOrder;
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,1000.);
    
    // Malha computacional
    fFluxCmesh = new TPZCompMesh(fGeoMesh);
    
    TRMMixedDarcy * mat = new TRMMixedDarcy(_ReservMatId);
    fFluxCmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, _ConfinementReservBCbottom, typePressure, val1, val2Pressure);
    fFluxCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, _ConfinementReservBCtop, typePressure, val1, val2Pressure);
    fFluxCmesh->InsertMaterialObject(bcS);
    
    // Bc E
    TPZBndCond * bcE = mat->CreateBC(mat, _LateralReservBC, typePressure, val1, val2Pressure);
    fFluxCmesh->InsertMaterialObject(bcE);
    
    // Bc W
//    TPZBndCond * bcW = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
//    fFluxCmesh->InsertMaterialObject(bcW);
    
    // Bc B
//    TPZBndCond * bcB = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
//    fFluxCmesh->InsertMaterialObject(bcB);
    
    // Bc T
//    TPZBndCond * bcT = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
//    fFluxCmesh->InsertMaterialObject(bcT);
    

    TPZBndCond * bcToe = mat->CreateBC(mat, _WellToeMatId, typeFlux, val1, val2Pressure);
    fFluxCmesh->InsertMaterialObject(bcToe);
    
    TPZBndCond * bcHeel = mat->CreateBC(mat, _WellHeelMatId, typePressure, val1, val2Pressure);
    fFluxCmesh->InsertMaterialObject(bcHeel);
    
    /*
    TPZBndCond * bcWellRes = mat->CreateBC(mat, _WellFacesMatId, typePressure, val1, val2Pressure);
    fFluxCmesh->InsertMaterialObject(bcWellRes);
     */
    
    TPZBndCond * bcWellFaces = mat->CreateBC(mat, _Well3DReservoirFaces, typePressure, val1, val2Pressure);
    fFluxCmesh->InsertMaterialObject(bcWellFaces);
     

    

    
    // Setando Hdiv
    fFluxCmesh->SetDimModel(dim);
    fFluxCmesh->SetDefaultOrder(qorder);
    fFluxCmesh->SetAllCreateFunctionsHDiv();
    fFluxCmesh->AutoBuild();
    
    long nel = fFluxCmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fFluxCmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if (!gel || (matid != _Well3DReservoirFaces && matid != _WellHeelMatId && matid != _WellToeMatId)) {
            continue;
        }
        TPZCompElHDivBound2<pzshape::TPZShapeQuad> *bound = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(cel);
        if (!bound) {
            DebugStop();
        }
        bound->SetSideOrient(8, -1);
    }
    
#ifdef PZDEBUG
    std::ofstream out("CmeshFlux.txt");
    fFluxCmesh->Print(out);
#endif
    
}

void PressFunc(const TPZVec<REAL> &x, TPZVec<STATE> &func)
{
    func[0] = 0.;
}


/** @brief Create a Discontinuous computational mesh L2 */
void TRMSpaceOdissey::CreatePressureCmesh(){
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int porder = fPOrder;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2Pressure(1,1,0.), val2Flux(1,1,1000.);
    
    // Malha computacional
    fPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    TRMMixedDarcy * mat = new TRMMixedDarcy(_ReservMatId);
    fPressureCmesh->InsertMaterialObject(mat);

    // Setando L2
    fPressureCmesh->SetDimModel(dim);
    fPressureCmesh->SetDefaultOrder(porder);
    
    fPressureCmesh->SetAllCreateFunctionsContinuous();
    fPressureCmesh->ApproxSpace().CreateDisconnectedElements(true);
    fPressureCmesh->AutoBuild();
    
    fPressureCmesh->AdjustBoundaryElements();
    fPressureCmesh->CleanUpUnconnectedNodes();
    
    int ncon = fPressureCmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = fPressureCmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    TPZDummyFunction<STATE> dummy(PressFunc);
    TPZCompMeshTools::LoadSolution(fPressureCmesh.operator->(), dummy);
#ifdef PZDEBUG
    std::ofstream out("../CmeshPress.txt");
    fPressureCmesh->Print(out);
#endif
    
}

void One(const TPZVec<REAL> &x, TPZVec<STATE> &f)
{
    f[0] = 3.*M_PI*M_PI*sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}

/** @brief Create a Mixed computational mesh Hdiv-L2 */
void TRMSpaceOdissey::CreateMixedCmesh(){
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    bool StaticCondensation = false;
    int dim = 3;
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(2,2,0.), val2Flux(2,1,0.), val2Pressure(2,1,0.);
    val2Pressure(0,0) = 1000.;
    
    // Malha computacional
    fMixedFluxPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    // Material medio poroso
    TRMMixedDarcy * mat = new TRMMixedDarcy(_ReservMatId);
//    mat->SetForcingFunction(One);
    fMixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, _ConfinementReservBCbottom, typePressure, val1, val2Pressure);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(Forcing);
    bcN->SetForcingFunction(0,force);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, _ConfinementReservBCtop, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    // Bc E
    val2Flux(0,0) = 0.0;
    TPZBndCond * bcE = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2Flux);
    
    // Bc W
    val2Flux(0,0) = 0.0;
    TPZBndCond * bcW = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2Flux);
    
    // Bc B
    TPZBndCond * bcB = mat->CreateBC(mat, _LateralReservBC, typePressure, val1, val2Pressure);
    bcB->SetForcingFunction(0, force);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcB);
    
    // Bc T
    val2Flux(0,0) = 0.0;
    TPZBndCond * bcT = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2Flux);
    
    
    TPZBndCond * bcToe = mat->CreateBC(mat, _WellToeMatId, typePressure, val1, val2Pressure);
    bcToe->SetForcingFunction(0, force);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcToe);
    
    TPZBndCond * bcHeel = mat->CreateBC(mat, _WellHeelMatId, typePressure, val1, val2Pressure);
    bcHeel->SetForcingFunction(0, force);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcHeel);
    
    /*
    TPZBndCond * bcWellRes = mat->CreateBC(mat, _WellFacesMatId, typePressure, val1, val2Pressure);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcWellRes);
    */
    
    TPZBndCond * bcWellFaces = mat->CreateBC(mat, _Well3DReservoirFaces, typePressure, val1, val2Pressure);
    bcWellFaces->SetForcingFunction(0, force);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcWellFaces);

    
    fMixedFluxPressureCmesh->SetDimModel(dim);
    fMixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fMixedFluxPressureCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    
    
    meshvector[0] = fFluxCmesh.operator->();
    meshvector[1] = fPressureCmesh.operator->();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fMixedFluxPressureCmesh.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fMixedFluxPressureCmesh.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fMixedFluxPressureCmesh.operator->());
    
    long nel = fMixedFluxPressureCmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fMixedFluxPressureCmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
}

/** @brief Statically condense the internal equations of the elements */
void TRMSpaceOdissey::StaticallyCondenseEquations()
{
    if (!fMixedFluxPressureCmesh.operator->()) {
        std::cout<< "No multiphysic computational mesh " << std::endl;
        DebugStop();
    }
    
    
    fMixedFluxPressureCmesh.operator->()->Reference()->ResetReference();
    fMixedFluxPressureCmesh.operator->()->LoadReferences();
    
    fMixedFluxPressureCmesh.operator->()->ComputeNodElCon();
    // create condensed elements
    // increase the NumElConnected of one pressure connects in order to prevent condensation
    for (long icel=0; icel < fMixedFluxPressureCmesh.operator->()->NElements(); icel++) {
        TPZCompEl  * cel = fMixedFluxPressureCmesh.operator->()->Element(icel);
        
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
            
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    meshvector[0] = fFluxCmesh.operator->();
    meshvector[1] = fPressureCmesh.operator->();
    
    int DOF = meshvector[0]->NEquations() + meshvector[1]->NEquations();
    REAL PercentCondensedDOF = 100.0*(1.0 - REAL(fMixedFluxPressureCmesh.operator->()->NEquations())/REAL(DOF));
    std::cout << "Degrees of freedom: " << DOF << std::endl;
    std::cout << "Percent of condensed Degrees of freedom: " << PercentCondensedDOF << std::endl;
}



/** @brief Create a H1 computational mesh */
void TRMSpaceOdissey::CreateH1Cmesh()
{
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int porder  = fPOrder;
    fH1Cmesh = new TPZCompMesh(fGeoMesh);
    fH1Cmesh->SetDimModel(3);
    
    TPZMatLaplacian *material = new TPZMatLaplacian(_ReservMatId,3);
    TPZAutoPointer<TPZFunction<STATE> > one = new TPZDummyFunction<STATE>(One);
    fH1Cmesh->InsertMaterialObject(material);

    TPZFNMatrix<1> val1(1,1,0.),val2(1,1,0);
    TPZBndCond *inflow = new TPZBndCond(material,_ConfinementReservBCbottom,0,val1,val2);
    val2(0,0) = 0.;
    TPZBndCond *outflow = new TPZBndCond(material,_ConfinementReservBCtop,0,val1,val2);
    
    // Bc B
    TPZBndCond * bcB = material->CreateBC(material, _LateralReservBC, 0, val1, val2);
//    bcB->SetForcingFunction(0, force);
    fH1Cmesh->InsertMaterialObject(bcB);

    fH1Cmesh->InsertMaterialObject(inflow);
    fH1Cmesh->InsertMaterialObject(outflow);
    fH1Cmesh->SetDefaultOrder(porder);
    
    TPZCreateApproximationSpace space;
    space.SetAllCreateFunctionsContinuous();    
    fH1Cmesh->ApproxSpace() = space;
    
    fH1Cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("../CmeshPressH1.txt");
    fH1Cmesh->Print(out);
#endif
    
}

/** @brief Create a computational mesh L2 */
void TRMSpaceOdissey::CreateTransportMesh()
{
    // criamos elementos descontinuos e de interface com memoria...
    if (fSaturationMesh) {
        DebugStop();
    }
    fSaturationMesh = new TPZCompMesh(fGeoMesh);
    fSaturationMesh->SetDimModel(3);
    fSaturationMesh->SetDefaultOrder(0);
    fSaturationMesh->SetAllCreateFunctionsDiscontinuous();
    
    TRMPhaseTransport *mat = new TRMPhaseTransport(_ReservMatId);
    fSaturationMesh->InsertMaterialObject(mat);
    
    TRMPhaseInterfaceTransport *matint = new TRMPhaseInterfaceTransport(_ReservoirInterface);
    fGeoMesh->AddInterfaceMaterial(_ReservMatId, _ReservMatId,_ReservoirInterface);
    
    // WE NEED TO ADD THE BOUNDARY CONDITION MATERIALS
    DebugStop();
    
    fSaturationMesh->ApproxSpace().CreateInterfaces(fSaturationMesh);
}



void TRMSpaceOdissey::PrintGeometry()
{
    //  Print Geometrical Base Mesh
    std::ofstream argument("../GeometricMesh.txt");
    fGeoMesh->Print(argument);
    std::ofstream Dummyfile("../GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh,Dummyfile, true);
}

/** @brief Create a reservoir-box geometry */
void TRMSpaceOdissey::CreateGeometricBoxMesh(TPZManVector<int,2> dx, TPZManVector<int,2> dy, TPZManVector<int,2> dz){
    
    REAL t=0.0;
    REAL dt;
    int n;
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh0D = new TPZGeoMesh;
    GeoMesh0D->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh0D->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,_ReservMatId,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom0D(GeoMesh0D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncX = new TPZDummyFunction<STATE>(ParametricfunctionX);
    CreateGridFrom0D.SetParametricFunction(ParFuncX);
    CreateGridFrom0D.SetFrontBackMatId(_ConfinementReservBCbottom,_ConfinementReservBCtop);
    
    dt = dx[0];
    n = dx[1];
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh1D = CreateGridFrom0D.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom1D(GeoMesh1D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncY = new TPZDummyFunction<STATE>(ParametricfunctionY);
    CreateGridFrom1D.SetParametricFunction(ParFuncY);
    CreateGridFrom1D.SetFrontBackMatId(_LateralReservBC,_LateralReservBC);
//    CreateGridFrom1D.SetTriangleExtrusion();
    
    dt = dy[0];
    n = dy[1];
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh2D = CreateGridFrom1D.ComputeExtrusion(t, dt, n);
        
    TPZHierarquicalGrid CreateGridFrom2D(GeoMesh2D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncZ = new TPZDummyFunction<STATE>(ParametricfunctionZ);
    CreateGridFrom2D.SetParametricFunction(ParFuncZ);
    CreateGridFrom2D.SetFrontBackMatId(_LateralReservBC,_LateralReservBC);
//    CreateGridFrom2D.SetPrismExtrusion();
    
    dt = dz[0];
    n = dz[1];
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fGeoMesh = CreateGridFrom2D.ComputeExtrusion(t, dt, n);
    
}
void TRMSpaceOdissey::ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0];
    X[1] = 0.0;
    X[2] = 0.0*25.0*sin(0.1*par[0]);
}

void TRMSpaceOdissey::ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void TRMSpaceOdissey::ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}



/** @brief Create the reservoir geometry */
void TRMSpaceOdissey::CreateGeometricReservoirMesh(){
    
    gRefDBase.ReadRefPatternDBase("../RefPatterns.rpt");
    TRMRawData rawdata;
    CreateExampleRawData(rawdata);
    TRMSimworxMeshGenerator meshGen;
    const bool withwellbc = true;
    fGeoMesh = meshGen.CreateSimworxGeoMesh(rawdata,withwellbc);
}

/** @brief Configure the boundary conditions of a well with reservoir boundary conditions */
void TRMSpaceOdissey::ConfigureWellConstantPressure(STATE wellpressure, STATE farfieldpressure)
{
    
    TPZMaterial *mat = fMixedFluxPressureCmesh->FindMaterial(_ReservMatId);
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2PressureFarField(1,1,farfieldpressure), val2PressureWell(1,1,wellpressure);

    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, _ConfinementReservBCbottom, typeFlux, val1, val2Flux);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, _ConfinementReservBCtop, typeFlux, val1, val2Flux);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    // Bc E
    //    val2(0,0) = 0.0;
    //    TPZBndCond * bcE = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
    
    // Bc W
    //    val2(0,0) = 0.0;
    //    TPZBndCond * bcW = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
    
    // Bc B
    TPZBndCond * bcB = mat->CreateBC(mat, _LateralReservBC, typePressure, val1, val2PressureFarField);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcB);
    // Bc T
    //    val2(0,0) = 0.0;
    //    TPZBndCond * bcT = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
    
    
    TPZBndCond * bcToe = mat->CreateBC(mat, _WellToeMatId, typeFlux, val1, val2Flux);
//    TPZBndCond * bcToe = mat->CreateBC(mat, _WellToeMatId, typePressure, val1, val2PressureWell);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcToe);
    
    TPZBndCond * bcHeel = mat->CreateBC(mat, _WellHeelMatId, typeFlux, val1, val2Flux);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcHeel);
    
    /*
     TPZBndCond * bcWellRes = mat->CreateBC(mat, _WellFacesMatId, typePressure, val1, val2Pressure);
     fMixedFluxPressureCmesh->InsertMaterialObject(bcWellRes);
     */
    
    TPZBndCond * bcWellFaces = mat->CreateBC(mat, _Well3DReservoirFaces, typePressure, val1, val2PressureWell);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcWellFaces);
    
}

static void IncludeNeighbours(TPZAutoPointer<TPZCompMesh> cmesh, long index, std::map<long,int> &extended)
{
    TPZCompEl *cel = cmesh->Element(index);
    TPZGeoEl *gel = cel->Reference();
    int nsides = gel->NSides();
    for (int is=0; is<nsides; is++) {
        if (gel->SideDimension(is) != 2) {
            continue;
        }
        TPZGeoElSide gelside(gel,is);
        TPZGeoElSide neighbour = gelside.Neighbour();
        TPZCompElSide celside = neighbour.Reference();
        if (!celside) {
            continue;
        }
        long celindex = celside.Element()->Index();
        if (extended.find(celindex) == extended.end()) {
            std::cout << "Including index " << celindex << " order " << extended[index] << std::endl;
            extended[celindex] = extended[index];
            // if celindex has neighbours of dimension 2, include them also
            TPZGeoEl *gelindex = neighbour.Element();
            // loop over all sides of dimension 2
            for (int is=0; is<gelindex->NSides(); is++) {
                if (gelindex->SideDimension(is) != 2) {
                    continue;
                }
                TPZGeoElSide gelside(gelindex,is);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->Dimension() == 2) {
                        long neighindex = neighbour.Element()->Reference()->Index();
                        std::cout << "Including index " << neighindex << " order " << extended[index] << std::endl;
                        extended[neighindex] = extended[index];
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
        }
        
    }

}

void TRMSpaceOdissey::ModifyElementOrders(std::map<long,int> &elorders)
{
    // settle the orders of the pressure elements first
    this->fMixedFluxPressureCmesh->Reference()->ResetReference();
    for (std::map<long,int>::iterator it = elorders.begin(); it != elorders.end(); it++) {
        long elindex = it->first;
        TPZCompEl *cel = fMixedFluxPressureCmesh->Element(elindex);
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mcel) {
            DebugStop();
        }
        TPZCompEl *press = mcel->Element(1);
        if (!press) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(press);
        intel->PRefine(it->second);
    }
    fPressureCmesh->ExpandSolution();
    fFluxCmesh->LoadReferences();
    for (std::map<long,int>::iterator it = elorders.begin(); it != elorders.end(); it++) {
        long elindex = it->first;
        TPZCompEl *cel = fMixedFluxPressureCmesh->Element(elindex);
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mcel) {
            DebugStop();
        }
        TPZCompEl *press = mcel->Element(0);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(press);
        intel->PRefine(it->second);
    }
    fFluxCmesh->ExpandSolution();
    
    CreateMixedCmesh();

}

/// Adjust the polinomial order of the elements
void TRMSpaceOdissey::IncreaseOrderAroundWell(int numlayers)
{
    fGeoMesh->ResetReference();
    this->fMixedFluxPressureCmesh->LoadReferences();
    
    // find well toe element
    //TPZGeoEl *geltoe = 0;
    TPZManVector<TPZGeoEl*,2> vecHellToe(2,NULL);
    long nelem = fGeoMesh->NElements();
    int ilocal = 0;
    for (long el=0; el < nelem; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if (gel && gel->MaterialId() == _WellToeMatId) {
            vecHellToe[ilocal] = gel;
            //geltoe = gel;
            ilocal++;
        }
        if (gel && gel->MaterialId() == _WellHeelMatId) {
            vecHellToe[ilocal] = gel;
            //geltoe = gel;
            ilocal++;
        }

        if (ilocal == 2){
            break;
        }
    }
    
    if(ilocal != 2)
    {
        DebugStop();
    }
    
    /// find a neighbouring element of type reservoir matid
    for (int i = 0; i < vecHellToe.NElements(); i++) {
        if (vecHellToe[i]->Type() != EQuadrilateral) {
            DebugStop();
        }

    }
    
    TPZManVector<TPZGeoElSide,2> vecGeoSide(2);
    TPZManVector<TPZGeoElSide,2> vecGeoNeigh(2);
    
    for (int i = 0; i < vecGeoSide.size(); i++) {
        vecGeoSide[i] = TPZGeoElSide(vecHellToe[i],vecHellToe[i]->NSides()-1);
        vecGeoNeigh[i] = vecGeoSide[i].Neighbour();
    }
    for (int i = 0; i < vecGeoSide.size(); i++) {
        while (vecGeoNeigh[i] != vecGeoSide[i]) {
            if (vecGeoNeigh[i].Element()->MaterialId() == _ReservMatId) {
                break;
            }
            vecGeoNeigh[i] = vecGeoNeigh[i].Neighbour();
        }
    }

    int quadside = -1;
    for (int i = 0; i < vecGeoSide.size(); i++) {
        if (vecGeoSide[i].Element()->MaterialId() == _WellToeMatId) {
            quadside = 20;
        }
        else if (vecGeoSide[i].Element()->MaterialId() == _WellHeelMatId){
            quadside = 25;
        }
        else{
            DebugStop();
        }
        if (vecGeoNeigh[i] == vecGeoSide[i] || vecGeoNeigh[i].Side() != quadside) {
            DebugStop();
        }
    }
    // go up the refinement tree and assign an order to the included computational elements
    TPZManVector<TPZGeoEl*,2> vecGelBase(2);

    for (int i = 0; i < vecGelBase.size(); i++) {
        vecGelBase[i] = vecGeoNeigh[i].Element();
    }
    
    TPZManVector<std::map<long,int>,2> contemplated(2);
//    std::cout << "Including index " << gelbase->Reference()->Index() << " order " << fPOrder+numlayers << std::endl;

    
    contemplated[0][vecGelBase[0]->Reference()->Index()]= fPOrder+numlayers;
    contemplated[1][vecGelBase[1]->Reference()->Index()]= fPOrder+numlayers;
    int refside = -1;
    for (int ic = 0; ic < vecGelBase.size(); ic++) {
        if (vecGeoSide[ic].Element()->MaterialId() == _WellToeMatId) {
            refside = 25;
        }
        else if (vecGeoSide[ic].Element()->MaterialId() == _WellHeelMatId) {
            refside = 20;
        }
        else{
            DebugStop();
        }
        for (int i=1; i< numlayers; i++) {
            vecGelBase[ic] = vecGelBase[ic]->Neighbour(refside).Element();
            if (!vecGelBase[ic]) {
                DebugStop();
            }
            std::cout << "Including index " << vecGelBase[ic]->Reference()->Index() << " order " << fPOrder+numlayers-i << std::endl;
            contemplated[ic][vecGelBase[ic]->Reference()->Index()] = fPOrder+numlayers-i;
        }

    }
    // include the neighbours of the elements within contemplated along faces
    
    TPZManVector<std::map<long,int>,2> original(2);
    //TPZManVector<std::map<long,int>,2> original(contemplated[0]);
    //std::map<long,int> original2(contemplated[1]);
    
    for (int i = 0; i < contemplated.size(); i++) {
        original[i] = contemplated[i];
    }
    
    for (int i = 0; i < original.size(); i++) {
        for (std::map<long,int>::iterator it = original[i].begin(); it != original[i].end(); it++) {
            IncludeNeighbours(this->fMixedFluxPressureCmesh, it->first, contemplated[i]);
        }
    }
    for (int i = 0 ; i < contemplated.size() ; i++){
        ModifyElementOrders(contemplated[i]);
    }
}


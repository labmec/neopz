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
TRMSpaceOdissey::TRMSpaceOdissey() : fMeshType(TRMSpaceOdissey::EBox)
{
    
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
    int qorder = 1;
    
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
    
#ifdef DEBUG
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
    int porder = 1;
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Pressure(1,1,0.), val2Flux(1,1,1000.);
    
    // Malha computacional
    fPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    TRMMixedDarcy * mat = new TRMMixedDarcy(_ReservMatId);
    fPressureCmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, _ConfinementReservBCbottom, typePressure, val1, val2Pressure);
    fPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, _ConfinementReservBCtop, typePressure, val1, val2Pressure);
    fPressureCmesh->InsertMaterialObject(bcS);
    
    // Bc E
    TPZBndCond * bcE = mat->CreateBC(mat, _LateralReservBC, typePressure, val1, val2Pressure);
    fPressureCmesh->InsertMaterialObject(bcE);
    
    // Bc W
//    TPZBndCond * bcW = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
//    fPressureCmesh->InsertMaterialObject(bcW);
    
    // Bc B
//    TPZBndCond * bcB = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
//    fPressureCmesh->InsertMaterialObject(bcB);
    
    // Bc T
//    TPZBndCond * bcT = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
//    fPressureCmesh->InsertMaterialObject(bcT);
    
    TPZBndCond * bcToe = mat->CreateBC(mat, _WellToeMatId, typeFlux, val1, val2Pressure);
    fPressureCmesh->InsertMaterialObject(bcToe);
    
    TPZBndCond * bcHeel = mat->CreateBC(mat, _WellHeelMatId, typePressure, val1, val2Pressure);
    fPressureCmesh->InsertMaterialObject(bcHeel);
    
    TPZBndCond * bcWellRes = mat->CreateBC(mat, _Well3DReservoirFaces, typePressure, val1, val2Pressure);
    fPressureCmesh->InsertMaterialObject(bcWellRes);
    

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
#ifdef DEBUG
    std::ofstream out("../CmeshPress.txt");
    fPressureCmesh->Print(out);
#endif
    
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
//    val2(0,0) = 0.0;
//    TPZBndCond * bcE = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
    
    // Bc W
//    val2(0,0) = 0.0;
//    TPZBndCond * bcW = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
    
    // Bc B
    TPZBndCond * bcB = mat->CreateBC(mat, _LateralReservBC, typePressure, val1, val2Pressure);
    bcB->SetForcingFunction(0, force);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcB);
    // Bc T
//    val2(0,0) = 0.0;
//    TPZBndCond * bcT = mat->CreateBC(mat, _LateralReservBC, typeFlux, val1, val2);
    
    
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
    fMixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    fMixedFluxPressureCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    
    this->CreateFluxCmesh();
    this->CreatePressureCmesh();
    
    meshvector[0] = fFluxCmesh.operator->();
    meshvector[1] = fPressureCmesh.operator->();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fMixedFluxPressureCmesh.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fMixedFluxPressureCmesh.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fMixedFluxPressureCmesh.operator->());
    
    if (StaticCondensation) {
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
        
        int DOF = meshvector[0]->NEquations() + meshvector[1]->NEquations();
        REAL PercentCondensedDOF = 100.0*(1.0 - REAL(fMixedFluxPressureCmesh.operator->()->NEquations())/REAL(DOF));
        std::cout << "Degrees of freedom: " << DOF << std::endl;
        std::cout << "Percent of condensed Degrees of freedom: " << PercentCondensedDOF << std::endl;

    }
    
    
    
}



/** @brief Create a H1 computational mesh */
void TRMSpaceOdissey::CreateH1Cmesh()
{
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int porder  = 1;
    fH1Cmesh = new TPZCompMesh(fGeoMesh);
    fH1Cmesh->SetDimModel(3);
    
    TPZMatLaplacian *material = new TPZMatLaplacian(_ReservMatId,3);
    fH1Cmesh->InsertMaterialObject(material);

    TPZFNMatrix<1> val1(1,1,0.),val2(1,1,20);
    TPZBndCond *inflow = new TPZBndCond(material,_ConfinementReservBCbottom,0,val1,val2);
    val2(0,0) = 10.;
    TPZBndCond *outflow = new TPZBndCond(material,_ConfinementReservBCtop,0,val1,val2);
    
    fH1Cmesh->InsertMaterialObject(inflow);
    fH1Cmesh->InsertMaterialObject(outflow);
    fH1Cmesh->SetDefaultOrder(porder);
    
    TPZCreateApproximationSpace space;
    space.SetAllCreateFunctionsContinuous();    
    fH1Cmesh->ApproxSpace() = space;
    
    fH1Cmesh->AutoBuild();
    
#ifdef DEBUG
    std::ofstream out("CmeshPressH1.txt");
    fH1Cmesh->Print(out);
#endif
    
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
    
    GeoMesh2D->Print();
    
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


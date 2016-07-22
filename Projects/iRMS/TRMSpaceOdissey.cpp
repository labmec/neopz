//
//  TRMSpaceOdissey.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpaceOdissey.h"
#include "TRMFlowConstants.h"

#include "TRMMultiphase.h"
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

#include "pzbndcond.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzcompelwithmem.h"


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
   
    fPOrder = 1;
    fSOrder = 0;
    fGeoMesh                    = NULL;
    fSimulationData             = NULL;
    fH1Cmesh                    = NULL;
    fFluxCmesh                  = NULL;
    fPressureCmesh              = NULL;
    fAlphaSaturationMesh        = NULL;
    fBetaSaturationMesh         = NULL;
    fGeoMechanicsCmesh          = NULL;
    fTransportMesh              = NULL;
    fMixedFluxPressureCmesh     = NULL;
    fPressureSaturationCmesh    = NULL;
    fMonolithicMultiphaseCmesh  = NULL;
//    fTransferGenerator          = new TRMBuildTransfers;
    
}

/** @brief Default desconstructor */
TRMSpaceOdissey::~TRMSpaceOdissey(){
    
    if(fH1Cmesh)                    fH1Cmesh->CleanUp();
    if(fFluxCmesh)                  fFluxCmesh->CleanUp();
    if(fPressureCmesh)              fPressureCmesh->CleanUp();
    if(fAlphaSaturationMesh)        fAlphaSaturationMesh->CleanUp();
    if(fBetaSaturationMesh)         fBetaSaturationMesh->CleanUp();
    if(fGeoMechanicsCmesh)          fGeoMechanicsCmesh->CleanUp();
    if(fTransportMesh)              fTransportMesh->CleanUp();
    if(fMixedFluxPressureCmesh)     fMixedFluxPressureCmesh->CleanUp();
    if(fPressureSaturationCmesh)    fPressureSaturationCmesh->CleanUp();
    if(fMonolithicMultiphaseCmesh)  fMonolithicMultiphaseCmesh->CleanUp();
    
}


/** @brief Create a Hdiv computational mesh Hdiv */
void TRMSpaceOdissey::CreateFluxCmesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int flux_or_pressure = 0;
    int qorder = fPOrder;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > bc_item;
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > bc;
    
    // Malha computacional
    fFluxCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id);
        fFluxCmesh->InsertMaterialObject(mat);
        
        // Inserting volumetric materials
        int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
        int bc_id = 0;

        for (int j = 0; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[flux_or_pressure];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            boundary_c->SetTimedependentBCForcingFunction(bc_item.second);
            fFluxCmesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    // Setando Hdiv
    fFluxCmesh->SetDimModel(dim);
    fFluxCmesh->SetDefaultOrder(qorder);
    fFluxCmesh->SetAllCreateFunctionsHDiv();
    fFluxCmesh->AutoBuild();
    
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
    
    // Malha computacional
    fPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id);
        fPressureCmesh->InsertMaterialObject(mat);
        
    }

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
    std::ofstream out("CmeshPress.txt");
    fPressureCmesh->Print(out);
#endif
    
}

void One(const TPZVec<REAL> &x, TPZVec<STATE> &f)
{
    f[0] = 3.*M_PI*M_PI*sin(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
}

/** @brief exact pressure */
void ExactPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &pressure)
{
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    STATE p = (1. - x)*x + (1. - y)*y + (1. - z)*z;//(1. - x)*x*(1. - y)*y*(1. - z)*z;
    pressure[0] = p;
}

/** @brief exact flux */
void ExactFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &flux)
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
void ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &f)
{
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    REAL rhs = 6;//2.*(-1. + x)*x*(-1. + y)*y + 2.*(-1. + x)*x*(-1. + z)*z + 2.*(-1. + y)*y*(-1. + z)*z;
    f[0] = rhs;
}

/** @brief Create a Mixed computational mesh Hdiv-L2 */
void TRMSpaceOdissey::CreateMixedCmesh(){
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }

    int dim = 3;
    int flux_or_pressure = 0;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > bc_item;
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > bc;
    
    // Malha computacional
    fMixedFluxPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id);
        mat->SetSimulationData(fSimulationData);        
        fMixedFluxPressureCmesh->InsertMaterialObject(mat);
        
        // Inserting boundary materials
        int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
        int bc_id = 0;

        for (int j = 0; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[flux_or_pressure];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            boundary_c->SetTimedependentBCForcingFunction(bc_item.second); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fMixedFluxPressureCmesh->InsertMaterialObject(boundary_c);

        }
        
    }

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
    
#ifdef PZDEBUG
    std::ofstream out("CmeshMixed.txt");
    fMixedFluxPressureCmesh->Print(out);
#endif
    
}

/** @brief Create a Mixed-Transport muliphase computational mesh Hdiv-L2-L2-L2 */
void TRMSpaceOdissey::CreateMultiphaseCmesh(){
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int flux_or_pressure = 0;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > bc_item;
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > bc;
    
    // Malha computacional
    fMonolithicMultiphaseCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMultiphase * mat = new TRMMultiphase(rock_id);
        mat->SetSimulationData(fSimulationData);
        fMonolithicMultiphaseCmesh->InsertMaterialObject(mat);
        
        // Inserting boundary materials
        int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
        int bc_id = 0;

        for (int j = 0; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[flux_or_pressure];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            boundary_c->SetTimedependentBCForcingFunction(bc_item.second); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fMonolithicMultiphaseCmesh->InsertMaterialObject(boundary_c);
            
        }
        
    }
    
    
    fMonolithicMultiphaseCmesh->SetDimModel(dim);
    fMonolithicMultiphaseCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fMonolithicMultiphaseCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    
    if (fSimulationData->IsOnePhaseQ()) {
        meshvector[0] = fFluxCmesh.operator->();
        meshvector[1] = fPressureCmesh.operator->();
    }
    
    if (fSimulationData->IsTwoPhaseQ()) {
        meshvector.Resize(3);
        meshvector[0] = fFluxCmesh.operator->();
        meshvector[1] = fPressureCmesh.operator->();
        meshvector[2] = fAlphaSaturationMesh.operator->();
    }

    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fMonolithicMultiphaseCmesh.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fMonolithicMultiphaseCmesh.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fMonolithicMultiphaseCmesh.operator->());
    
    long nel = fMonolithicMultiphaseCmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fMonolithicMultiphaseCmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CmeshMultiphase.txt");
    fMonolithicMultiphaseCmesh->Print(out);
#endif
}

/** @brief Create computational interfaces for jumps  */
void TRMSpaceOdissey::CreateInterfacesInside(TPZAutoPointer<TPZCompMesh> cmesh){ //@Omar:: It is required to robust for several materials in the same hydraulic unit!!!
    
    fGeoMesh->ResetReference();
    cmesh->LoadReferences();

    // Creation of interface elements
    int nel = cmesh->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        TPZGeoEl * gel = compEl->Reference();
        if(!gel) {continue;}
        if(gel->HasSubElement()) {continue;}
        int index = compEl ->Index();
        if(compEl->Dimension() == cmesh->Dimension())
        {
            TPZMultiphysicsElement * inter_el_mult = dynamic_cast<TPZMultiphysicsElement *>(cmesh->ElementVec()[index]);
            if(!inter_el_mult) continue;
            inter_el_mult->CreateInterfaces();
        }
    }
    
#ifdef PZDEBUG
    std::ofstream out("CmeshWithInterfaces.txt");
    cmesh->Print(out);
#endif
    
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
    material->SetForcingFunction(One,fPOrder);
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
    std::ofstream out("CmeshPressH1.txt");
    fH1Cmesh->Print(out);
#endif
    
}

/** @brief Create a computational mesh L2 */
void TRMSpaceOdissey::CreateAlphaTransportMesh()
{
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int saturation = 0;
    int sorder = fSOrder;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
    // Malha computacional
    fAlphaSaturationMesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMPhaseTransport * mat = new TRMPhaseTransport(rock_id);
        fAlphaSaturationMesh->InsertMaterialObject(mat);
        
        // Inserting volumetric materials
        int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
        int bc_id = 0;
        std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > bc_item;
        TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > bc;
        for (int j = 0; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[saturation];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            boundary_c->SetTimedependentBCForcingFunction(bc_item.second);
            fAlphaSaturationMesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    fAlphaSaturationMesh->SetDimModel(dim);
    fAlphaSaturationMesh->SetDefaultOrder(sorder);
    fAlphaSaturationMesh->SetAllCreateFunctionsDiscontinuous();
    fAlphaSaturationMesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshS_alpha.txt");
    fAlphaSaturationMesh->Print(out);
#endif
    
}

/** @brief Create a computational mesh L2 */
void TRMSpaceOdissey::CreateBetaTransportMesh()
{
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int saturation = 0;
    int sorder = fSOrder;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
    // Malha computacional
    fBetaSaturationMesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMPhaseTransport * mat = new TRMPhaseTransport(rock_id);
        fBetaSaturationMesh->InsertMaterialObject(mat);
        
        // Inserting volumetric materials
        int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
        int bc_id = 0;
        std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > bc_item;
        TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > bc;
        for (int j = 0; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[saturation];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            boundary_c->SetTimedependentBCForcingFunction(bc_item.second);
            fBetaSaturationMesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    fBetaSaturationMesh->SetDimModel(dim);
    fBetaSaturationMesh->SetDefaultOrder(sorder);
    fBetaSaturationMesh->SetAllCreateFunctionsDiscontinuous();
    fBetaSaturationMesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshS_beta.txt");
    fBetaSaturationMesh->Print(out);
#endif
    
//    if (fBetaSaturationMesh) {
//        DebugStop();
//    }
//    fBetaSaturationMesh = new TPZCompMesh(fGeoMesh);
//    fBetaSaturationMesh->SetDimModel(3);
//    fBetaSaturationMesh->SetDefaultOrder(0);
//    fBetaSaturationMesh->SetAllCreateFunctionsDiscontinuous();
//    
//    TRMPhaseTransport *mat = new TRMPhaseTransport(_ReservMatId);
//    fBetaSaturationMesh->InsertMaterialObject(mat);
//    
//    TRMPhaseInterfaceTransport *matint = new TRMPhaseInterfaceTransport(_ReservoirInterface);
//    fGeoMesh->AddInterfaceMaterial(_ReservMatId, _ReservMatId,_ReservoirInterface);
//    
//    // WE NEED TO ADD THE BOUNDARY CONDITION MATERIALS
//    DebugStop();
//    
//    fBetaSaturationMesh->ApproxSpace().CreateInterfaces(fBetaSaturationMesh);
    
}

/** @brief Create a multiphysics computational mesh L2 */
void TRMSpaceOdissey::CreateTransportMesh(){
    // Second option put all the transpor meshes inside a multiphysics mesh

    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int saturation = 0;
    int sorder = fSOrder;
    int interface_id = fSimulationData->InterfacesMatId();
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > bc_item;
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > bc;
    
    // Malha computacional
    fTransportMesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMPhaseTransport * mat = new TRMPhaseTransport(rock_id);
        mat->SetSimulationData(fSimulationData);
        fTransportMesh->InsertMaterialObject(mat);
        
        TRMPhaseInterfaceTransport * matint = new TRMPhaseInterfaceTransport(interface_id);
        matint->SetSimulationData(fSimulationData);
        fTransportMesh->InsertMaterialObject(matint);
        fGeoMesh->AddInterfaceMaterial(rock_id, rock_id,interface_id);

        
        // Inserting volumetric materials
        int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
        int bc_id = 0;

        for (int j = 0; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[saturation];
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> * boundary_bc = new TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>;
            boundary_bc->SetMaterial(mat);
            boundary_bc->SetId(bc_id);
            boundary_bc->SetType(bc_item.first);
            TPZMaterial * material_bc = dynamic_cast<TPZMaterial * >(boundary_bc);
            material_bc->SetTimedependentBCForcingFunction(bc_item.second);
            fTransportMesh->InsertMaterialObject(boundary_bc);
        }
        
    }

    fTransportMesh->SetDimModel(dim);
    fTransportMesh->SetDefaultOrder(sorder);
    fTransportMesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fTransportMesh->ApproxSpace().CreateWithMemory(true);// Force the creating of interfaces with memory.
    fTransportMesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector;
    
    if(this->SimulationData()->IsTwoPhaseQ()){
        
        meshvector.Resize(1);
        meshvector[0] = fAlphaSaturationMesh.operator->();
        
        // Transferindo para a multifisica
        TPZBuildMultiphysicsMesh::AddElements(meshvector, fTransportMesh.operator->());
        TPZBuildMultiphysicsMesh::AddConnects(meshvector, fTransportMesh.operator->());
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fTransportMesh.operator->());
        
    }
    
    if(this->SimulationData()->IsThreePhaseQ()){

        meshvector.Resize(2);
        meshvector[0] = fAlphaSaturationMesh.operator->();
        meshvector[1] = fBetaSaturationMesh.operator->();
        
        // Transferindo para a multifisica
        TPZBuildMultiphysicsMesh::AddElements(meshvector, fTransportMesh.operator->());
        TPZBuildMultiphysicsMesh::AddConnects(meshvector, fTransportMesh.operator->());
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fTransportMesh.operator->());
    }

    
    fGeoMesh->ResetReference();
    fTransportMesh->LoadReferences();
    long nel = fTransportMesh->ElementVec().NElements();
    // Creation of interface elements
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = fTransportMesh->ElementVec()[el];
        if(!compEl) continue;
        TPZGeoEl * gel = compEl->Reference();
        if(!gel) {continue;}
        if(gel->HasSubElement()) {continue;}
        int index = compEl ->Index();
        if(compEl->Dimension() == fTransportMesh->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fTransportMesh->ElementVec()[index]);
            if(!InterpEl) {
                continue;
            }
            InterpEl->CreateInterfaces();
        }
    }

    fTransportMesh->CleanUpUnconnectedNodes();
    fTransportMesh->AdjustBoundaryElements();
    fTransportMesh->AutoBuild();
    
    nel = fTransportMesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fTransportMesh->Element(el);
        
        if(cel->Dimension() != 3){
            continue;
        }
        
        TPZMultiphysicsElement  * mf_cel = dynamic_cast<TPZMultiphysicsElement *>(cel);
//        TPZCompElWithMem<TPZMultiphysicsInterfaceElement> * bc_cel = dynamic_cast<TPZCompElWithMem<TPZMultiphysicsInterfaceElement> *>(mf_cel);
//        TPZBndCond * bc = dynamic_cast<TPZBndCond *>(face_cel);
        if (!mf_cel) {
            
//            if(!bc){
////                bc->PrepareIntPtIndices();
//            }
            
            continue;
        }
        
        mf_cel->InitializeIntegrationRule();
        mf_cel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CmeshTransport.txt");
    fTransportMesh->Print(out);
#endif
    
    
}


void TRMSpaceOdissey::PrintGeometry()
{
    //  Print Geometrical Base Mesh
    std::ofstream planefile("GeometricMesh.txt");
    fGeoMesh->Print(planefile);
    std::ofstream file("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh,file, true);
}

/** @brief Create a reservoir-box geometry */
void TRMSpaceOdissey::CreateGeometricGIDMesh(std::string &grid){
    
    TPZReadGIDGrid GeometryInfo;
    REAL s = 1.0;
    GeometryInfo.SetfDimensionlessL(s);
    fGeoMesh = GeometryInfo.GeometricGIDMesh(grid);
    fGeoMesh->SetDimension(3);
    
}

/** @brief Create a reservoir-box geometry */
void TRMSpaceOdissey::CreateGeometricBoxMesh(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy, TPZManVector<REAL,2> dz){
    
    REAL t=0.0;
    REAL dt;
    int n;
    

    int rock =  this->SimulationData()->RawData()->fOmegaIds[0];
    
    int bc_W =  this->SimulationData()->RawData()->fGammaIds[0];
    int bc_E =  this->SimulationData()->RawData()->fGammaIds[1];
    int bc_S =  this->SimulationData()->RawData()->fGammaIds[2];
    int bc_N =  this->SimulationData()->RawData()->fGammaIds[3];
    int bc_B =  this->SimulationData()->RawData()->fGammaIds[4];
    int bc_T =  this->SimulationData()->RawData()->fGammaIds[5];
    
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
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,rock,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom0D(GeoMesh0D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncX = new TPZDummyFunction<STATE>(ParametricfunctionX);
    CreateGridFrom0D.SetParametricFunction(ParFuncX);
    CreateGridFrom0D.SetFrontBackMatId(bc_W,bc_E);
    
    dt = dx[0];
    n = int(dx[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh1D = CreateGridFrom0D.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom1D(GeoMesh1D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncY = new TPZDummyFunction<STATE>(ParametricfunctionY);
    CreateGridFrom1D.SetParametricFunction(ParFuncY);
    CreateGridFrom1D.SetFrontBackMatId(bc_S,bc_N);
    
    dt = dy[0];
    n = int(dy[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh2D = CreateGridFrom1D.ComputeExtrusion(t, dt, n);
        
    TPZHierarquicalGrid CreateGridFrom2D(GeoMesh2D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncZ = new TPZDummyFunction<STATE>(ParametricfunctionZ);
    CreateGridFrom2D.SetParametricFunction(ParFuncZ);
    CreateGridFrom2D.SetFrontBackMatId(bc_B,bc_T);
    
    dt = dz[0];
    n = int(dz[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fGeoMesh = CreateGridFrom2D.ComputeExtrusion(t, dt, n);

    const std::string name("Reservoir box");
    fGeoMesh->SetName(name);
    
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


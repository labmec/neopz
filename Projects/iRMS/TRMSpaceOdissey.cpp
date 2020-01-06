
;//
//  TRMSpaceOdissey.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpaceOdissey.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
    f.Resize(1,0.);
    f[0] = pt[2];
    return;
}

static void CreateExampleRawData(TRMRawData &data)
{
    data.fLw = 500.;
    data.fHasLiner = false;
    data.fHasCasing = false;
    
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
    fIsHexaDominatedQ           = false;
    fIsTetraDominatedQ          = false;
    fIsPrismDominatedQ          = false;
    fGeoMesh                    = NULL;
    fSimulationData             = NULL;
    fH1Cmesh                    = NULL;
    fBiotCmesh                  = NULL;
    fGP_BiotCmesh               = NULL;
    fGeoPressureCmesh           = NULL;
    fFluxCmesh                  = NULL;
    fPressureCmesh              = NULL;
    fAlphaSaturationMesh        = NULL;
    fBetaSaturationMesh         = NULL;
    fGalerkinProjectionsCmesh   = NULL;
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
    if(fBiotCmesh)                  fBiotCmesh->CleanUp();
    if(fGP_BiotCmesh)               fGP_BiotCmesh->CleanUp();
    if(fGeoPressureCmesh)           fGeoPressureCmesh->CleanUp();
    if(fFluxCmesh)                  fFluxCmesh->CleanUp();
    if(fPressureCmesh)              fPressureCmesh->CleanUp();
    if(fAlphaSaturationMesh)        fAlphaSaturationMesh->CleanUp();
    if(fBetaSaturationMesh)         fBetaSaturationMesh->CleanUp();
    if(fGalerkinProjectionsCmesh)   fGalerkinProjectionsCmesh->CleanUp();
    if(fGeoMechanicsCmesh)          fGeoMechanicsCmesh->CleanUp();
    if(fTransportMesh)              fTransportMesh->CleanUp();
    if(fMixedFluxPressureCmesh)     fMixedFluxPressureCmesh->CleanUp();
    if(fPressureSaturationCmesh)    fPressureSaturationCmesh->CleanUp();
    if(fMonolithicMultiphaseCmesh)  fMonolithicMultiphaseCmesh->CleanUp();
    
}

/** @brief Create a Biot H1 computational mesh */
void TRMSpaceOdissey::CreateBiotCmesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    int Sigma_or_displacement = 0;
    int uorder = fUOrder;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fBiotCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        
        TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(rock_id,dim);
        fBiotCmesh->InsertMaterialObject(mat);
        
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 0;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;
        
        for (int j = initial_bc; j < n_boundauries; j++) {
            
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[Sigma_or_displacement];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fBiotCmesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    
    // Sideburden
    int side_burden_rock = 14;
    int bc_W = 20;
    int bc_N = 19;
    int bc_E = 18;
    int bc_S = 17;
    int bc_T = 16;
    int bc_B = 15;
    
    if (dim == 2) {
        side_burden_rock = 12;
        bc_W = 14;
        bc_N = 15;
        bc_E = 16;
        bc_S = 13;
        bc_T = 1000;
        bc_B = 1000;
    }

    
    TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(side_burden_rock,dim);
    fBiotCmesh->InsertMaterialObject(mat);

    TPZMaterial * W_bndc = mat->CreateBC(mat, bc_W, Sigma_or_displacement, val1, val2);
    fBiotCmesh->InsertMaterialObject(W_bndc);
    
    TPZMaterial * N_bndc = mat->CreateBC(mat, bc_N, Sigma_or_displacement, val1, val2);
    fBiotCmesh->InsertMaterialObject(N_bndc);
    
    TPZMaterial * E_bndc = mat->CreateBC(mat, bc_E, Sigma_or_displacement, val1, val2);
    fBiotCmesh->InsertMaterialObject(E_bndc);
    
    TPZMaterial * S_bndc = mat->CreateBC(mat, bc_S, Sigma_or_displacement, val1, val2);
    fBiotCmesh->InsertMaterialObject(S_bndc);
    
    TPZMaterial * T_bndc = mat->CreateBC(mat, bc_T, Sigma_or_displacement, val1, val2);
    fBiotCmesh->InsertMaterialObject(T_bndc);
    
    TPZMaterial * B_bndc = mat->CreateBC(mat, bc_B, Sigma_or_displacement, val1, val2);
    fBiotCmesh->InsertMaterialObject(B_bndc);
    
    fBiotCmesh->SetDimModel(dim);
    fBiotCmesh->SetDefaultOrder(uorder);
    fBiotCmesh->SetAllCreateFunctionsContinuous();
    fBiotCmesh->AutoBuild();
    
    
#ifdef PZDEBUG
    std::ofstream out("CmeshBiot.txt");
    fBiotCmesh->Print(out);
#endif
    
}

/** @brief Create a Biot poroelastic elastic space for RB generation computational mesh */
void TRMSpaceOdissey::CreateGeoPressureBiotCmesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    int Sigma_or_displacement = 0;
    int uorder = fUOrder;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fGP_BiotCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        
        TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(rock_id,dim);
        fGP_BiotCmesh->InsertMaterialObject(mat);
        
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 0;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;
        
        for (int j = initial_bc; j < n_boundauries; j++) {
            
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[Sigma_or_displacement];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fGP_BiotCmesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    
    // Sideburden
    int side_burden_rock = 14;
    int bc_W = 20;
    int bc_N = 19;
    int bc_E = 18;
    int bc_S = 17;
    int bc_T = 16;
    int bc_B = 15;
    
    if (dim == 2) {
        side_burden_rock = 12;
        bc_W = 14;
        bc_N = 15;
        bc_E = 16;
        bc_S = 13;
        bc_T = 1000;
        bc_B = 1000;
    }
    
    
    TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(side_burden_rock,dim);
    fGP_BiotCmesh->InsertMaterialObject(mat);
    
    TPZMaterial * W_bndc = mat->CreateBC(mat, bc_W, Sigma_or_displacement, val1, val2);
    fGP_BiotCmesh->InsertMaterialObject(W_bndc);
    
    TPZMaterial * N_bndc = mat->CreateBC(mat, bc_N, Sigma_or_displacement, val1, val2);
    fGP_BiotCmesh->InsertMaterialObject(N_bndc);
    
    TPZMaterial * E_bndc = mat->CreateBC(mat, bc_E, Sigma_or_displacement, val1, val2);
    fGP_BiotCmesh->InsertMaterialObject(E_bndc);
    
    TPZMaterial * S_bndc = mat->CreateBC(mat, bc_S, Sigma_or_displacement, val1, val2);
    fGP_BiotCmesh->InsertMaterialObject(S_bndc);
    
    TPZMaterial * T_bndc = mat->CreateBC(mat, bc_T, Sigma_or_displacement, val1, val2);
    fGP_BiotCmesh->InsertMaterialObject(T_bndc);
    
    TPZMaterial * B_bndc = mat->CreateBC(mat, bc_B, Sigma_or_displacement, val1, val2);
    fGP_BiotCmesh->InsertMaterialObject(B_bndc);
    
    fGP_BiotCmesh->SetDimModel(dim);
    fGP_BiotCmesh->SetDefaultOrder(uorder);
    fGP_BiotCmesh->SetAllCreateFunctionsContinuous();
    fGP_BiotCmesh->AutoBuild();
    
    
#ifdef PZDEBUG
    std::ofstream out("CmeshGProjections.txt");
    fGP_BiotCmesh->Print(out);
#endif
    
}

/** @brief Create a Biot RB computational mesh */
TPZCompMesh * TRMSpaceOdissey::CreateRB_BiotCmesh(){
    
    if(!fGeoMesh || !fGP_BiotCmesh)
    {
        std::cout<< "Geometric mesh or GP mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    int Sigma_or_displacement = 0;
    int uorder = 0;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    TPZCompMeshReferred * cmesh_rb = new TPZCompMeshReferred(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        
        TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(rock_id,dim);
        if (fSimulationData->ReducedBasisResolution().first) {
            mat->Set_NStateVariables(1);
        }
        cmesh_rb->InsertMaterialObject(mat);
        
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 0;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;
        
        for (int j = initial_bc; j < n_boundauries; j++) {
            
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[Sigma_or_displacement];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            cmesh_rb->InsertMaterialObject(boundary_c);
        }
        
    }
    
    
    // Sideburden
    int side_burden_rock = 14;
    int bc_W = 20;
    int bc_N = 19;
    int bc_E = 18;
    int bc_S = 17;
    int bc_T = 16;
    int bc_B = 15;
    
    if (dim == 2) {
        side_burden_rock = 12;
        bc_W = 14;
        bc_N = 15;
        bc_E = 16;
        bc_S = 13;
        bc_T = 1000;
        bc_B = 1000;
    }
    
    
    TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(side_burden_rock,dim);
    if (fSimulationData->ReducedBasisResolution().first) {
        mat->Set_NStateVariables(1);
    }
    cmesh_rb->InsertMaterialObject(mat);
    
    TPZMaterial * W_bndc = mat->CreateBC(mat, bc_W, Sigma_or_displacement, val1, val2);
    cmesh_rb->InsertMaterialObject(W_bndc);
    
    TPZMaterial * N_bndc = mat->CreateBC(mat, bc_N, Sigma_or_displacement, val1, val2);
    cmesh_rb->InsertMaterialObject(N_bndc);
    
    TPZMaterial * E_bndc = mat->CreateBC(mat, bc_E, Sigma_or_displacement, val1, val2);
    cmesh_rb->InsertMaterialObject(E_bndc);
    
    TPZMaterial * S_bndc = mat->CreateBC(mat, bc_S, Sigma_or_displacement, val1, val2);
    cmesh_rb->InsertMaterialObject(S_bndc);
    
    TPZMaterial * T_bndc = mat->CreateBC(mat, bc_T, Sigma_or_displacement, val1, val2);
    cmesh_rb->InsertMaterialObject(T_bndc);
    
    TPZMaterial * B_bndc = mat->CreateBC(mat, bc_B, Sigma_or_displacement, val1, val2);
    cmesh_rb->InsertMaterialObject(B_bndc);
    
    // Setting RB approximation space
    cmesh_rb->SetDimModel(dim);
    int numsol = fGP_BiotCmesh->Solution().Cols();
    cmesh_rb->AllocateNewConnect(numsol, 1, uorder);
    TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmesh_rb);
    cmesh_rb->AutoBuild();
    
    cmesh_rb->AdjustBoundaryElements();
    cmesh_rb->CleanUpUnconnectedNodes();
    cmesh_rb->LoadReferred(fGP_BiotCmesh);
    
#ifdef PZDEBUG
    std::ofstream out("CmeshRB_Biot.txt");
    cmesh_rb->Print(out);
#endif
    
    return cmesh_rb;
    
}


/** @brief Create a constant pressure computational mesh */
void TRMSpaceOdissey::CreateGeoPressuresCmesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    int porder = 0;
    
    // Malha computacional
    fGeoPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMPoroelasticModes * mat = new TRMPoroelasticModes(rock_id,dim,1);
        fGeoPressureCmesh->InsertMaterialObject(mat);
        
    }
    
    // Setando L2
    fGeoPressureCmesh->SetDimModel(dim);
    fGeoPressureCmesh->SetDefaultOrder(porder);
    
    fGeoPressureCmesh->SetAllCreateFunctionsDiscontinuous();
    fGeoPressureCmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshGeoPress.txt");
    fGeoPressureCmesh->Print(out);
#endif
    
}

/** @brief Create a Hdiv computational mesh Hdiv */
void TRMSpaceOdissey::CreateFluxCmesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    int flux_or_pressure = 0;
    int qorder = fPOrder;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fFluxCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        
        // Inserting volumetric materials
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id,dim);
        fFluxCmesh->InsertMaterialObject(mat);
        
        TRMMixedDarcy * mat_skeleton = new TRMMixedDarcy(fSimulationData->Skeleton_material_Id(),dim-1);
        fFluxCmesh->InsertMaterialObject(mat_skeleton); // @omar::  skeleton material inserted
        

        if (rock_id == 5) { // Reservoir
            n_boundauries = 6;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;

        // Inserting boundary materials associated to the current mat
        for (int j = initial_bc; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[flux_or_pressure];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fFluxCmesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    fFluxCmesh->SetDimModel(dim);
    fFluxCmesh->SetDefaultOrder(qorder);
    fFluxCmesh->SetAllCreateFunctionsHDiv();
    fFluxCmesh->AutoBuild();
    
    bool IncreaseAccQ = fSimulationData->IsAdataptedQ();
    int wellbore_order = 2;
    
    if(IncreaseAccQ){
        
        int n_el = fFluxCmesh->NElements();
        for (int icel = 0; icel < n_el; icel++) {
            TPZCompEl * cel = fFluxCmesh->Element(icel);
            if(!cel) continue;
            
            TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
            if(!sp) continue;
            if (sp->Reference()->Dimension() == dim - 1) {
                continue;
            }
            if (sp->Reference()->MaterialId() == 6) {
                sp->PRefine(wellbore_order);
            }
            if (sp->Reference()->MaterialId() == 7) {
                sp->PRefine(wellbore_order);
            }
            
        }
        
        fFluxCmesh->AdjustBoundaryElements();
        fFluxCmesh->CleanUpUnconnectedNodes();
        fFluxCmesh->ExpandSolution();
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
    
    int dim = fGeoMesh->Dimension();
    int porder = fPOrder;
    
    // Malha computacional
    fPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id,dim);
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
    
    bool IncreaseAccQ = fSimulationData->IsAdataptedQ();
    int wellbore_order = 2;
    
    if(IncreaseAccQ){

        int n_el = fPressureCmesh->NElements();
        for (int icel = 0; icel < n_el; icel++) {
            TPZCompEl * cel = fPressureCmesh->Element(icel);
            if(!cel) continue;
            
            TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
            if(!sp) continue;
            if (sp->Reference()->Dimension() == dim - 1) {
                continue;
            }
            if (sp->Reference()->MaterialId() == 6) {
                sp->PRefine(wellbore_order);
            }
            if (sp->Reference()->MaterialId() == 7) {
                sp->PRefine(wellbore_order);
            }

        }

        fPressureCmesh->AdjustBoundaryElements();
        fPressureCmesh->CleanUpUnconnectedNodes();
        fPressureCmesh->ExpandSolution();
    }

    int ncon = fPressureCmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = fPressureCmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::ofstream out("CmeshPress.txt");
    fPressureCmesh->Print(out);
#endif
    
}

/// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
void TRMSpaceOdissey::AdjustFluxPolynomialOrders(int n_acc_terms)
{
    int dim = fFluxCmesh->Dimension();
    /// loop over all the elements
    long nel = fFluxCmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fFluxCmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        // compute the maxorder
        int maxorder = -1;
        int ncon = intel->NConnects();
        for (int i=0; i<ncon-1; i++) {
            int conorder = intel->Connect(i).Order();
            maxorder = maxorder < conorder ? conorder : maxorder;
        }
        int nsides = gel->NSides();
        int nconside = intel->NSideConnects(nsides-1);
        long cindex = intel->SideConnectIndex(nconside-1, nsides-1);
        TPZConnect &c = fFluxCmesh->ConnectVec()[cindex];
        if (c.NElConnected() != 1) {
            DebugStop();
        }
        if (c.Order()+n_acc_terms != maxorder) {
            intel->SetSideOrder(nsides-1, maxorder+n_acc_terms);
        }
    }
    fFluxCmesh->ExpandSolution();
}

void TRMSpaceOdissey::SetPressureOrders()
{
    // build a vector with the required order of each element in the pressuremesh
    // if an element of the mesh dimension of the fluxmesh does not have a corresponding element in the pressuremesh DebugStop is called
    int meshdim = fPressureCmesh->Dimension();
    fPressureCmesh->Reference()->ResetReference();
    fPressureCmesh->LoadReferences();
    TPZManVector<long> pressorder(fPressureCmesh->NElements(),-1);
    long nel = fFluxCmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fFluxCmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != meshdim) {
            continue;
        }
        int nsides = gel->NSides();
        long cindex = intel->SideConnectIndex(0, nsides-1);
        TPZConnect &c = fFluxCmesh->ConnectVec()[cindex];
        int order = c.Order();
        TPZCompEl *pressureel = gel->Reference();
        TPZInterpolatedElement *pintel = dynamic_cast<TPZInterpolatedElement *>(pressureel);
        if (!pintel) {
            DebugStop();
        }
        pressorder[pintel->Index()] = order;
    }
    fPressureCmesh->Reference()->ResetReference();
    nel = pressorder.size();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureCmesh->Element(el);
        TPZInterpolatedElement *pintel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!pintel) {
            continue;
        }
        if (pressorder[el] == -1) {
            continue;
        }
        pintel->PRefine(pressorder[el]);
    }
    
    fPressureCmesh->ExpandSolution();
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

/** @brief Build Geomechanic mesh form the current H1/RB mesh */
void TRMSpaceOdissey::BuildGeomechanic_Mesh(){
    
    this->CreateBiotCmesh();
    this->CreateGeoMechanicMesh();
}

/** @brief Create a RB computational mesh for Maurice Biot Linear Poroelasticity */
void TRMSpaceOdissey::BuildRBGeomechanic_Mesh(){
    
    
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    // offline split
    fSimulationData->ReducedBasisResolution().second.first = true;    // offline split
    this->CreateGeoPressureBiotCmesh();
    this->CreateGeoPressuresCmesh();
    this->CreateGeoModesCmesh();
    this->RB_Generator();
    fBiotCmesh = this->CreateRB_BiotCmesh();
    this->CreateGeoMechanicMesh();
    fSimulationData->ReducedBasisResolution().second.first = false;    // offline split
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    std::cout  << "iRMS::RB:: Time for offline process " << (t2-t1) << std::endl;
#endif

}

/** @brief Create The reduced basis */
void TRMSpaceOdissey::RB_Generator(){

    TRMBuildTransfers * transfer = new TRMBuildTransfers;
    transfer->SetSimulationData(fSimulationData);    
    bool mustOptimizeBandwidth_e = true;
    int numofThreads_e = 6;
    TRMGeomechanicAnalysis  * RB_generator      = new TRMGeomechanicAnalysis;

    RB_generator->Meshvec().Resize(2);
    // Analysis for elliptic part
    RB_generator->Meshvec()[0] = this->GP_BiotCmesh();
    RB_generator->Meshvec()[1] = this->GeoPressureCMesh();
    RB_generator->SetCompMesh(this->GalerkinProjectionsCMesh(),mustOptimizeBandwidth_e);
    
    
    if(fSimulationData->UsePardisoQ()){
        TPZSymetricSpStructMatrix strmat_e(this->GalerkinProjectionsCMesh());
        
        TPZStepSolver<STATE> step_e;
        step_e.SetDirect(ELDLt);
        strmat_e.SetNumThreads(numofThreads_e);
        RB_generator->SetStructuralMatrix(strmat_e);
        RB_generator->SetSolver(step_e);
        RB_generator->AdjustVectors();
        RB_generator->SetSimulationData(fSimulationData);
    }
    else{
        
//        TPZSymetricSpStructMatrix strmat_e(this->GalerkinProjectionsCMesh());
        
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat_e(this->GalerkinProjectionsCMesh());
        strmat_e.SetDecomposeType(ELDLt);
        
//        TPZSkylineStructMatrix strmat_e(this->GalerkinProjectionsCMesh());
        
        TPZStepSolver<STATE> step_e;
        step_e.SetDirect(ELDLt);
        strmat_e.SetNumThreads(numofThreads_e);
        RB_generator->SetStructuralMatrix(strmat_e);
        RB_generator->SetSolver(step_e);
        RB_generator->AdjustVectors();
        RB_generator->SetSimulationData(fSimulationData);
    }
        
    

    
    transfer->spatial_props_To_elliptic(RB_generator->Mesh()); // load properties inside memory

    bool state = RB_generator->SimulationData()->IsInitialStateQ();
    RB_generator->SimulationData()->SetInitialStateQ(false);
    RB_generator->SimulationData()->SetCurrentStateQ(true);
    RB_generator->Assemble();
    
    // Setting up the empirical interpolation based on unitary pressures
    
    REAL unit_p = 1.0e6;
    TPZStack<TPZVec<long> > cts_pressures;
    
    int n_blocks = DrawingPressureBlocks(RB_generator->Meshvec()[1], cts_pressures,5); // reservoir partition
    n_blocks += DrawingPressureBlocksByID(RB_generator->Meshvec()[1], cts_pressures); // wellbore region
    
    int ndof_elastic = RB_generator->Meshvec()[0]->NEquations();
    TPZFMatrix<REAL> galerkin_projts(ndof_elastic,n_blocks);
    galerkin_projts.Zero();
    
    int progress = (n_blocks/10) + 1;
    REAL percent = -10.0;
    fSimulationData->Set_m_RB_functions(n_blocks);
    std::cout<< "RB:: number of geomodes = " << n_blocks << std::endl;
    std::cout << "RB:: ndof rb generator = " << fGP_BiotCmesh->Solution().Rows() << std::endl;
    for (int ip = 0; ip < n_blocks; ip++) {
        
        
        RB_generator->Meshvec()[0]->Solution().Zero();
        RB_generator->Meshvec()[1]->Solution().Zero();
        
        for (int jp = 0; jp < cts_pressures[ip].size(); jp++) {
            RB_generator->Meshvec()[1]->Solution()(cts_pressures[ip][jp],0) = unit_p;
        }
        
        TPZBuildMultiphysicsMesh::TransferFromMeshes(RB_generator->Meshvec(), RB_generator->Mesh());
        RB_generator->X_n() = RB_generator->Mesh()->Solution();
        RB_generator->AssembleResidual();
        RB_generator->Rhs() *= -1.0;
        RB_generator->Solve();
        RB_generator->Solution() += RB_generator->X_n();
        RB_generator->LoadSolution();
        RB_generator->PostProcessStep();
        
        if(ip%progress == 0){
            percent += 10.0;
            std::cout << " Progress on offline stage " << percent  << " % " <<std::endl;
        }
        
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(RB_generator->Meshvec(), RB_generator->Mesh());
        galerkin_projts.AddSub(0, ip, RB_generator->Meshvec()[0]->Solution());
    }
    
    RB_generator->Meshvec()[0]->LoadSolution(galerkin_projts);
    
//    fGP_BiotCmesh = RB_generator->Meshvec()[0];
    
#ifdef PZDEBUG
    std::ofstream out("CmeshGeoModesLoaded.txt");
    fGP_BiotCmesh->Print(out);
#endif
    
    RB_generator->SimulationData()->SetInitialStateQ(state);

}

/** @brief Select regions for pressure fields */
int TRMSpaceOdissey::DrawingPressureBlocks(TPZCompMesh * cmesh, TPZStack<TPZVec<long> > & constant_pressures, int target_id){
    
#ifdef PZDEBUG
    if(!cmesh){
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = cmesh->Reference();
    
#ifdef PZDEBUG
    if(!geometry){
        DebugStop();
    }
#endif
    
    int dim = geometry->Dimension();
    cmesh->LoadReferences();
    
    TPZStack<REAL> min_x;
    TPZStack<REAL> max_x;
    bool Outline_is_doneQ = DrawingGeometryOutline(min_x, max_x,5);
    
    if(!Outline_is_doneQ){
        DebugStop();
    }
    
#ifdef PZDEBUG
    if(fSimulationData->ReducedBasisResolution().second.second.size() == 0){
        DebugStop();
    }
#endif
    
    int ni,nj,nk;
    
    if (dim == 2) {
        ni = fSimulationData->ReducedBasisResolution().second.second[0];
        nj = fSimulationData->ReducedBasisResolution().second.second[1];
        nk = 1;
    }
    else{
        ni = fSimulationData->ReducedBasisResolution().second.second[0];
        nj = fSimulationData->ReducedBasisResolution().second.second[1];
        nk = fSimulationData->ReducedBasisResolution().second.second[2];
    }
 
    
    TPZManVector<REAL,3> x0(3,0.0);
    
    TPZStack<REAL> rule_x;
    REAL dx = (max_x[0]-min_x[0])/ni;
    REAL xv;
    for (int i = 0; i < ni + 1; i++) {
        xv = REAL(i)*dx + min_x[0];
        rule_x.Push(xv);
    }
    
    TPZStack<REAL> rule_y;
    REAL dy = (max_x[1]-min_x[1])/nj;
    REAL yv;
    for (int j = 0; j < nj + 1; j++) {
        yv = REAL(j)*dy + min_x[1];
        rule_y.Push(yv);
    }
    
    TPZStack<REAL> rule_z;
    REAL dz = (max_x[2]-min_x[2])/nk;
    REAL zv;
    for (int k = 0; k < nk + 1; k++) {
        zv = REAL(k)*dz + min_x[2];
        rule_z.Push(zv);
    }
    
    if (rule_z[0] == rule_z[1] && dim == 2) {
        std::cout << "RB:: Drawing Pressure Blocks for 2D geometry " <<std::endl;
    }
    
    if (dim == 3) {
        std::cout << "RB:: Drawing Pressure Blocks for 3D geometry " <<std::endl;
    }
    
    // counting volumetric elements
    int nel = geometry->NElements();
    int n_volumes = 0;
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if(!gel){
            DebugStop();
        }
#endif
        
        if(gel->HasSubElement() || gel->Dimension() != dim){
            continue;
        }
        
//        bool target_regionQ = gel->MaterialId() == target_id;
//        bool target_regionQ = gel->MaterialId() == 5 || gel->MaterialId()== 6 || gel->MaterialId() == 7;
        bool target_regionQ = gel->MaterialId() == 12 || gel->MaterialId()== 14 || gel->MaterialId()== 6 || gel->MaterialId() == 7;
        if(target_regionQ){
            continue;
        }
        
        n_volumes++;
    }
    
    std::cout << "RB:: Number of volumetric elements =  " << n_volumes << std::endl;
    
    // goup elements by Cartesian Grid
    TPZStack< TPZStack<long> > geo_groups;
    
    if (dim == 2) {
        // Check if the element belong to the cartesian block
        for (int i = 0; i < rule_x.size() - 1; i++) {
            for (int j = 0; j < rule_y.size() - 1; j++) {
                for (int k = 0; k < rule_z.size() - 1; k++) {
                    
                    // for each i,j,k box
                    TPZStack<long> box_group;
                    
                    TPZManVector<REAL,3> x_c(3,0.0);
                    TPZManVector<REAL,3> par_c(dim,0.0);
                    int nel = geometry->NElements();
                    for (int iel = 0; iel < nel; iel++) {
                        TPZGeoEl * gel = geometry->Element(iel);
                        
#ifdef PZDEBUG
                        if(!gel){
                            DebugStop();
                        }
#endif
                        
                        if(gel->Level() != 0 || gel->Dimension() != dim){
                            continue;
                        }
                        
                        bool target_regionQ = gel->MaterialId() == 12 || gel->MaterialId()== 14 || gel->MaterialId()== 6 || gel->MaterialId() == 7;
                        if(target_regionQ){
                            continue;
                        }
                        
                        gel->CenterPoint(gel->NSides()-1, par_c);
                        gel->X(par_c, x_c);
                        
                        // It is inside x
                        if(rule_x[i] <= x_c[0] && x_c[0] < rule_x[i+1]){
                            // It is inside y
                            if(rule_y[j] <= x_c[1] && x_c[1] < rule_y[j+1]){
                                // It is inside y
                                box_group.Push(gel->Index());
                            }
                        }
                        
                    }
                    
                    if (box_group.size() == 0) {
                        continue;
                    }
                    
                    geo_groups.Push(box_group);
                    //            std::cout << " group of geo elements with indexes =  " << box_group << std::endl;
                }
                
            }
        }
    }
    else{
        // Check if the element belong to the cartesian block
        for (int i = 0; i < rule_x.size() - 1; i++) {
            for (int j = 0; j < rule_y.size() - 1; j++) {
                for (int k = 0; k < rule_z.size() - 1; k++) {
                    
                    // for each i,j,k box
                    TPZStack<long> box_group;
                    
                    TPZManVector<REAL,3> x_c(3,0.0);
                    TPZManVector<REAL,3> par_c(dim,0.0);
                    int nel = geometry->NElements();
                    for (int iel = 0; iel < nel; iel++) {
                        TPZGeoEl * gel = geometry->Element(iel);
                        
#ifdef PZDEBUG
                        if(!gel){
                            DebugStop();
                        }
#endif
                        
                        if(gel->Level() != 0 || gel->Dimension() != dim){
                            continue;
                        }
                        
                        bool target_regionQ = gel->MaterialId() == 12 || gel->MaterialId()== 14 || gel->MaterialId()== 6 || gel->MaterialId() == 7;
                        if(target_regionQ){
                            continue;
                        }
                        
                        gel->CenterPoint(gel->NSides()-1, par_c);
                        gel->X(par_c, x_c);
                        
                        // It is inside x
                        if(rule_x[i] <= x_c[0] && x_c[0] < rule_x[i+1]){
                            // It is inside y
                            if(rule_y[j] <= x_c[1] && x_c[1] < rule_y[j+1]){
                                // It is inside z

                                if(rule_z[k] <= x_c[2] && x_c[2] < rule_z[k+1]){
                                    
                                    box_group.Push(gel->Index());
                                }
                            }
                        }
                        
                    }
                    
                    if (box_group.size() == 0) {
                        continue;
                    }
                    
                    geo_groups.Push(box_group);
                }
                
            }
        }
    }
    
#ifdef PZDEBUG
    if(geo_groups.size() == 0){
        DebugStop();
    }
#endif
    
    // drained response mode
    TPZVec<long> dofs(0);
//    constant_pressures.Push(dofs);
    
    // Pick blocks dofs
    TPZVec<long> dof_indexes;
    TPZVec<long> igroup;
    int n_groups = geo_groups.size();
    
    // For all groups
    int group_n_vol = 0;
    for (int ig = 0; ig < n_groups; ig++) {
        igroup = geo_groups[ig];
        
        // For a selecteced group of elements
        TPZStack<long> dof_stack;
        for(int igel = 0; igel < igroup.size(); igel++){
            
            TPZGeoEl * gel = geometry->Element(igroup[igel]);
            
#ifdef PZDEBUG
            if(!gel || gel->Dimension() != dim){
                DebugStop();
            }
#endif
            
            bool target_regionQ = gel->MaterialId() == 12 || gel->MaterialId()== 14 || gel->MaterialId()== 6 || gel->MaterialId() == 7;
            if(target_regionQ){
                continue;
            }
            
            TPZVec<TPZGeoEl *> unrefined_sons;
            gel->GetHigherSubElements(unrefined_sons);
            int nsub = unrefined_sons.size();
            for (int isub = 0; isub < nsub; isub++) {
                TPZGeoEl * subgel = unrefined_sons[isub];
                
                if (subgel->Dimension() != dim) {
                    continue;
                }
                
                TPZCompEl *cel = subgel->Reference();
#ifdef PZDEBUG
                if(!cel){
                    DebugStop();
                }
#endif
                
                TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(cel);
#ifdef PZDEBUG
                if(!intel){
                    DebugStop();
                }
#endif
                group_n_vol++;
                ElementDofIndexes(intel, dof_indexes);
                for (int i = 0;  i < dof_indexes.size(); i++) {
                    dof_stack.Push(dof_indexes[i]);
                }
                
            }
            
            if(nsub == 0){
                TPZCompEl *cel = gel->Reference();
#ifdef PZDEBUG
                if(!cel){
                    DebugStop();
                }
#endif
                
                TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(cel);
#ifdef PZDEBUG
                if(!intel){
                    DebugStop();
                }
#endif
                
                group_n_vol++;
                ElementDofIndexes(intel, dof_indexes);
                for (int i = 0;  i < dof_indexes.size(); i++) {
                    dof_stack.Push(dof_indexes[i]);
                }
            }
        }
        
        TPZVec<long> dofs(dof_stack);
        constant_pressures.Push(dofs);
        
    }
    
    
    if(group_n_vol != n_volumes){
        std::cout << "RB:: Drawing Pressure Blocks left some elements out! " <<std::endl;
        //        DebugStop();
    }
    
    int n_pressure_blocks = constant_pressures.size();
    if(n_pressure_blocks == 0){
        DebugStop();
    }
    
    return n_pressure_blocks;
    
    
}

/** @brief Select regions for pressure fields */
int TRMSpaceOdissey::DrawingPressureBlocksByID(TPZCompMesh * cmesh, TPZStack<TPZVec<long> > & constant_pressures){
    
#ifdef PZDEBUG
    if(!cmesh){
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = cmesh->Reference();
    
#ifdef PZDEBUG
    if(!geometry){
        DebugStop();
    }
#endif
    
    int dim = geometry->Dimension();
    cmesh->LoadReferences();
    
    // goup elements by Cartesian Grid
    TPZStack< TPZStack<long> > geo_groups;
    
//    TPZStack<long> box_group_res;
    
    // counting volumetric elements
    int nel = geometry->NElements();
    int n_volumes = 0;
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if(!gel){
            DebugStop();
        }
#endif
        
        if(gel->HasSubElement() || gel->Dimension() != dim){
            continue;
        }
        
        bool skip_regionQ = gel->MaterialId() == 12 || gel->MaterialId()== 14;
        if(skip_regionQ){
            continue;
        }
        n_volumes++;
        
//        bool target_regionQ = gel->MaterialId() == 5;
//        if(target_regionQ){
//            box_group_res.Push(gel->Index());
//            continue;
//        }
    
        {
            TPZStack<long> box_group_wb;
            bool target_regionQ = gel->MaterialId() == 6 || gel->MaterialId() == 7;
            if(target_regionQ){
                box_group_wb.Push(gel->Index());
            }
            geo_groups.Push(box_group_wb);
        }
    }
//    geo_groups.Push(box_group_res);
    
    std::cout << "RB:: Number of volumetric elements =  " << n_volumes << std::endl;
    
#ifdef PZDEBUG
    if(geo_groups.size() == 0){
        DebugStop();
    }
#endif
    
    // drained response mode
    TPZVec<long> dofs(0);
    //    constant_pressures.Push(dofs);
    
    // Pick blocks dofs
    TPZVec<long> dof_indexes;
    TPZVec<long> igroup;
    int n_groups = geo_groups.size();
    
    // For all groups
    int group_n_vol = 0;
    for (int ig = 0; ig < n_groups; ig++) {
        igroup = geo_groups[ig];
        
        // For a selecteced group of elements
        TPZStack<long> dof_stack;
        for(int igel = 0; igel < igroup.size(); igel++){
            
            TPZGeoEl * gel = geometry->Element(igroup[igel]);
            
#ifdef PZDEBUG
            if(!gel || gel->Dimension() != dim){
                DebugStop();
            }
#endif
            
            bool target_regionQ = gel->MaterialId() == 12 || gel->MaterialId()== 14;
            if(target_regionQ){
                continue;
            }
            
            TPZVec<TPZGeoEl *> unrefined_sons;
            gel->GetHigherSubElements(unrefined_sons);
            int nsub = unrefined_sons.size();
            for (int isub = 0; isub < nsub; isub++) {
                TPZGeoEl * subgel = unrefined_sons[isub];
                
                if (subgel->Dimension() != dim) {
                    continue;
                }
                
                TPZCompEl *cel = subgel->Reference();
#ifdef PZDEBUG
                if(!cel){
                    DebugStop();
                }
#endif
                
                TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(cel);
#ifdef PZDEBUG
                if(!intel){
                    DebugStop();
                }
#endif
                group_n_vol++;
                ElementDofIndexes(intel, dof_indexes);
                for (int i = 0;  i < dof_indexes.size(); i++) {
                    dof_stack.Push(dof_indexes[i]);
                }
                
            }
            
            if(nsub == 0){
                TPZCompEl *cel = gel->Reference();
#ifdef PZDEBUG
                if(!cel){
                    DebugStop();
                }
#endif
                
                TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(cel);
#ifdef PZDEBUG
                if(!intel){
                    DebugStop();
                }
#endif
                
                group_n_vol++;
                ElementDofIndexes(intel, dof_indexes);
                for (int i = 0;  i < dof_indexes.size(); i++) {
                    dof_stack.Push(dof_indexes[i]);
                }
            }
        }
        
        TPZVec<long> dofs(dof_stack);
        constant_pressures.Push(dofs);
        
    }
    
    
    if(group_n_vol != n_volumes){
        std::cout << "RB:: Drawing Pressure Blocks left some elements out! " <<std::endl;
        //        DebugStop();
    }
    
    int n_pressure_blocks = constant_pressures.size();
    if(n_pressure_blocks == 0){
        DebugStop();
    }
    
    return group_n_vol;
    
}

/** @brief Compute reservoir outline */
bool TRMSpaceOdissey::DrawingGeometryOutline(TPZStack< REAL > & min_x, TPZStack< REAL > & max_x, int target_id){
    
#ifdef PZDEBUG
    if(!fGeoMesh){
        DebugStop();
    }
#endif
    
    long n_elemtens = fGeoMesh->NElements();
    int dim = fGeoMesh->Dimension();
    
    REAL min_xc = +1.0e12;
    REAL max_xc = -1.0e12;
    REAL min_yc = +1.0e12;
    REAL max_yc = -1.0e12;
    REAL min_zc = +1.0e12;
    REAL max_zc = -1.0e12;
    
    TPZManVector<REAL,3> co(3,0.0);
    for (long iel = 0; iel < n_elemtens; iel++) {
        
        TPZGeoEl * gel = fGeoMesh->Element(iel);
        
#ifdef PZDEBUG
        if(!gel){
            DebugStop();
        }
#endif
        
        if(gel->HasSubElement() || gel->Dimension() != dim){// apply RB process on lowest level
            continue;
        }
        
//        bool target_regionQ = gel->MaterialId() == target_id;                
//        bool target_regionQ = gel->MaterialId() == 5 || gel->MaterialId()== 6 || gel->MaterialId() == 7;
        bool target_regionQ = gel->MaterialId() == 12 || gel->MaterialId()== 14;
        if(target_regionQ){
            continue;
        }
        
        int n_nodes = gel->NNodes();
        
        for (long inode = 0; inode < n_nodes; inode++) {

            TPZGeoNode node = gel->Node(inode);
            node.GetCoordinates(co);

            // x limits
            if (min_xc > co[0]) {
            min_xc = co[0];
            }

            if (max_xc < co[0]) {
            max_xc = co[0];
            }

            // y limits
            if (min_yc > co[1]) {
            min_yc = co[1];
            }

            if (max_yc < co[1]) {
            max_yc = co[1];
            }

            // z limits
            if (min_zc > co[2]) {
            min_zc = co[2];
            }

            if (max_zc < co[2]) {
            max_zc = co[2];
            }
        }
    }
    
    min_x.Push(min_xc);
    min_x.Push(min_yc);
    min_x.Push(min_zc);
    
    max_x.Push(max_xc);
    max_x.Push(max_yc);
    max_x.Push(max_zc);
    
    return true;
    
}

void TRMSpaceOdissey::ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    int nconnect = intel->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = intel->Mesh()->Block().Position(seqnumber);
        int nshape = con.NShape();
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
        }
    }
    
    dof_indexes = index;
    return;
}


/** @brief Build MHM form the current hdvi mesh */
void TRMSpaceOdissey::BuildMixed_Mesh(){
    
    this->CreateFluxCmesh();
    this->CreatePressureCmesh();
    
    if (fSimulationData->IsEnhancedPressureQ()) {
        int n_acc_terms = 1;
        this->AdjustFluxPolynomialOrders(n_acc_terms);
        this->SetPressureOrders();
    }
    
    this->CreateMixedCmesh();
    
}

/** @brief Build MHM form the current hdvi mesh */
void TRMSpaceOdissey::BuildMHM_Mesh(){
    
    this->CreateFluxCmesh();
    this->CreatePressureCmesh();

    if (fSimulationData->IsEnhancedPressureQ()) {
        int n_acc_terms = 1;
        this->AdjustFluxPolynomialOrders(n_acc_terms);
        this->SetPressureOrders();
    }
    
    SeparateConnectsByNeighborhood();
    
    this->CreateMixedCmesh();
    std::cout << "ndof parabolic MHM = " << fMixedFluxPressureCmesh->Solution().Rows() << std::endl;
    
    
    this->CreateMixedCmeshMHM();
    this->BuildMacroElements(); // @omar:: require the destruction and construction of the substrutucture mhm mesh
#ifdef PZDEBUG
    std::ofstream out_mhm("CmeshMixedMHM.txt");
    this->MixedFluxPressureCmeshMHM()->Print(out_mhm);
#endif
    std::cout << "ndof parabolic MHM substructures = " << fMixedFluxPressureCmeshMHM->Solution().Rows() << std::endl;
    
//    this->UnwrapMacroElements();
    
#ifdef PZDEBUG
    std::ofstream out_unwrap("CmeshMixedMHMUnWrap.txt");
    this->MixedFluxPressureCmeshMHM()->Print(out_unwrap);
#endif
    
    
}

/** @brief Sparated connects by given selected skeleton ids */
void TRMSpaceOdissey::SeparateConnectsBySkeletonIds(TPZVec<long> skeleton_ids){
    DebugStop();
}

/** @brief Sparated connects by hdiv connect neighborhood */
void TRMSpaceOdissey::SeparateConnectsByNeighborhood(){
    
#ifdef PZDEBUG
    if(!fFluxCmesh){
        DebugStop();
    }
#endif
    
    TPZGeoMesh *gmesh = fFluxCmesh->Reference();
    gmesh->ResetReference();
    fFluxCmesh->LoadReferences();
    fFluxCmesh->ComputeNodElCon();
    long nel = fFluxCmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fFluxCmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (!gel || (gel->Dimension() != gmesh->Dimension() && gel->MaterialId() != fSimulationData->Skeleton_material_Id()) ) {
            continue;
        }
        int nc = cel->NConnects();
        for (int ic =0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if ((c.HasDependency() && c.NElConnected() == 2) || (gel->MaterialId() == fSimulationData->Skeleton_material_Id() && c.NElConnected() == 2)) // @omar:: Hdiv connects have this invariant characteristic
            {
                // duplicate the connect
                long cindex = fFluxCmesh->AllocateNewConnect(c);
                TPZConnect &newc = fFluxCmesh->ConnectVec()[cindex];
                newc = c;
                c.DecrementElConnected();
                newc.DecrementElConnected();
                cel->SetConnectIndex(ic, cindex);
            }
        }
    }
    fFluxCmesh->ExpandSolution();
}

/** @brief Build MHM form the current hdvi mesh */
void TRMSpaceOdissey::InsertSkeletonInterfaces(int skeleton_id){
    
#ifdef PZDEBUG
    if(!fGeoMesh){
        DebugStop();
    }
#endif
    
    if(skeleton_id != 0){
        fSimulationData->SetSkeleton_material_Id(skeleton_id);
    }
    int level = fSimulationData->MHMResolution().second.first;
    long nel = fGeoMesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if (!gel || gel->Level() != level || gel->Dimension() != fGeoMesh->Dimension()) {
            continue;
        }
        
        if (gel->MaterialId() == 12 || gel->MaterialId() == 14) { // sideburden material
            continue;
        }
        
        int nsides = gel->NSides();
        for (int is = gel->NCornerNodes(); is<nsides; is++) {
            if (gel->SideDimension(is) != fGeoMesh->Dimension() - 1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == fGeoMesh->Dimension() - 1) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                TPZGeoElBC(gelside, fSimulationData->Skeleton_material_Id());
            }
        }
    }
    fGeoMesh->BuildConnectivity();
    this->PrintGeometry();
}

/** @brief Construc computational macro elements */
void TRMSpaceOdissey::BuildMacroElements()
{
    
    std::cout << "ndof parabolic before MHM substructuring = " << fMixedFluxPressureCmeshMHM->Solution().Rows() << std::endl;
    
#ifdef PZDEBUG
    if(!fMixedFluxPressureCmeshMHM){
        DebugStop();
    }
#endif
    
    long n_total_dof_avg = fMixedFluxPressureCmeshMHM->Solution().Rows();
    
    bool KeepOneLagrangian = true;
    typedef std::set<long> TCompIndexes;
    std::map<long, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = fMixedFluxPressureCmeshMHM->Reference();
    gmesh->ResetReference();
    fMixedFluxPressureCmeshMHM->LoadReferences();
    long nelg = gmesh->NElements();
    for (long el=0; el<nelg; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Father() != NULL) {
            continue;
        }
        if (gel->Dimension() == gmesh->Dimension() - 1 && gel->MaterialId() == fSimulationData->Skeleton_material_Id()) {
            continue;
        }
        
        bool IsSideBurden_Rock = gel->MaterialId() == 12 || gel->MaterialId() == 14;
        if (IsSideBurden_Rock) {
            continue;
        }
        
        long mapindex = gel->Index();
        if (gel->Dimension() == gmesh->Dimension() - 1) {
            TPZGeoElSide neighbour = gel->Neighbour(gel->NSides()-1);
            if (neighbour.Element()->Dimension() != gmesh->Dimension()) {
                DebugStop();
            }
            mapindex= neighbour.Element()->Index();
        }
        TPZStack<TPZCompElSide> highlevel;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        gelside.HigherLevelCompElementList3(highlevel, 0, 0);
        long nelst = highlevel.size();
        for (long elst=0; elst<nelst; elst++) {
            ElementGroups[mapindex].insert(highlevel[elst].Element()->Index());
        }
        if (gel->Reference()) {
            if (nelst) {
                DebugStop();
            }
            ElementGroups[mapindex].insert(gel->Reference()->Index());
        }
    }

    // omar:: Volumetric methodology for the creation macroblocks/macroelements
    
#ifdef PZDEBUG
    std::map<long,TCompIndexes>::iterator it;
    for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
        std::cout << "Group " << it->first << " with size " << it->second.size() << std::endl;
        std::cout << " elements ";
        std::set<long>::iterator its;
        for (its = it->second.begin(); its != it->second.end(); its++) {
            std::cout << *its << " ";
        }
        std::cout << std::endl;
    }
#endif
    
    std::set<long> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(fMixedFluxPressureCmeshMHM, ElementGroups, submeshindices, KeepOneLagrangian);
 

    fMixedFluxPressureCmeshMHM->ComputeNodElCon();
    fMixedFluxPressureCmeshMHM->CleanUpUnconnectedNodes();
    long n_micro = 0;
    long n_dof = 0;
    int count = 0;
    for (std::set<long>::iterator it=submeshindices.begin(); it != submeshindices.end(); it++) {
        TPZCompEl *cel = fMixedFluxPressureCmeshMHM->Element(*it);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, KeepOneLagrangian);
        subcmesh->SetAnalysisSkyline(0, 0, 0);
        n_micro += subcmesh->NElements();
        n_dof += cel->NEquations();
        count++;
    }
    
    fMixedFluxPressureCmeshMHM->ComputeNodElCon();
    fMixedFluxPressureCmeshMHM->CleanUpUnconnectedNodes();

    {
        std::string file = "mhm_report";
        
        if (fSimulationData->ReducedBasisResolution().first && !fSimulationData->ReducedBasisResolution().second.first) {
            file += "_RB_" + std::to_string(fSimulationData->m_RB_functions());
        }
        
        if (fSimulationData->IsAdataptedQ()) {
            file += "_A";
        }
        
        if (fSimulationData->IsEnhancedPressureQ()) {
            file += "_E";
        }
        
        if (fSimulationData->MHMResolution().first) {
            file +=  "_MHM_Hdiv_l_" + std::to_string(fSimulationData->MHMResolution().second.first);
        }
        
        if (fSimulationData->TransporResolution().first && !fSimulationData->IsOnePhaseQ()) {
            file += "_T_res_" + std::to_string(fSimulationData->TransporResolution().second);
        }
        
        file += ".txt";
        std::ofstream file_mhm(file.c_str());
        
        
        int n_micro_avg = int(REAL(n_micro)/REAL(ElementGroups.size()));
        int n_micro_dof_avg = int(REAL(n_total_dof_avg - fMixedFluxPressureCmeshMHM->Solution().Rows())/REAL(ElementGroups.size()*n_micro_avg));
        int n_macro_dof_avg = int(REAL(n_total_dof_avg - fMixedFluxPressureCmeshMHM->Solution().Rows())/REAL(ElementGroups.size()));
        file_mhm << "MHM-Hdiv:: Total number of dof = " << n_total_dof_avg << std::endl;
        file_mhm << "MHM-Hdiv:: Number of effective external dof = " << fMixedFluxPressureCmeshMHM->Solution().Rows() << std::endl;
        file_mhm << "MHM-Hdiv:: Number of condensed internal dof = " << n_total_dof_avg - fMixedFluxPressureCmeshMHM->Solution().Rows()  << std::endl;
        file_mhm << "MHM-Hdiv:: Number of macro elements " << ElementGroups.size()  << std::endl;
        file_mhm << "MHM-Hdiv:: Avg number of micro elements per macro element = " << n_micro_avg << std::endl;
        file_mhm << "MHM-Hdiv:: Avg number of condensed internal dof per macro element = " << n_macro_dof_avg << std::endl;
        file_mhm << "MHM-Hdiv:: Avg number of condensed internal dof per micro element = " << n_micro_dof_avg << std::endl;
        file_mhm.flush();
    }
}

/** @brief Destruct computational macro elements */
void TRMSpaceOdissey::UnwrapMacroElements()
{
    
    std::cout << " ndof parabolic with MHM substructuring = " << fMixedFluxPressureCmeshMHM->Solution().Rows() << std::endl;
    
#ifdef PZDEBUG
    if(!fMixedFluxPressureCmeshMHM){
        DebugStop();
    }
#endif
    
    long nelem = fMixedFluxPressureCmeshMHM->NElements();
    //deleting subcompmesh
    for(long i=0; i<nelem; i++){
        TPZCompEl *el = fMixedFluxPressureCmeshMHM->ElementVec()[i];
        TPZSubCompMesh * subc = dynamic_cast<TPZSubCompMesh*>(el);
        if(subc){
            long sub_nel = subc->NElements();
            long elindex;
            for (long isubcel = 0; isubcel < sub_nel; isubcel++) {
                
                TPZCompEl * cel =  subc->Element(isubcel);
                elindex = cel->Index();
                subc->TPZCompMesh::TransferElementTo(fMixedFluxPressureCmeshMHM, elindex);
            }

            delete subc;
        }
    }
    
    TPZCompMeshTools::UnCondensedElements(fMixedFluxPressureCmeshMHM);
    TPZCompMeshTools::UnGroupElements(fMixedFluxPressureCmeshMHM);
    
    std::cout << "ndof parabolic without MHM substructuring = " << fMixedFluxPressureCmeshMHM->Solution().Rows() << std::endl;
    
}

/** @brief Create a H1 computational mesh for Maurice Biot Linear Poroelasticity */
void TRMSpaceOdissey::CreateGeoMechanicMesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(dim,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fGeoMechanicsCmesh = new TPZCompMesh(fGeoMesh);
    
    // bc types
    int u_fixed   = 0;
    int u_h_fixed = 1;
    int s_v_free  = 3;
    int p_normal  = 6;
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        
        TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(rock_id,dim);
        if (fSimulationData->ReducedBasisResolution().first) {
            mat->Set_NStateVariables(1);
        }
        mat->SetSimulationData(this->SimulationData());
        fGeoMechanicsCmesh->InsertMaterialObject(mat);
        
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 0;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;
        
        for (int j = initial_bc; j < n_boundauries; j++) {
            
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }

            bc_item = bc[0];
            TPZMatWithMem<TRMMemory,TPZBndCond> * boundary_c = new TPZMatWithMem<TRMMemory,TPZBndCond>;
            boundary_c->SetNumLoadCases(1);
            boundary_c->SetMaterial(mat);
            boundary_c->SetId(bc_id);
            boundary_c->SetType(p_normal);
            boundary_c->SetValues(val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(0,boundary_data);
            fGeoMechanicsCmesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    // Sideburden
    int side_burden_rock = 14;
    int bc_W = 20;
    int bc_N = 19;
    int bc_E = 18;
    int bc_S = 17;
    int bc_T = 16;
    int bc_B = 15;
    
    if (dim == 2) {
        side_burden_rock = 12;
        bc_W = 14;
        bc_N = 15;
        bc_E = 16;
        bc_S = 13;
        bc_T = 1000;
        bc_B = 1000;
        
        TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(side_burden_rock,dim);
        if (fSimulationData->ReducedBasisResolution().first) {
            mat->Set_NStateVariables(1);
        }
        mat->SetSimulationData(this->SimulationData());
        fGeoMechanicsCmesh->InsertMaterialObject(mat);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * W_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        W_bndc->SetNumLoadCases(1);
        W_bndc->SetMaterial(mat);
        W_bndc->SetId(bc_W);
        W_bndc->SetType(u_h_fixed);
        W_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(W_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * N_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        N_bndc->SetNumLoadCases(1);
        N_bndc->SetMaterial(mat);
        N_bndc->SetId(bc_N);
        N_bndc->SetType(s_v_free);
        N_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(N_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * E_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        E_bndc->SetNumLoadCases(1);
        E_bndc->SetMaterial(mat);
        E_bndc->SetId(bc_E);
        E_bndc->SetType(u_h_fixed);
        E_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(E_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * S_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        S_bndc->SetNumLoadCases(1);
        S_bndc->SetMaterial(mat);
        S_bndc->SetId(bc_S);
        S_bndc->SetType(u_fixed);
        S_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(S_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * T_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        T_bndc->SetNumLoadCases(1);
        T_bndc->SetMaterial(mat);
        T_bndc->SetId(bc_T);
        T_bndc->SetType(u_fixed);
        T_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(T_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * B_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        B_bndc->SetNumLoadCases(1);
        B_bndc->SetMaterial(mat);
        B_bndc->SetId(bc_B);
        B_bndc->SetType(u_fixed);
        B_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(B_bndc);
        
    }
    else{
        
        TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(side_burden_rock,dim);
        mat->SetSimulationData(this->SimulationData());
        fGeoMechanicsCmesh->InsertMaterialObject(mat);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * W_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        W_bndc->SetNumLoadCases(1);
        W_bndc->SetMaterial(mat);
        W_bndc->SetId(bc_W);
        W_bndc->SetType(u_h_fixed);
        W_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(W_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * N_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        N_bndc->SetNumLoadCases(1);
        N_bndc->SetMaterial(mat);
        N_bndc->SetId(bc_N);
        N_bndc->SetType(u_h_fixed);
        N_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(N_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * E_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        E_bndc->SetNumLoadCases(1);
        E_bndc->SetMaterial(mat);
        E_bndc->SetId(bc_E);
        E_bndc->SetType(u_h_fixed);
        E_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(E_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * S_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        S_bndc->SetNumLoadCases(1);
        S_bndc->SetMaterial(mat);
        S_bndc->SetId(bc_S);
        S_bndc->SetType(u_h_fixed);
        S_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(S_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * T_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        T_bndc->SetNumLoadCases(1);
        T_bndc->SetMaterial(mat);
        T_bndc->SetId(bc_T);
        T_bndc->SetType(s_v_free);
        T_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(T_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * B_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        B_bndc->SetNumLoadCases(1);
        B_bndc->SetMaterial(mat);
        B_bndc->SetId(bc_B);
        B_bndc->SetType(u_fixed);
        B_bndc->SetValues(val1, val2);
        fGeoMechanicsCmesh->InsertMaterialObject(B_bndc);
        
    }
    
    
    fGeoMechanicsCmesh->SetDimModel(dim);
    fGeoMechanicsCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fGeoMechanicsCmesh->ApproxSpace().CreateWithMemory(true);
    fGeoMechanicsCmesh->AutoBuild();
    
    fGeoMechanicsCmesh->AdjustBoundaryElements();
    fGeoMechanicsCmesh->CleanUpUnconnectedNodes();
    
    TPZManVector<TPZCompMesh * ,1> meshvector(1);
    meshvector[0] = fBiotCmesh;
    
    // Transfer information
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fGeoMechanicsCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fGeoMechanicsCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fGeoMechanicsCmesh);
    
    long nel = fGeoMechanicsCmesh->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fGeoMechanicsCmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CmeshGeomechanic.txt");
    fGeoMechanicsCmesh->Print(out);
#endif
    
}

/** @brief Create Geomodes computational mesh */
void TRMSpaceOdissey::CreateGeoModesCmesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(dim,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fGalerkinProjectionsCmesh = new TPZCompMesh(fGeoMesh);
    
    // bc types
    int u_fixed   = 0;
    int u_h_fixed = 1;
    int s_v_free  = 3;
    int p_normal  = 6;
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        
        TRMPoroelasticModes * mat = new TRMPoroelasticModes(rock_id,dim,1);
        mat->SetSimulationData(this->SimulationData());
        fGalerkinProjectionsCmesh->InsertMaterialObject(mat);
        
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 0;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;
        
        for (int j = initial_bc; j < n_boundauries; j++) {
            
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[0];
            TPZMatWithMem<TRMMemory,TPZBndCond> * boundary_c = new TPZMatWithMem<TRMMemory,TPZBndCond>;
            boundary_c->SetNumLoadCases(1);
            boundary_c->SetMaterial(mat);
            boundary_c->SetId(bc_id);
            boundary_c->SetType(p_normal);
            boundary_c->SetValues(val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(0,boundary_data);
            fGalerkinProjectionsCmesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    // Sideburden
    int side_burden_rock = 14;
    int bc_W = 20;
    int bc_N = 19;
    int bc_E = 18;
    int bc_S = 17;
    int bc_T = 16;
    int bc_B = 15;
    
    if (dim == 2) {
        side_burden_rock = 12;
        bc_W = 14;
        bc_N = 15;
        bc_E = 16;
        bc_S = 13;
        bc_T = 1000;
        bc_B = 1000;
        
        TRMPoroelasticModes * mat = new TRMPoroelasticModes(side_burden_rock,dim,1);
        mat->SetSimulationData(this->SimulationData());
        fGalerkinProjectionsCmesh->InsertMaterialObject(mat);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * W_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        W_bndc->SetNumLoadCases(1);
        W_bndc->SetMaterial(mat);
        W_bndc->SetId(bc_W);
        W_bndc->SetType(u_h_fixed);
        W_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(W_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * N_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        N_bndc->SetNumLoadCases(1);
        N_bndc->SetMaterial(mat);
        N_bndc->SetId(bc_N);
        N_bndc->SetType(s_v_free);
        N_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(N_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * E_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        E_bndc->SetNumLoadCases(1);
        E_bndc->SetMaterial(mat);
        E_bndc->SetId(bc_E);
        E_bndc->SetType(u_h_fixed);
        E_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(E_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * S_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        S_bndc->SetNumLoadCases(1);
        S_bndc->SetMaterial(mat);
        S_bndc->SetId(bc_S);
        S_bndc->SetType(u_fixed);
        S_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(S_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * T_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        T_bndc->SetNumLoadCases(1);
        T_bndc->SetMaterial(mat);
        T_bndc->SetId(bc_T);
        T_bndc->SetType(u_fixed);
        T_bndc->SetValues(val1, val2);
//        fGalerkinProjectionsCmesh->InsertMaterialObject(T_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * B_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        B_bndc->SetNumLoadCases(1);
        B_bndc->SetMaterial(mat);
        B_bndc->SetId(bc_B);
        B_bndc->SetType(u_fixed);
        B_bndc->SetValues(val1, val2);
//        fGalerkinProjectionsCmesh->InsertMaterialObject(B_bndc);
        
    }
    else{
        
        TRMPoroelasticModes * mat = new TRMPoroelasticModes(side_burden_rock,dim,1);
        mat->SetSimulationData(this->SimulationData());
        fGalerkinProjectionsCmesh->InsertMaterialObject(mat);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * W_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        W_bndc->SetNumLoadCases(1);
        W_bndc->SetMaterial(mat);
        W_bndc->SetId(bc_W);
        W_bndc->SetType(u_h_fixed);
        W_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(W_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * N_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        N_bndc->SetNumLoadCases(1);
        N_bndc->SetMaterial(mat);
        N_bndc->SetId(bc_N);
        N_bndc->SetType(u_h_fixed);
        N_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(N_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * E_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        E_bndc->SetNumLoadCases(1);
        E_bndc->SetMaterial(mat);
        E_bndc->SetId(bc_E);
        E_bndc->SetType(u_h_fixed);
        E_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(E_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * S_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        S_bndc->SetNumLoadCases(1);
        S_bndc->SetMaterial(mat);
        S_bndc->SetId(bc_S);
        S_bndc->SetType(u_h_fixed);
        S_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(S_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * T_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        T_bndc->SetNumLoadCases(1);
        T_bndc->SetMaterial(mat);
        T_bndc->SetId(bc_T);
        T_bndc->SetType(s_v_free);
        T_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(T_bndc);
        
        TPZMatWithMem<TRMMemory,TPZBndCond> * B_bndc = new TPZMatWithMem<TRMMemory,TPZBndCond>;
        B_bndc->SetNumLoadCases(1);
        B_bndc->SetMaterial(mat);
        B_bndc->SetId(bc_B);
        B_bndc->SetType(u_fixed);
        B_bndc->SetValues(val1, val2);
        fGalerkinProjectionsCmesh->InsertMaterialObject(B_bndc);
        
    }
    
    
    fGalerkinProjectionsCmesh->SetDimModel(dim);
    fGalerkinProjectionsCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fGalerkinProjectionsCmesh->ApproxSpace().CreateWithMemory(true);
    fGalerkinProjectionsCmesh->AutoBuild();
    
    fGalerkinProjectionsCmesh->AdjustBoundaryElements();
    fGalerkinProjectionsCmesh->CleanUpUnconnectedNodes();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    meshvector[0] = fGP_BiotCmesh;
    meshvector[1] = fGeoPressureCmesh;
    
    // Transfer information
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fGalerkinProjectionsCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fGalerkinProjectionsCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fGalerkinProjectionsCmesh);
    
    long nel = fGalerkinProjectionsCmesh->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fGalerkinProjectionsCmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CmeshGalerkinProjections.txt");
    fGalerkinProjectionsCmesh->Print(out);
#endif
    
}


/** @brief Create a Mixed computational mesh Hdiv-L2 */
void TRMSpaceOdissey::CreateMixedCmesh(){
    
#ifdef PZDEBUG
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
#endif
    
    int dim = fGeoMesh->Dimension();
    int flux_or_pressure = 0;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fMixedFluxPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id,dim);
        mat->SetSimulationData(this->SimulationData());        
        fMixedFluxPressureCmesh->InsertMaterialObject(mat);
        
        // Inserting boundary materials
        if (rock_id == 5) { // Reservoir
            n_boundauries = 6;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;

        for (int j = initial_bc; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
//            bc_item = bc[flux_or_pressure];
            
//            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
//            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
//            boundary_c->SetTimedependentBCForcingFunction(boundary_data);
//            fMixedFluxPressureCmesh->InsertMaterialObject(boundary_c);
            
            bc_item = bc[flux_or_pressure];
            TPZMatWithMem<TRMMemory,TPZBndCond> * boundary_c = new TPZMatWithMem<TRMMemory,TPZBndCond>;
            boundary_c->SetNumLoadCases(1);
            boundary_c->SetMaterial(mat);
            boundary_c->SetId(bc_id);
            boundary_c->SetType(bc_item.first);
            boundary_c->SetValues(val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(0,boundary_data);
            fMixedFluxPressureCmesh->InsertMaterialObject(boundary_c);

        }
        
    }

    fMixedFluxPressureCmesh->SetDimModel(dim);
    fMixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fMixedFluxPressureCmesh->ApproxSpace().CreateWithMemory(true);
    fMixedFluxPressureCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    meshvector[0] = fFluxCmesh;
    meshvector[1] = fPressureCmesh;
    
    // Trensfer information
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fMixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fMixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fMixedFluxPressureCmesh);
    
    long nel = fMixedFluxPressureCmesh->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fMixedFluxPressureCmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
    fMixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshMixed.txt");
    fMixedFluxPressureCmesh->Print(out);
#endif
    
}

/** @brief Create a Mixed MHM computational mesh Hdiv-L2 */
void TRMSpaceOdissey::CreateMixedCmeshMHM(){
    
#ifdef PZDEBUG
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
#endif
    
    int dim = fGeoMesh->Dimension();
    int flux_or_pressure = 0;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fMixedFluxPressureCmeshMHM = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id,dim);
        mat->SetSimulationData(this->SimulationData());
        fMixedFluxPressureCmeshMHM->InsertMaterialObject(mat);
        
        // Inserting boundary materials
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 6;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
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
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fMixedFluxPressureCmeshMHM->InsertMaterialObject(boundary_c);
            
        }
        
    }
    
    fMixedFluxPressureCmeshMHM->SetDimModel(dim);
    fMixedFluxPressureCmeshMHM->SetAllCreateFunctionsMultiphysicElemWithMem();
    fMixedFluxPressureCmeshMHM->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    meshvector[0] = fFluxCmesh;
    meshvector[1] = fPressureCmesh;
    
    // Trensfer information
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fMixedFluxPressureCmeshMHM);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fMixedFluxPressureCmeshMHM);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fMixedFluxPressureCmeshMHM);
    
    long nel = fMixedFluxPressureCmeshMHM->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fMixedFluxPressureCmeshMHM->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
    fMixedFluxPressureCmeshMHM->CleanUpUnconnectedNodes();
    
}

/** @brief Create a Mixed-Transport muliphase computational mesh Hdiv-L2-L2-L2 */
void TRMSpaceOdissey::CreateMultiphaseCmesh(){
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    int flux_or_pressure = 0;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fMonolithicMultiphaseCmesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMultiphase * mat = new TRMMultiphase(rock_id,dim);
        mat->SetSimulationData(fSimulationData);
        fMonolithicMultiphaseCmesh->InsertMaterialObject(mat);
        
        // Inserting boundary materials
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 6;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;

        for (int j = initial_bc; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[flux_or_pressure];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fMonolithicMultiphaseCmesh->InsertMaterialObject(boundary_c);
            
        }
        
    }

    fMonolithicMultiphaseCmesh->SetDimModel(dim);
    fMonolithicMultiphaseCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fMonolithicMultiphaseCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(3);
    
    if (fSimulationData->IsOnePhaseQ()) {
        meshvector[0] = fBiotCmesh;
        meshvector[1] = fFluxCmesh;
        meshvector[2] = fPressureCmesh;
    }
    
    if (fSimulationData->IsTwoPhaseQ()) {
        meshvector.Resize(4);
        meshvector[0] = fBiotCmesh;
        meshvector[1] = fFluxCmesh;
        meshvector[2] = fPressureCmesh;
        meshvector[3] = fAlphaSaturationMesh;
    }
    
    if (fSimulationData->IsThreePhaseQ()) {
        meshvector.Resize(5);
        meshvector[0] = fBiotCmesh;
        meshvector[2] = fFluxCmesh;
        meshvector[3] = fPressureCmesh;
        meshvector[4] = fAlphaSaturationMesh;
        meshvector[5] = fBetaSaturationMesh;
    }

    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fMonolithicMultiphaseCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fMonolithicMultiphaseCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fMonolithicMultiphaseCmesh);
    
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
void TRMSpaceOdissey::CreateInterfacesInside(TPZCompMesh * cmesh){ //@Omar:: It is required to robust for several materials in the same hydraulic unit!!!
    
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
    if (!fMixedFluxPressureCmesh) {
        std::cout<< "No multiphysic computational mesh " << std::endl;
        DebugStop();
    }
    
    
    fMixedFluxPressureCmesh->Reference()->ResetReference();
    fMixedFluxPressureCmesh->LoadReferences();
    
    fMixedFluxPressureCmesh->ComputeNodElCon();
    // create condensed elements
    // increase the NumElConnected of one pressure connects in order to prevent condensation
    for (long icel=0; icel < fMixedFluxPressureCmesh->NElements(); icel++) {
        TPZCompEl  * cel = fMixedFluxPressureCmesh->Element(icel);
        
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
    meshvector[0] = fFluxCmesh;
    meshvector[1] = fPressureCmesh;
    
    int DOF = meshvector[0]->NEquations() + meshvector[1]->NEquations();
    REAL PercentCondensedDOF = 100.0*(1.0 - REAL(fMixedFluxPressureCmesh->NEquations())/REAL(DOF));
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
    
    int dim = fGeoMesh->Dimension();
    int saturation = 0;
    int sorder = fSOrder;
    
    if(sorder > 0 ){
        fSimulationData->SetUseGradientR(true);
    }
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
    // Malha computacional
    fAlphaSaturationMesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TPZMatPoisson3d * mat = new TPZMatPoisson3d(rock_id,dim);
        fAlphaSaturationMesh->InsertMaterialObject(mat);
        
        if(fSimulationData->UseGradientR()){
            int L2Proj_material = fSimulationData->L2_Projection_material_Id();
            TPZVec<STATE> sol(1,0.);
            TPZL2Projection *matl2proj = new TPZL2Projection(L2Proj_material,dim,mat->NStateVariables(),sol);
            fAlphaSaturationMesh->InsertMaterialObject(matl2proj);
        }
        
        // Inserting volumetric materials
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 6;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;
        std::pair< int, TPZFunction<REAL> * > bc_item;
        TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
        for (int j = initial_bc; j < n_boundauries; j++) {
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
    PrintGeometry();
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
    
    int dim = fGeoMesh->Dimension();
    int saturation = 0;
    int sorder = fSOrder;
    
    if(sorder > 0 ){
        fSimulationData->SetUseGradientR(true);
    }
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
    // Malha computacional
    fBetaSaturationMesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TPZMatPoisson3d * mat = new TPZMatPoisson3d(rock_id,dim);
        fBetaSaturationMesh->InsertMaterialObject(mat);
        
        if(fSimulationData->UseGradientR()){
            int L2Proj_material = fSimulationData->L2_Projection_material_Id();
            TPZVec<STATE> sol(1,0.);
            TPZL2Projection *matl2proj = new TPZL2Projection(L2Proj_material,dim,mat->NStateVariables(),sol);
            fBetaSaturationMesh->InsertMaterialObject(matl2proj);
        }
        
        // Inserting volumetric materials
        
        if (rock_id == 5) { // Reservoir
            n_boundauries = 6;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;
        std::pair< int, TPZFunction<REAL> * > bc_item;
        TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
        for (int j = initial_bc; j < n_boundauries; j++) {
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
    PrintGeometry();
#endif
    
}

/** @brief Create a multiphysics computational mesh L2 */
void TRMSpaceOdissey::CreateTransportMesh(){
    // Second option put all the transpor meshes inside a multiphysics mesh

    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = fGeoMesh->Dimension();
    int saturation = 0;
    int sorder = fSOrder;
    int interface_id = fSimulationData->InterfacesMatId();
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    std::pair< int, TPZFunction<REAL> * > bc_item;
    TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
    
    // Malha computacional
    fTransportMesh = new TPZCompMesh(fGeoMesh);
    
    // Inserting volumetric materials
    int n_rocks = this->SimulationData()->RawData()->fOmegaIds.size();
    int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
    
    int initial_bc = 0;
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMPhaseTransport * mat = new TRMPhaseTransport(rock_id,dim);
        mat->SetSimulationData(fSimulationData);
        fTransportMesh->InsertMaterialObject(mat);
        
        TRMPhaseInterfaceTransport * matint = new TRMPhaseInterfaceTransport(interface_id);
        matint->SetSimulationData(fSimulationData);
        fTransportMesh->InsertMaterialObject(matint);
        
        // Inserting volumetric materials
        if (rock_id == 5) { // Reservoir
            n_boundauries = 6;
            initial_bc = 0;
        }
        
        if (rock_id == 6) { // Wellbore productors
            n_boundauries = 8;
            initial_bc = 6;
        }
        
        if (rock_id == 7) { // Wellbore injectors
            n_boundauries = 10;
            initial_bc = 8;
        }
        
        int bc_id = 0;

        for (int j = initial_bc; j < n_boundauries; j++) {
            bc_id   = this->SimulationData()->RawData()->fGammaIds[j];
            
            if (fSimulationData->IsInitialStateQ()) {
                bc      = this->SimulationData()->RawData()->fIntial_bc_data[j];
            }
            else{
                bc      = this->SimulationData()->RawData()->fRecurrent_bc_data[j];
            }
            
            bc_item = bc[saturation];
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> * boundary_bc = new TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>;
            boundary_bc->SetNumLoadCases(1);
            boundary_bc->SetMaterial(matint);
            boundary_bc->SetId(bc_id);
            boundary_bc->SetType(bc_item.first);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_bc->SetTimedependentBCForcingFunction(0,boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fTransportMesh->InsertMaterialObject(boundary_bc);
        }
        
    }
    
    { // It Works but doesn't have sense
        int res_id = 5;
        int pro_id = 6;
        int inj_id = 7;
        
        fGeoMesh->AddInterfaceMaterial(res_id, pro_id,interface_id);
        fGeoMesh->AddInterfaceMaterial(pro_id, res_id,interface_id);
        fGeoMesh->AddInterfaceMaterial(res_id, inj_id,interface_id);
        fGeoMesh->AddInterfaceMaterial(inj_id, res_id,interface_id);
        
        fGeoMesh->AddInterfaceMaterial(res_id, res_id,interface_id);
        fGeoMesh->AddInterfaceMaterial(pro_id, pro_id,interface_id);
        fGeoMesh->AddInterfaceMaterial(inj_id, inj_id,interface_id);
    }
    
    fTransportMesh->SetDimModel(dim);
    fTransportMesh->SetDefaultOrder(sorder);
    fTransportMesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fTransportMesh->ApproxSpace().CreateWithMemory(true);// Force the creation of interfaces with memory.
    fTransportMesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector;
    
    if(this->SimulationData()->IsTwoPhaseQ()){
        
        meshvector.Resize(1);
        meshvector[0] = fAlphaSaturationMesh;
        
        // Transferindo para a multifisica
        TPZBuildMultiphysicsMesh::AddElements(meshvector, fTransportMesh);
        TPZBuildMultiphysicsMesh::AddConnects(meshvector, fTransportMesh);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fTransportMesh);
        
    }
    
    if(this->SimulationData()->IsThreePhaseQ()){

        meshvector.Resize(2);
        meshvector[0] = fAlphaSaturationMesh;
        meshvector[1] = fBetaSaturationMesh;
        
        // Transferindo para a multifisica
        TPZBuildMultiphysicsMesh::AddElements(meshvector, fTransportMesh);
        TPZBuildMultiphysicsMesh::AddConnects(meshvector, fTransportMesh);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fTransportMesh);
    }

    
    fGeoMesh->ResetReference();
    fTransportMesh->LoadReferences();
    long nel = fTransportMesh->ElementVec().NElements();
    // Creation of interface elements
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * cel = fTransportMesh->ElementVec()[el];
        if(!cel) continue;
        TPZGeoEl * gel = cel->Reference();
        if(!gel) {continue;}
        if(gel->HasSubElement()) {continue;}
        int index = cel ->Index();
        if(cel->Dimension() == fTransportMesh->Dimension())
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
        
        if(cel->Dimension() != dim){
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
    PrintGeometry();
#endif
    
    
}

void TRMSpaceOdissey::UniformRefinement_cmesh(TPZCompMesh  *cmesh, int n_ref)
{
    TPZVec<long > subindex;
    if (n_ref == 0) {
        return;
    }
    for (long iref = 0; iref < n_ref; iref++) {
        TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
        long nel = elvec.NElements();
        for(long el=0; el < nel; el++){
            TPZCompEl * cel = elvec[el];
            if(!cel) continue;
            long ind = cel->Index();
            cel->Divide(ind, subindex);
        }
    }
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
    
    long last_node = fGeoMesh->NNodes() - 1;
    long last_element = fGeoMesh->NElements() - 1;
    long node_id = fGeoMesh->NodeVec()[last_node].Id();
    long element_id = fGeoMesh->Element(last_element)->Id();
    const std::string name("GID Reservoir geometry");
    fGeoMesh->SetName(name);
    fGeoMesh->SetMaxNodeId(node_id);
    fGeoMesh->SetMaxElementId(element_id);
    
}

/** @brief Create a reservoir-box geometry with cylindrical wells */
void TRMSpaceOdissey::CreateGeometricGmshMesh(std::string &grid){
    
    TRMGmshReader Geometry;
    fGeoMesh = Geometry.GeometricGmshMesh(grid);
    const std::string name("Reservoir with cylindrical wells");
    fGeoMesh->SetName(name);
    
    if (fGeoMesh->Dimension() == 3) {
        if (Geometry.NHexahedra() != 0 &&  Geometry.NTetrahera() == 0 && Geometry.NPrisms() == 0) {
            fIsHexaDominatedQ = true;
        }
        if (Geometry.NHexahedra() == 0 &&  Geometry.NTetrahera() != 0 && Geometry.NPrisms() == 0) {
            fIsTetraDominatedQ = true;
        }
        if (Geometry.NHexahedra() == 0 &&  Geometry.NTetrahera() == 0 && Geometry.NPrisms() != 0) {
            fIsPrismDominatedQ = true;
        }
    }
}

/** @brief Create a reservoir-box geometry */
void TRMSpaceOdissey::CreateGeometricExtrudedGIDMesh(std::string &grid, TPZManVector<REAL,2> dz){

    
    TPZReadGIDGrid GeometryInfo;
    REAL s = 1.0;
    GeometryInfo.SetfDimensionlessL(s);
    TPZGeoMesh * GeoMesh2D = GeometryInfo.GeometricGIDMesh(grid);
    GeoMesh2D->SetDimension(2);
    
    
    REAL t=0.0;
    REAL dt;
    int n;
    bool IsTetrahedronMeshQ = false;
    
    int bc_B =  this->SimulationData()->RawData()->fGammaIds[0];
    int bc_T =  this->SimulationData()->RawData()->fGammaIds[0];
    
//    int bc_B =  this->SimulationData()->RawData()->fGammaIds[3];
//    int bc_T =  this->SimulationData()->RawData()->fGammaIds[4];
    
    TPZHierarquicalGrid CreateGridFrom2D(GeoMesh2D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncZ = new TPZDummyFunction<STATE>(ParametricfunctionZ);
    CreateGridFrom2D.SetParametricFunction(ParFuncZ);
    CreateGridFrom2D.SetFrontBackMatId(bc_B,bc_T);
    if(IsTetrahedronMeshQ){
        CreateGridFrom2D.SetTriangleExtrusion();
        CreateGridFrom2D.SetTetrahedonExtrusion();
    }
    
    
    dt = dz[0];
    n = int(dz[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fGeoMesh = CreateGridFrom2D.ComputeExtrusion(t, dt, n);
    
    long last_node = fGeoMesh->NNodes() - 1;
    long last_element = fGeoMesh->NElements() - 1;
    long node_id = fGeoMesh->NodeVec()[last_node].Id();
    long element_id = fGeoMesh->Element(last_element)->Id();
    const std::string name("Reservoir with vertical extrusion");
    fGeoMesh->SetName(name);
    fGeoMesh->SetMaxNodeId(node_id);
    fGeoMesh->SetMaxElementId(element_id);
    
}

/** @brief Create a reservoir-box geometry */
void TRMSpaceOdissey::CreateGeometricBoxMesh2D(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy){
    
    REAL t=0.0;
    REAL dt;
    int n;
    bool IsTetrahedronMeshQ = false;
    
    int rock =  this->SimulationData()->RawData()->fOmegaIds[0];
    
    int bc_W =  this->SimulationData()->RawData()->fGammaIds[0];
    int bc_E =  this->SimulationData()->RawData()->fGammaIds[1];
    int bc_S =  this->SimulationData()->RawData()->fGammaIds[2];
    int bc_N =  this->SimulationData()->RawData()->fGammaIds[3];
    
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
    if(IsTetrahedronMeshQ){
        CreateGridFrom1D.SetTriangleExtrusion();
    }
    
    
    dt = dy[0];
    n = int(dy[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fGeoMesh = CreateGridFrom1D.ComputeExtrusion(t, dt, n);
    

    long last_node = fGeoMesh->NNodes() - 1;
    long last_element = fGeoMesh->NElements() - 1;
    long node_id = fGeoMesh->NodeVec()[last_node].Id();
    long element_id = fGeoMesh->Element(last_element)->Id();
    const std::string name("Reservoir box 2D");
    fGeoMesh->SetName(name);
    fGeoMesh->SetMaxNodeId(node_id);
    fGeoMesh->SetMaxElementId(element_id);
    
}

/** @brief Create a reservoir-box geometry */
void TRMSpaceOdissey::CreateGeometricBoxMesh(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy, TPZManVector<REAL,2> dz){
    
    REAL t=0.0;
    REAL dt;
    int n;
    bool IsTetrahedronMeshQ = false;

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
    if(IsTetrahedronMeshQ){
        CreateGridFrom1D.SetTriangleExtrusion();
    }
    
    
    dt = dy[0];
    n = int(dy[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh2D = CreateGridFrom1D.ComputeExtrusion(t, dt, n);
        
    TPZHierarquicalGrid CreateGridFrom2D(GeoMesh2D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncZ = new TPZDummyFunction<STATE>(ParametricfunctionZ);
    CreateGridFrom2D.SetParametricFunction(ParFuncZ);
    CreateGridFrom2D.SetFrontBackMatId(bc_B,bc_T);
    if(IsTetrahedronMeshQ){
        CreateGridFrom2D.SetTriangleExtrusion();
        CreateGridFrom2D.SetTetrahedonExtrusion();
    }
    
    
    dt = dz[0];
    n = int(dz[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fGeoMesh = CreateGridFrom2D.ComputeExtrusion(t, dt, n);

    long last_node = fGeoMesh->NNodes() - 1;
    long last_element = fGeoMesh->NElements() - 1;
    long node_id = fGeoMesh->NodeVec()[last_node].Id();
    long element_id = fGeoMesh->Element(last_element)->Id();
    const std::string name("Reservoir box");
    fGeoMesh->SetName(name);
    fGeoMesh->SetMaxNodeId(node_id);
    fGeoMesh->SetMaxElementId(element_id);
    
}


void TRMSpaceOdissey::ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0];
    X[1] = 0.0*sin(0.01*par[0]);
    X[2] = 0.0;
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

void TRMSpaceOdissey::CElemtentRefinement(TPZCompMesh  *cmesh, int element_index)
{
    
    TPZVec<long > subindex;
    long nel = cmesh->ElementVec().NElements();
    for(long el=0; el < nel; el++){
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        long ind = compEl->Index();
        if(ind==element_index){
            compEl->Divide(element_index, subindex, 1);
        }
    }
}

void TRMSpaceOdissey::CMeshRefinement(TPZCompMesh  *cmesh, int ndiv)
{
    
    TPZVec<long > subindex;
    for (long iref = 0; iref < ndiv; iref++) {
        TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
        long nel = elvec.NElements();
        for(long el=0; el < nel; el++){
            TPZCompEl * compEl = elvec[el];
            if(!compEl) continue;
            long ind = compEl->Index();
            compEl->Divide(ind, subindex, 0);
        }
    }
}

/** @brief Apply uniform refinement on the Geometric mesh */
void TRMSpaceOdissey::UniformRefinement(int n_ref){
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        long n = fGeoMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fGeoMesh->ElementVec() [i];
            if (gel->Dimension() != 0) gel->Divide (sons);
        }//for i
    }//ref
    fGeoMesh->BuildConnectivity();
}

void TRMSpaceOdissey::UniformRefineTetrahedrons(int n_ref){
    
    TPZAutoPointer<TPZRefPattern> refp3D;
    
    {// needed!
        char buf[] =
        "10 9 "
        "-50 Tet0000111111111 "
        "0 0 0 "
        "1 0 0 "
        "0 1 0 "
        "0 0 1 "
        "0.5 0 0 "
        "0 0.5 0 "
        "0 0 0.5 "
        "0.5 0.5 0 "
        "0 0.5 0.5 "
        "0.5 0 0.5 "
        "4 4 0  1  2  3 "
        "4 4 0  4  5  6 "
        "4 4 4  1  7  9 "
        "4 4 7  2  5  8 "
        "4 4 6  9  8  3 "
        "4 4 4  9  6  5 "
        "4 4 5  8  6  9 "
        "4 4 7  8  9  5 "
        "4 4 4  7  5  9 ";
        std::istringstream str(buf);
        refp3D = new TPZRefPattern(str);
        refp3D->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp3D);
        if(!refp3D)
        {
            DebugStop();
        }
    }
    
    //    TPZAutoPointer<TPZRefPattern> refp3D = gRefDBase.FindRefPattern("UnifTet");
    
    if(!refp3D)
    {
        DebugStop();
    }
    
    TPZGeoEl * gel = NULL;
    for(int r = 0; r < n_ref; r++)
    {
        int nels = fGeoMesh->NElements();
        for(int iel = 0; iel < nels; iel++)
        {
            gel = fGeoMesh->ElementVec()[iel];
            if(!gel) DebugStop();
            if(gel->Dimension()==3)
            {
                gel->SetRefPattern(refp3D);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
            if(gel->Dimension()==2)
            {
                //                gel->SetRefPattern(refp3D);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
            
        }
    }
    fGeoMesh->BuildConnectivity();
}

/** @brief Apply uniform refinement at specific material id */
void TRMSpaceOdissey::UniformRefinement_at_MaterialId(int n_ref, int mat_id){
    
    int dim = fGeoMesh->Dimension();
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        long n = fGeoMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fGeoMesh->ElementVec() [i];
            if (gel->Dimension() == dim || gel->Dimension() == dim - 1){
                if (gel->MaterialId()== mat_id){
                    gel->Divide(sons);
                }
            }
        }//for i
    }//ref
}

/** @brief Apply uniform refinement around at specific material id */
void TRMSpaceOdissey::UniformRefinement_Around_MaterialId(int n_ref, int mat_id){
    
    int dim = fGeoMesh->Dimension();
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        long n = fGeoMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fGeoMesh->ElementVec() [i];
            if (gel->Dimension() == dim || gel->Dimension() == dim - 1){
                
                if(gel->HasSubElement()){
                    continue;
                }
                
                if (gel->MaterialId()== mat_id){
                    
                    int gel_nsides = gel->NSides();
                    TPZGeoElSide gel_itself =  gel->Neighbour(gel_nsides-1);
                    
                    if(!gel_itself.Element()){
                        continue;
                    }
                    
                    if(gel_itself.Element()->HasSubElement()){
                        continue;
                    }
                    
                    int high_gel_nsides = gel_itself.Element()->NSides() - 1;
                    for (int is = gel_itself.Element()->NNodes() ; is < high_gel_nsides; is++) {
                        TPZGeoElSide side = gel_itself.Element()->Neighbour(is);
                        if(!side.Element()){
                            continue;
                        }
                        
                        if(side.Element()->Dimension() != dim-1 || side.Element()->HasSubElement()){
                            continue;
                        }
                        
                        side.Element()->Divide(sons);
                    }
                    
                    gel_itself.Element()->Divide(sons);
                    
                }
            }
        }//for i
    }//ref
    
    fGeoMesh->BuildConnectivity();
}

/** @brief Apply uniform refinement on the Geometric mesh */
void TRMSpaceOdissey::UniformRefinement_at_Father(int n_ref, int father_index){
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        long n = fGeoMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fGeoMesh->ElementVec() [i];
            if(gel->HasSubElement()){
                continue;
            }
            if(gel->Index() == father_index ){
                if (gel->Dimension() != 0) gel->Divide (sons);
            }

        }//for i
    }//ref
    
    fGeoMesh->BuildConnectivity();
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


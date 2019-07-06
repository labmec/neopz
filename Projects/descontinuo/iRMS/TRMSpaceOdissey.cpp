  //
//  TRMSpaceOdissey.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpaceOdissey.h"



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
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        
        TRMBiotPoroelasticity * mat = new TRMBiotPoroelasticity(rock_id,dim);
        fBiotCmesh->InsertMaterialObject(mat);
        
        
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
            
            bc_item = bc[Sigma_or_displacement];
            TPZMaterial * boundary_c = mat->CreateBC(mat, bc_id, bc_item.first, val1, val2);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fBiotCmesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
    fBiotCmesh->SetDimModel(dim);
    fBiotCmesh->SetDefaultOrder(uorder);
    fBiotCmesh->SetAllCreateFunctionsContinuous();
    fBiotCmesh->AutoBuild();
    
    fBiotCmesh->AdjustBoundaryElements();
    fBiotCmesh->CleanUpUnconnectedNodes();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshBiot.txt");
    fBiotCmesh->Print(out);
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
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id,dim);
        fFluxCmesh->InsertMaterialObject(mat);
        
        TRMMixedDarcy * mat_skeleton = new TRMMixedDarcy(fSimulationData->Skeleton_material_Id(),dim-1);
        fFluxCmesh->InsertMaterialObject(mat_skeleton); // @omar::  skeleton material inserted
        
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
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fFluxCmesh->InsertMaterialObject(boundary_c);
        }
        
    }
    
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
    
    int ncon = fPressureCmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = fPressureCmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
//    TPZDummyFunction<STATE> dummy(PressFunc);
//    TPZCompMeshTools::LoadSolution(fPressureCmesh, dummy);
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

/** @brief Build MHM form the current hdvi mesh */
void TRMSpaceOdissey::BuildMixed_Mesh(){
    
    this->CreateFluxCmesh();
    this->CreatePressureCmesh();
    this->CreateMixedCmesh();
    
}

/** @brief Build MHM form the current hdvi mesh */
void TRMSpaceOdissey::BuildMHM_Mesh(){
    
    this->CreateFluxCmesh();
    this->CreatePressureCmesh();
    
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
    
    this->UnwrapMacroElements();
    
#ifdef PZDEBUG
    std::ofstream out_unwrap("CmeshMixedMHMUnWrap.txt");
    this->MixedFluxPressureCmeshMHM()->Print(out_unwrap);
#endif
    
    
}

/** @brief Sparated connects by given selected skeleton ids */
void TRMSpaceOdissey::SeparateConnectsBySkeletonIds(TPZVec<int64_t> skeleton_ids){
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
    int64_t nel = fFluxCmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
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
                int64_t cindex = fFluxCmesh->AllocateNewConnect(c);
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
    
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if (!gel || gel->Level() != 0 || gel->Dimension() != fGeoMesh->Dimension()) {
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
    
    TPZMatLaplacian *material = new TPZMatLaplacian(_ReservMatId,3);
    material->SetForcingFunction(One,fPOrder);
    fH1Cmesh->InsertMaterialObject(material);

    TPZFNMatrix<1> val1(1,1,0.),val2(1,1,0);
    TPZBndCond *inflow = new TPZBndCond(material,_ConfinementReservBCbottom,0,val1,val2);
    val2(0,0) = 0.;
    TPZBndCond *outflow = new TPZBndCond(material,_ConfinementReservBCtop,0,val1,val2);
    
#ifdef PZDEBUG
    if(!fMixedFluxPressureCmeshMHM){
        DebugStop();
    }
#endif
    
    bool KeepOneLagrangian = true;
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = fMixedFluxPressureCmeshMHM->Reference();
    gmesh->ResetReference();
    fMixedFluxPressureCmeshMHM->LoadReferences();
    int64_t nelg = gmesh->NElements();
    for (int64_t el=0; el<nelg; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Father() != NULL) {
            continue;
        }
        if (gel->Dimension() == gmesh->Dimension() - 1 && gel->MaterialId() == fSimulationData->Skeleton_material_Id()) {
            continue;
        }
        int64_t mapindex = gel->Index();
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
        int64_t nelst = highlevel.size();
        for (int64_t elst=0; elst<nelst; elst++) {
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
    std::map<int64_t,TCompIndexes>::iterator it;
    for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
        std::cout << "Group " << it->first << " with size " << it->second.size() << std::endl;
        std::cout << " elements ";
        std::set<int64_t>::iterator its;
        for (its = it->second.begin(); its != it->second.end(); its++) {
            std::cout << *its << " ";
        }
        std::cout << std::endl;
    }
#endif
    
    std::set<int64_t> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(fMixedFluxPressureCmeshMHM, ElementGroups, submeshindices, KeepOneLagrangian);
 
    std::cout << "Inserting " << ElementGroups.size()  <<  " macro elements into MHM substructures" << std::endl;
    fMixedFluxPressureCmeshMHM->ComputeNodElCon();
    fMixedFluxPressureCmeshMHM->CleanUpUnconnectedNodes();
    for (std::set<int64_t>::iterator it=submeshindices.begin(); it != submeshindices.end(); it++) {
        TPZCompEl *cel = fMixedFluxPressureCmeshMHM->Element(*it);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, KeepOneLagrangian);
        subcmesh->SetAnalysisSkyline(0, 0, 0);
    }
    fMixedFluxPressureCmeshMHM->ComputeNodElCon();
    fMixedFluxPressureCmeshMHM->CleanUpUnconnectedNodes();

}

/** @brief Destruct computational macro elements */
void TRMSpaceOdissey::UnwrapMacroElements()
{
    
    std::cout << "ndof parabolic with MHM substructuring = " << fMixedFluxPressureCmeshMHM->Solution().Rows() << std::endl;
    
#ifdef PZDEBUG
    if(!fMixedFluxPressureCmeshMHM){
        DebugStop();
    }
#endif
    
    int64_t nelem = fMixedFluxPressureCmeshMHM->NElements();
    //deleting subcompmesh
    for(int64_t i=0; i<nelem; i++){
        TPZCompEl *el = fMixedFluxPressureCmeshMHM->ElementVec()[i];
        TPZSubCompMesh * subc = dynamic_cast<TPZSubCompMesh*>(el);
        if(subc){
            int64_t sub_nel = subc->NElements();
            int64_t elindex;
            for (int64_t isubcel = 0; isubcel < sub_nel; isubcel++) {
                
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
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id,dim);
        mat->SetSimulationData(this->SimulationData());        
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
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_c->SetTimedependentBCForcingFunction(boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
            fMixedFluxPressureCmesh->InsertMaterialObject(boundary_c);

        }
        
    }

    fMixedFluxPressureCmesh->SetDimModel(dim);
    fMixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    fMixedFluxPressureCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    meshvector[0] = fFluxCmesh;
    meshvector[1] = fPressureCmesh;
    
    // Trensfer information
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fMixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fMixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fMixedFluxPressureCmesh);
    
    int64_t nel = fMixedFluxPressureCmesh->NElements();
    TPZVec<int64_t> indices;
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fMixedFluxPressureCmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
    fMixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
//    TPZCompMeshTools::OptimizeBandwidth(fMixedFluxPressureCmesh);
    
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
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMixedDarcy * mat = new TRMMixedDarcy(rock_id,dim);
        mat->SetSimulationData(this->SimulationData());
        fMixedFluxPressureCmeshMHM->InsertMaterialObject(mat);
        
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
    
    int64_t nel = fMixedFluxPressureCmeshMHM->NElements();
    TPZVec<int64_t> indices;
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fMixedFluxPressureCmeshMHM->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
    fMixedFluxPressureCmeshMHM->CleanUpUnconnectedNodes();
    
    //    TPZCompMeshTools::OptimizeBandwidth(fMixedFluxPressureCmesh);
    
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
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMMultiphase * mat = new TRMMultiphase(rock_id,dim);
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
    
    int64_t nel = fMonolithicMultiphaseCmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
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
    for (int64_t icel=0; icel < fMixedFluxPressureCmesh->NElements(); icel++) {
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
        int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
        int bc_id = 0;
        std::pair< int, TPZFunction<REAL> * > bc_item;
        TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
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
        int n_boundauries = this->SimulationData()->RawData()->fGammaIds.size();
        int bc_id = 0;
        std::pair< int, TPZFunction<REAL> * > bc_item;
        TPZVec< std::pair< int, TPZFunction<REAL> * > > bc;
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
    int rock_id = 0;
    for (int i = 0; i < n_rocks; i++) {
        rock_id = this->SimulationData()->RawData()->fOmegaIds[i];
        TRMPhaseTransport * mat = new TRMPhaseTransport(rock_id,dim);
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
            boundary_bc->SetNumLoadCases(1);
            boundary_bc->SetMaterial(matint);
            boundary_bc->SetId(bc_id);
            boundary_bc->SetType(bc_item.first);
            TPZAutoPointer<TPZFunction<STATE> > boundary_data = bc_item.second;
            boundary_bc->SetTimedependentBCForcingFunction(0,boundary_data); // @Omar:: Modified for multiple rock materials and set the polynomial order of the functions
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
    int64_t nel = fTransportMesh->ElementVec().NElements();
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
    for (int64_t el = 0; el<nel; el++) {
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
    
    int64_t last_node = fGeoMesh->NNodes() - 1;
    int64_t last_element = fGeoMesh->NElements() - 1;
    int64_t node_id = fGeoMesh->NodeVec()[last_node].Id();
    int64_t element_id = fGeoMesh->Element(last_element)->Id();
    const std::string name("GID Reservoir geometry");
    fGeoMesh->SetName(name);
    fGeoMesh->SetMaxNodeId(node_id);
    fGeoMesh->SetMaxElementId(element_id);
    
}

/** @brief Create a reservoir-box geometry with cylindrical wells */
void TRMSpaceOdissey::CreateGeometricGmshMesh(std::string &grid){
    
    TRMGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    fGeoMesh = Geometry.GeometricGmshMesh(grid);
    
    const std::string name("Reservoir with cylindrical wells");
    fGeoMesh->SetName(name);
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
    TPZAutoPointer<TPZFunction<REAL> > ParFuncZ = new TPZDummyFunction<REAL>(ParametricfunctionZ);
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
    
    int64_t last_node = fGeoMesh->NNodes() - 1;
    int64_t last_element = fGeoMesh->NElements() - 1;
    int64_t node_id = fGeoMesh->NodeVec()[last_node].Id();
    int64_t element_id = fGeoMesh->Element(last_element)->Id();
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
    
    TPZVec<int64_t> Topology(1,0);
    int elid=0;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,rock,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom0D(GeoMesh0D);
    TPZAutoPointer<TPZFunction<REAL> > ParFuncX = new TPZDummyFunction<REAL>(ParametricfunctionX);
    CreateGridFrom0D.SetParametricFunction(ParFuncX);
    CreateGridFrom0D.SetFrontBackMatId(bc_W,bc_E);
    
    dt = dx[0];
    n = int(dx[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh1D = CreateGridFrom0D.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom1D(GeoMesh1D);
    TPZAutoPointer<TPZFunction<REAL> > ParFuncY = new TPZDummyFunction<REAL>(ParametricfunctionY);
    CreateGridFrom1D.SetParametricFunction(ParFuncY);
    CreateGridFrom1D.SetFrontBackMatId(bc_S,bc_N);
    if(IsTetrahedronMeshQ){
        CreateGridFrom1D.SetTriangleExtrusion();
    }
    
    
    dt = dy[0];
    n = int(dy[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fGeoMesh = CreateGridFrom1D.ComputeExtrusion(t, dt, n);
    

    int64_t last_node = fGeoMesh->NNodes() - 1;
    int64_t last_element = fGeoMesh->NElements() - 1;
    int64_t node_id = fGeoMesh->NodeVec()[last_node].Id();
    int64_t element_id = fGeoMesh->Element(last_element)->Id();
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
    
    TPZVec<int64_t> Topology(1,0);
    int elid=0;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,rock,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom0D(GeoMesh0D);
    TPZAutoPointer<TPZFunction<REAL> > ParFuncX = new TPZDummyFunction<REAL>(ParametricfunctionX);
    CreateGridFrom0D.SetParametricFunction(ParFuncX);
    CreateGridFrom0D.SetFrontBackMatId(bc_W,bc_E);
    
    dt = dx[0];
    n = int(dx[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh1D = CreateGridFrom0D.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom1D(GeoMesh1D);
    TPZAutoPointer<TPZFunction<REAL> > ParFuncY = new TPZDummyFunction<REAL>(ParametricfunctionY);
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
    TPZAutoPointer<TPZFunction<REAL> > ParFuncZ = new TPZDummyFunction<REAL>(ParametricfunctionZ);
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

    int64_t last_node = fGeoMesh->NNodes() - 1;
    int64_t last_element = fGeoMesh->NElements() - 1;
    int64_t node_id = fGeoMesh->NodeVec()[last_node].Id();
    int64_t element_id = fGeoMesh->Element(last_element)->Id();
    const std::string name("Reservoir box");
    fGeoMesh->SetName(name);
    fGeoMesh->SetMaxNodeId(node_id);
    fGeoMesh->SetMaxElementId(element_id);
    
}


void TRMSpaceOdissey::ParametricfunctionX(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = par[0];
    X[1] = 0.0*sin(0.01*par[0]);
    X[2] = 0.0;
}

void TRMSpaceOdissey::ParametricfunctionY(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void TRMSpaceOdissey::ParametricfunctionZ(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}

void TRMSpaceOdissey::CElemtentRefinement(TPZCompMesh  *cmesh, int element_index)
{
    
    TPZVec<int64_t > subindex;
    int64_t nel = cmesh->ElementVec().NElements();
    for(int64_t el=0; el < nel; el++){
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        int64_t ind = compEl->Index();
        if(ind==element_index){
            compEl->Divide(element_index, subindex, 1);
        }
    }
}

void TRMSpaceOdissey::CMeshRefinement(TPZCompMesh  *cmesh, int ndiv)
{
    
    TPZVec<int64_t > subindex;
    for (int64_t iref = 0; iref < ndiv; iref++) {
        TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
        int64_t nel = elvec.NElements();
        for(int64_t el=0; el < nel; el++){
            TPZCompEl * compEl = elvec[el];
            if(!compEl) continue;
            int64_t ind = compEl->Index();
            compEl->Divide(ind, subindex, 0);
        }
    }
}

/** @brief Apply uniform refinement on the Geometric mesh */
void TRMSpaceOdissey::UniformRefinement(int n_ref){
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        int64_t n = fGeoMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = fGeoMesh->ElementVec() [i];
            if (gel->Dimension() != 0) gel->Divide (sons);
        }//for i
    }//ref
    fGeoMesh->BuildConnectivity();
}

/** @brief Apply uniform refinement at specific material id */
void TRMSpaceOdissey::UniformRefinement_at_MaterialId(int n_ref, int mat_id){
    
    int dim = fGeoMesh->Dimension();
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        int64_t n = fGeoMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
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
        int64_t n = fGeoMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
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
        int64_t n = fGeoMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
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

/** @brief Apply uniform refinement on the Geometric mesh */
void TRMSpaceOdissey::UniformRefinement_at_Father(int n_ref, int father_index){
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = fGeoMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = fGeoMesh->ElementVec() [i];
            if(gel->HasSubElement()){
                continue;
            }
            if(gel->Index() == father_index ){
                if (gel->Dimension() != 0) gel->Divide (filhos);
            }

        }//for i
    }//ref
    
//    int side;
//    TPZStack<TPZGeoElSide> subel;
//    gel->GetSubElements2(<#int side#>, subel)
    
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

static void IncludeNeighbours(TPZAutoPointer<TPZCompMesh> cmesh, int64_t index, std::map<int64_t,int> &extended)
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
        int64_t celindex = celside.Element()->Index();
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
                        int64_t neighindex = neighbour.Element()->Reference()->Index();
                        std::cout << "Including index " << neighindex << " order " << extended[index] << std::endl;
                        extended[neighindex] = extended[index];
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
        }
        
    }

}

void TRMSpaceOdissey::ModifyElementOrders(std::map<int64_t,int> &elorders)
{
    // settle the orders of the pressure elements first
    this->fMixedFluxPressureCmesh->Reference()->ResetReference();
    for (std::map<int64_t,int>::iterator it = elorders.begin(); it != elorders.end(); it++) {
        int64_t elindex = it->first;
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
    for (std::map<int64_t,int>::iterator it = elorders.begin(); it != elorders.end(); it++) {
        int64_t elindex = it->first;
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
    int64_t nelem = fGeoMesh->NElements();
    int ilocal = 0;
    for (int64_t el=0; el < nelem; el++) {
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
    
    TPZManVector<std::map<int64_t,int>,2> contemplated(2);
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
    
    TPZManVector<std::map<int64_t,int>,2> original(2);
    //TPZManVector<std::map<int64_t,int>,2> original(contemplated[0]);
    //std::map<int64_t,int> original2(contemplated[1]);
    
    for (int i = 0; i < contemplated.size(); i++) {
        original[i] = contemplated[i];
    }
    
    for (int i = 0; i < original.size(); i++) {
        for (std::map<int64_t,int>::iterator it = original[i].begin(); it != original[i].end(); it++) {
            IncludeNeighbours(this->fMixedFluxPressureCmesh, it->first, contemplated[i]);
        }
    }
    for (int i = 0 ; i < contemplated.size() ; i++){
        ModifyElementOrders(contemplated[i]);
    }
}


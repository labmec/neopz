//
//  TRMBuildTransfers.cpp
//  PZ
//
//  Created by Omar on 10/27/15.
//
//


#include "TRMBuildTransfers.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


/** @brief Default constructor */
TRMBuildTransfers::TRMBuildTransfers(){
    
    fSimulationData = NULL;

}

/** @brief Default desconstructor */
TRMBuildTransfers::~TRMBuildTransfers(){

}

/** @brief Copy constructor $ */
TRMBuildTransfers::TRMBuildTransfers(const TRMBuildTransfers &copy)
{
    fSimulationData = copy.fSimulationData;
}

/** @brief Copy assignemnt operator $ */
TRMBuildTransfers & TRMBuildTransfers::operator=(const TRMBuildTransfers &other)
{
    if (this != & other) // prevent self-assignment
    {
        fSimulationData = other.fSimulationData;
    }
    return *this;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Matrices Initialization Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TRMBuildTransfers::Initialize_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = cmesh_multiphysics->Reference()->Dimension(); // vectorial
    int64_t element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fu_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(nel);
    
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(intel->Dimension() < n_var_dim){
            // there is boundary elements for normal flux where it is a scalar variable
//            mf_cel->GetMemoryIndices(int_point_indexes);
//            this->ElementDofIndexes(intel, dof_indexes);
//            fu_dof_scatter[element_index] = dof_indexes;
            blocks_dimensions[element_index].first = 0;
            blocks_dimensions[element_index].second = 0;
            fu_dof_scatter[element_index] = dof_indexes;
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);
        fu_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
        fu_dof_scatter[element_index] = dof_indexes;
    }
    
    // Initialize the matrix
    fu_To_Mixed.Initialize(blocks_dimensions);
    
}


void TRMBuildTransfers::Initialize_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    int64_t element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fp_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(nel);
    
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(!intel){
            // there is no boundary elements for pressure
            blocks_dimensions[element_index].first = 0*n_var_dim;
            blocks_dimensions[element_index].second = 0;
            fp_dof_scatter[element_index] = dof_indexes;
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);
        fp_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
        fp_dof_scatter[element_index] = dof_indexes;
    }
    
    // Initialize the matrix
    fp_To_Mixed.Initialize(blocks_dimensions);
    
    
}


void TRMBuildTransfers::Initialize_s_To_Transport(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    int64_t element_index = 0;
    int dimension = cmesh_multiphysics->Reference()->Dimension();
    // Compute destination index scatter by element (Omega and Gamma)
    if (!mesh_index) {
        fsa_dof_scatter.Resize(nel);
    }
    else{
        fsb_dof_scatter.Resize(nel);
    }

    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(nel);
    
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
        if(gel->HasSubElement()){
            continue;
        }
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
       
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);        
        
#ifdef PZDEBUG
        if(!mf_cel)
        {
            if(!mf_int_cel){
                DebugStop();
            }
        }
#endif
        if (mf_int_cel) {
            continue;
        }
        
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(intel->Dimension() < dimension){
            // there is no boundary elements for saturation
            blocks_dimensions[element_index].first = 0*n_var_dim;
            blocks_dimensions[element_index].second = 0;
            if (!mesh_index) {
                fsa_dof_scatter[element_index] = dof_indexes;
            }
            else{
                fsb_dof_scatter[element_index] = dof_indexes;
            }
            
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);

        if (!mesh_index) {
            fsa_dof_scatter[element_index] = dof_indexes;
        }
        else{
            fsb_dof_scatter[element_index] = dof_indexes;
        }
        
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
    }
    
    // Initialize the matrix
    fs_To_Transport.Initialize(blocks_dimensions);
    
    
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Matrices Filling Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TRMBuildTransfers::Fill_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_u_To_Mixed(cmesh_multiphysics, mesh_index);
    
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = cmesh_multiphysics->Reference()->Dimension();; // vector
    int64_t element_index = 0;
    
    TPZMaterialData data;
    
    std::pair<int64_t, int64_t> block_dim;
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        dof_indexes = fu_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block;
        
        if(intel->Dimension() < n_var_dim){ // lower dimensional elements
            block.Resize(block_dim.first,block_dim.second);        }
        else{
            block.Resize(block_dim.first*n_var_dim,block_dim.second);
        }
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            // Get the vectorial phi
            intel->Shape(qsi, phi, dphidxi);
            intel->InitMaterialData(data);
            intel->ComputeRequiredData(data,qsi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    int vector_index = data.fVecShapeIndex[jp].first;
                    int shape_index = data.fVecShapeIndex[jp].second;
                    block(ip*n_var_dim+id,jp) = phi(shape_index,0)*data.fDeformedDirections(id,vector_index);
                }
            }
            
        }
        
        fu_To_Mixed.SetBlock(element_index, block);
        
    }
    
    return;
}


void TRMBuildTransfers::Fill_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_p_To_Mixed(cmesh_multiphysics, mesh_index);
    
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    int64_t element_index = 0;
    
    std::pair<int64_t, int64_t> block_dim;
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(!intel){
            continue;
        }
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        dof_indexes = fp_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block(block_dim.first*n_var_dim,block_dim.second);
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            intel->Shape(qsi, phi, dphidxi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    block(ip+id*n_var_dim,jp) = phi(jp,0);
                }
            }
            
        }
        
        fp_To_Mixed.SetBlock(element_index, block);
        
    }
    
    return;
}



void TRMBuildTransfers::Fill_s_To_Transport(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_s_To_Transport(cmesh_multiphysics, mesh_index);
    
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    int64_t element_index = 0;
    int dimension = cmesh_multiphysics->Reference()->Dimension();
    std::pair<int64_t, int64_t> block_dim;
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        
        TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (mf_int_cel) {
            continue;
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(intel->Dimension() < dimension){
            continue;
        }
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        if (!mesh_index) {
            dof_indexes = fsa_dof_scatter[element_index];
        }
        else{
            dof_indexes = fsb_dof_scatter[element_index];
        }
        
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block(block_dim.first*n_var_dim,block_dim.second);
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            intel->Shape(qsi, phi, dphidxi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    block(ip+id*n_var_dim,jp) = phi(jp,0);
                }
            }
            
        }
        
        fs_To_Transport.SetBlock(element_index, block);
        
    }
    return;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Transfer Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TRMBuildTransfers::u_To_Mixed_Memory(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_multiphysics){
    

#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    int dim = cmesh_flux->Dimension();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZMaterial> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZMaterial> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterFlux(fu_To_Mixed.Cols(),1,0.0);
    int64_t pos = 0;
    for (int el = 0; el < nel; el++) {
        for(int ip = 0; ip < fu_dof_scatter[el].size(); ip++) {
            ScatterFlux(pos,0) = cmesh_flux->Solution()(fu_dof_scatter[el][ip],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> Flux_at_intpoints;
    fu_To_Mixed.Multiply(ScatterFlux,Flux_at_intpoints);
    // Trasnfering values
    TPZManVector<STATE,3> u(dim,0.0);
    for(int64_t i = 0; i <  np_cmesh; i++){
        for (int id = 0; id < dim ; id++) {
            u[id]= Flux_at_intpoints(i*dim+id,0);
        }

        if(fSimulationData->IsCurrentStateQ()){
            associated_material->GetMemory()[i].Set_u(u);
        }
        else{
            associated_material->GetMemory()[i].Set_u_n(u);
        }
        
    }
    
}

void TRMBuildTransfers::p_To_Mixed_Memory(TPZCompMesh * cmesh_pressure, TPZCompMesh * cmesh_multiphysics){

    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];

    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZMaterial> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZMaterial> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterPressure(fp_To_Mixed.Cols(),1,0.0);
    int64_t pos = 0;
    for (int el = 0; el < nel; el++) {
        for(int ip = 0; ip < fp_dof_scatter[el].size(); ip++) {
            ScatterPressure(pos,0) = cmesh_pressure->Solution()(fp_dof_scatter[el][ip],0);
            pos++;
        }
    }
    
    // Step two
    TPZFNMatrix<30,STATE> Pressure_at_intpoints;
    fp_To_Mixed.Multiply(ScatterPressure,Pressure_at_intpoints);
    // Trasnfering values
    for(int64_t i = 0; i <  np_cmesh; i++){
        if(fSimulationData->IsCurrentStateQ()){
            associated_material->GetMemory()[i].Set_p_n(Pressure_at_intpoints(i,0));
        }
        else{
            associated_material->GetMemory()[i].Set_p(Pressure_at_intpoints(i,0));
        }
    }
    
}

void TRMBuildTransfers::s_To_Transport_Memory(TPZCompMesh * cmesh_saturation, TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
//    if(fSimulationData->IsThreePhaseQ()){
//        DebugStop();
//    }
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZMaterial> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZMaterial> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterSaturation(fs_To_Transport.Cols(),1,0.0);
    int64_t pos = 0;
    for (int el = 0; el < nel; el++) {
        
        if (!mesh_index) {
            for(int ip = 0; ip < fsa_dof_scatter[el].size(); ip++) {
                ScatterSaturation(pos,0) = cmesh_saturation->Solution()(fsa_dof_scatter[el][ip],0);
                pos++;
            }
        }
        else{
            for(int ip = 0; ip < fsb_dof_scatter[el].size(); ip++) {
                ScatterSaturation(pos,0) = cmesh_saturation->Solution()(fsb_dof_scatter[el][ip],0);
                pos++;
            }
        }
        
    }
    
    // Step two
    TPZFMatrix<STATE> Saturation_at_intpoints;
    fs_To_Transport.Multiply(ScatterSaturation,Saturation_at_intpoints);
    // Trasnfering values
    for(int64_t i = 0; i <  np_cmesh; i++){
        
        if(!mesh_index){
            if(fSimulationData->IsCurrentStateQ()){
                associated_material->GetMemory()[i].Set_sa_n(Saturation_at_intpoints(i,0));
            }
            else{
                associated_material->GetMemory()[i].Set_sa(Saturation_at_intpoints(i,0));
            }
        }
        else{
            if(fSimulationData->IsCurrentStateQ()){
                associated_material->GetMemory()[i].Set_sb_n(Saturation_at_intpoints(i,0));
            }
            else{
                associated_material->GetMemory()[i].Set_sb(Saturation_at_intpoints(i,0));
            }
        }
    
    }
    
}


/** @brief Reciprocal (mixed <-> transpor) transfer average quantities to integration points of multiphysics meshes over volumetric elements */
void TRMBuildTransfers::Reciprocal_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_trans){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_cindexes.size() == 0 ) {
        DebugStop();
    }

    if (!cmesh_mf_mixed || !cmesh_mf_trans) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZMaterial> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZMaterial> *>(mixed_material);
    
    TPZMaterial * trans_material = cmesh_mf_trans->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZMaterial> * trans_memory = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZMaterial> *>(trans_material);
    
    TPZManVector<int64_t,30> p_point_indexes;
    TPZManVector<int64_t,30> s_point_indexes;
    int64_t nvolumes = fmixed_transport_cindexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
    TPZGeoEl  * gel = geometry->Element(fmixed_transport_cindexes[ivol].first);
    TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_cindexes[ivol].second.first);
    TPZCompEl * trans_cel = cmesh_mf_trans->Element(fmixed_transport_cindexes[ivol].second.second);
        
#ifdef PZDEBUG
        if (!mixed_cel || !trans_cel || !gel) {
            DebugStop();
        }
#endif

        REAL element_measure = DimensionalMeasure(mixed_cel->Reference());
    
        GlobalPointIndexes(mixed_cel, p_point_indexes);
        GlobalPointIndexes(trans_cel, s_point_indexes);
        
        TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        TPZMultiphysicsElement * mf_trans_cel = dynamic_cast<TPZMultiphysicsElement * >(trans_cel);
        
#ifdef PZDEBUG
        if (!mf_mixed_cel || !mf_trans_cel) {
            DebugStop();
        }
#endif
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
        int np_mixed_cel = int_points_mixed.NPoints();
        
        const TPZIntPoints & int_points_trans = mf_trans_cel->GetIntegrationRule();
        int np_trans_cel = int_points_trans.NPoints();
        
#ifdef PZDEBUG
        if (np_mixed_cel != p_point_indexes.size() || np_trans_cel != s_point_indexes.size()) {
            DebugStop();
        }
#endif
        
        REAL w;
        TPZManVector<REAL,3> triplet(3,0.0);
        
        REAL detjac;
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;

        REAL p_avg      = 0.0;
        REAL p_avg_n    = 0.0;

        // Integrating pressure
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            int_points_mixed.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            p_avg_n += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure;
            p_avg +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure;
            
        }

        REAL sa      = 0.0;
        REAL sa_n    = 0.0;
        REAL sb      = 0.0;
        REAL sb_n    = 0.0;
        
        // Integrating Saturation
        for (int ip = 0; ip < np_trans_cel; ip++) {
            
            int_points_trans.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            sa_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa_n()/element_measure;
            sb_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb_n()/element_measure;
            sa +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa()/element_measure;
            sb +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb()/element_measure;

        }
        
//        std::cout << "p_avg_n = "<< p_avg_n << std::endl;
//        std::cout << "p_avg = "<< p_avg << std::endl;
//        
//        std::cout << "sa_n = "<< sa_n << std::endl;
//        std::cout << "sa = "<< sa << std::endl;
        
        // Inserting average pressure and saturation in mixed memory
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            if (fSimulationData->IsCurrentStateQ()) {
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa_n(sa_n);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb_n(sb_n);
            }
            else{
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa(sa);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb(sb);
            }

        }

        // Inserting average pressure in transport memory
        for (int ip = 0; ip < np_trans_cel; ip++) {
            
            if (fSimulationData->IsCurrentStateQ()) {
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sa_n);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sb_n);
            }
            else{
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg(p_avg);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sa);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sb);
            }

        }
        
    }
    
    return;
    
}

/** @brief Transfer average pressure to integration points of multiphysics mixed meshes over volumetric elements */
void TRMBuildTransfers::p_avg_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_cindexes.size() == 0 ) {
        DebugStop();
    }
    
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZMaterial> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZMaterial> *>(mixed_material);
    
    TPZManVector<int64_t,30> p_point_indexes;
    int64_t nvolumes = fmixed_transport_cindexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
        TPZGeoEl  * gel = geometry->Element(fmixed_transport_cindexes[ivol].first);
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_cindexes[ivol].second.first);
        
#ifdef PZDEBUG
        if (!mixed_cel || !gel) {
            DebugStop();
        }
#endif
        
        REAL element_measure = DimensionalMeasure(mixed_cel->Reference());
        
        GlobalPointIndexes(mixed_cel, p_point_indexes);
        TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        
#ifdef PZDEBUG
        if (!mf_mixed_cel) {
            DebugStop();
        }
#endif
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
        int np_mixed_cel = int_points_mixed.NPoints();
        
        
#ifdef PZDEBUG
        if (np_mixed_cel != p_point_indexes.size()) {
            DebugStop();
        }
#endif
        
        REAL w;
        TPZManVector<REAL,3> triplet(3,0.0);
        
        REAL detjac;
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;
        
        REAL p_avg      = 0.0;
        REAL p_avg_n    = 0.0;
        
        // Integrating pressure
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            int_points_mixed.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            p_avg_n +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure;
            p_avg += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure;
            
        }
        
        // Inserting average pressure and saturation in mixed memory
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            if (fSimulationData->IsCurrentStateQ()) {
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
            }
            else{
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);

            }
            
        }
        
    }
        
}



/** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Initialize_un_To_Transport(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    
#ifdef PZDEBUG
    if (!flux_mesh || !transport_mesh) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    int mesh_index = 0;
    
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<int64_t,10> dof_indexes;
    
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;

    
    TPZManVector<int64_t> indices;
    std::pair<int64_t, int64_t> duplet;
    TPZManVector<int,10> face_sides;
    int64_t face_index;
    int64_t n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
        fun_dof_scatter_Gamma.Resize(n_interfaces);
        fun_To_Transport_Gamma.Resize(0, 0);
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
        fun_dof_scatter_gamma.Resize(n_interfaces);
        fun_To_Transport_gamma.Resize(0, 0);
    }

    
    // Block size structue (Gamma or gamma (Inner element interfaces))
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(n_interfaces);
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {


        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }
        
        face_gel = geometry->Element(face_index);
        
        if (!face_gel) {
            DebugStop();
        }
        
        int64_t left_geo_index     = duplet.first;
        int64_t right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        this->ComputeFaceIndex(left_gel,face_sides);
        
//        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
//        
//        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace *> (left_cel);
        
        
        int face_side     = -1;
        int64_t connect_index = -1;

        if(!IdentifyFace(face_side,left_gel,face_gel)){
            std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
            DebugStop();
        }
        else{
            
            for(int ic = 0; ic < intel_vol->NConnects(); ic++){
                connect_index = intel_vol->SideConnectLocId(ic, face_side);
                if (connect_index != -1) {
                    break;
                }
            }
            
            if (connect_index == -1) {
                std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
                DebugStop();
            }
        }

        this->ElementDofFaceIndexes(connect_index, mf_cel, dof_indexes);
        
        int nshapes = left_cel->Connect(connect_index).NShape();
        
        blocks_dimensions[k_face].first = 1;
        blocks_dimensions[k_face].second = nshapes;
        if (IsBoundaryQ) {
            fun_dof_scatter_Gamma[k_face] = dof_indexes;
        }
        else{
            fun_dof_scatter_gamma[k_face] = dof_indexes;
        }

        
    }
    
    // Initialize the matrix
    
    if (IsBoundaryQ) {
        fun_To_Transport_Gamma.Initialize(blocks_dimensions);
    }
    else{
        fun_To_Transport_gamma.Initialize(blocks_dimensions);
    }
    
    
}

/** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Fill_un_To_Transport(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    

#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif

    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_un_To_Transport(flux_mesh,transport_mesh,IsBoundaryQ);
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif

#ifdef USING_BOOST
    std::cout  << "Time for Matrix Initialization " << (t2-t1) << std::endl;
#endif
    int mesh_index = 0;
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<int64_t,10> dof_indexes;
    
    TPZCompEl * cel_face;
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    
    TPZManVector<int64_t> indices;
    std::pair<int64_t, int64_t> duplet;
    TPZManVector<int,10> face_sides;
    TPZFMatrix<REAL> normals;
    int64_t face_index;
    int64_t n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    TPZFNMatrix<100,double> block;
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }
        
//        cel_face    = transport_mesh->Element(face_index);
//        
//        if (!cel_face) {
//            DebugStop();
//        }

        face_gel = geometry->Element(face_index);

        
        if (!face_gel) {
            DebugStop();
        }
        
        int64_t left_geo_index     = duplet.first;
        int64_t right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        if(left_gel->HasSubElement() && right_gel->HasSubElement() && face_gel->HasSubElement()){
            DebugStop();
        }
        
//        // Computing face connect index associated to the element face being integrated
//        // The method is really simple the inteface is detected by connect shared between hdiv volumetric elements
//        this->ComputeFaceNormals(left_gel,face_sides,normals);

        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
//        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace *> (left_cel);
        
        int nfaces = face_sides.size();
        
        int face_side       = -1;
        int connect_index   = -1;
        
        if(!IdentifyFace(face_side,left_gel,face_gel)){
            std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
            DebugStop();
        }
        else{
            
            for(int ic = 0; ic < intel_vol->NConnects(); ic++){
                connect_index = intel_vol->SideConnectLocId(ic, face_side);
                if (connect_index != -1) {
                    break;
                }
            }
            
            if (connect_index == -1) {
                std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
                DebugStop();
            }
        }
        
        this->ElementDofFaceIndexes(connect_index, mf_cel, dof_indexes);
        TPZIntPoints *int_points   = left_gel->CreateSideIntegrationRule(face_side, left_cel->GetgOrder());
        
        int npoints = int_points->NPoints();
        int nshapes = left_cel->Connect(connect_index).NShape();
  
#ifdef PZDEBUG
        if (IsBoundaryQ) {
            if (1 != fun_To_Transport_Gamma.GetSizeofBlock(k_face).first || nshapes != fun_To_Transport_Gamma.GetSizeofBlock(k_face).second){
                DebugStop();
            }
        }
        else{
            if (1 != fun_To_Transport_gamma.GetSizeofBlock(k_face).first || nshapes != fun_To_Transport_gamma.GetSizeofBlock(k_face).second){
                DebugStop();
            }
            
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(nshapes,1,0.0);
        int el_dim = face_gel->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,nshapes,0.0);
        TPZFNMatrix<50,double> block(npoints,nshapes);
        TPZFNMatrix<50,double> block_integral(1,nshapes,0.0);
        
        REAL w;
        TPZManVector<STATE,2> par_duplet(el_dim,0.0);
        REAL ElementMeasure   = DimensionalMeasure(face_gel);

        REAL detjac;
        TPZFMatrix<REAL> jac,axes,jacinv;

        for (int ip = 0; ip < npoints; ip++) {
         
            // Get the vectorial phi
            int_points->Point(ip, par_duplet, w);
            face_gel->Jacobian(par_duplet, jac, axes, detjac, jacinv);
            intel_vol->SideShapeFunction(face_side, par_duplet, phi, dphidxi);
            
            
            for (int jp = 0; jp < nshapes; jp++) {
                block_integral(0,jp) +=  w * detjac * phi(jp,0)/ElementMeasure;
            }
            
        }
        
        if (IsBoundaryQ) {
            fun_To_Transport_Gamma.SetBlock(k_face, block_integral);
        }
        else{
            fun_To_Transport_gamma.SetBlock(k_face, block_integral);
        }
        

        
    }
    
    
}



/** @brief Transfer normal fluxes to integration points of transport meshes */
void TRMBuildTransfers::un_To_Transport_Mesh(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_transport, bool IsBoundaryQ){
  
#ifdef PZDEBUG
    if (!cmesh_flux || !cmesh_transport) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
#ifdef PZDEBUG
    if (this->SimulationData()->RawData()->fOmegaIds[0] != 4) {
        DebugStop();
    }
    
#endif
    
    TPZManVector<int64_t,30>  point_index_trans;
    TPZManVector<int64_t,30>  point_index_l;
    TPZManVector<int64_t,30>  point_index_r;
    
    REAL p_avg_n_l = -1.0;
    REAL p_avg_n_r = -1.0;
    
    cmesh_transport->LoadReferences();
    TPZGeoMesh * geometry = cmesh_transport->Reference();
    int64_t n_interfaces;
    int dimension = geometry->Dimension();
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    int rock_id = this->SimulationData()->RawData()->fOmegaIds[0];
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * rock_material = cmesh_flux->FindMaterial(rock_id);
    TPZMatWithMem<TRMMemory,TPZMaterial>  * material_mixe_mem = dynamic_cast<TPZMatWithMem<TRMMemory,TPZMaterial> *>(rock_material);

    
    if (IsBoundaryQ) {
        
        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int nbc = this->SimulationData()->RawData()->fGammaIds.size();
        
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            int material_id = this->SimulationData()->RawData()->fGammaIds[ibc];
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
            
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  * material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> *>(material);
            
            if (!material_bc_mem) {
                DebugStop();
            }
            
            // Step one
            TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_Gamma.Cols(),1,0.0);
            int64_t pos = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                for(int iflux = 0; iflux < fun_dof_scatter_Gamma[iface].size(); iflux++) {
                    ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_Gamma[iface][iflux],0);
                    pos++;
                }
            }
            
            // Step two
            TPZFMatrix<STATE> un_at_intpoints;
            fun_To_Transport_Gamma.Multiply(ScatterFluxes,un_at_intpoints);
            
//            un_at_intpoints.Print("u_Gamma = ");
            
            // Step three
            // Trasnfering integrated normal fluxes values
            int counter = 0;
            int i = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                TPZGeoEl *gel    = geometry->Element(finterface_g_indexes_Gamma[iface]);
                TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_Gamma[iface].first);
                TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_Gamma[iface].second);

                
#ifdef PZDEBUG
                if (!gel || !gel_l || !gel_r) {
                    DebugStop();
                }
                
                if (gel_l->Dimension() != dimension) {
                    DebugStop();
                }
#endif
                TPZCompEl * mixed_cel_l = gel_l->Reference();
                TPZCompEl * mixed_cel_r = gel_r->Reference();

                
                if(gel->MaterialId() != material_bc_mem->Id()){
                    counter++;
                    continue;
                }
                
                
#ifdef PZDEBUG
                if (!mixed_cel_l || !mixed_cel_r) {
                    DebugStop();
                }
                
#endif
                
                GlobalPointIndexes(mixed_cel_l, point_index_l);
                
                p_avg_n_l = material_mixe_mem->GetMemory()[point_index_l[0]].p_avg_n();
                material_bc_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
                material_bc_mem->GetMemory()[i].Set_un(un_at_intpoints(counter,0));
                i++;
                counter++;
                
            }
        }

        
    }
    else{
        
        
        int material_id = this->SimulationData()->InterfacesMatId();        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
        
        TPZMatWithMem<TRMPhaseInterfaceMemory,TPZMaterial>  * material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZMaterial> *>(material);
        
        if (!material_mem) {
            DebugStop();
        }

        int np_cmesh = material_mem->GetMemory().NElements();

        // Step one
        TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_gamma.Cols(),1,0.0);
        int64_t pos = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {
            for(int iflux = 0; iflux < fun_dof_scatter_gamma[iface].size(); iflux++) {
                ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_gamma[iface][iflux],0);
                pos++;
            }
        }

        // Step two
        TPZFMatrix<STATE> un_at_intpoints;
        fun_To_Transport_gamma.Multiply(ScatterFluxes,un_at_intpoints);
        
//        un_at_intpoints.Print("u_gamma = ");        
        
        // Step three
        // Trasnfering integrated normal fluxes values
        for(int64_t i = 0; i < np_cmesh; i++){
            material_mem->GetMemory()[i].Set_un(un_at_intpoints(i,0));
        }

        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int i = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {

            TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_gamma[iface].first);
            TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_gamma[iface].second);
            
            
#ifdef PZDEBUG
            if (!gel_l || !gel_r) {
                DebugStop();
            }
            
            if (gel_l->Dimension() != dimension || gel_r->Dimension() != dimension) {
                DebugStop();
            }
#endif
            
            TPZCompEl * mixed_cel_l = gel_l->Reference();
            TPZCompEl * mixed_cel_r = gel_r->Reference();
            
            
#ifdef PZDEBUG
            if (!mixed_cel_l || !mixed_cel_r) {
                DebugStop();
            }
            
#endif

            GlobalPointIndexes(mixed_cel_l, point_index_l);
            GlobalPointIndexes(mixed_cel_r, point_index_r);
            
            p_avg_n_l = material_mixe_mem->GetMemory()[point_index_l[0]].p_avg_n();
            p_avg_n_r = material_mixe_mem->GetMemory()[point_index_r[0]].p_avg_n();
            
            material_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
            material_mem->GetMemory()[i].Set_p_avg_n_r(p_avg_n_r);
            i++;
            
        }
       
    }
    
    
    return;
    
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Utility Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief Get Global integration point indexes associaded  */
void TRMBuildTransfers::GlobalPointIndexes(TPZCompEl * cel, TPZManVector<int64_t,30> &int_point_indexes){
    
    TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
    
#ifdef PZDEBUG
    if(!mf_cel)
    {
        DebugStop();
    }
#endif
    
    mf_cel->GetMemoryIndices(int_point_indexes);
    
}

/** @brief Get Global integration point indexes associaded  */
void TRMBuildTransfers::GlobalPointIndexesInterface(TPZCompEl * int_cel, TPZManVector<int64_t,30> &int_point_indexes){
    
    TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(int_cel);
    
#ifdef PZDEBUG
    if(!mf_int_cel)
    {
        DebugStop();
    }
#endif
    
    mf_int_cel->GetMemoryIndices(int_point_indexes);
    
}

bool TRMBuildTransfers::IdentifyFace(int &side, TPZGeoEl * vol, TPZGeoEl * face){
    
    int volu_nsides = vol->NSides();
    int face_nsides = face->NSides();
    side = -1;
    TPZGeoElSide face_itself =  face->Neighbour(face_nsides-1);
    bool IsMembershipQ = false;
    
    for (int iside = 0; iside < volu_nsides; iside++) {
        IsMembershipQ = bool(vol->NeighbourExists(iside, face_itself));
        if (IsMembershipQ) {
            side = iside;
            break;
        }
    }
    
    TPZGeoElSide vol_itself =  vol->Neighbour(volu_nsides-1);
    
    if(!IsMembershipQ){
        TPZGeoElSide neigh = face->Neighbour(face_nsides-1);
        
        if(!neigh.Element()){
            DebugStop();
        }
        
        TPZGeoElSide neigh_father = neigh.Father2();
        
        if(!neigh_father.Element()){
            DebugStop();
        }
        bool IsNeighQ = false;
        for (int iside = vol_itself.Element()->NNodes(); iside < volu_nsides; iside++) {
 
            IsNeighQ = bool(vol_itself.Element()->NeighbourExists(iside, neigh_father));
            
            if (IsNeighQ) {
                side = iside;
                IsMembershipQ = IsNeighQ;
                break;
            }
        }
    }
    
    return IsMembershipQ;
}

/** @brief Compute indices associated to faces on 3D topologies */
void TRMBuildTransfers::ComputeFaceIndex(TPZGeoEl * gel , TPZVec<int> &sides){
    
    
    switch (gel->Type()) {
        case ECube:
        {
            int nfaces = 6;
            sides.Resize(nfaces);
            sides[0] = 20;
            sides[1] = 21;
            sides[2] = 22;
            sides[3] = 23;
            sides[4] = 24;
            sides[5] = 25;
            
        }
            break;
        case ETetraedro:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 10;
            sides[1] = 11;
            sides[2] = 12;
            sides[3] = 13;
            
        }
            break;
        case EQuadrilateral:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 4;
            sides[1] = 5;
            sides[2] = 6;
            sides[3] = 7;
            
        }
            break;
        case ETriangle:
        {
            int nfaces = 3;
            sides.Resize(nfaces);
            sides[0] = 3;
            sides[1] = 4;
            sides[2] = 5;
            
        }
            break;
        default:
        {
            std::cout << "Element not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compute sides associated to faces on 3D topologies */
void TRMBuildTransfers::ComputeFaceNormals(TPZGeoEl * gel , TPZVec<int> &sides, TPZFMatrix<STATE> &normals){
    
    //  @omar:: Just for linear mapping
    
    TPZFMatrix<REAL> mat_normals;
    TPZVec<int> v_sides;
    gel->ComputeNormals(mat_normals, v_sides);
    
    switch (gel->Type()) {
        case ECube:
        {
            int nfaces = 6;
            sides.Resize(nfaces);
            sides[0] = 20;
            sides[1] = 21;
            sides[2] = 22;
            sides[3] = 23;
            sides[4] = 24;
            sides[5] = 25;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
        break;
        case ETetraedro:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 10;
            sides[1] = 11;
            sides[2] = 12;
            sides[3] = 13;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
        break;
        case EQuadrilateral:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 4;
            sides[1] = 5;
            sides[2] = 6;
            sides[3] = 7;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
        case ETriangle:
        {
            int nfaces = 3;
            sides.Resize(nfaces);
            sides[0] = 3;
            sides[1] = 4;
            sides[2] = 5;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
            
        default:
        {
            std::cout << "Element not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compute left and right geometric element indexes */
void TRMBuildTransfers::ComputeLeftRight(TPZCompMesh * transport_mesh){
    
#ifdef PZDEBUG
    if (!transport_mesh) {
        std::cout << "There is no computational transport mesh, transport_mesh = Null." << std::endl;
        DebugStop();
    }
#endif
    
    int64_t nel = transport_mesh->NElements();
    int64_t face_index;
    std::pair <int64_t,int64_t> duplet;
    transport_mesh->LoadReferences();
    int dimension = transport_mesh->Reference()->Dimension();
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = transport_mesh->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsInterfaceElement * interface = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (!interface) {
            continue;
        }
        
        TPZCompEl * left_cel = interface->LeftElement();
        TPZCompEl * right_cel = interface->RightElement();
        
        if(!left_cel || !right_cel){
            DebugStop();
        }
        
        face_index  = interface->Reference()->Index();
        duplet      = std::make_pair(left_cel->Reference()->Index(), right_cel->Reference()->Index());
        
        if(left_cel->Dimension() != dimension ||  right_cel->Dimension() != dimension){
            
            fleft_right_g_indexes_Gamma.Push(duplet);
            finterface_g_indexes_Gamma.Push(face_index);
            continue;
        }
        
        fleft_right_g_indexes_gamma.Push(duplet);
        finterface_g_indexes_gamma.Push(face_index);
        
    }
    
#ifdef PZDEBUG
    if (finterface_g_indexes_Gamma.size() == 0) {
        DebugStop();
    }
    if (finterface_g_indexes_gamma.size() == 0) {
        std::cout << "Warning:: No inner interfaces were found" << std::endl;
    }
#endif
    
}

/** @brief Dimensionla Measure of the elemnt */
REAL TRMBuildTransfers::DimensionalMeasure(TPZGeoEl * gel){
    
#ifdef PZDEBUG
    if (!gel) {
        DebugStop();
    }
#endif
    REAL measure = 0.0;
    int order = 5;
    int element_itself  = gel->NSides() - 1;
    TPZIntPoints * int_points = gel->CreateSideIntegrationRule(element_itself, order);
    REAL detjac, w;
    TPZVec<REAL> par(gel->Dimension(),0.0);
    TPZFMatrix<REAL> jac;
    TPZFMatrix<REAL> axes;
    TPZFMatrix<REAL> jacinv;
    for (int i = 0; i < int_points->NPoints(); i++) {
        int_points->Point(i, par, w);
        gel->Jacobian(par, jac, axes, detjac, jacinv);
        measure += w * detjac;
    }
    
    return measure;
    
}


void TRMBuildTransfers::ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<int64_t> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<int64_t> index(0,0);
    int nconnect = intel->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        int64_t seqnumber = con.SequenceNumber();
        int64_t position = intel->Mesh()->Block().Position(seqnumber);
        int nshape = con.NShape();
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
        }
    }
    
    dof_indexes = index;
    return;
}

void TRMBuildTransfers::ElementDofFaceIndexes(int connect_index,TPZInterpolationSpace * &intel, TPZVec<int64_t> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<int64_t> index(0,0);
    TPZConnect  & con = intel->Connect(connect_index);
    int64_t seqnumber = con.SequenceNumber();
    int64_t position = intel->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}

void TRMBuildTransfers::ElementDofFaceIndexes(int connect_index, TPZMultiphysicsElement * &m_el, TPZVec<int64_t> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!m_el && connect_index > 4) {
        DebugStop();
    }
#endif
    

    TPZStack<int64_t> index(0,0);
    TPZConnect  & con = m_el->Connect(connect_index);
    int64_t seqnumber = con.SequenceNumber();
    int64_t position = m_el->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}


/** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
void TRMBuildTransfers::FillComputationalElPairs(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport){

    fmixed_transport_cindexes.Resize(0);
    
#ifdef PZDEBUG
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
    
    if (!fSimulationData->IsOnePhaseQ() && !cmesh_mf_transport) {
        DebugStop();
    }
    
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<int64_t, std::pair <int64_t,int64_t> > gel_indexes;
    
    for (int64_t i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() != dimension || gel->HasSubElement()) {
            continue;
        }
        gel_indexes.first = gel->Index();
        gel_indexes.second.first = -1;
        gel_indexes.second.second = -1;
        fmixed_transport_cindexes.Push(gel_indexes);

    }
    
    // counting volumetric elements
    int64_t nvol_elements = fmixed_transport_cindexes.size();
    fmixed_transport_cindexes.Resize(nvol_elements);
    
    // inserting mixed indexes
    cmesh_mf_mixed->LoadReferences();
    for(int64_t ivol = 0; ivol < nvol_elements; ivol++){

        TPZCompEl * mixed_cel = geometry->Element(fmixed_transport_cindexes[ivol].first)->Reference();
        
#ifdef PZDEBUG
        if (!mixed_cel) {
            //continue;
            DebugStop();
        }
#endif
        fmixed_transport_cindexes[ivol].second.first = mixed_cel->Index();
        
    }
    
    if(fSimulationData->IsOnePhaseQ()){
        return;
    }
    
    // inserting transport indexes
    cmesh_mf_transport->LoadReferences();
    for(int64_t ivol = 0; ivol < nvol_elements; ivol++){
        
        TPZCompEl * trans_cel = geometry->Element(fmixed_transport_cindexes[ivol].first)->Reference();
        
#ifdef PZDEBUG
        if (!trans_cel) {
            DebugStop();
        }
#endif
        fmixed_transport_cindexes[ivol].second.second = trans_cel->Index();
        
    }
    
}

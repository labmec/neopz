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


void TRMBuildTransfers::Initialize_u_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = cmesh_multiphysics->Reference()->Dimension(); // vectorial
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fu_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions(nel);
    
    
    for (long icel = 0; icel < nel; icel++) {
        
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
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
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


void TRMBuildTransfers::Initialize_p_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fp_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions(nel);
    
    
    for (long icel = 0; icel < nel; icel++) {
        
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
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
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


void TRMBuildTransfers::Initialize_s_To_Transport(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    long element_index = 0;
    int dimension = cmesh_multiphysics->Reference()->Dimension();
    // Compute destination index scatter by element (Omega and Gamma)
    if (!mesh_index) {
        fsa_dof_scatter.Resize(nel);
    }
    else{
        fsb_dof_scatter.Resize(nel);
    }

    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions(nel);
    
    
    for (long icel = 0; icel < nel; icel++) {
        
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
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
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

void TRMBuildTransfers::Fill_u_To_Mixed(TPZAutoPointer< TPZCompMesh > cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_u_To_Mixed(cmesh_multiphysics, mesh_index);
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = cmesh_multiphysics->Reference()->Dimension();; // vector
    long element_index = 0;
    
    TPZMaterialData data;
    
    std::pair<long, long> block_dim;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
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
                    block(ip*n_var_dim+id,jp) = phi(shape_index,0)*data.fNormalVec(id,vector_index);
                }
            }
            
        }
        
        fu_To_Mixed.SetBlock(element_index, block);
        
    }
    
    return;
}


void TRMBuildTransfers::Fill_p_To_Mixed(TPZAutoPointer< TPZCompMesh > cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_p_To_Mixed(cmesh_multiphysics, mesh_index);
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    long element_index = 0;
    
    std::pair<long, long> block_dim;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
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



void TRMBuildTransfers::Fill_s_To_Transport(TPZAutoPointer< TPZCompMesh > cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_s_To_Transport(cmesh_multiphysics, mesh_index);
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    long element_index = 0;
    int dimension = cmesh_multiphysics->Reference()->Dimension();
    std::pair<long, long> block_dim;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        
        TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (mf_int_cel) {
            continue;
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
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
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterFlux(fu_To_Mixed.Cols(),1,0.0);
    long pos = 0;
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
    for(long i = 0; i <  np_cmesh; i++){
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
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterPressure(fp_To_Mixed.Cols(),1,0.0);
    long pos = 0;
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
    for(long i = 0; i <  np_cmesh; i++){
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
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterSaturation(fs_To_Transport.Cols(),1,0.0);
    long pos = 0;
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
    for(long i = 0; i <  np_cmesh; i++){
        
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
    if ( fmixed_transport_indexes.size() == 0 ) {
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
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
    
    TPZMaterial * trans_material = cmesh_mf_trans->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * trans_memory = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(trans_material);
    
    TPZManVector<long,30> p_point_indexes;
    TPZManVector<long,30> s_point_indexes;
    long nvolumes = fmixed_transport_indexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
    TPZGeoEl  * gel = geometry->Element(fmixed_transport_indexes[ivol].first);
    TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_indexes[ivol].second.first);
    TPZCompEl * trans_cel = cmesh_mf_trans->Element(fmixed_transport_indexes[ivol].second.second);
        
#ifdef PZDEBUG
        if (!mixed_cel || !trans_cel || !gel) {
            DebugStop();
        }
#endif

        REAL element_measure = DimensionalMeasure(mixed_cel);
    
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
    if ( fmixed_transport_indexes.size() == 0 ) {
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
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
    
    TPZManVector<long,30> p_point_indexes;
    long nvolumes = fmixed_transport_indexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
        TPZGeoEl  * gel = geometry->Element(fmixed_transport_indexes[ivol].first);
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_indexes[ivol].second.first);
        
#ifdef PZDEBUG
        if (!mixed_cel || !gel) {
            DebugStop();
        }
#endif
        
        REAL element_measure = DimensionalMeasure(mixed_cel);
        
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
void TRMBuildTransfers::Initialize_un_To_Transport(TPZAutoPointer< TPZCompMesh> flux_mesh, TPZAutoPointer< TPZCompMesh> transport_mesh, bool IsBoundaryQ){
    
    
#ifdef PZDEBUG
    if (!flux_mesh || !transport_mesh) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZVec<long> dof_indexes;
    
    TPZCompEl * cel_face;
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;

    
    TPZManVector<long> indices;
    std::pair<long, long> duplet;
    TPZManVector<int,10> face_sides;
    long face_index;
    long n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_indexes_Gamma.size();
        fun_dof_scatter_Gamma.Resize(n_interfaces);
        fun_To_Transport_Gamma.Resize(0, 0);
    }
    else{
        n_interfaces = fleft_right_indexes_gamma.size();
        fun_dof_scatter_gamma.Resize(n_interfaces);
        fun_To_Transport_gamma.Resize(0, 0);
    }

    
    // Block size structue (Gamma or gamma (Inner element interfaces))
    TPZVec< std::pair<long, long> > blocks_dimensions(n_interfaces);
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {


        if (IsBoundaryQ) {
            face_index  = finterface_indexes_Gamma[k_face];
            duplet      = fleft_right_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_indexes_gamma[k_face];
            duplet      = fleft_right_indexes_gamma[k_face];
        }
        
        cel_face    = transport_mesh->Element(face_index);
        
        
        if (!cel_face) {
            DebugStop();
        }
        face_gel = cel_face->Reference();
        
        if (!face_gel) {
            DebugStop();
        }
        
        long left_geo_index     = duplet.first;
        long right_geo_index    = duplet.second;
        
        left_gel    = flux_mesh->Reference()->Element(left_geo_index);
        right_gel   = flux_mesh->Reference()->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        // Computing face connec index associated to the element face being integrated
        this->ComputeFaceIndex(left_gel,face_sides);
        
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace *> (left_cel);
        
        int nfaces = face_sides.size();
        
        int face_side   = -1;
        int i_face      = -1;
        
        if(!IdentifyFace(face_side,left_gel,face_gel)){
            std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
            DebugStop();
        }
        else{
            for (int i = 0; i < nfaces; i++) {
                if (face_sides[i] == face_side) {
                    i_face = i;
                    break;
                }
            }
            if (i_face == -1) {
                std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
                DebugStop();
            }
        }
        
        this->ElementDofFaceIndexes(i_face,intel_vol, dof_indexes);
        
        TPZIntPoints *face_int_points = left_gel->CreateSideIntegrationRule(face_sides[i_face], cel_face->GetgOrder());
        int npoints = face_int_points->NPoints();
        int nshapes = left_cel->Connect(i_face).NShape();
        
        blocks_dimensions[k_face].first = 1;//npoints;
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
void TRMBuildTransfers::Fill_un_To_Transport(TPZAutoPointer< TPZCompMesh> flux_mesh, TPZAutoPointer< TPZCompMesh> transport_mesh, bool IsBoundaryQ){
    

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
    
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<long,10> dof_indexes;
    
    TPZCompEl * cel_face;
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    
    TPZManVector<long> indices;
    std::pair<long, long> duplet;
    TPZManVector<int,10> face_sides;
    TPZFMatrix<REAL> normals;
    long face_index;
    long n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_indexes_gamma.size();
    }
    
    TPZFNMatrix<100,double> block;
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (IsBoundaryQ) {
            face_index  = finterface_indexes_Gamma[k_face];
            duplet      = fleft_right_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_indexes_gamma[k_face];
            duplet      = fleft_right_indexes_gamma[k_face];
        }
        
        cel_face    = transport_mesh->Element(face_index);
        
        if (!cel_face) {
            DebugStop();
        }

        face_gel = cel_face->Reference();
        
        if (!face_gel) {
            DebugStop();
        }
        
        long left_geo_index     = duplet.first;
        long right_geo_index    = duplet.second;
        
        left_gel    = flux_mesh->Reference()->Element(left_geo_index);
        right_gel   = flux_mesh->Reference()->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        // Computing face connect index associated to the element face being integrated
        // The method is really simple the inteface is detected by connect shared between hdiv volumetric elements
        this->ComputeFaceNormals(left_gel,face_sides,normals);
        
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace *> (left_cel);
        
        int nfaces = face_sides.size();
        
        int face_side   = -1;
        int i_face      = -1;
        
        if(!IdentifyFace(face_side,left_gel,face_gel)){
            std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
            DebugStop();
        }
        else{
            for (int i = 0; i < nfaces; i++) {
                if (face_sides[i] == face_side) {
                    i_face = i;
                    break;
                }
            }
            if (i_face == -1) {
                std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
                DebugStop();
            }
        }
        
        this->ElementDofFaceIndexes(i_face,intel_vol, dof_indexes);
        TPZIntPoints *face_int_points   = left_gel->CreateSideIntegrationRule(face_sides[i_face], cel_face->GetgOrder());
        
        int npoints = face_int_points->NPoints();
        int nshapes = left_cel->Connect(i_face).NShape();
  
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
        REAL ElementMeasure   = DimensionalMeasure(cel_face);

        REAL detjac;
        TPZFMatrix<REAL> jac,axes,jacinv;

        for (int ip = 0; ip < npoints; ip++) {
         
            // Get the vectorial phi
            face_int_points->Point(ip, par_duplet, w);
            face_gel->Jacobian(par_duplet, jac, axes, detjac, jacinv);
            intel_vol->SideShapeFunction(face_sides[i_face], par_duplet, phi, dphidxi);
            
            
            for (int jp = 0; jp < nshapes; jp++) {
                block(ip,jp) =  w * detjac * phi(jp,0)/ElementMeasure;
            }
            
        }
        
        for (int ip = 0; ip < npoints; ip++) {
        
            for (int jp = 0; jp < nshapes; jp++) {
                block_integral(0,jp) += block(ip,jp);
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
    if (this->SimulationData()->RawData()->fOmegaIds[0] != 1) {
        DebugStop();
    }
    
#endif
    
    TPZManVector<long,30>  point_index_trans;
    TPZManVector<long,30>  point_index_l;
    TPZManVector<long,30>  point_index_r;
    
    REAL p_avg_n_l = -1.0;
    REAL p_avg_n_r = -1.0;
    
    cmesh_transport->LoadReferences();
    TPZGeoMesh * geometry = cmesh_transport->Reference();
    long n_interfaces;
    int dimension = geometry->Dimension();
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_indexes_gamma.size();
    }
    
    int rock_id = this->SimulationData()->RawData()->fOmegaIds[0];
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * rock_material = cmesh_flux->FindMaterial(rock_id);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material);

    
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
            long pos = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                for(int iflux = 0; iflux < fun_dof_scatter_Gamma[iface].size(); iflux++) {
                    ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_Gamma[iface][iflux],0);
                    pos++;
                }
            }
            
            // Step two
            TPZFMatrix<STATE> un_at_intpoints;
            fun_To_Transport_Gamma.Multiply(ScatterFluxes,un_at_intpoints);
            
            // Step three
            // Trasnfering integrated normal fluxes values
            int counter = 0;
            int i = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                TPZGeoEl *gel    = geometry->Element(finterface_indexes_Gamma[iface]);
                TPZGeoEl *gel_l  = geometry->Element(fleft_right_indexes_Gamma[iface].first);
                TPZGeoEl *gel_r  = geometry->Element(fleft_right_indexes_Gamma[iface].second);

                
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
        
        TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin>  * material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin> *>(material);
        
        if (!material_mem) {
            DebugStop();
        }

        int np_cmesh = material_mem->GetMemory().NElements();

        // Step one
        TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_gamma.Cols(),1,0.0);
        long pos = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {
            for(int iflux = 0; iflux < fun_dof_scatter_gamma[iface].size(); iflux++) {
                ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_gamma[iface][iflux],0);
                pos++;
            }
        }

        // Step two
        TPZFMatrix<STATE> un_at_intpoints;
        fun_To_Transport_gamma.Multiply(ScatterFluxes,un_at_intpoints);

        // Step three
        // Trasnfering integrated normal fluxes values
        for(long i = 0; i < np_cmesh; i++){
            material_mem->GetMemory()[i].Set_un(un_at_intpoints(i,0));
        }
        
        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int i = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {

            TPZGeoEl *gel_l  = geometry->Element(fleft_right_indexes_gamma[iface].first);
            TPZGeoEl *gel_r  = geometry->Element(fleft_right_indexes_gamma[iface].second);
            
            
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
void TRMBuildTransfers::GlobalPointIndexes(TPZCompEl * cel, TPZManVector<long,30> &int_point_indexes){
    
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
void TRMBuildTransfers::GlobalPointIndexesInterface(TPZCompEl * int_cel, TPZManVector<long,30> &int_point_indexes){
    
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
    bool IsMembershipQ;
    for (int iside = 0; iside < volu_nsides; iside++) {
        IsMembershipQ = bool(vol->NeighbourExists(iside, face_itself));
        if (IsMembershipQ) {
            side = iside;
            break;
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
            sides[2] = 4;
            
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
void TRMBuildTransfers::ComputeLeftRight(TPZAutoPointer< TPZCompMesh> transport_mesh){
    
#ifdef PZDEBUG
    if (!transport_mesh) {
        std::cout << "There is no computational transport mesh, transport_mesh = Null." << std::endl;
        DebugStop();
    }
#endif
    
    long nel = transport_mesh->NElements();
    long face_index;
    std::pair <long,long> duplet;
    int dimension = transport_mesh->Reference()->Dimension();
    for (long icel = 0; icel < nel; icel++) {
        
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
        
        face_index  = interface->Index();
        duplet      = std::make_pair(left_cel->Reference()->Index(), right_cel->Reference()->Index());
        
        if(left_cel->Dimension() != dimension ||  right_cel->Dimension() != dimension){
            
            fleft_right_indexes_Gamma.Push(duplet);
            finterface_indexes_Gamma.Push(face_index);
            continue;
        }
        
        fleft_right_indexes_gamma.Push(duplet);
        finterface_indexes_gamma.Push(face_index);
        
    }
    
#ifdef PZDEBUG
    if (finterface_indexes_Gamma.size() == 0) {
        DebugStop();
    }
    if (finterface_indexes_gamma.size() == 0) {
        std::cout << "Warning:: No inner interfaces were found" << std::endl;
    }
#endif
    
}

/** @brief Dimensionla Measure of the elemnt */
REAL TRMBuildTransfers::DimensionalMeasure(TPZCompEl * cel){
    
#ifdef PZDEBUG
    if (!cel) {
        DebugStop();
    }
#endif
    REAL measure = 0.0;
    int order = 5;
    TPZGeoEl * gel = cel->Reference();
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


void TRMBuildTransfers::ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
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

void TRMBuildTransfers::ElementDofFaceIndexes(int connect_index,TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    TPZConnect  & con = intel->Connect(connect_index);
    long seqnumber = con.SequenceNumber();
    long position = intel->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}


/** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
void TRMBuildTransfers::FillComputationalElPairs(TPZAutoPointer< TPZCompMesh>  cmesh_mf_mixed, TPZAutoPointer< TPZCompMesh>  cmesh_mf_transport){

    fmixed_transport_indexes.Resize(0);
    
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
    
    std::pair<long, std::pair <long,long> > gel_indexes;
    
    for (long i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() != dimension) {
            continue;
        }
        gel_indexes.first = gel->Index();
        gel_indexes.second.first = -1;
        gel_indexes.second.second = -1;
        fmixed_transport_indexes.Push(gel_indexes);

    }
    
    // counting volumetric elements
    long nvol_elements = fmixed_transport_indexes.size();
    fmixed_transport_indexes.Resize(nvol_elements);
    
    // inserting mixed indexes
    cmesh_mf_mixed->LoadReferences();
    for(long ivol = 0; ivol < nvol_elements; ivol++){

        TPZCompEl * mixed_cel = geometry->Element(fmixed_transport_indexes[ivol].first)->Reference();
        
#ifdef PZDEBUG
        if (!mixed_cel) {
            DebugStop();
        }
#endif
        fmixed_transport_indexes[ivol].second.first = mixed_cel->Index();
        
    }
    
    if(fSimulationData->IsOnePhaseQ()){
        return;
    }
    
    // inserting transport indexes
    cmesh_mf_transport->LoadReferences();
    for(long ivol = 0; ivol < nvol_elements; ivol++){
        
        TPZCompEl * trans_cel = geometry->Element(fmixed_transport_indexes[ivol].first)->Reference();
        
#ifdef PZDEBUG
        if (!trans_cel) {
            DebugStop();
        }
#endif
        fmixed_transport_indexes[ivol].second.second = trans_cel->Index();
        
    }
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Old Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//
//void TRMBuildTransfers::TransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_pressure, TPZAutoPointer< TPZCompMesh> cmesh_multiphysics){
//
//#ifdef PZDEBUG
//    if (!cmesh_multiphysics) {
//        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
//        DebugStop();
//    }
//#endif
//
//    // Getting the integration points
//    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
//    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
//    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
//    int np_cmesh_des = material_memory.NElements();
//    
//    TPZFMatrix<STATE> pressure_at_intpoints;
//    fTransferPressure_To_Mixed.Multiply(cmesh_pressure->Solution(),pressure_at_intpoints);
//    
////    TRMMemory data;
//    for(long i = 0; i <  np_cmesh_des; i++){
//        STATE pressure = pressure_at_intpoints(i,0);
//        material_memory[i].SetPressure(pressure);
////        material_memory[i] = data;
//    }
//    
//    associated_material->SetMemory(material_memory);
//    
//}
//
//void TRMBuildTransfers::TransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer< TPZCompMesh> cmesh_multiphysics){
//    
//#ifdef PZDEBUG
//    if (!cmesh_multiphysics) {
//        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
//        DebugStop();
//    }
//#endif
//    
//    // Getting the integration points
//    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
//    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
//    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
//    int n_int_points = material_memory.NElements();
//    
//    TPZFMatrix<STATE> x_flux_at_intpoints, y_flux_at_intpoints, z_flux_at_intpoints, divflux_at_points;
//    fTransfer_X_Flux_To_Mixed.Multiply(cmesh_flux->Solution(),x_flux_at_intpoints);
//    fTransfer_Y_Flux_To_Mixed.Multiply(cmesh_flux->Solution(),y_flux_at_intpoints);
//    fTransfer_Z_Flux_To_Mixed.Multiply(cmesh_flux->Solution(),z_flux_at_intpoints);
//    fTransferDivergenceTo_Mixed.Multiply(cmesh_flux->Solution(), divflux_at_points);
//    
////    TRMMemory data;
//    TPZManVector<STATE> flux(3);
//    REAL divu;
//    for(long i = 0; i <  n_int_points; i++)
//    {
//        flux[0] = x_flux_at_intpoints(i,0);
//        flux[1] = y_flux_at_intpoints(i,0);
//        flux[2] = z_flux_at_intpoints(i,0);
//        divu = divflux_at_points(i,0);
//        material_memory[i].SetTotal_Flux(flux);
//        material_memory[i].SetDiv_Flux(divu);
////        material_memory[i] = data;
//    }
//    
//    associated_material->SetMemory(material_memory);
//    
//}

//void TRMBuildTransfers::ComputeTransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
//
//#ifdef PZDEBUG
//    if (!cmesh_multiphysics) {
//        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
//        DebugStop();
//    }
//#endif
//
//    long nel = cmesh_multiphysics->NElements();
//
//    // Getting the total integration point of the destination cmesh
//    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
//    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
//    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
//    int n_int_points = material_memory.NElements();
//
//
//    // Resizing IA JA and A
//    TPZAutoPointer<TPZVec<long> > IA = new TPZVec<long>;
//    TPZAutoPointer<TPZVec<long> > JA = new TPZVec<long>;
//    TPZAutoPointer<TPZVec<double> > Ax = new TPZVec<double>;
//    TPZAutoPointer<TPZVec<double> > Ay = new TPZVec<double>;
//    TPZAutoPointer<TPZVec<double> > Az = new TPZVec<double>;
//    TPZAutoPointer<TPZVec<double> > Ad = new TPZVec<double>;
//
//    this->CreateTransferFlux_To_Mixed_V(cmesh_multiphysics, mesh_index, IA, JA, Ax, Ay, Az, Ad);
//
//    TPZAutoPointer<TPZVec<long> > IAg = new TPZVec<long>;
//    TPZAutoPointer<TPZVec<long> > JAg = new TPZVec<long>;
//    TPZAutoPointer<TPZVec<double> >Aw = new TPZVec<double>;
//    TPZAutoPointer<TPZVec<double> >Adet = new TPZVec<double>;
//    TPZAutoPointer<TPZVec<double> >Arhs = new TPZVec<double>;
//
//    this->CreateTransferGeometricData_To_Mixed_V(cmesh_multiphysics, mesh_index, IAg, JAg, Aw, Adet, Arhs);
//
//    int origin = mesh_index;
//    int destination = mesh_index;
//    int volumetric_elements = 0;
//    int iblock = 0;
//    TPZMaterialData data;
//
//    for (long icel = 0; icel < nel; icel++) {
//
//        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
//        if (!cel) {
//            DebugStop();
//        }
//
//        //        TPZCondensedCompEl * mf_cel_condensed = dynamic_cast<TPZCondensedCompEl *> (cel);
//        //        if(!mf_cel_condensed){
//        //            DebugStop();
//        //        }
//        //        TPZCompEl * mf_cel_cond_ref = mf_cel_condensed->ReferenceCompEl();
//
//        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
//        if(!mf_cel)
//        {
//            DebugStop();
//        }
//
//        // Avoiding al the materials that don't correspond to the volume
//        if(mf_cel->Material()->Id() != _ReservMatId)
//        {
//            continue;
//        }
//
//        TPZInterpolationSpace * intel_o = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(origin));
//        TPZInterpolationSpace * intel_d = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(destination));
//        if (!intel_o || !intel_d) {
//            DebugStop();
//        }
//
//        TPZCompMesh * cmesh_o = intel_o->Mesh();
//
//        // Computing the global integration points indexes
//        TPZManVector<long> globindexes;
//        mf_cel->GetMemoryIndices(globindexes);
//
//        // Computing the local integration points indexes
//        const TPZIntPoints & int_points_cel_d = mf_cel->GetIntegrationRule();
//        int np_cel_d = int_points_cel_d.NPoints();
//
//        if (globindexes.size() != np_cel_d) {
//            DebugStop();
//        }
//
//        volumetric_elements++;
//
//        // Computing over all integration points of the compuational element cel
//        TPZFNMatrix<100,REAL> phi(intel_o->NShapeF(),1,0.0);
//        int el_dim = mf_cel->Reference()->Dimension();
//        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel_o->NShapeF(),0.0);
//
//        int nconnect = intel_o->NConnects();
//        int iphicount = 0;
//        // Compute all the phi values and push inside corresponding j-destination;
//        int el_nshape = 0;
//        for (int icon = 0; icon < nconnect; icon++) {
//
//            TPZConnect  & con = intel_o->Connect(icon);
//            long seqnumber = con.SequenceNumber();
//            long position = cmesh_o->Block().Position(seqnumber);
//            el_nshape += con.NShape();
//
//        }
//
//        // Creating the block with the right size
//        TPZFMatrix<double> blockx(np_cel_d,el_nshape);
//        TPZFMatrix<double> blocky(np_cel_d,el_nshape);
//        TPZFMatrix<double> blockz(np_cel_d,el_nshape);
//        TPZFMatrix<double> blockd(np_cel_d,el_nshape);
//
//        for (int ip = 0; ip < np_cel_d ; ip++)
//        {
//
//            TPZManVector<REAL,3> qsi(3,0.0);
//            STATE w;
//            int_points_cel_d.Point(ip, qsi, w);
//
//            // Indetifying the right global index
//            long globindex = globindexes[ip];
//
//            long JAcount = IA->operator[](globindex);
//            long JgAcount = IAg->operator[](globindex);
//
//            // Get the vectorial phi
//            intel_o->Shape(qsi, phi, dphidxi);
//            intel_o->InitMaterialData(data);
//            intel_o->ComputeRequiredData(data,qsi);
//
//            // Hdiv space
//            REAL JacobianDet  = data.detjac;
//            TPZFNMatrix<100,REAL> phiuH1         = data.phi;   // For H1  test functions Q
//            TPZFNMatrix<300,STATE> dphiuH1       = data.dphi; // Derivative For H1  test functions
//            TPZFNMatrix<300,STATE> dphiuH1axes   = data.dphix; // Derivative For H1  test functions
//
//            TPZFNMatrix<300,STATE> Qaxes = data.axes;
//            TPZFNMatrix<300,STATE> QaxesT;
//            TPZFNMatrix<300,STATE> Jacobian = data.jacobian;
//            TPZFNMatrix<300,STATE> JacobianInverse = data.jacinv;
//
////            TPZFMatrix<STATE> GradOfX;
//            TPZFNMatrix<300,STATE> GradOfXInverse;
//            TPZFNMatrix<3,STATE> VectorOnMaster;
//            TPZFNMatrix<3,STATE> VectorOnXYZ(3,1,0.0);
//
//            Qaxes.Transpose(&QaxesT);
////                                QaxesT.Multiply(Jacobian, GradOfX);
//            JacobianInverse.Multiply(Qaxes, GradOfXInverse);
//
////            TPZManVector<STATE,1> fval(1,0.0);
////            ExactLaplacian(data.x,fval); // defined as static member
////            REAL rhs = fval[0];
////
//////            // Feed material memory
//////            material_memory[globindex].SetWeight(w);
//////            material_memory[globindex].SetDetJac(JacobianDet);
//////            material_memory[globindex].SetX(data.x);
//////            material_memory[globindex].SetRhs(rhs);
////////            std::cout << "globindex = " <<  globindex << std::endl;
//////            // Feed vectors of geometric data for each volumetric element index
//////            JAg->operator[](JgAcount) = volumetric_elements-1;
//////            Aw->operator[](JgAcount) = w;
//////            Adet->operator[](JgAcount) = JacobianDet;
//////            Arhs->operator[](JgAcount) = rhs;
//////            JgAcount++;
//
//            int nconnect = intel_o->NConnects();
//            int iphicount = 0;
//            // Compute all the phi values and push inside corresponding j-destination;
//
//            for (int icon = 0; icon < nconnect; icon++) {
//
//                TPZConnect  & con = intel_o->Connect(icon);
//                long seqnumber = con.SequenceNumber();
////                long position = cmesh_o->Block().Position(seqnumber);
//                int nshape = con.NShape();
//
//                for (int ish=0; ish < nshape; ish++) {
//
//                    int vector_index = data.fVecShapeIndex[iphicount].first;
//                    int shape_index = data.fVecShapeIndex[iphicount].second;
//
//                    VectorOnXYZ(0,0) = data.fNormalVec(0,vector_index);
//                    VectorOnXYZ(1,0) = data.fNormalVec(1,vector_index);
//                    VectorOnXYZ(2,0) = data.fNormalVec(2,vector_index);
//
//                    GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
//                    VectorOnMaster *= JacobianDet;
//
////                    std::cout << "VectorOnMaster = " << VectorOnMaster << std::endl;
//
//                    /* Contravariant Piola mapping preserves the divergence */
//                    int dim =VectorOnMaster.Rows();
//                    TPZFNMatrix<9,REAL> GradOfPhi(dim,dim);
//                    for (int ir = 0; ir < dim; ir++) {
//
//                        //  Compute grad_{hat}(PhiH1)
//
//                        GradOfPhi(ir,0) = VectorOnMaster(ir,0)*dphiuH1(0,shape_index);
//                        GradOfPhi(ir,1) = VectorOnMaster(ir,0)*dphiuH1(1,shape_index);
//                        GradOfPhi(ir,2) = VectorOnMaster(ir,0)*dphiuH1(2,shape_index);
//
//                    }
//                    GradOfPhi *= (1.0/JacobianDet);
//                    REAL divPhi_Hdiv = ( GradOfPhi(0,0) + GradOfPhi(1,1) + GradOfPhi(2,2));
//
//                    TPZManVector<STATE> phi_Hdiv(3);
//                    phi_Hdiv[0] = phi(shape_index,0)*VectorOnXYZ(0,0);
//                    phi_Hdiv[1] = phi(shape_index,0)*VectorOnXYZ(1,0);
//                    phi_Hdiv[2] = phi(shape_index,0)*VectorOnXYZ(2,0);
//
////                    JA->operator[](JAcount) = position + ish;
////                    Ax->operator[](JAcount) = phi_Hdiv[0];
////                    Ay->operator[](JAcount) = phi_Hdiv[1];
////                    Az->operator[](JAcount) = phi_Hdiv[2];
////                    Ad->operator[](JAcount) = divPhi_Hdiv;
//                    blockx(ip,iphicount) = phi_Hdiv[0];
//                    blocky(ip,iphicount) = phi_Hdiv[1];
//                    blockz(ip,iphicount) = phi_Hdiv[2];
//                    blockd(ip,iphicount) = divPhi_Hdiv;
//                    iphicount++;
////                    JAcount++;
//                }
//
//            }
//        }
//
//        fTransfer_X_Flux_To_Mixed_V.SetBlock(iblock, blockx);
//        fTransfer_Y_Flux_To_Mixed_V.SetBlock(iblock, blocky);
//        fTransfer_Z_Flux_To_Mixed_V.SetBlock(iblock, blockz);
//        fTransferDivergenceTo_Mixed_V.SetBlock(iblock, blockd);
//        iblock++;
//    }
//
//
////    fTransfer_X_Flux_To_Mixed_V.
//
////    fTransfer_X_Flux_To_Mixed_V.SetData(IA, JA, Ax);
////    fTransfer_Y_Flux_To_Mixed_V.SetData(IA, JA, Ay);
////    fTransfer_Z_Flux_To_Mixed_V.SetData(IA, JA, Az);
////    fTransferDivergenceTo_Mixed_V.SetData(IA, JA, Ad);
////
////    fIntegrationWeights_V.SetData(IAg,JAg,Aw);
////    fJacobianDet_V.SetData(IAg,JAg,Adet);
////    fRhs_V.SetData(IAg,JAg,Arhs);
//    associated_material->SetMemory(material_memory);
//}
//
///** @brief exact laplacian */
//void TRMBuildTransfers::ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &f)
//{
//    REAL x,y,z;
//    x = pt[0];
//    y = pt[1];
//    z = pt[2];
//
//    f[0] = 6;//2.*(-1. + x)*x*(-1. + y)*y + 2.*(-1. + x)*x*(-1. + z)*z + 2.*(-1. + y)*y*(-1. + z)*z;
//}
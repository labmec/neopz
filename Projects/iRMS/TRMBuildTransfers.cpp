//
//  TRMBuildTransfers.cpp
//  PZ
//
//  Created by Omar on 10/27/15.
//
//


#include "TRMBuildTransfers.h"

#include <stdio.h>
#include "pzgmesh.h"
#include "pzcmesh.h"

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

void TRMBuildTransfers::Initialize_u_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 3; // vectorial
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

void TRMBuildTransfers::Fill_u_To_Mixed(TPZAutoPointer< TPZCompMesh > cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_u_To_Mixed(cmesh_multiphysics, mesh_index);
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 3; // vector
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
        
        if(intel->Dimension() < n_var_dim){ // two dimensional elements
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

void TRMBuildTransfers::Transfer_u_To_Mixed_Memory(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_multiphysics){
    

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
    TPZManVector<STATE> u(dim,0.0);
    for(long i = 0; i <  np_cmesh; i++){
        for (int id = 0; id < dim ; id++) {
            u[id]= Flux_at_intpoints(i*dim+id,0);
        }
        associated_material->GetMemory()[i].SetTotal_Flux(u);
    }
    
}

void TRMBuildTransfers::Initialize_p_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
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
            // there is no boundary elements for pressure
            block_dim.first = 0;
            block_dim.second = 0;
            fp_dof_scatter[element_index] = dof_indexes;
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


void TRMBuildTransfers::Transfer_p_To_Mixed_Memory(TPZCompMesh * cmesh_pressure, TPZCompMesh * cmesh_multiphysics){

    
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
    TPZFMatrix<STATE> Pressure_at_intpoints;
    fp_To_Mixed.Multiply(ScatterPressure,Pressure_at_intpoints);
    // Trasnfering values
    for(long i = 0; i <  np_cmesh; i++){
        associated_material->GetMemory()[i].SetPressure(Pressure_at_intpoints(i,0));
    }
    
}

void TRMBuildTransfers::ElementDofIndexes(TPZInterpolationSpace * intel, TPZVec<long> &dof_indexes){

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

/** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Initialize_un_To_Transport_a(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    //* seeking for total blocks */
    
    
    long nel = cmesh_multiphysics->NElements();
    int n_var_dim = 3; // vectorial
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fun_dof_scatter.Resize(nel);
    
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

void TRMBuildTransfers::ElementDofFaceIndexes(TPZInterpolationSpace * intel, TPZVec<long> &dof_indexes){
    
    
    DebugStop(); // method not implemented!
    
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

/** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void Initialize_un_To_Transport_a(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
}

/** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void Fill_un_To_Transport_a(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_un_To_Transport_a(cmesh_multiphysics, mesh_index);
    
//    long nel = cmesh_multiphysics->NElements();
//    int n_var_dim = 3; // vector
//    long element_index = 0;
//    
//    TPZMaterialData data;
//    
//    std::pair<long, long> block_dim;
//    
//    for (long icel = 0; icel < nel; icel++) {
//        
//        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
//        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
//        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
//        element_index = mf_cel->Index();
//        
//        // Getting local integration index
//        TPZManVector<long> int_point_indexes(0,0);
//        TPZManVector<long> dof_indexes(0,0);
//        
//        mf_cel->GetMemoryIndices(int_point_indexes);
//        dof_indexes = fu_dof_scatter[element_index];
//        
//        block_dim.first = int_point_indexes.size();
//        block_dim.second = dof_indexes.size();
//        
//        
//        // Computing the local integration points indexes
//        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
//        int np_cel = int_points_mixed.NPoints();
//        
//#ifdef PZDEBUG
//        if (int_point_indexes.size() != np_cel) {
//            DebugStop();
//        }
//#endif
//        
//        // Computing over all integration points of the compuational element cel
//        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
//        int el_dim = mf_cel->Reference()->Dimension();
//        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
//        TPZFMatrix<double> block;
//        
//        if(intel->Dimension() < n_var_dim){ // two dimensional elements
//            block.Resize(block_dim.first,block_dim.second);        }
//        else{
//            block.Resize(block_dim.first*n_var_dim,block_dim.second);
//        }
//        
//        for (int ip = 0; ip < block_dim.first ; ip++)
//        {
//            TPZManVector<REAL,3> qsi(el_dim,0.0);
//            STATE w;
//            int_points_mixed.Point(ip, qsi, w);
//            // Get the vectorial phi
//            intel->Shape(qsi, phi, dphidxi);
//            intel->InitMaterialData(data);
//            intel->ComputeRequiredData(data,qsi);
//            
//            for (int id = 0; id < n_var_dim; id++) {
//                for (int jp = 0; jp < block_dim.second; jp++) {
//                    int vector_index = data.fVecShapeIndex[jp].first;
//                    int shape_index = data.fVecShapeIndex[jp].second;
//                    block(ip*n_var_dim+id,jp) = phi(shape_index,0)*data.fNormalVec(id,vector_index);
//                }
//            }
//            
//        }
//        
//        fu_To_Mixed.SetBlock(element_index, block);
//        
//    }
    
    return;
    
    
}


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
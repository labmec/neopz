//
//  TRMBuildTransfers.cpp
//  PZ
//
//  Created by Omar on 10/27/15.
//
//


#include "TRMBuildTransfers.h"

#include "tpzintpoints.h"
#include "pzmatwithmem.h"
#include "TRMMemory.h"
#include "TRMMixedDarcy.h"
#include "pzmaterial.h"
#include "TRMFlowConstants.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"
#include "pzcondensedcompel.h"

/** @brief Default constructor */
TRMBuildTransfers::TRMBuildTransfers(){

}

/** @brief Default desconstructor */
TRMBuildTransfers::~TRMBuildTransfers(){

}

void TRMBuildTransfers::CreateTransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZManVector<long> &IA, TPZManVector<long> &JA, TPZManVector<STATE> &Ax, TPZManVector<STATE> &Ay, TPZManVector<STATE> &Az, TPZManVector<STATE> Ad){
    
    long nel = cmesh_multiphysics->NElements();
    
    int dimension = cmesh_multiphysics->Dimension();
    
    // Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
    int n_int_points = material_memory.NElements();
    
    IA.Resize((n_int_points+1), 0);
    
    int origin      = mesh_index;
    int destination = mesh_index;
    TPZCompMesh * cmesh_o;
    TPZCompMesh * cmesh_d;
    int mesh_o_nequ = 0;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        if (!cel) {
            DebugStop();
        }
        
        //         TPZCondensedCompEl * mf_cel_condensed = dynamic_cast<TPZCondensedCompEl *> (cel);
        //         if(!mf_cel_condensed){
        //             DebugStop();
        //         }
        //        TPZCompEl * mf_cel_cond_ref = mf_cel_condensed->ReferenceCompEl();
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        if(!mf_cel)
        {
            DebugStop();
        }
        
        // Avoiding al the materials that don't correspond to the volume
        if(mf_cel->Material()->Id() != _ReservMatId)
        {
            continue;
        }
        
        TPZInterpolationSpace * intel_o = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(origin));
        TPZInterpolationSpace * intel_d = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(destination));
        if (!intel_o || !intel_d) {
            DebugStop();
        }
        
        cmesh_o = intel_o->Mesh();
        cmesh_d = intel_d->Mesh();
        if (!cmesh_o || !cmesh_d) {
            DebugStop();
        }
        
        mesh_o_nequ = cmesh_o->NEquations();
        
        // Computing the global integration points indexes
        TPZManVector<long> globindexes;
        mf_cel->GetMemoryIndices(globindexes);
        int nshapes = intel_d->NShapeF();
        
        for(long i=0; i< globindexes.size(); i++)
        {
            long glob = globindexes[i];
            IA[(glob+1)] = nshapes;
        }
        
    }
    
    for(long i=0; i< n_int_points; i++ )
    {
        IA[(i+1)] += IA[i];
    }

    JA.Resize(IA[n_int_points]);
    Ax.Resize(IA[n_int_points]);
    Ay.Resize(IA[n_int_points]);
    Az.Resize(IA[n_int_points]);
    Ad.Resize(IA[n_int_points]);
    fTransfer_X_Flux_To_Mixed.Resize(n_int_points, mesh_o_nequ);
    fTransfer_Y_Flux_To_Mixed.Resize(n_int_points, mesh_o_nequ);
    fTransfer_Z_Flux_To_Mixed.Resize(n_int_points, mesh_o_nequ);
    fTransferDivergenceTo_Mixed.Resize(n_int_points, mesh_o_nequ);
    
}

void TRMBuildTransfers::CreateTransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZManVector<long> &IA, TPZManVector<long> &JA, TPZManVector<STATE> &A){
    
    long nel = cmesh_multiphysics->NElements();
    
    // Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
    int n_int_points = material_memory.NElements();
    
    IA.Resize(n_int_points+1, 0);
    
    int origin      = mesh_index;
    int destination = mesh_index;
    TPZCompMesh * cmesh_o;
    TPZCompMesh * cmesh_d;
    int mesh_o_nequ = 0;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        if (!cel) {
            DebugStop();
        }
        
        //         TPZCondensedCompEl * mf_cel_condensed = dynamic_cast<TPZCondensedCompEl *> (cel);
        //         if(!mf_cel_condensed){
        //             DebugStop();
        //         }
        //        TPZCompEl * mf_cel_cond_ref = mf_cel_condensed->ReferenceCompEl();
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        if(!mf_cel)
        {
            DebugStop();
        }
        
        // Avoiding all materials that don't correspond to volumetric ones
        if(mf_cel->Material()->Id() != _ReservMatId)
        {
            continue;
        }
        
        TPZInterpolationSpace * intel_o = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(origin));
        TPZInterpolationSpace * intel_d = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(destination));
        if (!intel_o || !intel_d) {
            DebugStop();
        }
        
        cmesh_o = intel_o->Mesh();
        cmesh_d = intel_d->Mesh();
        if (!cmesh_o || !cmesh_d) {
            DebugStop();
        }
        
        mesh_o_nequ = cmesh_o->NEquations();
        
        // Computing the global integration points indexes
        TPZManVector<long> globindexes;
        mf_cel->GetMemoryIndices(globindexes);
        int nshapes = intel_d->NShapeF();
        long nglob_indexes = globindexes.size();
        for(long i=0; i< nglob_indexes; i++ )
        {
            long glob = globindexes[i];
            IA[glob+1] = nshapes;
        }
        
    }
    
    for(long i=0; i< n_int_points; i++ )
    {
        IA[i+1] += IA[i];
    }
    
    JA.Resize(IA[n_int_points]);
    A.Resize(IA[n_int_points]);
    
    fTransferPressure_To_Mixed.Resize(n_int_points, mesh_o_nequ);
    
}

void TRMBuildTransfers::ComputeTransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    long nel = cmesh_multiphysics->NElements();
    
    // Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
    int n_int_points = material_memory.NElements();
    
    // Resizing IA JA and A
    TPZManVector<long> IA;
    TPZManVector<long> JA;
    TPZManVector<double> Ax,Ay,Az,Ad;
    
    TPZFMatrix<STATE> VectorOnXYZ(3,1,0.0);
    
    this->CreateTransferFlux_To_Mixed(cmesh_multiphysics, mesh_index, IA, JA, Ax,Ay,Az,Ad);
    
    int origin = mesh_index;
    int destination = mesh_index;
    TPZMaterialData data;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        if (!cel) {
            DebugStop();
        }
        
        //        TPZCondensedCompEl * mf_cel_condensed = dynamic_cast<TPZCondensedCompEl *> (cel);
        //        if(!mf_cel_condensed){
        //            DebugStop();
        //        }
        //        TPZCompEl * mf_cel_cond_ref = mf_cel_condensed->ReferenceCompEl();
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        if(!mf_cel)
        {
            DebugStop();
        }
        
        // Avoiding al the materials that don't correspond to the volume
        if(mf_cel->Material()->Id() != _ReservMatId)
        {
            continue;
        }
        
        TPZInterpolationSpace * intel_o = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(origin));
        TPZInterpolationSpace * intel_d = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(destination));
        if (!intel_o || !intel_d) {
            DebugStop();
        }
        
        TPZCompMesh * cmesh_o = intel_o->Mesh();
        
        // Computing the global integration points indexes
        TPZManVector<long> globindexes;
        mf_cel->GetMemoryIndices(globindexes);
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_cel_d = mf_cel->GetIntegrationRule();
        int np_cel_d = int_points_cel_d.NPoints();
        
        if (globindexes.size() != np_cel_d) {
            DebugStop();
        }
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel_o->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel_o->NShapeF(),0.0);
        for (int ip = 0; ip < np_cel_d ; ip++)
        {
            
            // Here for each integration point
            TPZManVector<REAL,3> qsi(3,0.0);
            STATE w;
            int_points_cel_d.Point(ip, qsi, w);
            
            // Indetifying the right global index
            long globindex = globindexes[ip];
            long JAcount = IA[globindex];
            
            // 
            
            // Get the vectorial phi
            intel_o->Shape(qsi, phi, dphidxi);
            

            intel_o->InitMaterialData(data);
            intel_o->ComputeRequiredData(data,qsi);
        
            int nconnect = intel_o->NConnects();
            int iphicount = 0;
            // Compute all the phi values and push inside corresponding j-destination;
            for (int icon = 0; icon < nconnect; icon++) {
                
                TPZConnect  & con = intel_o->Connect(icon);
                long seqnumber = con.SequenceNumber();
                long position = cmesh_o->Block().Position(seqnumber);
                int nshape = con.NShape();
                
                for (int ish=0; ish < nshape; ish++) {
                    
                    int vector_index = data.fVecShapeIndex[iphicount].first;
                    int shape_index = data.fVecShapeIndex[iphicount].second;
                    
                    VectorOnXYZ(0,0) = data.fNormalVec(0,vector_index);
                    VectorOnXYZ(1,0) = data.fNormalVec(0,vector_index);
                    VectorOnXYZ(2,0) = data.fNormalVec(0,vector_index);
                    
                    TPZManVector<STATE> phi_Hdiv(3);
                    phi_Hdiv[0] = phi(shape_index,0)*VectorOnXYZ(0,0);
                    phi_Hdiv[1] = phi(shape_index,0)*VectorOnXYZ(1,0);
                    phi_Hdiv[2] = phi(shape_index,0)*VectorOnXYZ(2,0);
                    
                    JA[JAcount] = position + ish;
                    Ax[JAcount] = phi_Hdiv[0];
                    Ay[JAcount] = phi_Hdiv[1];
                    Az[JAcount] = phi_Hdiv[2];
                    iphicount++;
                    JAcount++;
                }
                
            }
        }
        
    }
    
    fTransfer_X_Flux_To_Mixed.SetData(IA, JA, Ax);
    fTransfer_Y_Flux_To_Mixed.SetData(IA, JA, Ay);
    fTransfer_Z_Flux_To_Mixed.SetData(IA, JA, Az);
    fTransferDivergenceTo_Mixed.SetData(IA, JA, Ad);
    
}


void TRMBuildTransfers::ComputeTransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index){

#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    long nel = cmesh_multiphysics->NElements();
    
    // Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
    int n_int_points = material_memory.NElements();
    
    // Resizing IA JA and A
    TPZManVector<long> IA;
    TPZManVector<long> JA;
    TPZManVector<double> A;
    
    this->CreateTransferPressure_To_Mixed(cmesh_multiphysics, mesh_index, IA, JA, A);
    
    int origin = mesh_index;
    int destination = mesh_index;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        if (!cel) {
            DebugStop();
        }
        
//        TPZCondensedCompEl * mf_cel_condensed = dynamic_cast<TPZCondensedCompEl *> (cel);
//        if(!mf_cel_condensed){
//            DebugStop();
//        }
//        TPZCompEl * mf_cel_cond_ref = mf_cel_condensed->ReferenceCompEl();
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        if(!mf_cel)
        {
            DebugStop();
        }
        
        // Avoiding al the materials that don't correspond to the volume
        if(mf_cel->Material()->Id() != _ReservMatId)
        {
            continue;
        }
        
        TPZInterpolationSpace * intel_o = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(origin));
        TPZInterpolationSpace * intel_d = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(destination));
        if (!intel_o || !intel_d) {
            DebugStop();
        }
        TPZCompMesh * cmesh_o = intel_o->Mesh();
        
        // Computing the global integration points indexes
        TPZManVector<long> globindexes;
        mf_cel->GetMemoryIndices(globindexes);
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_cel_d = mf_cel->GetIntegrationRule();
        int np_cel_d = int_points_cel_d.NPoints();
        
        if (globindexes.size() != np_cel_d) {
            DebugStop();
        }
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel_o->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel_o->NShapeF(),0.0);
        for (int ip = 0; ip < np_cel_d ; ip++)
        {
            TPZManVector<REAL,3> qsi(3,0.0);
            STATE w;
            int_points_cel_d.Point(ip, qsi, w);
            
            // Indetifying the right global index
           long globindex = globindexes[ip];
           long JAcount = IA[globindex];

            intel_o->Shape(qsi, phi, dphidxi);

            int nconnect = intel_o->NConnects();
            int iphicount = 0;
            // Compute all the phi values and push inside corresponding j-destination;
            for (int icon = 0; icon < nconnect; icon++) {
                TPZConnect  & con = intel_o->Connect(icon);
                long seqnumber = con.SequenceNumber();
                long position = cmesh_o->Block().Position(seqnumber);
                int nshape = con.NShape();

                for (int ish=0; ish < nshape; ish++) {
                    JA[JAcount] = position + ish;
                    A[JAcount] = phi(iphicount,0);
                    iphicount++;
                    JAcount++;
                }
            }
        }
        
    }
    
    fTransferPressure_To_Mixed.SetData(IA, JA, A);
    
}


void TRMBuildTransfers::TransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_pressure, TPZAutoPointer< TPZCompMesh> cmesh_multiphysics){

#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif

    // Getting the integration points
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
    int np_cmesh_des = material_memory.NElements();
    
    TPZFMatrix<STATE> pressure_at_intpoints;
    fTransferPressure_To_Mixed.Multiply(cmesh_pressure->Solution(),pressure_at_intpoints);
    
//    TRMMemory data;
    for(long i = 0; i <  np_cmesh_des; i++){
        STATE pressure = pressure_at_intpoints(i,0);
        material_memory[i].SetPressure(pressure);
//        material_memory[i] = data;
    }
    
    associated_material->SetMemory(material_memory);
    
}

void TRMBuildTransfers::TransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer< TPZCompMesh> cmesh_multiphysics){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    // Getting the integration points
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
    int n_int_points = material_memory.NElements();
    
    TPZFMatrix<STATE> x_flux_at_intpoints, y_flux_at_intpoints, z_flux_at_intpoints;
    fTransfer_X_Flux_To_Mixed.Multiply(cmesh_flux->Solution(),x_flux_at_intpoints);
    fTransfer_Y_Flux_To_Mixed.Multiply(cmesh_flux->Solution(),y_flux_at_intpoints);
    fTransfer_Z_Flux_To_Mixed.Multiply(cmesh_flux->Solution(),z_flux_at_intpoints);
    
//    TRMMemory data;
    TPZManVector<STATE> flux(3);
    for(long i = 0; i <  n_int_points; i++)
    {
        flux[0] = x_flux_at_intpoints(i,0);
        flux[1] = y_flux_at_intpoints(i,0);
        flux[2] = z_flux_at_intpoints(i,0);
        material_memory[i].SetTotal_Flux(flux);
//        material_memory[i] = data;
    }
    
    associated_material->SetMemory(material_memory);
    
}
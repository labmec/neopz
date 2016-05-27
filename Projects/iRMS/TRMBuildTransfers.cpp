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
    
    fSimulationData == NULL;

}

/** @brief Default desconstructor */
TRMBuildTransfers::~TRMBuildTransfers(){

}

void TRMBuildTransfers::CreateTransferFlux_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &Ax, TPZVec<STATE> &Ay, TPZVec<STATE> &Az, TPZVec<STATE> &Ad){
    
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
    fTransfer_X_Flux_To_Mixed_V.Resize(n_int_points, mesh_o_nequ);
    fTransfer_Y_Flux_To_Mixed_V.Resize(n_int_points, mesh_o_nequ);
    fTransfer_Z_Flux_To_Mixed_V.Resize(n_int_points, mesh_o_nequ);
    fTransferDivergenceTo_Mixed_V.Resize(n_int_points, mesh_o_nequ);
    
}

void TRMBuildTransfers::CreateTransferPressure_To_Mixed_V(TPZAutoPointer< TPZCompMesh> &cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &A){
    
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
    fp_To_Mixed.Resize(n_int_points, mesh_o_nequ);


    
}


void TRMBuildTransfers::CreateTransferGeometricData_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &Aweights, TPZVec<STATE> &Adet, TPZVec<STATE> &Arhs){
    
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
    int volumetric_elements = 0;
    
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
        volumetric_elements++;
        // Computing the global integration points indexes
        TPZManVector<long> globindexes;
        mf_cel->GetMemoryIndices(globindexes);
        int ndata = 1;
        long nglob_indexes = globindexes.size();
        for(long i=0; i< nglob_indexes; i++ )
        {
            long glob = globindexes[i];
            IA[glob+1] = ndata;
        }
        
    }
    
    for(long i=0; i< n_int_points; i++ )
    {
        IA[i+1] += IA[i];
    }
    
    JA.Resize(IA[n_int_points]);
    Aweights.Resize(IA[n_int_points]);
    Adet.Resize(IA[n_int_points]);
    Arhs.Resize(IA[n_int_points]);

    
    fIntegrationWeights_V.Resize(n_int_points, volumetric_elements);
    fJacobianDet_V.Resize(n_int_points, volumetric_elements);
    fRhs_V.Resize(n_int_points, volumetric_elements);
    
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
    TPZAutoPointer<TPZVec<long> > IA = new TPZVec<long>;
    TPZAutoPointer<TPZVec<long> > JA = new TPZVec<long>;
    TPZAutoPointer<TPZVec<double> > Ax = new TPZVec<double>;
    TPZAutoPointer<TPZVec<double> > Ay = new TPZVec<double>;
    TPZAutoPointer<TPZVec<double> > Az = new TPZVec<double>;
    TPZAutoPointer<TPZVec<double> > Ad = new TPZVec<double>;
    
    this->CreateTransferFlux_To_Mixed_V(cmesh_multiphysics, mesh_index, IA, JA, Ax, Ay, Az, Ad);
    
    TPZAutoPointer<TPZVec<long> > IAg = new TPZVec<long>;
    TPZAutoPointer<TPZVec<long> > JAg = new TPZVec<long>;
    TPZAutoPointer<TPZVec<double> >Aw = new TPZVec<double>;
    TPZAutoPointer<TPZVec<double> >Adet = new TPZVec<double>;
    TPZAutoPointer<TPZVec<double> >Arhs = new TPZVec<double>;

    this->CreateTransferGeometricData_To_Mixed_V(cmesh_multiphysics, mesh_index, IAg, JAg, Aw, Adet, Arhs);
    
    int origin = mesh_index;
    int destination = mesh_index;
    int volumetric_elements = 0;
    int iblock = 0;
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
        
        volumetric_elements++;

        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel_o->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel_o->NShapeF(),0.0);
        
        int nconnect = intel_o->NConnects();
        int iphicount = 0;
        // Compute all the phi values and push inside corresponding j-destination;
        int el_nshape = 0;
        for (int icon = 0; icon < nconnect; icon++) {
            
            TPZConnect  & con = intel_o->Connect(icon);
            long seqnumber = con.SequenceNumber();
            long position = cmesh_o->Block().Position(seqnumber);
            el_nshape += con.NShape();
            
        }
        
        // Creating the block with the right size
        TPZFMatrix<double> blockx(np_cel_d,el_nshape);
        TPZFMatrix<double> blocky(np_cel_d,el_nshape);
        TPZFMatrix<double> blockz(np_cel_d,el_nshape);
        TPZFMatrix<double> blockd(np_cel_d,el_nshape);
        
        for (int ip = 0; ip < np_cel_d ; ip++)
        {
            
            TPZManVector<REAL,3> qsi(3,0.0);
            STATE w;
            int_points_cel_d.Point(ip, qsi, w);
            
            // Indetifying the right global index
            long globindex = globindexes[ip];

            long JAcount = IA->operator[](globindex);
            long JgAcount = IAg->operator[](globindex);
            
            // Get the vectorial phi
            intel_o->Shape(qsi, phi, dphidxi);
            intel_o->InitMaterialData(data);
            intel_o->ComputeRequiredData(data,qsi);
            
            // Hdiv space
            REAL JacobianDet  = data.detjac;
            TPZFNMatrix<100,REAL> phiuH1         = data.phi;   // For H1  test functions Q
            TPZFNMatrix<300,STATE> dphiuH1       = data.dphi; // Derivative For H1  test functions
            TPZFNMatrix<300,STATE> dphiuH1axes   = data.dphix; // Derivative For H1  test functions
            
            TPZFNMatrix<300,STATE> Qaxes = data.axes;
            TPZFNMatrix<300,STATE> QaxesT;
            TPZFNMatrix<300,STATE> Jacobian = data.jacobian;
            TPZFNMatrix<300,STATE> JacobianInverse = data.jacinv;
            
//            TPZFMatrix<STATE> GradOfX;
            TPZFNMatrix<300,STATE> GradOfXInverse;
            TPZFNMatrix<3,STATE> VectorOnMaster;
            TPZFNMatrix<3,STATE> VectorOnXYZ(3,1,0.0);
            
            Qaxes.Transpose(&QaxesT);
//                                QaxesT.Multiply(Jacobian, GradOfX);
            JacobianInverse.Multiply(Qaxes, GradOfXInverse);
            
//            TPZManVector<STATE,1> fval(1,0.0);
//            ExactLaplacian(data.x,fval); // defined as static member
//            REAL rhs = fval[0];
//            
////            // Feed material memory
////            material_memory[globindex].SetWeight(w);
////            material_memory[globindex].SetDetJac(JacobianDet);
////            material_memory[globindex].SetX(data.x);
////            material_memory[globindex].SetRhs(rhs);
//////            std::cout << "globindex = " <<  globindex << std::endl;
////            // Feed vectors of geometric data for each volumetric element index
////            JAg->operator[](JgAcount) = volumetric_elements-1;
////            Aw->operator[](JgAcount) = w;
////            Adet->operator[](JgAcount) = JacobianDet;
////            Arhs->operator[](JgAcount) = rhs;
////            JgAcount++;
            
            int nconnect = intel_o->NConnects();
            int iphicount = 0;
            // Compute all the phi values and push inside corresponding j-destination;
            
            for (int icon = 0; icon < nconnect; icon++) {
                
                TPZConnect  & con = intel_o->Connect(icon);
                long seqnumber = con.SequenceNumber();
//                long position = cmesh_o->Block().Position(seqnumber);
                int nshape = con.NShape();
                
                for (int ish=0; ish < nshape; ish++) {
                    
                    int vector_index = data.fVecShapeIndex[iphicount].first;
                    int shape_index = data.fVecShapeIndex[iphicount].second;
                    
                    VectorOnXYZ(0,0) = data.fNormalVec(0,vector_index);
                    VectorOnXYZ(1,0) = data.fNormalVec(1,vector_index);
                    VectorOnXYZ(2,0) = data.fNormalVec(2,vector_index);
                    
                    GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
                    VectorOnMaster *= JacobianDet;
                    
//                    std::cout << "VectorOnMaster = " << VectorOnMaster << std::endl;
    
                    /* Contravariant Piola mapping preserves the divergence */
                    int dim =VectorOnMaster.Rows();
                    TPZFNMatrix<9,REAL> GradOfPhi(dim,dim);
                    for (int ir = 0; ir < dim; ir++) {
                        
                        //  Compute grad_{hat}(PhiH1)
                        
                        GradOfPhi(ir,0) = VectorOnMaster(ir,0)*dphiuH1(0,shape_index);
                        GradOfPhi(ir,1) = VectorOnMaster(ir,0)*dphiuH1(1,shape_index);
                        GradOfPhi(ir,2) = VectorOnMaster(ir,0)*dphiuH1(2,shape_index);
                        
                    }
                    GradOfPhi *= (1.0/JacobianDet);
                    REAL divPhi_Hdiv = ( GradOfPhi(0,0) + GradOfPhi(1,1) + GradOfPhi(2,2));
                    
                    TPZManVector<STATE> phi_Hdiv(3);
                    phi_Hdiv[0] = phi(shape_index,0)*VectorOnXYZ(0,0);
                    phi_Hdiv[1] = phi(shape_index,0)*VectorOnXYZ(1,0);
                    phi_Hdiv[2] = phi(shape_index,0)*VectorOnXYZ(2,0);
                    
//                    JA->operator[](JAcount) = position + ish;
//                    Ax->operator[](JAcount) = phi_Hdiv[0];
//                    Ay->operator[](JAcount) = phi_Hdiv[1];
//                    Az->operator[](JAcount) = phi_Hdiv[2];
//                    Ad->operator[](JAcount) = divPhi_Hdiv;
                    blockx(ip,iphicount) = phi_Hdiv[0];
                    blocky(ip,iphicount) = phi_Hdiv[1];
                    blockz(ip,iphicount) = phi_Hdiv[2];
                    blockd(ip,iphicount) = divPhi_Hdiv;
                    iphicount++;
//                    JAcount++;
                }
                
            }
        }

        fTransfer_X_Flux_To_Mixed_V.SetBlock(iblock, blockx);
        fTransfer_Y_Flux_To_Mixed_V.SetBlock(iblock, blocky);
        fTransfer_Z_Flux_To_Mixed_V.SetBlock(iblock, blockz);
        fTransferDivergenceTo_Mixed_V.SetBlock(iblock, blockd);
        iblock++;
    }
    
    
//    fTransfer_X_Flux_To_Mixed_V.
    
//    fTransfer_X_Flux_To_Mixed_V.SetData(IA, JA, Ax);
//    fTransfer_Y_Flux_To_Mixed_V.SetData(IA, JA, Ay);
//    fTransfer_Z_Flux_To_Mixed_V.SetData(IA, JA, Az);
//    fTransferDivergenceTo_Mixed_V.SetData(IA, JA, Ad);
//    
//    fIntegrationWeights_V.SetData(IAg,JAg,Aw);
//    fJacobianDet_V.SetData(IAg,JAg,Adet);
//    fRhs_V.SetData(IAg,JAg,Arhs);
    associated_material->SetMemory(material_memory);
}

/** @brief exact laplacian */
void TRMBuildTransfers::ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &f)
{
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    f[0] = 6;//2.*(-1. + x)*x*(-1. + y)*y + 2.*(-1. + x)*x*(-1. + z)*z + 2.*(-1. + y)*y*(-1. + z)*z;
}

void TRMBuildTransfers::ComputeTransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh > cmesh_multiphysics, int mesh_index){

#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    long nel = cmesh_multiphysics->NElements();
    //int n_data = 1; // scalar
    mesh_index = 0;
    long element_index = 0;
    int imat = 0;
    
    // For the imat
    
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    // Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> material_memory =  associated_material->GetMemory();
    //int n_int_points = material_memory.NElements();
    
    // Compute destination index scatter by element (Omega and Gamma)
    TPZManVector< TPZVec<long>, 1000000 > dof_scatter(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions(nel);
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        if (!cel) {
            DebugStop();
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        if(!mf_cel)
        {
            DebugStop();
        }
        
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Increassing new block
        element_index++;
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes;
        TPZManVector<long> dof_indexes;
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);
        dof_scatter[element_index] = dof_indexes;
        
        
        
    }
    
    
    
    // Initialize block matrix
    TRMIrregularBlockDiagonal<REAL> p_to_Mixed;
    

    
    
    
    // Resizing IA JA and A
    TPZAutoPointer<TPZVec<long> > IA = new TPZVec<long>;
    TPZAutoPointer<TPZVec<long> > JA = new TPZVec<long>;
    TPZAutoPointer<TPZVec<double> > A = new TPZVec<double>;
    
    
    

//    
//    this->CreateTransferPressure_To_Mixed_V(cmesh_multiphysics, mesh_index, IA, JA, A);
    
    int origin = mesh_index;
    int destination = mesh_index;
    
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        if (!cel) {
            DebugStop();
        }
        
        
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
           long JAcount = IA->operator[](globindex);

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
                    JA->operator[](JAcount) = position + ish;
                    A->operator[](JAcount) = phi(iphicount,0);
                    iphicount++;
                    JAcount++;
                }
            }
        }
        
    }
    
//    fTransferPressure_To_Mixed_V.SetData(IA, JA, A);
}

void TRMBuildTransfers::ElementDofIndexes(TPZInterpolationSpace * intel, TPZVec<long> &dof_indexes){

#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index;
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
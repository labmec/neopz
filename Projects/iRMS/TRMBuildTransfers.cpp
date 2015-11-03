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
    fTransferScalar_Vol.Resize(0, 0);
}

/** @brief Default desconstructor */
TRMBuildTransfers::~TRMBuildTransfers(){

}

/** @brief Compute the sparse matrix to transfer information between two cmesh that belongs to multiphysics mesh  */
void TRMBuildTransfers::ComputeTransferScalar_Vol(TPZCompMesh *cmesh_multiphysics, int origin, int destination){

    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
    
    long nel = cmesh_multiphysics->NElements();
    
    // Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(_ReservMatId);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> destination_memory =  associated_material->GetMemory();
    int np_cmesh_des = destination_memory.NElements();
    
    // Resizing IA JA and A
    TPZManVector<long> IA(np_cmesh_des+1,0);
    TPZManVector<long> JA;
    TPZManVector<double> A;
    
    TPZCompMesh * cmesh_o;
    TPZCompMesh * cmesh_d;
    int mesh_o_nequ = 0;

    for (long icel = 0; icel < nel; icel++) {

        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        if (!cel) {
            DebugStop();
        }

         TPZCondensedCompEl * mf_cel_condensed = dynamic_cast<TPZCondensedCompEl *> (cel);
 
         if(!mf_cel_condensed){
             DebugStop();
         }
    
        TPZCompEl * mf_cel_cond_ref = mf_cel_condensed->ReferenceCompEl();
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(mf_cel_cond_ref);
        if(!mf_cel)
        {
            DebugStop();
        }
        
        // filtrar os elementos que nao sao volumetricos.
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
        for(long i=0; i< globindexes.size(); i++ )
        {
            long glob = globindexes[i];
            IA[glob+1] = nshapes;
        }
        
    }    

    for(long i=0; i< np_cmesh_des; i++ )
    {
       IA[i+1] += IA[i];
    }
    
    JA.Resize(IA[np_cmesh_des]);
    A.Resize(IA[np_cmesh_des]);

    TPZFYsmpMatrix<REAL> TransferScalar(np_cmesh_des,mesh_o_nequ);
    
    long row = 0;
    for (long icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        if (!cel) {
            DebugStop();
        }
        
        TPZCondensedCompEl * mf_cel_condensed = dynamic_cast<TPZCondensedCompEl *> (cel);
        
        if(!mf_cel_condensed){
            DebugStop();
        }
        
        TPZCompEl * mf_cel_cond_ref = mf_cel_condensed->ReferenceCompEl();
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(mf_cel_cond_ref);
        if(!mf_cel)
        {
            DebugStop();
        }
        
        // filtrar os elementos que nao sao volumetricos.
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
        TPZCompMesh * cmesh_d = intel_d->Mesh();
        
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
        int el_dim =mf_cel->Reference()->Dimension();
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
//            row++;
//           IA[row] = JAcount;
        }
        
    }
    
    TransferScalar.SetData(IA, JA, A);
    fTransferScalar_Vol = TransferScalar;
    
    
    
}
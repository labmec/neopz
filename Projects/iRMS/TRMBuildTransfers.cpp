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

/** @brief Default constructor */
TRMBuildTransfers::TRMBuildTransfers(){
    fTransferScalar_Vol.Resize(0, 0);
}

/** @brief Default desconstructor */
TRMBuildTransfers::~TRMBuildTransfers(){
    
}

/** @brief Compute the  */
void TRMBuildTransfers::ComputeTransferScalar_Vol(TPZCompMesh *cmesh_origin, TPZCompMesh *cmesh_destination){
    
    TPZFYsmpMatrix<STATE> TransferScalar;
    if (!cmesh_origin || !cmesh_destination) {
        std::cout << "There is no computational mesh for cmesh_origin or cmesh_destination, cmesh = Null." << std::endl;
        DebugStop();
    }
    
    long nel_origin         = cmesh_origin->NElements();
    long nel_destination    = cmesh_destination->NElements();
    
    // Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_origin->FindMaterial(_ReservMatId);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    TPZAdmChunkVector<TRMMemory> destination_memory =  associated_material->GetMemory();
//    int np_cmesh_des = destination_memory.NElements();
    int np_cmesh_des = 8;
    
    // Resizing IA JA and A
    TPZVec<long> IA(np_cmesh_des,0.0);
    TPZVec<long> JA;
    TPZVec<double> A;

    long row = 0;
    for (long icel_d = 0; icel_d < nel_destination; icel_d++) {
        TPZCompEl * cel_d = cmesh_destination->Element(icel_d);
        if (!cel_d) {
            DebugStop();
        }
        
        TPZInterpolationSpace * intel_d = dynamic_cast<TPZInterpolationSpace * >(cel_d);
        if (!intel_d) {
            DebugStop();
        }
        
        // filtrar os elementos que nao sao volumetricos.
        
        
        // Computing the global integration points indexes
        TPZManVector<long> globindexes;
        cel_d->GetMemoryIndices(globindexes);
        
        // Computing the local integration points indexes
        TPZIntPoints & int_points_cel_d = intel_d->GetIntegrationRule();
        int np_cel_d = int_points_cel_d.NPoints();
        
        // Computing over all integration points of the compuational element cel
        for (int ip = 0; ip < np_cel_d ; ip++) {
            TPZVec<REAL> qsi(3,0.0);
            STATE w;
            int_points_cel_d.Point(ip, qsi, w);
            
            // Indetifying the right global index
//            long globindex = globindexes[ip];
//            long JAcount = IA[globindex];
            long globindex = ip;
            
            IA[globindex] = ip;
            
            TPZFMatrix<REAL> phi;
            TPZFMatrix<REAL> dphidxi;
            intel_d->Shape(qsi, phi, dphidxi);
            
            for (long icel_o = 0; icel_o < nel_origin; icel_o++) {
                TPZCompEl * cel_o = cmesh_origin->Element(icel_o);
                if (!cel_o) {
                    DebugStop();
                }
            
                int nconnect = cel_o->NConnects();
                int iphicount = 0;
                // Compute all the phi values and push inside corresponding j-destination;
                for (int icon = 0; icon < nconnect; icon++) {
                    TPZConnect  & con = cel_o->Connect(icon);
                    long seqnumber = con.SequenceNumber();
                    long position = cmesh_origin->Block().Position(seqnumber);
                    int nshape = con.NShape();
                    
                    for (int ish=0; ish < nshape; ish++) {
//                        JA[JAcount] = position + ish;
//                        A[JAcount] = phi(iphicount,0);
//                        iphicount++;
//                        JAcount++;
                    }
                }
            }
            row++;
//            IA[row] = JAcount;
        }
        
    }
    
    TransferScalar.SetData(IA, JA, A);
    fTransferScalar_Vol = TransferScalar;
    
}
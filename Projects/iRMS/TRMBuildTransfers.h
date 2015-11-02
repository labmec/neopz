//
//  TRMBuildTransfers.h
//  PZ
//
//  Created by Omar on 10/27/15.
//
//

#ifndef TRMBuildTransfers_h
#define TRMBuildTransfers_h

#include "pzysmp.h"
#include <stdio.h>
#include "pzgmesh.h"
#include "pzcmesh.h"



class TRMBuildTransfers{
    
public:
    
    /** @brief Default constructor */
    TRMBuildTransfers();
    
    /** @brief Default desconstructor */
    ~TRMBuildTransfers();

    /** @brief Compute the sparse matrix to transfer information between two cmesh that belongs to multiphysics mesh  */
    void ComputeTransferScalar_Vol(TPZCompMesh *cmesh_multiphysics, int origin, int destination);
    
    /** @brief Get the sparse matrix to transfer scalar solution to integrations points other mesh */
    TPZFYsmpMatrix<STATE> GetTransferScalar_Vol(){
        return fTransferScalar_Vol;
    }
    
private:
    
    /** @brief Sparse matrix to transfer scalar solution to integrations points other mesh */
    TPZFYsmpMatrix<STATE> fTransferScalar_Vol;
    
    
};


#endif /* TRMBuildTransfers_h */

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

    /**
     * @defgroup Do the transfers
     * @{
     */
    
    /** @brief Transfer Pressure to integration points of multiphysics mesh  */
    void TransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_pressure, TPZAutoPointer< TPZCompMesh > cmesh_multiphysics);
    
    /** @brief Transfer Flux to integration points of multiphysics mesh  */
    void TransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer< TPZCompMesh > cmesh_multiphysics);
    
    /** @brief Transfer Sw to integration points of multiphysics mesh  */
    void TransferSw_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_sw, TPZAutoPointer< TPZCompMesh > cmesh_multiphysics);
    
    /** @brief Transfer Sw to integration points of multiphysics mesh  */
    void TransferSo_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_so, TPZAutoPointer< TPZCompMesh > cmesh_multiphysics);
    
    // @}
    
    /**
     * @defgroup Create, compute and get matrices
     * @{
     */
    
    /** @brief Create the sparse matrix to transfer Pressure to multiphysics mesh  */
    void CreateTransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZManVector<long> &IA, TPZManVector<long> &JA, TPZManVector<STATE> &A);
    
    /** @brief Create the sparse matrix to transfer the total flux to multiphysics mesh  */
    void CreateTransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZManVector<long> &IA, TPZManVector<long> &JA, TPZManVector<STATE> &Ax, TPZManVector<STATE> &Ay, TPZManVector<STATE> &Az,  TPZManVector<STATE> Ad);
    
    /** @brief Compute the sparse matrix to transfer Pressure to multiphysics mesh  */
    void ComputeTransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Compute the sparse matrix to transfer the total flux to multiphysics mesh  */
    void ComputeTransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZFYsmpMatrix<STATE> GetTransfer_X_Flux_To_Mixed(){
        return fTransfer_X_Flux_To_Mixed;
    }

    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZFYsmpMatrix<STATE> GetTransfer_Y_Flux_To_Mixed(){
        return fTransfer_Y_Flux_To_Mixed;
    }
    
    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZFYsmpMatrix<STATE> GetTransfer_Z_Flux_To_Mixed(){
        return fTransfer_Z_Flux_To_Mixed;
    }
    
    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZFYsmpMatrix<STATE> GetTransferDivergenceTo_Mixed(){
        return fTransferDivergenceTo_Mixed;
    }
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TPZFYsmpMatrix<STATE> GetTransferPressure_To_Mixed(){
        return fTransferPressure_To_Mixed;
    }

    
    // @}
    
private:
    
    /**
     * @defgroup Sparses Matrices to transfer information to the mixed problem
     * @{
     */
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of the mixed mesh */
    TPZFYsmpMatrix<REAL> fTransferPressure_To_Mixed;
    
    /** @brief Sparse matrix to transfer x-Flux solution to integrations points of the mixed mesh */
    TPZFYsmpMatrix<REAL> fTransfer_X_Flux_To_Mixed;
    
    /** @brief Sparse matrix to transfer y-Flux solution to integrations points of the mixed mesh */
    TPZFYsmpMatrix<REAL> fTransfer_Y_Flux_To_Mixed;
    
    /** @brief Sparse matrix to transfer z-Flux solution to integrations points of the mixed mesh */
    TPZFYsmpMatrix<REAL> fTransfer_Z_Flux_To_Mixed;
    
    /** @brief Sparse matrix to transfer divergence of flux solution to integrations points of the mixed mesh */
    TPZFYsmpMatrix<REAL> fTransferDivergenceTo_Mixed;
    
    /** @brief Sparse matrix to transfer Sw solution to integrations points of the mixed mesh */
    TPZFYsmpMatrix<REAL> fTransferSw_To_Mixed;
    
    /** @brief Sparse matrix to transfer So solution to integrations points of the mixed mesh */
    TPZFYsmpMatrix<REAL> fTransferSo_To_Mixed;
    
    // @}
    
    /**
     * @defgroup Sparses Matrices ot transfer information to saturations meshes
     * @{
     */
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    TPZFYsmpMatrix<REAL> fTransferPressure_To_Sw;
    
    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    TPZFYsmpMatrix<REAL> fTransferFlux_To_Sw;
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    TPZFYsmpMatrix<REAL> fTransferPressure_To_So;
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    TPZFYsmpMatrix<REAL> fTransferFlux_To_So;
    
    // @}
    
    /** @brief exact laplacian */
    void ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &f);
    
};


#endif /* TRMBuildTransfers_h */

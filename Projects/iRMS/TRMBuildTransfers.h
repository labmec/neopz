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
    void CreateTransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZManVector<long> &IA, TPZManVector<long> &JA, TPZManVector<STATE> &A);
    
    /** @brief Compute the sparse matrix to transfer Pressure to multiphysics mesh  */
    void ComputeTransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Compute the sparse matrix to transfer the total flux to multiphysics mesh  */
    void ComputeTransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer the total flux to multiphysics mesh  */
    TPZFYsmpMatrix<STATE> GetTransferFlux_To_Mixed(){
        return fTransferFlux_To_Mixed;
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
    
    /** @brief Sparse matrix to transfer Flux solution to integrations points of the mixed mesh */
    TPZFYsmpMatrix<REAL> fTransferFlux_To_Mixed;
    
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
    
};


#endif /* TRMBuildTransfers_h */

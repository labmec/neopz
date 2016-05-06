//
//  TRMBuildTransfers.h
//  PZ
//
//  Created by Omar on 10/27/15.
//  This class storage approximation space in global integration points arrays for volumetric and boundary elements
//

#ifndef TRMBuildTransfers_h
#define TRMBuildTransfers_h

#include "pzysmp.h"
#include "pzblockdiag.h"
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
    
    /** @brief Transfer Pressure to integration points of multiphysics mesh over volumetric elements */
    void TransferPressure_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_pressure, TPZAutoPointer< TPZCompMesh > cmesh_multiphysics);
    
    /** @brief Transfer Flux to integration points of multiphysics mesh  */
    void TransferFlux_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_flux, TPZAutoPointer< TPZCompMesh > cmesh_multiphysics);
    
    /** @brief Transfer Sw to integration points of multiphysics mesh  */
    void TransferSw_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_sw, TPZAutoPointer< TPZCompMesh > cmesh_multiphysics);
    
    /** @brief Transfer Sw to integration points of multiphysics mesh  */
    void TransferSo_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_so, TPZAutoPointer< TPZCompMesh > cmesh_multiphysics);
    
    // @}
    
    /**
     * @defgroup Create, compute and get matrices
     * @{
     */
    
    /** @brief Create the sparse matrix to transfer Pressure to multiphysics mesh over volumetric elements */
    void CreateTransferPressure_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &A);
    
    /** @brief Create the sparse matrix to transfer Geometric information to integration points over volumetric elements */
    void CreateTransferGeometricData_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &Aweights, TPZVec<STATE> &Adet, TPZVec<STATE> &Arhs);
    
    /** @brief Create the sparse matrix to transfer the total flux to multiphysics mesh over volumetric elements */
    void CreateTransferFlux_To_Mixed_V(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &Ax, TPZVec<STATE> &Ay, TPZVec<STATE> &Az,  TPZVec<STATE> &Ad);
    
    /** @brief Create the sparse matrix to transfer Pressure to multiphysics mesh over surface elements */
    void CreateTransferPressure_To_Mixed_S(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &A);
    
    /** @brief Create the sparse matrix to transfer Geometric information to integration points  */
    void CreateTransferGeometricData_To_Mixed_S(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &Aweights, TPZVec<STATE> &Adet, TPZVec<STATE> &Arhs);
    
    /** @brief Create the sparse matrix to transfer the total flux to multiphysics mesh  */
    void CreateTransferFlux_To_Mixed_S(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index, TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<STATE> &Ax, TPZVec<STATE> &Ay, TPZVec<STATE> &Az,  TPZVec<STATE> &Ad);
    
    /** @brief Compute the sparse matrix to transfer Pressure to multiphysics mesh  */
    void ComputeTransferPressure_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Compute the sparse matrix to transfer the total flux to multiphysics mesh  */
    void ComputeTransferFlux_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZBlockDiagonal<STATE> GetTransfer_X_Flux_To_Mixed_V(){
        return fTransfer_X_Flux_To_Mixed_V;
    }

    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZBlockDiagonal<STATE> GetTransfer_Y_Flux_To_Mixed_V(){
        return fTransfer_Y_Flux_To_Mixed_V;
    }
    
    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZBlockDiagonal<STATE> GetTransfer_Z_Flux_To_Mixed_V(){
        return fTransfer_Z_Flux_To_Mixed_V;
    }
    
    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZBlockDiagonal<STATE> GetTransferDivergenceTo_Mixed_V(){
        return fTransferDivergenceTo_Mixed_V;
    }
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TPZBlockDiagonal<STATE> GetTransferPressure_To_Mixed_V(){
        return fTransferPressure_To_Mixed_V;
    }
    
    /** @brief Get the sparse matrix to transfer the total x-flux to multiphysics mesh  */
    TPZBlockDiagonal<STATE> GetWeightsTo_Mixed_V(){
        return fIntegrationWeights_V;
    }
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TPZBlockDiagonal<STATE> GetJacobianDet_To_Mixed_V(){
        return fJacobianDet_V;
    }
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TPZBlockDiagonal<STATE> GetRhs_To_Mixed_V(){
        return fRhs_V;
    }

    
    // @}
    
private:
    
    /**
     * @defgroup Sparses Matrices to transfer information to the mixed problem on Omega
     * @{
     */
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransferPressure_To_Mixed_V;
    
    /** @brief Sparse matrix to transfer x-Flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransfer_X_Flux_To_Mixed_V;
    
    /** @brief Sparse matrix to transfer y-Flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransfer_Y_Flux_To_Mixed_V;
    
    /** @brief Sparse matrix to transfer z-Flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransfer_Z_Flux_To_Mixed_V;
    
    /** @brief Sparse matrix to transfer divergence of flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransferDivergenceTo_Mixed_V;
    
    /** @brief Sparse matrix to transfer Sw solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransferSw_To_Mixed_V;
    
    /** @brief Sparse matrix to transfer So solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransferSo_To_Mixed_V;
    
    // @}
    
    /**
     * @defgroup Sparses Matrices to transfer information to saturations meshes on Omega
     * @{
     */
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    TPZBlockDiagonal<REAL> fTransferPressure_To_Sw_V;
    
    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    TPZBlockDiagonal<REAL> fTransferFlux_To_Sw_V;
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    TPZBlockDiagonal<REAL> fTransferPressure_To_So_V;
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    TPZBlockDiagonal<REAL> fTransferFlux_To_So_V;
    
    // @}
    
    /**
     * @defgroup Sparses Matrices to transfer information to the mixed problem on Omega
     * @{
     */
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransferPressure_To_Mixed_S;
    
    /** @brief Sparse matrix to transfer x-Flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransfer_X_Flux_To_Mixed_S;
    
    /** @brief Sparse matrix to transfer y-Flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransfer_Y_Flux_To_Mixed_S;
    
    /** @brief Sparse matrix to transfer z-Flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransfer_Z_Flux_To_Mixed_S;
    
    /** @brief Sparse matrix to transfer divergence of flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransferDivergenceTo_Mixed_S;
    
    /** @brief Sparse matrix to transfer Sw solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransferSw_To_Mixed_S;
    
    /** @brief Sparse matrix to transfer So solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<REAL> fTransferSo_To_Mixed_S;
    
    // @}
    
    /**
     * @defgroup Sparses Matrices to transfer information to saturations meshes on Omega
     * @{
     */
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    TPZBlockDiagonal<REAL> fTransferPressure_To_Sw_S;
    
    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    TPZBlockDiagonal<REAL> fTransferFlux_To_Sw_S;
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    TPZBlockDiagonal<REAL> fTransferPressure_To_So_S;
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    TPZBlockDiagonal<REAL> fTransferFlux_To_So_S;
    
    // @}
    
    /**
     * @defgroup Vectors to transfer geometric information to integration points at volumes on Omega
     * @{
     */
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    TPZBlockDiagonal<REAL> fIntegrationWeights_V;
    
    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    TPZBlockDiagonal<REAL> fJacobianDet_V;
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    TPZBlockDiagonal<REAL> fRhs_V;
    
    
    // @}
    
    /**
     * @defgroup Vectors to transfer geometric information to integration points at volumes on Gamma
     * @{
     */
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    TPZBlockDiagonal<REAL> fIntegrationWeights_S;
    
    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    TPZBlockDiagonal<REAL> fJacobianDet_S;
    
    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    TPZBlockDiagonal<REAL> fRhs_S;
    
    
    // @}
    
    /** @brief exact laplacian */
    void ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &f);
    
};


#endif /* TRMBuildTransfers_h */

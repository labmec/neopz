//
//  TRMBuildTransfers.h
//  PZ
//
//  Created by Omar on 10/27/15.
//  This class storage approximation space in global integration points arrays for volumetric and boundary elements
//

#ifndef TRMBuildTransfers_h
#define TRMBuildTransfers_h

#include <stdio.h>

#include "tpzintpoints.h"
#include "pzmatwithmem.h"
#include "TRMMemory.h"
#include "TRMMixedDarcy.h"
#include "pzmaterial.h"
#include "TRMFlowConstants.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"
#include "pzcondensedcompel.h"


#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZInterfaceEl.h"

#include "pzysmp.h"
#include "pzblockdiag.h"
#include "TRMSimulationData.h"
#include "TRMIrregularBlockDiagonal.h"


class TRMBuildTransfers{
    
    
private:
    
    
    
    /**
     * @defgroup Attributes
     * @{
     */
    
    /** @brief Autopointer of simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    // @}
    
    
    /**
     * @defgroup Methods for change the attributes
     * @{
     */
    
    /** @brief Compute element dof indexes */
    void ElementDofIndexes(TPZInterpolationSpace * &intel,  TPZVec<long> &dof_indexes);

    /** @brief Compute element dof indexes at faces */    
    void ElementDofFaceIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes);
    
    // @}
    
    /** @brief Autopointer of simulation data */
    
    /**
     * @defgroup Sparses Matrices to transfer information to the mixed problem on Omega
     * @{
     */
    
    /** @brief Diagonal block matrix to transfer u flux solution to integrations points of the mixed mesh */
    TRMIrregularBlockDiagonal<STATE> fu_To_Mixed;
    
    /** @brief flux dof indexes per element */
    TPZVec< TPZVec<long> > fu_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer Pressure solution to integrations points of the mixed mesh */
    TRMIrregularBlockDiagonal<STATE> fp_To_Mixed;
    
    /** @brief Diagonal block matrix to transfer Average pressure solution to integrations points of the mixed mesh */
    TRMIrregularBlockDiagonal<STATE> fap_To_Mixed;
    
    /** @brief pressure dof indexes per element */
    TPZVec< TPZVec<long> > fp_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer Average alpha saturation solution to integrations points of the mixed mesh */
    TRMIrregularBlockDiagonal<STATE> fs_a_To_Mixed;
    
    /** @brief pressure dof indexes per element */
    TPZVec< TPZVec<long> > fs_a_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer Average normal flux solution to integrations points of the transport mesh */
    TRMIrregularBlockDiagonal<STATE> fun_To_Transport_a;
    
    /** @brief pressure dof indexes per element */
    TPZVec< TPZVec<long> > fun_dof_scatter;

    /** @brief left and right geometric element indexes */
    TPZStack < std::pair<long, long> > fleft_right_indexes;
    
    /** @brief left and right geometric element indexes */
    TPZStack < long > finterface_indexes;
    
    //    /** @brief Sparse matrix to transfer x-Flux solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransfer_X_Flux_To_Mixed_V;
    //
    //    /** @brief Sparse matrix to transfer y-Flux solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransfer_Y_Flux_To_Mixed_V;
    //
    //    /** @brief Sparse matrix to transfer z-Flux solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransfer_Z_Flux_To_Mixed_V;
    //
    //    /** @brief Sparse matrix to transfer divergence of flux solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransferDivergenceTo_Mixed_V;
    //
    //    /** @brief Sparse matrix to transfer Sw solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransferSw_To_Mixed_V;
    //
    //    /** @brief Sparse matrix to transfer So solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransferSo_To_Mixed_V;
    
    // @}
    
    /**
     * @defgroup Sparses Matrices to transfer information to saturations meshes on Omega
     * @{
     */
    
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    //    TPZBlockDiagonal<REAL> fTransferPressure_To_Sw_V;
    //
    //    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    //    TPZBlockDiagonal<REAL> fTransferFlux_To_Sw_V;
    //
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    //    TPZBlockDiagonal<REAL> fTransferPressure_To_So_V;
    //
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    //    TPZBlockDiagonal<REAL> fTransferFlux_To_So_V;
    
    // @}
    
    /**
     * @defgroup Sparses Matrices to transfer information to the mixed problem on Omega
     * @{
     */
    
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransferPressure_To_Mixed_S;
    //
    //    /** @brief Sparse matrix to transfer x-Flux solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransfer_X_Flux_To_Mixed_S;
    //
    //    /** @brief Sparse matrix to transfer y-Flux solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransfer_Y_Flux_To_Mixed_S;
    //
    //    /** @brief Sparse matrix to transfer z-Flux solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransfer_Z_Flux_To_Mixed_S;
    //
    //    /** @brief Sparse matrix to transfer divergence of flux solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransferDivergenceTo_Mixed_S;
    //
    //    /** @brief Sparse matrix to transfer Sw solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransferSw_To_Mixed_S;
    //
    //    /** @brief Sparse matrix to transfer So solution to integrations points of the mixed mesh */
    //    TPZBlockDiagonal<REAL> fTransferSo_To_Mixed_S;
    
    // @}
    
    /**
     * @defgroup Sparses Matrices to transfer information to saturations meshes on Omega
     * @{
     */
    
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    //    TPZBlockDiagonal<REAL> fTransferPressure_To_Sw_S;
    //
    //    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    //    TPZBlockDiagonal<REAL> fTransferFlux_To_Sw_S;
    //
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    //    TPZBlockDiagonal<REAL> fTransferPressure_To_So_S;
    //
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    //    TPZBlockDiagonal<REAL> fTransferFlux_To_So_S;
    
    // @}
    
    /**
     * @defgroup Vectors to transfer geometric information to integration points at volumes on Omega
     * @{
     */
    
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    //    TPZBlockDiagonal<REAL> fIntegrationWeights_V;
    //
    //    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    //    TPZBlockDiagonal<REAL> fJacobianDet_V;
    //
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    //    TPZBlockDiagonal<REAL> fRhs_V;
    
    
    // @}
    
    /**
     * @defgroup Vectors to transfer geometric information to integration points at volumes on Gamma
     * @{
     */
    
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of water saturation mesh */
    //    TPZBlockDiagonal<REAL> fIntegrationWeights_S;
    //
    //    /** @brief Sparse matrix to transfer flux solution to integrations points of the water saturation mesh */
    //    TPZBlockDiagonal<REAL> fJacobianDet_S;
    //    
    //    /** @brief Sparse matrix to transfer Pressure solution to integrations points of oil saturation mesh */
    //    TPZBlockDiagonal<REAL> fRhs_S;
    
    
    // @}
    
    //    /** @brief exact laplacian */
    //    void ExactLaplacian(const TPZVec<REAL> &pt, TPZVec<STATE> &f);
    
    
public:
    
    /** @brief Default constructor */
    TRMBuildTransfers();
    
    /** @brief Default desconstructor */
    ~TRMBuildTransfers();
    
    /** @brief Copy constructor $ */
    TRMBuildTransfers(const TRMBuildTransfers &copy);
    
    /** @brief Copy assignemnt operator $ */
    TRMBuildTransfers &operator=(const TRMBuildTransfers &other);
    

    /**
     * @defgroup Apply transfers to different meshes
     * @{
     */

    /** @brief Transfer Flux to integration points of multiphysics mesh over volumetric elements */
    void Transfer_u_To_Mixed_Memory(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer Pressure to integration points of multiphysics mesh over volumetric elements */
    void Transfer_p_To_Mixed_Memory(TPZCompMesh * cmesh_pressure, TPZCompMesh * cmesh_multiphysics);
    
    // @}
    
    /**
     * @defgroup Create, compute and get transfer matrices
     * @{
     */

    /** @brief Initializate diagonal block matrix to transfer flux to multiphysics mesh  */
    void Initialize_u_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Initializate diagonal block matrix to transfer flux to multiphysics mesh  */
    void Fill_u_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_u_To_Mixed(){
        return fu_To_Mixed;
    }
    
    /** @brief Initializate  diagonal block matrix to transfer Pressure to multiphysics mesh  */
    void Initialize_p_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Initializate diagonal block matrix to transfer Pressure to multiphysics mesh  */
    void Fill_p_To_Mixed(TPZAutoPointer< TPZCompMesh> cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_p_To_Mixed(){
        return fp_To_Mixed;
    }
        
    /** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
    void Initialize_un_To_Transport_a(TPZAutoPointer< TPZCompMesh> flux_mesh, TPZAutoPointer< TPZCompMesh> transport_mesh);
    
    /** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
    void Fill_un_To_Transport_a(TPZAutoPointer< TPZCompMesh> flux_mesh, TPZAutoPointer< TPZCompMesh> transport_mesh);
    
    /** @brief Get the sparse matrix to transfer average normal flux solution to integrations points of the transport mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_un_To_Transport_a(){
        return fun_To_Transport_a;
    }
        

    
    // @}
    
    
    /**
     * @defgroup Fill class attributes
     * @{
     */
    
    /** @brief Compute left and right geometric element indexes associated with the transport mesh */
    void ComputeLeftRight(TPZAutoPointer< TPZCompMesh> transport_mesh);
    
    
    
    // @}
    

    
    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TPZAutoPointer<TRMSimulationData> &SimulationData){
        fSimulationData = SimulationData;
    }
    
    /** @brief Get autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData>  SimulationData(){
        return fSimulationData;
    }
    
};


#endif /* TRMBuildTransfers_h */

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
#include "TPZMatWithMem.h"
#include "TRMMemory.h"
#include "TRMPhaseMemory.h"
#include "TRMPhaseInterfaceMemory.h"
#include "TRMMixedDarcy.h"
#include "TPZMaterial.h"
#include "TRMFlowConstants.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"
#include "pzcondensedcompel.h"


#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

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
    TRMSimulationData * fSimulationData;
    
    // @}

    
    /** @brief Autopointer of simulation data */
    
    /**
     * @defgroup Sparses Matrices to transfer information to the mixed problem on Omega
     * @{
     */
    
    /** @brief Diagonal block matrix to transfer u flux solution to integrations points of the mixed mesh */
    TRMIrregularBlockDiagonal<STATE> fu_To_Mixed;
    
    /** @brief flux dof indexes per element */
    TPZVec< TPZVec<int64_t> > fu_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer Pressure solution to integrations points of the mixed mesh */
    TRMIrregularBlockDiagonal<STATE> fp_To_Mixed;
    
    /** @brief pressure dof indexes per element */
    TPZVec< TPZVec<int64_t> > fp_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer saturation solution to integrations points of the transport mesh */
    TRMIrregularBlockDiagonal<STATE> fs_To_Transport;
    
//    /** @brief saturations dof indexes per element */
//    TPZVec< TPZVec<int64_t> > fs_dof_scatter;
    
    /** @brief saturation a dof indexes per element */
    TPZVec< TPZVec<int64_t> > fsa_dof_scatter;
    
    /** @brief staruation b dof indexes per element */
    TPZVec< TPZVec<int64_t> > fsb_dof_scatter;
    
//    /** @brief Diagonal block matrix to transfer Average alpha saturation solution to integrations points of the mixed mesh */
//    TRMIrregularBlockDiagonal<STATE> fs_a_To_Mixed;
    
//    /** @brief pressure dof indexes per element */
//    TPZVec< TPZVec<int64_t> > fs_a_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer Average normal flux solution to integrations points of the transport mesh over gamma */
    TRMIrregularBlockDiagonal<STATE> fun_To_Transport_gamma;
    
    /** @brief Diagonal block matrix to transfer Average normal flux solution to integrations points of the transport mesh over Gamma */
    TRMIrregularBlockDiagonal<STATE> fun_To_Transport_Gamma;
    
    /** @brief normal flux dof indexes per interface element on gamma (inner interfaces)*/
    TPZVec< TPZVec<int64_t> > fun_dof_scatter_gamma;
    
    /** @brief normal flux dof indexes per interface element on Gamma (boundary interfaces) */
    TPZVec< TPZVec<int64_t> > fun_dof_scatter_Gamma;

    /** @brief mixed and transpor computational multiphysics element indexes, every element is indexed by geometric element */
    TPZStack< std::pair<int64_t, std::pair<int64_t, int64_t> >  > fmixed_transport_cindexes;
    
    /** @brief left and right geometric element indexes on gamma */
    TPZStack < std::pair<int64_t, int64_t> > fleft_right_g_indexes_gamma;
    
    /** @brief geometric interface element indexes on Gamma */
    TPZStack < int64_t > finterface_g_indexes_gamma;
    
    /** @brief left and right geometric element indexes on Gamma */
    TPZStack < std::pair<int64_t, int64_t> > fleft_right_g_indexes_Gamma;
    
    /** @brief geometric interface element indexes on Gamma */
    TPZStack < int64_t > finterface_g_indexes_Gamma;
    
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
    void u_To_Mixed_Memory(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer Pressure to integration points of multiphysics mesh over volumetric elements */
    void p_To_Mixed_Memory(TPZCompMesh * cmesh_pressure, TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer saturations to integration points of multiphysics transport mesh over volumetric elements */
    void s_To_Transport_Memory(TPZCompMesh * cmesh_saturation, TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Transfer average pressure to integration points of multiphysics mixed meshes over volumetric elements */
    void p_avg_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed);
    
    /** @brief Transfer average quantities to integration points of multiphysics mixed/ transpor meshes over volumetric elements */
    void Reciprocal_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_trans);
    
    /** @brief Transfer normal fluxes to integration points of transport meshes */
    void un_To_Transport_Mesh(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_transport, bool IsBoundaryQ);
    
    // @}
    
    /**
     * @defgroup Create, compute and get transfer matrices
     * @{
     */

    /** @brief Initializate diagonal block matrix to transfer flux to multiphysics mesh  */
    void Initialize_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Initializate diagonal block matrix to transfer flux to multiphysics mesh  */
    void Fill_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_u_To_Mixed(){
        return fu_To_Mixed;
    }
    
    /** @brief Initializate  diagonal block matrix to transfer Pressure to multiphysics mesh  */
    void Initialize_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Initializate diagonal block matrix to transfer Pressure to multiphysics mesh  */
    void Fill_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_p_To_Mixed(){
        return fp_To_Mixed;
    }
    
    /** @brief Initializate  diagonal block matrix to transfer Pressure to multiphysics mesh  */
    void Initialize_s_To_Transport(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Initializate diagonal block matrix to transfer Pressure to multiphysics mesh  */
    void Fill_s_To_Transport(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer saturation to multiphysics mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_s_To_Transport(){
        return fs_To_Transport;
    }
        
    /** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh over Gamma or gamma */
    void Initialize_un_To_Transport(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ);
    
    /** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh over Gamma or gamma */
    void Fill_un_To_Transport(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ);
    
    
    /** @brief Get the sparse matrix to transfer average normal flux solution to integrations points of the transport mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_un_To_Transport_gamma(){
        return fun_To_Transport_gamma;
    }
    
    /** @brief Get the sparse matrix to transfer average normal flux solution to integrations points of the transport mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_un_To_Transport_Gamma(){
        return fun_To_Transport_Gamma;
    }
    
    // @}
    
    
    /**
     * @defgroup Fill class attributes
     * @{
     */

    /** @brief Compute left and right geometric element indexes associated with the transport mesh */
    void ComputeLeftRight(TPZCompMesh * transport_mesh);
    
    
    // @}
    
    
    
    /**
     * @defgroup Utility Methods
     * @{
     */
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Memory operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Get Global integration point indexes associaded  */
    void GlobalPointIndexes(TPZCompEl * cel, TPZManVector<int64_t,30> &int_point_indexes);

    /** @brief Get Global integration point indexes associaded with interfaces */
    void GlobalPointIndexesInterface(TPZCompEl * int_cel, TPZManVector<int64_t,30> &int_point_indexes);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Computational mesh operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
    void FillComputationalElPairs(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Computational element operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Compute element dof indexes */
    void ElementDofIndexes(TPZInterpolationSpace * &intel,  TPZVec<int64_t> &dof_indexes);
    
    /** @brief Compute element dof indexes at given connect */
    void ElementDofFaceIndexes(int connect,TPZInterpolationSpace * &intel, TPZVec<int64_t> &dof_indexes);
    
    /** @brief Compute element dof indexes at given connect */
    void ElementDofFaceIndexes(int connect,TPZMultiphysicsElement * &m_el, TPZVec<int64_t> &dof_indexes);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Geometry Operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Ifdentify the which side of the volume is associated with the face */
    bool IdentifyFace(int &side, TPZGeoEl * vol, TPZGeoEl * face);
    
    /** @brief Dimensionla Measure of the elemnt */
    REAL DimensionalMeasure(TPZGeoEl * gel);
    
    /** @brief Compute indices associated to faces on 3D topologies */
    void ComputeFaceIndex(TPZGeoEl * gel , TPZVec<int> &sides);
    
    /** @brief Compute sides associated to faces on 3D topologies */
    void ComputeFaceNormals(TPZGeoEl * gel , TPZVec<int> &sides, TPZFMatrix<STATE> &normals);
    
    
    // @}
    

    
    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData){
        fSimulationData = SimulationData;
    }
    
    /** @brief Get autopointer of Simulation data */
    TRMSimulationData * SimulationData(){
        return fSimulationData;
    }
    
};


#endif /* TRMBuildTransfers_h */

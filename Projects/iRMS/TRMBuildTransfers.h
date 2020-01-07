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
#include "TRMPhaseMemory.h"
#include "TRMPhaseInterfaceMemory.h"
#include "TRMMixedDarcy.h"
#include "pzmaterial.h"
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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Transfers:: Iterative Coupling by Operator Splitting //////////////////////////////

    /** @brief Autopointer of simulation data */
    TRMSimulationData * fSimulationData;

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Attributes Elliptic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief full order u dof indexes per element */
    TPZVec< TPZVec<long> > fu_dof_scatter;
    
    /** @brief full order u dof indexes per element to parabolic */
    TPZVec< TPZVec<long> > fu_p_dof_scatter;

    /** @brief integration point indexes geo_intp_o_intp_t */
    TPZStack< std::pair< long, long >  > fe_e_cindexes;
    
    /** @brief integration point indexes geo_cel_o_cel_t */
    TPZStack< std::pair<long, TPZVec<long> >  > fe_e_intp_indexes;
    
    /** @brief integration point indexes geo_intp_o_intp_t */
    TPZStack< std::pair<long, std::pair<long, long> >  > fe_p_cindexes;
    
    /** @brief integration point indexes geo_intp_o_intp_t */
    TPZStack< std::pair<long, std::pair<long, long> >  > fe_h_cindexes;
    
    /** @brief integration point indexes geo_cel_o_cel_t */
    TPZStack< std::pair<long, std::pair< TPZVec<long>, TPZVec<long> > >  > fe_p_intp_indexes;
    
    /** @brief linear application u to elliptic mesh */
    TRMIrregularBlockDiagonal<STATE> fu_To_elliptic;
    
    /** @brief linear application grad_u to elliptic mesh */
    TRMIrregularBlockDiagonal<STATE> fgrad_u_To_elliptic;

    /** @brief linear application u to parabolic mesh */
    TRMIrregularBlockDiagonal<STATE> fu_To_parabolic;
    
    /** @brief linear application grad_u to parabolic mesh */
    TRMIrregularBlockDiagonal<STATE> fgrad_u_To_parabolic;
    
    /** @brief linear application grad_u to hyperbolic mesh */
    TRMIrregularBlockDiagonal<STATE> fgrad_u_To_hyperbolic;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Attributes Parabolic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief pressure dof indexes per element */
    TPZVec< TPZVec<long> > fp_dof_scatter;
    
    /** @brief pressure dof indexes per element to elliptic */
    TPZVec< TPZVec<long> > fp_e_dof_scatter;
    
    /** @brief total velocity dof indexes per element */
    TPZVec< TPZVec<long> > fq_dof_scatter;

    /** @brief integration point indexes geo_intp_o_intp_t */
    TPZStack< std::pair<long, long >  > fp_p_cindexes;
    
    /** @brief integration point indexes geo_cel_o_cel_t */
    TPZStack< std::pair<long, std::pair< TPZVec<long>, TPZVec<long> > >  > fp_p_intp_indexes;
    
    /** @brief integration point indexes geo_intp_o_intp_t */
    TPZStack< std::pair<long, std::pair<long, long> >  > fp_e_cindexes;

    /** @brief integration point indexes geo_cel_o_cel_t */
    TPZStack< std::pair<long, std::pair< TPZVec<long>, TPZVec<long> > >  > fp_e_intp_indexes;
    
    /** @brief linear application p to parabolic mesh */
    TRMIrregularBlockDiagonal<STATE> fp_To_parabolic;
    
    /** @brief linear application q to parabolic mesh */
    TRMIrregularBlockDiagonal<STATE> fq_To_parabolic;
    
    /** @brief linear application div_q to parabolic mesh */
    TRMIrregularBlockDiagonal<STATE> fdiv_q_To_parabolic;
    
    /** @brief linear application p to elliptic mesh */
    TRMIrregularBlockDiagonal<STATE> fp_To_elliptic;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Attributes Hyperbolic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /** @brief s_w dof indexes per element */
    TPZVec< TPZVec<long> > fsw_dof_scatter;
    
    /** @brief s_o dof indexes per element */
    TPZVec< TPZVec<long> > fso_dof_scatter;
    
    /** @brief p_avg dof indexes per element */
    TPZVec< TPZVec<long> > fp_avg_dof_scatter;
    
    /** @brief p_avg dof indexes per element */
    TPZVec< TPZVec<long> > fu_avg_dof_scatter;
    
    /** @brief p_avg dof indexes per element */
    TPZVec< TPZVec<long> > fqn_avg_dof_scatter_Gamma;
    
    /** @brief p_avg dof indexes per element */
    TPZVec< TPZVec<long> > fqn_avg_dof_scatter_gamma;
    
    /** @brief sw_avg dof indexes per element */
    TPZVec< TPZVec<long> > fsw_avg_dof_scatter;

    /** @brief integration point indexes geo_intp_o_intp_t */
    TPZStack< std::pair<long, long >  > fh_h_cindexes;
    
    /** @brief integration point indexes geo_cel_o_cel_t */
    TPZStack< std::pair<long, TPZVec<long> >  > fh_h_intp_indexes;
    
    /** @brief integration point indexes geo_intp_o_intp_t */
    TPZStack< std::pair<long, std::pair<long, long> >  > fh_p_cindexes;
    
    /** @brief integration point indexes geo_cel_o_cel_t */
    TPZStack< std::pair<long, std::pair< TPZVec<long>, TPZVec<long> > >  > fh_p_intp_indexes;
    
    /** @brief integration point indexes geo_cel_o_cel_t */
    TPZStack< std::pair<long, std::pair<long, std::vector<long> > >  > fparabolic_hyperbolic_cel_pairs;
    
    /** @brief integration point indexes geo_cel_o_cel_t */
    TPZStack< std::pair<long, std::pair<long, std::vector<long> > >  > felliptic_hyperbolic_cel_pairs;
    
    /** @brief left and right geo/cel element pairs indexes by inner interfaces gamma */
    TPZStack < std::pair< TPZVec<long>, std::pair< TPZVec<long>, TPZVec<long> > >  > fleft_right_g_c_indexes_gamma;
    
    /** @brief left and right geo/cel element pairs indexes by boundaries interfaces partial Omega */
    TPZStack < std::pair< TPZVec<long>, std::pair< TPZVec<long>, TPZVec<long> > >  > fleft_right_g_c_indexes_Gamma;
    
    /** @brief computational interface element and associated mixed computational element */
    TPZStack < std::pair<long, std::pair< std::pair<long, long> , std::pair<long, long> > >   > fcelint_celh_celp_gamma;
    
    /** @brief computational interface element and associated mixed computational element */
    TPZStack < std::pair<long, std::pair< std::pair<long, long> , std::pair<long, long> > >   > fcelint_celh_celp_Gamma;
    
    /** @brief linear application sw to hyperbolic mesh */
    TRMIrregularBlockDiagonal<STATE> fsw_To_hyperbolic;
    
    /** @brief linear application p_avg to hyperbolic mesh */
    TRMIrregularBlockDiagonal<STATE> fp_avg_To_hyperbolic;
    
    /** @brief linear application grad_u_avg to hyperbolic mesh */
    TRMIrregularBlockDiagonal<STATE> fgrad_u_avg_To_hyperbolic;
    
    /** @brief linear application p_avg to hyperbolic mesh */
    TRMIrregularBlockDiagonal<STATE> fqn_avg_To_hyperbolic_gamma;
    
    /** @brief linear application p_avg to hyperbolic mesh */
    TRMIrregularBlockDiagonal<STATE> fqn_avg_To_hyperbolic_Gamma;
    
    /** @brief linear application sw_avg to parabolic mesh */
    TRMIrregularBlockDiagonal<STATE> fsw_avg_To_parabolic;
    
    /** @brief linear application sw_avg to elliptic mesh */
    TRMIrregularBlockDiagonal<STATE> fsw_avg_To_elliptic;
    

    ////////////////////////// Transfers:: Iterative Coupling by Operator Splitting //////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /** @brief Autopointer of simulation data */
    
    /** @brief Diagonal block matrix to transfer u flux solution to integrations points of the mixed mesh */
    TRMIrregularBlockDiagonal<STATE> fu_To_Mixed;
    
    /** @brief flux dof indexes per element */
//    TPZVec< TPZVec<long> > fu_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer Pressure solution to integrations points of the mixed mesh */
    TRMIrregularBlockDiagonal<STATE> fp_To_Mixed;
    
    /** @brief pressure dof indexes per element */
//    TPZVec< TPZVec<long> > fp_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer saturation solution to integrations points of the transport mesh */
    TRMIrregularBlockDiagonal<STATE> fs_To_Transport;
    
//    /** @brief saturations dof indexes per element */
//    TPZVec< TPZVec<long> > fs_dof_scatter;
    
    /** @brief saturation a dof indexes per element */
    TPZVec< TPZVec<long> > fsa_dof_scatter;
    
    /** @brief staruation b dof indexes per element */
    TPZVec< TPZVec<long> > fsb_dof_scatter;
    
//    /** @brief Diagonal block matrix to transfer Average alpha saturation solution to integrations points of the mixed mesh */
//    TRMIrregularBlockDiagonal<STATE> fs_a_To_Mixed;
    
//    /** @brief pressure dof indexes per element */
//    TPZVec< TPZVec<long> > fs_a_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer Average normal flux solution to integrations points of the transport mesh over gamma */
    TRMIrregularBlockDiagonal<STATE> fun_To_Transport_gamma;
    
    /** @brief Diagonal block matrix to transfer Average normal flux solution to integrations points of the transport mesh over Gamma */
    TRMIrregularBlockDiagonal<STATE> fun_To_Transport_Gamma;
    
    /** @brief Diagonal block matrix to transfer Average normal flux solution to integrations points of the transport mesh over interfaces */
    TRMIrregularBlockDiagonal<STATE> fun_To_Transport;
    
    /** @brief normal flux dof indexes per interface element on gamma (inner interfaces)*/
    TPZVec< TPZVec<long> > fun_dof_scatter_gamma;
    
    /** @brief normal flux dof indexes per interface element on Gamma (boundary interfaces) */
    TPZVec< TPZVec<long> > fun_dof_scatter_Gamma;
    
    /** @brief normal flux dof indexes per interface on inner and boundary interfaces */
    TPZVec< TPZVec<long> > fun_dof_scatter;

    /** @brief mixed and transpor computational multiphysics element indexes, every element is indexed by geometric element */
    TPZStack< std::pair<long, std::pair<long, long> >  > fmixed_transport_cindexes;

    /** @brief mixed and transpor computational multiphysics element indexes, every element is indexed by correponding geometric element index */
    TPZStack< std::pair<long, std::pair<long, std::vector<long> > >  > fmixed_transport_comp_indexes;
    
    /** @brief left and right geometric element indexes on gamma */
    TPZStack < std::pair<long, long> > fleft_right_g_indexes_gamma;
    
    /** @brief geometric interface element indexes on Gamma */
    TPZStack < long > finterface_g_indexes_gamma;
    
    /** @brief left and right geometric element indexes on gamma */
    TPZStack < std::pair<long, long> > fleft_right_g_indexes_Gamma;
    
    /** @brief geometric interface element indexes on Gamma */
    TPZStack < long > finterface_g_indexes_Gamma;

    /** @brief computational interface element and associated mixed computational element */
    TPZStack < std::pair<long, std::pair< std::pair<long, long> , std::pair<long, long> > >   > fcinterface_ctransport_cmixed_indexes_gamma;
    
    /** @brief computational interface element and associated mixed computational element */
    TPZStack < std::pair<long, std::pair< std::pair<long, long> , std::pair<long, long> > >   > fcinterface_ctransport_cmixed_indexes_Gamma;
    

    
    
public:
    
    /** @brief Default constructor */
    TRMBuildTransfers();
    
    /** @brief Default desconstructor */
    ~TRMBuildTransfers();
    
    /** @brief Copy constructor $ */
    TRMBuildTransfers(const TRMBuildTransfers &copy);
    
    /** @brief Copy assignemnt operator $ */
    TRMBuildTransfers &operator=(const TRMBuildTransfers &other);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Transfers:: Iterative Coupling by Operator Splitting //////////////////////////////
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Build methods Elliptic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief bluid linear applications: u and grad_u to elliptic $ */
    void Build_elliptic_To_elliptic(TPZCompMesh * elliptic);
    
    void space_To_elliptic(TPZCompMesh * elliptic);
    
    void spatial_props_To_elliptic(TPZCompMesh * elliptic);
    
    void phi_To_elliptic(TPZCompMesh * elliptic);
    
    void elliptic_To_elliptic(TPZCompMesh * elliptic);
    
    /** @brief bluid linear applications: u and grad_u to parabolic $ */
    void Build_elliptic_To_parabolic(TPZCompMesh * elliptic, TPZCompMesh * parabolic);
    
    void elliptic_To_parabolic(TPZCompMesh * elliptic, TPZCompMesh * parabolic);
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Build methods Parabolic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void Build_parabolic_To_parabolic(TPZCompMesh * parabolic);
    
    void space_To_parabolic(TPZCompMesh * parabolic);
    
    void spatial_props_To_parabolic(TPZCompMesh * parabolic);
    
    void parabolic_To_parabolic(TPZCompMesh * parabolic);
    
    void Build_parabolic_To_elliptic(TPZCompMesh * parabolic, TPZCompMesh * elliptic);
    
    void parabolic_To_elliptic(TPZCompMesh * parabolic, TPZCompMesh * elliptic);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Build methods Hyperbolic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void Build_hyperbolic_To_hyperbolic(TPZCompMesh * hyperbolic);
    
    void spatial_props_To_hyperbolic(TPZCompMesh * hyperbolic);
    
    void hyperbolic_To_hyperbolic(TPZCompMesh * hyperbolic);
    
    void Build_parabolic_hyperbolic_cel_pairs(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic);
    
    void Build_parabolic_hyperbolic_volumetric(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic);
    
    void Build_parabolic_hyperbolic_interfaces(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic, bool bc_interfaceQ);
    
    void Build_parabolic_hyperbolic_left_right_pairs(TPZCompMesh * hyperbolic);
    
    void Build_hyperbolic_parabolic_volumetric(TPZCompMesh * hyperbolic, TPZCompMesh * parabolic);
    
    void parabolic_To_hyperbolic_volumetric(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic);
    
    void parabolic_To_hyperbolic_interfaces(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic, bool bc_interfaceQ);
    
    void hyperbolic_To_parabolic_volumetric(TPZCompMesh * hyperbolic, TPZCompMesh * parabolic);
    
    
    void Build_elliptic_hyperbolic_cel_pairs(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic);
    
    /** @brief bluid linear applications: u and grad_u to hyperbolic $ */
    void Build_elliptic_hyperbolic_volumetric(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic);
    
    void Build_hyperbolic_elliptic_volumetric(TPZCompMesh * hyperbolic, TPZCompMesh * elliptic);
    
    void elliptic_To_hyperbolic(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic);
    
    void hyperbolic_To_elliptic(TPZCompMesh * hyperbolic,TPZCompMesh * elliptic);
    
    void elliptic_To_hyperbolic_volumetric(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic);
    
    void elliptic_To_hyperbolic_volumetricII(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic);
    
    void hyperbolic_To_elliptic_volumetric(TPZCompMesh * hyperbolic, TPZCompMesh * elliptic);
    
    ////////////////////////// Transfers:: Iterative Coupling by Operator Splitting //////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    

    /** @brief Transfer Flux to integration points of multiphysics mesh over volumetric elements */
    void kappa_phi_To_Mixed_Memory(TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer Flux to integration points of multiphysics mesh over volumetric elements */
    void u_To_Mixed_Memory(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer Pressure to integration points of multiphysics mesh over volumetric elements */
    void p_To_Mixed_Memory(TPZCompMesh * cmesh_pressure, TPZCompMesh * cmesh_multiphysics);
    
    
    /** @brief Transfer saturations to integration points of multiphysics transport mesh over volumetric elements */
    void kappa_phi_To_Transport_Memory(TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer saturations to integration points of multiphysics transport mesh over volumetric elements */
    void s_To_Transport_Memory(TPZCompMesh * cmesh_saturation, TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Transfer average pressure to integration points of multiphysics mixed meshes over volumetric elements */
    void p_avg_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed);
    
    /** @brief Transfer average pressure to integration points of multiphysics mixed meshes over volumetric elements */
    void p_avg_Memory_TransferII(TPZCompMesh * cmesh_mf_mixed);
    
    /** @brief Transfer average quantities to integration points of multiphysics mixed/ transpor meshes over volumetric elements */
    void Reciprocal_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_trans);
    
    /** @brief Transfer average quantities to integration points of multiphysics mixed/ transpor meshes over volumetric elements */
    void Reciprocal_Memory_TransferII(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_trans);
    
    /** @brief Transfer normal fluxes to integration points of transport meshes */
    void un_To_Transport_Mesh(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_transport, bool IsBoundaryQ);
    
    /** @brief Transfer normal fluxes to integration points of transport meshes */
    void un_To_Transport_MeshII(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_transport, bool IsBoundaryQ);

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
    
    /** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh over Gamma or gamma */
    void Initialize_un_To_TransportII(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ);
    
    /** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh over Gamma or gamma */
    void Fill_un_To_TransportII(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ);
    
    /** @brief Get the sparse matrix to transfer average normal flux solution to integrations points of the transport mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_un_To_Transport_gamma(){
        return fun_To_Transport_gamma;
    }
    
    
    /** @brief Get the sparse matrix to transfer average normal flux solution to integrations points of the transport mesh  */
    TRMIrregularBlockDiagonal<STATE> Transfer_un_To_Transport_Gamma(){
        return fun_To_Transport_Gamma;
    }
    
    
    

    /** @brief Compute left and right geometric element indexes associated with the transport mesh */
    void ComputeLeftRight(TPZCompMesh * transport_mesh);
    
    /** @brief Compute left and right geometric element indexes associated with the transport mesh */
    void ComputeLeftRightII(TPZCompMesh * transport_mesh);
    

    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Memory operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Get Global integration point indexes associaded  */
    void GlobalPointIndexes(TPZCompEl * cel, TPZManVector<long,30> &int_point_indexes);

    /** @brief Get Global integration point indexes associaded with interfaces */
    void GlobalPointIndexesInterface(TPZCompEl * int_cel, TPZManVector<long,30> &int_point_indexes);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Computational mesh operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Compute geometric mesh pair (mixed, transport) indexed by geometric volumetic element index */
    void FillComputationalElPairs(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport);

    /** @brief Compute computational mesh pair (mixed, transport) indexed by geometric volumetic element index */
    void FillComputationalElPairsII(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport);
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Computational element operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Compute element dof indexes */
    void ElementDofIndexes(TPZInterpolationSpace * &intel,  TPZVec<long> &dof_indexes);
    
    /** @brief Compute element dof indexes */    
    void ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes);

    /** @brief Compute element dof indexes */        
    void ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes, int el_index);
    
    /** @brief Compute element dof indexes at given connect */
    void ElementDofFaceIndexes(int connect,TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes);
    
    /** @brief Compute element dof indexes at given connect */
    void ElementDofFaceIndexes(int connect,TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes);
    
    
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
    
    /** @brief Compute indices associated to faces on 3D topologies */
    void ComputeTransformation(TPZGeoEl * face_gel_origin, TPZGeoEl * gel_origin , TPZGeoEl * gel_target, TPZVec<REAL> & origin, TPZVec<REAL> & target);
    

    
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
//
//  TPZTransferFunctions.h
//  PZ
//
//  Created by Omar on 3/11/16.
//
//

#ifndef TPZTransferFunctions_h
#define TPZTransferFunctions_h

#include <stdio.h>

#include "tpzintpoints.h"
#include "pzmatwithmem.h"
#include "TPZPoroPermMemory.h"
#include "pzmaterial.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"
#include "pzcondensedcompel.h"
#include "TPZBiotPoroelasticity.h"


#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "pzysmp.h"
#include "pzblockdiag.h"
#include "TPZSimulationData.h"
#include "TPZBiotIrregularBlockDiagonal.h"


class TPZTransferFunctions{
    
    
private:
    
    
    /** @brief Autopointer of simulation data */
    TPZSimulationData * fSimulationData;

    /** @brief Galerkin projections from offline process */
    TPZFMatrix<REAL> fgalerkin_projections;
    
    /** @brief Computacional mesh for Galerkin projections from offline process */
    TPZCompMesh * fCmeshRB_projections;
    
    /** @brief Diagonal block matrix to transfer phi_u displacement solution to integrations points of the coupled mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fphi_u_To_Geomechanic;
    
    /** @brief RB dof indexes per element */
    TPZVec< TPZVec<long> > fphi_u_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer grad of phi_u displacement solution to integrations points of the coupled mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgrad_phi_u_To_Geomechanic;
    
    /** @brief Diagonal block matrix to transfer u displacement solution to integrations points of the coupled mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fu_To_Geomechanic;
    
    /** @brief displacement dof indexes per element */
    TPZVec< TPZVec<long> > fu_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer grad of u displacement solution to integrations points of the coupled mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgrad_u_To_Geomechanic;
    

    
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
    
    //    /** @brief left and right geometric element indexes */
    //    TPZStack < std::pair<long, long> > fleft_right_g_indexes;
    //
    //    /** @brief geometric interface element indexes */
    //    TPZStack < long > finterface_g_indexes;
    
    /** @brief computational interface element and associated mixed computational element */
    TPZStack < std::pair<long, std::pair< std::pair<long, long> , std::pair<long, long> > >   > fcinterface_ctransport_cmixed_indexes_gamma;
    
    /** @brief computational interface element and associated mixed computational element */
    TPZStack < std::pair<long, std::pair< std::pair<long, long> , std::pair<long, long> > >   > fcinterface_ctransport_cmixed_indexes_Gamma;
    

    
    
public:
    
    /** @brief Default constructor */
    TPZTransferFunctions();
    
    /** @brief Default desconstructor */
    ~TPZTransferFunctions();
    
    /** @brief Copy constructor $ */
    TPZTransferFunctions(const TPZTransferFunctions &copy);
    
    /** @brief Copy assignemnt operator $ */
    TPZTransferFunctions &operator=(const TPZTransferFunctions &other);
    
    /** @brief Set Galerkin projections from offline process */
    void SetGalerkingProjections(TPZFMatrix<REAL> & galerkin_projections){
        fgalerkin_projections = galerkin_projections;
    }
    
    /** @brief Computacional mesh for Galerkin projections from offline process */
    void SetCmeshGalerkingProjections(TPZCompMesh * CmeshRB_projections){
        fCmeshRB_projections = CmeshRB_projections;
        this->SetGalerkingProjections(fCmeshRB_projections->Solution());
    }
    
    /** @brief Transfer RB basis to integration points of multiphysics mesh over volumetric elements */
    void RB_basis_To_Geomechanic_Memory(TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer u to integration points of multiphysics mesh over volumetric elements */
    void u_To_Geomechanic_Memory(TPZCompMesh * cmesh_elastic, TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer u to integration points of multiphysics mesh over volumetric elements */
    void grad_u_To_Geomechanic_Memory(TPZCompMesh * cmesh_elastic, TPZCompMesh * cmesh_multiphysics);
    
    
    /** @brief Initializate diagonal block matrix to transfer u to multiphysics mesh  */
    void Initialize_RB_basis_To_Geomechanic(TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Initializate diagonal block matrix to transfer u to multiphysics mesh  */
    void Fill_RB_basis_To_Geomechanic(TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TPZBiotIrregularBlockDiagonal<STATE> Transfer_phi_u_To_Geomechanic(){
        return fphi_u_To_Geomechanic;
    }
    
    /** @brief Transfer the RB Solution to multiphysics mesh  */
    void RB_Solution_To_Geomechanic(TPZCompMesh * cmesh_multiphysics, TPZFMatrix<STATE> & rb_solution);
    
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
    void SetSimulationData(TPZSimulationData * SimulationData){
        fSimulationData = SimulationData;
    }
    
    /** @brief Get autopointer of Simulation data */
    TPZSimulationData * SimulationData(){
        return fSimulationData;
    }
    
};


#endif /* TPZTransferFunctions_h */

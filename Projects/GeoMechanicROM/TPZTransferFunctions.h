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
#include "pzmaterial.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"
#include "pzcondensedcompel.h"
#include "TPZBiotPoroelasticity.h"

// Memory
#include "TPZPoroPermMemory.h"
#include "TPZElasticBiotMemory.h"
#include "TPZDarcyFlowMemory.h"

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
    
    /** @brief Diagonal block matrix to transfer grad of u displacement solution to integrations points of the coupled mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgrad_u_To_Geomechanic;

    
    /** @brief geomechanic and galerkin projection mesh computational multiphysics element indexes, every element is indexed by geometric element */
    TPZStack< std::pair<long, std::pair<long, long> >  > fgeomechanic_galerkinp_cindexes;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Attributes Elliptic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief full order dof indexes per element */
    TPZVec< TPZVec<long> > fu_dof_scatter;

    /** @brief integration point indexes */
    TPZStack< std::pair<long, std::pair< TPZVec<long>, TPZVec<long> > >  > fe_p_intp_indexes;
    
    TPZStack< std::pair<long, std::pair<long, long> >  > fe_p_cindexes;
    
    /** @brief linear application u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fu_To_elliptic;
    
    /** @brief linear application grad_u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgrad_u_To_elliptic;
    
    /** @brief linear application u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fu_To_parabolic;
    
    /** @brief linear application grad_u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgrad_u_To_parabolic;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Attributes Elliptic RB
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief full order dof indexes per element */
    TPZVec< TPZVec<long> > fgp_u_dof_scatter;
    
    TPZStack< std::pair<long, std::pair<long, long> >  > fgp_rb_e_cindexes;
    
    TPZStack< std::pair<long, std::pair<long, long> >  > fgp_rb_p_cindexes;
    
    /** @brief integration point indexes */
    TPZStack< std::pair<long, std::pair< TPZVec<long>, TPZVec<long> > >  > frb_e_p_intp_indexes;
    
    TPZStack< std::pair<long, std::pair<long, long> >  > frb_e_p_cindexes;
    
    /** @brief linear application u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgp_u_To_rb_elliptic;
    
    /** @brief linear application grad_u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgp_grad_u_To_rb_elliptic;
    
    /** @brief linear application u to parabolic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgp_u_To_parabolic;
    
    /** @brief linear application grad_u to parabolic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgp_grad_u_To_parabolic;
    
    
    /** @brief linear application rb u to elliptic mesh */
    TPZFMatrix<STATE> frb_u_To_rb_elliptic;
    
    /** @brief linear application rb grad_u to elliptic mesh */
    TPZFMatrix<STATE> frb_grad_u_To_rb_elliptic;
    
    
    /** @brief linear application u to parabolic mesh */
    TPZFMatrix<STATE> frb_u_To_parabolic;
    
    /** @brief linear application grad_u to parabolic mesh */
    TPZFMatrix<STATE> frb_grad_u_To_parabolic;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) :: Attributes Parabolic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief  p dof indexes per element */
    TPZVec< TPZVec<long> > fp_dof_scatter;
    
    /** @brief  flux dof indexes per element */
    TPZVec< TPZVec<long> > fq_dof_scatter;
    
    /** @brief integration point indexes */
    TPZStack< std::pair<long, std::pair< TPZVec<long>, TPZVec<long> > >  > fp_e_intp_indexes;
    
    TPZStack< std::pair<long, std::pair<long, long> >  > fp_e_cindexes;
    
    /** @brief linear application u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fp_To_elliptic;
    
    /** @brief linear application grad_u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgrad_p_To_elliptic;
    
    /** @brief linear application u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fp_To_parabolic;
    
    /** @brief linear application grad_u to elliptic mesh */
    TPZBiotIrregularBlockDiagonal<STATE> fgrad_p_To_parabolic;
    
    
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
    // Segregated Transfer methods (Gamma and Omega) Elliptic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void Fill_elliptic_To_elliptic(TPZCompMesh * elliptic);
    
    void Fill_elliptic_To_parabolic(TPZCompMesh * elliptic, TPZCompMesh * parabolic);
    
    void space_To_elliptic(TPZCompMesh * elliptic);
    
    void elliptic_To_elliptic(TPZCompMesh * elliptic);
    
    void elliptic_To_parabolic(TPZCompMesh * elliptic, TPZCompMesh * parabolic);
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) Elliptic RB
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void Fill_gp_elliptic_To_rb_elliptic(TPZCompMesh * gp_mesh, TPZCompMesh * elliptic);
    
    void Fill_gp_elliptic_To_parabolic(TPZCompMesh * gp_mesh, TPZCompMesh * parabolic);
    
    void rb_space_To_rb_elliptic(TPZCompMesh * elliptic);
    
    void rb_elliptic_To_rb_elliptic(TPZCompMesh * elliptic);
    
    void rb_elliptic_To_parabolic(TPZCompMesh * elliptic, TPZCompMesh * parabolic);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Transfer methods (Gamma and Omega) Parabolic
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void Fill_parabolic_To_parabolic(TPZCompMesh * parabolic);
    
    void Fill_M_parabolic_To_parabolic(TPZCompMesh * parabolic);
    
    void Fill_parabolic_To_elliptic(TPZCompMesh * parabolic, TPZCompMesh * elliptic);

    void space_To_parabolic(TPZCompMesh * parabolic);
    
    void space_M_To_parabolic(TPZCompMesh * parabolic);
    
    void parabolic_To_parabolic(TPZCompMesh * parabolic);
    
    void M_parabolic_To_parabolic(TPZCompMesh * parabolic);
    
    void parabolic_To_elliptic(TPZCompMesh * parabolic, TPZCompMesh * elliptic);
    
    
    
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
    
    /** @brief Compute geometric mesh pair (geomechanic, gp cmesh) indexed by geometric volumetic element index */
    void FillGeomechanicElPairs(TPZCompMesh * cmesh_mphysics);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Computational element operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Compute element dof indexes */
    void ElementDofIndexes(TPZInterpolationSpace * &intel,  TPZVec<long> &dof_indexes);
    
    /** @brief Compute element dof indexes */
    void ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes);
    
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
    void SetSimulationData(TPZSimulationData * SimulationData){
        fSimulationData = SimulationData;
    }
    
    /** @brief Get autopointer of Simulation data */
    TPZSimulationData * SimulationData(){
        return fSimulationData;
    }
    
};


#endif /* TPZTransferFunctions_h */

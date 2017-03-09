//
//  TPZTransfer.cpp
//  PZ
//
//  Created by Omar on 3/1/17.
//
//

#include "TPZTransferFunctions.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


/** @brief Default constructor */
TPZTransferFunctions::TPZTransferFunctions(){
    
    fSimulationData = NULL;
    
}

/** @brief Default desconstructor */
TPZTransferFunctions::~TPZTransferFunctions(){
    
}

/** @brief Copy constructor $ */
TPZTransferFunctions::TPZTransferFunctions(const TPZTransferFunctions &copy)
{
    fSimulationData = copy.fSimulationData;
}

/** @brief Copy assignemnt operator $ */
TPZTransferFunctions & TPZTransferFunctions::operator=(const TPZTransferFunctions &other)
{
    if (this != & other) // prevent self-assignment
    {
        fSimulationData = other.fSimulationData;
    }
    return *this;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Matrices Initialization Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TPZTransferFunctions::Initialize_RB_basis_To_Geomechanic(TPZCompMesh * cmesh_multiphysics){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics || !fCmeshRB_projections || fgeomechanic_galerkinp_cindexes.size() == 0) {
        std::cout << "There is no computational meshes, cmesh_multiphysics||fCmeshRB_projections = Null." << std::endl;
        DebugStop();
    }
#endif
    
    
    fCmeshRB_projections->LoadReferences();
    
    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    long nel = fgeomechanic_galerkinp_cindexes.size();
    int n_var_dim = fCmeshRB_projections->Reference()->Dimension(); // vectorial
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fphi_u_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions_phi(nel);
    TPZVec< std::pair<long, long> > blocks_dimensions_grad_phi(nel);
    
    long gel_index, geomech_index, gp_index;
    for (long iel = 0; iel < nel; iel++) {
        
        gel_index = fgeomechanic_galerkinp_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        geomech_index   = fgeomechanic_galerkinp_cindexes[iel].second.first;
        gp_index        = fgeomechanic_galerkinp_cindexes[iel].second.second;

        TPZCompEl * geomec_cel = cmesh_multiphysics->Element(geomech_index);
        TPZCompEl * gp_cel = fCmeshRB_projections->Element(gp_index);
        
#ifdef PZDEBUG
        if (!geomec_cel || !gp_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(geomec_cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = gp_cel->Index();
        
        TPZInterpolationSpace * gp_intel = dynamic_cast<TPZInterpolationSpace * >(gp_cel);
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        int gel_dim = gel->Dimension();
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(gp_intel, dof_indexes);
        fphi_u_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions_phi[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions_phi[element_index].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[element_index].first = int_point_indexes.size()*n_var_dim*gel_dim;
        blocks_dimensions_grad_phi[element_index].second = dof_indexes.size();
    }
    
    // Initialize the matrix
    fphi_u_To_Geomechanic.Initialize(blocks_dimensions_phi);
    fgrad_phi_u_To_Geomechanic.Initialize(blocks_dimensions_grad_phi);
    
}


void TPZTransferFunctions::Fill_RB_basis_To_Geomechanic(TPZCompMesh * cmesh_multiphysics){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_RB_basis_To_Geomechanic(cmesh_multiphysics);
    
    
    fCmeshRB_projections->LoadReferences();
    
    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    long nel = fgeomechanic_galerkinp_cindexes.size();
    int n_var_dim = fCmeshRB_projections->Reference()->Dimension(); // vectorial
    long element_index = 0;
    
    std::pair<long, long> block_dim;
    
    long gel_index, geomech_index, gp_index;
    for (long iel = 0; iel < nel; iel++) {
        
        gel_index = fgeomechanic_galerkinp_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        geomech_index   = fgeomechanic_galerkinp_cindexes[iel].second.first;
        gp_index        = fgeomechanic_galerkinp_cindexes[iel].second.second;
        
        TPZCompEl * geomec_cel = cmesh_multiphysics->Element(geomech_index);
        TPZCompEl * gp_cel = fCmeshRB_projections->Element(gp_index);
        
#ifdef PZDEBUG
        if (!geomec_cel || !gp_cel) {
            DebugStop();
        }
#endif
        

        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(geomec_cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = gp_cel->Index();
        TPZInterpolationSpace * gp_intel = dynamic_cast<TPZInterpolationSpace * >(gp_cel);
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        dof_indexes = fphi_u_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_geomechanic = mf_cel->GetIntegrationRule();
        int np_cel = int_points_geomechanic.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi, block_grad_phi;
        
        block_phi.Resize(block_dim.first*n_var_dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*n_var_dim*gel_dim,block_dim.second);

        // for derivatives in real space
        int nshape = gp_intel->NShapeF();
        TPZFNMatrix<220> phi(nshape,1);
        TPZFNMatrix<660> dphi(gel_dim,nshape),dphix_axes(gel_dim,nshape);
        TPZFMatrix<double> dphidx;
        TPZFNMatrix<9,STATE> jacobian(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> jacinv(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> axes;
        REAL detjac;
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(gel_dim,0.0);
            STATE w;
            int_points_geomechanic.Point(ip, qsi, w);
            
            // Get the phi and dphix for H1 elasticity
            gp_intel->Shape(qsi, phi, dphi);
            gel->Jacobian( qsi, jacobian, axes, detjac , jacinv);
            
            switch(gel_dim) {
                case 0:
                    break;
                case 1:
                    dphix_axes = dphi;
                    dphix_axes *= (1./detjac);
                    break;
                case 2:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
                    }
                    break;
                case 3:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
                        dphix_axes(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
                    }
                    break;
                default:
                    std::stringstream sout;
                    sout << "pzintel.c please implement the " << gel_dim << "d Jacobian and inverse\n";
                    LOGPZ_ERROR(logger,sout.str());
            }
            
            TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            
#ifdef PZDEBUG
            if(block_dim.second != phi.Rows() * n_var_dim){
                DebugStop();
            }
#endif
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    block_phi(ip*n_var_dim+id,jp*n_var_dim+id) = phi(jp,0);
                }
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    for (int jd = 0; jd < gel_dim; jd++) {
                        if(gel_dim == n_var_dim){
                            block_grad_phi(ip*n_var_dim*gel_dim + id*gel_dim + jd,jp*n_var_dim+id) = dphidx(jd,jp);
                        }
                        else{
                            block_grad_phi(ip*n_var_dim*gel_dim+id*gel_dim + jd,jp*n_var_dim+id) = 0.0;
                        }
                    }
                }
            }
            
        }
        
        fphi_u_To_Geomechanic.SetBlock(element_index, block_phi);
        fgrad_phi_u_To_Geomechanic.SetBlock(element_index, block_grad_phi);
    }
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segregated Transfer methods (Gamma and Omega) Elliptic
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void TPZTransferFunctions::Fill_elliptic_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
#endif
    
    fe_p_cindexes.resize(0);
    
    
    elliptic->LoadReferences();
    TPZGeoMesh * geometry = elliptic->Reference();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<long, std::pair <long,long> > chunk_indexes;
    long valid_elements = 0;
    // inserting geometric indexes
    for (long i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        valid_elements++;
    }
    
    
    long nel = geometry->NElements();
    int n_var_dim = geometry->Dimension();
    
    // Compute destination index scatter by element (Omega and Gamma)
    
    fu_dof_scatter.resize(valid_elements);
    std::pair<long, std::pair <TPZVec<long>, TPZVec<long> > > chunk_intp_indexes;
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions_phi(valid_elements);
    TPZVec< std::pair<long, long> > blocks_dimensions_grad_phi(valid_elements);
    
    long element_index = 0;
    for (long iel = 0; iel < nel; iel++) {
        
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
    
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        // Getting local integration index
        TPZManVector<long> e_int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        int gel_dim = gel->Dimension();
        
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        
        this->ElementDofIndexes(mf_e_cel, dof_indexes);
        fu_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions_phi[element_index].first = e_int_point_indexes.size()*n_var_dim;
        blocks_dimensions_phi[element_index].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[element_index].first = e_int_point_indexes.size()*n_var_dim*gel_dim;
        blocks_dimensions_grad_phi[element_index].second = dof_indexes.size();
        
        element_index++;
    }
    
    // Initialize the matrix
    fu_To_elliptic.Initialize(blocks_dimensions_phi);
    fgrad_u_To_elliptic.Initialize(blocks_dimensions_grad_phi);
    
    
    element_index = 0;
    
    TPZManVector<long> e_int_point_indexes(0,0);
    TPZManVector<long> dof_indexes(0,0);
    
    std::pair<long, long> block_dim;
    for (long iel = 0; iel < nel; iel++) {
        
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        
        
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);

        
#ifdef PZDEBUG
        if(!mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        
        TPZInterpolationSpace * e_intel = dynamic_cast<TPZInterpolationSpace * >(mf_e_cel->Element(0));
        
        // Getting local integration index
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        dof_indexes = fu_dof_scatter[element_index];
        
        block_dim.first = e_int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_geomechanic = mf_e_cel->GetIntegrationRule();
        int np_cel = int_points_geomechanic.NPoints();
        
#ifdef PZDEBUG
        if (e_int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi, block_grad_phi;
        
        block_phi.Resize(block_dim.first*n_var_dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*n_var_dim*gel_dim,block_dim.second);
        
        // for derivatives in real space
        int nshape = e_intel->NShapeF();
        TPZFNMatrix<220> phi(nshape,1);
        TPZFNMatrix<660> dphi(gel_dim,nshape),dphix_axes(gel_dim,nshape);
        TPZFMatrix<double> dphidx;
        TPZFNMatrix<9,STATE> jacobian(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> jacinv(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> axes;
        REAL detjac;
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(gel_dim,0.0);
            STATE w;
            int_points_geomechanic.Point(ip, qsi, w);
            
            // Get the phi and dphix for H1 elasticity
            e_intel->Shape(qsi, phi, dphi);
            gel->Jacobian( qsi, jacobian, axes, detjac , jacinv);
            
            switch(gel_dim) {
                case 0:
                    break;
                case 1:
                    dphix_axes = dphi;
                    dphix_axes *= (1./detjac);
                    break;
                case 2:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
                    }
                    break;
                case 3:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
                        dphix_axes(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
                    }
                    break;
                default:
                    std::stringstream sout;
                    sout << "pzintel.c please implement the " << gel_dim << "d Jacobian and inverse\n";
                    LOGPZ_ERROR(logger,sout.str());
            }
            
            TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            
#ifdef PZDEBUG
            if(block_dim.second != phi.Rows() * n_var_dim){
                DebugStop();
            }
#endif
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    block_phi(ip*n_var_dim+id,jp*n_var_dim+id) = phi(jp,0);
                }
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    for (int jd = 0; jd < gel_dim; jd++) {
                        if(gel_dim == n_var_dim){
                            block_grad_phi(ip*n_var_dim*gel_dim + id*gel_dim + jd,jp*n_var_dim+id) = dphidx(jd,jp);
                        }
                        else{
                            block_grad_phi(ip*n_var_dim*gel_dim+id*gel_dim + jd,jp*n_var_dim+id) = 0.0;
                        }
                    }
                }
            }
            
        }
        
        fu_To_elliptic.SetBlock(element_index, block_phi);
        fgrad_u_To_elliptic.SetBlock(element_index, block_grad_phi);
        
        element_index++;
    }
    
//    fu_To_elliptic.Print(" u_to_e ");
//    fgrad_u_To_elliptic.Print(" grad_u_to_e ");
    
    return;
    
}

void TPZTransferFunctions::Fill_elliptic_To_parabolic(TPZCompMesh * origin, TPZCompMesh * target){
    
#ifdef PZDEBUG
    if (!origin || !target) {
        DebugStop();
    }
#endif
    
    fe_p_cindexes.resize(0);
    
    
    origin->LoadReferences();
    TPZGeoMesh * geometry = origin->Reference();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<long, std::pair <long,long> > chunk_indexes;

    // inserting geometric indexes
    for (long i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        chunk_indexes.first = gel->Index();
        chunk_indexes.second.first = -1;
        chunk_indexes.second.second = -1;
        fe_p_cindexes.Push(chunk_indexes);
        
    }
    
    // counting volumetric elements
    long n_elements = fe_p_cindexes.size();
    
    // inserting origin indexes
    origin->LoadReferences();
    for(long iel = 0; iel < n_elements; iel++){
        
        TPZCompEl * cel = geometry->Element(fe_p_cindexes[iel].first)->Reference();
        
#ifdef PZDEBUG
        if (!cel) {
            //continue;
            DebugStop();
        }
#endif
        fe_p_cindexes[iel].second.first = cel->Index();
        
    }
    
    // inserting target indexes
    target->LoadReferences();
    for(long iel = 0; iel < n_elements; iel++){
        
        TPZCompEl * cel = geometry->Element(fe_p_cindexes[iel].first)->Reference();
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        fe_p_cindexes[iel].second.second = cel->Index();
        
    }
    
//    for (int k = 0; k < comp_indexes.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << fe_p_cindexes[k].first <<std::endl;
//        std::cout << " volume c_o : " << fe_p_cindexes[k].second.first <<std::endl;
//        std::cout << " volume c_t : " << fe_p_cindexes[k].second.second <<std::endl;
//    }
    
    
    origin->LoadReferences();
    geometry = origin->Reference();
    
    long nel = fe_p_cindexes.size();
    int n_var_dim = origin->Reference()->Dimension();
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    
    fu_dof_scatter.resize(nel);
    std::pair<long, std::pair <TPZVec<long>, TPZVec<long> > > chunk_intp_indexes;
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions_phi(nel);
    TPZVec< std::pair<long, long> > blocks_dimensions_grad_phi(nel);
    
    long gel_index, o_index, t_index;
    for (long iel = 0; iel < nel; iel++) {
        
        gel_index = fe_p_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        o_index   = fe_p_cindexes[iel].second.first;
        t_index   = fe_p_cindexes[iel].second.second;
        
        TPZCompEl * o_cel = origin->Element(o_index);
        TPZCompEl * t_cel = target->Element(t_index);
        
#ifdef PZDEBUG
        if (!o_cel || !t_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_o_cel = dynamic_cast<TPZMultiphysicsElement * >(o_cel);
        TPZMultiphysicsElement * mf_t_cel = dynamic_cast<TPZMultiphysicsElement * >(t_cel);
        
#ifdef PZDEBUG
        if(!mf_o_cel || !mf_t_cel)
        {
            DebugStop();
        }
#endif
        
        // Getting local integration index
        TPZManVector<long> o_int_point_indexes(0,0);
        TPZManVector<long> t_int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        int gel_dim = gel->Dimension();
        
        mf_o_cel->GetMemoryIndices(o_int_point_indexes);
        mf_t_cel->GetMemoryIndices(t_int_point_indexes);
        
        chunk_intp_indexes.first = gel->Index();
        chunk_intp_indexes.second.first  = o_int_point_indexes;
        chunk_intp_indexes.second.second = t_int_point_indexes;
        fe_p_intp_indexes.Push(chunk_intp_indexes);
        
        this->ElementDofIndexes(mf_o_cel, dof_indexes);
        fu_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions_phi[element_index].first = t_int_point_indexes.size()*n_var_dim;
        blocks_dimensions_phi[element_index].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[element_index].first = t_int_point_indexes.size()*n_var_dim*gel_dim;
        blocks_dimensions_grad_phi[element_index].second = dof_indexes.size();
        
        element_index++;
    }
    
    // Initialize the matrix
    fu_To_parabolic.Initialize(blocks_dimensions_phi);
    fgrad_u_To_parabolic.Initialize(blocks_dimensions_grad_phi);
    
    
    origin->LoadReferences();
    element_index = 0;
    
    TPZManVector<long> int_point_indexes(0,0);
    TPZManVector<long> dof_indexes(0,0);
    
    std::pair<long, long> block_dim;
    for (long iel = 0; iel < nel; iel++) {
        
        gel_index = fe_p_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        o_index   = fe_p_cindexes[iel].second.first;
        t_index   = fe_p_cindexes[iel].second.second;
        
        TPZCompEl * o_cel = origin->Element(o_index);
        TPZCompEl * t_cel = target->Element(t_index);
        
#ifdef PZDEBUG
        if (!o_cel || !t_cel) {
            DebugStop();
        }
#endif
        
        
        TPZMultiphysicsElement * mf_o_cel = dynamic_cast<TPZMultiphysicsElement * >(o_cel);
        TPZMultiphysicsElement * mf_t_cel = dynamic_cast<TPZMultiphysicsElement * >(t_cel);
        
#ifdef PZDEBUG
        if(!mf_o_cel || !mf_t_cel)
        {
            DebugStop();
        }
#endif
      
        
        TPZInterpolationSpace * o_intel = dynamic_cast<TPZInterpolationSpace * >(mf_o_cel->Element(0));

        // Getting local integration index
        int_point_indexes = fe_p_intp_indexes[element_index].second.second;
        dof_indexes = fu_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_geomechanic = mf_t_cel->GetIntegrationRule();
        int np_cel = int_points_geomechanic.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi, block_grad_phi;
        
        block_phi.Resize(block_dim.first*n_var_dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*n_var_dim*gel_dim,block_dim.second);
        
        // for derivatives in real space
        int nshape = o_intel->NShapeF();
        TPZFNMatrix<220> phi(nshape,1);
        TPZFNMatrix<660> dphi(gel_dim,nshape),dphix_axes(gel_dim,nshape);
        TPZFMatrix<double> dphidx;
        TPZFNMatrix<9,STATE> jacobian(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> jacinv(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> axes;
        REAL detjac;
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(gel_dim,0.0);
            STATE w;
            int_points_geomechanic.Point(ip, qsi, w);
            
            // Get the phi and dphix for H1 elasticity
            o_intel->Shape(qsi, phi, dphi);
            gel->Jacobian( qsi, jacobian, axes, detjac , jacinv);
            
            switch(gel_dim) {
                case 0:
                    break;
                case 1:
                    dphix_axes = dphi;
                    dphix_axes *= (1./detjac);
                    break;
                case 2:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
                    }
                    break;
                case 3:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
                        dphix_axes(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
                    }
                    break;
                default:
                    std::stringstream sout;
                    sout << "pzintel.c please implement the " << gel_dim << "d Jacobian and inverse\n";
                    LOGPZ_ERROR(logger,sout.str());
            }
            
            TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            
#ifdef PZDEBUG
            if(block_dim.second != phi.Rows() * n_var_dim){
                DebugStop();
            }
#endif
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    block_phi(ip*n_var_dim+id,jp*n_var_dim+id) = phi(jp,0);
                }
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    for (int jd = 0; jd < gel_dim; jd++) {
                        if(gel_dim == n_var_dim){
                            block_grad_phi(ip*n_var_dim*gel_dim + id*gel_dim + jd,jp*n_var_dim+id) = dphidx(jd,jp);
                        }
                        else{
                            block_grad_phi(ip*n_var_dim*gel_dim+id*gel_dim + jd,jp*n_var_dim+id) = 0.0;
                        }
                    }
                }
            }
            
        }
        
        fu_To_parabolic.SetBlock(element_index, block_phi);
        fgrad_u_To_parabolic.SetBlock(element_index, block_grad_phi);
        
        element_index++;
    }
    
//    fu_To_parabolic.Print(" u_to_p ");
//    fgrad_u_To_parabolic.Print(" grad_u_to_p ");

    return;
    
}

void TPZTransferFunctions::space_To_elliptic(TPZCompMesh * elliptic){
    
    
#ifdef PZDEBUG
    if (!elliptic || fu_To_elliptic.Rows() == 0 || fgrad_u_To_elliptic.Rows() == 0) {
        DebugStop();
    }
#endif
    
    elliptic->LoadReferences();
    TPZGeoMesh * geometry = elliptic->Reference();
    long nel = geometry->NElements();
    int dim = elliptic->Dimension();

    
    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    TPZFMatrix<double> block_phi;
    TPZFMatrix<double> block_grad_phi;
    
    for (int iel = 0; iel < nel; iel++) {
        
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fu_To_elliptic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_u_To_elliptic.GetSizeofBlock(iblock);
        fu_To_elliptic.GetBlock(iblock, block_phi);
        fgrad_u_To_elliptic.GetBlock(iblock, block_grad_phi);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TPZElasticBiotMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZElasticBiotMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            int gel_dim = gel->Dimension();
            int nshapes = b_size_phi.second/dim;
            TPZFNMatrix<3,STATE> phi_u(nshapes,1,0.0);
            TPZFNMatrix<9,STATE> grad_phi_u(dim,nshapes);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < nshapes; is++) {
                    phi_u(is,0) = block_phi(ip*dim,is*dim);
                    
                }
                
                for (int is = 0; is < nshapes; is++) {
                    for (int id = 0 ; id < dim; id++) {
                        grad_phi_u(id,is) = block_grad_phi(ip*dim*gel_dim + id,is*dim);
                    }
                    
                }
                
                associated_material->GetMemory()[ipos].Set_phi_u(phi_u);
                associated_material->GetMemory()[ipos].Set_grad_phi_u(grad_phi_u);
            }
            
            
        }
        else{
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            int gel_dim = gel->Dimension();
            int nshapes = b_size_phi.second/dim;
            TPZFNMatrix<3,STATE> phi_u(nshapes,1,0.0);
            TPZFNMatrix<9,STATE> grad_phi_u(dim,nshapes);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < nshapes; is++) {
                    phi_u(is,0) = block_phi(ip*dim,is*dim);
                    
                }
                
                for (int is = 0; is < nshapes; is++) {
                    for (int id = 0 ; id < dim; id++) {
                            grad_phi_u(id,is) = block_grad_phi(ip*dim*gel_dim + id,is*dim);
                    }
                    
                }
                
                associated_material->GetMemory()[ipos].Set_phi_u(phi_u);
                associated_material->GetMemory()[ipos].Set_grad_phi_u(grad_phi_u);
            }
            
        }
        
        iblock++;
        
    }
    
    
}


void TPZTransferFunctions::elliptic_To_elliptic(TPZCompMesh * elliptic){
    
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
#endif
    
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_u(fu_To_parabolic.Cols(),1,0.0);
    int n = fe_p_cindexes.size();
    long pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fu_dof_scatter[i].size(); iequ++) {
            Scatter_u(pos,0) = elliptic->Solution()(fu_dof_scatter[i][iequ],0);
            pos++;
        }
    }

    // Step two
    TPZFMatrix<STATE> u_at_elliptic,grad_u_at_elliptic;
    fu_To_elliptic.Multiply(Scatter_u,u_at_elliptic);
    fgrad_u_To_elliptic.Multiply(Scatter_u, grad_u_at_elliptic);

    elliptic->LoadReferences();
    TPZGeoMesh * geometry = elliptic->Reference();
    long nel = geometry->NElements();
    int dim = elliptic->Dimension();

    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    for (int iel = 0; iel < nel; iel++) {
    
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fu_To_elliptic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_u_To_elliptic.GetSizeofBlock(iblock);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TPZElasticBiotMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZElasticBiotMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,dim,0.0);
            TPZFNMatrix<9,STATE> grad_u(dim,dim);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_elliptic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        grad_u(jd,id)= grad_u_at_elliptic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                    }
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_u_n(u);
                    associated_material->GetMemory()[ipos].Set_grad_u_n(grad_u);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_u(u);
                    associated_material->GetMemory()[ipos].Set_grad_u(grad_u);
                }
                
            }
            
            
        }
        else{
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,dim,0.0);
            TPZFNMatrix<9,STATE> grad_u(dim,dim);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_elliptic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        grad_u(jd,id)= grad_u_at_elliptic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                    }
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_u_n(u);
                    associated_material->GetMemory()[ipos].Set_grad_u_n(grad_u);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_u(u);
                    associated_material->GetMemory()[ipos].Set_grad_u(grad_u);
                }
                
            }
            
        }
        
        iblock++;
        
    }
    
    
}

void TPZTransferFunctions::elliptic_To_parabolic(TPZCompMesh * elliptic, TPZCompMesh * parabolic){
    
    
#ifdef PZDEBUG
    if (!elliptic || !parabolic) {
        DebugStop();
    }
#endif
    
    
    TPZGeoMesh * geometry = elliptic->Reference();
    long nel = fe_p_cindexes.size();
    int dim = elliptic->Dimension();
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_u(fu_To_parabolic.Cols(),1,0.0);
    
    long pos = 0;
    for (int el = 0; el < nel; el++) {
        for(int iequ = 0; iequ < fu_dof_scatter[el].size(); iequ++) {
            Scatter_u(pos,0) = elliptic->Solution()(fu_dof_scatter[el][iequ],0);
            pos++;
        }
    }
    
    
    
    // Step two
    TPZFMatrix<STATE> u_at_parabolic,grad_u_at_parabolic;
    fu_To_parabolic.Multiply(Scatter_u,u_at_parabolic);
    fgrad_u_To_parabolic.Multiply(Scatter_u, grad_u_at_parabolic);
    
    
    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    long gel_index, e_index, p_index;
    for (int iel = 0; iel < nel; iel++) {
        
        
        gel_index = fe_p_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        e_index   = fe_p_cindexes[iel].second.first;
        p_index   = fe_p_cindexes[iel].second.second;
        
        TPZCompEl * e_cel = elliptic->Element(e_index);
        TPZCompEl * p_cel = parabolic->Element(p_index);
        
#ifdef PZDEBUG
        if (!e_cel || !p_cel) {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fu_To_parabolic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_u_To_parabolic.GetSizeofBlock(iblock);

        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = parabolic->FindMaterial(matd_id);
            TPZMatWithMem<TPZDarcyFlowMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZDarcyFlowMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            int_point_indexes = fe_p_intp_indexes[iblock].second.second;
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,dim,0.0);
            TPZFNMatrix<9,STATE> grad_u(dim,dim);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_parabolic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        grad_u(id,jd)= grad_u_at_parabolic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                    }
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_u_n(u);
                    associated_material->GetMemory()[ipos].Set_grad_u_n(grad_u);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_u(u);
                    associated_material->GetMemory()[ipos].Set_grad_u(grad_u);
                }

            }
            
            
        }
        else{
            
            TPZMaterial * material = parabolic->FindMaterial(matd_id);
            TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            int_point_indexes = fe_p_intp_indexes[iblock].second.second;
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,dim,0.0);
            TPZFNMatrix<9,STATE> grad_u(dim,dim);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_parabolic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        grad_u(id,jd)= grad_u_at_parabolic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                    }
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_u_n(u);
                    associated_material->GetMemory()[ipos].Set_grad_u_n(grad_u);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_u(u);
                    associated_material->GetMemory()[ipos].Set_grad_u(grad_u);
                }
                
            }
            
        }
        
        iblock++;
        
    }
    
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segregated Transfer methods (Gamma and Omega) Parabolic
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void TPZTransferFunctions::Fill_parabolic_To_parabolic(TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    fp_e_cindexes.resize(0);
    
    
    parabolic->LoadReferences();
    TPZGeoMesh * geometry = parabolic->Reference();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<long, std::pair <long,long> > chunk_indexes;
    long valid_elements = 0;
    // inserting geometric indexes
    for (long i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        valid_elements++;
    }
    
    
    long nel = geometry->NElements();
    int n_var_dim = geometry->Dimension();
    
    // Compute destination index scatter by element (Omega and Gamma)
    
    fp_dof_scatter.resize(valid_elements);
    std::pair<long, std::pair <TPZVec<long>, TPZVec<long> > > chunk_intp_indexes;
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions_phi(valid_elements);
    TPZVec< std::pair<long, long> > blocks_dimensions_grad_phi(valid_elements);
    
    long element_index = 0;
    for (long iel = 0; iel < nel; iel++) {
        
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        TPZCompEl * p_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        // Getting local integration index
        TPZManVector<long> p_int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);
        
        this->ElementDofIndexes(mf_p_cel, dof_indexes);
        fp_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions_phi[element_index].first = p_int_point_indexes.size();
        blocks_dimensions_phi[element_index].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[element_index].first = p_int_point_indexes.size()*n_var_dim;
        blocks_dimensions_grad_phi[element_index].second = dof_indexes.size();
        
        element_index++;
    }
    
    // Initialize the matrix
    fp_To_parabolic.Initialize(blocks_dimensions_phi);
    fgrad_p_To_parabolic.Initialize(blocks_dimensions_grad_phi);
    
    
    element_index = 0;
    
    TPZManVector<long> p_int_point_indexes(0,0);
    TPZManVector<long> dof_indexes(0,0);
    
    std::pair<long, long> block_dim;
    for (long iel = 0; iel < nel; iel++) {
        
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        
        
        TPZCompEl * p_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
        
#ifdef PZDEBUG
        if(!mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        
        TPZInterpolationSpace * p_intel = dynamic_cast<TPZInterpolationSpace * >(mf_p_cel->Element(0));
        
        // Getting local integration index
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);
        dof_indexes = fp_dof_scatter[element_index];
        
        block_dim.first = p_int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points = mf_p_cel->GetIntegrationRule();
        int np_cel = int_points.NPoints();
        
#ifdef PZDEBUG
        if (p_int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi, block_grad_phi;
        
        block_phi.Resize(block_dim.first*n_var_dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*n_var_dim*gel_dim,block_dim.second);
        
        // for derivatives in real space
        int nshape = p_intel->NShapeF();
        TPZFNMatrix<220> phi(nshape,1);
        TPZFNMatrix<660> dphi(gel_dim,nshape),dphix_axes(gel_dim,nshape);
        TPZFMatrix<double> dphidx;
        TPZFNMatrix<9,STATE> jacobian(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> jacinv(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> axes;
        REAL detjac;
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(gel_dim,0.0);
            STATE w;
            int_points.Point(ip, qsi, w);
            
            // Get the phi and dphix for H1 elasticity
            p_intel->Shape(qsi, phi, dphi);
            gel->Jacobian( qsi, jacobian, axes, detjac , jacinv);
            
            switch(gel_dim) {
                case 0:
                    break;
                case 1:
                    dphix_axes = dphi;
                    dphix_axes *= (1./detjac);
                    break;
                case 2:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
                    }
                    break;
                case 3:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
                        dphix_axes(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
                    }
                    break;
                default:
                    std::stringstream sout;
                    sout << "pzintel.c please implement the " << gel_dim << "d Jacobian and inverse\n";
                    LOGPZ_ERROR(logger,sout.str());
            }
            
            TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            
#ifdef PZDEBUG
            if(block_dim.second != phi.Rows()){
                DebugStop();
            }
#endif
            for (int jp = 0; jp < phi.Rows(); jp++) {
                block_phi(ip,jp) = phi(jp,0);
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    block_grad_phi(ip*n_var_dim + id,jp) = dphidx(id,jp);
                }
            }
            
        }
        
        fp_To_parabolic.SetBlock(element_index, block_phi);
        fgrad_p_To_parabolic.SetBlock(element_index, block_grad_phi);
        
        element_index++;
    }
    
//    fp_To_parabolic.Print(" p_to_p ");
//    fgrad_p_To_parabolic.Print(" grad_p_to_p ");
    
    return;
    
}

void TPZTransferFunctions::Fill_parabolic_To_elliptic(TPZCompMesh * parabolic, TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!parabolic || !elliptic) {
        DebugStop();
    }
#endif
    
    fp_e_cindexes.resize(0);
    
    
    parabolic->LoadReferences();
    TPZGeoMesh * geometry = parabolic->Reference();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<long, std::pair <long,long> > chunk_indexes;
    
    // inserting geometric indexes
    for (long i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        chunk_indexes.first = gel->Index();
        chunk_indexes.second.first = -1;
        chunk_indexes.second.second = -1;
        fp_e_cindexes.Push(chunk_indexes);
        
    }
    
    // counting volumetric elements
    long n_elements = fp_e_cindexes.size();
    
    // inserting origin indexes
    parabolic->LoadReferences();
    for(long iel = 0; iel < n_elements; iel++){
        
        TPZCompEl * cel = geometry->Element(fp_e_cindexes[iel].first)->Reference();
        
#ifdef PZDEBUG
        if (!cel) {
            //continue;
            DebugStop();
        }
#endif
        fp_e_cindexes[iel].second.first = cel->Index();
        
    }
    
    // inserting target indexes
    elliptic->LoadReferences();
    for(long iel = 0; iel < n_elements; iel++){
        
        TPZCompEl * cel = geometry->Element(fp_e_cindexes[iel].first)->Reference();
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        fp_e_cindexes[iel].second.second = cel->Index();
        
    }
    
    //    for (int k = 0; k < fp_e_cindexes.size(); k++) {
    //        std::cout << " volume k : " << k <<std::endl;
    //        std::cout << " volume gel : " << fp_e_cindexes[k].first <<std::endl;
    //        std::cout << " volume c_o : " << fp_e_cindexes[k].second.first <<std::endl;
    //        std::cout << " volume c_t : " << fp_e_cindexes[k].second.second <<std::endl;
    //    }
    
    
    parabolic->LoadReferences();
    geometry = parabolic->Reference();
    
    long nel = fp_e_cindexes.size();
    int n_var_dim = parabolic->Reference()->Dimension();
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    
    fp_dof_scatter.resize(nel);
    std::pair<long, std::pair <TPZVec<long>, TPZVec<long> > > chunk_intp_indexes;
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions_phi(nel);
    TPZVec< std::pair<long, long> > blocks_dimensions_grad_phi(nel);
    
    long gel_index, p_index, e_index;
    for (long iel = 0; iel < nel; iel++) {
        
        gel_index = fp_e_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        p_index   = fp_e_cindexes[iel].second.first;
        e_index   = fp_e_cindexes[iel].second.second;
        
        TPZCompEl * p_cel = parabolic->Element(p_index);
        TPZCompEl * e_cel = elliptic->Element(e_index);
        
#ifdef PZDEBUG
        if (!p_cel || !e_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel || !mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        // Getting local integration index
        TPZManVector<long> p_int_point_indexes(0,0);
        TPZManVector<long> e_int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
//        int gel_dim = gel->Dimension();
        
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        
        chunk_intp_indexes.first = gel->Index();
        chunk_intp_indexes.second.first  = p_int_point_indexes;
        chunk_intp_indexes.second.second = e_int_point_indexes;
        fp_e_intp_indexes.Push(chunk_intp_indexes);
        
        this->ElementDofIndexes(mf_p_cel, dof_indexes);
        fp_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions_phi[element_index].first = e_int_point_indexes.size();
        blocks_dimensions_phi[element_index].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[element_index].first = e_int_point_indexes.size()*n_var_dim;
        blocks_dimensions_grad_phi[element_index].second = dof_indexes.size();
        
        element_index++;
    }
    
    // Initialize the matrix
    fp_To_elliptic.Initialize(blocks_dimensions_phi);
    fgrad_p_To_elliptic.Initialize(blocks_dimensions_grad_phi);
    
    
    parabolic->LoadReferences();
    element_index = 0;
    
    TPZManVector<long> int_point_indexes(0,0);
    TPZManVector<long> dof_indexes(0,0);
    
    std::pair<long, long> block_dim;
    for (long iel = 0; iel < nel; iel++) {
        
        gel_index = fp_e_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        p_index   = fp_e_cindexes[iel].second.first;
        e_index   = fp_e_cindexes[iel].second.second;
        
        TPZCompEl * p_cel = parabolic->Element(p_index);
        TPZCompEl * e_cel = elliptic->Element(e_index);
        
#ifdef PZDEBUG
        if (!p_cel || !e_cel) {
            DebugStop();
        }
#endif
        
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel || !mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        
        TPZInterpolationSpace * p_intel = dynamic_cast<TPZInterpolationSpace * >(mf_p_cel->Element(0));
        
        // Getting local integration index
        int_point_indexes = fp_e_intp_indexes[element_index].second.second;
        dof_indexes = fp_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_geomechanic = mf_e_cel->GetIntegrationRule();
        int np_cel = int_points_geomechanic.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi, block_grad_phi;
        
        block_phi.Resize(block_dim.first*n_var_dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*n_var_dim*gel_dim,block_dim.second);
        
        // for derivatives in real space
        int nshape = p_intel->NShapeF();
        TPZFNMatrix<220> phi(nshape,1);
        TPZFNMatrix<660> dphi(gel_dim,nshape),dphix_axes(gel_dim,nshape);
        TPZFMatrix<double> dphidx;
        TPZFNMatrix<9,STATE> jacobian(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> jacinv(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> axes;
        REAL detjac;
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(gel_dim,0.0);
            STATE w;
            int_points_geomechanic.Point(ip, qsi, w);
            
            // Get the phi and dphix for H1 elasticity
            p_intel->Shape(qsi, phi, dphi);
            gel->Jacobian( qsi, jacobian, axes, detjac , jacinv);
            
            switch(gel_dim) {
                case 0:
                    break;
                case 1:
                    dphix_axes = dphi;
                    dphix_axes *= (1./detjac);
                    break;
                case 2:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
                    }
                    break;
                case 3:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
                        dphix_axes(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
                    }
                    break;
                default:
                    std::stringstream sout;
                    sout << "pzintel.c please implement the " << gel_dim << "d Jacobian and inverse\n";
                    LOGPZ_ERROR(logger,sout.str());
            }
            
            TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            
#ifdef PZDEBUG
            if(block_dim.second != phi.Rows()){
                DebugStop();
            }
#endif
            for (int jp = 0; jp < phi.Rows(); jp++) {
                block_phi(ip,jp) = phi(jp,0);
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    block_grad_phi(ip*n_var_dim + id,jp) = dphidx(id,jp);
                }
            }
            
        }
        
        fp_To_elliptic.SetBlock(element_index, block_phi);
        fgrad_p_To_elliptic.SetBlock(element_index, block_grad_phi);
        
        element_index++;
    }
    
//    fp_To_elliptic.Print(" p_to_e ");
//    fgrad_p_To_elliptic.Print(" grad_p_to_e ");
    
    return;
    
}

void TPZTransferFunctions::space_To_parabolic(TPZCompMesh * parabolic){
    
    
#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    parabolic->LoadReferences();
    TPZGeoMesh * geometry = parabolic->Reference();
    long nel = geometry->NElements();
    int dim = parabolic->Dimension();
    
    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    TPZFMatrix<double> block_phi;
    TPZFMatrix<double> block_grad_phi;
    
    for (int iel = 0; iel < nel; iel++) {
        
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        TPZCompEl * p_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fp_To_parabolic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_p_To_parabolic.GetSizeofBlock(iblock);
        
        fp_To_parabolic.GetBlock(iblock, block_phi);
        fgrad_p_To_parabolic.GetBlock(iblock, block_grad_phi);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = parabolic->FindMaterial(matd_id);
            TPZMatWithMem<TPZDarcyFlowMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZDarcyFlowMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            int gel_dim = gel->Dimension();
            int nshapes = b_size_phi.second;
            TPZFNMatrix<3,STATE> phi_p(nshapes,1,0.0);
            TPZFNMatrix<9,STATE> grad_phi_p(dim,nshapes);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < nshapes; is++) {
                    phi_p(is,0) = block_phi(ip,is);
                    
                }
                
                for (int is = 0; is < nshapes; is++) {
                    for (int id = 0 ; id < dim; id++) {
                        grad_phi_p(id,is) = block_grad_phi(ip*gel_dim + id,is);
                    }
                    
                }
                
                associated_material->GetMemory()[ipos].Set_phi_p(phi_p);
                associated_material->GetMemory()[ipos].Set_grad_phi_p(grad_phi_p);
            }
            
            
        }
        else{
            
            TPZMaterial * material = parabolic->FindMaterial(matd_id);
            TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            int gel_dim = gel->Dimension();
            int nshapes = b_size_phi.second;
            TPZFNMatrix<3,STATE> phi_p(nshapes,1,0.0);
            TPZFNMatrix<9,STATE> grad_phi_p(dim,nshapes);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < nshapes; is++) {
                    phi_p(is,0) = block_phi(ip,is);
                    
                }
                
                for (int is = 0; is < nshapes; is++) {
                    for (int id = 0 ; id < dim; id++) {
                        grad_phi_p(id,is) = block_grad_phi(ip*gel_dim + id,is);
                    }
                    
                }
                
                associated_material->GetMemory()[ipos].Set_phi_p(phi_p);
                associated_material->GetMemory()[ipos].Set_grad_phi_p(grad_phi_p);
            }
            
        }
        
        iblock++;
        
    }
    
    
}


void TPZTransferFunctions::parabolic_To_parabolic(TPZCompMesh * parabolic){
    
    
#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_p(fp_To_parabolic.Cols(),1,0.0);
    
    long pos = 0;
    int n = fp_e_cindexes.size();
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fp_dof_scatter[i].size(); iequ++) {
            Scatter_p(pos,0) = parabolic->Solution()(fp_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    
    // Step two
    TPZFMatrix<STATE> p_at_parabolic,grad_p_at_parabolic;
    fp_To_parabolic.Multiply(Scatter_p,p_at_parabolic);
    fgrad_p_To_parabolic.Multiply(Scatter_p, grad_p_at_parabolic);
    
    parabolic->LoadReferences();
    TPZGeoMesh * geometry = parabolic->Reference();
    long nel = geometry->NElements();
    int dim = parabolic->Dimension();
    
    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    for (int iel = 0; iel < nel; iel++) {
        
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        TPZCompEl * p_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fp_To_parabolic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_p_To_parabolic.GetSizeofBlock(iblock);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = parabolic->FindMaterial(matd_id);
            TPZMatWithMem<TPZDarcyFlowMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZDarcyFlowMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            REAL p;
            TPZFNMatrix<9,STATE> grad_p(dim,1);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                p = p_at_parabolic(first_point_phi + ip,0);
                
                for (int id = 0; id < dim ; id++) {
                    grad_p(id,0)= grad_p_at_parabolic(first_point_dphi + ip*dim + id,0);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_p_n(p);
                    associated_material->GetMemory()[ipos].Set_grad_p_n(grad_p);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_p(p);
                    associated_material->GetMemory()[ipos].Set_grad_p(grad_p);
                }
            }
            
            
        }
        else{
            
            TPZMaterial * material = parabolic->FindMaterial(matd_id);
            TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            REAL p;
            TPZFNMatrix<9,STATE> grad_p(dim,1);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                p = p_at_parabolic(first_point_phi + ip,0);
                
                for (int id = 0; id < dim ; id++) {
                    grad_p(id,0)= grad_p_at_parabolic(first_point_dphi + ip*dim + id,0);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_p_n(p);
                    associated_material->GetMemory()[ipos].Set_grad_p_n(grad_p);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_p(p);
                    associated_material->GetMemory()[ipos].Set_grad_p(grad_p);
                }
                

            }
            
        }
        
        iblock++;
        
    }
    
    
}

void TPZTransferFunctions::parabolic_To_elliptic(TPZCompMesh * parabolic, TPZCompMesh * elliptic){
    
    
#ifdef PZDEBUG
    if (!parabolic || !elliptic) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = parabolic->Reference();
    long nel = fp_e_cindexes.size();
    int dim = parabolic->Dimension();
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_p(fp_To_elliptic.Cols(),1,0.0);
    
    long pos = 0;
    for (int el = 0; el < nel; el++) {
        for(int iequ = 0; iequ < fp_dof_scatter[el].size(); iequ++) {
            Scatter_p(pos,0) = parabolic->Solution()(fp_dof_scatter[el][iequ],0);
            pos++;
        }
    }
    
    
    // Step two
    TPZFMatrix<STATE> p_at_elliptic,grad_p_at_elliptic;
    fp_To_elliptic.Multiply(Scatter_p,p_at_elliptic);
    fgrad_p_To_elliptic.Multiply(Scatter_p, grad_p_at_elliptic);
    
    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    long gel_index, e_index, p_index;
    for (int iel = 0; iel < nel; iel++) {
        
        
        gel_index = fp_e_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        p_index   = fp_e_cindexes[iel].second.first;
        e_index   = fp_e_cindexes[iel].second.second;
        
        TPZCompEl * p_cel = parabolic->Element(p_index);
        TPZCompEl * e_cel = elliptic->Element(e_index);
        
#ifdef PZDEBUG
        if (!p_cel || !e_cel) {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fp_To_elliptic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_p_To_elliptic.GetSizeofBlock(iblock);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TPZElasticBiotMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZElasticBiotMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            int_point_indexes = fp_e_intp_indexes[iblock].second.second;
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            REAL p;
            TPZFNMatrix<9,STATE> grad_p(dim,1);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                p = p_at_elliptic(first_point_phi + ip,0);
                
                for (int id = 0; id < dim ; id++) {
                    grad_p(id,0)= grad_p_at_elliptic(first_point_dphi + ip*dim + id,0);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_p_n(p);
                    associated_material->GetMemory()[ipos].Set_grad_p_n(grad_p);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_p(p);
                    associated_material->GetMemory()[ipos].Set_grad_p(grad_p);
                }
            }
            
            
        }
        else{
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            int_point_indexes = fp_e_intp_indexes[iblock].second.second;
            int n_points = int_point_indexes.size();
            long ipos;
            
            
            REAL p;
            TPZFNMatrix<9,STATE> grad_p(dim,1);
            for(long ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                p = p_at_elliptic(first_point_phi + ip,0);
                
                for (int id = 0; id < dim ; id++) {
                    grad_p(id,0)= grad_p_at_elliptic(first_point_dphi + ip*dim + id,0);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_p_n(p);
                    associated_material->GetMemory()[ipos].Set_grad_p_n(grad_p);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_p(p);
                    associated_material->GetMemory()[ipos].Set_grad_p(grad_p);
                }
            }
            
        }
        
        iblock++;
        
    }
    
    
}

void TPZTransferFunctions::RB_basis_To_Geomechanic_Memory(TPZCompMesh * cmesh_multiphysics){
    
    
    this->Fill_RB_basis_To_Geomechanic(cmesh_multiphysics);
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics || fgeomechanic_galerkinp_cindexes.size() == 0) {
        DebugStop();
    }
#endif
    

    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    long nel = fgeomechanic_galerkinp_cindexes.size();
    int dim = cmesh_multiphysics->Dimension();
    int n_rb = fCmeshRB_projections->Solution().Cols();
    

    // Step zero scatter
    TPZFMatrix<STATE> ScatterDisplacements(fphi_u_To_Geomechanic.Cols(),n_rb,0.0);

    //    for all RB functions
    for(int i_rb = 0; i_rb < n_rb; i_rb++){
        long pos = 0;
        for (int el = 0; el < nel; el++) {
            for(int ip = 0; ip < fphi_u_dof_scatter[el].size(); ip++) {
                ScatterDisplacements(pos,i_rb) = fCmeshRB_projections->Solution()(fphi_u_dof_scatter[el][ip],i_rb);
                pos++;
            }
        }
        
    }
    
    
    // Step two
    TPZFMatrix<STATE> u_at_intpoints,grad_u_at_intpoints;
    fphi_u_To_Geomechanic.Multiply(ScatterDisplacements,u_at_intpoints);
    fgrad_phi_u_To_Geomechanic.Multiply(ScatterDisplacements, grad_u_at_intpoints);
    
    
    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    long gel_index, geomech_index, gp_index;
    for (int iel = 0; iel < nel; iel++) {
        
        
        gel_index = fgeomechanic_galerkinp_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        geomech_index   = fgeomechanic_galerkinp_cindexes[iel].second.first;
        gp_index        = fgeomechanic_galerkinp_cindexes[iel].second.second;
        
        TPZCompEl * geomec_cel = cmesh_multiphysics->Element(geomech_index);
        TPZCompEl * gp_cel = fCmeshRB_projections->Element(gp_index);
        
#ifdef PZDEBUG
        if (!geomec_cel || !gp_cel) {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fphi_u_To_Geomechanic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_phi_u_To_Geomechanic.GetSizeofBlock(iblock);
        iblock++;
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = cmesh_multiphysics->FindMaterial(matd_id);
            TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> *>(material);
            
            long np_cmesh = associated_material->GetMemory().NElements();

            for(long i = 0; i <  np_cmesh; i++) {
                associated_material->GetMemory()[i].phi_u().Resize(n_rb,dim);
                associated_material->GetMemory()[i].grad_phi_u().Resize(n_rb,dim*dim);
            }
            
            TPZManVector<long, 30> int_point_indexes;
            GlobalPointIndexes(geomec_cel, int_point_indexes);

            // Transfering values
            int n_points = int_point_indexes.size();
            TPZFMatrix<STATE> phi_u(n_rb,dim,0.0);
            TPZFMatrix<STATE> grad_phi_u(n_rb,dim*dim);
            long i_pos;
            for(long ip = 0; ip <  n_points; ip++){
                i_pos = int_point_indexes[ip];
                
                //    for all RB functions
                for(int i_rb = 0; i_rb < n_rb; i_rb++){
                    for (int id = 0; id < dim ; id++) {
                        phi_u(i_rb,id) = u_at_intpoints(first_point_phi + ip*dim + id,i_rb);
                    }
                    
                    int c_pos = 0;
                    for (int id = 0; id < dim ; id++) {
                        for (int jd = 0; jd < dim ; jd++) {
                            grad_phi_u(i_rb,c_pos)= grad_u_at_intpoints(first_point_dphi + ip*dim*dim + id*dim + jd,i_rb);
                            c_pos++;
                        }
                    }
                }
                associated_material->GetMemory()[int_point_indexes[ip]].Set_phi_u(phi_u);
                associated_material->GetMemory()[int_point_indexes[ip]].Set_grad_phi_u(grad_phi_u);
            }
            
        }
        else{
            
            TPZMaterial * material = cmesh_multiphysics->FindMaterial(matd_id);
            TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> *>(material);
            
            long np_cmesh = associated_material->GetMemory().NElements();
            
            for(long i = 0; i <  np_cmesh; i++) {
                associated_material->GetMemory()[i].phi_u().Resize(n_rb,dim);
            }
            
            TPZManVector<long, 30> int_point_indexes;
            cmesh_multiphysics->LoadReferences();
            TPZCompEl * mf_cel = gel->Reference();
            GlobalPointIndexes(mf_cel, int_point_indexes);
            
            // Transfering values
            int n_points = int_point_indexes.size();
            TPZFMatrix<STATE> phi_u(n_rb,dim,0.0);
            long i_pos;
            for(long ip = 0; ip <  n_points; ip++){
                i_pos = int_point_indexes[ip];
                
                //    for all RB functions
                for(int i_rb = 0; i_rb < n_rb; i_rb++){
                    for (int id = 0; id < dim ; id++) {
                        phi_u(i_rb,id) = u_at_intpoints(first_point_phi + ip*dim + id,i_rb);
                    }
                }
                associated_material->GetMemory()[i_pos].Set_phi_u(phi_u);
            }
            
        }

    }
    
}

/** @brief Transfer the RB Solution to multiphysics mesh  */
void TPZTransferFunctions::RB_Solution_To_Geomechanic(TPZCompMesh * cmesh_multiphysics, TPZFMatrix<STATE> & rb_solution){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics || fgeomechanic_galerkinp_cindexes.size() == 0) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    long nel = fgeomechanic_galerkinp_cindexes.size();
    int dim = cmesh_multiphysics->Dimension();
    int n_gp = fCmeshRB_projections->Solution().Cols();
    
#ifdef PZDEBUG
    if (n_gp != rb_solution.Rows()) {
        DebugStop();
    }
#endif
    
    
    std::pair<long, long> block_size;
    block_size.first = 0;
    block_size.second = 0;
    
    
    long gel_index, geomech_index, gp_index;
    for (int iel = 0; iel < nel; iel++) {
        
        
        gel_index = fgeomechanic_galerkinp_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        geomech_index   = fgeomechanic_galerkinp_cindexes[iel].second.first;
        gp_index        = fgeomechanic_galerkinp_cindexes[iel].second.second;
        
        TPZCompEl * geomec_cel = cmesh_multiphysics->Element(geomech_index);
        TPZCompEl * gp_cel = fCmeshRB_projections->Element(gp_index);
        
#ifdef PZDEBUG
        if (!geomec_cel || !gp_cel) {
            DebugStop();
        }
#endif
        
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = cmesh_multiphysics->FindMaterial(matd_id);
            TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> *>(material);
            
            
            TPZManVector<long, 30> int_point_indexes;
            GlobalPointIndexes(geomec_cel, int_point_indexes);
            
            // Transfering values
            long n_points = int_point_indexes.size();
            TPZFNMatrix<3,REAL> u(dim,1);
            TPZFNMatrix<9,REAL> grad_u(dim,dim);
            
            long i_pos;
            for(long ip = 0; ip <  n_points; ip++){
                i_pos = int_point_indexes[ip];
                u.Zero();
                grad_u.Zero();

                 for (int igp = 0; igp <n_gp; igp++) {
                     int c_pos = 0;
                     for (int id = 0; id <dim; id++) {
                        u(id,0) += associated_material->GetMemory()[i_pos].phi_u()(igp,id) * rb_solution(igp,0);
                        for (int jd = 0; jd <dim; jd++) {
                            grad_u(id,jd) += associated_material->GetMemory()[i_pos].grad_phi_u()(igp,c_pos) * rb_solution(igp,0);
                            c_pos++;
                        }
                    }
                }
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[i_pos].Set_u_n(u);
                    associated_material->GetMemory()[i_pos].Set_grad_u_n(grad_u);
                }
                else{
                    associated_material->GetMemory()[i_pos].Set_u(u);
                    associated_material->GetMemory()[i_pos].Set_grad_u(grad_u);                    
                }
            }
            
        }
        else{
            
            TPZMaterial * material = cmesh_multiphysics->FindMaterial(matd_id);
            TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            GlobalPointIndexes(geomec_cel, int_point_indexes);
            
            // Transfering values
            long n_points = int_point_indexes.size();
            TPZFNMatrix<3,REAL> u(dim,1);
            
            long i_pos;
            for(long ip = 0; ip <  n_points; ip++){
                i_pos = int_point_indexes[ip];
                u.Zero();
                for (int igp = 0; igp <n_gp; igp++) {
                    int c_pos = 0;
                    for (int id = 0; id <dim; id++) {
                        u(id,0) += associated_material->GetMemory()[i_pos].phi_u()(igp,id) * rb_solution(igp,0);
                    }
                }
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[i_pos].Set_u_n(u);
                }
                else{
                    associated_material->GetMemory()[i_pos].Set_u(u);
                }
            }
            
        }
        
    }
    
    
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Utility Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief Get Global integration point indexes associaded  */
void TPZTransferFunctions::GlobalPointIndexes(TPZCompEl * cel, TPZManVector<long,30> &int_point_indexes){
    
    TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
    
#ifdef PZDEBUG
    if(!mf_cel)
    {
        DebugStop();
    }
#endif
    
    mf_cel->GetMemoryIndices(int_point_indexes);
    
}

/** @brief Get Global integration point indexes associaded  */
void TPZTransferFunctions::GlobalPointIndexesInterface(TPZCompEl * int_cel, TPZManVector<long,30> &int_point_indexes){
    
    TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(int_cel);
    
#ifdef PZDEBUG
    if(!mf_int_cel)
    {
        DebugStop();
    }
#endif
    
    mf_int_cel->GetMemoryIndices(int_point_indexes);
    
}

bool TPZTransferFunctions::IdentifyFace(int &side, TPZGeoEl * vol, TPZGeoEl * face){
    
    int volu_nsides = vol->NSides();
    int face_nsides = face->NSides();
    side = -1;
    TPZGeoElSide face_itself =  face->Neighbour(face_nsides-1);
    bool IsMembershipQ = false;
    
    for (int iside = 0; iside < volu_nsides; iside++) {
        IsMembershipQ = bool(vol->NeighbourExists(iside, face_itself));
        if (IsMembershipQ) {
            side = iside;
            break;
        }
    }
    
    TPZGeoElSide vol_itself =  vol->Neighbour(volu_nsides-1);
    
    if(!IsMembershipQ){
        TPZGeoElSide neigh = face->Neighbour(face_nsides-1);
        
        if(!neigh.Element()){
            DebugStop();
        }
        
        TPZGeoElSide neigh_father = neigh.Father2();
        
        if(!neigh_father.Element()){
            DebugStop();
        }
        bool IsNeighQ = false;
        for (int iside = vol_itself.Element()->NNodes(); iside < volu_nsides; iside++) {
            
            IsNeighQ = bool(vol_itself.Element()->NeighbourExists(iside, neigh_father));
            
            if (IsNeighQ) {
                side = iside;
                IsMembershipQ = IsNeighQ;
                break;
            }
        }
    }
    
    return IsMembershipQ;
}

/** @brief Compute parametric transformation form origin to tarfet (xinverse based) */
void TPZTransferFunctions::ComputeTransformation(TPZGeoEl * face_gel_origin, TPZGeoEl * gel_origin , TPZGeoEl * gel_target, TPZVec<REAL> & origin, TPZVec<REAL> & target){
    
#ifdef PZDEBUG
    if (!face_gel_origin || !gel_origin || !gel_target) {
        DebugStop();
    }
    
    if (gel_origin->Dimension() != gel_target->Dimension()) {
        DebugStop();
    }
    
#endif
    
    bool IsmemberQ;
    int face_side;
    TPZManVector<REAL,3> origin_vol(gel_origin->Dimension());
    target.Resize(gel_origin->Dimension(),0.0);
    
    // Transformation step one
    TPZTransform<REAL> tr_face_to_vol(gel_origin->Dimension());
    IsmemberQ =  IdentifyFace(face_side, gel_origin, face_gel_origin);
    
    
#ifdef PZDEBUG
    if (!IsmemberQ) {
        DebugStop();
    }
#endif
    
    tr_face_to_vol = gel_origin->SideToSideTransform(face_side, gel_origin->NSides()-1);
    tr_face_to_vol.Apply(origin, origin_vol);
    gel_origin->TransformSonToFather(gel_target, origin_vol, target);
    
    //    std::cout << "origin = "        << origin << std::endl;
    //    std::cout << "origin_vol = "    << origin_vol << std::endl;
    //    std::cout << "target = "        << target << std::endl;
    
}

/** @brief Compute indices associated to faces on 3D topologies */
void TPZTransferFunctions::ComputeFaceIndex(TPZGeoEl * gel , TPZVec<int> &sides){
    
    
    switch (gel->Type()) {
        case ECube:
        {
            int nfaces = 6;
            sides.Resize(nfaces);
            sides[0] = 20;
            sides[1] = 21;
            sides[2] = 22;
            sides[3] = 23;
            sides[4] = 24;
            sides[5] = 25;
            
        }
            break;
        case ETetraedro:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 10;
            sides[1] = 11;
            sides[2] = 12;
            sides[3] = 13;
            
        }
            break;
        case EQuadrilateral:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 4;
            sides[1] = 5;
            sides[2] = 6;
            sides[3] = 7;
            
        }
            break;
        case ETriangle:
        {
            int nfaces = 3;
            sides.Resize(nfaces);
            sides[0] = 3;
            sides[1] = 4;
            sides[2] = 5;
            
        }
            break;
        default:
        {
            std::cout << "Element not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compute sides associated to faces on 3D topologies */
void TPZTransferFunctions::ComputeFaceNormals(TPZGeoEl * gel , TPZVec<int> &sides, TPZFMatrix<STATE> &normals){
    
    //  @omar:: Just for linear mapping
    
    TPZFMatrix<REAL> mat_normals;
    TPZVec<int> v_sides;
    gel->ComputeNormals(mat_normals, v_sides);
    
    switch (gel->Type()) {
        case ECube:
        {
            int nfaces = 6;
            sides.Resize(nfaces);
            sides[0] = 20;
            sides[1] = 21;
            sides[2] = 22;
            sides[3] = 23;
            sides[4] = 24;
            sides[5] = 25;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
        case ETetraedro:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 10;
            sides[1] = 11;
            sides[2] = 12;
            sides[3] = 13;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
        case EQuadrilateral:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 4;
            sides[1] = 5;
            sides[2] = 6;
            sides[3] = 7;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
        case ETriangle:
        {
            int nfaces = 3;
            sides.Resize(nfaces);
            sides[0] = 3;
            sides[1] = 4;
            sides[2] = 5;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
            
        default:
        {
            std::cout << "Element not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compute left and right geometric element indexes */
void TPZTransferFunctions::ComputeLeftRight(TPZCompMesh * transport_mesh){
    
    DebugStop();
    
}


/** @brief Compute left and right geometric element indexes */
void TPZTransferFunctions::ComputeLeftRightII(TPZCompMesh * transport_mesh){
    
    DebugStop();
    
}


/** @brief Dimensionla Measure of the elemnt */
REAL TPZTransferFunctions::DimensionalMeasure(TPZGeoEl * gel){
    
#ifdef PZDEBUG
    if (!gel) {
        DebugStop();
    }
#endif
    REAL measure = 0.0;
    int order = 10;
    int element_itself  = gel->NSides() - 1;
    TPZIntPoints * int_points = gel->CreateSideIntegrationRule(element_itself, order);
    REAL detjac, w;
    TPZVec<REAL> par(gel->Dimension(),0.0);
    TPZFMatrix<REAL> jac;
    TPZFMatrix<REAL> axes;
    TPZFMatrix<REAL> jacinv;
    for (int i = 0; i < int_points->NPoints(); i++) {
        int_points->Point(i, par, w);
        gel->Jacobian(par, jac, axes, detjac, jacinv);
        measure += w * detjac;
    }
    
    return measure;
    
}


void TPZTransferFunctions::ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
//    int dfvar = block.Size(dfseq);
//    long pos = block.Position(dfseq);
//    for(int jn=0; jn<dfvar; jn++) {
//        for (long is=0; is<numbersol; is++) {
//            sol[is][iv%nstate] += (STATE)phi(iv/nstate,0)*MeshSol(pos+jn,is);
//            for(d=0; d<dphix.Rows(); d++){
//                dsol[is](d,iv%nstate) += (STATE)dphix(d,iv/nstate)*MeshSol(pos+jn,is);
//            }
//        }
//        iv++;
//    }
    
    TPZStack<long> index(0,0);
    int nconnect = intel->NConnects();
    
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = intel->Mesh()->Block().Position(seqnumber);
        long b_size = intel->Mesh()->Block().Size(seqnumber);
            for(int jb=0; jb<b_size; jb++) {
                index.Push(position + jb);
            }
    }
    
    dof_indexes = index;
    return;
}

void TPZTransferFunctions::ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!m_el) {
        DebugStop();
    }
#endif
    
    TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(m_el->Element(0));
    
#ifdef PZDEBUG
    if (!intel_vol) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    int nconnect = intel_vol->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = m_el->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = m_el->Mesh()->Block().Position(seqnumber);
        int b_size = m_el->Mesh()->Block().Size(seqnumber);
        for (int ib=0; ib < b_size; ib++) {
            index.Push(position+ ib);
        }
    }
    
    dof_indexes = index;
    return;
    
}

void TPZTransferFunctions::ElementDofFaceIndexes(int connect_index,TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    TPZConnect  & con = intel->Connect(connect_index);
    long seqnumber = con.SequenceNumber();
    long position = intel->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}

void TPZTransferFunctions::ElementDofFaceIndexes(int connect_index, TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!m_el && connect_index > 4) {
        DebugStop();
    }
#endif
    
    
    TPZStack<long> index(0,0);
    TPZConnect  & con = m_el->Connect(connect_index);
    long seqnumber = con.SequenceNumber();
    long position = m_el->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}


/** @brief Compute compuational mesh pair (geomechanic, gp mesh) indexed by geometric element index */
void TPZTransferFunctions::FillGeomechanicElPairs(TPZCompMesh * cmesh_mphysics){
    

    fgeomechanic_galerkinp_cindexes.Resize(0);
    
#ifdef PZDEBUG
    if (!cmesh_mphysics || !fCmeshRB_projections) {
        DebugStop();
    }
    
#endif

    fCmeshRB_projections->LoadReferences();
    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif

    std::pair<long, std::pair <long,long> > gel_indexes;
    
    for (long i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        gel_indexes.first = gel->Index();
        gel_indexes.second.first = -1;
        gel_indexes.second.second = -1;
        fgeomechanic_galerkinp_cindexes.Push(gel_indexes);
        
    }
    
    // counting volumetric elements
    long n_elements = fgeomechanic_galerkinp_cindexes.size();
    fgeomechanic_galerkinp_cindexes.Resize(n_elements);
    
    // inserting geomechanic indexes
    cmesh_mphysics->LoadReferences();
    for(long iel = 0; iel < n_elements; iel++){
        
        TPZCompEl * cel = geometry->Element(fgeomechanic_galerkinp_cindexes[iel].first)->Reference();
        
#ifdef PZDEBUG
        if (!cel) {
            //continue;
            DebugStop();
        }
#endif
        fgeomechanic_galerkinp_cindexes[iel].second.first = cel->Index();
        
    }
    
    // inserting gp indexes
    fCmeshRB_projections->LoadReferences();
    for(long iel = 0; iel < n_elements; iel++){
        
        TPZCompEl * gp_cel = geometry->Element(fgeomechanic_galerkinp_cindexes[iel].first)->Reference();
        
#ifdef PZDEBUG
        if (!gp_cel) {
            DebugStop();
        }
#endif
        fgeomechanic_galerkinp_cindexes[iel].second.second = gp_cel->Index();
        
    }
    
//    for (int k = 0; k < fgeomechanic_galerkinp_cindexes.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << fgeomechanic_galerkinp_cindexes[k].first <<std::endl;
//        std::cout << " volume cgeo : " << fgeomechanic_galerkinp_cindexes[k].second.first <<std::endl;
//        std::cout << " volume cgp : " << fgeomechanic_galerkinp_cindexes[k].second.second <<std::endl;
//    }
    
}


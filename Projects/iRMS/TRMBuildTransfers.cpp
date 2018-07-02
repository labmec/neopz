//
//  TRMBuildTransfers.cpp
//  PZ
//
//  Created by Omar on 10/27/15.
//
//


#include "TRMBuildTransfers.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


/** @brief Default constructor */
TRMBuildTransfers::TRMBuildTransfers(){
    
    fSimulationData = NULL;

}

/** @brief Default desconstructor */
TRMBuildTransfers::~TRMBuildTransfers(){

}

/** @brief Copy constructor $ */
TRMBuildTransfers::TRMBuildTransfers(const TRMBuildTransfers &copy)
{
    fSimulationData = copy.fSimulationData;
}

/** @brief Copy assignemnt operator $ */
TRMBuildTransfers & TRMBuildTransfers::operator=(const TRMBuildTransfers &other)
{
    if (this != & other) // prevent self-assignment
    {
        fSimulationData = other.fSimulationData;
    }
    return *this;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Transfers:: Iterative Coupling by Operator Splitting //////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segregated Transfer methods (Gamma and Omega) :: Build methods Elliptic
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief bluid linear applications: u and grad_u to elliptic $ */
void TRMBuildTransfers::Build_elliptic_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    elliptic->LoadReferences();    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    fe_e_cindexes.resize(0);
    std::pair< int64_t, int64_t > chunk_geo_cel_indexes;
    
    // Step 1 :: Counting for valid elements (apply all the needed filters in this step)
    for (int64_t i = 0; i < geometry->NElements(); i++) {
        
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        int parabolic_res = fSimulationData->MHMResolution().second.second;
        if (gel->Level() != parabolic_res) {
            continue;
        }
        
        int mat_id = gel->MaterialId();
        if ((dim == 2 && (mat_id >= 8 && mat_id <= 11)) && (gel->Dimension() == dim-1)) { // Filtering bc reservoir elements
            continue;
        }
        
        if ((gel->Dimension() == dim-1) && (dim == 3 && (mat_id >= 8 && mat_id <= 13))) { // Filtering bc reservoir elements
            continue;
        }
        
        if (mat_id == fSimulationData->Skeleton_material_Id() || mat_id == fSimulationData->InterfacesMatId()) { // Filtering skeleton reservoir elements
            continue;
        }
        
        if ((gel->Dimension() == dim-1) && gel->NumInterfaces() !=0) { // Filtering interface reservoir elements for transport

            continue;
        }
        
        
        chunk_geo_cel_indexes.first = gel->Index();
        chunk_geo_cel_indexes.second = -1;
        fe_e_cindexes.push_back(chunk_geo_cel_indexes);
    }
    
    int64_t n_el = fe_e_cindexes.size();
    fu_dof_scatter.resize(n_el);
    fe_e_intp_indexes.resize(n_el);
    
    // Block size structue including (Omega and Gamma)
    std::pair<int64_t, TPZVec<int64_t> > chunk_intp_indexes;
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi(n_el);
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_grad_phi(n_el);
    
    // Step 2 :: filling linking vectors
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
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
        
        int gel_dim = gel->Dimension();
        
        // Geometry and cel link
        fe_e_cindexes[iel].second = e_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> e_int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        
        chunk_intp_indexes.first = gel->Index();
        chunk_intp_indexes.second  = e_int_point_indexes;
        fe_e_intp_indexes.push_back(chunk_intp_indexes);
        
        this->ElementDofIndexes(mf_e_cel, dof_indexes);
        fu_dof_scatter[iel] = dof_indexes;
        
        
        blocks_dimensions_phi[iel].first = e_int_point_indexes.size()*dim;
        blocks_dimensions_phi[iel].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[iel].first = e_int_point_indexes.size()*dim*gel_dim;
        blocks_dimensions_grad_phi[iel].second = dof_indexes.size();
    }
    
    // Step 3 :: Initialize the matrix
    fu_To_elliptic.Initialize(blocks_dimensions_phi);
    fgrad_u_To_elliptic.Initialize(blocks_dimensions_grad_phi);
    
    TPZManVector<int64_t> e_int_point_indexes(0,0);
    TPZManVector<int64_t> dof_indexes(0,0);
    
    // Step 4 :: Filling the matrix
    std::pair<int64_t, int64_t> block_dim;
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
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
        dof_indexes = fu_dof_scatter[iel];
        
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
        
        block_phi.Resize(block_dim.first*dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*dim*gel_dim,block_dim.second);
        
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
            
            // Get the phi and dphix for H1/rb elasticity
            e_intel->Shape(qsi, phi, dphi);
            
            if (!fSimulationData->ReducedBasisResolution().first){
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
//                        LOGPZ_ERROR(logger,sout.str());
                }
                
                TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            }
            
#ifdef PZDEBUG
            
            if(block_dim.second != phi.Rows() * dim){
                if (fSimulationData->ReducedBasisResolution().first) {
                    if (block_dim.second != phi.Rows()) {
                        DebugStop();
                    }
                }
                else{
                    DebugStop();
                }

            }
#endif
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < dim; id++) {
                    if (!fSimulationData->ReducedBasisResolution().first){
                        block_phi(ip*dim+id,jp*dim+id) = phi(jp,0);
                    }
                    else{
                        block_phi(ip*dim+id,jp) = phi(jp,id);
                    }
                }
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                int d_count = 0;
                for (int id = 0; id < dim; id++) {
                    for (int jd = 0; jd < gel_dim; jd++) {
                        if(gel_dim == dim){

                            if (!fSimulationData->ReducedBasisResolution().first){
                                block_grad_phi(ip*dim*gel_dim + id*gel_dim + jd,jp*dim+id) = dphidx(jd,jp);
                            }
                            else{
                                block_grad_phi(ip*dim*gel_dim + id*gel_dim + jd,jp) = dphi(jp,d_count);
                                d_count++;
                            }
                        }
                        else{
                            block_grad_phi(ip*dim*gel_dim+id*gel_dim + jd,jp) = 0.0;
                        }
                    }
                }
            }
            
        }
        
        fu_To_elliptic.SetBlock(iel, block_phi);
        fgrad_u_To_elliptic.SetBlock(iel, block_grad_phi);

    }
    
//    fu_To_elliptic.Print(" u_to_e ");
//    fgrad_u_To_elliptic.Print(" grad_u_to_e ");
    
    return;
    
}

void TRMBuildTransfers::space_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic || fu_To_elliptic.Rows() == 0 || fgrad_u_To_elliptic.Rows() == 0) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = elliptic->Dimension();
    int n_el = fe_e_cindexes.size();
    
//    int64_t first_point_phi = 0;
//    int64_t first_point_dphi = 0;
    std::pair<int64_t, int64_t> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    TPZFMatrix<double> block_phi;
    TPZFMatrix<double> block_grad_phi;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_e_cindexes[iel].first);
    
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif

        TPZCompEl * e_cel = elliptic->Element(fe_e_cindexes[iel].second);
        
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
        
//        first_point_phi += b_size_phi.first;
//        first_point_dphi += b_size_dphi.first;
        b_size_phi = fu_To_elliptic.GetSizeofBlock(iel);
        b_size_dphi = fgrad_u_To_elliptic.GetSizeofBlock(iel);
        fu_To_elliptic.GetBlock(iel, block_phi);
        fgrad_u_To_elliptic.GetBlock(iel, block_grad_phi);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = elliptic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            int gel_dim = gel->Dimension();
            int nshapes = b_size_phi.second/dim;
            
            if(fSimulationData->ReducedBasisResolution().first){

                nshapes = b_size_phi.second;
                TPZFNMatrix<3,STATE> phi_u(nshapes,dim,0.0);
                TPZFNMatrix<9,STATE> grad_phi_u(dim*dim,nshapes);
                for(int64_t ip = 0; ip <  n_points; ip++) {
                    ipos  = int_point_indexes[ip];
                    
                    for (int is = 0; is < nshapes; is++) {
                        for (int id = 0 ; id < dim; id++) {
                            phi_u(is,id) = block_phi(ip*dim+id,is);
                        }
                    }
                    
                    for (int is = 0; is < nshapes; is++) {
                        int dcout = 0;
                        for (int id = 0 ; id < dim; id++) {
                            for (int jd = 0 ; jd < gel_dim; jd++) {
                                grad_phi_u(dcout,is) = block_grad_phi(ip*dim*gel_dim + id*gel_dim + jd,is);
                                dcout++;
                            }
                        }
                        
                    }
                    
                    associated_material->GetMemory()[ipos].Set_phi_u(phi_u);
                    associated_material->GetMemory()[ipos].Set_grad_phi_u(grad_phi_u);
                }
                
            }
            else{
                TPZFNMatrix<3,STATE> phi_u(nshapes,1,0.0);
                TPZFNMatrix<9,STATE> grad_phi_u(dim,nshapes);
                for(int64_t ip = 0; ip <  n_points; ip++) {
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

            
            
        }
        else{
            
            TPZMatWithMem<TRMMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            
            int gel_dim = gel->Dimension();
            int nshapes = b_size_phi.second/dim;
            
            if(fSimulationData->ReducedBasisResolution().first){
                
                nshapes = b_size_phi.second;
                TPZFNMatrix<3,STATE> phi_u(nshapes,dim,0.0);
                TPZFNMatrix<9,STATE> grad_phi_u(dim*dim,nshapes);
                for(int64_t ip = 0; ip <  n_points; ip++) {
                    ipos  = int_point_indexes[ip];
                    
                    for (int is = 0; is < nshapes; is++) {
                        for (int id = 0 ; id < dim; id++) {
                            phi_u(is,id) = block_phi(ip*dim+id,is);
                        }
                    }
                    
                    for (int is = 0; is < nshapes; is++) {
                        int dcout = 0;
                        for (int id = 0 ; id < dim; id++) {
                            for (int jd = 0 ; jd < gel_dim; jd++) {
                                grad_phi_u(dcout,is) = block_grad_phi(ip*dim*gel_dim + id*gel_dim + jd,is);
                                dcout++;
                            }
                        }
                        
                    }
                    
                    associated_material->GetMemory()[ipos].Set_phi_u(phi_u);
                    associated_material->GetMemory()[ipos].Set_grad_phi_u(grad_phi_u);
                }
            }
            else{

                TPZFNMatrix<3,STATE> phi_u(nshapes,1,0.0);
                TPZFNMatrix<9,STATE> grad_phi_u(dim,nshapes);
                for(int64_t ip = 0; ip <  n_points; ip++) {
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
            
        }

        
    }
    
}

void TRMBuildTransfers::spatial_props_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
#endif
    
    int dim = elliptic->Dimension();
    
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;
    TPZManVector<STATE, 10> lambda,mu,S_e,alpha;
    
    // Step one
    int n_elements = elliptic->NElements();
    TPZManVector<int64_t, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = elliptic->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = elliptic->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);
            int map_model = fSimulationData->Map()->MapModel();
            if (map_model != 0 && (rockid == 12 || rockid == 14)) {
                fSimulationData->Map()->SetMapModel(0);
                fSimulationData->Map()->phi(x, porosity, vars);
                fSimulationData->Map()->lambda(x, lambda, vars);
                fSimulationData->Map()->mu(x, mu, vars);
                fSimulationData->Map()->S_e(x, S_e, vars);
                fSimulationData->Map()->alpha(x, alpha, vars);
                fSimulationData->Map()->SetMapModel(map_model);
            }
            else{

                fSimulationData->Map()->phi(x, porosity, vars);
                fSimulationData->Map()->lambda(x, lambda, vars);
                fSimulationData->Map()->mu(x, mu, vars);
                fSimulationData->Map()->S_e(x, S_e, vars);
                fSimulationData->Map()->alpha(x, alpha, vars);
                
                if (map_model != 0) {
                    lambda[0]   *= cos(porosity[0]);
                    mu[0]       *= cos(porosity[0]);
                    S_e[0]      *= sin(porosity[0]);
                    alpha[0]    *= sin(porosity[0]);
                }
                
            }
            
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
            associated_material->GetMemory()[indexes[ip]].Set_lambda(lambda[0]);
            associated_material->GetMemory()[indexes[ip]].Set_mu(mu[0]);
            associated_material->GetMemory()[indexes[ip]].Set_S_e(S_e[0]);
            associated_material->GetMemory()[indexes[ip]].Set_alpha(alpha[0]);
        }
    }
    
    
}

void TRMBuildTransfers::phi_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
#endif
    
    int dim = elliptic->Dimension();
    
    
    TPZFMatrix<STATE> kappa, kappa_inv;
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;
    
    // Step one
    int n_elements = elliptic->NElements();
    TPZManVector<int64_t, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = elliptic->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = elliptic->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);
            //            fSimulationData->Map()->ComputePropertieSPE10Map(index, x, kappa, kappa_inv, porosity);
//            fSimulationData->Map()->Kappa(x, kappa, kappa_inv, vars);
            fSimulationData->Map()->phi(x, porosity, vars);
//            associated_material->GetMemory()[indexes[ip]].Set_K_0(kappa);
//            associated_material->GetMemory()[indexes[ip]].Set_Kinv_0(kappa_inv);
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
        }
    }
    
}


void TRMBuildTransfers::elliptic_To_elliptic(TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_u(fu_To_elliptic.Cols(),1,0.0);
    int n = fe_e_cindexes.size();
    int64_t pos = 0;
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
    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = elliptic->Dimension();
    int n_el = fe_e_cindexes.size();
    
    int64_t first_point_phi = 0;
    int64_t first_point_dphi = 0;
    std::pair<int64_t, int64_t> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(fe_e_cindexes[iel].second);
        
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
        b_size_phi = fu_To_elliptic.GetSizeofBlock(iel);
        b_size_dphi = fgrad_u_To_elliptic.GetSizeofBlock(iel);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMaterial * material = elliptic->FindMaterial(matd_id);
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,3,0.0);
            TPZFNMatrix<9,STATE> grad_u(3,3,0.0);
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_elliptic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        if(fSimulationData->ReducedBasisResolution().first){
                            grad_u(id,jd) = grad_u_at_elliptic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                        }
                        else{
                            grad_u(jd,id) = grad_u_at_elliptic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                        }

                    }
                }
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_grad_u_0(grad_u);
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
            TPZMatWithMem<TRMMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,3,0.0);
            TPZFNMatrix<9,STATE> grad_u(3,3,0.0);
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_elliptic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        grad_u(jd,id)= grad_u_at_elliptic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                    }
                }
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_grad_u_0(grad_u);
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
        
    }
    
}

/** @brief bluid linear applications: u and grad_u to parabolic $ */
void TRMBuildTransfers::Build_elliptic_To_parabolic(TPZCompMesh * elliptic, TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!elliptic || !parabolic) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    elliptic->LoadReferences();
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    
    fe_p_cindexes.resize(0);
    std::pair<int64_t, std::pair<int64_t, int64_t> > chunk_geo_cel_indexes;
    
    // Step 1 :: Counting for valid elements (apply all the needed filters in this step)
    for (int64_t i = 0; i < geometry->NElements(); i++) {
        
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        int parabolic_res = fSimulationData->MHMResolution().second.second;
        if (gel->Level() != parabolic_res) {
            continue;
        }
        
        int mat_id = gel->MaterialId();
        if ((dim == 2 && mat_id == 12) || (gel->Dimension() == dim-1)) { // Filtering bc and sideburden elements
            continue;
        }
        
        if ((gel->Dimension() == dim-1) || (dim == 3 && mat_id == 14)) { // Filtering bc and sideburden elements
            continue;
        }
        
        if (mat_id == fSimulationData->Skeleton_material_Id() || mat_id == fSimulationData->InterfacesMatId()) { // Filtering skeleton reservoir elements
            continue;
        }
        
        if ((gel->Dimension() == dim-1) && gel->NumInterfaces() !=0) { // Filtering interface reservoir elements for transport
            continue;
        }
        
        chunk_geo_cel_indexes.first = gel->Index();
        chunk_geo_cel_indexes.second.first  = -1;
        chunk_geo_cel_indexes.second.second = -1;
        fe_p_cindexes.push_back(chunk_geo_cel_indexes);
    }
    
    int64_t n_el = fe_p_cindexes.size();
    fu_p_dof_scatter.resize(n_el);
    fe_p_intp_indexes.resize(n_el);
    
    
    // Inserting elliptic elements
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        fe_p_cindexes[iel].second.first = e_cel->Index();
    
    }
    
    
    // Inserting parabolic elements
    parabolic->LoadReferences();
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        fe_p_cindexes[iel].second.second = p_cel->Index();
        
    }
    
    // Block size structue including (Omega and Gamma)
    std::pair<int64_t, std::pair< TPZVec<int64_t>, TPZVec<int64_t> > > chunk_intp_indexes;
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi(n_el);
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_grad_phi(n_el);
    
    // Step 2 :: filling linking vectors
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(fe_p_cindexes[iel].second.first);
        TPZCompEl * p_cel = parabolic->Element(fe_p_cindexes[iel].second.second);
        
#ifdef PZDEBUG
        if (!e_cel || !p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel || !mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        int gel_dim = gel->Dimension();
        
        // Getting local integration index
        TPZManVector<int64_t> e_int_point_indexes(0,0);
        TPZManVector<int64_t> p_int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);
        
        chunk_intp_indexes.first = gel->Index();
        chunk_intp_indexes.second.first   = e_int_point_indexes;
        chunk_intp_indexes.second.second  = e_int_point_indexes;
        fe_p_intp_indexes.push_back(chunk_intp_indexes);
        
        this->ElementDofIndexes(mf_e_cel, dof_indexes);
        fu_p_dof_scatter[iel] = dof_indexes;
        
        
        blocks_dimensions_phi[iel].first = p_int_point_indexes.size()*dim;
        blocks_dimensions_phi[iel].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[iel].first = p_int_point_indexes.size()*dim*gel_dim;
        blocks_dimensions_grad_phi[iel].second = dof_indexes.size();
    }
    
    // Step 3 :: Initialize the matrix
    fu_To_parabolic.Initialize(blocks_dimensions_phi);
    fgrad_u_To_parabolic.Initialize(blocks_dimensions_grad_phi);
    
    
    TPZManVector<int64_t> p_int_point_indexes(0,0);
    TPZManVector<int64_t> dof_indexes(0,0);
    
    // Step 4 :: Filling the matrix
    std::pair<int64_t, int64_t> block_dim;
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(fe_p_cindexes[iel].second.first);
        TPZCompEl * p_cel = parabolic->Element(fe_p_cindexes[iel].second.second);
        
#ifdef PZDEBUG
        if (!e_cel || !p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel || !mf_p_cel)
        {
            DebugStop();
        }
#endif
        
        
        TPZInterpolationSpace * e_intel = dynamic_cast<TPZInterpolationSpace * >(mf_e_cel->Element(0));
        
        // Getting local integration index
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);
        dof_indexes = fu_p_dof_scatter[iel];
        
        block_dim.first = p_int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_geomechanic = mf_p_cel->GetIntegrationRule();
        int np_cel = int_points_geomechanic.NPoints();
        
#ifdef PZDEBUG
        if (p_int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi, block_grad_phi;
        
        block_phi.Resize(block_dim.first*dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*dim*gel_dim,block_dim.second);
        
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
            
            if (!fSimulationData->ReducedBasisResolution().first){
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
//                        LOGPZ_ERROR(logger,sout.str());
                }
                
                TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            }
            
#ifdef PZDEBUG
            
            if(block_dim.second != phi.Rows() * dim){
                if (fSimulationData->ReducedBasisResolution().first) {
                    if (block_dim.second != phi.Rows()) {
                        DebugStop();
                    }
                }
                else{
                    DebugStop();
                }
                
            }
#endif
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < dim; id++) {
                    if (!fSimulationData->ReducedBasisResolution().first){
                        block_phi(ip*dim+id,jp*dim+id) = phi(jp,0);
                    }
                    else{
                        block_phi(ip*dim+id,jp) = phi(jp,id);
                    }
                }
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                int d_count = 0;
                for (int id = 0; id < dim; id++) {
                    for (int jd = 0; jd < gel_dim; jd++) {
                        if(gel_dim == dim){
                            
                            if (!fSimulationData->ReducedBasisResolution().first){
                                block_grad_phi(ip*dim*gel_dim + id*gel_dim + jd,jp*dim+id) = dphidx(jd,jp);
                            }
                            else{
                                block_grad_phi(ip*dim*gel_dim + id*gel_dim + jd,jp) = dphi(jp,d_count);
                                d_count++;
                            }
                        }
                        else{
                            block_grad_phi(ip*dim*gel_dim+id*gel_dim + jd,jp) = 0.0;
                        }
                    }
                }
            }
        }
        
        fu_To_parabolic.SetBlock(iel, block_phi);
        fgrad_u_To_parabolic.SetBlock(iel, block_grad_phi);
        
    }
    
//    fu_To_parabolic.Print(" u_to_p ");
//    fgrad_u_To_parabolic.Print(" grad_u_to_p ");
    
    return;
    
    
}

void TRMBuildTransfers::elliptic_To_parabolic(TPZCompMesh * elliptic, TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!elliptic || !parabolic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_u(fu_To_parabolic.Cols(),1,0.0);
    int n = fe_p_cindexes.size();
    int64_t pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fu_p_dof_scatter[i].size(); iequ++) {
            Scatter_u(pos,0) = elliptic->Solution()(fu_p_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> u_at_parabolic,grad_u_at_parabolic;
    fu_To_parabolic.Multiply(Scatter_u,u_at_parabolic);
    fgrad_u_To_parabolic.Multiply(Scatter_u, grad_u_at_parabolic);
    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = elliptic->Dimension();
    int n_el = fe_p_cindexes.size();
    
    int64_t first_point_phi = 0;
    int64_t first_point_dphi = 0;
    std::pair<int64_t, int64_t> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fe_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(fe_p_cindexes[iel].second.first);
        TPZCompEl * p_cel = parabolic->Element(fe_p_cindexes[iel].second.second);
        
#ifdef PZDEBUG
        if (!e_cel || !p_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_e_cel = dynamic_cast<TPZMultiphysicsElement * >(e_cel);
        TPZMultiphysicsElement * mf_p_cel = dynamic_cast<TPZMultiphysicsElement * >(p_cel);
        
#ifdef PZDEBUG
        if(!mf_e_cel || !mf_e_cel)
        {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fu_To_parabolic.GetSizeofBlock(iel);
        b_size_dphi = fgrad_u_To_parabolic.GetSizeofBlock(iel);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMaterial * material = parabolic->FindMaterial(matd_id);
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            
            TPZFNMatrix<3,STATE> u(1,3,0.0);
            TPZFNMatrix<9,STATE> grad_u(3,3,0.0);
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    u(0,id) = u_at_parabolic(first_point_phi + ip*dim + id,0);
                }
                
                for (int id = 0; id < dim ; id++) {
                    for (int jd = 0; jd < dim ; jd++) {
                        if(fSimulationData->ReducedBasisResolution().first){
                            grad_u(id,jd) = grad_u_at_parabolic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                        }
                        else{
                            grad_u(jd,id) = grad_u_at_parabolic(first_point_dphi + ip*dim*dim + id*dim + jd,0);
                        }
                    }
                }
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_grad_u_0(grad_u);
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
            
            DebugStop(); // Just volumetric coupling
            
        }
        
    }
    
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segregated Transfer methods (Gamma and Omega) :: Build methods Parabolic
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void TRMBuildTransfers::Build_parabolic_To_parabolic(TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    parabolic->LoadReferences();
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    
    fp_p_cindexes.resize(0);
    std::pair< int64_t, int64_t > chunk_geo_cel_indexes;
    
    // Step 1 :: Counting for valid elements (apply all the needed filters in this step)
    for (int64_t i = 0; i < geometry->NElements(); i++) {
        
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        int parabolic_res = fSimulationData->MHMResolution().second.second;
        if (gel->Level() != parabolic_res) {
            continue;
        }
        
        
        int mat_id = gel->MaterialId();
        if ( (dim == 2 && mat_id > 11) || (dim == 3 && mat_id > 13) ) { // Filtering bc reservoir elements
            continue;
        }
        
        if (mat_id == fSimulationData->Skeleton_material_Id() || mat_id == fSimulationData->InterfacesMatId()) { // Filtering skeleton reservoir elements
            continue;
        }
        
        if ((gel->Dimension() == dim-1) && gel->NumInterfaces() !=0) { // Filtering interface reservoir elements for transport
            continue;
        }
        
        chunk_geo_cel_indexes.first = gel->Index();
        chunk_geo_cel_indexes.second = -1;
        fp_p_cindexes.push_back(chunk_geo_cel_indexes);

    }

    
    int64_t n_el = fp_p_cindexes.size();
    fq_dof_scatter.resize(n_el);
    fp_dof_scatter.resize(n_el);
    fp_p_intp_indexes.resize(n_el);
    
    std::pair<int64_t, TPZVec<int64_t>  > chunk_intp_indexes;
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi_q(n_el);
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_div_phi_q(n_el);
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi_p(n_el);
    
    int q_index = 0;
    int p_index = 1;
    
    int q_points = 0;
    int div_q_points = 0;
    int p_points = 0;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
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
        
        // Geometry and cel link
        fp_p_cindexes[iel].second = p_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> p_int_point_indexes(0,0);
        
        TPZManVector<int64_t> q_dof_indexes(0,0);
        TPZManVector<int64_t> p_dof_indexes(0,0);
        
        int vec_dim = dim;
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);
        q_points        = p_int_point_indexes.size();
        div_q_points    = p_int_point_indexes.size();
        p_points        = p_int_point_indexes.size();
        
        this->ElementDofIndexes(mf_p_cel, q_dof_indexes,q_index);

        if (gel->Dimension() == dim) {
            this->ElementDofIndexes(mf_p_cel, p_dof_indexes,p_index);

        }
        else{
            div_q_points = 0;
            p_points     = 0;
            vec_dim      = 1;
        }

        fq_dof_scatter[iel] = q_dof_indexes;
        fp_dof_scatter[iel] = p_dof_indexes;
        
        blocks_dimensions_phi_q[iel].first = q_points*vec_dim;
        blocks_dimensions_phi_q[iel].second = q_dof_indexes.size();
        
        blocks_dimensions_div_phi_q[iel].first = div_q_points;
        blocks_dimensions_div_phi_q[iel].second = q_dof_indexes.size();
        
        blocks_dimensions_phi_p[iel].first = p_points;
        blocks_dimensions_phi_p[iel].second = p_dof_indexes.size();
    
        
    }
    
    // Initialize the matrix
    fp_To_parabolic.Initialize(blocks_dimensions_phi_p);
    fq_To_parabolic.Initialize(blocks_dimensions_phi_q);
    fdiv_q_To_parabolic.Initialize(blocks_dimensions_div_phi_q);

    
    TPZManVector<int64_t> p_int_point_indexes(0,0);
    
    std::pair<int64_t, int64_t> q_block_dim;
    std::pair<int64_t, int64_t> div_q_block_dim;
    std::pair<int64_t, int64_t> p_block_dim;
    
    // for velocity functions
    TPZMaterialData data;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_p_cindexes[iel].second);
        
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
        mf_p_cel->GetMemoryIndices(p_int_point_indexes);

        q_block_dim     = fq_To_parabolic.GetSizeofBlock(iel);
        div_q_block_dim = fdiv_q_To_parabolic.GetSizeofBlock(iel);
        p_block_dim     = fp_To_parabolic.GetSizeofBlock(iel);
        
        
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
        TPZFMatrix<double> block_phi_q, block_div_phi_q, block_phi_p;
        
        block_phi_q.Resize(q_block_dim.first,q_block_dim.second);
        block_div_phi_q.Resize(div_q_block_dim.first, div_q_block_dim.second);
        block_phi_p.Resize(p_block_dim.first,p_block_dim.second);
        
        // Velocity functions
        TPZInterpolationSpace * q_intel = dynamic_cast<TPZInterpolationSpace * >(mf_p_cel->Element(0));
        if(q_intel)
        {
            
            // Computing over all integration points of the compuational element cel
            TPZFNMatrix<100,REAL> phi(q_intel->NShapeF(),1,0.0);
            int el_dim = gel->Dimension();
            TPZFNMatrix<300,REAL> dphidxi(el_dim,q_intel->NShapeF(),0.0);

            for (int ip = 0; ip < np_cel; ip++)
            {
                TPZManVector<REAL,3> qsi(el_dim,0.0);
                STATE w;
                int_points.Point(ip, qsi, w);
                // Get the vectorial phi
                q_intel->Shape(qsi, phi, dphidxi);
                q_intel->InitMaterialData(data);
                q_intel->ComputeRequiredData(data,qsi);
                
                TPZFMatrix<STATE> dphi       = data.dphi;
                REAL JacobianDet = data.detjac;
                
                TPZFMatrix<STATE> Qaxes = data.axes;
                TPZFMatrix<STATE> QaxesT;
                TPZFMatrix<STATE> Jacobian = data.jacobian;
                TPZFMatrix<STATE> JacobianInverse = data.jacinv;
                
                TPZFMatrix<STATE> GradOfX;
                TPZFMatrix<STATE> GradOfXInverse;
                TPZFMatrix<STATE> VectorOnMaster;
                TPZFMatrix<STATE> VectorOnXYZ(3,1,0.0);
                Qaxes.Transpose(&QaxesT);
                QaxesT.Multiply(Jacobian, GradOfX);
                JacobianInverse.Multiply(Qaxes, GradOfXInverse);

                if (el_dim == dim) {
                    for (int jp = 0; jp < q_block_dim.second; jp++) {
                        int vector_index = data.fVecShapeIndex[jp].first;
                        int shape_index = data.fVecShapeIndex[jp].second;
                        
                        for (int k = 0; k < dim; k++) {
                            VectorOnXYZ(k,0) = data.fNormalVec(k,vector_index);
                        }

                        GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
                        VectorOnMaster *= JacobianDet;

                        for (int id = 0; id < dim; id++) {
                            block_phi_q(ip*dim+id,jp) = phi(shape_index,0)*VectorOnXYZ(id,0);
                            block_div_phi_q(ip,jp) += dphi(id,shape_index)*VectorOnMaster(id,0)/JacobianDet;
                        }
                    }
                }
                else{
                    for (int jp = 0; jp < q_block_dim.second; jp++) {
                        block_phi_q(ip,jp) = phi(jp,0);
                    }
                }
            }
            
        }
        
        
        // Pressure functions
        TPZInterpolationSpace * p_intel = dynamic_cast<TPZInterpolationSpace * >(mf_p_cel->Element(1));
        
        if(p_intel)
        {
            // for derivatives in real space
            int nshape = p_intel->NShapeF();
            TPZFNMatrix<220> phi(nshape,1);
            TPZFNMatrix<660> dphi(gel_dim,nshape);
            
            for (int ip = 0; ip < np_cel ; ip++)
            {
                TPZManVector<REAL,3> qsi(gel_dim,0.0);
                STATE w;
                int_points.Point(ip, qsi, w);
                p_intel->Shape(qsi, phi, dphi);
                
#ifdef PZDEBUG
                if(p_block_dim.second != phi.Rows()){
                    DebugStop();
                }
#endif
                for (int jp = 0; jp < phi.Rows(); jp++) {
                    block_phi_p(ip,jp) = phi(jp,0);
                }
                
            }
        }
        
        fq_To_parabolic.SetBlock(iel, block_phi_q);
        fdiv_q_To_parabolic.SetBlock(iel, block_div_phi_q);
        fp_To_parabolic.SetBlock(iel, block_phi_p);
        
    }

    
//    fq_To_parabolic.Print(" q_to_p ");
//    fdiv_q_To_parabolic.Print(" div_q_to_p ");
//    fp_To_parabolic.Print(" p_to_p ");
    
    return;

}

void TRMBuildTransfers::space_To_parabolic(TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!parabolic || fq_To_parabolic.Rows() == 0 || fdiv_q_To_parabolic.Rows() == 0 || fp_To_parabolic.Rows() == 0) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = parabolic->Dimension();
    int n_el = fp_p_cindexes.size();
    
    std::pair<int64_t, int64_t> b_size_phi_q, b_size_div_phi_q, b_size_phi_p;

    b_size_phi_q.first = 0;
    b_size_phi_q.second = 0;
    
    b_size_div_phi_q.first = 0;
    b_size_div_phi_q.second = 0;
    
    b_size_phi_p.first = 0;
    b_size_phi_p.second = 0;
    
    TPZFMatrix<STATE> block_phi_q, block_div_phi_q, block_phi_p;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_p_cindexes[iel].second);
        
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
        
        b_size_phi_q        = fq_To_parabolic.GetSizeofBlock(iel);
        b_size_div_phi_q    = fdiv_q_To_parabolic.GetSizeofBlock(iel);
        b_size_phi_p        = fp_To_parabolic.GetSizeofBlock(iel);

        fq_To_parabolic.GetBlock(iel, block_phi_q);
        fdiv_q_To_parabolic.GetBlock(iel, block_div_phi_q);
        fp_To_parabolic.GetBlock(iel, block_phi_p);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = parabolic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            int n_phi_q = b_size_phi_q.second;
            int n_phi_p = b_size_phi_p.second;
            
            TPZFNMatrix<3,STATE> phi_q(n_phi_q,dim,0.0);
            TPZFNMatrix<3,STATE> div_phi_q(n_phi_q,1,0.0);
            TPZFNMatrix<3,STATE> phi_p(n_phi_p,1,0.0);
            
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < n_phi_q; is++) {
                    for (int d = 0; d < dim; d++) {
                        phi_q(is,d) = block_phi_q(ip*dim+d,is);
                    }
                    div_phi_q(is,0) = block_div_phi_q(ip,is);
                }
                
                for (int is = 0; is < n_phi_p; is++) {
                    phi_p(is,0) = block_phi_p(ip,is);
                }
                
                associated_material->GetMemory()[ipos].Set_phi_q(phi_q);
                associated_material->GetMemory()[ipos].Set_div_phi_q(div_phi_q);
                associated_material->GetMemory()[ipos].Set_phi_p(phi_p);
                
            }
            
        }
        else{
            
            TPZMatWithMem<TRMMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            int n_phi_q = b_size_phi_q.second;
            
            TPZFNMatrix<3,STATE> phi_q(n_phi_q,1,0.0);
            
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int is = 0; is < n_phi_q; is++) {
                    phi_q(is,0) = block_phi_q(ip,is);
                }
            
                associated_material->GetMemory()[ipos].Set_phi_q(phi_q);
                
            }
            
        }
        
        
    }
    
}

void TRMBuildTransfers::spatial_props_To_parabolic(TPZCompMesh * parabolic){
    
    
#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    int dim = parabolic->Dimension();
    
    
    TPZFMatrix<STATE> kappa, kappa_inv;
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;
    TPZManVector<STATE, 10> lambda,mu,S_e,alpha;
    
    // Step one
    int n_elements = parabolic->NElements();
    TPZManVector<int64_t, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = parabolic->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = parabolic->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);

            int map_model = fSimulationData->Map()->MapModel();
            fSimulationData->Map()->Kappa(x, kappa, kappa_inv, vars);
            fSimulationData->Map()->phi(x, porosity, vars);
            fSimulationData->Map()->lambda(x, lambda, vars);
            fSimulationData->Map()->mu(x, mu, vars);
            fSimulationData->Map()->S_e(x, S_e, vars);
            fSimulationData->Map()->alpha(x, alpha, vars);
            
            if (map_model != 0) {
                lambda[0]   *= cos(porosity[0]);
                mu[0]       *= cos(porosity[0]);
                S_e[0]      *= sin(porosity[0]);
                alpha[0]    *= sin(porosity[0]);
            }
            
            associated_material->GetMemory()[indexes[ip]].Set_K_0(kappa);
            associated_material->GetMemory()[indexes[ip]].Set_Kinv_0(kappa_inv);
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
            associated_material->GetMemory()[indexes[ip]].Set_lambda(lambda[0]);
            associated_material->GetMemory()[indexes[ip]].Set_mu(mu[0]);
            associated_material->GetMemory()[indexes[ip]].Set_S_e(S_e[0]);
            associated_material->GetMemory()[indexes[ip]].Set_alpha(alpha[0]);
        }
    }
    
}


void TRMBuildTransfers::parabolic_To_parabolic(TPZCompMesh * parabolic){

#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_q(fq_To_parabolic.Cols(),1,0.0);
    TPZFMatrix<STATE> Scatter_p(fp_To_parabolic.Cols(),1,0.0);
    
    int n = fp_p_cindexes.size();
    int64_t pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fq_dof_scatter[i].size(); iequ++) {
            Scatter_q(pos,0) = parabolic->Solution()(fq_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fp_dof_scatter[i].size(); iequ++) {
            Scatter_p(pos,0) = parabolic->Solution()(fp_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> q_at_parabolic,div_q_at_parabolic,p_at_parabolic;
    fq_To_parabolic.Multiply(Scatter_q,q_at_parabolic);
    fdiv_q_To_parabolic.Multiply(Scatter_q,div_q_at_parabolic);
    fp_To_parabolic.Multiply(Scatter_p,p_at_parabolic);

    
    
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = parabolic->Dimension();
    int n_el = fp_p_cindexes.size();
    
    int64_t first_point_phi_q = 0;
    int64_t first_point_div_phi_q = 0;
    int64_t first_point_phi_p = 0;
    
    std::pair<int64_t, int64_t> b_size_phi_q, b_size_div_phi_q, b_size_phi_p;
    
    b_size_phi_q.first = 0;
    b_size_phi_q.second = 0;
    
    b_size_div_phi_q.first = 0;
    b_size_div_phi_q.second = 0;
    
    b_size_phi_p.first = 0;
    b_size_phi_p.second = 0;
    
    TPZFMatrix<STATE> block_phi_q, block_div_phi_q, block_phi_p;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_p_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_p_cindexes[iel].second);
        
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
        
        first_point_phi_q     += b_size_phi_q.first;
        first_point_div_phi_q += b_size_div_phi_q.first;
        first_point_phi_p     += b_size_phi_p.first;
        
        b_size_phi_q        = fq_To_parabolic.GetSizeofBlock(iel);
        b_size_div_phi_q    = fdiv_q_To_parabolic.GetSizeofBlock(iel);
        b_size_phi_p        = fp_To_parabolic.GetSizeofBlock(iel);
        
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = parabolic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
        
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            
            TPZManVector<REAL,3> q(3,0.0);
            STATE div_q, p;
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                for (int id = 0; id < dim ; id++) {
                    q[id] = q_at_parabolic(first_point_phi_q + ip*dim + id,0);
                }
                
                div_q   = div_q_at_parabolic(first_point_div_phi_q + ip,0);
                p       = p_at_parabolic(first_point_phi_p + ip,0);
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_p_0(p);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_q_n(q);
                    associated_material->GetMemory()[ipos].Set_div_q_n(div_q);
                    associated_material->GetMemory()[ipos].Set_p_n(p);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_q(q);
                    associated_material->GetMemory()[ipos].Set_div_q(div_q);
                    associated_material->GetMemory()[ipos].Set_p(p);
                }
                
            }
            
            
        }
        else{
            

            TPZMatWithMem<TRMMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            
            TPZManVector<REAL,3> q(1,0.0);
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                q[0]       = q_at_parabolic(first_point_phi_q + ip,0);
                
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_q_n(q);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_q(q);
                }
                
            }
            
        }
        
    }


}

void TRMBuildTransfers::Build_parabolic_To_elliptic(TPZCompMesh * parabolic, TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!parabolic || !elliptic) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    parabolic->LoadReferences();
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    
    fp_e_cindexes.resize(0);
    std::pair<int64_t, std::pair<int64_t, int64_t> > chunk_geo_cel_indexes;
    
    // Step 1 :: Counting for valid elements (apply all the needed filters in this step)
    for (int64_t i = 0; i < geometry->NElements(); i++) {
        
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        int parabolic_res = fSimulationData->MHMResolution().second.second;
        if (gel->Level() != parabolic_res) {
            continue;
        }
        
        
        int mat_id = gel->MaterialId();
        if ( (gel->Dimension() == dim-1) || ( (dim == 2 && mat_id > 11) || (dim == 3 && mat_id > 13) ) ) { // Filtering bc side burden elements
            continue;
        }
        
        if (mat_id == fSimulationData->Skeleton_material_Id() || mat_id == fSimulationData->InterfacesMatId()) { // Filtering skeleton reservoir elements
            continue;
        }
        
        if ((gel->Dimension() == dim-1) && gel->NumInterfaces() !=0) { // Filtering interface reservoir elements for transport
            continue;
        }
        
        chunk_geo_cel_indexes.first = gel->Index();
        chunk_geo_cel_indexes.second.first  = -1;
        chunk_geo_cel_indexes.second.second = -1;
        fp_e_cindexes.push_back(chunk_geo_cel_indexes);
        
    }
    
    
    int64_t n_el = fp_e_cindexes.size();
    fp_e_dof_scatter.Resize(n_el);
    fp_e_intp_indexes.resize(n_el);
    
    
    // Inserting parabolic
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!p_cel) {
            DebugStop();
        }
#endif
        
        // Geometry and cel link
        fp_e_cindexes[iel].second.first = p_cel->Index();
    }
    
    // Inserting elliptic
    elliptic->LoadReferences();
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        // Geometry and cel link
        fp_e_cindexes[iel].second.second = e_cel->Index();
    }
    
    
    std::pair<int64_t, std::pair< TPZVec<int64_t>, TPZVec<int64_t> > > chunk_intp_indexes;
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi_p(n_el);
    
    int p_index = 1;
    int p_points = 0;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_e_cindexes[iel].second.first);
        TPZCompEl * e_cel = elliptic->Element(fp_e_cindexes[iel].second.second);
        
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
        TPZManVector<int64_t> e_int_point_indexes(0,0);
        TPZManVector<int64_t> p_dof_indexes(0,0);
        
        int vec_dim = dim;
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        p_points        = e_int_point_indexes.size();
        
        if (gel->Dimension() == dim) {
            this->ElementDofIndexes(mf_p_cel, p_dof_indexes, p_index);
        }
        else{
            p_points     = 0;
            vec_dim      = 1;
        }
        
        fp_e_dof_scatter[iel] = p_dof_indexes;
        
        blocks_dimensions_phi_p[iel].first = p_points;
        blocks_dimensions_phi_p[iel].second = p_dof_indexes.size();
        
        
    }
    
    // Initialize the matrix
    fp_To_elliptic.Initialize(blocks_dimensions_phi_p);
    
    
    TPZManVector<int64_t> e_int_point_indexes(0,0);
    std::pair<int64_t, int64_t> p_block_dim;
    
    // for velocity functions
    TPZMaterialData data;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_e_cindexes[iel].second.first);
        TPZCompEl * e_cel = elliptic->Element(fp_e_cindexes[iel].second.second);
        
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
        mf_e_cel->GetMemoryIndices(e_int_point_indexes);
        p_block_dim     = fp_To_elliptic.GetSizeofBlock(iel);
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points = mf_e_cel->GetIntegrationRule();
        int np_cel = int_points.NPoints();
        
#ifdef PZDEBUG
        if (e_int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi_p;
        block_phi_p.Resize(p_block_dim.first,p_block_dim.second);
        
        
        // Pressure functions
        TPZInterpolationSpace * p_intel = dynamic_cast<TPZInterpolationSpace * >(mf_p_cel->Element(1));
        
        if(p_intel)
        {
            // for derivatives in real space
            int nshape = p_intel->NShapeF();
            TPZFNMatrix<220> phi(nshape,1);
            TPZFNMatrix<660> dphi(gel_dim,nshape);
            
            for (int ip = 0; ip < np_cel ; ip++)
            {
                TPZManVector<REAL,3> qsi(gel_dim,0.0);
                STATE w;
                int_points.Point(ip, qsi, w);
                p_intel->Shape(qsi, phi, dphi);
                
#ifdef PZDEBUG
                if(p_block_dim.second != phi.Rows()){
                    DebugStop();
                }
#endif
                for (int jp = 0; jp < phi.Rows(); jp++) {
                    block_phi_p(ip,jp) = phi(jp,0);
                }
                
            }
        }
        
        fp_To_elliptic.SetBlock(iel, block_phi_p);
        
    }
    
//    fp_To_elliptic.Print(" p_to_e ");
    
    return;
    
}

void TRMBuildTransfers::parabolic_To_elliptic(TPZCompMesh * parabolic, TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!parabolic || !elliptic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_p(fp_To_elliptic.Cols(),1,0.0);
    
    int n = fp_e_cindexes.size();
    int64_t pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fp_e_dof_scatter[i].size(); iequ++) {
            Scatter_p(pos,0) = parabolic->Solution()(fp_e_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> p_at_elliptic;
    fp_To_elliptic.Multiply(Scatter_p,p_at_elliptic);

    
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = parabolic->Dimension();
    int n_el = fp_e_cindexes.size();
    
    int64_t first_point_phi_p = 0;
    
    std::pair<int64_t, int64_t> b_size_phi_p;
    
    b_size_phi_p.first = 0;
    b_size_phi_p.second = 0;
    
    TPZFMatrix<STATE> block_phi_p;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fp_e_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fp_e_cindexes[iel].second.first);
        TPZCompEl * e_cel = elliptic->Element(fp_e_cindexes[iel].second.second);
        
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
        
        first_point_phi_p     += b_size_phi_p.first;
        b_size_phi_p        = fp_To_elliptic.GetSizeofBlock(iel);
        
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = elliptic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            STATE p;
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                p       = p_at_elliptic(first_point_phi_p + ip,0);
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_p_0(p);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_p_n(p);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_p(p);
                }
                
            }
            
            
        }
        else{
            DebugStop(); // Just volumetric coupling
        }
        
    }
    
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segregated Transfer methods (Gamma and Omega) :: Build methods Hyperbolic
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void TRMBuildTransfers::Build_hyperbolic_To_hyperbolic(TPZCompMesh * hyperbolic){
    
#ifdef PZDEBUG
    if (!hyperbolic) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    hyperbolic->LoadReferences();
    TPZGeoMesh * geometry = hyperbolic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    
    fh_h_cindexes.resize(0);
    std::pair< int64_t, int64_t > chunk_geo_cel_indexes;
    
    // Step 1 :: Counting for valid elements (apply all the needed filters in this step)
    for (int64_t i = 0; i < geometry->NElements(); i++) {
        
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        
        
        int mat_id = gel->MaterialId();
        if ( (dim == 2 && mat_id > 11) || (dim == 3 && mat_id > 13) ) { // Filtering bc sideburden elements
            continue;
        }
        
        if ((gel->Dimension() == dim-1) && (dim == 3 && (mat_id >= 8 && mat_id <= 13))) { // Filtering bc reservoir elements
            continue;
        }
        
        if (mat_id == fSimulationData->Skeleton_material_Id() || mat_id == fSimulationData->InterfacesMatId()) { // Filtering skeleton reservoir elements
            continue;
        }
        
        if ((gel->Dimension() == dim-1)) { // Filtering all dim-1 elements for transport
            continue;
        }
        
        
        
        chunk_geo_cel_indexes.first = gel->Index();
        chunk_geo_cel_indexes.second = -1;
        fh_h_cindexes.push_back(chunk_geo_cel_indexes);
        
    }
    
    
    int64_t n_el = fh_h_cindexes.size();
    fsw_dof_scatter.resize(n_el);
    
    std::pair<int64_t, TPZVec<int64_t>  > chunk_intp_indexes;
    
    // Block size structue including (Omega)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi_sw(n_el);
    
    int sw_index = 0;
    int sw_points = 0;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fh_h_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * h_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!h_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
        
#ifdef PZDEBUG
        if(!mf_h_cel)
        {
            DebugStop();
        }
#endif
        
        // Geometry and cel link
        fh_h_cindexes[iel].second = h_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> sw_int_point_indexes(0,0);
        TPZManVector<int64_t> sw_dof_indexes(0,0);
        
        mf_h_cel->GetMemoryIndices(sw_int_point_indexes);
        sw_points        = sw_int_point_indexes.size();
        
        this->ElementDofIndexes(mf_h_cel, sw_dof_indexes,sw_index);
        
        if (gel->Dimension() != dim) {
            DebugStop();
            
        }
        
        fsw_dof_scatter[iel] = sw_dof_indexes;
        
        blocks_dimensions_phi_sw[iel].first = sw_points;
        blocks_dimensions_phi_sw[iel].second = sw_dof_indexes.size();
        
        
    }
    
    // Initialize the matrix
    fsw_To_hyperbolic.Initialize(blocks_dimensions_phi_sw);
    
    
    TPZManVector<int64_t> sw_int_point_indexes(0,0);
    std::pair<int64_t, int64_t> sw_block_dim;
    
    // for velocity functions
    TPZMaterialData data;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fh_h_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * h_cel = hyperbolic->Element(fh_h_cindexes[iel].second);
        
#ifdef PZDEBUG
        if (!h_cel) {
            DebugStop();
        }
#endif
        
        
        TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
        
#ifdef PZDEBUG
        if(!mf_h_cel)
        {
            DebugStop();
        }
#endif
        
        
        // Getting local integration index
        mf_h_cel->GetMemoryIndices(sw_int_point_indexes);
        
        sw_block_dim     = fsw_To_hyperbolic.GetSizeofBlock(iel);
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points = mf_h_cel->GetIntegrationRule();
        int np_cel = int_points.NPoints();
        
#ifdef PZDEBUG
        if (sw_int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi_sw;
        
        block_phi_sw.Resize(sw_block_dim.first,sw_block_dim.second);
        
        
        // Water Saturations functions
        TPZInterpolationSpace * h_intel = dynamic_cast<TPZInterpolationSpace * >(mf_h_cel->Element(0));
        
        if(h_intel)
        {
            // for derivatives in real space
            int nshape = h_intel->NShapeF();
            TPZFNMatrix<220> phi(nshape,1);
            TPZFNMatrix<660> dphi(gel_dim,nshape);
            
            for (int ip = 0; ip < np_cel ; ip++)
            {
                TPZManVector<REAL,3> qsi(gel_dim,0.0);
                STATE w;
                int_points.Point(ip, qsi, w);
                h_intel->Shape(qsi, phi, dphi);
                
#ifdef PZDEBUG
                if(sw_block_dim.second != phi.Rows()){
                    DebugStop();
                }
#endif
                for (int jp = 0; jp < phi.Rows(); jp++) {
                    block_phi_sw(ip,jp) = phi(jp,0);
                }
                
            }
        }
        
        fsw_To_hyperbolic.SetBlock(iel, block_phi_sw);
        
    }
    
    
//    fsw_To_hyperbolic.Print(" sw_to_h ");
    
    return;
    
}

void TRMBuildTransfers::spatial_props_To_hyperbolic(TPZCompMesh * hyperbolic){
    
#ifdef PZDEBUG
    if (!hyperbolic) {
        DebugStop();
    }
#endif
    
    int dim = hyperbolic->Dimension();
    
    
    TPZFMatrix<STATE> kappa, kappa_inv;
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;
    TPZManVector<STATE, 10> lambda,mu,S_e,alpha;
    
    // Step one
    int n_elements = hyperbolic->NElements();
    TPZManVector<int64_t, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = hyperbolic->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = hyperbolic->FindMaterial(rockid);
        TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);
            
            int map_model = fSimulationData->Map()->MapModel();
            fSimulationData->Map()->Kappa(x, kappa, kappa_inv, vars);
            fSimulationData->Map()->phi(x, porosity, vars);
            fSimulationData->Map()->lambda(x, lambda, vars);
            fSimulationData->Map()->mu(x, mu, vars);
            fSimulationData->Map()->S_e(x, S_e, vars);
            fSimulationData->Map()->alpha(x, alpha, vars);
            
            if (map_model != 0) {
                lambda[0]   *= cos(porosity[0]);
                mu[0]       *= cos(porosity[0]);
                S_e[0]      *= sin(porosity[0]);
                alpha[0]    *= sin(porosity[0]);
            }
            
            associated_material->GetMemory()[indexes[ip]].Set_K_0(kappa);
            associated_material->GetMemory()[indexes[ip]].Set_Kinv_0(kappa_inv);
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
            associated_material->GetMemory()[indexes[ip]].Set_lambda(lambda[0]);
            associated_material->GetMemory()[indexes[ip]].Set_mu(mu[0]);
            associated_material->GetMemory()[indexes[ip]].Set_S_e(S_e[0]);
            associated_material->GetMemory()[indexes[ip]].Set_alpha(alpha[0]);
        }
    }
    
}

void TRMBuildTransfers::Build_parabolic_hyperbolic_cel_pairs(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic){
    
    fparabolic_hyperbolic_cel_pairs.resize(0);
    
#ifdef PZDEBUG
    if (!parabolic) {
        DebugStop();
    }
    
    if (!fSimulationData->IsOnePhaseQ() && !hyperbolic) {
        DebugStop();
    }
    
    if (!fSimulationData->TransporResolution().first) {
        DebugStop();
    }
    
#endif
    
    parabolic->LoadReferences();
    TPZGeoMesh * geometry = parabolic->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<int64_t, std::pair <int64_t, std::vector<int64_t> > > gel_cel_indexes;
    
    for (int64_t igel = 0; igel < geometry->NElements(); igel++) {
        
        TPZGeoEl * gel = geometry->Element(igel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() != dimension) {
            continue;
        }
        int mat_id = gel->MaterialId();
        bool filter_check = (mat_id == 5) || (mat_id == 6) || (mat_id == 7);
        if (!filter_check) {
            continue;
        }
        
        int parabolic_res = fSimulationData->MHMResolution().second.second;
        if (gel->Level() != parabolic_res) {
            continue;
        }
        
        TPZCompEl * mixed_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!mixed_cel) {
            DebugStop();
        }
#endif
    
        gel_cel_indexes.first = gel->Index();
        gel_cel_indexes.second.first = mixed_cel->Index();
        gel_cel_indexes.second.second.resize(0);
        fparabolic_hyperbolic_cel_pairs.push_back(gel_cel_indexes);
        
    }
    
    // counting volumetric elements
    int nvol_elements = fparabolic_hyperbolic_cel_pairs.size();
    
    if(fSimulationData->IsOnePhaseQ()){
        
        //        for (int k = 0; k < fmixed_transport_comp_indexes.size(); k++) {
        //            std::cout << " volume k : " << k <<std::endl;
        //            std::cout << " volume gel : " << fmixed_transport_comp_indexes[k].first <<std::endl;
        //            std::cout << " volume cmixed : " << fmixed_transport_comp_indexes[k].second.first <<std::endl;
        //        }
        
        return;
    }
    
    // inserting transport indexes
    hyperbolic->LoadReferences();
    TPZVec<TPZGeoEl *> n_refined_sons;
    int cel_index;
    for(int64_t ivol = 0; ivol < nvol_elements; ivol++){
        
        TPZGeoEl * father_gel = geometry->Element(fparabolic_hyperbolic_cel_pairs[ivol].first);
        
#ifdef PZDEBUG
        if (!father_gel) {
            DebugStop();
        }
#endif
        n_refined_sons.resize(0);
        father_gel->GetHigherSubElements(n_refined_sons);
        
        if (n_refined_sons.size() == 0) {
            n_refined_sons.resize(1);
            n_refined_sons[0] = father_gel;
        }
        
        int n_sons = n_refined_sons.size();
        for (int igel = 0; igel < n_sons; igel++) {
            cel_index = geometry->Element(n_refined_sons[igel]->Index())->Reference()->Index();
            fparabolic_hyperbolic_cel_pairs[ivol].second.second.push_back(cel_index);
        }
        
    }
    
//    for (int k = 0; k < fparabolic_hyperbolic_cel_pairs.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << fparabolic_hyperbolic_cel_pairs[k].first <<std::endl;
//        std::cout << " volume cmixed : " << fparabolic_hyperbolic_cel_pairs[k].second.first <<std::endl;
//        int n_sons = n_refined_sons.size();
//        for (int igel = 0; igel < n_sons; igel++) {
//            int index = fparabolic_hyperbolic_cel_pairs[k].second.second[igel]; ;
//            std::cout << " volume ctransport : " << index <<std::endl;
//        }
//    }
    
}

void TRMBuildTransfers::Build_elliptic_hyperbolic_cel_pairs(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic){
    
    felliptic_hyperbolic_cel_pairs.resize(0);
    
#ifdef PZDEBUG
    if (!elliptic) {
        DebugStop();
    }
    
    if (!fSimulationData->IsOnePhaseQ() && !hyperbolic) {
        DebugStop();
    }
    
    if (!fSimulationData->TransporResolution().first) {
        DebugStop();
    }
    
#endif
    
    elliptic->LoadReferences();
    TPZGeoMesh * geometry = elliptic->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<int64_t, std::pair <int64_t, std::vector<int64_t> > > gel_cel_indexes;
    
    for (int64_t igel = 0; igel < geometry->NElements(); igel++) {
        
        TPZGeoEl * gel = geometry->Element(igel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() != dimension) {
            continue;
        }
        
        int mat_id = gel->MaterialId();
        bool filter_check = (mat_id == 5) || (mat_id == 6) || (mat_id == 7);
        if (!filter_check) {
            continue;
        }
        
        int parabolic_res = fSimulationData->MHMResolution().second.second;
        if (gel->Level() != parabolic_res) {
            continue;
        }
        
        TPZCompEl * e_cel = gel->Reference();
        
#ifdef PZDEBUG
        if (!e_cel) {
            DebugStop();
        }
#endif
        
        gel_cel_indexes.first = gel->Index();
        gel_cel_indexes.second.first = e_cel->Index();
        gel_cel_indexes.second.second.resize(0);
        felliptic_hyperbolic_cel_pairs.push_back(gel_cel_indexes);
        
    }
    
    // counting volumetric elements
    int nvol_elements = felliptic_hyperbolic_cel_pairs.size();
    
    if(fSimulationData->IsOnePhaseQ()){
        
//        for (int k = 0; k < felliptic_hyperbolic_cel_pairs.size(); k++) {
//            std::cout << " volume k : " << k <<std::endl;
//            std::cout << " volume gel : " << felliptic_hyperbolic_cel_pairs[k].first <<std::endl;
//            std::cout << " volume celliptic : " << felliptic_hyperbolic_cel_pairs[k].second.first <<std::endl;
//        }
        
        return;
    }
    
    // inserting transport indexes
    hyperbolic->LoadReferences();
    TPZVec<TPZGeoEl *> n_refined_sons;
    int cel_index;
    for(int64_t ivol = 0; ivol < nvol_elements; ivol++){
        
        TPZGeoEl * father_gel = geometry->Element(felliptic_hyperbolic_cel_pairs[ivol].first);
        
#ifdef PZDEBUG
        if (!father_gel) {
            DebugStop();
        }
#endif
        n_refined_sons.resize(0);
        father_gel->GetHigherSubElements(n_refined_sons);
        
        if (n_refined_sons.size() == 0) {
            n_refined_sons.resize(1);
            n_refined_sons[0] = father_gel;
        }
        
        int n_sons = n_refined_sons.size();
        for (int igel = 0; igel < n_sons; igel++) {
            cel_index = geometry->Element(n_refined_sons[igel]->Index())->Reference()->Index();
            felliptic_hyperbolic_cel_pairs[ivol].second.second.push_back(cel_index);
        }
        
    }
    
//    for (int k = 0; k < felliptic_hyperbolic_cel_pairs.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << felliptic_hyperbolic_cel_pairs[k].first <<std::endl;
//        std::cout << " volume celliptic : " << felliptic_hyperbolic_cel_pairs[k].second.first <<std::endl;
//        int n_sons = n_refined_sons.size();
//        for (int igel = 0; igel < n_sons; igel++) {
//            int index = felliptic_hyperbolic_cel_pairs[k].second.second[igel]; ;
//            std::cout << " volume ctransport : " << index <<std::endl;
//        }
//    }
    
}

void TRMBuildTransfers::hyperbolic_To_hyperbolic(TPZCompMesh * hyperbolic){
    
#ifdef PZDEBUG
    if (!hyperbolic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_sw(fsw_To_hyperbolic.Cols(),1,0.0);
    
    int n = fh_h_cindexes.size();
    int64_t pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fsw_dof_scatter[i].size(); iequ++) {
            Scatter_sw(pos,0) = hyperbolic->Solution()(fsw_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step one
    TPZFMatrix<STATE> sw_at_hyperbolic;
    fsw_To_hyperbolic.Multiply(Scatter_sw,sw_at_hyperbolic);
    
    
    
    TPZGeoMesh * geometry = hyperbolic->Reference();
    int dim = hyperbolic->Dimension();
    int n_el = fh_h_cindexes.size();
    
    int64_t first_point_phi_sw = 0;
    
    std::pair<int64_t, int64_t> b_size_phi_sw;
    
    b_size_phi_sw.first = 0;
    b_size_phi_sw.second = 0;
    
    TPZFMatrix<STATE> block_phi_sw;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fh_h_cindexes[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * h_cel = hyperbolic->Element(fh_h_cindexes[iel].second);
        
#ifdef PZDEBUG
        if (!h_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
        
#ifdef PZDEBUG
        if(!mf_h_cel)
        {
            DebugStop();
        }
#endif
        
        first_point_phi_sw     += b_size_phi_sw.first;
        
        b_size_phi_sw        = fsw_To_hyperbolic.GetSizeofBlock(iel);
        
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = hyperbolic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_h_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            
            TPZManVector<REAL,3> q(3,0.0);
            STATE sw;
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                sw       = sw_at_hyperbolic(first_point_phi_sw + ip,0);
                
                if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                    associated_material->GetMemory()[ipos].Set_sa_0(sw);
                }
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_sa_n(sw);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_sa(sw);
                }
                
            }
            
            
        }
        else{
            
            DebugStop();// Volumetric transfer
        }
    }

    
}

void TRMBuildTransfers::Build_parabolic_hyperbolic_volumetric(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic){
 
#ifdef PZDEBUG
    if (!parabolic || !hyperbolic) {
        DebugStop();
    }
#endif
    
#ifdef PZDEBUG
    if (fparabolic_hyperbolic_cel_pairs.size() == 0) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    int64_t n_el = fparabolic_hyperbolic_cel_pairs.size();
    fp_avg_dof_scatter.Resize(n_el);
    
    std::pair<int64_t, std::pair< TPZVec<int64_t>, TPZVec<int64_t> > > chunk_intp_indexes;
    
    // Block size structue for Omega
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi_p_avg(n_el);
    
    int p_avg_index = 1;
    int p_avg_points = 0;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fparabolic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.first);
        
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
        
        p_avg_points = fparabolic_hyperbolic_cel_pairs[iel].second.second.size();
        
        // Getting local integration index
        TPZManVector<int64_t> p_avg_dof_indexes(0,0);
        
        if (gel->Dimension() == dim) {
            this->ElementDofIndexes(mf_p_cel, p_avg_dof_indexes, p_avg_index);
        }
        else{
            DebugStop();
        }
        
        fp_avg_dof_scatter[iel] = p_avg_dof_indexes;
        
        blocks_dimensions_phi_p_avg[iel].first = p_avg_points;
        blocks_dimensions_phi_p_avg[iel].second = p_avg_dof_indexes.size();
        
        
    }
    
    // Initialize the matrix
    fp_avg_To_hyperbolic.Initialize(blocks_dimensions_phi_p_avg);
    
    
    TPZManVector<int64_t> sw_int_point_indexes(0,0);
    std::pair<int64_t, int64_t> p_avg_block_dim;
    
    // for velocity functions
    TPZMaterialData data;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fparabolic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.first);

        
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
        
        int n_sons = fparabolic_hyperbolic_cel_pairs[iel].second.second.size();
        p_avg_block_dim     = fp_avg_To_hyperbolic.GetSizeofBlock(iel);
        
        TPZFMatrix<double> block_phi_p_avg(p_avg_block_dim.first,p_avg_block_dim.second,0.0);
        
        for (int ison = 0; ison < n_sons ; ison++) {
            TPZCompEl * h_cel = hyperbolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            
            // Getting local integration index
            mf_h_cel->GetMemoryIndices(sw_int_point_indexes);
            
            // Pressure functions
            TPZInterpolationSpace * p_intel = dynamic_cast<TPZInterpolationSpace * >(mf_p_cel->Element(1));
            int gel_dim = sub_gel->Dimension();
            TPZManVector<REAL,3> xi_origin_vol(gel_dim),xi_target(gel_dim);
            REAL sub_volume = this->DimensionalMeasure(sub_gel);
            
            REAL detjac,w;
            TPZFMatrix<REAL> jac;
            TPZFMatrix<REAL> axes;
            TPZFMatrix<REAL> jacinv;
            
            if(p_intel)
            {
                // for derivatives in real space
                int nshape = p_intel->NShapeF();
                TPZFNMatrix<220> phi(nshape,1);
                TPZFNMatrix<660> dphi(gel_dim,nshape);
                
                int order = 4; // increase to 4 o 5
                TPZIntPoints * int_points = sub_gel->CreateSideIntegrationRule(sub_gel->NSides()-1, order);
                int n_points = int_points->NPoints();
                for (int ip = 0; ip < n_points; ip++) {

  
                    int_points->Point(ip, xi_origin_vol, w);
                    sub_gel->Jacobian(xi_origin_vol, jac, axes, detjac, jacinv);
                    
                    if (sub_gel->FatherIndex() == -1) {
                        xi_target = xi_origin_vol;
                    }
                    else{
                        sub_gel->TransformSonToFather(sub_gel->Father(), xi_origin_vol, xi_target);
                    }
        
                    p_intel->Shape(xi_target, phi, dphi);
                    
                    
                    for (int jp = 0; jp < phi.Rows(); jp++) {
                        block_phi_p_avg(ison,jp) += w * detjac * phi(jp,0)/sub_volume;
                    }
                    
                    
                }

            }
            
        }
        
        fp_avg_To_hyperbolic.SetBlock(iel, block_phi_p_avg);
        
    }
    
//    fp_avg_To_hyperbolic.Print(" p_avg_to_h ");
    
    return;
    
}

void TRMBuildTransfers::Build_elliptic_hyperbolic_volumetric(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic){
 
#ifdef PZDEBUG
    if (!elliptic || !hyperbolic) {
        DebugStop();
    }
#endif
    
#ifdef PZDEBUG
    if (felliptic_hyperbolic_cel_pairs.size() == 0) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    int64_t n_el = felliptic_hyperbolic_cel_pairs.size();
    fu_avg_dof_scatter.Resize(n_el);
    
    std::pair<int64_t, std::pair< TPZVec<int64_t>, TPZVec<int64_t> > > chunk_intp_indexes;
    
    // Block size structue for Omega
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi_p_avg(n_el);
    
    int u_avg_index = 0;
    int u_avg_points = 0;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(felliptic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(felliptic_hyperbolic_cel_pairs[iel].second.first);
        
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
        
        u_avg_points = felliptic_hyperbolic_cel_pairs[iel].second.second.size()*dim*dim;
        
        
        // Getting local integration index
        TPZManVector<int64_t> u_avg_dof_indexes(0,0);
        
        if (gel->Dimension() == dim) {
            this->ElementDofIndexes(mf_e_cel, u_avg_dof_indexes, u_avg_index);
        }
        else{
            DebugStop();
        }
        
        fu_avg_dof_scatter[iel] = u_avg_dof_indexes;
        
        blocks_dimensions_phi_p_avg[iel].first = u_avg_points;
        blocks_dimensions_phi_p_avg[iel].second = u_avg_dof_indexes.size();
        
        
    }
    
    // Initialize the matrix
    fgrad_u_avg_To_hyperbolic.Initialize(blocks_dimensions_phi_p_avg);
    
    
    TPZManVector<int64_t> sw_int_point_indexes(0,0);
    std::pair<int64_t, int64_t> gra_u_avg_block_dim;
    
    // for velocity functions
    TPZMaterialData data;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(felliptic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(felliptic_hyperbolic_cel_pairs[iel].second.first);
        
        
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
        
        int n_sons = felliptic_hyperbolic_cel_pairs[iel].second.second.size();
        gra_u_avg_block_dim     = fgrad_u_avg_To_hyperbolic.GetSizeofBlock(iel);
        
        TPZFMatrix<double> block_grad_phi_u_avg(gra_u_avg_block_dim.first,gra_u_avg_block_dim.second,0.0);
        
        for (int ison = 0; ison < n_sons ; ison++) {
            TPZCompEl * h_cel = hyperbolic->Element(felliptic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            
            // Getting local integration index
            mf_h_cel->GetMemoryIndices(sw_int_point_indexes);
            
            // Pressure functions
            TPZInterpolationSpace * e_intel = dynamic_cast<TPZInterpolationSpace * >(mf_e_cel->Element(0));
            int gel_dim = sub_gel->Dimension();
            TPZManVector<REAL,3> xi_origin_vol(gel_dim),xi_target(gel_dim);
            REAL sub_volume = this->DimensionalMeasure(sub_gel);
            
            REAL detjac_sub_cel,detjac,w;
            TPZFMatrix<REAL> jac;
            TPZFMatrix<REAL> axes;
            TPZFMatrix<REAL> jacinv;
            TPZFMatrix<REAL> jacobian;
            
            if(e_intel) // Compute gradient
            {
                
                // for derivatives in real space
                int nshape = e_intel->NShapeF();
                TPZFNMatrix<220> phi(nshape,1);
                TPZFNMatrix<660> dphi(gel_dim,nshape),dphix_axes(gel_dim,nshape);
                TPZFMatrix<double> dphidx;
                
                
                int order = 4; // increase to 4 o 5
                TPZIntPoints * int_points = sub_gel->CreateSideIntegrationRule(sub_gel->NSides()-1, order);
                int n_points = int_points->NPoints();
                for (int ip = 0; ip < n_points; ip++)
                {
                    
                    int_points->Point(ip, xi_origin_vol, w);
                    sub_gel->Jacobian(xi_origin_vol, jac, axes, detjac_sub_cel, jacinv);

                    if (sub_gel->FatherIndex() == -1) {
                        xi_target = xi_origin_vol;
                    }
                    else{
                        sub_gel->TransformSonToFather(sub_gel->Father(), xi_origin_vol, xi_target);
                    }
                    
                    // Get the phi and dphix for H1 elasticity
                    e_intel->Shape(xi_target, phi, dphi);
                    
                    if (!fSimulationData->ReducedBasisResolution().first){
                        gel->Jacobian( xi_origin_vol, jacobian, axes, detjac , jacinv);
                        
                        switch(gel->Dimension()) {
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
//                                LOGPZ_ERROR(logger,sout.str());
                        }
                        
                        TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
                    }
                    
                    dphidx *= w/sub_volume; // avegare gradient
                    
                    for (int jp = 0; jp < phi.Rows(); jp++) {
                        int d_count = 0;
                        for (int id = 0; id < dim; id++) {
                            for (int jd = 0; jd < gel_dim; jd++) {
                                if(gel_dim == dim){
                                    
                                    if (!fSimulationData->ReducedBasisResolution().first){
                                        block_grad_phi_u_avg(ison*dim*dim + id*gel_dim + jd,jp*dim+id) += dphidx(jd,jp);
                                    }
                                    else{
                                        block_grad_phi_u_avg(ison*dim*dim + id*gel_dim + jd,jp) += w*dphi(jp,d_count)/sub_volume;
                                        d_count++;
                                    }
                                }
                                else{
                                    block_grad_phi_u_avg(ison*dim*dim + ip*dim*gel_dim+id*gel_dim + jd,jp) = 0.0;
                                }
                            }
                        }
                    }
                    
                }
                
            }
            
        }
        
        fgrad_u_avg_To_hyperbolic.SetBlock(iel, block_grad_phi_u_avg);
        
    }
    
//    fgrad_u_avg_To_hyperbolic.Print(" grad_u_avg_to_h ");
    
    return;
    
}

void TRMBuildTransfers::Build_parabolic_hyperbolic_left_right_pairs(TPZCompMesh * hyperbolic){
 
    fleft_right_g_c_indexes_Gamma.resize(0);
    fleft_right_g_c_indexes_gamma.resize(0);
    
    fcelint_celh_celp_Gamma.resize(0);
    fcelint_celh_celp_gamma.resize(0);
    
#ifdef PZDEBUG
    if (!hyperbolic) {
        std::cout << "There is no computational transport mesh, transport_mesh = Null." << std::endl;
        DebugStop();
    }
#endif
    
    hyperbolic->LoadReferences();
    TPZGeoMesh * geometry = hyperbolic->Reference();
    int dim = geometry->Dimension();
    
    int64_t n_el = hyperbolic->NElements();
    
    std::pair< TPZVec<int64_t> , std::pair< TPZVec<int64_t>, TPZVec<int64_t> > >  chunk;
    
    for (int64_t icel = 0; icel < n_el; icel++) {
        
        TPZCompEl * h_cel = hyperbolic->Element(icel);
        
#ifdef PZDEBUG
        if (!h_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsInterfaceElement * interface = dynamic_cast<TPZMultiphysicsInterfaceElement * >(h_cel);
        
        if (!interface) {
            continue;
        }
        
        TPZCompEl * left_cel = interface->LeftElement();
        TPZCompEl * right_cel = interface->RightElement();
        
#ifdef PZDEBUG
        
        if(!left_cel || !right_cel){
            DebugStop();
        }
#endif
        
        if(interface->Reference()->HasSubElement()) {
            continue;
        }
        
//        if( dim == 2 && (left_cel->Reference()->MaterialId() == 12 || right_cel->Reference()->MaterialId() == 12)) {
//            continue;
//        }
//        
//        if( dim == 3 && (left_cel->Reference()->MaterialId() == 14 || right_cel->Reference()->MaterialId() == 14)) {
//            continue;
//        }
        
        chunk.first.Resize(2);
        chunk.second.first.Resize(2);
        chunk.second.second.Resize(2);
        
        chunk.first[0] = h_cel->Reference()->Index();
        chunk.first[1] = h_cel->Index();
        
        chunk.second.first[0] = left_cel->Reference()->Index();
        chunk.second.first[1] = left_cel->Index();
        
        chunk.second.second[0] = right_cel->Reference()->Index();
        chunk.second.second[1] = right_cel->Index();
        
//#ifdef PZDEBUG
//        
//        std::cout << "h:: face index  = " << chunk.first <<  std::endl;
//        std::cout << "h:: volume index  = " << chunk.second.first <<  std::endl;
//        std::cout << "p:: volume index  = " << chunk.second.second <<  std::endl;
//        
//#endif
        
        if(left_cel->Dimension() != dim ||  right_cel->Dimension() != dim){
            
            fleft_right_g_c_indexes_Gamma.push_back(chunk);
            continue;
        }
        
        fleft_right_g_c_indexes_gamma.push_back(chunk);
        
    }
    
//#ifdef PZDEBUG
//    
//    std::cout << " on Gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_c_indexes_Gamma.size(); k++) {
//        std::cout << " index k : " << k << std::endl;
//        std::cout << " face : " << fleft_right_g_c_indexes_Gamma[k].first << std::endl;
//        std::cout << " volume left : " << fleft_right_g_c_indexes_Gamma[k].second.first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_c_indexes_Gamma[k].second.second <<std::endl;
//    }
//
//    std::cout << " on gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_c_indexes_gamma.size(); k++) {
//        std::cout << " index k : " << k << std::endl;
//        std::cout << " face : " << fleft_right_g_c_indexes_gamma[k].first << std::endl;
//        std::cout << " volume left : " << fleft_right_g_c_indexes_gamma[k].second.first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_c_indexes_gamma[k].second.second <<std::endl;
//    }
//    
//#endif
    
#ifdef PZDEBUG
    if (fleft_right_g_c_indexes_Gamma.size() == 0) {
        DebugStop();
    }
    if (fleft_right_g_c_indexes_gamma.size() == 0) {
        std::cout << "Warning:: No inner interfaces were found" << std::endl;
    }
#endif
    
}

void TRMBuildTransfers::Build_parabolic_hyperbolic_interfaces(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic, bool bc_interfaceQ){
 
    
#ifdef PZDEBUG
    if (!parabolic || !hyperbolic) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = parabolic->Reference();
    parabolic->LoadReferences();
    
    TPZManVector<int64_t,10> qn_dof_indexes;
    int n_shapes;
    
    
    TPZManVector<int64_t> indices;
    std::pair<TPZVec<int64_t> , std::pair< TPZVec<int64_t>, TPZVec<int64_t> > > h_chunk;
    std::pair<int64_t, std::pair< std::pair<int64_t, int64_t> , std::pair<int64_t, int64_t> > >   celint_celh_cp_chunk;
    TPZManVector<int,10> face_sides;
    
    int64_t n_interfaces;
    
    if (bc_interfaceQ) {
        n_interfaces = fleft_right_g_c_indexes_Gamma.size();
        fqn_avg_dof_scatter_Gamma.Resize(n_interfaces);
        fqn_avg_To_hyperbolic_Gamma.Resize(0, 0);
        fcelint_celh_celp_Gamma.resize(n_interfaces);
    }
    else{
        n_interfaces = fleft_right_g_c_indexes_gamma.size();
        fqn_avg_dof_scatter_gamma.Resize(n_interfaces);
        fqn_avg_To_hyperbolic_gamma.Resize(0, 0);
        fcelint_celh_celp_gamma.resize(n_interfaces);
    }
    
    
    // Block size structue Gamma and gamma (inner element interfaces)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(n_interfaces);
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (bc_interfaceQ) {
            h_chunk      = fleft_right_g_c_indexes_Gamma[k_face];
        }
        else{
            h_chunk      = fleft_right_g_c_indexes_gamma[k_face];
        }
        
        
        TPZGeoEl * left_gel    = geometry->Element(h_chunk.second.first[0]);
        TPZGeoEl * right_gel   = geometry->Element(h_chunk.second.second[0]);
        
        TPZCompEl *left_cel = hyperbolic->Element(h_chunk.second.first[1]);
        TPZCompEl *right_cel = hyperbolic->Element(h_chunk.second.second[1]);
        
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        // New method
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        // Identify the father element in mixed mesh that touch or contain the left element in transport mesh
        
        int n_l = fSimulationData->TransporResolution().second;
        
        TPZCompEl * left_mixed_cel = NULL;
        TPZCompEl * right_mixed_cel= NULL;
        
        TPZGeoEl * father_left_gel  = left_gel->Father();
        TPZGeoEl * father_right_gel = right_gel->Father();
        
        int mhm_div_l = fSimulationData->MHMResolution().second.second;
        bool mhm_div_level_check = (left_gel->Level() == mhm_div_l && right_gel->Level() == mhm_div_l);
        
        if ((!father_left_gel && !father_right_gel) || mhm_div_level_check) {
            
#ifdef PZDEBUG
            
            if (fSimulationData->TransporResolution().second !=0) {
                DebugStop();
            }
            
#endif

            left_mixed_cel = left_gel->Reference();
            right_mixed_cel = right_gel->Reference();
        }
        else{
            
        bool hyperbolic_level_check = (left_gel->Level() == mhm_div_l + n_l && right_gel->Level() == mhm_div_l + n_l);

            if (hyperbolic_level_check) {
                for (int l = 0; l < n_l - 1; l++) {
                    father_left_gel     = father_left_gel->Father();
                    father_right_gel    = father_right_gel->Father();
                }
                
                left_mixed_cel = father_left_gel->Reference();
                right_mixed_cel = father_right_gel->Reference();
            }
            else{
                DebugStop();
            }
            
        }

        
#ifdef PZDEBUG
        
        if (bc_interfaceQ) {
            
            if(!left_mixed_cel){
                DebugStop();
            }
        }
        else{
            
            if(!left_mixed_cel || !right_mixed_cel){
                DebugStop();
            }
        }
        
#endif
        
        
        // End
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        
        // Left based element
        this->ComputeFaceIndex(left_mixed_cel->Reference(),face_sides);
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_mixed_cel);
        
        this->ElementDofIndexes(mf_cel, qn_dof_indexes);
        n_shapes = qn_dof_indexes.size();
#ifdef PZDEBUG
        if (qn_dof_indexes.size()==0) {
            DebugStop();
        }
#endif
        
        
        blocks_dimensions[k_face].first = 1;
        blocks_dimensions[k_face].second = n_shapes;
        
        // hyperbolic mesh intercaes and pairs
        celint_celh_cp_chunk.first = h_chunk.first[1];
        celint_celh_cp_chunk.second.first.first     = left_cel->Index();
        celint_celh_cp_chunk.second.first.second    = right_cel->Index();
        
        
        if (bc_interfaceQ) {
            celint_celh_cp_chunk.second.second.first = left_mixed_cel->Index();
            celint_celh_cp_chunk.second.second.second = left_mixed_cel->Index(); // @omar:: left
            fqn_avg_dof_scatter_Gamma[k_face] = qn_dof_indexes;
            fcelint_celh_celp_Gamma[k_face] = celint_celh_cp_chunk;
        }
        else{
            celint_celh_cp_chunk.second.second.first = left_mixed_cel->Index();
            celint_celh_cp_chunk.second.second.second = right_mixed_cel->Index();
            fqn_avg_dof_scatter_gamma[k_face] = qn_dof_indexes;
            fcelint_celh_celp_gamma[k_face] = celint_celh_cp_chunk;
        }
        
    }
    
    // Initialize the matrix
    
    if (bc_interfaceQ) {
        fqn_avg_To_hyperbolic_Gamma.Initialize(blocks_dimensions);
    }
    else{
        fqn_avg_To_hyperbolic_gamma.Initialize(blocks_dimensions);
    }
    
    
    /////////////////////////////////////////////// linear application entries /////////////////////////////////////////
    
    int mesh_index = 0;
    TPZManVector<int64_t,10> dof_indexes;
    
    TPZCompEl * face_cel;
    TPZCompEl * mixed_cel;
    TPZGeoEl * left_gel;
//    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    TPZGeoEl * mixed_gel;
    
    int int_order_interfaces = 1;
    
//    TPZManVector<int64_t> indices;
    std::pair<int64_t, int64_t> duplet;
//    TPZManVector<int,10> face_sides;
    TPZFMatrix<REAL> normals;
//    int64_t face_index;
    
    TPZFNMatrix<100,double> block;
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (bc_interfaceQ) {
            celint_celh_cp_chunk      = fcelint_celh_celp_Gamma[k_face];
        }
        else{
            celint_celh_cp_chunk      = fcelint_celh_celp_gamma[k_face];
        }
        
        face_cel = hyperbolic->Element(celint_celh_cp_chunk.first);
        
#ifdef PZDEBUG
        if (!face_cel) {
            DebugStop();
        }
#endif
        
        face_gel = face_cel->Reference();
      
#ifdef PZDEBUG
        if (!face_gel) {
            DebugStop();
        }
#endif
        
        left_gel =  hyperbolic->Element(celint_celh_cp_chunk.second.first.first)->Reference(); // left
        
#ifdef PZDEBUG
        if (!left_gel) {
            DebugStop();
        }
#endif
        
        if (bc_interfaceQ) {
            mixed_cel = parabolic->Element(celint_celh_cp_chunk.second.second.first); // left
            mixed_gel = mixed_cel->Reference();
        }
        else{
            mixed_cel = parabolic->Element(celint_celh_cp_chunk.second.second.first); // left
            mixed_gel = mixed_cel->Reference();
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        TPZMultiphysicsInterfaceElement * mf_face_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(face_cel);
        
        
        this->ElementDofIndexes(intel_vol, dof_indexes);
        TPZIntPoints *int_points   = face_gel->CreateSideIntegrationRule(face_gel->NSides()-1, int_order_interfaces);
        
        int npoints = int_points->NPoints();
        int nshapes = dof_indexes.size();
        
        // Computing over all integration points of the compuational element mf_cel
        TPZFNMatrix<100,REAL> phi_dot_n(nshapes,1,0.0);
        int face_gel_dim = face_gel->Dimension();
        int intel_vol_dim = intel_vol->Dimension();
        TPZFNMatrix<50,double> block(npoints,nshapes);
        TPZFNMatrix<50,double> block_integral(1,nshapes,0.0);
        
        REAL w;
        TPZManVector<STATE,2> par_duplet(face_gel_dim,0.0);
        TPZManVector<STATE,3> par_mixed_duplet(intel_vol_dim,0.0);
        
        REAL ElementMeasure   = DimensionalMeasure(face_gel);
        
        REAL detjac;
        TPZFMatrix<REAL> jac,axes,jacinv;
        TPZFNMatrix<3,REAL> n;
        TPZVec<int> vectorsides;
        TPZMaterialData data, face_data;
        face_data.fNeedsNormal = true;
        TPZFNMatrix<100,STATE> phi_qs;
        int nphiu,s_i,v_i;
        
        for (int ip = 0; ip < npoints; ip++) {
            
            // Get the vectorial phi
            int_points->Point(ip, par_duplet, w);
            face_gel->Jacobian(par_duplet, jac, axes, detjac, jacinv);
            
            mf_face_cel->ComputeRequiredData(face_data, par_duplet);
            ComputeTransformation(face_gel, left_gel, mixed_gel, par_duplet, par_mixed_duplet);
            
            //            std::cout << "normal = " << face_data.normal <<  std::endl;
            
            intel_vol->InitMaterialData(data);
            intel_vol->ComputeRequiredData(data, par_mixed_duplet);
            
            phi_qs       = data.phi;
            nphiu       = data.fVecShapeIndex.NElements();
            
            for (int iu = 0; iu < nphiu; iu++)
            {
                
                v_i = data.fVecShapeIndex[iu].first;
                s_i = data.fVecShapeIndex[iu].second;
                
                for (int k = 0; k < intel_vol_dim; k++) {
                    phi_dot_n(iu,0) += 1.0 * phi_qs(s_i,0) * data.fNormalVec(k,v_i) * face_data.normal[k];
                }
                
            }
            
            for (int j = 0; j < nshapes; j++) {
                block_integral(0,j) +=  w * detjac * phi_dot_n(j,0)/ElementMeasure;
            }
            
        }
        
//        std::cout << "face index = " << face_gel->Index() <<  std::endl;
//        std::cout << "mat id = " << face_gel->MaterialId() <<  std::endl;
//        std::cout << "k_face = " << k_face <<  std::endl;
//        std::cout << "dof_indexes = " << dof_indexes <<  std::endl;
//        block_integral.Print(std::cout);
        
        if (bc_interfaceQ) {
            fqn_avg_To_hyperbolic_Gamma.SetBlock(k_face, block_integral);
        }
        else{
            fqn_avg_To_hyperbolic_gamma.SetBlock(k_face, block_integral);
        }
        
    }
    
//    fqn_avg_To_hyperbolic_Gamma.Print(" qn to h G = ");
//    fqn_avg_To_hyperbolic_gamma.Print(" qn to h g = ");
    
    return;
}

void TRMBuildTransfers::Build_hyperbolic_parabolic_volumetric(TPZCompMesh * hyperbolic, TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!hyperbolic || !parabolic) {
        DebugStop();
    }
#endif
    
#ifdef PZDEBUG
    if (fparabolic_hyperbolic_cel_pairs.size() == 0) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    int64_t n_el = fparabolic_hyperbolic_cel_pairs.size();
    fsw_avg_dof_scatter.Resize(n_el);
    
    std::pair<int64_t, std::pair< TPZVec<int64_t>, TPZVec<int64_t> > > chunk_intp_indexes;
    
    // Block size structue for Omega
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi_sw_avg(n_el);
    
    int sw_avg_index = 0;
    int sw_avg_points = 0;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fparabolic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.first);
        
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
        
        sw_avg_points = fparabolic_hyperbolic_cel_pairs[iel].second.second.size();
        
        // Getting local integration index
        TPZManVector<int64_t> sw_avg_dof_index(0,0);
        TPZStack<int64_t> sw_avg_dof_indexes;
        
        for (int ison = 0; ison < sw_avg_points ; ison++) {
            
            TPZCompEl * h_cel = hyperbolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            if (sub_gel->Dimension() == dim) {
                this->ElementDofIndexes(mf_h_cel, sw_avg_dof_index, sw_avg_index);
            }
            else{
                DebugStop();
            }
            
            sw_avg_dof_indexes.Push(sw_avg_dof_index[0]);
            
        }

        
        
        fsw_avg_dof_scatter[iel] = sw_avg_dof_indexes;
        
        blocks_dimensions_phi_sw_avg[iel].first = 1;
        blocks_dimensions_phi_sw_avg[iel].second = sw_avg_dof_indexes.size();
        
        
    }
    
    // Initialize the matrix
    fsw_avg_To_parabolic.Initialize(blocks_dimensions_phi_sw_avg);
    
    
    TPZManVector<int64_t> sw_int_point_indexes(0,0);
    std::pair<int64_t, int64_t> sw_avg_block_dim;
    
    // for velocity functions
    TPZMaterialData data;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fparabolic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.first);
        
        
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
        
        int n_sons = fparabolic_hyperbolic_cel_pairs[iel].second.second.size();
        sw_avg_block_dim     = fsw_avg_To_parabolic.GetSizeofBlock(iel);
        
        TPZFMatrix<double> block_phi_sw_avg(sw_avg_block_dim.first,sw_avg_block_dim.second,0.0);
        
        for (int ison = 0; ison < n_sons ; ison++) {
            TPZCompEl * h_cel = hyperbolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            
            REAL volume = this->DimensionalMeasure(gel);
            REAL sub_volume = this->DimensionalMeasure(sub_gel);
            
            block_phi_sw_avg(0,ison) = sub_volume/volume;
            
        }
        
        fsw_avg_To_parabolic.SetBlock(iel, block_phi_sw_avg);
        
    }
    
//    fsw_avg_To_parabolic.Print(" p_avg_to_h ");
    
    return;
    
}

void TRMBuildTransfers::Build_hyperbolic_elliptic_volumetric(TPZCompMesh * hyperbolic, TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!hyperbolic || !elliptic) {
        DebugStop();
    }
#endif
    
#ifdef PZDEBUG
    if (felliptic_hyperbolic_cel_pairs.size() == 0) {
        DebugStop();
    }
#endif
    
    // Loading the links to the geometry (expensive for big geometric meshes)
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = geometry ->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    int64_t n_el = felliptic_hyperbolic_cel_pairs.size();
    fsw_avg_dof_scatter.Resize(n_el);
    
    std::pair<int64_t, std::pair< TPZVec<int64_t>, TPZVec<int64_t> > > chunk_intp_indexes;
    
    // Block size structue for Omega
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions_phi_sw_avg(n_el);
    
    int sw_avg_index = 0;
    int sw_avg_points = 0;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(felliptic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.first);
        
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
        
        sw_avg_points = felliptic_hyperbolic_cel_pairs[iel].second.second.size();
        
        // Getting local integration index
        TPZManVector<int64_t> sw_avg_dof_index(0,0);
        TPZStack<int64_t> sw_avg_dof_indexes;
        
        for (int ison = 0; ison < sw_avg_points ; ison++) {
            
            TPZCompEl * h_cel = hyperbolic->Element(felliptic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            if (sub_gel->Dimension() == dim) {
                this->ElementDofIndexes(mf_h_cel, sw_avg_dof_index, sw_avg_index);
            }
            else{
                DebugStop();
            }
            
            sw_avg_dof_indexes.Push(sw_avg_dof_index[0]);
            
        }
        
        
        
        fsw_avg_dof_scatter[iel] = sw_avg_dof_indexes;
        
        blocks_dimensions_phi_sw_avg[iel].first = 1;
        blocks_dimensions_phi_sw_avg[iel].second = sw_avg_dof_indexes.size();
        
        
    }
    
    // Initialize the matrix
    fsw_avg_To_elliptic.Initialize(blocks_dimensions_phi_sw_avg);
    
    
    TPZManVector<int64_t> sw_int_point_indexes(0,0);
    std::pair<int64_t, int64_t> sw_avg_block_dim;
    
    // for velocity functions
    TPZMaterialData data;
    
    for (int64_t iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(felliptic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(felliptic_hyperbolic_cel_pairs[iel].second.first);
        
        
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
        
        int n_sons = felliptic_hyperbolic_cel_pairs[iel].second.second.size();
        sw_avg_block_dim     = fsw_avg_To_elliptic.GetSizeofBlock(iel);
        
        TPZFMatrix<double> block_phi_sw_avg(sw_avg_block_dim.first,sw_avg_block_dim.second,0.0);
        
        for (int ison = 0; ison < n_sons ; ison++) {
            TPZCompEl * h_cel = hyperbolic->Element(felliptic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            
            REAL volume = this->DimensionalMeasure(gel);
            REAL sub_volume = this->DimensionalMeasure(sub_gel);
            
            block_phi_sw_avg(0,ison) = sub_volume/volume;
            
        }
        
        fsw_avg_To_elliptic.SetBlock(iel, block_phi_sw_avg);
        
    }
    
//    fsw_avg_To_elliptic.Print(" sw_avg_to_e ");
    
    return;
    
}


void TRMBuildTransfers::parabolic_To_hyperbolic_volumetric(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic){
    
#ifdef PZDEBUG
    if (!parabolic || !hyperbolic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_p(fp_avg_To_hyperbolic.Cols(),1,0.0);
    
    int n = fparabolic_hyperbolic_cel_pairs.size();
    int64_t pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fp_avg_dof_scatter[i].size(); iequ++) {
            Scatter_p(pos,0) = parabolic->Solution()(fp_avg_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> p_avg_at_hyperbolic;
    fp_avg_To_hyperbolic.Multiply(Scatter_p,p_avg_at_hyperbolic);
    
    
    TPZGeoMesh * geometry = parabolic->Reference();
    int dim = parabolic->Dimension();
    int n_el = fparabolic_hyperbolic_cel_pairs.size();
    
    int64_t first_point_phi_p_avg = 0;
    
    std::pair<int64_t, int64_t> b_size_phi_p_avg;
    
    b_size_phi_p_avg.first = 0;
    b_size_phi_p_avg.second = 0;
    
    TPZFMatrix<STATE> block_phi_p_avg;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fparabolic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.first);

        
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
        
        first_point_phi_p_avg     += b_size_phi_p_avg.first;
        b_size_phi_p_avg        = fp_avg_To_hyperbolic.GetSizeofBlock(iel);
        
        
        int n_sons = fparabolic_hyperbolic_cel_pairs[iel].second.second.size();
        
        for (int ison = 0; ison < n_sons ; ison++) {
            
            TPZCompEl * h_cel = hyperbolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            
            //  Getting the total integration point of the destination cmesh
            int matd_id = sub_gel->MaterialId();
            TPZMaterial * material = hyperbolic->FindMaterial(matd_id);
            
            if(gel->Dimension() == dim){ // The volumetric ones!
                
                TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
                
                TPZManVector<int64_t, 30> int_point_indexes;
                mf_h_cel->GetMemoryIndices(int_point_indexes);
                int n_points = int_point_indexes.size();
                int64_t ipos;
                
                REAL p_avg;
                for(int64_t ip = 0; ip <  n_points; ip++) {
                    ipos  = int_point_indexes[ip];
                    
                    p_avg       = p_avg_at_hyperbolic(first_point_phi_p_avg + ison,0); // ip -> ison
                    
                    if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                        associated_material->GetMemory()[ipos].Set_p_0(p_avg);
                    }
                    
                    if (fSimulationData->IsCurrentStateQ()) {
                        associated_material->GetMemory()[ipos].Set_p_avg_n(p_avg);
                    }
                    else{
                        associated_material->GetMemory()[ipos].Set_p_avg(p_avg);
                    }
                    
                }
                
                
            }
            else{
                DebugStop(); // Just volumetric coupling
            }
        }
        
    }
    
    
}

void TRMBuildTransfers::parabolic_To_hyperbolic_interfaces(TPZCompMesh * parabolic, TPZCompMesh * hyperbolic, bool bc_interfaceQ){
    
#ifdef PZDEBUG
    if (!parabolic || !hyperbolic) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    TPZManVector<int64_t,30>  point_index_trans;
    TPZManVector<int64_t,30>  point_index_l;
    TPZManVector<int64_t,30>  point_index_r;
    
    
    TPZGeoMesh * geometry = hyperbolic->Reference();
    int64_t n_interfaces;
//    int dim = geometry->Dimension();
    
    if (bc_interfaceQ) {
        n_interfaces = fleft_right_g_c_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_c_indexes_gamma.size();
    }
    
    if (bc_interfaceQ) {
        
        // Step one
        TPZFMatrix<STATE> ScatterFluxes(fqn_avg_To_hyperbolic_Gamma.Cols(),1,0.0);
        int64_t pos = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {
            for(int iflux = 0; iflux < fqn_avg_dof_scatter_Gamma[iface].size(); iflux++) {
                ScatterFluxes(pos,0) = parabolic->Solution()(fqn_avg_dof_scatter_Gamma[iface][iflux],0);
                pos++;
            }
        }
        
        // Step two
        TPZFMatrix<STATE> qn_at_intpoints;
        fqn_avg_To_hyperbolic_Gamma.Multiply(ScatterFluxes,qn_at_intpoints);
    
        
        // Step three
        // Trasnfering integrated normal fluxes values
        TPZVec<int64_t> p_point_indices;
        TPZVec<int64_t> left_p_point_indices;
        
        for (int iface = 0; iface < n_interfaces; iface++) {
            
            int material_id = geometry->Element(fleft_right_g_c_indexes_Gamma[iface].first[0])->MaterialId();
            TPZMaterial * material = hyperbolic->FindMaterial(material_id);
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  * material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> *>(material);
            
#ifdef PZDEBUG
            if (!material_bc_mem) {
                DebugStop();
            }
            
#endif
            TPZCompEl * face_cel = hyperbolic->Element(fcelint_celh_celp_Gamma[iface].first);
            TPZCompEl * left_cel = hyperbolic->Element(fcelint_celh_celp_Gamma[iface].second.first.first); // left cel
            
            
#ifdef PZDEBUG
            if (!face_cel || !left_cel) {
                DebugStop();
            }
            
#endif
            
            TPZMultiphysicsInterfaceElement * mf_face_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(face_cel);
            TPZMultiphysicsElement * mf_left_cel = dynamic_cast<TPZMultiphysicsElement * >(left_cel);
            
#ifdef PZDEBUG
            if (!mf_face_cel || !mf_left_cel) {
                DebugStop();
            }
            
#endif
            
            mf_face_cel->GetMemoryIndices(p_point_indices);
            mf_left_cel->GetMemoryIndices(left_p_point_indices);
            
            int n_points = p_point_indices.size();
            
#ifdef PZDEBUG
            if (n_points != 1 || left_p_point_indices.size() == 0) {
                DebugStop();
            }
            
#endif
            
            int left_material_id = geometry->Element(fleft_right_g_c_indexes_Gamma[iface].second.first[0])->MaterialId();
            TPZMaterial * left_material = hyperbolic->FindMaterial(left_material_id);
            TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin>  * left_material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(left_material);
            
#ifdef PZDEBUG
            if (!left_material_mem) {
                DebugStop();
            }
            
#endif
            int left_pos = left_p_point_indices[0];
            REAL & p_avg_n = left_material_mem->GetMemory()[left_pos].p_avg_n();
            REAL & p_avg = left_material_mem->GetMemory()[left_pos].p_avg();
            REAL s_a_n = left_material_mem->GetMemory()[left_pos].sa_n();
            REAL s_a = left_material_mem->GetMemory()[left_pos].sa();
            
            int pos = p_point_indices[0];
            material_bc_mem->GetMemory()[pos].Set_un(qn_at_intpoints(iface,0));
            material_bc_mem->GetMemory()[pos].Set_p_avg_n_l(p_avg_n);
            material_bc_mem->GetMemory()[pos].Set_sa_n_l(s_a_n);
            material_bc_mem->GetMemory()[pos].Set_p_avg_l(p_avg);
            material_bc_mem->GetMemory()[pos].Set_sa_l(s_a);
            
        }
        
        
    }
    else{
        
        // Step one
        TPZFMatrix<STATE> ScatterFluxes(fqn_avg_To_hyperbolic_gamma.Cols(),1,0.0);
        int64_t pos = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {
            for(int iflux = 0; iflux < fqn_avg_dof_scatter_gamma[iface].size(); iflux++) {
                ScatterFluxes(pos,0) = parabolic->Solution()(fqn_avg_dof_scatter_gamma[iface][iflux],0);
                pos++;
            }
        }
        
        // Step two
        TPZFMatrix<STATE> qn_at_intpoints;
        fqn_avg_To_hyperbolic_gamma.Multiply(ScatterFluxes,qn_at_intpoints);
        
        
        // Step three
        // Trasnfering integrated normal fluxes values
        TPZVec<int64_t> p_point_indices;
        TPZVec<int64_t> left_p_point_indices,right_p_point_indices;
        
        for (int iface = 0; iface < n_interfaces; iface++) {
            
            int material_id = geometry->Element(fleft_right_g_c_indexes_gamma[iface].first[0])->MaterialId();
            TPZMaterial * material = hyperbolic->FindMaterial(material_id);
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin>  * material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin> *>(material);
            
#ifdef PZDEBUG
            if (!material_mem) {
                DebugStop();
            }
            
#endif
            TPZCompEl * face_cel = hyperbolic->Element(fcelint_celh_celp_gamma[iface].first);
            TPZCompEl * left_cel = hyperbolic->Element(fcelint_celh_celp_gamma[iface].second.first.first); // left cel
            TPZCompEl * right_cel = hyperbolic->Element(fcelint_celh_celp_gamma[iface].second.first.second); // right cel
            
            
#ifdef PZDEBUG
            if (!face_cel || !left_cel || !right_cel) {
                DebugStop();
            }
            
#endif
            
            TPZMultiphysicsInterfaceElement * mf_face_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(face_cel);
            TPZMultiphysicsElement * mf_left_cel = dynamic_cast<TPZMultiphysicsElement * >(left_cel);
            TPZMultiphysicsElement * mf_right_cel = dynamic_cast<TPZMultiphysicsElement * >(right_cel);
            
#ifdef PZDEBUG
            if (!mf_face_cel || !mf_left_cel || !mf_right_cel) {
                DebugStop();
            }
            
#endif
            
            mf_face_cel->GetMemoryIndices(p_point_indices);
            mf_left_cel->GetMemoryIndices(left_p_point_indices);
            mf_right_cel->GetMemoryIndices(right_p_point_indices);
            
            int n_points = p_point_indices.size();
            
#ifdef PZDEBUG
            if (n_points != 1 || left_p_point_indices.size() == 0 || right_p_point_indices.size() == 0) {
                DebugStop();
            }
            
#endif
            
            int left_material_id = geometry->Element(fleft_right_g_c_indexes_gamma[iface].second.first[0])->MaterialId();
            TPZMaterial * left_material = hyperbolic->FindMaterial(left_material_id);
            TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin>  * left_material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(left_material);
            
            int right_material_id = geometry->Element(fleft_right_g_c_indexes_gamma[iface].second.second[0])->MaterialId();
            TPZMaterial * right_material = hyperbolic->FindMaterial(right_material_id);
            TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin>  * right_material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(right_material);
            
#ifdef PZDEBUG
            if (!left_material_mem || !right_material_mem) {
                DebugStop();
            }
            
#endif
            int left_pos = left_p_point_indices[0];
            int right_pos = right_p_point_indices[0];
            
            REAL & p_avg_n_l = left_material_mem->GetMemory()[left_pos].p_avg_n();
            REAL & p_avg_n_r = right_material_mem->GetMemory()[right_pos].p_avg_n();
            REAL & p_avg_l = left_material_mem->GetMemory()[left_pos].p_avg();
            REAL & p_avg_r = right_material_mem->GetMemory()[right_pos].p_avg();
            
            REAL s_a_n_l = left_material_mem->GetMemory()[left_pos].sa_n();
            REAL s_a_n_r = right_material_mem->GetMemory()[right_pos].sa_n();
            REAL s_a_l = left_material_mem->GetMemory()[left_pos].sa();
            REAL s_a_r = right_material_mem->GetMemory()[right_pos].sa();
            
            TPZFMatrix<REAL> & k_l =  left_material_mem->GetMemory()[left_pos].K_0();
            TPZFMatrix<REAL> & k_r =  right_material_mem->GetMemory()[right_pos].K_0();
            TPZFMatrix<REAL>  k_avg(3,3);
            REAL epsilon = 1.0e-25;
            
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    k_avg(i,j) = 2.0*(k_l(i,j)*k_r(i,j))/(k_l(i,j) + k_r(i,j) + epsilon);
                }
            }
            
            int pos = p_point_indices[0];
            material_mem->GetMemory()[pos].Set_un(qn_at_intpoints(iface,0));
            material_mem->GetMemory()[pos].Set_p_avg_n_l(p_avg_n_l);
            material_mem->GetMemory()[pos].Set_p_avg_n_r(p_avg_n_r);
            material_mem->GetMemory()[pos].Set_sa_n_l(s_a_n_l);
            material_mem->GetMemory()[pos].Set_sa_n_r(s_a_n_r);
            
            material_mem->GetMemory()[pos].Set_p_avg_l(p_avg_l);
            material_mem->GetMemory()[pos].Set_p_avg_r(p_avg_r);
            material_mem->GetMemory()[pos].Set_sa_l(s_a_l);
            material_mem->GetMemory()[pos].Set_sa_r(s_a_r);
            
            material_mem->GetMemory()[pos].Set_K_0(k_avg);
        }
        
    }
    
    
    return;

}

void TRMBuildTransfers::hyperbolic_To_parabolic_volumetric(TPZCompMesh * hyperbolic, TPZCompMesh * parabolic){
    
#ifdef PZDEBUG
    if (!hyperbolic || !parabolic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_sw(fsw_avg_To_parabolic.Cols(),1,0.0);
    
    int n = fparabolic_hyperbolic_cel_pairs.size();
    int64_t pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fsw_avg_dof_scatter[i].size(); iequ++) {
            Scatter_sw(pos,0) = hyperbolic->Solution()(fsw_avg_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> sw_avg_at_parabolic;
    fsw_avg_To_parabolic.Multiply(Scatter_sw,sw_avg_at_parabolic);
    
    
    TPZGeoMesh * geometry = hyperbolic->Reference();
    int dim = hyperbolic->Dimension();
    int n_el = fparabolic_hyperbolic_cel_pairs.size();
    
    int64_t first_point_phi_sw_avg = 0;
    
    std::pair<int64_t, int64_t> b_size_phi_sw_avg;
    
    b_size_phi_sw_avg.first = 0;
    b_size_phi_sw_avg.second = 0;
    
    TPZFMatrix<STATE> block_phi_sw_avg;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(fparabolic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * p_cel = parabolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.first);
        
        
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
        
        first_point_phi_sw_avg     += b_size_phi_sw_avg.first;
        b_size_phi_sw_avg        = fsw_avg_To_parabolic.GetSizeofBlock(iel);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = parabolic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_p_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            REAL sw_avg;
            sw_avg       = sw_avg_at_parabolic(first_point_phi_sw_avg,0);
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_sa_n(sw_avg);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_sa(sw_avg);
                }
                
            }
            
            
        }
        else{
            DebugStop(); // Just volumetric coupling
        }
    }
    
    
}
void TRMBuildTransfers::elliptic_To_hyperbolic(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic){
    this->elliptic_To_hyperbolic_volumetricII(elliptic,hyperbolic);
    return;
}

void TRMBuildTransfers::hyperbolic_To_elliptic(TPZCompMesh * hyperbolic,TPZCompMesh * elliptic){
    this->hyperbolic_To_elliptic_volumetric(hyperbolic,elliptic);
    return;
}

void TRMBuildTransfers::elliptic_To_hyperbolic_volumetric(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic){
    
#ifdef PZDEBUG
    if (!elliptic || !hyperbolic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_u(fgrad_u_avg_To_hyperbolic.Cols(),1,0.0);
    
    int n = felliptic_hyperbolic_cel_pairs.size();
    int64_t pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fu_avg_dof_scatter[i].size(); iequ++) {
            Scatter_u(pos,0) = elliptic->Solution()(fu_avg_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> grad_u_avg_at_hyperbolic;
    fgrad_u_avg_To_hyperbolic.Multiply(Scatter_u,grad_u_avg_at_hyperbolic);
    
    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = elliptic->Dimension();
    int n_el = felliptic_hyperbolic_cel_pairs.size();
    
    int64_t first_point_grad_phi_u_avg = 0;
    
    std::pair<int64_t, int64_t> b_size_grad_phi_u_avg;
    
    b_size_grad_phi_u_avg.first = 0;
    b_size_grad_phi_u_avg.second = 0;
    
    TPZFMatrix<STATE> block_grad_phi_u_avg;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(felliptic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(felliptic_hyperbolic_cel_pairs[iel].second.first);
        
        
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
        
        first_point_grad_phi_u_avg     += b_size_grad_phi_u_avg.first;
        b_size_grad_phi_u_avg        = fgrad_u_avg_To_hyperbolic.GetSizeofBlock(iel);
        
        
        int n_sons = felliptic_hyperbolic_cel_pairs[iel].second.second.size();
        
        for (int ison = 0; ison < n_sons ; ison++) {
            
            TPZCompEl * h_cel = hyperbolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            
            //  Getting the total integration point of the destination cmesh
            int matd_id = sub_gel->MaterialId();
            TPZMaterial * material = hyperbolic->FindMaterial(matd_id);
            
            if(gel->Dimension() == dim){ // The volumetric ones!
                
                TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
                
                TPZManVector<int64_t, 30> int_point_indexes;
                mf_h_cel->GetMemoryIndices(int_point_indexes);
                int n_points = int_point_indexes.size();
                int64_t ipos;
                
                TPZFMatrix<REAL> grad_u(dim,dim,0.0);
                for(int64_t ip = 0; ip <  n_points; ip++) {
                    ipos  = int_point_indexes[ip];
                    
                    
                    TPZFNMatrix<3,STATE> u(1,3,0.0);
                    TPZFNMatrix<9,STATE> grad_u(3,3,0.0);
                    for(int64_t ip = 0; ip <  n_points; ip++) {
                        ipos  = int_point_indexes[ip];
                    
                        for (int id = 0; id < dim ; id++) {
                            for (int jd = 0; jd < dim ; jd++) {
                                if(fSimulationData->ReducedBasisResolution().first){
                                    grad_u(id,jd) = grad_u_avg_at_hyperbolic(first_point_grad_phi_u_avg + ison*dim*dim + id*dim + jd,0);
                                }
                                else{
                                    grad_u(jd,id) = grad_u_avg_at_hyperbolic(first_point_grad_phi_u_avg + ison*dim*dim + id*dim + jd,0);
                                }
                            }
                        }
                            
                        if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                            associated_material->GetMemory()[ipos].Set_grad_u_0(grad_u);
                        }
                        
                        if (fSimulationData->IsCurrentStateQ()) {
                            associated_material->GetMemory()[ipos].Set_grad_u_n(grad_u);
                        }
                        else{
                            associated_material->GetMemory()[ipos].Set_grad_u(grad_u);
                        }
                    }
                }
                
                
            }
            else{
                DebugStop(); // Just volumetric coupling
            }
        }
        
    }
    
    
}

void TRMBuildTransfers::elliptic_To_hyperbolic_volumetricII(TPZCompMesh * elliptic, TPZCompMesh * hyperbolic){

#ifdef PZDEBUG
    if (!elliptic || !hyperbolic) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = elliptic->Reference();
    int dim = elliptic->Dimension();
    int n_el = felliptic_hyperbolic_cel_pairs.size();
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(felliptic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(felliptic_hyperbolic_cel_pairs[iel].second.first);
        
        
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
        
        REAL vol_e = DimensionalMeasure(gel);
        int matd_id_e = gel->MaterialId();
        TPZMaterial * material_e = elliptic->FindMaterial(matd_id_e);
        
        TPZManVector<int64_t, 30> int_point_indexes_e;
        mf_e_cel->GetMemoryIndices(int_point_indexes_e);
        int n_points_e = int_point_indexes_e.size();
        int64_t ipos_e;
        
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material_e = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material_e);
        
        TPZFMatrix<REAL> grad_u_avg(3,3,0.0);
        for(int64_t ip = 0; ip <  n_points_e; ip++) {
            ipos_e = int_point_indexes_e[ip];
            if(fSimulationData->IsCurrentStateQ()){
                TPZFMatrix<REAL> & grad_u = associated_material_e->GetMemory()[ipos_e].grad_u_n();
                grad_u *= 1.0/vol_e;
                grad_u_avg += grad_u;
            }
            else{
                TPZFMatrix<REAL> & grad_u = associated_material_e->GetMemory()[ipos_e].grad_u();
                grad_u *= 1.0/vol_e;
                grad_u_avg += grad_u;
            }
            
        }
        
        
        int n_sons = felliptic_hyperbolic_cel_pairs[iel].second.second.size();
        
        for (int ison = 0; ison < n_sons ; ison++) {
            
            TPZCompEl * h_cel = hyperbolic->Element(fparabolic_hyperbolic_cel_pairs[iel].second.second[ison]);
            
#ifdef PZDEBUG
            if (!h_cel) {
                DebugStop();
            }
#endif
            
            TPZMultiphysicsElement * mf_h_cel = dynamic_cast<TPZMultiphysicsElement * >(h_cel);
            
#ifdef PZDEBUG
            if(!mf_h_cel)
            {
                DebugStop();
            }
#endif
            
            TPZGeoEl * sub_gel = h_cel->Reference();
            
#ifdef PZDEBUG
            if(!sub_gel)
            {
                DebugStop();
            }
#endif
            
            //  Getting the total integration point of the destination cmesh
            int matd_id = sub_gel->MaterialId();
            TPZMaterial * material = hyperbolic->FindMaterial(matd_id);
            
            if(gel->Dimension() == dim){ // The volumetric ones!
                
                TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
                
                TPZManVector<int64_t, 30> int_point_indexes;
                mf_h_cel->GetMemoryIndices(int_point_indexes);
                int n_points = int_point_indexes.size();
                int64_t ipos;
                
                TPZFMatrix<REAL> grad_u(3,3,0.0);
                for(int64_t ip = 0; ip <  n_points; ip++) {
                    ipos  = int_point_indexes[ip];
                    
                    TPZFNMatrix<9,STATE> grad_u(3,3,0.0);
                    for(int64_t ip = 0; ip <  n_points; ip++) {
                        ipos  = int_point_indexes[ip];
                        
                        REAL vol_h = DimensionalMeasure(sub_gel);
                        grad_u = (vol_h)*grad_u_avg;
                        
                        if(fSimulationData->IsInitialStateQ() && fSimulationData->IsCurrentStateQ()){
                            associated_material->GetMemory()[ipos].Set_grad_u_0(grad_u);
                        }
                        
                        if (fSimulationData->IsCurrentStateQ()) {
                            associated_material->GetMemory()[ipos].Set_grad_u_n(grad_u);
                        }
                        else{
                            associated_material->GetMemory()[ipos].Set_grad_u(grad_u);
                        }
                    }
                }
                
                
            }
            else{
                DebugStop(); // Just volumetric coupling
            }
        }
        
    }
    
}



void TRMBuildTransfers::hyperbolic_To_elliptic_volumetric(TPZCompMesh * hyperbolic, TPZCompMesh * elliptic){
    
#ifdef PZDEBUG
    if (!hyperbolic || !elliptic) {
        DebugStop();
    }
#endif
    
    
    // Step zero scatter
    TPZFMatrix<STATE> Scatter_sw(fsw_avg_To_elliptic.Cols(),1,0.0);
    
    int n = felliptic_hyperbolic_cel_pairs.size();
    int64_t pos = 0;
    for (int i = 0; i < n; i++) {
        for(int iequ = 0; iequ < fsw_avg_dof_scatter[i].size(); iequ++) {
            Scatter_sw(pos,0) = hyperbolic->Solution()(fsw_avg_dof_scatter[i][iequ],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> sw_avg_at_elliptic;
    fsw_avg_To_elliptic.Multiply(Scatter_sw,sw_avg_at_elliptic);
    
    
    TPZGeoMesh * geometry = hyperbolic->Reference();
    int dim = hyperbolic->Dimension();
    int n_el = felliptic_hyperbolic_cel_pairs.size();
    
    int64_t first_point_phi_sw_avg = 0;
    
    std::pair<int64_t, int64_t> b_size_phi_sw_avg;
    
    b_size_phi_sw_avg.first = 0;
    b_size_phi_sw_avg.second = 0;
    
    TPZFMatrix<STATE> block_phi_sw_avg;
    
    for (int iel = 0; iel < n_el; iel++) {
        
        TPZGeoEl * gel = geometry->Element(felliptic_hyperbolic_cel_pairs[iel].first);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZCompEl * e_cel = elliptic->Element(felliptic_hyperbolic_cel_pairs[iel].second.first);
        
        
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
        
        first_point_phi_sw_avg     += b_size_phi_sw_avg.first;
        b_size_phi_sw_avg        = fsw_avg_To_elliptic.GetSizeofBlock(iel);
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        TPZMaterial * material = elliptic->FindMaterial(matd_id);
        
        if(gel->Dimension() == dim){ // The volumetric ones!
            
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
            
            TPZManVector<int64_t, 30> int_point_indexes;
            mf_e_cel->GetMemoryIndices(int_point_indexes);
            int n_points = int_point_indexes.size();
            int64_t ipos;
            
            REAL sw_avg;
            sw_avg       = sw_avg_at_elliptic(first_point_phi_sw_avg,0);
            for(int64_t ip = 0; ip <  n_points; ip++) {
                ipos  = int_point_indexes[ip];
                
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[ipos].Set_sa_n(sw_avg);
                }
                else{
                    associated_material->GetMemory()[ipos].Set_sa(sw_avg);
                }
                
            }
            
            
        }
        else{
            DebugStop(); // Just volumetric coupling
        }
    }
    
    
}

////////////////////////// Transfers:: Iterative Coupling by Operator Splitting //////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////









































/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Matrices Initialization Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TRMBuildTransfers::Initialize_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = cmesh_multiphysics->Reference()->Dimension(); // vectorial
    int64_t element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fu_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(nel);
    
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(intel->Dimension() < n_var_dim){
            // there is boundary elements for normal flux where it is a scalar variable
//            mf_cel->GetMemoryIndices(int_point_indexes);
//            this->ElementDofIndexes(intel, dof_indexes);
//            fu_dof_scatter[element_index] = dof_indexes;
            blocks_dimensions[element_index].first = 0;
            blocks_dimensions[element_index].second = 0;
            fu_dof_scatter[element_index] = dof_indexes;
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);
        fu_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
        fu_dof_scatter[element_index] = dof_indexes;
    }
    
    // Initialize the matrix
    fu_To_Mixed.Initialize(blocks_dimensions);
    
}


void TRMBuildTransfers::Initialize_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    int64_t element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fp_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(nel);
    
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(!intel){
            // there is no boundary elements for pressure
            blocks_dimensions[element_index].first = 0*n_var_dim;
            blocks_dimensions[element_index].second = 0;
            fp_dof_scatter[element_index] = dof_indexes;
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);
        fp_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
        fp_dof_scatter[element_index] = dof_indexes;
    }
    
    // Initialize the matrix
    fp_To_Mixed.Initialize(blocks_dimensions);
    
    
}


void TRMBuildTransfers::Initialize_s_To_Transport(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    cmesh_multiphysics->LoadReferences();
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    int64_t element_index = 0;
    int dimension = cmesh_multiphysics->Reference()->Dimension();
    // Compute destination index scatter by element (Omega and Gamma)
    if (!mesh_index) {
        fsa_dof_scatter.Resize(nel);
    }
    else{
        fsb_dof_scatter.Resize(nel);
    }

    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(nel);
    
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
        if(gel->HasSubElement()){
            continue;
        }
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
       
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);        
        
#ifdef PZDEBUG
        if(!mf_cel)
        {
            if(!mf_int_cel){
                DebugStop();
            }
        }
#endif
        if (mf_int_cel) {
            continue;
        }
        
        element_index = mf_cel->Index();
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(intel->Dimension() < dimension){
            // there is no boundary elements for saturation
            blocks_dimensions[element_index].first = 0*n_var_dim;
            blocks_dimensions[element_index].second = 0;
            if (!mesh_index) {
                fsa_dof_scatter[element_index] = dof_indexes;
            }
            else{
                fsb_dof_scatter[element_index] = dof_indexes;
            }
            
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel, dof_indexes);

        if (!mesh_index) {
            fsa_dof_scatter[element_index] = dof_indexes;
        }
        else{
            fsb_dof_scatter[element_index] = dof_indexes;
        }
        
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
    }
    
    // Initialize the matrix
    fs_To_Transport.Initialize(blocks_dimensions);
    
    
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Matrices Filling Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TRMBuildTransfers::Fill_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_u_To_Mixed(cmesh_multiphysics, mesh_index);
    
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = cmesh_multiphysics->Reference()->Dimension();; // vector
    int64_t element_index = 0;
    
    TPZMaterialData data;
    
    std::pair<int64_t, int64_t> block_dim;
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        dof_indexes = fu_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block;
        
        if(intel->Dimension() < n_var_dim){ // lower dimensional elements
            block.Resize(block_dim.first,block_dim.second);        }
        else{
            block.Resize(block_dim.first*n_var_dim,block_dim.second);
        }
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            // Get the vectorial phi
            intel->Shape(qsi, phi, dphidxi);
            intel->InitMaterialData(data);
            intel->ComputeRequiredData(data,qsi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    int vector_index = data.fVecShapeIndex[jp].first;
                    int shape_index = data.fVecShapeIndex[jp].second;
                    block(ip*n_var_dim+id,jp) = phi(shape_index,0)*data.fNormalVec(id,vector_index);
                }
            }
            
        }
        
        fu_To_Mixed.SetBlock(element_index, block);
        
    }
    
    return;
}


void TRMBuildTransfers::Fill_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_p_To_Mixed(cmesh_multiphysics, mesh_index);
    
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    int64_t element_index = 0;
    
    std::pair<int64_t, int64_t> block_dim;
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(!intel){
            continue;
        }
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        dof_indexes = fp_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block(block_dim.first*n_var_dim,block_dim.second);
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            intel->Shape(qsi, phi, dphidxi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    block(ip+id*n_var_dim,jp) = phi(jp,0);
                }
            }
            
        }
        
        fp_To_Mixed.SetBlock(element_index, block);
        
    }
    
    return;
}



void TRMBuildTransfers::Fill_s_To_Transport(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_s_To_Transport(cmesh_multiphysics, mesh_index);
    
    int64_t nel = cmesh_multiphysics->NElements();
    int n_var_dim = 1; // scalar
    int64_t element_index = 0;
    int dimension = cmesh_multiphysics->Reference()->Dimension();
    std::pair<int64_t, int64_t> block_dim;
    
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        
        TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (mf_int_cel) {
            continue;
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        element_index = mf_cel->Index();
        
        // Getting local integration index
        TPZManVector<int64_t> int_point_indexes(0,0);
        TPZManVector<int64_t> dof_indexes(0,0);
        
        if(intel->Dimension() < dimension){
            continue;
        }
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        if (!mesh_index) {
            dof_indexes = fsa_dof_scatter[element_index];
        }
        else{
            dof_indexes = fsb_dof_scatter[element_index];
        }
        
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_cel->GetIntegrationRule();
        int np_cel = int_points_mixed.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(intel->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel->NShapeF(),0.0);
        TPZFMatrix<double> block(block_dim.first*n_var_dim,block_dim.second);
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_mixed.Point(ip, qsi, w);
            intel->Shape(qsi, phi, dphidxi);
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < block_dim.second; jp++) {
                    block(ip+id*n_var_dim,jp) = phi(jp,0);
                }
            }
            
        }
        
        fs_To_Transport.SetBlock(element_index, block);
        
    }
    return;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Transfer Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TRMBuildTransfers::kappa_phi_To_Mixed_Memory(TPZCompMesh * cmesh_multiphysics){
    
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int dim = cmesh_multiphysics->Dimension();
    
    
    TPZFMatrix<STATE> kappa, kappa_inv;
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;

    // Step one
    int n_elements = cmesh_multiphysics->NElements();
    TPZManVector<int64_t, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);
//            fSimulationData->Map()->ComputePropertieSPE10Map(index, x, kappa, kappa_inv, porosity);
            fSimulationData->Map()->Kappa(x, kappa, kappa_inv, vars);
            fSimulationData->Map()->phi(x, porosity, vars);
            associated_material->GetMemory()[indexes[ip]].Set_K_0(kappa);
            associated_material->GetMemory()[indexes[ip]].Set_Kinv_0(kappa_inv);
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
        }
    }
    
}


void TRMBuildTransfers::u_To_Mixed_Memory(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_multiphysics){
    

#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    int dim = cmesh_flux->Dimension();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterFlux(fu_To_Mixed.Cols(),1,0.0);
    int64_t pos = 0;
    for (int el = 0; el < nel; el++) {
        for(int ip = 0; ip < fu_dof_scatter[el].size(); ip++) {
            ScatterFlux(pos,0) = cmesh_flux->Solution()(fu_dof_scatter[el][ip],0);
            pos++;
        }
    }
    
    // Step two
    TPZFMatrix<STATE> Flux_at_intpoints;
    fu_To_Mixed.Multiply(ScatterFlux,Flux_at_intpoints);
    // Trasnfering values
    TPZManVector<STATE,3> u(dim,0.0);
    for(int64_t i = 0; i <  np_cmesh; i++){
        for (int id = 0; id < dim ; id++) {
            u[id]= Flux_at_intpoints(i*dim+id,0);
        }

        if(fSimulationData->IsCurrentStateQ()){
            associated_material->GetMemory()[i].Set_u(u);
        }
        else{
            associated_material->GetMemory()[i].Set_u_n(u);
        }
        
    }
    
}

void TRMBuildTransfers::p_To_Mixed_Memory(TPZCompMesh * cmesh_pressure, TPZCompMesh * cmesh_multiphysics){

    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];

    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterPressure(fp_To_Mixed.Cols(),1,0.0);
    int64_t pos = 0;
    for (int el = 0; el < nel; el++) {
        for(int ip = 0; ip < fp_dof_scatter[el].size(); ip++) {
            ScatterPressure(pos,0) = cmesh_pressure->Solution()(fp_dof_scatter[el][ip],0);
            pos++;
        }
    }
    
    // Step two
    TPZFNMatrix<30,STATE> Pressure_at_intpoints;
    fp_To_Mixed.Multiply(ScatterPressure,Pressure_at_intpoints);
    // Trasnfering values
    for(int64_t i = 0; i <  np_cmesh; i++){
        if(fSimulationData->IsCurrentStateQ()){
            associated_material->GetMemory()[i].Set_p_n(Pressure_at_intpoints(i,0));
        }
        else{
            associated_material->GetMemory()[i].Set_p(Pressure_at_intpoints(i,0));
        }
    }
    
}

void TRMBuildTransfers::kappa_phi_To_Transport_Memory(TPZCompMesh * cmesh_multiphysics){
    
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int dim = cmesh_multiphysics->Dimension();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
    

    TPZFMatrix<STATE> kappa, kappa_inv;
    TPZManVector<STATE, 10> vars;
    TPZManVector<STATE, 10> porosity;
//    REAL porosity;
//    int64_t index = 0;
    // Step one
    int n_elements = cmesh_multiphysics->NElements();
    TPZManVector<int64_t, 30> indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension()!= dim) {
            continue;
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
        
#ifdef PZDEBUG
        if (!mf_cel) {
            DebugStop();
        }
#endif
        
        
        const TPZIntPoints & int_points = mf_cel->GetIntegrationRule();
        int np = int_points.NPoints();
        GlobalPointIndexes(cel, indexes);
        
#ifdef PZDEBUG
        if (indexes.size() != np) {
            DebugStop();
        }
#endif
        
        TPZManVector<REAL,3> par_triplet(3,0.0);
        TPZManVector<REAL,3> x(3,0.0);
        REAL w;
        for (int ip = 0; ip<np; ip++) {
            int_points.Point(ip, par_triplet, w);
            gel->X(par_triplet, x);
            
            associated_material->GetMemory()[indexes[ip]].Set_x(x);
            //            fSimulationData->Map()->ComputePropertieSPE10Map(index, x, kappa, kappa_inv, porosity);
            fSimulationData->Map()->Kappa(x, kappa, kappa_inv, vars);
            fSimulationData->Map()->phi(x, porosity, vars);
            associated_material->GetMemory()[indexes[ip]].Set_K_0(kappa);
            associated_material->GetMemory()[indexes[ip]].Set_Kinv_0(kappa_inv);
            associated_material->GetMemory()[indexes[ip]].Set_phi_0(porosity[0]);
        }
    }
    
}


void TRMBuildTransfers::s_To_Transport_Memory(TPZCompMesh * cmesh_saturation, TPZCompMesh * cmesh_multiphysics, int mesh_index){

    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = cmesh_multiphysics->NElements();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step one
    TPZFMatrix<STATE> ScatterSaturation(fs_To_Transport.Cols(),1,0.0);
    long pos = 0;
    for (int el = 0; el < nel; el++) {
        
        if (!mesh_index) {
            for(int ip = 0; ip < fsa_dof_scatter[el].size(); ip++) {
                ScatterSaturation(pos,0) = cmesh_saturation->Solution()(fsa_dof_scatter[el][ip],0);
                pos++;
            }
        }
        else{
            for(int ip = 0; ip < fsb_dof_scatter[el].size(); ip++) {
                ScatterSaturation(pos,0) = cmesh_saturation->Solution()(fsb_dof_scatter[el][ip],0);
                pos++;
            }
        }
        
    }
    
    // Step two
    TPZFMatrix<STATE> Saturation_at_intpoints;
    fs_To_Transport.Multiply(ScatterSaturation,Saturation_at_intpoints);
    // Trasnfering values
    for(long i = 0; i <  np_cmesh; i++){
        
        if(!mesh_index){
            if(fSimulationData->IsCurrentStateQ()){
                associated_material->GetMemory()[i].Set_sa_n(Saturation_at_intpoints(i,0));
            }
            else{
                associated_material->GetMemory()[i].Set_sa(Saturation_at_intpoints(i,0));
            }
        }
        else{
            if(fSimulationData->IsCurrentStateQ()){
                associated_material->GetMemory()[i].Set_sb_n(Saturation_at_intpoints(i,0));
            }
            else{
                associated_material->GetMemory()[i].Set_sb(Saturation_at_intpoints(i,0));
            }
        }
        
    }
    
}

/** @brief Reciprocal (mixed <-> transpor) transfer average quantities to integration points of multiphysics meshes over volumetric elements */
void TRMBuildTransfers::Reciprocal_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_trans){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_cindexes.size() == 0 ) {
        DebugStop();
    }

    if (!cmesh_mf_mixed || !cmesh_mf_trans) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    TPZManVector<int64_t,30> p_point_indexes;
    TPZManVector<int64_t,30> s_point_indexes;
    int64_t nvolumes = fmixed_transport_cindexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
    TPZGeoEl  * gel = geometry->Element(fmixed_transport_cindexes[ivol].first);
    TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_cindexes[ivol].second.first);
    TPZCompEl * trans_cel = cmesh_mf_trans->Element(fmixed_transport_cindexes[ivol].second.second);
        
#ifdef PZDEBUG
        if (!mixed_cel || !trans_cel || !gel) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
        
        TPZMaterial * trans_material = cmesh_mf_trans->FindMaterial(rockid);
        TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * trans_memory = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(trans_material);

        REAL element_measure = DimensionalMeasure(mixed_cel->Reference());
    
        GlobalPointIndexes(mixed_cel, p_point_indexes);
        GlobalPointIndexes(trans_cel, s_point_indexes);
        
        TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        TPZMultiphysicsElement * mf_trans_cel = dynamic_cast<TPZMultiphysicsElement * >(trans_cel);
        
#ifdef PZDEBUG
        if (!mf_mixed_cel || !mf_trans_cel) {
            DebugStop();
        }
#endif
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
        int np_mixed_cel = int_points_mixed.NPoints();
        
        const TPZIntPoints & int_points_trans = mf_trans_cel->GetIntegrationRule();
        int np_trans_cel = int_points_trans.NPoints();
        
#ifdef PZDEBUG
        if (np_mixed_cel != p_point_indexes.size() || np_trans_cel != s_point_indexes.size()) {
            DebugStop();
        }
#endif
        
        REAL w;
        TPZManVector<REAL,3> triplet(3,0.0);
        
        REAL detjac;
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;

        REAL p_avg      = 0.0;
        REAL p_avg_n    = 0.0;

        // Integrating pressure
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            int_points_mixed.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            p_avg_n += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure;
            p_avg +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure;
            
        }

        REAL sa      = 0.0;
        REAL sa_n    = 0.0;
        REAL sb      = 0.0;
        REAL sb_n    = 0.0;
        
        // Integrating Saturation
        for (int ip = 0; ip < np_trans_cel; ip++) {
            
            int_points_trans.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            sa_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa_n()/element_measure;
            sb_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb_n()/element_measure;
            sa +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa()/element_measure;
            sb +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb()/element_measure;

        }
        
//        std::cout << "p_avg_n = "<< p_avg_n << std::endl;
//        std::cout << "p_avg = "<< p_avg << std::endl;
//        
//        std::cout << "sa_n = "<< sa_n << std::endl;
//        std::cout << "sa = "<< sa << std::endl;
        
        // Inserting average pressure and saturation in mixed memory
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            if (fSimulationData->IsCurrentStateQ()) {
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa_n(sa_n);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb_n(sb_n);
            }
            else{
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa(sa);
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb(sb);
            }

        }

        // Inserting average pressure in transport memory
        for (int ip = 0; ip < np_trans_cel; ip++) {
            
            if (fSimulationData->IsCurrentStateQ()) {
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sa_n);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sb_n);
            }
            else{
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg(p_avg);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sa);
                trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sb);
            }

        }
        
    }
    
    return;
    
}

/** @brief Reciprocal (mixed <-> transpor) transfer average quantities to integration points of multiphysics meshes over volumetric elements */
void TRMBuildTransfers::Reciprocal_Memory_TransferII(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_trans){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_comp_indexes.size() == 0 ) {
        DebugStop();
    }
    
    if (!cmesh_mf_mixed || !cmesh_mf_trans) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
    
    TPZMaterial * trans_material = cmesh_mf_trans->FindMaterial(rockid);
    TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> * trans_memory = dynamic_cast<TPZMatWithMem<TRMPhaseMemory,TPZDiscontinuousGalerkin> *>(trans_material);
    
    TPZManVector<int64_t,30> p_point_indexes;
    TPZManVector<int64_t,30> s_point_indexes;
    int64_t nvolumes = fmixed_transport_comp_indexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
        TPZGeoEl  * gel = geometry->Element(fmixed_transport_comp_indexes[ivol].first);
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_comp_indexes[ivol].second.first);
        
        // for each transport subelement
        int n_subcels = fmixed_transport_comp_indexes[ivol].second.second.size();
        
        REAL avg_sa      = 0.0;
        REAL avg_sa_n    = 0.0;
        REAL avg_sb      = 0.0;
        REAL avg_sb_n    = 0.0;
        
        for (int icel = 0; icel < n_subcels; icel++) {
            TPZCompEl * trans_cel = cmesh_mf_trans->Element( fmixed_transport_comp_indexes[ivol].second.second[icel]);
            
#ifdef PZDEBUG
            if (!mixed_cel || !trans_cel || !gel) {
                DebugStop();
            }
#endif
            
            REAL element_measure_mixed = DimensionalMeasure(mixed_cel->Reference());
            REAL element_measure_transport = DimensionalMeasure(trans_cel->Reference());
            
            GlobalPointIndexes(mixed_cel, p_point_indexes);
            GlobalPointIndexes(trans_cel, s_point_indexes);
            
            TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
            TPZMultiphysicsElement * mf_trans_cel = dynamic_cast<TPZMultiphysicsElement * >(trans_cel);
            
#ifdef PZDEBUG
            if (!mf_mixed_cel || !mf_trans_cel) {
                DebugStop();
            }
#endif
            
            // Computing the local integration points indexes
            const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
            int np_mixed_cel = int_points_mixed.NPoints();
            
            const TPZIntPoints & int_points_trans = mf_trans_cel->GetIntegrationRule();
            int np_trans_cel = int_points_trans.NPoints();
            
#ifdef PZDEBUG
            if (np_mixed_cel != p_point_indexes.size() || np_trans_cel != s_point_indexes.size()) {
                DebugStop();
            }
#endif
            
            REAL w;
            TPZManVector<REAL,3> triplet(3,0.0);
            
            REAL detjac;
            TPZFMatrix<REAL> jac;
            TPZFMatrix<REAL> axes;
            TPZFMatrix<REAL> jacinv;
            
            REAL p_avg      = 0.0;
            REAL p_avg_n    = 0.0;
            
            // Integrating pressure
            for (int ip = 0; ip < np_mixed_cel; ip++) {
                int_points_mixed.Point(ip, triplet, w);
                gel->Jacobian(triplet, jac, axes, detjac, jacinv);
                
                p_avg_n += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure_mixed;
                p_avg +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure_mixed;
                
            }
            
            TPZGeoEl  * gel_transport = mf_trans_cel->Reference();
            
            REAL sa      = 0.0;
            REAL sa_n    = 0.0;
            REAL sb      = 0.0;
            REAL sb_n    = 0.0;
            
            // Integrating Saturation
            for (int ip = 0; ip < np_trans_cel; ip++) {
                
                int_points_trans.Point(ip, triplet, w);
                gel_transport->Jacobian(triplet, jac, axes, detjac, jacinv);
                
                sa_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa_n()/element_measure_transport;
                sb_n += w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb_n()/element_measure_transport;
                sa +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sa()/element_measure_transport;
                sb +=  w * detjac * trans_memory->GetMemory()[s_point_indexes[ip]].sb()/element_measure_transport;
                
            }

            if (fSimulationData->IsCurrentStateQ()) {
                avg_sa_n += sa_n*element_measure_transport/element_measure_mixed;
                avg_sb_n += sb_n*element_measure_transport/element_measure_mixed;
            }
            else{
                avg_sa += sa*element_measure_transport/element_measure_mixed;
                avg_sb += sb*element_measure_transport/element_measure_mixed;
            }
            
            // Inserting average pressure and saturation in mixed memory
            for (int ip = 0; ip < np_mixed_cel; ip++) {
                if (fSimulationData->IsCurrentStateQ()) {
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa_n(avg_sa_n);
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb_n(avg_sb_n);
                }
                else{
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sa(avg_sa);
                    mixed_memory->GetMemory()[p_point_indexes[ip]].Set_sb(avg_sb);
                }
                
            }
            
            // Inserting average pressure in transport memory
            for (int ip = 0; ip < np_trans_cel; ip++) {
                
                if (fSimulationData->IsCurrentStateQ()) {
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg_n(p_avg_n);
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sa_n);
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa_n(sb_n);
                }
                else{
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_p_avg(p_avg);
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sa);
                    trans_memory->GetMemory()[s_point_indexes[ip]].Set_sa(sb);
                }
                
            }
        }
        
    }
    
    return;
    
}

/** @brief Transfer average pressure to integration points of multiphysics mixed meshes over volumetric elements */
void TRMBuildTransfers::p_avg_Memory_Transfer(TPZCompMesh * cmesh_mf_mixed){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_cindexes.size() == 0 ) {
        DebugStop();
    }
    
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    TPZManVector<int64_t,30> p_point_indexes;
    int64_t nvolumes = fmixed_transport_cindexes.size();
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
        TPZGeoEl  * gel = geometry->Element(fmixed_transport_cindexes[ivol].first);
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_cindexes[ivol].second.first);
        
#ifdef PZDEBUG
        if (!mixed_cel || !gel) {
            DebugStop();
        }
#endif
        
        REAL element_measure = DimensionalMeasure(mixed_cel->Reference());
        
        GlobalPointIndexes(mixed_cel, p_point_indexes);
        TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        
#ifdef PZDEBUG
        if (!mf_mixed_cel) {
            DebugStop();
        }
#endif
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
        int np_mixed_cel = int_points_mixed.NPoints();
        
        
#ifdef PZDEBUG
        if (np_mixed_cel != p_point_indexes.size()) {
            DebugStop();
        }
#endif
        
        int rockid = gel->MaterialId();
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
        TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
        
        REAL w;
        TPZManVector<REAL,3> triplet(3,0.0);
        
        REAL detjac;
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;
        
        REAL p_avg      = 0.0;
        REAL p_avg_n    = 0.0;
        
        // Integrating pressure
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            int_points_mixed.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            p_avg_n +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure;
            p_avg += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure;
            
        }
        
        // Inserting average pressure in mixed memory
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            if (fSimulationData->IsCurrentStateQ()) {
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
            }
            else{
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);

            }
            
        }
        
    }
        
}

/** @brief Transfer average pressure to integration points of multiphysics mixed meshes over volumetric elements */
void TRMBuildTransfers::p_avg_Memory_TransferII(TPZCompMesh * cmesh_mf_mixed){
    
    
#ifdef PZDEBUG
    if ( fmixed_transport_comp_indexes.size() == 0 ) {
        DebugStop();
    }
    
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    
    // For the imat
    int imat = 0;
    int rockid = this->SimulationData()->RawData()->fOmegaIds[imat];
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * mixed_material = cmesh_mf_mixed->FindMaterial(rockid);
    TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> * mixed_memory = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(mixed_material);
    
    TPZManVector<int64_t,30> p_point_indexes;
    int64_t nvolumes = fmixed_transport_comp_indexes.size();
    
    for (int ivol = 0; ivol < nvolumes; ivol++) {
        
        TPZGeoEl  * gel = geometry->Element(fmixed_transport_comp_indexes[ivol].first);
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(fmixed_transport_comp_indexes[ivol].second.first);
        
#ifdef PZDEBUG
        if (!mixed_cel || !gel) {
            DebugStop();
        }
#endif
        
        REAL element_measure = DimensionalMeasure(mixed_cel->Reference());
        
        GlobalPointIndexes(mixed_cel, p_point_indexes);
        TPZMultiphysicsElement * mf_mixed_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        
#ifdef PZDEBUG
        if (!mf_mixed_cel) {
            DebugStop();
        }
#endif
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_mixed = mf_mixed_cel->GetIntegrationRule();
        int np_mixed_cel = int_points_mixed.NPoints();
        
        
#ifdef PZDEBUG
        if (np_mixed_cel != p_point_indexes.size()) {
            DebugStop();
        }
#endif
        
        REAL w;
        TPZManVector<REAL,3> triplet(3,0.0);
        
        REAL detjac;
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;
        
        REAL p_avg      = 0.0;
        REAL p_avg_n    = 0.0;
        
        // Integrating pressure
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            int_points_mixed.Point(ip, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            
            p_avg_n +=  w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p_n()/element_measure;
            p_avg += w * detjac * mixed_memory->GetMemory()[p_point_indexes[ip]].p()/element_measure;
            
        }
        
        // Inserting average pressure in mixed memory
        for (int ip = 0; ip < np_mixed_cel; ip++) {
            if (fSimulationData->IsCurrentStateQ()) {
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg_n(p_avg_n);
            }
            else{
                mixed_memory->GetMemory()[p_point_indexes[ip]].Set_p_avg(p_avg);
                
            }
            
        }
        
    }
    
}


/** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Initialize_un_To_Transport(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    
#ifdef PZDEBUG
    if (!flux_mesh || !transport_mesh) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    int mesh_index = 0;
    
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<int64_t,10> dof_indexes;
    
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;

    
    TPZManVector<int64_t> indices;
    std::pair<int64_t, int64_t> duplet;
    TPZManVector<int,10> face_sides;
    int64_t face_index;
    int64_t n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
        fun_dof_scatter_Gamma.Resize(n_interfaces);
        fun_To_Transport_Gamma.Resize(0, 0);
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
        fun_dof_scatter_gamma.Resize(n_interfaces);
        fun_To_Transport_gamma.Resize(0, 0);
    }

    
    // Block size structue (Gamma or gamma (Inner element interfaces))
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(n_interfaces);
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {


        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }
        
        face_gel = geometry->Element(face_index);
        
        if (!face_gel) {
            DebugStop();
        }
        
        int64_t left_geo_index     = duplet.first;
        int64_t right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        this->ComputeFaceIndex(left_gel,face_sides);
        
//        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
//        
//        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace *> (left_cel);
        
        
        int face_side     = -1;
        int64_t connect_index = -1;

        if(!IdentifyFace(face_side,left_gel,face_gel)){
            std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
            DebugStop();
        }
        else{
            
            for(int ic = 0; ic < intel_vol->NConnects(); ic++){
                connect_index = intel_vol->SideConnectLocId(ic, face_side);
                if (connect_index != -1) {
                    break;
                }
            }
            
            if (connect_index == -1) {
                std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
                DebugStop();
            }
        }

        this->ElementDofFaceIndexes(connect_index, mf_cel, dof_indexes);
        
        int nshapes = left_cel->Connect(connect_index).NShape();
        
        blocks_dimensions[k_face].first = 1;
        blocks_dimensions[k_face].second = nshapes;
        if (IsBoundaryQ) {
            fun_dof_scatter_Gamma[k_face] = dof_indexes;
        }
        else{
            fun_dof_scatter_gamma[k_face] = dof_indexes;
        }

        
    }
    
    // Initialize the matrix
    
    if (IsBoundaryQ) {
        fun_To_Transport_Gamma.Initialize(blocks_dimensions);
    }
    else{
        fun_To_Transport_gamma.Initialize(blocks_dimensions);
    }
    
    
}

/** @brief Initializate  diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Initialize_un_To_TransportII(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    
#ifdef PZDEBUG
    if (!flux_mesh || !transport_mesh) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
//    int mesh_index = 0;
    
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<int64_t,10> dof_indexes;
    int n_shapes;
    
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    
    TPZManVector<int64_t> indices;
    std::pair<int64_t, int64_t> duplet;
    std::pair<int64_t, std::pair< std::pair<int64_t, int64_t> , std::pair<int64_t, int64_t> > >   cint_ctransport_cmixed_duplet;
    TPZManVector<int,10> face_sides;
    int64_t face_index;
    int64_t n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
        fun_dof_scatter_Gamma.Resize(n_interfaces);
        fun_To_Transport_Gamma.Resize(0, 0);
        fcinterface_ctransport_cmixed_indexes_Gamma.resize(n_interfaces);
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
        fun_dof_scatter_gamma.Resize(n_interfaces);
        fun_To_Transport_gamma.Resize(0, 0);
        fcinterface_ctransport_cmixed_indexes_gamma.resize(n_interfaces);
    }

    // Block size structue Gamma and gamma (inner element interfaces)
    TPZVec< std::pair<int64_t, int64_t> > blocks_dimensions(n_interfaces);
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }
        
        face_gel = geometry->Element(face_index);
        
#ifdef PZDEBUG
        if (!face_gel) {
            DebugStop();
        }
#endif
        
        int64_t left_geo_index     = duplet.first;
        int64_t right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        // Identify the father element in mixed mesh that touch or contain the left element in transport mesh
        
        int left_cel_index, father_left_cel_index_index;
        left_cel_index = left_cel->Index();
        TPZCompEl * left_mixed_cel = NULL;
        
        int64_t n_data = fmixed_transport_comp_indexes.size();
        int n_cels;
        for (int i = 0; i < n_data; i++) {
            n_cels = fmixed_transport_comp_indexes[i].second.second.size();
            for (int j = 0; j < n_cels; j++) {
                if (fmixed_transport_comp_indexes[i].second.second[j] == left_cel_index) {
                    father_left_cel_index_index = fmixed_transport_comp_indexes[i].second.first;
                    left_mixed_cel = flux_mesh->Element(father_left_cel_index_index);
                    break;
                }
            }
            if(left_mixed_cel){
                break;
            }
        }
        
        int right_cel_index, father_right_cel_index_index;
        right_cel_index = right_cel->Index();
        TPZCompEl * right_mixed_cel = NULL;
        

        for (int i = 0; i < n_data; i++) {
            n_cels = fmixed_transport_comp_indexes[i].second.second.size();
            for (int j = 0; j < n_cels; j++) {
                if (fmixed_transport_comp_indexes[i].second.second[j] == right_cel_index) {
                    father_right_cel_index_index = fmixed_transport_comp_indexes[i].second.first;
                    right_mixed_cel = flux_mesh->Element(father_right_cel_index_index);
                    break;
                }
            }
            if(right_mixed_cel){
                break;
            }
        }
        
#ifdef PZDEBUG
        
        if (IsBoundaryQ) {
            
            if(!left_mixed_cel){
                DebugStop();
            }
            
        }
        else{
            
            if(!left_mixed_cel || !right_mixed_cel){
                DebugStop();
            }
            
        }

#endif
        
        // End
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // Left based element
        this->ComputeFaceIndex(left_mixed_cel->Reference(),face_sides);
        
        //        TPZCompEl * cel = cmesh_multiphysics->Element(icel);
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_mixed_cel);
//        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
    
        this->ElementDofIndexes(mf_cel, dof_indexes);
        n_shapes = dof_indexes.size();
#ifdef PZDEBUG
        if (dof_indexes.size()==0) {
            DebugStop();
        }
#endif

        
        blocks_dimensions[k_face].first = 1;
        blocks_dimensions[k_face].second = n_shapes;
        
        cint_ctransport_cmixed_duplet.first = face_gel->Reference()->Index();
        cint_ctransport_cmixed_duplet.second.first.first = left_cel->Index();
        cint_ctransport_cmixed_duplet.second.first.second = right_cel->Index();

        
        if (IsBoundaryQ) {
            cint_ctransport_cmixed_duplet.second.second.first = left_mixed_cel->Index();
            cint_ctransport_cmixed_duplet.second.second.second = left_mixed_cel->Index(); // @omar:: left
            fun_dof_scatter_Gamma[k_face] = dof_indexes;
            fcinterface_ctransport_cmixed_indexes_Gamma[k_face] = cint_ctransport_cmixed_duplet;
        }
        else{
            cint_ctransport_cmixed_duplet.second.second.first = left_mixed_cel->Index();
            cint_ctransport_cmixed_duplet.second.second.second = right_mixed_cel->Index();
            fun_dof_scatter_gamma[k_face] = dof_indexes;            
            fcinterface_ctransport_cmixed_indexes_gamma[k_face] = cint_ctransport_cmixed_duplet;
        }
        
    }
    
    // Initialize the matrix
    // Initialize the matrix
    
    if (IsBoundaryQ) {
        fun_To_Transport_Gamma.Initialize(blocks_dimensions);
    }
    else{
        fun_To_Transport_gamma.Initialize(blocks_dimensions);
    }
    
}

/** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Fill_un_To_Transport(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    
    if (fSimulationData->TransporResolution().first) {
        this->Fill_un_To_TransportII(flux_mesh,transport_mesh,IsBoundaryQ);
        return;
    }
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_un_To_Transport(flux_mesh,transport_mesh,IsBoundaryQ);
    
    int mesh_index = 0;
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<int64_t,10> dof_indexes;
    
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    
    TPZManVector<int64_t> indices;
    std::pair<int64_t, int64_t> duplet;
    TPZManVector<int,10> face_sides;
    TPZFMatrix<REAL> normals;
    int64_t face_index;
    int64_t n_interfaces;
    
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    TPZFNMatrix<100,double> block;
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }

        face_gel = geometry->Element(face_index);

        
        if (!face_gel) {
            DebugStop();
        }
        
        int64_t left_geo_index     = duplet.first;
        int64_t right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        
        if(left_gel->HasSubElement() && right_gel->HasSubElement() && face_gel->HasSubElement()){
            DebugStop();
        }
        

        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(left_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        int face_side       = -1;
        int connect_index   = -1;
        
        if(!IdentifyFace(face_side,left_gel,face_gel)){
            std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
            DebugStop();
        }
        else{
            
            for(int ic = 0; ic < intel_vol->NConnects(); ic++){
                connect_index = intel_vol->SideConnectLocId(ic, face_side);
                if (connect_index != -1) {
                    break;
                }
            }
            
            if (connect_index == -1) {
                std::cout << "iRMS Error:: Given Face is not part of the volume element" << std::endl;
                DebugStop();
            }
        }
        
        this->ElementDofFaceIndexes(connect_index, mf_cel, dof_indexes);
        TPZIntPoints *int_points   = left_gel->CreateSideIntegrationRule(face_side, left_cel->GetgOrder());
        
        int npoints = int_points->NPoints();
        int nshapes = left_cel->Connect(connect_index).NShape();
  
#ifdef PZDEBUG
        if (IsBoundaryQ) {
            if (1 != fun_To_Transport_Gamma.GetSizeofBlock(k_face).first || nshapes != fun_To_Transport_Gamma.GetSizeofBlock(k_face).second){
                DebugStop();
            }
        }
        else{
            if (1 != fun_To_Transport_gamma.GetSizeofBlock(k_face).first || nshapes != fun_To_Transport_gamma.GetSizeofBlock(k_face).second){
                DebugStop();
            }
            
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        TPZFNMatrix<100,REAL> phi(nshapes,1,0.0);
        int el_dim = face_gel->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,nshapes,0.0);
        TPZFNMatrix<50,double> block(npoints,nshapes);
        TPZFNMatrix<50,double> block_integral(1,nshapes,0.0);
        
        REAL w;
        TPZManVector<STATE,2> par_duplet(el_dim,0.0);
        REAL ElementMeasure   = DimensionalMeasure(face_gel);

        REAL detjac;
        TPZFMatrix<REAL> jac,axes,jacinv;

        for (int ip = 0; ip < npoints; ip++) {
         
            // Get the vectorial phi
            int_points->Point(ip, par_duplet, w);
            face_gel->Jacobian(par_duplet, jac, axes, detjac, jacinv);
            intel_vol->SideShapeFunction(face_side, par_duplet, phi, dphidxi);
            
            for (int jp = 0; jp < nshapes; jp++) {
                block_integral(0,jp) +=  w * detjac * phi(jp,0)/ElementMeasure;
            }
            
        }
        
        if (IsBoundaryQ) {
            fun_To_Transport_Gamma.SetBlock(k_face, block_integral);
        }
        else{
            fun_To_Transport_gamma.SetBlock(k_face, block_integral);
        }
        

        
    }
    
    
}


/** @brief Initializate diagonal block matrix to transfer average normal flux solution to integrations points of the transport mesh  */
void TRMBuildTransfers::Fill_un_To_TransportII(TPZCompMesh * flux_mesh, TPZCompMesh * transport_mesh, bool IsBoundaryQ){
    
    // It verify the consistency of dynamic_cast and mesh structure and at the end Initialize diagonal matrix blocks
    Initialize_un_To_TransportII(flux_mesh,transport_mesh,IsBoundaryQ);
    
    
    int mesh_index = 0;
    TPZGeoMesh * geometry = flux_mesh->Reference();
    
    //* seeking for total blocks */
    flux_mesh->LoadReferences();
    TPZManVector<int64_t,10> dof_indexes;
    
    TPZCompEl * face_cel;
    TPZCompEl * mixed_cel;
    TPZGeoEl * left_gel;
    TPZGeoEl * right_gel;
    TPZGeoEl * face_gel;
    TPZGeoEl * mixed_gel;
    
    int int_order_interfaces = 1;
    
    TPZManVector<int64_t> indices;
    std::pair<int64_t, int64_t> duplet;
    TPZManVector<int,10> face_sides;
    TPZFMatrix<REAL> normals;
    int64_t face_index;
    int64_t n_interfaces;

    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    TPZFNMatrix<100,double> block;
    
    for (int k_face = 0; k_face < n_interfaces; k_face++) {
        
        if (IsBoundaryQ) {
            face_index  = finterface_g_indexes_Gamma[k_face];
            duplet      = fleft_right_g_indexes_Gamma[k_face];
        }
        else{
            face_index  = finterface_g_indexes_gamma[k_face];
            duplet      = fleft_right_g_indexes_gamma[k_face];
        }
        
        face_gel    = geometry->Element(face_index);
        
#ifdef PZDEBUG
        if (!face_gel) {
            DebugStop();
        }
#endif
        
        face_cel = face_gel->Reference();
        
#ifdef PZDEBUG
        if (!face_cel) {
            DebugStop();
        }
#endif
        
        int64_t left_geo_index     = duplet.first;
        int64_t right_geo_index    = duplet.second;
        
        left_gel    = geometry->Element(left_geo_index);
        right_gel   = geometry->Element(right_geo_index);
        
        TPZCompEl *left_cel = left_gel->Reference();
        TPZCompEl *right_cel = right_gel->Reference();
        
        if (!left_cel || !right_cel) {
            DebugStop();
        }
        

        
        if (IsBoundaryQ) {
            mixed_cel = flux_mesh->Element(fcinterface_ctransport_cmixed_indexes_Gamma[k_face].second.second.first);
            mixed_gel = mixed_cel->Reference();
        }
        else{
            mixed_cel = flux_mesh->Element(fcinterface_ctransport_cmixed_indexes_gamma[k_face].second.second.first);
            mixed_gel = mixed_cel->Reference();
        }
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(mixed_cel);
        TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(mf_cel->Element(mesh_index));
        
        TPZMultiphysicsInterfaceElement * mf_face_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(face_cel);
        
        
        this->ElementDofIndexes(intel_vol, dof_indexes);
        TPZIntPoints *int_points   = face_gel->CreateSideIntegrationRule(face_gel->NSides()-1, int_order_interfaces);
        
        int npoints = int_points->NPoints();
        int nshapes = dof_indexes.size();
        
        // Computing over all integration points of the compuational element mf_cel
        TPZFNMatrix<100,REAL> phi_dot_n(nshapes,1,0.0);
        int face_gel_dim = face_gel->Dimension();
        int intel_vol_dim = intel_vol->Dimension();
        TPZFNMatrix<50,double> block(npoints,nshapes);
        TPZFNMatrix<50,double> block_integral(1,nshapes,0.0);
        
        REAL w;
        TPZManVector<STATE,2> par_duplet(face_gel_dim,0.0);
        TPZManVector<STATE,3> par_mixed_duplet(intel_vol_dim,0.0);
        
        REAL ElementMeasure   = DimensionalMeasure(face_gel);
        
        REAL detjac;
        TPZFMatrix<REAL> jac,axes,jacinv;
        TPZFNMatrix<3,REAL> n;
        TPZVec<int> vectorsides;
        TPZMaterialData data, face_data;
        face_data.fNeedsNormal = true;
        TPZFNMatrix<100,STATE> phi_qs;
        int nphiu,s_i,v_i;
        
        for (int ip = 0; ip < npoints; ip++) {
            
            // Get the vectorial phi
            int_points->Point(ip, par_duplet, w);
            face_gel->Jacobian(par_duplet, jac, axes, detjac, jacinv);

            mf_face_cel->ComputeRequiredData(face_data, par_duplet);
            ComputeTransformation(face_gel, left_gel, mixed_gel, par_duplet, par_mixed_duplet);
            
//            std::cout << "normal = " << face_data.normal <<  std::endl;
            
            intel_vol->InitMaterialData(data);
            intel_vol->ComputeRequiredData(data, par_mixed_duplet);
            
            phi_qs       = data.phi;
            nphiu       = data.fVecShapeIndex.NElements();

            for (int iu = 0; iu < nphiu; iu++)
            {
                
                v_i = data.fVecShapeIndex[iu].first;
                s_i = data.fVecShapeIndex[iu].second;
                
                for (int k = 0; k < intel_vol_dim; k++) {
                        phi_dot_n(iu,0) += 1.0 * phi_qs(s_i,0) * data.fNormalVec(k,v_i) * face_data.normal[k];
                }
                
            }
        
            for (int j = 0; j < nshapes; j++) {
                block_integral(0,j) +=  w * detjac * phi_dot_n(j,0)/ElementMeasure;
            }
            
        }
        
//        std::cout << "face index = " << face_gel->Index() <<  std::endl;
//        std::cout << "mat id = " << face_gel->MaterialId() <<  std::endl;
//        std::cout << "k_face = " << k_face <<  std::endl;
//        std::cout << "dof_indexes = " << dof_indexes <<  std::endl;
//        block_integral.Print(std::cout);
        
        if (IsBoundaryQ) {
            fun_To_Transport_Gamma.SetBlock(k_face, block_integral);
        }
        else{
            fun_To_Transport_gamma.SetBlock(k_face, block_integral);
        }
        
    }
    
    fun_To_Transport_Gamma.Print(" qn to h G = ");
    fun_To_Transport_gamma.Print(" qn to h g = ");
    
}

/** @brief Transfer normal fluxes to integration points of transport meshes */
void TRMBuildTransfers::un_To_Transport_Mesh(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_transport, bool IsBoundaryQ){
  
#ifdef PZDEBUG
    if (!cmesh_flux || !cmesh_transport) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
//#ifdef PZDEBUG
//    if (this->SimulationData()->RawData()->fOmegaIds[0] != 4) {
//        DebugStop();
//    }
//    
//#endif
    
    TPZManVector<int64_t,30>  point_index_trans;
    TPZManVector<int64_t,30>  point_index_l;
    TPZManVector<int64_t,30>  point_index_r;
    
    REAL p_avg_n_l = -1.0;
    REAL p_avg_n_r = -1.0;
    
    cmesh_transport->LoadReferences();
    TPZGeoMesh * geometry = cmesh_transport->Reference();
    int64_t n_interfaces;
    int dimension = geometry->Dimension();
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    if (IsBoundaryQ) {
        
        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int nbc = this->SimulationData()->RawData()->fGammaIds.size();
        
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            int material_id = this->SimulationData()->RawData()->fGammaIds[ibc];
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
            
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  * material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> *>(material);
            
            if (!material_bc_mem) {
                DebugStop();
            }
            
            // Step one
            TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_Gamma.Cols(),1,0.0);
            int64_t pos = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                for(int iflux = 0; iflux < fun_dof_scatter_Gamma[iface].size(); iflux++) {
                    ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_Gamma[iface][iflux],0);
                    pos++;
                }
            }
            
            // Step two
            TPZFMatrix<STATE> un_at_intpoints;
            fun_To_Transport_Gamma.Multiply(ScatterFluxes,un_at_intpoints);
            
            // Step three
            // Trasnfering integrated normal fluxes values
            int counter = 0;
            int i = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                
                TPZGeoEl *gel    = geometry->Element(finterface_g_indexes_Gamma[iface]);
                TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_Gamma[iface].first);
                TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_Gamma[iface].second);

                
#ifdef PZDEBUG
                if (!gel || !gel_l || !gel_r) {
                    DebugStop();
                }
                
                if (gel_l->Dimension() != dimension) {
                    DebugStop();
                }
#endif
                TPZCompEl * mixed_cel_l = gel_l->Reference();
                TPZCompEl * mixed_cel_r = gel_r->Reference();

                
                if(gel->MaterialId() != material_bc_mem->Id()){
                    counter++;
                    continue;
                }
                
                
#ifdef PZDEBUG
                if (!mixed_cel_l || !mixed_cel_r) {
                    DebugStop();
                }
                
#endif
                
                int rock_id = gel_l->MaterialId();
                //  Getting the total integration point of the destination cmesh
                TPZMaterial * rock_material = cmesh_flux->FindMaterial(rock_id);
                TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material);
                
                GlobalPointIndexes(mixed_cel_l, point_index_l);
                
                p_avg_n_l = material_mixe_mem->GetMemory()[point_index_l[0]].p_avg_n();
                material_bc_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
                material_bc_mem->GetMemory()[i].Set_un(un_at_intpoints(counter,0));
                i++;
                counter++;
                
            }
        }

        
    }
    else{
        
        
        int material_id = this->SimulationData()->InterfacesMatId();        
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
        
        TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin>  * material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin> *>(material);
        
        if (!material_mem) {
            DebugStop();
        }

        int np_cmesh = material_mem->GetMemory().NElements();

        // Step one
        TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_gamma.Cols(),1,0.0);
        int64_t pos = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {
            for(int iflux = 0; iflux < fun_dof_scatter_gamma[iface].size(); iflux++) {
                ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_gamma[iface][iflux],0);
                pos++;
            }
        }

        // Step two
        TPZFMatrix<STATE> un_at_intpoints;
        fun_To_Transport_gamma.Multiply(ScatterFluxes,un_at_intpoints);
        
        // Step three
        // Trasnfering integrated normal fluxes values
        for(int64_t i = 0; i < np_cmesh; i++){
            material_mem->GetMemory()[i].Set_un(un_at_intpoints(i,0));
        }

        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int i = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {

            TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_gamma[iface].first);
            TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_gamma[iface].second);
            
#ifdef PZDEBUG
            if (!gel_l || !gel_r) {
                DebugStop();
            }
            
            if (gel_l->Dimension() != dimension || gel_r->Dimension() != dimension) {
                DebugStop();
            }
#endif
            
            TPZCompEl * mixed_cel_l = gel_l->Reference();
            TPZCompEl * mixed_cel_r = gel_r->Reference();
            
#ifdef PZDEBUG
            if (!mixed_cel_l || !mixed_cel_r) {
                DebugStop();
            }
            
#endif
            
            int rock_id_l = gel_l->MaterialId();
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * rock_material_l = cmesh_flux->FindMaterial(rock_id_l);
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem_l = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material_l);
            
            int rock_id_r = gel_r->MaterialId();
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * rock_material_r = cmesh_flux->FindMaterial(rock_id_r);
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem_r = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material_r);

            GlobalPointIndexes(mixed_cel_l, point_index_l);
            GlobalPointIndexes(mixed_cel_r, point_index_r);
            
            p_avg_n_l = material_mixe_mem_l->GetMemory()[point_index_l[0]].p_avg_n();
            p_avg_n_r = material_mixe_mem_r->GetMemory()[point_index_r[0]].p_avg_n();
            
            material_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
            material_mem->GetMemory()[i].Set_p_avg_n_r(p_avg_n_r);
            i++;
            
        }
       
    }

    return;
}


/** @brief Transfer normal fluxes to integration points of transport meshes */
void TRMBuildTransfers::un_To_Transport_MeshII(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_transport, bool IsBoundaryQ){
    
#ifdef PZDEBUG
    if (!cmesh_flux || !cmesh_transport) {
        std::cout << "There is no computational mesh cmesh_multiphysics, cmesh_multiphysics = Null." << std::endl;
        DebugStop();
    }
#endif
    
    TPZManVector<int64_t,30>  point_index_trans;
    TPZManVector<int64_t,30>  point_index_l;
    TPZManVector<int64_t,30>  point_index_r;
    
    REAL p_avg_n_l = -1.0;
    REAL p_avg_n_r = -1.0;
    
    cmesh_transport->LoadReferences();
    TPZGeoMesh * geometry = cmesh_transport->Reference();
    int64_t n_interfaces;
    int dimension = geometry->Dimension();
    if (IsBoundaryQ) {
        n_interfaces = fleft_right_g_indexes_Gamma.size();
    }
    else{
        n_interfaces = fleft_right_g_indexes_gamma.size();
    }
    
    if (IsBoundaryQ) {
        
        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int nbc = this->SimulationData()->RawData()->fGammaIds.size();
        
        for (int ibc = 0; ibc < nbc; ibc++) {
            
            int material_id = this->SimulationData()->RawData()->fGammaIds[ibc];
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
            
            TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  * material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> *>(material);
            
            if (!material_bc_mem) {
                DebugStop();
            }
            
            // Step one
            TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_Gamma.Cols(),1,0.0);
            int64_t pos = 0;
            for (int iface = 0; iface < n_interfaces; iface++) {
                for(int iflux = 0; iflux < fun_dof_scatter_Gamma[iface].size(); iflux++) {
                    ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_Gamma[iface][iflux],0);
                    pos++;
                }
            }
            
            // Step two
            TPZFMatrix<STATE> un_at_intpoints;
            fun_To_Transport_Gamma.Multiply(ScatterFluxes,un_at_intpoints);
            
            // Step three
            // Trasnfering integrated normal fluxes values
            int counter = 0;
            int i = 0;
            int face_g_index, left_mixed_g_index, right_mixed_g_index;
            for (int iface = 0; iface < n_interfaces; iface++) {
                
                //                TPZGeoEl *gel    = geometry->Element(finterface_g_indexes_Gamma[iface]);
                //                TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_Gamma[iface].first);
                //                TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_Gamma[iface].second);
                
                face_g_index = fcinterface_ctransport_cmixed_indexes_Gamma[iface].first;
                left_mixed_g_index = fcinterface_ctransport_cmixed_indexes_Gamma[iface].second.second.first;
                right_mixed_g_index = fcinterface_ctransport_cmixed_indexes_Gamma[iface].second.second.second;
                
                TPZGeoEl *gel    = cmesh_transport->Element(face_g_index)->Reference();
                TPZGeoEl *gel_l  = cmesh_flux->Element(left_mixed_g_index)->Reference();
                TPZGeoEl *gel_r  = cmesh_flux->Element(right_mixed_g_index)->Reference();
                
                
#ifdef PZDEBUG
                if (!gel || !gel_l || !gel_r) {
                    DebugStop();
                }
                
                if (gel_l->Dimension() != dimension) {
                    DebugStop();
                }
#endif
                TPZCompEl * mixed_cel_l = gel_l->Reference();
                TPZCompEl * mixed_cel_r = gel_r->Reference();
                
                
                if(gel->MaterialId() != material_bc_mem->Id()){
                    counter++;
                    continue;
                }
                
                
#ifdef PZDEBUG
                if (!mixed_cel_l || !mixed_cel_r) {
                    DebugStop();
                }
                
#endif
                
                int rock_id = gel_l->MaterialId();
                //  Getting the total integration point of the destination cmesh
                TPZMaterial * rock_material = cmesh_flux->FindMaterial(rock_id);
                TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material);
                
                GlobalPointIndexes(mixed_cel_l, point_index_l);
                
                p_avg_n_l = material_mixe_mem->GetMemory()[point_index_l[0]].p_avg_n();
                material_bc_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
                material_bc_mem->GetMemory()[i].Set_un(un_at_intpoints(counter,0));
                i++;
                counter++;
                
            }
        }
        
        
    }
    else{
        
        
        int material_id = this->SimulationData()->InterfacesMatId();
        //  Getting the total integration point of the destination cmesh
        TPZMaterial * material = cmesh_transport->FindMaterial(material_id);
        
        TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin>  * material_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZDiscontinuousGalerkin> *>(material);
        
        if (!material_mem) {
            DebugStop();
        }
        
        int np_cmesh = material_mem->GetMemory().NElements();
        
        // Step one
        TPZFMatrix<STATE> ScatterFluxes(fun_To_Transport_gamma.Cols(),1,0.0);
        int64_t pos = 0;
        for (int iface = 0; iface < n_interfaces; iface++) {
            for(int iflux = 0; iflux < fun_dof_scatter_gamma[iface].size(); iflux++) {
                ScatterFluxes(pos,0) = cmesh_flux->Solution()(fun_dof_scatter_gamma[iface][iflux],0);
                pos++;
            }
        }
        
        // Step two
        TPZFMatrix<STATE> un_at_intpoints;
        fun_To_Transport_gamma.Multiply(ScatterFluxes,un_at_intpoints);
        
        // Step three
        // Trasnfering integrated normal fluxes values
        for(int64_t i = 0; i < np_cmesh; i++){
            material_mem->GetMemory()[i].Set_un(un_at_intpoints(i,0));
        }
        
        geometry->ResetReference();
        cmesh_flux->LoadReferences();
        int i = 0;
        int left_mixed_g_index, right_mixed_g_index;
        for (int iface = 0; iface < n_interfaces; iface++) {
            
            //            TPZGeoEl *gel_l  = geometry->Element(fleft_right_g_indexes_gamma[iface].first);
            //            TPZGeoEl *gel_r  = geometry->Element(fleft_right_g_indexes_gamma[iface].second);
            
            left_mixed_g_index = fcinterface_ctransport_cmixed_indexes_gamma[iface].second.second.first;
            right_mixed_g_index = fcinterface_ctransport_cmixed_indexes_gamma[iface].second.second.second;
            
            TPZGeoEl *gel_l  = cmesh_flux->Element(left_mixed_g_index)->Reference();
            TPZGeoEl *gel_r  = cmesh_flux->Element(right_mixed_g_index)->Reference();
            
#ifdef PZDEBUG
            if (!gel_l || !gel_r) {
                DebugStop();
            }
            
            if (gel_l->Dimension() != dimension || gel_r->Dimension() != dimension) {
                DebugStop();
            }
#endif
            
            TPZCompEl * mixed_cel_l = gel_l->Reference();
            TPZCompEl * mixed_cel_r = gel_r->Reference();
            
            
#ifdef PZDEBUG
            if (!mixed_cel_l || !mixed_cel_r) {
                DebugStop();
            }
            
#endif
            
            int rock_id_l = gel_l->MaterialId();
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * rock_material_l = cmesh_flux->FindMaterial(rock_id_l);
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem_l = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material_l);
            
            int rock_id_r = gel_r->MaterialId();
            //  Getting the total integration point of the destination cmesh
            TPZMaterial * rock_material_r = cmesh_flux->FindMaterial(rock_id_r);
            TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin>  * material_mixe_mem_r = dynamic_cast<TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> *>(rock_material_r);
            
            GlobalPointIndexes(mixed_cel_l, point_index_l);
            GlobalPointIndexes(mixed_cel_r, point_index_r);
            
            p_avg_n_l = material_mixe_mem_l->GetMemory()[point_index_l[0]].p_avg_n();
            p_avg_n_r = material_mixe_mem_r->GetMemory()[point_index_r[0]].p_avg_n();
            
            material_mem->GetMemory()[i].Set_p_avg_n_l(p_avg_n_l);
            material_mem->GetMemory()[i].Set_p_avg_n_r(p_avg_n_r);
            i++;
            
        }
        
    }
    
    
    return;
    
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Utility Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief Get Global integration point indexes associaded  */
void TRMBuildTransfers::GlobalPointIndexes(TPZCompEl * cel, TPZManVector<int64_t,30> &int_point_indexes){
    
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
void TRMBuildTransfers::GlobalPointIndexesInterface(TPZCompEl * int_cel, TPZManVector<int64_t,30> &int_point_indexes){
    
    TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(int_cel);
    
#ifdef PZDEBUG
    if(!mf_int_cel)
    {
        DebugStop();
    }
#endif
    
    mf_int_cel->GetMemoryIndices(int_point_indexes);
    
}

bool TRMBuildTransfers::IdentifyFace(int &side, TPZGeoEl * vol, TPZGeoEl * face){
    
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
void TRMBuildTransfers::ComputeTransformation(TPZGeoEl * face_gel_origin, TPZGeoEl * gel_origin , TPZGeoEl * gel_target, TPZVec<REAL> & origin, TPZVec<REAL> & target){
    
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
void TRMBuildTransfers::ComputeFaceIndex(TPZGeoEl * gel , TPZVec<int> &sides){
    
    
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
void TRMBuildTransfers::ComputeFaceNormals(TPZGeoEl * gel , TPZVec<int> &sides, TPZFMatrix<STATE> &normals){
    
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
void TRMBuildTransfers::ComputeLeftRight(TPZCompMesh * transport_mesh){
    
    if (fSimulationData->TransporResolution().first) {
        this->ComputeLeftRightII(transport_mesh);
        return;
    }
    
    fleft_right_g_indexes_Gamma.resize(0);
    fleft_right_g_indexes_gamma.resize(0);
    
#ifdef PZDEBUG
    if (!transport_mesh) {
        std::cout << "There is no computational transport mesh, transport_mesh = Null." << std::endl;
        DebugStop();
    }
#endif
    
    int64_t nel = transport_mesh->NElements();
    int64_t face_index;
    std::pair <int64_t,int64_t> duplet;
    transport_mesh->LoadReferences();
    int dimension = transport_mesh->Reference()->Dimension();
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = transport_mesh->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsInterfaceElement * interface = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (!interface) {
            continue;
        }
        
        TPZCompEl * left_cel = interface->LeftElement();
        TPZCompEl * right_cel = interface->RightElement();
        
#ifdef PZDEBUG
        
        if(!left_cel || !right_cel){
            DebugStop();
        }
#endif
        
        if(interface->Reference()->HasSubElement()) {
            continue;
        }
        
        face_index  = interface->Reference()->Index();
        duplet      = std::make_pair(left_cel->Reference()->Index(), right_cel->Reference()->Index());
        
        if(left_cel->Dimension() != dimension ||  right_cel->Dimension() != dimension){
            
            fleft_right_g_indexes_Gamma.push_back(duplet);
            finterface_g_indexes_Gamma.push_back(face_index);
            continue;
        }
        
        fleft_right_g_indexes_gamma.push_back(duplet);
        finterface_g_indexes_gamma.push_back(face_index);
        
    }
    
//    std::cout << " on Gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_indexes_Gamma.size(); k++) {
//        std::cout << " volume k : " << k << std::endl;
//        std::cout << " volume left : " << fleft_right_g_indexes_Gamma[k].first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_indexes_Gamma[k].second <<std::endl;
//    }
//    
//    std::cout << " on gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_indexes_gamma.size(); k++) {
//        std::cout << " volume k : " << k << std::endl;
//        std::cout << " volume left : " << fleft_right_g_indexes_gamma[k].first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_indexes_gamma[k].second <<std::endl;
//    }
    
#ifdef PZDEBUG
    if (finterface_g_indexes_Gamma.size() == 0) {
        DebugStop();
    }
    if (finterface_g_indexes_gamma.size() == 0) {
        std::cout << "Warning:: No inner interfaces were found" << std::endl;
    }
#endif
    
}


/** @brief Compute left and right geometric element indexes */
void TRMBuildTransfers::ComputeLeftRightII(TPZCompMesh * transport_mesh){
    
    
    fleft_right_g_indexes_Gamma.resize(0);
    fleft_right_g_indexes_gamma.resize(0);
    
#ifdef PZDEBUG
    if (!transport_mesh) {
        std::cout << "There is no computational transport mesh, transport_mesh = Null." << std::endl;
        DebugStop();
    }
#endif
    
    int64_t nel = transport_mesh->NElements();
    int64_t face_index;
    std::pair <int64_t,int64_t> duplet;
    transport_mesh->LoadReferences();
    int dimension = transport_mesh->Reference()->Dimension();
    for (int64_t icel = 0; icel < nel; icel++) {
        
        TPZCompEl * cel = transport_mesh->Element(icel);
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsInterfaceElement * interface = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel);
        
        if (!interface) {
            continue;
        }
        
        TPZCompEl * left_cel = interface->LeftElement();
        TPZCompEl * right_cel = interface->RightElement();
        
#ifdef PZDEBUG
        
        if(!left_cel || !right_cel){
            DebugStop();
        }
#endif
        
        if(interface->Reference()->HasSubElement()) {
            continue;
        }
        
        face_index  = interface->Reference()->Index();
        duplet      = std::make_pair(left_cel->Reference()->Index(), right_cel->Reference()->Index());
        
        if(left_cel->Dimension() != dimension ||  right_cel->Dimension() != dimension){
            
            fleft_right_g_indexes_Gamma.push_back(duplet);
            finterface_g_indexes_Gamma.push_back(face_index);
            continue;
        }
        
        fleft_right_g_indexes_gamma.push_back(duplet);
        finterface_g_indexes_gamma.push_back(face_index);
        
    }
    
    
//    std::cout << " on Gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_indexes_Gamma.size(); k++) {
//        std::cout << " volume k : " << k << std::endl;
//        std::cout << " volume left : " << fleft_right_g_indexes_Gamma[k].first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_indexes_Gamma[k].second <<std::endl;
//    }
//    
//    std::cout << " on gamma " << std::endl;
//    for (int k = 0; k < fleft_right_g_indexes_gamma.size(); k++) {
//        std::cout << " volume k : " << k << std::endl;
//        std::cout << " volume left : " << fleft_right_g_indexes_gamma[k].first << std::endl;
//        std::cout << " volume ritgh : " << fleft_right_g_indexes_gamma[k].second <<std::endl;
//    }
    
#ifdef PZDEBUG
    if (finterface_g_indexes_Gamma.size() == 0) {
        DebugStop();
    }
    if (finterface_g_indexes_gamma.size() == 0) {
        std::cout << "Warning:: No inner interfaces were found" << std::endl;
    }
#endif
    
}


/** @brief Dimensionla Measure of the elemnt */
REAL TRMBuildTransfers::DimensionalMeasure(TPZGeoEl * gel){
    
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


void TRMBuildTransfers::ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<int64_t> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<int64_t> index(0,0);
    int nconnect = intel->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        int64_t seqnumber = con.SequenceNumber();
        int64_t position = intel->Mesh()->Block().Position(seqnumber);
        int nshape = con.NShape();
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
        }
    }
    
    dof_indexes = index;
    return;
}

void TRMBuildTransfers::ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<int64_t> &dof_indexes){
    
    
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
    
//    TPZStack<int64_t> index(0,0);
//    int nconnect = intel_vol->NConnects();
//    for (int icon = 0; icon < nconnect; icon++) {
//        TPZConnect  & con = m_el->Connect(icon);
//        int64_t seqnumber = con.SequenceNumber();
//        int64_t position = m_el->Mesh()->Block().Position(seqnumber);
//        int nshape = con.NShape();
//        for (int ish=0; ish < nshape; ish++) {
//            index.Push(position+ ish);
//        }
//    }
//    
//    dof_indexes = index;
//    return;
    
    
    TPZStack<int64_t> index(0,0);
    int nconnect = intel_vol->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = m_el->Connect(icon);
        int64_t seqnumber = con.SequenceNumber();
        int64_t position = m_el->Mesh()->Block().Position(seqnumber);
        int b_size = m_el->Mesh()->Block().Size(seqnumber);
        for (int ib=0; ib < b_size; ib++) {
            index.Push(position+ ib);
        }
    }
    
    dof_indexes = index;
    return;

}

void TRMBuildTransfers::ElementDofFaceIndexes(int connect_index,TPZInterpolationSpace * &intel, TPZVec<int64_t> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<int64_t> index(0,0);
    TPZConnect  & con = intel->Connect(connect_index);
    int64_t seqnumber = con.SequenceNumber();
    int64_t position = intel->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}

void TRMBuildTransfers::ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<int64_t> &dof_indexes, int el_index){
    
#ifdef PZDEBUG
    if (!m_el) {
        DebugStop();
    }
#endif
    
    TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(m_el->Element(el_index));
    
#ifdef PZDEBUG
    if (!intel_vol) {
        DebugStop();
    }
#endif
    
    int start = 0;
    int end = intel_vol->NConnects();
    if (el_index == 1) {
        
        TPZInterpolationSpace * intel_vol_q = dynamic_cast<TPZInterpolationSpace * >(m_el->Element(0));
        
#ifdef PZDEBUG
        if (!intel_vol_q) {
            DebugStop();
        }
#endif
        start = intel_vol_q->NConnects();
        end = m_el->NConnects();
    }
    
    TPZStack<int64_t> index(0,0);
    int nconnect = end;
    for (int icon = start; icon < nconnect; icon++) {
        TPZConnect  & con = m_el->Connect(icon);
        int64_t seqnumber = con.SequenceNumber();
        int64_t position = m_el->Mesh()->Block().Position(seqnumber);
        int b_size = m_el->Mesh()->Block().Size(seqnumber);
        for (int ib=0; ib < b_size; ib++) {
            index.Push(position+ ib);
        }
    }
    
    dof_indexes = index;
    return;    
}

void TRMBuildTransfers::ElementDofFaceIndexes(int connect_index, TPZMultiphysicsElement * &m_el, TPZVec<int64_t> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!m_el && connect_index > 4) {
        DebugStop();
    }
#endif
    

    TPZStack<int64_t> index(0,0);
    TPZConnect  & con = m_el->Connect(connect_index);
    int64_t seqnumber = con.SequenceNumber();
    int64_t position = m_el->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}


/** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
void TRMBuildTransfers::FillComputationalElPairs(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport){
    
    
    if (fSimulationData->TransporResolution().first) {
        this->FillComputationalElPairsII(cmesh_mf_mixed,cmesh_mf_transport);
        return;
    }
    else{
        DebugStop();
    }

    fmixed_transport_cindexes.resize(0);
    
#ifdef PZDEBUG
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
    
    if (!fSimulationData->IsOnePhaseQ() && !cmesh_mf_transport) {
        DebugStop();
    }
    
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<int64_t, std::pair <int64_t,int64_t> > gel_indexes;
    
    for (int64_t i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dimension) {
            continue;
        }
        int mat_id = gel->MaterialId();
        bool filter_check = (mat_id == 5) || (mat_id == 6) || (mat_id == 7);
        if (!filter_check) {
            continue;
        }
        
        int parabolic_res = fSimulationData->MHMResolution().second.second;
        if (gel->Level() != parabolic_res) {
            continue;
        }
        
        gel_indexes.first = gel->Index();
        gel_indexes.second.first = -1;
        gel_indexes.second.second = -1;
        fmixed_transport_cindexes.push_back(gel_indexes);

    }
    
    // counting volumetric elements
    int64_t nvol_elements = fmixed_transport_cindexes.size();
    fmixed_transport_cindexes.resize(nvol_elements);
    
    // inserting mixed indexes
    cmesh_mf_mixed->LoadReferences();
    for(int64_t ivol = 0; ivol < nvol_elements; ivol++){

        TPZCompEl * mixed_cel = geometry->Element(fmixed_transport_cindexes[ivol].first)->Reference();
        
#ifdef PZDEBUG
        if (!mixed_cel) {
            //continue;
            DebugStop();
        }
#endif
        fmixed_transport_cindexes[ivol].second.first = mixed_cel->Index();
        
    }
    
    if(fSimulationData->IsOnePhaseQ()){
        
//        for (int k = 0; k < fmixed_transport_cindexes.size(); k++) {
//            std::cout << " volume k : " << k <<std::endl;
//            std::cout << " volume gel : " << fmixed_transport_cindexes[k].first <<std::endl;
//            std::cout << " volume cmixed : " << fmixed_transport_cindexes[k].second.first <<std::endl;
//            std::cout << " volume ctransport : " << fmixed_transport_cindexes[k].second.second <<std::endl;
//        }
        
        return;
    }
    
    // inserting transport indexes
    cmesh_mf_transport->LoadReferences();
    for(int64_t ivol = 0; ivol < nvol_elements; ivol++){
        
        TPZCompEl * trans_cel = geometry->Element(fmixed_transport_cindexes[ivol].first)->Reference();
        
#ifdef PZDEBUG
        if (!trans_cel) {
            DebugStop();
        }
#endif
        fmixed_transport_cindexes[ivol].second.second = trans_cel->Index();
        
    }
    
//    for (int k = 0; k < fmixed_transport_cindexes.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << fmixed_transport_cindexes[k].first <<std::endl;
//        std::cout << " volume cmixed : " << fmixed_transport_cindexes[k].second.first <<std::endl;
//        std::cout << " volume ctransport : " << fmixed_transport_cindexes[k].second.second <<std::endl;
//    }
    
}


/** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
void TRMBuildTransfers::FillComputationalElPairsII(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport){
    
    fmixed_transport_comp_indexes.resize(0);
    
#ifdef PZDEBUG
    if (!cmesh_mf_mixed) {
        DebugStop();
    }
    
    if (!fSimulationData->IsOnePhaseQ() && !cmesh_mf_transport) {
        DebugStop();
    }
    
#endif
    
    cmesh_mf_mixed->LoadReferences();
    TPZGeoMesh * geometry = cmesh_mf_mixed->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    std::pair<int64_t, std::pair <int64_t, std::vector<int64_t> > > gel_cel_indexes;
    
    for (int64_t icel = 0; icel < cmesh_mf_mixed->NElements(); icel++) {
        
        TPZCompEl * mixed_cel = cmesh_mf_mixed->Element(icel);
        
#ifdef PZDEBUG
        if (!mixed_cel) {
            DebugStop();
        }
#endif
        TPZGeoEl * gel = mixed_cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dimension) {
            continue;
        }
        int mat_id = gel->MaterialId();
        bool filter_check = (mat_id == 5) || (mat_id == 6) || (mat_id == 7);
        if (!filter_check) {
            continue;
        }
        
        int parabolic_res = fSimulationData->MHMResolution().second.second;
        if (gel->Level() != parabolic_res) {
            continue;
        }
        
        gel_cel_indexes.first = gel->Index();
        gel_cel_indexes.second.first = mixed_cel->Index();
        gel_cel_indexes.second.second.resize(0);
        fmixed_transport_comp_indexes.push_back(gel_cel_indexes);
        
    }
    
    // counting volumetric elements
    int nvol_elements = fmixed_transport_comp_indexes.size();
    
    if(fSimulationData->IsOnePhaseQ()){
        
//        for (int k = 0; k < fmixed_transport_comp_indexes.size(); k++) {
//            std::cout << " volume k : " << k <<std::endl;
//            std::cout << " volume gel : " << fmixed_transport_comp_indexes[k].first <<std::endl;
//            std::cout << " volume cmixed : " << fmixed_transport_comp_indexes[k].second.first <<std::endl;
//        }
        
        return;
    }
    
    // inserting transport indexes
    cmesh_mf_transport->LoadReferences();
    TPZVec<TPZGeoEl *> n_refined_sons;
    int cel_index;
    for(int64_t ivol = 0; ivol < nvol_elements; ivol++){
        
        TPZGeoEl * father_gel = geometry->Element(fmixed_transport_comp_indexes[ivol].first);
        
#ifdef PZDEBUG
        if (!father_gel) {
            DebugStop();
        }
#endif
        n_refined_sons.resize(0);
        father_gel->GetHigherSubElements(n_refined_sons);
        
        if (n_refined_sons.size() == 0) {
            n_refined_sons.resize(1);
            n_refined_sons[0] = father_gel;
        }
        
        for (int igel = 0; igel < n_refined_sons.size(); igel++) {
            cel_index = geometry->Element(n_refined_sons[igel]->Index())->Reference()->Index();
            fmixed_transport_comp_indexes[ivol].second.second.push_back(cel_index);
        }
        
    }
    
//    for (int k = 0; k < fmixed_transport_comp_indexes.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << fmixed_transport_comp_indexes[k].first <<std::endl;
//        std::cout << " volume cmixed : " << fmixed_transport_comp_indexes[k].second.first <<std::endl;
//        for (int igel = 0; igel < n_refined_sons.size(); igel++) {
//            int index = fmixed_transport_comp_indexes[k].second.second[igel]; ;
//            std::cout << " volume ctransport : " << index <<std::endl;
//        }
//    }
    
}

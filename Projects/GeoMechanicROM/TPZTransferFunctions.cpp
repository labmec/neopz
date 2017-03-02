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


void TPZTransferFunctions::Initialize_phi_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics || !fCmeshRB_projections) {
        std::cout << "There is no computational meshes, cmesh_multiphysics||fCmeshRB_projections = Null." << std::endl;
        DebugStop();
    }
#endif
    
    
    int material_id = 1; // rock target
    
    fCmeshRB_projections->LoadReferences();
    
    long nel = fCmeshRB_projections->NElements();
    int n_var_dim = fCmeshRB_projections->Reference()->Dimension(); // vectorial
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fphi_u_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions(nel);
    
    
    for (long icel = 0; icel < nel; icel++) {
        
        fCmeshRB_projections->LoadReferences();
        TPZCompEl * cel_GP = fCmeshRB_projections->Element(icel);
        
#ifdef PZDEBUG
        if (!cel_GP) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel_GP->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
//        if (material_id != gel->MaterialId()) {
//            continue;
//        }
        
        gel->ResetReference();
        cmesh_multiphysics->LoadReferences();
        TPZCompEl * cel = gel->Reference();
        
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
        element_index = cel_GP->Index();
        
        TPZInterpolationSpace * intel_GP = dynamic_cast<TPZInterpolationSpace * >(cel_GP);
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        if(intel_GP->Dimension() < n_var_dim){

            mf_cel->GetMemoryIndices(int_point_indexes);
            this->ElementDofIndexes(intel_GP, dof_indexes);
            fphi_u_dof_scatter[element_index] = dof_indexes;
            blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
            blocks_dimensions[element_index].second = dof_indexes.size();
            continue;
        }
        
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(intel_GP, dof_indexes);
        fphi_u_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions[element_index].second = dof_indexes.size();
    }
    
    // Initialize the matrix
    fphi_u_To_Geomechanic.Initialize(blocks_dimensions);
    
}


void TPZTransferFunctions::Fill_phi_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_phi_u_To_Mixed(cmesh_multiphysics, mesh_index);
    

    int material_id = 1; // rock target
    
    fCmeshRB_projections->LoadReferences();
    
    long nel = fCmeshRB_projections->NElements();
    int n_var_dim = fCmeshRB_projections->Reference()->Dimension(); // vectorial
    long element_index = 0;
    
    std::pair<long, long> block_dim;
    TPZMaterialData data;
    
    for (long icel = 0; icel < nel; icel++) {
        
        fCmeshRB_projections->LoadReferences();
        TPZCompEl * cel_GP = fCmeshRB_projections->Element(icel);
        
#ifdef PZDEBUG
        if (!cel_GP) {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel_GP->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        gel->ResetReference();
        cmesh_multiphysics->LoadReferences();
        TPZCompEl * cel = gel->Reference();
        
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
        element_index = cel_GP->Index();
        TPZInterpolationSpace * intel_GP = dynamic_cast<TPZInterpolationSpace * >(cel_GP);
        
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
        TPZFNMatrix<100,REAL> phi(intel_GP->NShapeF(),1,0.0);
        int el_dim = mf_cel->Reference()->Dimension();
        TPZFNMatrix<300,REAL> dphidxi(el_dim,intel_GP->NShapeF(),0.0);
        TPZFMatrix<double> block;
        
        block.Resize(block_dim.first*n_var_dim,block_dim.second);
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(el_dim,0.0);
            STATE w;
            int_points_geomechanic.Point(ip, qsi, w);
            // Get the phi for H1 elasticity
            intel_GP->Shape(qsi, phi, dphidxi);
            intel_GP->InitMaterialData(data);
            intel_GP->ComputeRequiredData(data,qsi);
            
#ifdef PZDEBUG
            if(block_dim.second != phi.Rows() * n_var_dim){
                DebugStop();
            }
#endif
            
            for (int id = 0; id < n_var_dim; id++) {
                for (int jp = 0; jp < phi.Rows(); jp++) {
                    block(ip*n_var_dim+id,n_var_dim*jp+id) = phi(jp%n_var_dim,0); // or phi(jp,0); check it
                }
            }
            
        }
        
        fphi_u_To_Geomechanic.SetBlock(element_index, block);
    }
    
    fphi_u_To_Geomechanic.Print(" phi_u = ");
    
}

void TPZTransferFunctions::phi_u_To_Geomechanic_Memory(TPZCompMesh * cmesh_multiphysics){
    
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics) {
        DebugStop();
    }
#endif
    
    int nel = fCmeshRB_projections->NElements();
    int dim = cmesh_multiphysics->Dimension();
    int n_rb = fCmeshRB_projections->Solution().Cols();
    
    fCmeshRB_projections->Solution().Print("GP = ");
    
    // For the imat
    int imat = 0;
    int rockid = 1;
    
    //  Getting the total integration point of the destination cmesh
    TPZMaterial * material = cmesh_multiphysics->FindMaterial(rockid);
    TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> *>(material);
    int np_cmesh = associated_material->GetMemory().NElements();
    
    // Step zero
    for(long i = 0; i <  np_cmesh; i++){
        associated_material->GetMemory()[i].phi_u().Resize(n_rb);
    }
    
//    for all RB functions
    for(int i_rb = 0; i_rb < n_rb; i_rb++){
    
        // Step one
        TPZFMatrix<STATE> ScatterDisplacements(fphi_u_To_Geomechanic.Cols(),1,0.0);
        long pos = 0;
        for (int el = 0; el < nel; el++) {
            for(int ip = 0; ip < fphi_u_dof_scatter[el].size(); ip++) {
                ScatterDisplacements(pos,0) = fCmeshRB_projections->Solution()(fphi_u_dof_scatter[el][ip],i_rb);
                pos++;
            }
            std::cout << "el dof = " << fphi_u_dof_scatter[el] <<std::endl;
        }
        
        ScatterDisplacements.Print(" u  at points ");
        // Step two
        TPZFMatrix<STATE> u_at_intpoints;
        fphi_u_To_Geomechanic.Multiply(ScatterDisplacements,u_at_intpoints);
        
        for (int icel = 0; icel < nel; icel++) {
            fCmeshRB_projections->LoadReferences();
            TPZCompEl *cel_GP = fCmeshRB_projections->Element(icel);
#ifdef PZDEBUG
            if (!cel_GP) {
                DebugStop();
            }
#endif
            TPZGeoEl * gel = cel_GP->Reference();
            
#ifdef PZDEBUG
            if (!gel) {
                DebugStop();
            }
#endif
            
            if (gel->Dimension() != dim) {
                continue;
            }
            
            TPZManVector<long, 30> int_point_indexes;
            cmesh_multiphysics->LoadReferences();
            TPZCompEl * mf_cel = gel->Reference();
            GlobalPointIndexes(mf_cel, int_point_indexes);
            

            // Trasnfering values
            int n_points = int_point_indexes.size();
            
            TPZManVector<STATE,3> phi_u(dim,0.0);
            for(long ip = 0; ip <  n_points; ip++){
                for (int id = 0; id < dim ; id++) {
                    phi_u[id]= u_at_intpoints(int_point_indexes[ip]+id,0);
                }
                associated_material->GetMemory()[int_point_indexes[ip]].Set_phi_u_n(i_rb, phi_u);
            }
            
        }


        
//        // Step two
//        TPZFMatrix<STATE> u_at_intpoints;
//        fphi_u_To_Geomechanic.Multiply(ScatterDisplacements,u_at_intpoints);
//        // Trasnfering values
//        TPZManVector<STATE,3> phi_u(dim,0.0);
//        for(long i = 0; i <  np_cmesh; i++){
//            for (int id = 0; id < dim ; id++) {
//                phi_u[id]= u_at_intpoints(i*dim+id,0);
//            }
//            associated_material->GetMemory()[i].Set_phi_u_n(i_rb, phi_u);
//        }
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
        int nvars = intel->Mesh()->Block().Size(seqnumber); // @omar:: must Verify
        for (int ish=0; ish < nvars; ish++) {
                index.Push(position + ish);
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
        int nshape = con.NShape();
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
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


/** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
void TPZTransferFunctions::FillComputationalElPairs(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport){
    
    
    DebugStop();
    
}


/** @brief Compute compuational mesh pair (mixed, transport) indexed by geometric volumetic element index */
void TPZTransferFunctions::FillComputationalElPairsII(TPZCompMesh * cmesh_mf_mixed, TPZCompMesh * cmesh_mf_transport){
    
    DebugStop();
}

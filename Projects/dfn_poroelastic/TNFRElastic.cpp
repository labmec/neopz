//
//  TNFRElastic.cpp
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#include "TNFRElastic.h"
#include "TNRFElasticMemory.h"
#include "pzbndcond.h"
#include "pzgeoelbc.h"
#include "pzintel.h"
#include "TPZMatWithMem.h"
//#include "TPZVTKGeoMesh.h"

// Material for construction of fracture's shape functions
#include "pzmat1dlin.h"

// Constructor
TNFRElastic::TNFRElastic(){
    
    m_geometry  = nullptr;
    m_cmesh     = nullptr;
    m_Operator  = nullptr;
    m_boundary_data.clear();
    m_n_threads = 0;
}

// Descructor
TNFRElastic::~TNFRElastic(){
    
}

void TNFRElastic::SetBoundaryData(std::vector<TNFRBoundaryDescription> boundary_data){
    m_boundary_data = boundary_data;
}

void TNFRElastic::SetNumberOfThreads(unsigned int n_threads){
    m_n_threads = n_threads;
}

void TNFRElastic::SetEnableBandwidthReduction(bool bandwidth_reduction_Q){
    m_bandwidth_reduction_Q = bandwidth_reduction_Q;
}

bool TNFRElastic::BuildOperator(TPZGeoMesh *  geometry, int p_order){
    
    bool CoherenceQ = true;
    m_geometry = geometry;
    m_cmesh = new TPZCompMesh(m_geometry);
    m_cmesh->SetDefaultOrder(p_order);
    m_cmesh->SetDimModel(m_geometry->Dimension());
    m_cmesh->SetAllCreateFunctionsContinuous();
    
    // Inserting a single material and related boundary conditions
    int rock_id = 1;
    TNFRElasticMaterial * rock_material = new TNFRElasticMaterial(rock_id);
    
    // Lamé parameters
    REAL lambda = 32e3;
    REAL mu     = 32e3;
    rock_material->SetLameParameters(lambda, mu);
    
    m_cmesh->InsertMaterialObject(rock_material);

    unsigned int n_bc = m_boundary_data.size();
    
    if (n_bc == 0) {
        std::cout << "A_e:: No bounary data exists." << std::endl;
        CoherenceQ = false;
    }
    
    TPZFMatrix<REAL> val1(1,1,0.0),val2(1,1,0.0);
    for (unsigned int i = 0; i < n_bc; i++) {
        TNFRBoundaryDescription bc = m_boundary_data[i];
        unsigned int n_values = bc.GetBCValues().size();
        
        if (n_values == 0) {
            std::cout << "A_e:: Bounary without values exists." << std::endl;
            CoherenceQ = false;
        }
        
        val2.Resize(n_values, 1);
        for (unsigned int j = 0; j < n_values; j++) {
            val2(j,0) = bc.GetBCValues()[j];
        }
        
        TPZMatWithMem<TNRFElasticMemory,TPZBndCond> * boundary_c = new TPZMatWithMem<TNRFElasticMemory,TPZBndCond>;
        boundary_c->SetNumLoadCases(1);
        boundary_c->SetMaterial(rock_material);
        boundary_c->SetId(bc.GetBCId());
        boundary_c->SetType(bc.GetBCType());
        boundary_c->SetValues(val1, val2);
        m_cmesh->InsertMaterialObject(boundary_c);
    }
    
    m_cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("Cmesh_Deformation.txt");
    m_cmesh->Print(out);
#endif
    
    bool fracture_detected_Q = false;
    // Detecting fractures
    unsigned int n_elements = m_geometry->NElements();
    for (unsigned int igel = 0; igel < n_elements; igel++) {
        TPZGeoEl * gel = geometry->Element(igel);
        if (gel->MaterialId() != m_fracture_bc_id) { // Filtering by fracture material id
            continue;
        }
        fracture_detected_Q = true;
    }
    
    if (fracture_detected_Q) {
        
        // Auxiliary objects
        std::vector<REAL> tensor_values(6,0.0);
        REAL wb_p = -100.0;
        tensor_values[0] = wb_p;
        tensor_values[3] = wb_p;
        tensor_values[5] = wb_p;
        
        // Fracture bc data
        TNFRBoundaryDescription Fracture_bc;
        Fracture_bc.SetBCId(9);
        Fracture_bc.SetBCType(2);
        Fracture_bc.SetBCValues(tensor_values);
        
        unsigned int n_values = Fracture_bc.GetBCValues().size();
        TPZFMatrix<REAL> val1(1,1,0.0),val2(n_values,1,0.0);
        for (unsigned int j = 0; j < n_values; j++) {
            val2(j,0) = Fracture_bc.GetBCValues()[j];
        }
        
        TPZMatWithMem<TNRFElasticMemory,TPZBndCond> * boundary_c_fractures = new TPZMatWithMem<TNRFElasticMemory,TPZBndCond>;
        boundary_c_fractures->SetNumLoadCases(1);
        boundary_c_fractures->SetMaterial(rock_material);
        boundary_c_fractures->SetId(Fracture_bc.GetBCId());
        boundary_c_fractures->SetType(Fracture_bc.GetBCType());
        boundary_c_fractures->SetValues(val1, val2);
        m_cmesh->InsertMaterialObject(boundary_c_fractures);
        
        InsertFractureRepresentation();
        
#ifdef PZDEBUG
        std::ofstream out("Cmesh_Deformation_with_Fractures.txt");
        m_cmesh->Print(out);
#endif
    }
    
    // Creating the analysis object with pardiso solver
    m_Operator = new TPZAnalysis;
    m_Operator->SetCompMesh(m_cmesh, m_bandwidth_reduction_Q);
    TPZSymetricSpStructMatrix strmat(m_cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    strmat.SetNumThreads(m_n_threads);
    m_Operator->SetSolver(step);
    m_Operator->SetStructuralMatrix(strmat);
    
    if (CoherenceQ) {
        std::cout << "A_e:: Operator was constructed properly." << std::endl;
    }
    return CoherenceQ;
    
}


void TNFRElastic::InsertFractureRepresentation(){
    
    // Break mesh by inserting fractures as boundary elements
    unsigned int n_elements = m_geometry->NElements();
    for (unsigned int igel = 0; igel < n_elements; igel++) {
        TPZGeoEl * gel = m_geometry->Element(igel);
        if (gel->MaterialId() != 9) { // Filtering by fracture material id
            continue;
        }
        
        unsigned int nsides = gel->NSides();
        TPZStack<TPZCompElSide> neighbors;
        TPZGeoElSide gelside(gel,nsides-1);
        
        gelside.EqualLevelCompElementList(neighbors, 0, 0);
        
        if(neighbors.size()!=2){ // Fractures should have two volumetric elements
            DebugStop();
        }
        
        gel->ResetReference();
        neighbors[0].Element()->Reference()->ResetReference(); // +
        neighbors[1].Element()->Reference()->ResetReference(); // -
        
        // Working with + side
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(neighbors[0].Element());
            if(!intel){
                DebugStop();
            }
            intel->LoadElementReference();
            
            std::set<int64_t> conner_indexes;
            intel->BuildCornerConnectList(conner_indexes);
            
            // Midside connect
            int local_index = intel->MidSideConnectLocId(neighbors[0].Side());
            TPZConnect & inner_side_connect = intel->MidSideConnect(neighbors[0].Side());
            if(inner_side_connect.NElConnected() != 2)
            {
                DebugStop();
            }
            int64_t index = CreateAndDuplicateConnect(intel, local_index, inner_side_connect);
            
            
            // The fracture boundary with + side
            TPZGeoElBC bc(intel->Reference(),neighbors[0].Side(),m_fracture_bc_id);
            m_cmesh->CreateCompEl(bc.CreatedElement(), index);
            TPZCompEl *var = m_cmesh->Element(index);
            var->Reference()->ResetReference();
            intel->Reference()->ResetReference();
            
        }
        
        // Working with - side
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(neighbors[1].Element());
            
            if(!intel){
                DebugStop();
            }
            
            intel->LoadElementReference();
            //            intel->SetSideOrient(neighbors[1].Side(), 1);
            TPZGeoElBC bc(intel->Reference(),neighbors[1].Side(),m_fracture_bc_id);
            
            int64_t index;
            m_cmesh->CreateCompEl(bc.CreatedElement(), index);
            TPZCompEl *var = m_cmesh->Element(index);
            var->Reference()->ResetReference();
            intel->Reference()->ResetReference();
            
        }
        
    }
    
    m_cmesh->ExpandSolution();
    
//    //  Print Geometrical Base Mesh
//    std::ofstream argument("geometry_with_fractures.txt");
//    m_geometry->Print(argument);
//    std::ofstream Dummyfile("geometry_with_fractures.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(m_geometry,Dummyfile, true);
    
}

int64_t TNFRElastic::CreateAndDuplicateConnect(TPZInterpolatedElement *intel, int local_index, TPZConnect &connect){
    
    // Create a duplicated connect
    int64_t index = m_cmesh->AllocateNewConnect(connect.NShape(), connect.NState(), connect.Order());
    intel->SetConnectIndex(local_index, index);
    connect.DecrementElConnected();
    m_cmesh->ConnectVec()[index].IncrementElConnected();
    return index;
}

void TNFRElastic::ExecuteASingleTimeStep(){
    
    /// Getting solution at n
    m_x_n = m_Operator->Solution();
    
    m_r_tolerance = 0.001;
    unsigned int n_iterations = 10;
    bool stop_criterion_Q = false;
    
    for (m_k = 0; m_k < n_iterations; m_k++) {
        
        /// Compute R
        m_Operator->AssembleResidual();
        m_r_norm = Norm(m_Operator->Rhs());
        
        stop_criterion_Q = m_r_norm <= m_r_tolerance;
        if (stop_criterion_Q) {
            std::cout << "A_e:: Newton process reaches convergence tolerance." << std::endl;
            std::cout << "A_e:: Redisue norm        = " << m_r_norm << std::endl;
            std::cout << "A_e:: Correction norm     = " << m_delta_x_norm << std::endl;
            std::cout << "A_e:: Iterations          = " << m_k << std::endl;
            // UpdateGlobalSolution
            //        UpdateGlobalSolution();
            return;
        }
        
        /// Computes J and -R
        m_Operator->Assemble(); // J
        m_Operator->Rhs() *= -1.0; // -R
        
        /// Computes correction and updates it
        m_Operator->Solve();
        m_delta_x_norm = Norm(m_Operator->Solution());
        m_Operator->Solution() += m_x_n;
    }
    
    if (m_k == n_iterations) {
        std::cout << "A_e:: Newton process does not reach convergence tolerance." << std::endl;
    }
    
}

/// Post-Process operator variables
void TNFRElastic::PostProcess(){
    
    const int dim = m_geometry->Dimension();
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile;
    plotfile =  "elliptic.vtk";
//    scalnames.Push("s_xx");
//    scalnames.Push("s_xy");
//    scalnames.Push("s_xz");
//    scalnames.Push("s_yy");
//    scalnames.Push("s_yz");
//    scalnames.Push("s_zz");
    scalnames.Push("i_1");
    scalnames.Push("j_2");
    scalnames.Push("s_1");
    scalnames.Push("s_2");
    scalnames.Push("s_3");
    vecnames.Push("u");
    m_Operator->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    m_Operator->PostProcess(div);
}




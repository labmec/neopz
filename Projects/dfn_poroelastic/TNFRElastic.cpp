//
//  TNFRElastic.cpp
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#include "TNFRElastic.h"
#include "TNRFElasticMemory.h"
#include "pzbndcond.h"
#include "TPZMatWithMem.h"

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
    std::ofstream out("CmeshDeformation.txt");
    m_cmesh->Print(out);
#endif
    
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
    scalnames.Push("s_xx");
    scalnames.Push("s_xy");
    scalnames.Push("s_xz");
    scalnames.Push("s_yy");
    scalnames.Push("s_yz");
    scalnames.Push("s_zz");
//    scalnames.Push("i_1");
//    scalnames.Push("j_2");
//    scalnames.Push("s_1");
//    scalnames.Push("s_2");
//    scalnames.Push("s_3");
    vecnames.Push("u");
    m_Operator->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    m_Operator->PostProcess(div);
}




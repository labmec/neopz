//
//  TNFRElastic.h
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#ifndef TNFRElastic_h
#define TNFRElastic_h

#include <stdio.h>
#include <vector>

// neopz objects
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"


// materials objects
#include "TNFRElasticMaterial.h"

// boundary description
#include "TNFRBoundaryDescription.h"

class TNFRElastic {

public:
    
    /// Constructor
    TNFRElastic();
    
    /// Descructor
    ~TNFRElastic();
    
    /// Geometry description
    TPZGeoMesh * m_geometry;
    
    /// Computational mesh
    TPZCompMesh * m_cmesh;
    
    /// Construct the operator and the approximation
    TPZAnalysis * m_Operator;
    
    /// Material identifier for fracture boundaries
    unsigned int m_fracture_bc_id = 9;
    
    /// vector of objects describing the boundary data
    std::vector<TNFRBoundaryDescription> m_boundary_data;
    
    /// Number of threads for parallele execution
    unsigned int m_n_threads;
    
    /// Directive for active neopz matrix bandwidth reduction
    bool m_bandwidth_reduction_Q;
    
    /// Set vector of objects describing the boundary data
    void SetBoundaryData(std::vector<TNFRBoundaryDescription> boundary_data);
    
    /// Set the number of threads for parallele execution
    void SetNumberOfThreads(unsigned int n_threads);
    
    /// Set the directive for active neopz matrix bandwidth reduction
    void SetEnableBandwidthReduction(bool bandwidth_reduction_Q);
    
    /// Create and duplicate a connect associated to the provided interpolated element
    int64_t CreateAndDuplicateConnect(TPZInterpolatedElement *intel, int local_index, TPZConnect &connect);
    
    /// Insert one boundary elements to represent fractures
    void InsertFractureRepresentation();
    
    /// build the operator with polynomial order p_order
    bool BuildOperator(TPZGeoMesh *  geometry, int p_order);
    
    /// Execute a single time step
    void ExecuteASingleTimeStep();
    
    /// Post-Process operator variables
    void PostProcess();
    
private:

    /// last state unknowns vector
    TPZFMatrix<REAL> m_x_n;
    
    /// Residue tolerance
    REAL m_r_tolerance;
    
    /// Residue norm
    REAL m_r_norm;
    
    /// Delta x norm
    REAL m_delta_x_norm;
    
    /// Newton iteration counter
    unsigned int m_k;
    
};

#endif /* TNFRElastic_h */

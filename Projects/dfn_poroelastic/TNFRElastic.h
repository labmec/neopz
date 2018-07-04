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
    
    /// vector of objects describing the boundary data
    std::vector<TNFRBoundaryDescription> m_boundary_data;
    
    // build the operator with polynomial order p_order
    void BuildOperator(TPZGeoMesh *  geometry, int p_order);
    
private:
    
    
    
};

#endif /* TNFRElastic_h */

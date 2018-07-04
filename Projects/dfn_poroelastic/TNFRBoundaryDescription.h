//
//  TNFRBoundaryDescription.h
//  pz
//
//  Created by Omar Durán on 7/3/18.
//

#ifndef TNFRBoundaryDescription_h
#define TNFRBoundaryDescription_h

#include <stdio.h>
#include <vector>
#include "pzreal.h"

class TNFRBoundaryDescription {
    
public:
    
    /// Constructor
    TNFRBoundaryDescription();
    
    /// Destructor
    ~TNFRBoundaryDescription();
    
    /// Copy constructor
    TNFRBoundaryDescription(const TNFRBoundaryDescription & other);
    
    /// assignment operator
    const TNFRBoundaryDescription & operator=(const TNFRBoundaryDescription & other);
    
    /// Set boundary description
    void SetBC(int bc_id, int bc_type, std::vector<REAL> bc_values);
    
    /// Set boundary indentifier
    void SetBCId(int bc_id);
    
    /// Set boundary type
    void SetBCType(int bc_type);
    
    /// Set boundary values
    void SetBCValues(std::vector<REAL> bc_values);
    
    /// Get boundary indentifier
    int GetBCId();
    
    /// Get boundary type
    int GetBCType();
    
    /// Get boundary values
    std::vector<REAL> GetBCValues();
    
private:
    
    /// Boundary indentifier
    int m_bc_id;
    
    /// Boundary type
    int m_bc_type;
    
    /// Boundary values
    std::vector<REAL> m_bc_values;
    
};

#endif /* TNFRBoundaryDescription_h */

//
//  TNFRBoundaryDescription.cpp
//  pz
//
//  Created by Omar Durán on 7/3/18.
//

#include "TNFRBoundaryDescription.h"


TNFRBoundaryDescription::TNFRBoundaryDescription(){
    
    m_bc_id     = 0;
    m_bc_type   = 0;
    m_bc_values.resize(0);
    
}

TNFRBoundaryDescription::~TNFRBoundaryDescription(){
    
}

/// Copy constructor
TNFRBoundaryDescription::TNFRBoundaryDescription(const TNFRBoundaryDescription & other){
    m_bc_id     = other.m_bc_id;
    m_bc_type   = other.m_bc_type;
    m_bc_values = other.m_bc_values;
}

/// assignment operator
const TNFRBoundaryDescription & TNFRBoundaryDescription::operator=(const TNFRBoundaryDescription & other){
    
    if (this != & other) // prevent self-assignment
    {
        m_bc_id     = other.m_bc_id;
        m_bc_type   = other.m_bc_type;
        m_bc_values = other.m_bc_values;
    }
    return *this;
}

/// Set boundary description
void TNFRBoundaryDescription::SetBC(int bc_id, int bc_type, std::vector<REAL> bc_values){
    SetBCId(bc_id);
    SetBCType(bc_type);
    SetBCValues(bc_values);
}

/// Set boundary indentifier
void TNFRBoundaryDescription::SetBCId(int bc_id){
    m_bc_id     = bc_id;
}

/// Set boundary type
void TNFRBoundaryDescription::SetBCType(int bc_type){
    m_bc_type   = bc_type;
}

/// Set boundary values
void TNFRBoundaryDescription::SetBCValues(std::vector<REAL> bc_values){
    m_bc_values = bc_values;
}

/// Get boundary indentifier
int TNFRBoundaryDescription::GetBCId(){
    return m_bc_id;
}

/// Get boundary type
int TNFRBoundaryDescription::GetBCType(){
    return m_bc_type;
}

/// Get boundary values
std::vector<REAL> TNFRBoundaryDescription::GetBCValues(){
    return m_bc_values;
}



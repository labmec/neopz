//
//  TPZPMRSMemory.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#include "TPZPMRSMemory.h"

/** @brief Default constructor */
TPZPMRSMemory::TPZPMRSMemory(){
  
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // initial state items
    
    /** @brief gradient of u_n at intial state*/
    m_grad_u_0.Resize(3, 3);
    m_grad_u_0.Zero();
    
    
    // last time state items
    
    /** @brief displacements */
    m_u.Resize(3, 0);
    
    /** @brief gradient of u_n */
    m_grad_u.Resize(3, 3);
    m_grad_u.Zero();
    
    
    // current time state items
    
    /** @brief displacements */
    m_u_n.Resize(3, 0);
    
    /** @brief gradient of u_n */
    m_grad_u_n.Resize(3, 3);
    m_grad_u_n.Zero();
    
    
    /** @brief elastic strain at n */
    m_epsilon_e_n.Resize(3, 3);
    m_epsilon_e_n.Zero();
    
    
    /** @brief plastic strain at n */
    m_epsilon_p_n.Resize(3, 3);
    m_epsilon_p_n.Zero();
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // initial state items
    
    /** @brief weighted pressure at intial state */
    m_p_0 = 0.0;
    
    /** @brief weighted pressure at the previous timestep */
    m_p_n = 0.0;
    
    
}

/** @brief Default destructor */
TPZPMRSMemory::~TPZPMRSMemory(){
    
}

int TPZPMRSMemory::ClassId() const{
    return Hash("TPZPMRSMemory");
}


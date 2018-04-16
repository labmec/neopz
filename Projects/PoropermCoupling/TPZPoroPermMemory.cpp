//
//  TPZPoroPermMemory.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#include "TPZPoroPermMemory.h"

/** @brief Default constructor */
TPZPoroPermMemory::TPZPoroPermMemory(){
    
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
}

/** @brief Default destructor */
TPZPoroPermMemory::~TPZPoroPermMemory(){
    
}

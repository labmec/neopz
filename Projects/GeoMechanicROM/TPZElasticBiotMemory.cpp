//
//  TPZElasticBiotMemory.cpp
//  PZ
//
//  Created by omar on 3/6/17.
//
//

#include "TPZElasticBiotMemory.h"

/** @brief Default constructor */
TPZElasticBiotMemory::TPZElasticBiotMemory(){
    
    /** @brief RB functions */
    fphi_u.Resize(0,0);
    fgrad_phi_u.Resize(0,0);
    
    /** @brief displacements */
    fu_n.Resize(3,1);
    fu_n.Zero();
    
    /** @brief gradient of u_n */
    fgrad_u_n.Resize(3, 3);
    fgrad_u_n.Zero();
    
    // Cross terms
    /** @brief displacements */
    fu.Resize(3,1);
    fu.Zero();
    
    /** @brief gradient of u */
    fgrad_u.Resize(3, 3);
    fgrad_u.Zero();
    
    /** @brief displacements_n */
    fu_n.Resize(3,1);
    fu_n.Zero();
    
    /** @brief gradient of u_n */
    fgrad_u_n.Resize(3, 3);
    fgrad_u_n.Zero();
    
    /** @brief gradient of u_n */
    fporosity_0 = 0.0;
    
    /** @brief pressure at initial state */
    fp_0 = 0.0;
    
    /** @brief gradient of u_0 at intial state*/
    fgrad_u_0.Resize(3, 3);
    fgrad_u_0.Zero();
    
    /** @brief pressure at current state */
    fp_n = 0.0;
    
    /** @brief pressure at current state */
    fp = 0.0;
    
    /** @brief gradient of u_n at intial state*/
    fgrad_p_n.Resize(0, 0);
    
    /** @brief gradient of u_n at intial state*/
    fgrad_p.Resize(0, 0);
    
}

/** @brief Default destructor */
TPZElasticBiotMemory::~TPZElasticBiotMemory(){
    
}

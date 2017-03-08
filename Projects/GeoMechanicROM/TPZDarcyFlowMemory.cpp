//
//  TPZDarcyFlowMemory.cpp
//  PZ
//
//  Created by omar on 3/7/17.
//
//

#include "TPZDarcyFlowMemory.h"

/** @brief Default constructor */
TPZDarcyFlowMemory::TPZDarcyFlowMemory(){

    /** @brief flux space functions */
    fphi_q.Resize(0, 0);
    
    /** @brief flux space div functions */
    fdiv_phi_q.Resize(0, 0);
    
    /** @brief pressure space functions */
    fphi_p.Resize(0, 0);
    
    /** @brief pressure space functions */
    fgrad_phi_p.Resize(0, 0);
    
    /** @brief flux */
    fq_n.Resize(0,0);
    
    /** @brief displacements */
    fu_n.Resize(0, 0);
    
    /** @brief gradient of u_n */
    fgrad_u_n.Resize(0, 0);
    
    /** @brief displacements */
    fu.Resize(0, 0);
    
    /** @brief gradient of u_n */
    fgrad_u.Resize(0, 0);
    
    /** @brief gradient of u_n */
    fporosity_0 = 0.0;
    
    /** @brief pressure at initial state */
    fp_0 = 0.0;
    
    /** @brief gradient of u_n at intial state*/
    fgrad_u_0.Resize(0, 0);
    
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
TPZDarcyFlowMemory::~TPZDarcyFlowMemory(){
    
}

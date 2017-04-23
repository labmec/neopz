//
//  TRMMemory.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//
//

#include "TRMMemory.h"


/** @brief Default constructor */
TRMMemory::TRMMemory(){
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Elliptic memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // Basis functions
    
    /** @brief elliptic functions functions */
    f_e_phi_u.Resize(3, 0.0);
    
    /** @brief elliptic functions functions */
    f_e_grad_phi_u.Resize(3, 0.0);
    
    
    // initial state items
    
    /** @brief gradient of u_n at intial state*/
    f_e_grad_u_0.Resize(3, 0.0);
    
    /** @brief sigma at intial state*/
    f_e_sigma_0.Resize(3, 0.0);
    
    
    // last time state items
    
    /** @brief displacements */
    f_e_u.Resize(3, 0.0);
    
    /** @brief gradient of u_n */
    f_e_grad_u.Resize(3, 0.0);
    
    // current time state items
    
    /** @brief displacements */
    f_e_u_n.Resize(3, 0.0);
    
    /** @brief gradient of u_n */
    f_e_grad_u_n.Resize(3, 0.0);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Parabolic memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // Basis functions
    
    /** @brief parabolic q base functions */
    f_p_phi_q.Resize(3, 0.0);
    
    /** @brief parabolic div_q base functions */
    f_p_div_phi_q.Resize(3, 0.0);
    
    /** @brief parabolic p base functions */
    f_p_phi_p.Resize(3, 0.0);
    
    
    // initial state items
    
    /** @brief weighted pressure at intial state */
    f_p_p_0 = 0.0;
    
    /** @brief Rock Porosity */
    f_p_phi_0 = 0.0;
    
    /** @brief absolute permeability */
    f_p_K_0.Resize(3, 0.0);
    
    /** @brief absolute permeability inverse */
    f_p_Kinv_0.Resize(3, 0.0);
    
    // last time state items
    
    /** @brief total velocity */
    f_p_q.Resize(0);
    
    /** @brief divergence of velocity at last state */
    f_p_div_q = 0.0;
    
    /** @brief weighted pressure */
    f_p_p = 0.0;
    
    
    // current time state items
    
    /** @brief total velocity */
    f_p_q_n.Resize(0);
    
    /** @brief divergence of velocity at last state */
    f_p_div_q_n = 0.0;
    
    /** @brief weighted pressure at the previous timestep */
    f_p_p_n = 0.0;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    fu.Resize(3, 0.0);
    fu_n.Resize(3, 0.0);
    fdivu               = 0.0;
    fdivu_n             = 0.0;
    
    // Required
    fp              = 0.0;
    fp_n            = 0.0;
    fp_avg          = 0.0;
    fp_avg_n        = 0.0;
    fsa             = 0.0;
    fsa_n           = 0.0;
    fsb             = 0.0;
    fsb_n           = 0.0;

    
    fx.Resize(3,0.0);
    fw = 0.0;
    fdet = 0.0;
    frhs = 0.0;
    
}

/** @brief Default destructor */
TRMMemory::~TRMMemory(){
    
}
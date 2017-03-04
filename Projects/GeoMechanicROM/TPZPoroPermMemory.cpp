//
//  TPZPoroPermMemory.cpp
//  PZ
//
//  Created by Omar on 9/6/16.
//
//

#include "TPZPoroPermMemory.h"

/** @brief Default constructor */
TPZPoroPermMemory::TPZPoroPermMemory(){
    
    /** @brief RB functions */
    fphi_u.Resize(0,0);
    fgrad_phi_u.Resize(0,0);
    
    /** @brief displacements */
    fu_n.Resize(3,1);
    fu_n.Zero();
    
    /** @brief gradient of u_n */
    fgrad_u_n.Resize(3, 3);
    fgrad_u_n.Zero();
    
    /** @brief displacements */
    fu.Resize(3,1);
    fu.Zero();
    
    /** @brief gradient of u_n */
    fgrad_u.Resize(3, 3);
    fgrad_u.Zero();
}

/** @brief Default destructor */
TPZPoroPermMemory::~TPZPoroPermMemory(){
    
}

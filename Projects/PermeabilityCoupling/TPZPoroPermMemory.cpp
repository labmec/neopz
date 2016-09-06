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
    
    /** @brief displacements */
    fu_n.Resize(0, 0);
    
    /** @brief gradient of u_n */
    fgrad_u_n.Resize(0, 0);
    
    /** @brief elastic strain at n */
    fepsilon_e_n.Resize(0, 0);
    
    /** @brief plastic strain at n */
    fepsilon_p_n.Resize(0, 0);
    
}

/** @brief Default destructor */
TPZPoroPermMemory::~TPZPoroPermMemory(){
    
}

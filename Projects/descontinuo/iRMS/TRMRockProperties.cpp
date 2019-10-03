//
//  TRMRockProperties.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMRockProperties.h"


/** @brief Default constructor */
TRMRockProperties::TRMRockProperties(){
    
    /** @brief Characteristic length - m */
    fLref = 1.0;
    
    /** @brief Characteristic Permeability - m2 */
    fKref = 1.0;
    
    /** @brief Characteristic Pressure - Pa */
    fPref = 1.0;
    
    /** @brief Layer Top depth  - m */
    fDepth = 1.0;
    
    /** @brief Porosity at P of reference - */
    fPhiref = 1.0;
    
    /** @brief Rock Compressibility 1/pa - */
    fcrock = 1.0;
    
    /** @brief Is GID geometry - */
    fIsGIDGeometry = false;
    
    /** @brief absolute permeability */
    fKab.Resize(3,3);
    fKab.Zero();
    
    /** @brief absolute permeability inverse */
    fKabinv.Resize(3,3);
    fKabinv.Zero();
    
}

/** @brief Default destructor */
TRMRockProperties::~TRMRockProperties(){
    
}
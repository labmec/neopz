//
//  TRMSpatialPropertiesMap.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpatialPropertiesMap.h"

/** @brief default constructor */
TRMSpatialPropertiesMap::TRMSpatialPropertiesMap(){
    
}

/** @brief default destructor */
TRMSpatialPropertiesMap::~TRMSpatialPropertiesMap(){
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::Kappa(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MapModel()) {
        case 0:
        {
            this->Kappa_c(x, kappa, inv_kappa, state_vars);
        }
            break;
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        default:
        {
            std::cout << "Error: Model not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::phi(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MapModel()) {
        case 0:
        {
            this->phi_c(x, phi, state_vars);
        }
            break;
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        default:
        {
            std::cout << "Error: Model not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::Kappa_c(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars){

    kappa.Resize(3,3);
    kappa.Zero();
    STATE val = 1.0e-13;
    kappa(0,0) = val;
    kappa(1,1) = val;
    kappa(2,2) = val;
    
    kappa.Inverse(inv_kappa,ELU);
}




/** @brief Porosity fraction  $\phi$ */
void TRMSpatialPropertiesMap::phi_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars){
    
    phi.Resize(10, 0.0);
    STATE val = 0.1;
    phi[0] = val;
    
}
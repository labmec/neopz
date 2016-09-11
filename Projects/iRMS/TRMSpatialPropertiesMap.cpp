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

/** @brief Porosity fraction  $\phi$ */
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

/** @brief first lamé parameter $\lambda$ */
void TRMSpatialPropertiesMap::lambda(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MapModel()) {
        case 0:
        {
            this->lambda_c(x, lambda, state_vars);
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


/** @brief undrained first lamé parameter  $\lambda_{u}$ */
void TRMSpatialPropertiesMap::lambda_u(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda_u, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MapModel()) {
        case 0:
        {
            this->lambda_u_c(x, lambda_u, state_vars);
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

/** @brief second lamé parameter  $\mu$ */
void TRMSpatialPropertiesMap::mu(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MapModel()) {
        case 0:
        {
            this->mu_c(x, mu, state_vars);
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

/** @brief Biot's poroelastic parameter  $\alpha$ */
void TRMSpatialPropertiesMap::alpha(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &alpha, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MapModel()) {
        case 0:
        {
            this->alpha_c(x, alpha, state_vars);
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

/** @brief Geological Stress $\sigma_{0}$ */
void TRMSpatialPropertiesMap::S_0(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &s_0){
    
    s_0.Resize(3, 3);
    s_0.Zero();
    REAL MPa = 1.0e6;
    s_0(0,0) = -40.0*MPa;
    s_0(1,1) = -50.0*MPa;
    s_0(2,2) = -60.0*MPa;
    
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::Kappa_c(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars){

    kappa.Resize(3,3);
    kappa.Zero();
    STATE val = 1.0e-13;
    kappa(0,0) = val;
    kappa(1,1) = val;
    kappa(2,2) = val;
    
    inv_kappa.Resize(3,3);
    inv_kappa.Zero();
    inv_kappa(0,0) = 1.0/kappa(0,0);
    inv_kappa(1,1) = 1.0/kappa(1,1);
    inv_kappa(2,2) = 1.0/kappa(2,2);
    
}




/** @brief Porosity fraction  $\phi$ */
void TRMSpatialPropertiesMap::phi_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars){
    
    phi.Resize(10, 0.0);
    STATE val = 0.25;
    phi[0] = val;
    
}

/** @brief first lamé parameter $\lambda$ */
void TRMSpatialPropertiesMap::lambda_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda, TPZManVector<STATE,10> &state_vars){
    
    REAL GPa = 1.0e9;
    lambda.Resize(10, 0.0);
    STATE val = 9.4541*GPa;
    lambda[0] = val;
    
}


/** @brief undrained first lamé parameter  $\lambda_{u}$ */
void TRMSpatialPropertiesMap::lambda_u_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda_u, TPZManVector<STATE,10> &state_vars){
    
    REAL GPa = 1.0e9;
    lambda_u.Resize(10, 0.0);
    STATE val = 9.5541*GPa;
    lambda_u[0] = val;
    
}

/** @brief second lamé parameter  $\mu$ */
void TRMSpatialPropertiesMap::mu_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars){
    
    REAL GPa = 1.0e9;
    mu.Resize(10, 0.0);
    STATE val = 0.0965*GPa;
    mu[0] = val;
    
}

/** @brief Biot's poroelastic parameter  $\alpha$ */
void TRMSpatialPropertiesMap::alpha_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &alpha, TPZManVector<STATE,10> &state_vars){
    
    alpha.Resize(10, 0.0);
    STATE val = 1.0;
    alpha[0] = val;
    
}
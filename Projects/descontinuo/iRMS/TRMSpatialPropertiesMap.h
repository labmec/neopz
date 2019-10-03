//
//  TRMSpatialPropertiesMap.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMSpatialPropertiesMap__
#define __PZ__TRMSpatialPropertiesMap__

#include <stdio.h>
#include "pzmanvector.h"
#include "pzfmatrix.h"

class TRMSpatialPropertiesMap{
    
    
private:
    
    // Here, given a x, it will return the permeability, porosity and dphidp
    
    /** @brief spatial properties model */
    int fMap_model; // map_model = {0 -> constan map, 1 -> linear map, 2 -> kriged map}
    
    // Constant case
    
    void Kappa_c(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    void lambda_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda, TPZManVector<STATE,10> &state_vars);
    
    void lambda_u_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda_u, TPZManVector<STATE,10> &state_vars);
    
    void mu_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars);
    
    void alpha_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &alpha, TPZManVector<STATE,10> &state_vars);
    
    // CGAL interpoaltion
    
    void Kappa_f(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi_f(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    double phi_cdf(double x);
    
public:
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap();
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap(TRMSpatialPropertiesMap &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap &operator=(TRMSpatialPropertiesMap &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief default destructor */
    ~TRMSpatialPropertiesMap();
    
    /**
     * @defgroup Set and get methods
     * @{
     */
    
    /** @brief Set spatial properties model */
    void SetMapModel(int model){
        fMap_model = model;
    }
    
    /** @brief Get spatial properties model */
    int MapModel()
    {
        return fMap_model;
    }

    /** @brief Geological Stress $\sigma_{0}$ */    
    void S_0(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &s_0);
    
    void Kappa(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    void lambda(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda, TPZManVector<STATE,10> &state_vars);
    
    void lambda_u(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &lambda_u, TPZManVector<STATE,10> &state_vars);
    
    void mu(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars);
    
    void alpha(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &alpha, TPZManVector<STATE,10> &state_vars);
    
    

};

#endif /* defined(__PZ__TRMSpatialPropertiesMap__) */

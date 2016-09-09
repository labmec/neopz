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
    
    // @}
    
    
    /**
     * @defgroup Constant map models models
     * @{
     */
    
    void Kappa_c(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    // @}
    
    
    
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
    
    
    /**
     * @defgroup map models models
     * @{
     */
    
    void Kappa(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars);
    
    void phi(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars);
    
    // @}
    

    
    
};

#endif /* defined(__PZ__TRMSpatialPropertiesMap__) */

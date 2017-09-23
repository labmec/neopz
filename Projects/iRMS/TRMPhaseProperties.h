//
//  TRMPhaseProperties.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMPhaseProperties__
#define __PZ__TRMPhaseProperties__

#include <stdio.h>
#include "pzmanvector.h"

class TRMPhaseProperties{

protected:
    
    /** @brief Temperature @ reservoir conditions  - F */
    REAL fTRes;
    
    /** @brief Temperature for references values [K] */
    REAL fTRef;
    
    /** @brief Pressure for references values - Pa */
    REAL fPRef;
    
private:
    
    /** @brief Density - kg/m3  $\rho$ */
    TPZManVector<STATE,10> frho;
    
    /** @brief viscosity - Pa s  $\mu$ */
    TPZManVector<STATE,10> fmu;
    
    /** @brief Compressibility - 1/pa $c$ */
    TPZManVector<STATE,10> fc;
    
    /** @brief Density model = {0,1,2}  */
    int frho_model;
    
    /** @brief Viscosity model = {0,1,2} */
    int fmu_model;
    
    /** @brief Compressibility model = {0,1,2} */
    int fc_model;
    
public:
    
    /** @brief default constructor */
    TRMPhaseProperties();
    
    /** @brief default constructor */
    TRMPhaseProperties(TRMPhaseProperties &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMPhaseProperties &operator=(TRMPhaseProperties &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief default destructor */
    ~TRMPhaseProperties();

    /**
     * @brief Density - kg/m3  $\rho$. Compute the propertie f and partial derivatives
     * @param state_vars = {P,T,S1,S2,S3}.
     * @return f = {f,dfdP,dfdT,dfdS1,dfdS2,dfdS3}
     */
    virtual void Density(TPZManVector<STATE,10> &rho, TPZManVector<STATE,10> &state_vars) = 0;
    
    /**
     * @brief Viscosity - Pa s  $\mu$. Compute the propertie f and partial derivatives
     * @param state_vars = {P,T,S1,S2,S3}.
     * @return f = {f,dfdP,dfdT,dfdS1,dfdS2,dfdS3}
     */
    virtual void Viscosity(TPZManVector<STATE,10> &mu, TPZManVector<STATE,10> &state_vars) = 0;
    
    /**
     * @brief Compressibility - 1/pa $c$ . Compute the propertie f and partial derivatives
     * @param state_vars = {P,T,S1,S2,S3}.
     * @return f = {f,dfdP,dfdT,dfdS1,dfdS2,dfdS3}
     */
    virtual void Compressibility(TPZManVector<STATE,10> &c, TPZManVector<STATE,10> &state_vars) = 0;
    
    /** @brief Set density model = {0,1,2}  */
    void SetRhoModel(int model){
        frho_model = model;
    }
    
    /** @brief Get density model = {0,1,2}  */
    int RhoModel(){
        return frho_model;
    }
    
    /** @brief Set viscosity model = {0,1,2}  */
    void SetMuModel(int model){
        fmu_model = model;
    }
    
    /** @brief Get viscosity model = {0,1,2}  */
    int MuModel(){
        return fmu_model;
    }
    
    /** @brief Set compressibility model = {0,1,2}  */
    void SetCModel(int model){
        fc_model = model;
    }
    
    /** @brief Get compressibility model = {0,1,2}  */
    int CModel(){
        return fc_model;
    }
    
    
};

#endif /* defined(__PZ__TRMPhaseProperties__) */

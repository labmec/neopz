/* 
 * File:   TPZYCCamClayPV.h
 * Author: quinelato
 *
 * Created on August 14, 2017, 5:18 PM
 */

#ifndef TPZYCCAMCLAYPV_H
#define TPZYCCAMCLAYPV_H

#include "pzreal.h"
#include "pzvec.h"
#include "TPZElasticResponse.h"
#include "pzfmatrix.h"
#include "TPZPlasticCriterion.h"

class TPZYCCamClayPV : public TPZPlasticCriterion {
public:

    enum {
        NYield = 1
    };

    TPZYCCamClayPV();
    TPZYCCamClayPV(const TPZYCCamClayPV& other);
    void SetUp(const TPZElasticResponse &ER, REAL gamma, REAL m, REAL pt, REAL logHardening, REAL logBulkModulus, REAL a0, REAL e0);
    void SetElasticResponse(const TPZElasticResponse &ER);
    virtual int ClassId() const override;
    void Read(TPZStream& buf, void* context) override;
    void Write(TPZStream& buf, int withclassid) const override;
    REAL bFromP(const REAL p, const REAL a) const;
    REAL bFromTheta(REAL theta) const;
    void Phi(TPZVec<REAL> sigmaPV, REAL a, TPZVec<REAL> &phi) const;
    
    /**
     * Computes the cylindrical coordinates of the point on the yield surface with the given xi and beta.
     * @param xi Hydrostatic component
     * @param beta Lode angle
     * @param returnValue Cylindrical coordinates of the point on the yield surface
     */
    void SurfaceInCyl(const REAL theta, const REAL beta, const REAL a, TPZVec<REAL> &returnValue) const;
    
    REAL ResLFunc(const TPZVec<STATE> &sigma_trial_pv, STATE theta, STATE beta, REAL a, REAL aPrev) const;
//    REAL DResLFunc(const TPZVec<STATE> &sigma_trial_pv, STATE theta, STATE beta, REAL a, REAL aPrev) const;
    
    REAL DistanceToSurface(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta, const REAL a) const;
    
    /**
     * Compute initial damage variable from the given principal stress state
     * @param stress_p principal values
     * @return the damage variable
     */
    REAL InitialDamage(const TPZVec<REAL> &stress_p) const;
    
    /**
     * Computes the derivative of the distance function to the yield surface as a function of theta, beta and a
     * @param pt
     * @param xi
     * @param beta
     * @param fxn
     */
    void DDistanceToSurface(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, const REAL aPrev, TPZVec<STATE> &fxn) const;
    
    /**
     * Computes the second derivative of the distance as a function of theta, beta and a
     * @param pt
     * @param xi
     * @param beta
     * @param jac
     */
    void D2DistanceToSurface(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, TPZFNMatrix<9, STATE> &jac) const;
    
    void ProjectToSurfaceConstantBeta(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma_pv, REAL &aProj, const REAL tol) const;
    void ProjectToSurface(const TPZVec<REAL> &sigma_trial, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, const REAL tol) const;
    void ProjectSigma(const TPZVec<REAL> &sigma_trial, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, int &m_type, TPZFMatrix<REAL> * gradient = NULL) const;
    
    void SurfaceParam(const TPZVec<STATE> &sigma_pv, const STATE a, STATE &theta, STATE &beta) const;
    
    void GradSigmaTrial(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta, const REAL aProj, TPZFNMatrix<9, STATE> &ddist_dsigmatrial) const;
    
    /* \frac{\partial \sigma}{\partial (theta, beta, a)}*/
    void DFuncCart(STATE theta, STATE beta, STATE a, TPZFNMatrix<9, STATE> &DFunccart) const;
    
    void ProjectSigmaDep(const TPZVec<REAL> &sigma_trial, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, TPZFMatrix<REAL> &GradSigma) const;
    STATE PlasticVolumetricStrain(STATE a) const;
    virtual ~TPZYCCamClayPV();
    
    void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override{
        Phi(sigma, kprev, yield);
    }

    virtual int GetNYield() const override {
        return as_integer(NYield);
    }
    
    friend class TPZYCDruckerPragerPV;
private:
    TPZElasticResponse fER;

    REAL fGamma;
    REAL fM;
    REAL fPt;
    REAL fLogHardening; // Logarithmic hardening constant
    REAL fLogBulkModulus; // Logarithmic bulk modulus
    REAL fA0; // Initial size of the yield surface
    REAL fE0; // initial void ratio
    
};

#endif /* TPZYCCAMCLAYPV_H */


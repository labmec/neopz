/* 
 * File:   TPZYCDruckerPragerPV.h
 * Author: quinelato
 *
 * Created on September 12, 2017, 2:24 PM
 */

#ifndef TPZYCDRUCKERPRAGERPV_H
#define TPZYCDRUCKERPRAGERPV_H

#include "TPZYCCamClayPV.h"

class TPZYCDruckerPragerPV : public TPZPlasticCriterion {
public:

    enum {
        NYield = 2
    };

    TPZYCDruckerPragerPV();
    TPZYCDruckerPragerPV(const TPZYCDruckerPragerPV& other);
    TPZYCDruckerPragerPV& operator=(const TPZYCDruckerPragerPV& other);

    void SetUp(const TPZElasticResponse &ER, REAL gamma, REAL m, REAL pt, REAL logHardening, REAL logBulkModulus, REAL a0, REAL e0);
    void SetElasticResponse(const TPZElasticResponse &ER);
    virtual int ClassId() const override;
    void Read(TPZStream& buf, void* context) override;
    void Write(TPZStream& buf, int withclassid) const override;
    REAL bFromP(const REAL p, const REAL a) const;
    REAL bFromTheta(REAL theta) const;
    void Phi(TPZVec<REAL> sigmaPV, REAL a, TPZVec<REAL> &phi) const;

    /**
     * Compute initial damage variable from the given principal stress state
     * @param stress_p principal values
     * @return the damage variable
     */
    REAL InitialDamage(const TPZVec<REAL> &stress_p) const;
    
    /**
     * Computes the cylindrical coordinates of the point on the yield surface with the given xi and beta.
     * @param xi Hydrostatic component
     * @param beta Lode angle
     * @param returnValue Cylindrical coordinates of the point on the yield surface
     */
    void SurfaceF1InCyl(const REAL xi, const REAL beta, TPZVec<REAL> &returnValue) const;

    /**
     * Computes the cylindrical coordinates of the point on the yield surface with the given xi and beta.
     * @param theta Deviatoric component
     * @param beta Lode angle
     * @param returnValue Cylindrical coordinates of the point on the yield surface
     */
    void SurfaceF2InCyl(const REAL theta, const REAL beta, const REAL a, TPZVec<REAL> &returnValue) const;

    REAL ResLF1(const TPZVec<STATE> &sigma_trial_pv, const TPZVec<STATE> &sigma_proj_pv, const STATE a, const STATE aPrev) const;

    REAL ResLF2(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, const REAL aPrev) const;

    /**
     * Derivative of the residual equation with respect to a.
     * @param sigma_trial_pv
     * @param sigma_proj_pv
     * @param a
     * @param aPrev
     * @return 
     */
    STATE DResLF1(const TPZVec<STATE> &sigma_trial_pv, const TPZVec<STATE> &sigma_proj_pv, const STATE a, const STATE aPrev) const;

    REAL DistanceToSurfaceF1(const TPZVec<REAL> &sigma_trial_pv, const REAL xi, const REAL beta) const;

    REAL DistanceToSurfaceF2(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta, const REAL a) const;

    /**
     * Computes the derivative of the distance function to the yield surface as a function of xi and beta
     * @param sigma_trial_pv
     * @param xi
     * @param beta
     * @param fxn
     */
    void DDistanceToSurfaceF1(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, TPZVec<STATE> &fxn) const;

    /**
     * Computes the derivative of the distance function to the yield surface as a function of xi, beta and a
     * @param sigma_trial_pv
     * @param theta
     * @param beta
     * @param a
     * @param aPrev
     * @param fxn
     */
    void DDistanceToSurfaceF2(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, const REAL aPrev, TPZVec<STATE> &fxn) const;

    /**
     * Computes the second derivative of the distance as a function of xi and beta
     * @param sigma_trial_pv
     * @param xi
     * @param beta
     * @param jac
     */
    void D2DistanceToSurfaceF1(const TPZVec<STATE> &sigma_trial_pv, const STATE xi, const STATE beta, TPZFNMatrix<4, STATE> &jac) const;

    /**
     * Computes the second derivative of the distance as a function of theta, beta and a
     * @param sigma_trial_pv
     * @param theta
     * @param beta
     * @param a
     * @param jac
     */
    void D2DistanceToSurfaceF2(const TPZVec<STATE> &sigma_trial_pv, const STATE theta, const STATE beta, const REAL a, TPZFNMatrix<9, STATE> &jac) const;

    void ProjectToSurfaceF1(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, const REAL tol) const;

    void ProjectToSurfaceF2(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, const REAL tol) const;

    void ProjectSigma(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma_pv, REAL &aProj, int &m_type, TPZFMatrix<REAL> * gradient = NULL) const;

    void SurfaceParamF1(const TPZVec<STATE> &sigma_pv, STATE &xi, STATE &beta, const REAL tol=1e-5) const;

    void SurfaceParamF2(const TPZVec<STATE> &sigma_pv, const STATE a, STATE &theta, STATE &beta, const REAL tol = 1e-5) const;

    void GradF1SigmaTrial(const TPZVec<REAL> &sigma_trial_pv, const REAL xi, const REAL beta, const REAL aProj, TPZFNMatrix<6, STATE> &ddist_dsigmatrial) const;

    void GradF2SigmaTrial(const TPZVec<REAL> &sigma_trial_pv, const REAL theta, const REAL beta, const REAL aProj, TPZFNMatrix<9, STATE> &ddist_dsigmatrial) const;

    /* \frac{\partial \sigma}{\partial (xi, beta)}*/
    void DF1Cart(STATE xi, STATE beta, TPZFNMatrix<6, STATE> &DFunccart) const;

    /* \frac{\partial \sigma}{\partial (theta, beta, a)}*/
    void DF2Cart(STATE theta, STATE beta, STATE a, TPZFNMatrix<9, STATE> &DFunccart) const;

    void ProjectSigmaDep(const TPZVec<REAL> &sigma_trial_pv, const REAL aPrev, TPZVec<REAL> &sigma, REAL &aProj, TPZFMatrix<REAL> &GradSigma) const;
    
    STATE PlasticVolumetricStrain(STATE a) const;
    
    virtual ~TPZYCDruckerPragerPV();
    
    void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override{
        Phi(sigma, kprev, yield);
    }
    
    virtual int GetNYield() const override{
        return as_integer(NYield);
    }

private:

    TPZYCCamClayPV fCap;
    TPZElasticResponse &fER;

    REAL &fM;
    REAL &fPt;
    REAL &fLogHardening; // Logarithmic hardening constant
    REAL &fLogBulkModulus; // Logarithmic bulk modulus
    REAL &fA0; // Initial size of the yield surface
    REAL &fE0; // initial void ratio

};

#endif /* TPZYCDRUCKERPRAGERPV_H */


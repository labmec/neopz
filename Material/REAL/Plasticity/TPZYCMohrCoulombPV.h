/*
 *  TPZMohrCoulomb.h
 *  FEMPZ
 *
 *  Created by Nathan Shauer on 5/4/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZYCMOHRCOULOMBPV_H
#define TPZYCMOHRCOULOMBPV_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticCriterion.h"

#ifdef LOG4CXX
static LoggerPtr loggerMohrCoulombPV(Logger::getLogger("pz.plasticity.mohrcoulombpv"));
#endif

class TPZYCMohrCoulombPV : public TPZPlasticCriterion {  
public:

    enum {
        NYield = 3
    };

private:
    
    REAL fPhi;
    REAL fPsi;
    REAL fc;
    TPZElasticResponse fER;

protected:
    REAL fEpsPlasticBar;

public:

    /// structure which contains the decision tree of the return map
    // we can only expect a consistent tangent matrix if the decision tree remains the same

    struct TComputeSequence {

        TComputeSequence() : fWhichPlane(ENoPlane), fGamma(0) {

        }

        TComputeSequence(const TComputeSequence &copy) : fWhichPlane(copy.fWhichPlane), fGamma(copy.fGamma) {

        }

        TComputeSequence &operator=(const TComputeSequence &copy) {
            fWhichPlane = copy.fWhichPlane;
            fGamma = copy.fGamma;
            return *this;
        }

        enum MPlane {
            ENoPlane, EElastic, EMainPlane, ERightEdge, ELeftEdge, EApex
        };

        MPlane fWhichPlane;

        TPZManVector<REAL> fGamma;
    };

public:

    /**
     * @brief empty constructor
     */
    TPZYCMohrCoulombPV();

    /**
     * @brief Constructor seting yc parameters
     */
    TPZYCMohrCoulombPV(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER);

    /**
     * @brief Copy Constructor
     */
    TPZYCMohrCoulombPV(const TPZYCMohrCoulombPV &cp);

    /**
     * @brief Sets up the data
     */
    
    /**
     Setup attributes

     @param Phi Friction angle
     @param Psi Dilation angle
     @param c Cohesion yield stress
     @param ER <#ER description#>
     */
    void SetUp(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER) {
        fPhi = Phi;
        fPsi = Psi;
        fc = c;
        fER = ER;
    }

    /**
     * @brief Operator =
     */
    TPZYCMohrCoulombPV & operator=(const TPZYCMohrCoulombPV &cp);

    virtual int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;

    /**
     * @brief Sets epsbar
     */
    void SetEpsBar(REAL &epsbar) {
        fEpsPlasticBar = epsbar;
    }

    /**
     * @brief Print Method
     */
    virtual void Print(std::ostream &out) const override;

    void SetElasticResponse(const TPZElasticResponse &ER) {
        fER = ER;
    }

    virtual TPZElasticResponse GetElasticResponse() const {
        return fER;
    }

    /**
     * @brief Compute initial damage variable from the given principal stress state
     */
    REAL InitialDamage(const TPZVec<REAL> &stress_p) const;
    
    /**
     * @brief Calculates the value c(epsp) and its derivative
     */
    template <class T>
    void PlasticityFunction(const T epsp, T &c, T &H) const;

    /**
     * @brief sigma = lambda Tr(E)I + 2 mu E
     */
    template<class T>
    TPZVec<T> SigmaElastPV(const TPZVec<T> &deform) const;

    /**
     * @brief Calcula o valor da funcao criteiro de plastificacao
     */
    template<class T>
    T PhiPlane(const TPZVec<T> &sigma) const;

    /**
     * @brief Implements the return map in the plane of the surface
     */
    template<class T>
    bool ReturnMapPlane(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
            TComputeSequence &memory, REAL &epsbarnew) const;

    /**
     * @brief Computes dsigmapr/dsigmatr for the ReturnMapPlane
     */
    void ComputePlaneTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const;

    /**
     * @brief Implements the return map in the left edge of the surface
     */
    template<class T>
    bool ReturnMapLeftEdge(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
            TComputeSequence &memory, REAL &epsbarnew) const;

    /**
     * @brief Computes dsigmapr/dsigmatr for the ReturnMapLeftEdge
     */
    void ComputeLeftEdgeTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const;

    /**
     * @brief Implements the return map in the right edge of the surface
     */
    template<class T>
    bool ReturnMapRightEdge(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
            TComputeSequence &memory, REAL &epsbarnew) const;

    /**
     * @brief Computes dsigmapr/dsigmatr for the ReturnMapRightEdge
     */
    void ComputeRightEdgeTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const;

    /**
     * @brief Implements the return map in the apex
     */
    template<class T>
    bool ReturnMapApex(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
            TComputeSequence &memory, REAL &epsbarnew) const;
    
    /**
     Computes gradient of projected stress at sigma_trial for the ReturnMapApex

     @param gradient gradient of projected stress at sigma_trial
     @param eps_bar_p accumulated plastic strain
     */
    void ComputeApexGradient(TPZMatrix<REAL> & gradient, REAL & eps_bar_p) const;
    
    /**
     Execute the integration algorithm for the Mohr-Coulomb model.
     Source: Computational Methods for Plasticity: Theory and Applications. Eduardo A. de Souza Neto, Djordje Peric, David R. J. Owen (2008).

     @param sigma_trial principal values of trial elastic stress
     @param k_prev previous state of the damage variable
     @param sigma projected stress
     @param k_proj current state of the damage variable after projection
     @param m_type variable that indentify the material deformation behavior
     @param gradient gradient of projected stress at sigma_trial
     */
    void ProjectSigma(const TPZVec<STATE> & sigma_trial, STATE k_prev, TPZVec<STATE> & sigma, STATE &k_proj, int & m_type, TPZFMatrix<REAL> * gradient = NULL);

    // Deprecated
    void ProjectSigmaDep(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kpro, TPZFMatrix<STATE> &tang){
        std::cerr << "Deprecated gradient calculation is incorporated on ProjectSigma method." << std::endl;
        DebugStop();
    }

    /**
     Evaluates the yield criterion

     @param sig_vec principal stress
     @param alpha internal damage variable
     @param phi yield criterion function
     */
    void Phi(TPZVec<STATE> sig_vec, STATE alpha, TPZVec<STATE> &phi)const;

    /**
     Access to Friction angle
     
     @return Friction angle
     */
    STATE Phi() {
        return fPhi;
    }

    /**
     Access to Dilation angle

     @return Dilation angle
     */
    STATE Psi() {
        return fPsi;
    }

    /**
     Access to Cohesion yield stress

     @return Cohesion yield stress
     */
    STATE Cohesion() {
        return fc;
    }

    /**
     Access to Young's modulus

     @return Young's modulus
     */
    STATE E() {
        return fER.E();
    }

    /**
     Access to Poisson's ratio

     @return Poisson's ratio
     */
    STATE Poisson() {
        return fER.Poisson();
    }


    virtual void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override{
        Phi(sigma, kprev, yield);
    }
    
    virtual int GetNYield() const override{
        return as_integer(NYield);
    }
};


#endif //TPZYCMOHRCOULOMBPV_H

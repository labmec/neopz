//
//  TPZElasticResponse.h
//  pz
//
//  Created by Erick Slis Raggio Santos on 04/07/2009.
//


#ifndef TPZElasticResponse_H
#define TPZElasticResponse_H

#include "TPZTensor.h"
#include "pzreal.h"


class TPZElasticResponse : public TPZSavable {
    
    /// First Lamé parameter
    REAL m_lambda;
    
    /// Second Lamé parameter
    REAL m_mu;
    
    /// Reference strain at zero stress state.
    TPZTensor<REAL> m_epsilon_star;
    
    /// Reference stress at zero strain state.
    TPZTensor<REAL> m_sigma_star;
    
public:
    
    /**
     A unique class identifier
     
     @return The class identifier as integer
     */
    int ClassId() const override;
    
    /**
     Default constructor
     */
    TPZElasticResponse();
    
    /**
     Copy constructor
     */
    TPZElasticResponse(const TPZElasticResponse & other);
    
    /**
     Assignment constructor
     */
    TPZElasticResponse & operator=(const TPZElasticResponse & other);
    
    /**
     Write (persistency)
     
     @param buf The TPZStream object
     @param withclassid The class identifier
     */
    void Write(TPZStream &buf, int withclassid) const override;
    
    /**
     Read (persistency)
     
     @param buf The TPZStream object
     @param context pointer to the associated object
     */
    void Read(TPZStream &buf, void *context) override;
    
    
    /**
     Class name
     
     @return constant char with the class name
     */
    const char * Name() const;
    
    /**
     Print
     
     @param out ostream object to write the output
     */
    void Print(std::ostream & out) const;
    
    /**
     Computes the stress tensor
     
     @param epsilon The strain tensor
     @param sigma The stress tensor
     */
    template<class T>
    void ComputeStress(const TPZTensor<T> & epsilon, TPZTensor<T> & sigma) const {

        TPZTensor<T> delta_epsilon(epsilon);
        
        delta_epsilon.XX() -= m_epsilon_star.XX();
        delta_epsilon.YY() -= m_epsilon_star.YY();
        delta_epsilon.ZZ() -= m_epsilon_star.ZZ();
        delta_epsilon.XY() -= m_epsilon_star.XY();
        delta_epsilon.XZ() -= m_epsilon_star.XZ();
        delta_epsilon.YZ() -= m_epsilon_star.YZ();
        
        
        T trace = delta_epsilon.I1();
        sigma.Identity();
        sigma.Multiply(trace, m_lambda);
        sigma.Add(delta_epsilon, 2. * m_mu);
        
        sigma.XX() += m_sigma_star.XX();
        sigma.YY() += m_sigma_star.YY();
        sigma.ZZ() += m_sigma_star.ZZ();
        sigma.XY() += m_sigma_star.XY();
        sigma.XZ() += m_sigma_star.XZ();
        sigma.YZ() += m_sigma_star.YZ();

    }
    
    /**
     Computes the strain tensor
     
     @param sigma The stress tensor
     @param epsilon The strain tensor
     */
    template<class T>
    void ComputeStrain(const TPZTensor<T> & sigma, TPZTensor<T> & epsilon) const {
        const T fac = T((1 / 3.)*(1. / (3. * m_lambda + 2. * m_mu) - 1. / (2. * m_mu)));
        TPZTensor<T> delta_sigma(sigma);
        delta_sigma -= m_sigma_star;
        
        delta_sigma.XX() -= m_sigma_star.XX();
        delta_sigma.YY() -= m_sigma_star.YY();
        delta_sigma.ZZ() -= m_sigma_star.ZZ();
        delta_sigma.XY() -= m_sigma_star.XY();
        delta_sigma.XZ() -= m_sigma_star.XZ();
        delta_sigma.YZ() -= m_sigma_star.YZ();
        
        REAL trace = delta_sigma.I1();
        epsilon.Identity();
        epsilon.Multiply(trace, fac);
        epsilon.Add(delta_sigma, 1. / (2. * m_mu));
        
        epsilon.XX() += m_epsilon_star.XX();
        epsilon.YY() += m_epsilon_star.YY();
        epsilon.ZZ() += m_epsilon_star.ZZ();
        epsilon.XY() += m_epsilon_star.XY();
        epsilon.XZ() += m_epsilon_star.XZ();
        epsilon.YZ() += m_epsilon_star.YZ();
        
    }
    
    
    /**
     Incremental constitutive relation in Voigt notation
     
     @param De Return the De operator
     */
    void De(TPZFMatrix<REAL> & De);
    
    
    /**
     Set elastic parameters using engineering data, i.e. Young modulus and Poisson ratio
     
     @param Eyoung Young modulus
     @param Poisson Poisson ratio
     */
    void SetEngineeringData(REAL Eyoung, REAL Poisson);
    
    
    /**
     Set elastic parameters using Lamé data
     
     @param lambda First Lamé parameter
     @param mu Second Lamé parameter (Shear modulus)
     */
    void SetLameData(REAL lambda, REAL mu);
    
    /**
     Access to the first Lamé parameter
     
     @return The first Lamé parameter
     */
    REAL Lambda() const;
    
    /**
     Access to the bulk modulus
     
     @return The bulk modulus
     */
    REAL K() const;
    
    /**
     Access to the second Lamé parameter
     
     @return The Second Lamé parameter
     */
    REAL Mu() const;
    
    /**
     Access to the shear modulus
     
     @return The shear modulus
     */
    REAL G() const;
    
    /**
     Access to the Young modulus
     
     @return The Young modulus
     */
    REAL E() const;
    
    /**
     Access to the Poisson ratio
     
     @return The Poisson ratio
     */
    REAL Poisson() const;
    
    /// Set the reference strain
    void SetReferenceStrainData(TPZTensor<REAL> & eps_star);
    
    /// Get the reference strain
    TPZTensor<REAL> & ReferenceStrainData();
    
    /// Set the reference stress
    void SetReferenceStressData(TPZTensor<REAL> & sigma_star);
    
    /// Get the reference strain
    TPZTensor<REAL> & ReferenceStressData();
};

#endif /* TPZElasticResponse_h */

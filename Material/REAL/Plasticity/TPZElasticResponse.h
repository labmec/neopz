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
    
    /// Residual strain at zero stress state.
    TPZTensor<REAL> m_eps_star;
    
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

        TPZTensor<T> delta_epsilon;
        delta_epsilon.XX() = -m_eps_star.XX();
        delta_epsilon.XY() = -m_eps_star.XY();
        delta_epsilon.XZ() = -m_eps_star.XZ();
        delta_epsilon.YY() = -m_eps_star.YY();
        delta_epsilon.YZ() = -m_eps_star.YZ();
        delta_epsilon.ZZ() = -m_eps_star.ZZ();
        delta_epsilon += epsilon; // Substract residual strain
        T trace = delta_epsilon.I1();
        sigma.Identity();
        sigma.Multiply(trace, m_lambda);
        sigma.Add(delta_epsilon, 2. * m_mu);
    }
    
    /**
     Computes the strain tensor
     
     @param sigma The stress tensor
     @param epsilon The strain tensor
     */
    template<class T>
    void ComputeStrain(const TPZTensor<T> & sigma, TPZTensor<T> & epsilon) const {
        const T fac = T((1 / 3.)*(1. / (3. * m_lambda + 2. * m_mu) - 1. / (2. * m_mu)));
        REAL trace = sigma.I1();
        epsilon.Identity();
        epsilon.Multiply(trace, fac);
        epsilon.Add(sigma, 1. / (2. * m_mu));
        epsilon += m_eps_star;// Adding residual strain
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
    
    /**
     Set the residual strain
     
     @param lambda First Lamé parameter
     @param mu Second Lamé parameter (Shear modulus)
     */
    void SetResidualStrainData(TPZTensor<REAL> & eps_res);
    
    /**
     Get the residual strain
     
     @param lambda First Lamé parameter
     @param mu Second Lamé parameter (Shear modulus)
     */
    TPZTensor<REAL> & ResidualStrainData();
};

#endif /* TPZElasticResponse_h */

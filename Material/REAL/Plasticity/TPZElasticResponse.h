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
    
public:
    
    
    /**
     A unique class identifier
     
     @return The class identifier as integer
     */
    virtual int ClassId() const;
    
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
    void Write(TPZStream& buf, int withclassid) const;
    
    /**
     Read (persistency)
     
     @param buf The TPZStream object
     @param context pointer to the associated object
     */
    void Read(TPZStream& buf, void* context);
    
    
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
        T trace = epsilon.I1();
        sigma.Identity();
        sigma.Multiply(trace, m_lambda);
        sigma.Add(epsilon, 2. * m_mu);
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
};

#endif /* TPZElasticResponse_h */

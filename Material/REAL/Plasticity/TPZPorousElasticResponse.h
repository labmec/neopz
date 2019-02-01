//
//  TPZPorousElasticResponse.h
//  pz
//
//  Created by Omar Durán on 1/16/19.
//

#ifndef TPZPorousElasticResponse_h
#define TPZPorousElasticResponse_h

#include <stdio.h>
#include "TPZTensor.h"
#include "pzreal.h"
#include "TPZElasticResponse.h"

class TPZPorousElasticResponse : public TPZSavable {
    
private:
    
    /// Logarithmic bulk modulus
    STATE m_kappa;
    
    /// Elastic tensile strengh
    STATE m_pt_el;
    
    /// Initial void ratio
    STATE m_e_0;
    
    /// Initial equivalent pressure stress
    STATE m_p_0;
    
    /// Poisson ratio
    STATE m_nu;
    
    /// Second lamé parameter
    STATE m_mu;
    
    /// Directive for define constant shear modulus calculations (false means constant Poisson ratio)
    bool m_is_G_constant_Q;
    
    /// Directive for define Plain stress state or plane strain state
    bool m_plane_stress_Q;
    
public:

    //@TODO:: Document the Class

    virtual int ClassId() const;
    
    TPZPorousElasticResponse();
    
    TPZPorousElasticResponse(const TPZPorousElasticResponse & other);
    
    TPZPorousElasticResponse & operator=(const TPZPorousElasticResponse & other);
    
    void Write(TPZStream& buf, int withclassid) const;
    
    void Read(TPZStream& buf, void* context);
    
    void SetPorousElasticity(STATE kappa, STATE pt_el, STATE e_0, STATE p_0);
    
    void SetShearModulusConstant(STATE G);
    
    void SetPoissonRatioConstant(STATE nu);
    
    void SetPlaneStrain();
    
    void SetPlaneStress();
    
    const char * Name() const ;
    
    void Print(std::ostream & out) const ;
    
    void G(TPZTensor<STATE> &epsilon, STATE & G, STATE & dG_desp_vol);
    
    void K(TPZTensor<STATE> &epsilon, STATE & G, STATE & dGdesp_vol);
    
    void De_porous(TPZTensor<STATE> & epsilon, TPZFMatrix<STATE> & De);
    
    void Compute(TPZTensor<STATE> & epsilon, TPZTensor<STATE> & sigma);
    
    void ComputeDeformation(TPZTensor<STATE> & sigma, TPZTensor<STATE> & epsilon, TPZTensor<STATE> & sigma_n, TPZTensor<STATE> & epsilon_n);
    
    TPZElasticResponse LinearizedElasticResponse(TPZTensor<STATE> & epsilon);
    
};

#endif /* TPZPorousElasticResponse_h */

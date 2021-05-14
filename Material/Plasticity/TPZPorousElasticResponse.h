//
//  TPZPorousElasticResponse.h
//  pz
//
//  Created by Omar Durán on 1/16/19.
//

#ifndef TPZPorousElasticResponse_h
#define TPZPorousElasticResponse_h

#include <stdio.h>
#include <iostream>
#include "TPZTensor.h"
#include "pzreal.h"
#include "TPZElasticResponse.h"

class TPZPorousElasticResponse : public TPZElasticResponse {
    
private:
    
    /// Logarithmic bulk modulus
    STATE m_kappa;
    
    /// Elastic tensile strengh
    STATE m_pt_el;
    
    /// Initial void ratio
    STATE m_e_0;
    
    /// Initial equivalent pressure stress
    STATE m_p_0;
    
    /// Initial volumetric strain
    STATE m_eps_v_0;
    
    /// Poisson ratio
    STATE m_nu;
    
    /// Second lamé parameter
    STATE m_mu;
    
    /// Directive for define constant shear modulus calculations (false means constant Poisson ratio)
    bool m_is_G_constant_Q;
    
    /// Directive for define Plain stress state or plane strain state
    bool m_plane_stress_Q;
    
public:
    
    /// Unique class identifier
    int ClassId() const override;
    
    /// Default constructor
    TPZPorousElasticResponse();
    
    /// Copy constructor
    TPZPorousElasticResponse(const TPZPorousElasticResponse & other);
    
    /// Assingment constructor
    TPZPorousElasticResponse & operator=(const TPZPorousElasticResponse & other);
    
    /// Write class members
    void Write(TPZStream &buf, int withclassid) const override;
    
    /// Read class members
    void Read(TPZStream &buf, void *context) override;
    
    /// Set Porous elasticity data
    void SetPorousElasticity(STATE kappa, STATE pt_el, STATE e_0, STATE p_0);
    
    /// Set initial volumetric stress
    void Setp_0(STATE p_0);
    
    /// Set void ratio at initial volumetric stress
    void Sete_0(STATE e_0);

    /// Set initial volumetric strain
    void Seteps_v_0(STATE eps_v_0);
    
    /// Set directive to keep constant the Shear modulus
    void SetShearModulusConstant(STATE G);
    
    /// Set directive to keep constant Poisson ratio
    void SetPoissonRatioConstant(STATE nu);
    
    /// Set directive for plain strain state
    void SetPlaneStrain();
    
    /// set directive for plain stress state
    void SetPlaneStress();
    
    const char * Name() const ;
    
    void Print(std::ostream & out) const ;
    
    void p(const TPZTensor<STATE> &epsilon, STATE & p, STATE & dp_desp_vol) const;
    
    void G(const TPZTensor<STATE> &delta_eps, STATE & G, STATE & dG_desp_vol) const;
    
    void K(const TPZTensor<STATE> &delta_eps, STATE & K, STATE & dK_desp_vol) const;
    
    void De(const TPZTensor<STATE> & delta_eps, TPZFMatrix<STATE> & De) const;
    
    void De_G_constant(const TPZTensor<STATE> & delta_eps, TPZFMatrix<STATE> & De) const;
    
    void De_Poisson_constant(const TPZTensor<STATE> & delta_eps, TPZFMatrix<STATE> & De) const;
    
    /// Computes a linear elastic response from function evaluation of non linear expressions
    TPZElasticResponse EvaluateElasticResponse(const TPZTensor<STATE> & epsilon) const;
    
    /// Computes the linear elastic response corresponding to the initial state p_0
    TPZElasticResponse LEInitialState();
    
    template<class T>
    void ComputeStress(const TPZTensor<T> & epsilon, TPZTensor<T> & sigma) const {
        
        if (m_is_G_constant_Q) {
            T trace = T(epsilon.I1());
            TPZTensor<T> Se(epsilon);
            Se.XX() -= trace/3.0;
            Se.YY() -= trace/3.0;
            Se.ZZ() -= trace/3.0;
            sigma = Se;
            sigma *= m_mu;
        
            STATE p, dp_desp_vol;
            this->p(epsilon, p, dp_desp_vol);
            sigma.XX() -= p;
            sigma.YY() -= p;
            sigma.ZZ() -= p;
        }else{
            DebugStop();
//            T trace = T(delta_eps.I1());
//            REAL lambda, G, dG_desp_vol;
//            this->G(delta_eps,G,dG_desp_vol);
//            lambda = T((2.0*G*m_nu)/(1.0-2.0*m_nu));
//            sigma.Identity();
//            sigma.Multiply(trace, lambda);
//            sigma.Add(delta_eps, 2. * G);
        }

    }
    
    template<class T>
    void ComputeStrain(const TPZTensor<T> & sigma, TPZTensor<T> & epsilon) const {

        // Initial guess is for the deviatoric is obtained from the given epsilon
        // Computing the initial guess for the volumetric part
        REAL p_star;
        p_star = -sigma.I1()/3;
        int n_iterations = 20;
        bool stop_criterion = false;
        REAL res_tol = 1.0e-5;
        REAL eps_v = (m_kappa*log((m_p_0 + m_pt_el)/(p_star + m_pt_el)))/(1 + m_e_0);
        
        epsilon.XX() = eps_v/3.0;
        epsilon.YY() = eps_v/3.0;
        epsilon.ZZ() = eps_v/3.0;

        
        n_iterations = 40;
        stop_criterion = false;
        res_tol = 1.0e-8;
        REAL corr_tol = 1.0e-8;
        TPZFNMatrix<6,STATE> res_sigma(6,1), delta_eps(6,1), eps_k(6,1);
        TPZTensor<T> delta_sigma, sigma_i;
        epsilon.CopyTo(eps_k);

        REAL res_norm = 1.0;
        REAL corr_norm = 1.0;
        int i;
        for (i = 0; i < n_iterations; i++) {
            
            ComputeStress(epsilon, sigma_i);
            delta_sigma = sigma - sigma_i;
            delta_sigma.CopyTo(res_sigma);
            res_norm = Norm(res_sigma);
            stop_criterion = res_norm < res_tol && corr_norm < corr_tol;
            if (stop_criterion) {
                break;
            }
            TPZFNMatrix<36,STATE> De(6,6,0.0),De_inv;
            STATE det;
            this->De(epsilon, De);

            De.DeterminantInverse(det, De_inv);
            if (IsZero(det)) {
                std::cout << "TPZPorousElasticResponse:: De matrix does not have an inverse." << std::endl;
                DebugStop();
            }
            
            De_inv.Multiply(res_sigma, delta_eps);
            corr_norm = Norm(delta_eps);
            eps_k += delta_eps;
            epsilon.CopyFrom(eps_k);
        }
        if ( i == n_iterations) {
            std::cout << "TPZPorousElasticResponse:: Inversion process does not converge." << std::endl;
            std::cout << "TPZPorousElasticResponse:: Residual norm   = " << res_norm << std::endl;
            std::cout << "TPZPorousElasticResponse:: Correction norm = " << corr_norm << std::endl;
            DebugStop();
        }
    
    }
    
};

#endif /* TPZPorousElasticResponse_h */

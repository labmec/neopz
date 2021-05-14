/**
 * @file
 * @brief Contains Unit Tests for methods of the material classes.
 */

//#include "pzmatrix.h"
#include "Elasticity/TPZElasticity3D.h"
#include "iostream"
#include "fstream"
#include <chrono>
#include "TPZElasticResponse.h" // linear elastic (LE)
#include "TPZElasticCriterion.h"
#include "TPZPorousElasticResponse.h" // Porous elasticity (PE)
#include "TPZPorousElasticCriterion.h"
#include "TPZSandlerExtended.h" // LE with DiMaggio Sandler (LEDS)
#include "TPZYCMohrCoulombPV.h" // LE with Mohr Coulomb (LEMC)

#include "TPZMatElastoPlastic2D_impl.h"
#include "TPZPlasticStepPV.h" // Plastic Integrator
#include <catch2/catch.hpp>
/**
 Read DiMaggio Sandler data
 Original source: An algorithm and a modular subroutine for the CAP model (April 1979)
 DOI: 10.1002/nag.1610030206
 @param file_name file containg computed data for the original cap model
 @return stress_strain uniaxial loading strain stress data
 */
TPZFMatrix<STATE> readStressStrain(std::string &file_name) {
    
    std::ifstream in(file_name.c_str());
    int n_data = 26;
    TPZFMatrix<STATE> stress_strain(n_data, 5, 0.);
    
    REAL epsilon_1, sigma_1, sigma_3, alpha_damage, m_type;
    
    int count = 0;
    while(in)
    {
        in >> epsilon_1;
        in >> sigma_1;
        in >> sigma_3;
        in >> alpha_damage;
        in >> m_type;
        
        stress_strain(count,0) = epsilon_1;
        stress_strain(count,1) = sigma_1;
        stress_strain(count,2) = sigma_3;
        stress_strain(count,3) = alpha_damage;
        stress_strain(count,4) = m_type - 1; // Because m_type = 0 , means elastic behavior
        count++;
        if (count == n_data) {
            break;
        }
    }
    return stress_strain;
}

/**
 Read pre-computed strain path and projected stresses at plane strain conditions
 @param file_name file containg computed strain and projected stresses
 @return stress_strain elastoplastic response data
 */
TPZFMatrix<STATE> readPlaneStrainPath(std::string &file_name, int n_data) {
    
    std::ifstream in(file_name.c_str());
    TPZFMatrix<STATE> strain_pv(n_data, 13, 0.);
    
    REAL s;
    REAL epsilon_xx, epsilon_yy, epsilon_zz, epsilon_xy;
    REAL epsilon_p_xx, epsilon_p_yy, epsilon_p_zz, epsilon_p_xy;
    REAL sigma_xx, sigma_yy, sigma_zz, sigma_xy;
    
    int count = 0;
    while(in)
    {
        in >> s;
        in >> epsilon_xx;
        in >> epsilon_xy;
        in >> epsilon_yy;
        in >> epsilon_zz;

        
        in >> epsilon_p_xx;
        in >> epsilon_p_xy;
        in >> epsilon_p_yy;
        in >> epsilon_p_zz;

        
        in >> sigma_xx;
        in >> sigma_xy;
        in >> sigma_yy;
        in >> sigma_zz;

        
        strain_pv(count,0) = s;
        strain_pv(count,1) = epsilon_xx;
        strain_pv(count,2) = epsilon_xy;
        strain_pv(count,3) = epsilon_yy;
        strain_pv(count,4) = epsilon_zz;

        
        strain_pv(count,5) = epsilon_p_xx;
        strain_pv(count,6) = epsilon_p_xy;
        strain_pv(count,7) = epsilon_p_yy;
        strain_pv(count,8) = epsilon_p_zz;

        
        strain_pv(count,9)  = sigma_xx;
        strain_pv(count,10) = sigma_xy;
        strain_pv(count,11) = sigma_yy;
        strain_pv(count,12) = sigma_zz;
        
        count++;
        if (count == n_data) {
            break;
        }
    }
    return strain_pv;
}

/**
 Read pre-computed strain path and projected stresses expressed in principal values
 @param file_name file containg computed strain and projected stresses in principal values
 @return stress_strain elastoplastic response data
 */
TPZFMatrix<STATE> readStrainPVPath(std::string &file_name, int n_data) {
    
    std::ifstream in(file_name.c_str());
    TPZFMatrix<STATE> strain_pv(n_data, 6, 0.);
    
    REAL epsilon_1, epsilon_2, epsilon_3;
    REAL sigma_1, sigma_2, sigma_3;
    
    int count = 0;
    while(in)
    {
        in >> epsilon_1;
        in >> epsilon_2;
        in >> epsilon_3;
        in >> sigma_1;
        in >> sigma_2;
        in >> sigma_3;
        
        strain_pv(count,0) = epsilon_1;
        strain_pv(count,1) = epsilon_2;
        strain_pv(count,2) = epsilon_3;
        strain_pv(count,3) = sigma_1;
        strain_pv(count,4) = sigma_2;
        strain_pv(count,5) = sigma_3;
        count++;
        if (count == n_data) {
            break;
        }
    }
    return strain_pv;
}

//#define PlotDataQ

/**
 * @brief Compute and compare linear elastic response
 */
void LECompareStressStrainResponse() {
    
    // The reference data
    TPZTensor<REAL> eps_t,sigma_ref;
    REAL Eyoug = 29269.0;
    REAL nu = 0.203;
    sigma_ref.Zero(); // The reference stress
    sigma_ref.XX() = -19.771203724695;
    sigma_ref.XY() = -4.86600166251039;
    sigma_ref.XZ() = -21.897007481296754;
    sigma_ref.YY() = -29.503207049715776;
    sigma_ref.ZZ() = -24.6372053872054;

    eps_t.Zero(); // The reference strain
    eps_t.XX() = -0.0003;
    eps_t.XY() = -0.0002;
    eps_t.XZ() = -0.0009;
    eps_t.YY() = -0.0007;
    eps_t.ZZ() = -0.0005;
    
    TPZElasticResponse ER;
    TPZElasticCriterion EC;

    ER.SetEngineeringData(Eyoug, nu);
    EC.SetElasticResponse(ER);
    TPZTensor<REAL> sigma, delta_sigma;
    
    TPZFMatrix<REAL> Dep(6,6,0.0);
    EC.ApplyStrainComputeSigma(eps_t, sigma, &Dep);
    delta_sigma = sigma-sigma_ref;
    REAL norm = delta_sigma.Norm();
    bool check = IsZero(norm);
    REQUIRE(check);
    return;
}

/// Integer power
int power(int a, int n) {
    int res = 1;
    while (n) {
        if (n & 1)
            res *= a;
        a *= a;
        n >>= 1;
    }
    return res;
}

/**
 * @brief Compute and compare porous elastic response
 */
void PECompareStressStrainResponse() {
    
    TPZPorousElasticResponse PER;
    
    //Testing TPZPorousElasticResponse for constan Poisson ratio
    {
        // The reference data
        TPZTensor<REAL> eps_e,sigma_ref;
        sigma_ref.Zero(); // The reference stress
        sigma_ref.XX() = -0.54395711345815;
        sigma_ref.YY() = 0.564215469327881;
        sigma_ref.ZZ() = 0.00411244624155542;

        eps_e.Zero(); // The reference strain
        eps_e.XX() = -0.000113727;
        eps_e.YY() = 0.000116224;
        
        STATE nu = 0.203;
        STATE kappa = 0.0024;
        STATE pt_el = 5.835;
        STATE e_0 = 0.34;
        STATE p_0 = 0.0;
        PER.SetPorousElasticity(kappa, pt_el, e_0, p_0);
        PER.SetPoissonRatioConstant(nu);

        TPZTensor<STATE> sigma, epsilon_target;
        PER.ComputeStress(eps_e, sigma);
        REAL sigma_norm =(sigma-sigma_ref).Norm();
        bool sigma_check = IsZero(sigma_norm);
        REQUIRE(sigma_check);
        
        epsilon_target.Zero();
        PER.ComputeStrain(sigma,epsilon_target);
        REAL esp_norm = (epsilon_target-eps_e).Norm();
        bool esp_check = IsZero(esp_norm);
        REQUIRE(esp_check);
        
        // Check for convergence
        int n_data = 10;
        TPZManVector<REAL,10> alpha(10), error(10);
        REAL s = 0.00001;
        for (int i = 0; i < n_data; i++) {
            alpha[i] = s*(1.0/power(2,(i)));
        }
        
        TPZFNMatrix<36,STATE> De(6,6);
        TPZFNMatrix<6,STATE> epsilon_v(6,1), sigma_v(6,1), delta_sigma_v(6,1), delta_eps_v(6,1);
        De.Zero();
        
        PER.De(eps_e, De);
        for(int i = 0; i < n_data; i++){

            REAL alpha_val = alpha[i];
            
            TPZTensor<STATE> sigma_approx, epsilon, delta_sigma, delta_eps;
            eps_e.CopyTo(epsilon_v);
            epsilon_v += alpha_val;
            epsilon.CopyFrom(epsilon_v);
            PER.ComputeStress(epsilon, sigma);
            
            delta_eps_v.Zero();
            delta_eps_v += alpha_val;
            De.Multiply(delta_eps_v, delta_sigma_v);
            delta_sigma.CopyFrom(delta_sigma_v);
            sigma_approx = sigma_ref + delta_sigma;
            delta_sigma = sigma - sigma_approx;
            REAL error_norm = delta_sigma.Norm();
            error[i] = log(error_norm);
            alpha[i] = log(alpha_val);
        }
        for(int i = 0; i < n_data-1; i++){
            REAL n = (error[i+1] - error[i])/(alpha[i+1] - alpha[i]);
            REAL diff = n - 1.99;
            bool rate_check = fabs(diff) < 0.01;
            REQUIRE(rate_check);
        }
    }
    
    //Testing TPZPorousElasticResponse for constan Shear modulus
    {
        // The reference data
        TPZTensor<REAL> eps_e,sigma_ref;
        sigma_ref.Zero(); // The reference stress
        sigma_ref.XX() = -0.544354891592515;
        sigma_ref.YY() = 0.564579745407486;
        sigma_ref.ZZ() = -0.00671541759251463;

        eps_e.Zero(); // The reference strain
        eps_e.XX() = -2.20978e-5;
        eps_e.YY() = 2.34811e-5;

        STATE mu = 12165.0;
        STATE kappa = 0.0024;
        STATE pt_el = 5.835;
        STATE e_0 = 0.34;
        STATE p_0 = 0.0;
        PER.SetPorousElasticity(kappa, pt_el, e_0, p_0);
        PER.SetShearModulusConstant(mu);

        TPZTensor<STATE> sigma, epsilon_target;
        PER.ComputeStress(eps_e, sigma);
        REAL sigma_norm =(sigma-sigma_ref).Norm();
        bool sigma_check = IsZero(sigma_norm);
        REQUIRE(sigma_check);

        epsilon_target.Zero();
        PER.ComputeStrain(sigma,epsilon_target);
        REAL esp_norm = (epsilon_target-eps_e).Norm();
        bool esp_check = IsZero(esp_norm);
        REQUIRE(esp_check);

        // Check for convergence
        int n_data = 10;
        TPZManVector<REAL,10> alpha(10), error(10);
        REAL s = 0.00001;
        for (int i = 0; i < n_data; i++) {
            alpha[i] = s*(1.0/power(2,(i)));
        }
        
        TPZFNMatrix<36,STATE> De(6,6);
        TPZFNMatrix<6,STATE> epsilon_v(6,1), sigma_v(6,1), delta_sigma_v(6,1), delta_eps_v(6,1);
        De.Zero();
        
        PER.De(eps_e, De);
        for(int i = 0; i < n_data; i++){
            
            REAL alpha_val = alpha[i];
            
            TPZTensor<STATE> sigma_approx, epsilon, delta_sigma, delta_eps;
            eps_e.CopyTo(epsilon_v);
            epsilon_v += alpha_val;
            epsilon.CopyFrom(epsilon_v);
            PER.ComputeStress(epsilon, sigma);
            
            delta_eps_v.Zero();
            delta_eps_v += alpha_val;
            De.Multiply(delta_eps_v, delta_sigma_v);
            delta_sigma.CopyFrom(delta_sigma_v);
            sigma_approx = sigma_ref + delta_sigma;
            delta_sigma = sigma - sigma_approx;
            REAL error_norm = delta_sigma.Norm();
            error[i] = log(error_norm);
            alpha[i] = log(alpha_val);
        }
        for(int i = 0; i < n_data-1; i++){
            REAL n = (error[i+1] - error[i])/(alpha[i+1] - alpha[i]);
            REAL diff = n - 1.99;
            bool rate_check = fabs(diff) < 0.01;
            REQUIRE(rate_check);
        }

    }
    
    /// Testing linearized elastic response with constant nu
    {
        STATE nu = 0.203;
        STATE kappa = 0.0024;
        STATE pt_el = 5.835;
        STATE e_0 = 0.34;
        STATE p_0 = 0.0;
        PER.SetPorousElasticity(kappa, pt_el, e_0, p_0);
        PER.SetPoissonRatioConstant(nu);
        
        // The reference data
        TPZTensor<REAL> eps_e,sigma_ref, sigma_linear;
        sigma_ref.Zero(); // The reference stress
        sigma_ref.XX() = -0.54395711345815;
        sigma_ref.YY() = 0.564215469327881;
        sigma_ref.ZZ() = 0.00411244624155542;
        
        eps_e.Zero(); // The reference strain
        eps_e.XX() = -0.000113727;
        eps_e.YY() = 0.000116224;
        
        TPZElasticResponse LE_equivalent = PER.EvaluateElasticResponse(eps_e);
        LE_equivalent.ComputeStress(eps_e, sigma_linear);
        REAL sigma_norm = (sigma_linear-sigma_ref).Norm();
        bool sigma_check = IsZero(sigma_norm);
        REQUIRE(sigma_check);
        
    }
    
    /// Testing linearized elastic response with constant mu (Shear G)
    {
        STATE mu = 12165.0;
        STATE kappa = 0.0024;
        STATE pt_el = 5.835;
        STATE e_0 = 0.34;
        STATE p_0 = 0.0;
        PER.SetPorousElasticity(kappa, pt_el, e_0, p_0);
        PER.SetShearModulusConstant(mu);
        
        // The reference data
        TPZTensor<REAL> eps_e,sigma_ref, sigma_linear;
        sigma_ref.Zero(); // The reference stress
        sigma_ref.XX() = -0.544354891592515;
        sigma_ref.YY() = 0.564579745407486;
        sigma_ref.ZZ() = -0.00671541759251463;
        
        eps_e.Zero(); // The reference strain
        eps_e.XX() = -2.20978e-5;
        eps_e.YY() = 2.34811e-5;
        
        TPZElasticResponse LE_equivalent = PER.EvaluateElasticResponse(eps_e);
        LE_equivalent.ComputeStress(eps_e, sigma_linear);
        REAL sigma_norm = (sigma_linear-sigma_ref).Norm();
        bool sigma_check = IsZero(sigma_norm);
        REQUIRE(sigma_check);
        
    }
    
    return;
}

/**
 * @brief Compute and compare DiMaggio Sandler elastoplastic response with DiMaggio Sandler data
 */
void LEDSCompareStressStrainAlphaMType() {
   
    std::string dirname = PZSOURCEDIR;
    std::string file_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/Sandler_Rubin_data_1979.txt";
    
    TPZFNMatrix<80,STATE> ref_epsilon_stress;
    ref_epsilon_stress = readStressStrain(file_name);
    int n_data_to_compare = 15; // @omar:: fev/2018: 13 because we do not care of tensile states
    
    TPZFNMatrix<80,STATE> LEDS_epsilon_stress(n_data_to_compare,ref_epsilon_stress.Cols());
    TPZFNMatrix<80,int> comparison(n_data_to_compare,ref_epsilon_stress.Cols());
    
    // DS Dimaggio Sandler PV
    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    // MCormick Ranch sand data:
    REAL K = 66.67; // ksi
    REAL G = 40.00; // ksi
    
    REAL E       = (9.0*K*G)/(3.0*K+G);
    REAL nu      = (3.0*K - 2.0*G)/(2.0*(3.0*K+G));
    REAL CA      = 0.250;
    REAL CB      = 0.670;
    REAL CC      = 0.180;
    REAL CD      = 0.670;
    REAL CR      = 2.500;
    REAL CW      = 0.066;
    REAL phi = 0, psi = 1., N = 0.;
    
    ER.SetEngineeringData(E, nu);
    LEDS.SetElasticResponse(ER);
    LEDS.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
    
    TPZTensor<REAL> epsilon_t,sigma;
    TPZFMatrix<REAL> source(6,1,0.0);
    sigma.Zero();
    
    // Initial damage data
    REAL k_0;
    LEDS.InitialDamage(sigma, k_0); // resolve the initial damage when reach two roots
    LEDS.fN.m_hardening = k_0;
    LEDS.fYC.SetInitialDamage(k_0);
    
    TPZPlasticState<STATE> plastic_state;
    for (int i = 0; i < n_data_to_compare; i++) {
        
        source(3,0) = ref_epsilon_stress(i,0);
        epsilon_t.CopyFrom(source);
        LEDS.ApplyStrainComputeSigma(epsilon_t, sigma);
        
        LEDS_epsilon_stress(i,0) = epsilon_t.YY();
        LEDS_epsilon_stress(i,1) = sigma.YY();
        LEDS_epsilon_stress(i,2) = sigma.XX();
        LEDS_epsilon_stress(i,3) = LEDS.fN.VolHardening();
        LEDS_epsilon_stress(i,4) = LEDS.fN.MType();
        
       
    }
    
#ifdef PlotDataQ
    LEDS_epsilon_stress.Print("LEDSdata = ",std::cout,EMathematicaInput);
#endif
    
    REAL tolerance = 1.0e-2;
    
    // Force second point comparison with expdata to zero
    LEDS_epsilon_stress(1,1) = ref_epsilon_stress(1,1);
    LEDS_epsilon_stress(1,2) = ref_epsilon_stress(1,2);
    LEDS_epsilon_stress(1,3) = ref_epsilon_stress(1,3);
    for (int i = 0; i < n_data_to_compare; i++) {
        
        for (int j = 0; j < 5; j++) {
            comparison(i,j) = fabs(LEDS_epsilon_stress(i,j) - ref_epsilon_stress(i,j)) <= tolerance;
#ifdef PZDEBUG
            if(comparison(i,j) == 0){
                std::cout << "LEDS_epsilon_stress(i,j) = " << LEDS_epsilon_stress(i,j) << std::endl;
                std::cout << "ref_epsilon_stress(i,j) = " << ref_epsilon_stress(i,j) << std::endl;
                std::cout << "comparison(i,j) = " << comparison(i,j) << std::endl;
            }
#endif
        }
    }
    
    for (int i = 0; i < comparison.Rows(); i++) {
        for (int j = 0; j < comparison.Cols(); j++) {
            bool check = comparison(i,j);
            REQUIRE(check);
        }
    }
    return;
}

/**
 * @brief Compute and compare DiMaggio Sandler elastoplastic response with erick experimental
 */
void LEDSCompareStressStrainErickTest() {
    
    std::string dirname = PZSOURCEDIR;
    std::string file_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/DS_Path_and_Erick_Stress.txt";
    
    int n_data = 10067;
    TPZFNMatrix<18,STATE> epsilon_path_proj_sigma;
    epsilon_path_proj_sigma = readStrainPVPath(file_name,n_data);
    
    TPZFNMatrix<18,STATE> LEDS_stress(n_data,4);
    TPZFNMatrix<18,int> comparison(n_data,epsilon_path_proj_sigma.Cols());
    
    // DS Dimaggio Sandler PV
    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    // Erick's rock data:
    REAL E       = 29269.0; // MPa
    REAL nu      = 0.203;   // MPa
    
    REAL K = E/(3.0*(1.0-2.0*nu));
    REAL G = E/(2.0*(1.0+nu));
    
    REAL CA      = 2.0*152.54;
    REAL CB      = 0.00155;
    REAL CC      = 0.5*146.29;
    REAL CD      = 0.01877;
    REAL CR      = 0.91969;
    REAL CW      = 0.06605;
    REAL phi = 0, psi = 1., N = 0.;
    
    ER.SetEngineeringData(E, nu);
    LEDS.SetElasticResponse(ER);
    LEDS.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
    
    TPZTensor<REAL> epsilon_t,sigma;
    TPZFMatrix<REAL> source(6,1,0.0);
    sigma.Zero();
    
    sigma.XX() = -45.8051;
    sigma.YY() = -62.1051;
    sigma.ZZ() = -48.2052;
    
    // Initial damage data
    REAL k_0;
    LEDS.InitialDamage(sigma, k_0);
    LEDS.fN.m_hardening = k_0;
    LEDS.fYC.SetInitialDamage(k_0);
    
    
    TPZPlasticState<STATE> plastic_state;
    for (int i = 0; i < n_data; i++) {
        
        source(0,0) = epsilon_path_proj_sigma(i,0);
        source(3,0) = epsilon_path_proj_sigma(i,1);
        source(5,0) = epsilon_path_proj_sigma(i,2);
        epsilon_t.CopyFrom(source);
        
        std::cout << "i = " << i << std::endl;
        LEDS.ApplyStrainComputeSigma(epsilon_t, sigma);
        
        
        
        LEDS_stress(i,0) = epsilon_t.ZZ();
        LEDS_stress(i,1) = sigma.XX();
        LEDS_stress(i,2) = sigma.YY();
        LEDS_stress(i,3) = sigma.ZZ();
        
    }
    
#ifdef PlotDataQ
    LEDS_stress.Print("LEDSdata = ",std::cout,EMathematicaInput);
#endif
    
//    for (int i = 0; i < n_data; i++) {
//        for (int j = 0; j < 3; j++) {
//            bool check = IsZero(LEDS_stress(i,j) - epsilon_path_proj_sigma(i,3+j));
//            REQUIRE(check);
//        }
//    }
    return;
}

/**
 * @brief Compute and compare DiMaggio Sandler elastoplastic response
 */
void LEDSCompareStressStrainResponse() {
    
    std::string dirname = PZSOURCEDIR;
    std::string file_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/DS_Path_and_Projected_Stress.txt";
    
    int n_data = 72;
    TPZFNMatrix<18,STATE> epsilon_path_proj_sigma;
    epsilon_path_proj_sigma = readStrainPVPath(file_name,n_data);
    
    TPZFNMatrix<18,STATE> LEDS_stress(n_data,3);
    TPZFNMatrix<18,int> comparison(n_data,epsilon_path_proj_sigma.Cols());
    
    // DS Dimaggio Sandler PV
    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    // MCormick Ranch sand data:
    REAL K = 66.67; // ksi
    REAL G = 40.00; // ksi
    
    REAL E       = (9.0*K*G)/(3.0*K+G);
    REAL nu      = (3.0*K - 2.0*G)/(2.0*(3.0*K+G));
    REAL CA      = 0.250;
    REAL CB      = 0.670;
    REAL CC      = 0.180;
    REAL CD      = 0.670;
    REAL CR      = 2.500;
    REAL CW      = 0.066;
    REAL phi = 0, psi = 1., N = 0.;
    
    ER.SetEngineeringData(E, nu);
    LEDS.SetElasticResponse(ER);
    LEDS.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
    
    TPZTensor<REAL> epsilon_t,sigma;
    TPZFMatrix<REAL> source(6,1,0.0);
    sigma.Zero();
    
    // Initial damage data
    REAL k_0;
    LEDS.InitialDamage(sigma, k_0);

    for (int i = 0; i < n_data; i++) {
        
        LEDS.fN.m_hardening = k_0;
        source(0,0) = epsilon_path_proj_sigma(i,0);
        source(3,0) = epsilon_path_proj_sigma(i,1);
        source(5,0) = epsilon_path_proj_sigma(i,2);
        epsilon_t.CopyFrom(source);
        
        LEDS.ApplyStrainComputeSigma(epsilon_t, sigma);
        
        LEDS_stress(i,0) = sigma.XX();
        LEDS_stress(i,1) = sigma.YY();
        LEDS_stress(i,2) = sigma.ZZ();
        
        LEDS.fN.m_eps_p.Zero();
        LEDS.fN.m_eps_t.Zero();
        LEDS.fN.m_hardening = 0.0;
    }
    
//#ifdef PlotDataQ
//    LEDS_stress.Print("LEDSdata = ",std::cout,EMathematicaInput);
//#endif
    
    // @omar:: update precomputed stress path
    for (int i = 0; i < n_data; i++) {
        for (int j = 0; j < 3; j++) {
            bool check = IsZero(LEDS_stress(i,j) - epsilon_path_proj_sigma(i,3+j));
            REQUIRE(check);
        }
    }
    
    return;
}

//#define VerboseMode

/**
 * @brief Compute the tangent at several points of rubin experiment.
 */
void LEDSCompareStressStrainTangent() {
    
    std::string dirname = PZSOURCEDIR;
    std::string file_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/Sandler_Rubin_data_1979.txt";
    
    TPZFNMatrix<80,STATE> ref_epsilon_stress;
    ref_epsilon_stress = readStressStrain(file_name);
    int n_data_to_compare = 17; // @omar:: fev/2018: 13 because we do not care of tensile states
    
    TPZFNMatrix<80,STATE> LEDS_epsilon_stress(n_data_to_compare,ref_epsilon_stress.Cols());
    TPZFNMatrix<80,int> comparison(n_data_to_compare,ref_epsilon_stress.Cols());
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    // MCormick Ranch sand data:
    REAL K = 66.67; // ksi
    REAL G = 40.00; // ksi
    
    REAL E       = (9.0*K*G)/(3.0*K+G);
    REAL nu      = (3.0*K - 2.0*G)/(2.0*(3.0*K+G));
    REAL CA      = 0.250;
    REAL CB      = 0.670;
    REAL CC      = 0.180;
    REAL CD      = 0.670;
    REAL CR      = 2.500;
    REAL CW      = 0.066;
    REAL phi = 0, psi = 1., N = 0.;
    
    ER.SetEngineeringData(E, nu);
    
    // Testing Cap Vertex tangent and projected stress
    {
        TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS_cap_vertex;
        REAL delta_strain = -0.005;
        TPZFMatrix<REAL> Dep(6,6,0.0);
        TPZTensor<REAL> epsilon_t,sigma;
        TPZFMatrix<REAL> source(6,1,0.0);
        sigma.Zero();
        
        LEDS_cap_vertex.SetElasticResponse(ER);
        LEDS_cap_vertex.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
        
        // Initial damage data
        REAL k_0;
        LEDS_cap_vertex.InitialDamage(sigma, k_0); // resolve the initial damage when reach two roots
        LEDS_cap_vertex.fN.m_hardening = k_0;
        LEDS_cap_vertex.fYC.SetInitialDamage(k_0);
    
        TPZPlasticState<STATE> ref_state = LEDS_cap_vertex.fN;
        
        source(0,0) = delta_strain;
        source(3,0) = delta_strain;
        source(5,0) = delta_strain;
        epsilon_t.CopyFrom(source);
        LEDS_cap_vertex.ApplyStrainComputeSigma(epsilon_t, sigma, &Dep);
        
        TPZManVector<STATE,4> alpha(4);
        alpha[0] = -1.0e-3;
        alpha[1] = -1.0e-4;
        alpha[2] = -1.0e-5;
        alpha[3] = -1.0e-6;
        
        TPZTensor<REAL> epsilon_neigh, delta_epsilon;
        TPZTensor<REAL> sigma_approx,sigma_neigh, delta_sigma, sigma_error;
        TPZFNMatrix<6,REAL> source_t(6,1,0.0),origin(6,1,0.0);
        TPZFNMatrix<6,REAL> delta_epsilon_t(6,1,0.0),delta_sigma_t(6,1,0.0);
        
        TPZFNMatrix<4,REAL> errors(4,2,0.0),rates(3,1,0.0);
        
        for (int i = 0; i < 4; i++) {
            
            LEDS_cap_vertex.fN = ref_state;
            LEDS_cap_vertex.fYC.SetInitialDamage(k_0);
            source(0,0) = delta_strain + alpha[i];
            source(3,0) = delta_strain + alpha[i];
            source(5,0) = delta_strain + alpha[i];
            epsilon_neigh.CopyFrom(source);
            LEDS_cap_vertex.ApplyStrainComputeSigma(epsilon_neigh, sigma_neigh);
            
            delta_epsilon = epsilon_neigh - epsilon_t;
            delta_epsilon.CopyTo(delta_epsilon_t);
            
            Dep.Multiply(delta_epsilon_t, delta_sigma_t);
            delta_sigma.CopyFrom(delta_sigma_t);
            sigma_approx = delta_sigma + sigma;
            
            sigma_error = sigma_neigh - sigma_approx;
            errors(i,0) = delta_epsilon.Norm();
            errors(i,1) = sigma_error.Norm();
            delta_sigma_t.Zero();
            
        }
        

        
        STATE rate;
        for (int j = 0; j < 3; j++) {
            rate = (log(errors(j,1)) - log(errors(j+1,1)))/(log(errors(j,0)) - log(errors(j+1,0)));
            rates(j,0) = rate;
            bool check = fabs(rate-2.0) < 1.0e-1;
            REQUIRE(check);
        }
#ifdef VerboseMode
        std::cout << "Cap Vertex tangent " << std::endl;
        errors.Print(std::cout);
        rates.Print(std::cout);
#endif
    }
    
    
    // Testing Cap tangent and projected stress
    {
        TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS_cap;
        REAL delta_strain = -0.005;
        TPZFMatrix<REAL> Dep(6,6,0.0);
        TPZTensor<REAL> epsilon_t,sigma;
        TPZFMatrix<REAL> source(6,1,0.0);
        sigma.Zero();
        
        LEDS_cap.SetElasticResponse(ER);
        LEDS_cap.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
        
        // Initial damage data
        REAL k_0;
        LEDS_cap.InitialDamage(sigma, k_0); // resolve the initial damage when reach two roots
        LEDS_cap.fN.m_hardening = k_0;
        LEDS_cap.fYC.SetInitialDamage(k_0);
        
        TPZPlasticState<STATE> ref_state = LEDS_cap.fN;
        
        source(0,0) = 0.0;
        source(3,0) = delta_strain;
        source(5,0) = 0.0;
        epsilon_t.CopyFrom(source);
        LEDS_cap.ApplyStrainComputeSigma(epsilon_t, sigma, &Dep);
        
        TPZManVector<STATE,4> alpha(4);
        alpha[0] = -1.0e-3;
        alpha[1] = -1.0e-4;
        alpha[2] = -1.0e-5;
        alpha[3] = -1.0e-6;
        
        TPZTensor<REAL> epsilon_neigh, delta_epsilon;
        TPZTensor<REAL> sigma_approx,sigma_neigh, delta_sigma, sigma_error;
        TPZFNMatrix<6,REAL> source_t(6,1,0.0),origin(6,1,0.0);
        TPZFNMatrix<6,REAL> delta_epsilon_t(6,1,0.0),delta_sigma_t(6,1,0.0);
        
        TPZFNMatrix<4,REAL> errors(4,2,0.0),rates(3,1,0.0);
        
        for (int i = 0; i < 4; i++) {
            
            LEDS_cap.fN = ref_state;
            LEDS_cap.fYC.SetInitialDamage(k_0);
            source(0,0) = alpha[i];
            source(3,0) = delta_strain + alpha[i];
            source(5,0) = alpha[i];
            epsilon_neigh.CopyFrom(source);
            LEDS_cap.ApplyStrainComputeSigma(epsilon_neigh, sigma_neigh);
            
            delta_epsilon = epsilon_neigh - epsilon_t;
            delta_epsilon.CopyTo(delta_epsilon_t);
            
            Dep.Multiply(delta_epsilon_t, delta_sigma_t);
            delta_sigma.CopyFrom(delta_sigma_t);
            sigma_approx = delta_sigma + sigma;
            
            sigma_error = sigma_neigh - sigma_approx;
            errors(i,0) = delta_epsilon.Norm();
            errors(i,1) = sigma_error.Norm();
            delta_sigma_t.Zero();
            
        }
        
        STATE rate;
        for (int j = 0; j < 3; j++) {
            rate = (log(errors(j,1)) - log(errors(j+1,1)))/(log(errors(j,0)) - log(errors(j+1,0)));
            rates(j,0) = rate;
            bool check = fabs(rate-2.0) < 1.0e-1;
            REQUIRE(check);
        }
#ifdef VerboseMode
        std::cout << "Cap tangent " << std::endl;
        errors.Print(std::cout);
        rates.Print(std::cout);
#endif
    }
    
    // Testing Cap Covertex tangent and projected stress
    {
        TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS_cap_covertex;
        REAL strain_xx = -0.00206416;
        REAL strain_xy = -0.00579677;
        REAL strain_yy = 0.002147444999999999;
        
        TPZFMatrix<REAL> Dep(6,6,0.0);
        TPZTensor<REAL> epsilon_t,sigma;
        TPZFMatrix<REAL> source(6,1,0.0);
        sigma.Zero();
        
        LEDS_cap_covertex.SetElasticResponse(ER);
        LEDS_cap_covertex.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
        
        // Reference data
        TPZTensor<REAL> eps_ref,eps_p_ref;
        eps_ref.Zero();
        eps_p_ref.Zero();
        
        eps_ref.XX() = -0.002434284;
        eps_ref.XY() = -0.0035826;
        eps_ref.YY() = 0.0024967469999999997;
        eps_p_ref.XX() = -0.002428475298331567;
        eps_p_ref.XY() = -0.0029359125875911608;
        eps_p_ref.YY() = 0.002077105696021569;
        eps_p_ref.ZZ() = -0.00020310684196511735;
        REAL k_0 = 0.12357552315237577;
    
        LEDS_cap_covertex.fN.m_eps_t = eps_ref;
        LEDS_cap_covertex.fN.m_eps_p = eps_p_ref;
        LEDS_cap_covertex.fN.m_hardening = k_0;
        LEDS_cap_covertex.fYC.SetInitialDamage(k_0);

        TPZPlasticState<STATE> ref_state = LEDS_cap_covertex.fN;
        
        source(0,0) = strain_xx;
        source(1,0) = strain_xy;
        source(3,0) = strain_yy;
        epsilon_t.CopyFrom(source);
        LEDS_cap_covertex.ApplyStrainComputeSigma(epsilon_t, sigma, &Dep);
        
        TPZManVector<STATE,4> alpha(4);
        alpha[0] = -1.0e-4;
        alpha[1] = -1.0e-5;
        alpha[2] = -1.0e-6;
        alpha[3] = -1.0e-7;
        
        TPZTensor<REAL> epsilon_neigh, delta_epsilon;
        TPZTensor<REAL> sigma_approx,sigma_neigh, delta_sigma, sigma_error;
        TPZFNMatrix<6,REAL> source_t(6,1,0.0),origin(6,1,0.0);
        TPZFNMatrix<6,REAL> delta_epsilon_t(6,1,0.0),delta_sigma_t(6,1,0.0);
        
        TPZFNMatrix<4,REAL> errors(4,2,0.0),rates(3,1,0.0);
        
        for (int i = 0; i < 4; i++) {
            
            LEDS_cap_covertex.fN = ref_state;
            LEDS_cap_covertex.fYC.SetInitialDamage(k_0);
            source(0,0) = strain_xx + alpha[i];
            source(1,0) = strain_xy;
            source(3,0) = strain_yy + alpha[i];
            source(5,0) = alpha[i];
            
            epsilon_neigh.CopyFrom(source);
            LEDS_cap_covertex.ApplyStrainComputeSigma(epsilon_neigh, sigma_neigh);
        
            
            delta_epsilon = epsilon_neigh - epsilon_t;
            delta_epsilon.CopyTo(delta_epsilon_t);
            
            Dep.Multiply(delta_epsilon_t, delta_sigma_t);
            delta_sigma.CopyFrom(delta_sigma_t);
            sigma_approx = delta_sigma + sigma;
            
            sigma_error = sigma_neigh - sigma_approx;
            errors(i,0) = delta_epsilon.Norm();
            errors(i,1) = sigma_error.Norm();
            delta_sigma_t.Zero();
            
        }
        
        STATE rate;
        for (int j = 0; j < 3; j++) {
            rate = (log(errors(j,1)) - log(errors(j+1,1)))/(log(errors(j,0)) - log(errors(j+1,0)));
            rates(j,0) = rate;
            bool check = fabs(rate-2.0) < 1.0e-1;
//            REQUIRE(check);.
        }
#ifdef VerboseMode
        std::cout << "Cap CoVertex tangent " << std::endl;
        errors.Print(std::cout);
        rates.Print(std::cout);
#endif
    }
    
    // Testing Failure tangent and projected stress
    {
        TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS_failure;
        REAL delta_strain = -0.052;
        TPZFMatrix<REAL> Dep(6,6,0.0);
        TPZTensor<REAL> epsilon_t,sigma;
        TPZFMatrix<REAL> source(6,1,0.0);
        sigma.Zero();
        
        LEDS_failure.SetElasticResponse(ER);
        LEDS_failure.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
        
        // Initial damage data
        REAL k_0;
        LEDS_failure.InitialDamage(sigma, k_0); // resolve the initial damage when reach two roots
        LEDS_failure.fN.m_hardening = k_0;
        LEDS_failure.fYC.SetInitialDamage(k_0);
        
        for (int i = 0; i < 16; i++) {
            source(3,0) = ref_epsilon_stress(i,0);
            epsilon_t.CopyFrom(source);
            LEDS_failure.ApplyStrainComputeSigma(epsilon_t, sigma);
        }
        
        TPZPlasticState<STATE> ref_state = LEDS_failure.fN;
        
        source(0,0) = 0.0;
        source(3,0) = delta_strain;
        source(5,0) = 0.0;
        epsilon_t.CopyFrom(source);
        LEDS_failure.ApplyStrainComputeSigma(epsilon_t, sigma, &Dep);
        
        TPZManVector<STATE,4> alpha(4);
        alpha[0] = -1.0e-3;
        alpha[1] = -1.0e-4;
        alpha[2] = -1.0e-5;
        alpha[3] = -1.0e-6;
        
        TPZTensor<REAL> epsilon_neigh, delta_epsilon;
        TPZTensor<REAL> sigma_approx,sigma_neigh, delta_sigma, sigma_error;
        TPZFNMatrix<6,REAL> source_t(6,1,0.0),origin(6,1,0.0);
        TPZFNMatrix<6,REAL> delta_epsilon_t(6,1,0.0),delta_sigma_t(6,1,0.0);
        
        TPZFNMatrix<4,REAL> errors(4,2,0.0),rates(3,1,0.0);
        
        for (int i = 0; i < 4; i++) {
            
            LEDS_failure.fN = ref_state;
            LEDS_failure.fYC.SetInitialDamage(k_0);
            source(0,0) = alpha[i];
            source(3,0) = delta_strain + alpha[i];
            source(5,0) = alpha[i];
            epsilon_neigh.CopyFrom(source);
            LEDS_failure.ApplyStrainComputeSigma(epsilon_neigh, sigma_neigh);
            
            delta_epsilon = epsilon_neigh - epsilon_t;
            delta_epsilon.CopyTo(delta_epsilon_t);
            
            Dep.Multiply(delta_epsilon_t, delta_sigma_t);
            delta_sigma.CopyFrom(delta_sigma_t);
            sigma_approx = delta_sigma + sigma;
            
            sigma_error = sigma_neigh - sigma_approx;
            errors(i,0) = delta_epsilon.Norm();
            errors(i,1) = sigma_error.Norm();
            delta_sigma_t.Zero();
            
        }
        
        STATE rate;
        for (int j = 0; j < 3; j++) {
            rate = (log(errors(j,1)) - log(errors(j+1,1)))/(log(errors(j,0)) - log(errors(j+1,0)));
            rates(j,0) = rate;
            bool check = fabs(rate-2.0) < 1.0e-1;
            REQUIRE(check);
        }
#ifdef VerboseMode
        std::cout << "Failure tangent " << std::endl;
        errors.Print(std::cout);
        rates.Print(std::cout);
#endif
    }
    
    return;
}

/**
 * @brief Compute and compare Mohr-Coulomb elastoplastic response
 */
void LEMCCompareStressStrainResponseAbaqus() {
    
    std::string dirname = PZSOURCEDIR;
    std::string file_name,file_ref_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/MC_Abaqus_S1_Path_Projected_Stress.txt";
    
    int n_data = 16;
    TPZFNMatrix<80,STATE> epsilon_path_proj_sigma;
    epsilon_path_proj_sigma = readPlaneStrainPath(file_name,n_data);
    
    int n_data_to_compare = epsilon_path_proj_sigma.Rows();
    TPZFNMatrix<80,STATE> LEMC_epsilon_stress(n_data_to_compare,7);
    
    // MC Mohr Coloumb PV
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    REAL to_Rad = M_PI/180.0;
    REAL K =  16424.8; // MPa
    REAL G =  12165.0; // MPa
    
    REAL E       = (9.0*K*G)/(3.0*K+G);
    REAL nu      = (3.0*K - 2.0*G)/(2.0*(3.0*K+G));
    REAL c = 23.3*3.0; // MPa
    REAL phi = 0.5364, psi = 0.5364;
    
    ER.SetEngineeringData(E, nu);
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(phi, psi, c, ER);
    
    
    TPZTensor<REAL> epsilon_e,sigma;
    TPZFMatrix<REAL> source(6,1,0.0);
    TPZFNMatrix<80,REAL> Dep;
    REAL epsilon_p_xx, epsilon_p_yy, epsilon_p_zz, epsilon_p_xy;
//#define _XX_ 0
//#define _XY_ 1
//#define _XZ_ 2
//#define _YY_ 3
//#define _YZ_ 4
//#define _ZZ_ 5
//    TPZFNMatrix<36,REAL> dep(6,6,0.0);
    
    for (int i = 0; i < n_data_to_compare; i++) {
        
        epsilon_p_xx = epsilon_path_proj_sigma(i,5);
        epsilon_p_xy = epsilon_path_proj_sigma(i,6);
        epsilon_p_yy = epsilon_path_proj_sigma(i,7);
        epsilon_p_zz = epsilon_path_proj_sigma(i,8);
        
        source(0,0) = epsilon_path_proj_sigma(i,1) - epsilon_p_xx; // _XX_ 0 -> 1
        source(1,0) = epsilon_path_proj_sigma(i,2) - epsilon_p_xy; // _XY_ 1 -> 2
        source(3,0) = epsilon_path_proj_sigma(i,3) - epsilon_p_yy; // _YY_ 3 -> 4
        source(5,0) = epsilon_path_proj_sigma(i,4) - epsilon_p_zz; // _ZZ_ 5 -> 6
        
        epsilon_e.CopyFrom(source);
        LEMC.ApplyStrainComputeSigma(epsilon_e, sigma);
        epsilon_e.Print(std::cout);
        sigma.Print(std::cout);
        
        LEMC_epsilon_stress(i,0) = epsilon_path_proj_sigma(i,0);
        LEMC_epsilon_stress(i,1) = sigma.XX();
        LEMC_epsilon_stress(i,2) = sigma.XY();
        LEMC_epsilon_stress(i,3) = sigma.XZ();
        LEMC_epsilon_stress(i,4) = sigma.YY();
        LEMC_epsilon_stress(i,5) = sigma.YZ();
        LEMC_epsilon_stress(i,6) = sigma.ZZ();
        
        LEMC.fN.m_eps_p.Zero();
        LEMC.fN.m_eps_t.Zero();
        LEMC.fN.m_hardening = 0.0;
        
        
    }
    epsilon_path_proj_sigma.Print("Abaqus=",std::cout,EMathematicaInput);
    LEMC_epsilon_stress.Print("LEMC=",std::cout,EMathematicaInput);
    DebugStop();
    return;
}

/**
 * @brief Compute and compare Mohr-Coulomb elastoplastic response
 */
void LEMCCompareStressStrainResponse() {
    
    std::string dirname = PZSOURCEDIR;
    std::string file_name,file_ref_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/MC_Path_and_Projected_Stress.txt";
    
    int n_data = 66;
    TPZFNMatrix<80,STATE> epsilon_path_proj_sigma;
    epsilon_path_proj_sigma = readStrainPVPath(file_name,n_data);
    
    int n_data_to_compare = epsilon_path_proj_sigma.Rows();
    TPZFNMatrix<80,STATE> LEMC_epsilon_stress(n_data_to_compare,3);
    TPZFNMatrix<80,int> comparison(n_data_to_compare,3);
    
    // MC Mohr Coloumb PV
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    REAL to_Rad = M_PI/180.0;
    REAL K =  37029.0; // MPa
    REAL G =  17784.2; // MPa
    
    REAL E       = (9.0*K*G)/(3.0*K+G);
    REAL nu      = (3.0*K - 2.0*G)/(2.0*(3.0*K+G));
    REAL c = 23.3; // MPa
    REAL phi = 30.733456130817356*to_Rad, psi = 30.733456130817356*to_Rad;
    
    ER.SetEngineeringData(E, nu);
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(phi, psi, c, ER);

    
    TPZTensor<REAL> epsilon_t,sigma;
    TPZFMatrix<REAL> source(6,1,0.0);
    TPZFNMatrix<80,REAL> Dep;
    
    for (int i = 0; i < n_data_to_compare; i++) {
        
        source(0,0) = epsilon_path_proj_sigma(i,0);
        source(3,0) = epsilon_path_proj_sigma(i,1);
        source(5,0) = epsilon_path_proj_sigma(i,2);
        epsilon_t.CopyFrom(source);
        LEMC.ApplyStrainComputeSigma(epsilon_t, sigma);
        
        LEMC_epsilon_stress(i,0) = sigma.XX();
        LEMC_epsilon_stress(i,1) = sigma.YY();
        LEMC_epsilon_stress(i,2) = sigma.ZZ();
        
        LEMC.fN.m_eps_p.Zero();
        LEMC.fN.m_eps_t.Zero();
        LEMC.fN.m_hardening = 0.0;
        
        
    }
    
    for (int i = 0; i < n_data_to_compare; i++) {
        for (int j = 0; j < 3; j++) {
            bool check = IsZero(LEMC_epsilon_stress(i,j) - epsilon_path_proj_sigma(i,3+j));
            REQUIRE(check);
        }
    }
    // const auto t1 = std::chrono::high_resolution_clock::now();
//    std::cout << "Computing plastic steps ... " << std::endl;
//    int n = 5;
//    int n_points = pow(10, n);
//    for (int i = 0; i < n_points; i++) {
//        source(3,0) = -0.0150;
//        LEMC.fN = plastic_state;
//        epsilon_t.CopyFrom(source);
//        LEMC.ApplyStrainComputeSigma(epsilon_t, sigma);
//    }
    // const auto t2 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, milli> elapsed_time = t2-t1;
    // std::cout << "Absolute Time (ms) = " << elapsed_time.count() << std::endl;
    
    return;
}

void LEMCCompareStressStrainResponse_PlasticLab_Test() {
    
    std::string dirname = PZSOURCEDIR;
    std::string file_name,file_ref_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/MC_Path_and_Projected_Stress_PlasticLab.txt";
    
    int n_data = 16;
    TPZFNMatrix<80,STATE> epsilon_path_proj_sigma;
    epsilon_path_proj_sigma = readStrainPVPath(file_name,n_data);
    
    int n_data_to_compare = epsilon_path_proj_sigma.Rows();
    TPZFNMatrix<80,STATE> LEMC_epsilon_stress(n_data_to_compare,3);
    TPZFNMatrix<80,int> comparison(n_data_to_compare,3);
    
    // MC Mohr Coloumb PV
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    REAL to_Rad = M_PI/180.0;
    REAL K =  51074.3546; // MPa
    REAL G =  15960.8923; // MPa
    
    REAL E       = (9.0*K*G)/(3.0*K+G);
    REAL nu      = (3.0*K - 2.0*G)/(2.0*(3.0*K+G));
    REAL c = 30.0;        // MPa
    REAL phi = 10.0*to_Rad, psi = 10.0*to_Rad;
    
    ER.SetEngineeringData(E, nu);
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(phi, psi, c, ER);
    
    TPZTensor<REAL> epsilon_t,sigma;
    TPZFMatrix<REAL> source(6,1,0.0);
    TPZFNMatrix<80,REAL> Dep;
    
    for (int i = 0; i < n_data_to_compare; i++) {
        
        source(0,0) = epsilon_path_proj_sigma(i,0);
        source(3,0) = epsilon_path_proj_sigma(i,1);
        source(5,0) = epsilon_path_proj_sigma(i,2);
        epsilon_t.CopyFrom(source);
        LEMC.ApplyStrainComputeSigma(epsilon_t, sigma);
        
        LEMC_epsilon_stress(i,0) = sigma.XX();
        LEMC_epsilon_stress(i,1) = sigma.YY();
        LEMC_epsilon_stress(i,2) = sigma.ZZ();
        
    }
    
//    LEMC_epsilon_stress.Print("LEMCdata = ",std::cout,EMathematicaInput);
    
    for (int i = 0; i < n_data_to_compare; i++) {
        for (int j = 0; j < 3; j++) {
            bool check = IsZero(LEMC_epsilon_stress(i,j) - epsilon_path_proj_sigma(i,3+j));
            REQUIRE(check);
        }
    }
    
    return;
}

/**
 * @brief Compute approximation rate of Tangent operator
 */
void LEMCCompareStressStrainTangent() {
    
    std::string dirname = PZSOURCEDIR;
    std::string file_name,file_ref_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/MC_Tangent_Path_and_Projected_Stress.txt";
    
    int n_data = 126;
    TPZFNMatrix<80,STATE> epsilon_path_proj_sigma;
    epsilon_path_proj_sigma = readStrainPVPath(file_name,n_data);
    
    int n_data_to_compare = epsilon_path_proj_sigma.Rows();
    TPZFNMatrix<80,STATE> LEMC_epsilon_stress(n_data_to_compare,3);
    TPZFNMatrix<80,int> comparison(n_data_to_compare,3);
    
    // MC Mohr Coloumb PV
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    REAL to_Rad = M_PI/180.0;
    REAL K =  37029.0; // MPa
    REAL G =  17784.2; // MPa
    
    REAL E       = (9.0*K*G)/(3.0*K+G);
    REAL nu      = (3.0*K - 2.0*G)/(2.0*(3.0*K+G));
    REAL c = 23.3; // MPa
    REAL phi = 30.733456130817356*to_Rad, psi = 30.733456130817356*to_Rad;
    
    ER.SetEngineeringData(E, nu);
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(phi, psi, c, ER);
    
    TPZTensor<REAL> epsilon,epsilon_neigh, delta_epsilon;
    TPZTensor<REAL> sigma,sigma_approx,sigma_neigh, delta_sigma, sigma_error;
    TPZFNMatrix<6,REAL> source(6,1,0.0),source_t(6,1,0.0),origin(6,1,0.0);
    TPZFNMatrix<6,REAL> delta_epsilon_t(6,1,0.0),delta_sigma_t(6,1,0.0);
    TPZFMatrix<REAL> Dep(6,6,0.0);
    
    TPZFNMatrix<6,REAL> errors(6,2,0.0),alpha(6,1,0.0),rates(5,1,0.0);
    
    REAL s = 1.0;
    for (int i = 0; i < 18; i++) {
        int icp = i*7;
        
//        std::cout << "Cluster number  = " << i <<  std::endl;
        source(0,0) = s*epsilon_path_proj_sigma(icp,0);
        source(3,0) = s*epsilon_path_proj_sigma(icp,1);
        source(5,0) = s*epsilon_path_proj_sigma(icp,2);
        epsilon.CopyFrom(source);
        LEMC.ApplyStrainComputeSigma(epsilon, sigma, &Dep);
        
//        Dep.Print(std::cout);
        
        LEMC.fN.m_eps_p.Zero();
        LEMC.fN.m_eps_t.Zero();
        LEMC.fN.m_hardening = 0.0;
        
        for (int j = 1; j <= 6; j++) {

            source_t(0,0) = epsilon_path_proj_sigma(icp + j ,0);
            source_t(3,0) = epsilon_path_proj_sigma(icp + j ,1);
            source_t(5,0) = epsilon_path_proj_sigma(icp + j ,2);
            epsilon_neigh.CopyFrom(source_t);

            delta_epsilon = epsilon_neigh - epsilon;
            delta_epsilon.CopyTo(delta_epsilon_t);

            Dep.Multiply(delta_epsilon_t, delta_sigma_t);
            delta_sigma.CopyFrom(delta_sigma_t);
            sigma_approx = delta_sigma + sigma;
            
            LEMC.ApplyStrainComputeSigma(epsilon_neigh, sigma_neigh);
            
            sigma_error = sigma_approx - sigma_neigh;
            errors(j-1,0) = delta_epsilon.Norm();
            errors(j-1,1) = sigma_error.Norm();

            LEMC.fN.m_eps_p.Zero();
            LEMC.fN.m_eps_t.Zero();
            LEMC.fN.m_hardening = 0.0;

        }
        Dep.Zero();
        
//        errors.Print(std::cout);
        for (int k = 0; k < 6; k++) {
            bool check = IsZero(errors(k,1));
            
            
            if (!check) {
                for (int j = 0; j < 5; j++) {
                    rates(j,0) = (log(errors(j,1)) - log(errors(j+1,1)))/(log(errors(j,0)) - log(errors(j+1,0)));
                    check = fabs(rates(j,0)-2.0) < 1.0e-1;
                    REQUIRE(check);
                }
            }else{
                REQUIRE(check);
            }
            
            
        }
        
        
    }
    return;
}



TEST_CASE("test_sandler_dimaggio", "[plasticity_tests]") {
    
    // Elastic response
    LECompareStressStrainResponse();
    PECompareStressStrainResponse();

    LEDSCompareStressStrainAlphaMType();
//    LEDSCompareStressStrainResponse(); //  Improve this test with Ericks Data
//    LEDSCompareStressStrainErickTest();
    LEDSCompareStressStrainTangent();
    
    
    // Complete
//    LEMCCompareStressStrainResponseAbaqus(); // Test projection //  Improve this test with Abaqus data
    LEMCCompareStressStrainResponse(); // Test projection
    LEMCCompareStressStrainResponse_PlasticLab_Test(); // Test projection for PlasticLab Simulated experiment
    LEMCCompareStressStrainTangent(); //  Test Tangent
    
}

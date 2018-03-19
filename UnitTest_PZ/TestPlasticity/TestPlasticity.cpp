/**
 * @file
 * @brief Contains Unit Tests for methods of the material classes.
 */

//#include "pzmatrix.h"
#include "pzelast3d.h"
#include "iostream"
#include "fstream"

#include "TPZElasticResponse.h" // linear elastic (LE)
#include "TPZPlasticStepPV.h" // Plastic Integrator
#include "TPZSandlerExtended.h" // LE with DiMaggio Sandler (LEDS)
#include "TPZYCMohrCoulombPV.h" // LE with Mohr Coulomb (LEMC)


#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz material tests

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

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
 Read pre-computed strain path and projected stresses
 @param file_name file containg computed strain and projected stresses
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

#define PlotDataQ

/**
 * @brief Compute and compare DiMaggio Sandler elastoplastic response
 */
void LEDSCompareStressStrainAlphaMType() {
   
    std::string dirname = PZSOURCEDIR;
    std::string file_name;
    file_name = dirname + "/UnitTest_PZ/TestPlasticity/StressPaths/Sandler_Rubin_data_1979.txt";
    
    TPZFNMatrix<80,STATE> ref_epsilon_stress;
    ref_epsilon_stress = readStressStrain(file_name);
    int n_data_to_compare = 13; // @omar:: fev/2018: 19 because we do not care of tensile states
    
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
    
    ER.SetUp(E, nu);
    LEDS.SetElasticResponse(ER);
    LEDS.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
    
    TPZTensor<REAL> epsilon_t,sigma;
    TPZFMatrix<REAL> source(6,1,0.0);
    TPZFNMatrix<80,REAL> Dep;
    
    sigma.Zero();
    
    // Initial damage data
    REAL k_0;
    LEDS.InitialDamage(sigma, k_0);
    LEDS.fN.fAlpha = k_0;

    TPZPlasticState<STATE> plastic_state;
    for (int i = 0; i < n_data_to_compare; i++) {

        source(3,0) = ref_epsilon_stress(i,0);
        epsilon_t.CopyFrom(source);
        LEDS.ApplyStrainComputeSigma(epsilon_t, sigma);
        
        LEDS_epsilon_stress(i,0) = epsilon_t.YY();
        LEDS_epsilon_stress(i,1) = sigma.YY();
        LEDS_epsilon_stress(i,2) = sigma.XX();
        LEDS_epsilon_stress(i,3) = LEDS.fN.Alpha();
        LEDS_epsilon_stress(i,4) = LEDS.fN.MType();
        
        if (i == 2) {
            plastic_state = LEDS.fN;
        }
        
       
    }
    
#ifdef PlotDataQ
    LEDS_epsilon_stress.Print("LEDSdata = ",std::cout,EMathematicaInput);
#endif
    
    REAL tolerance = 1.0e-2;
    
    // Force second point comparison with expdata to zero
    LEDS_epsilon_stress(1,1) = ref_epsilon_stress(1,1);
    LEDS_epsilon_stress(1,2) = ref_epsilon_stress(1,2);
    for (int i = 0; i < n_data_to_compare; i++) {
        
        for (int j = 0; j < 5; j++) {
            comparison(i,j) = fabs(LEDS_epsilon_stress(i,j) - ref_epsilon_stress(i,j)) <= tolerance;
        }
    }
    
    for (int i = 0; i < comparison.Rows(); i++) {
        for (int j = 0; j < comparison.Cols(); j++) {
            bool check = comparison(i,j);
            BOOST_CHECK(check);
        }
    }
    
    //    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
    //    std::cout << "Computing plastic steps ... " << std::endl;
    //    int n = 5;
    //    int n_points = pow(10, n);
    //    for (int i = 0; i < n_points; i++) {
    //        source(3,0) = -0.0150;
    //        LEDS.fN = plastic_state;
    //        epsilon_t.CopyFrom(source);
    //        LEDS.ApplyStrainComputeSigma(epsilon_t, sigma);
    //    }
    //    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
    //    REAL absolute_time = boost::numeric_cast<double>((t2-t1).total_milliseconds());
    //    std::cout << "Absolute Time (seconds) = " << absolute_time/1000.0 << std::endl;
    
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
    
    ER.SetUp(E, nu);
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
        
        LEMC.fN.fEpsP.Zero();
        LEMC.fN.fEpsT.Zero();
        LEMC.fN.fAlpha = 0.0;
        
        
    }
    
    for (int i = 0; i < comparison.Rows(); i++) {
        for (int j = 0; j < comparison.Cols(); j++) {
            bool check = IsZero(LEMC_epsilon_stress(i,j) - epsilon_path_proj_sigma(i,3+j));
            BOOST_CHECK(check);
        }
    }
    
//    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
//    std::cout << "Computing plastic steps ... " << std::endl;
//    int n = 5;
//    int n_points = pow(10, n);
//    for (int i = 0; i < n_points; i++) {
//        source(3,0) = -0.0150;
//        LEMC.fN = plastic_state;
//        epsilon_t.CopyFrom(source);
//        LEMC.ApplyStrainComputeSigma(epsilon_t, sigma);
//    }
//    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
//    REAL absolute_time = boost::numeric_cast<double>((t2-t1).total_milliseconds());
//    std::cout << "Absolute Time (seconds) = " << absolute_time/1000.0 << std::endl;
    
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
    
    ER.SetUp(E, nu);
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
        
        LEMC.fN.fEpsP.Zero();
        LEMC.fN.fEpsT.Zero();
        LEMC.fN.fAlpha = 0.0;
        
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

            LEMC.fN.fEpsP.Zero();
            LEMC.fN.fEpsT.Zero();
            LEMC.fN.fAlpha = 0.0;

        }
        Dep.Zero();
        
//        errors.Print(std::cout);
        for (int k = 0; k < 6; k++) {
            bool check = IsZero(errors(k,1));
            
            
            if (!check) {
                for (int j = 0; j < 5; j++) {
                    rates(j,0) = (log(errors(j,1)) - log(errors(j+1,1)))/(log(errors(j,0)) - log(errors(j+1,0)));
                    check = fabs(rates(j,0)-2.0) < 1.0e-1;
                    BOOST_CHECK(check);
                }
            }else{
                BOOST_CHECK(check);
            }
            
            
        }
        
        
    }
    return;
}

#ifdef USING_BOOST

BOOST_AUTO_TEST_SUITE(plasticity_tests)

BOOST_AUTO_TEST_CASE(test_sandler_dimaggio) {
    
    LEDSCompareStressStrainAlphaMType();
    
    // Complete
//    LEMCCompareStressStrainResponse(); // Test projection
//    LEMCCompareStressStrainTangent(); //  Test Tangent
    
}

BOOST_AUTO_TEST_SUITE_END()



#endif

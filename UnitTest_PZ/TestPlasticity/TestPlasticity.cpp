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
#include "pzsandlerextPV.h" // LE with DiMaggio Sandler (LEDS)
#include "TPZMohrCoulomb.h" // LE with Mohr Coulomb (LEMC)


#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz material tests

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#endif

/**
 Read DiMaggio Sandler data
 Original source: An algorithm and a modular subroutine for the CAP model (April 1979)
 DOI: 10.1002/nag.1610030206
 @param file_name file containg experimental data
 @return stress_strain uniaxial loading strain stress data
 */
TPZFMatrix<STATE> readStressStrain(std::string &file_name) {
    std::ifstream in(file_name.c_str());
    TPZFMatrix<STATE> stress_strain(19, 2, 0.);
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 2; j++) {
            in >> stress_strain(i, j);
        }
    }
    return stress_strain;
}


/**
 * @brief Create the stiffness matrix of a cube from -1 and 1 in cartesian coordinates, with elastic material
 * @param
 * @note The stiffness is created with the following basis function: (x,0,0), (0,x,0), (0,0,x), (y,0,0), (0,y,0), (0,0,y), (z,0,0), (0,z,0), (0,0,z) 
 * 
 */
TPZFMatrix<STATE> computeStressStrain() {
   
    std::string file_name = "Sandler_experimental_data_1979.txt";
    TPZFNMatrix<100,STATE> ref_epsilon_stress(readStressStrain(file_name));
    
    // Sandler Dimaggio PV
    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
    
    // LE
    TPZElasticResponse ER;
    
    // MCormick Ranch sand data:
    REAL K = 66.67e6;
    REAL G = 40.00e6;
    
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
    LEDS.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
    LEDS.fER.SetUp(E, nu);
    
    TPZTensor<REAL> epsilon_t;
    TPZTensor<REAL> sigma;
    
    LEDS.ApplyStrainComputeSigma(epsilon_t, sigma);
    
    
//    TPZFMatrix<STATE> stressStrain(19, 2, 0);
//    TPZManVector<STATE, 3> epsPnext(3), // plastic strain
//            epsT(3), // total strain
//            deleps(3), // Increments of strain
//            epssol(3), // elastic strain
//            deltaepsP(3),
//            sigproj(3),
//            sigtrial(3), // eigenvalues of trial stress tensor
//            deltasigma(3);
//
//    TPZSandlerDimaggioPV materialmodel(0.25, 0.67, 0.18, 0.67, 66.67, 40., 0.066, 2.5);
//
//    deleps[0] = -0.005;
//
//    for (int k = 0; k < 3; k++) {
//        deltaepsP[k] = 0.;
//        epsT[k] = 0.;
//        epsPnext[k] = 0.;
//        epssol[k] = 0.;
//    }
//
//    STATE kproj, kprev;
//    TPZFMatrix<STATE> GradSigma(3, 3);
//
//    // The point where we change from f1 (Sandler-DiMaggio) to f2 (strain-hardening cap) in the yield function
//    kproj = 0.;
//    kprev = 0.13304;
//    for (int i = 0; i < 19; i++) {
//        for (int k = 0; k < 3; k++) {
//            epssol[k] = epsT[k] - epsPnext[k];
//        }
//
//        materialmodel.ApplyStrainComputeElasticStress(epssol, sigtrial); // Computes a "trial" stress tensor sigtrial
//        materialmodel.ProjectSigmaDep(sigtrial, kprev, sigproj, kproj, GradSigma); // Determines the regime (elastic or plastic) and projects the trial stress tensor (if needed)
//
//        stressStrain(i, 0) = -epsT[0];
//        stressStrain(i, 1) = -sigproj[0];
//
//        if (i == 12) {
//            deleps[0] = -0.002;
//            deleps[0] *= -1;
//        }
//
//        for (int k = 0; k < 3; k++) {
//            deltasigma[k] = sigtrial[k] - sigproj[k];
//        }
//
//        materialmodel.ApplyStressComputeElasticStrain(deltasigma, deltaepsP);
//
//        for (int k = 0; k < 3; k++) {
//            epsPnext[k] += deltaepsP[k];
//            epsT[k] += deleps[k];
//        }
//        kprev = kproj;
//    }
//    return stressStrain;
}


#ifdef USING_BOOST

BOOST_AUTO_TEST_SUITE(plasticity_tests)

BOOST_AUTO_TEST_CASE(test_tonto) {
    int a = 3;
    BOOST_CHECK_EQUAL(a, 3);
}

BOOST_AUTO_TEST_CASE(test_sandler_dimaggio) {
    
    InitializePZLOG();
    
    std::string name = "ExpectedSandler1979.txt";
    TPZFMatrix<STATE> rightStressStrain(readStressStrain(name)), computedStressStrain(computeStressStrain());
    REAL tol = 1e-6;
    REAL dif;
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 2; j++) {
            dif = fabs(computedStressStrain(i, j) - rightStressStrain(i, j));
            BOOST_CHECK_SMALL(dif, tol);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif

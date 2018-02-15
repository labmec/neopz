/**
 * @file
 * @brief Contains Unit Tests for methods of the material classes.
 */

//#include "pzmatrix.h"
#include "pzelast3d.h"
#include "iostream"
#include "fstream"
#include "pzSandlerDimaggioPV.h"
#include "TPZMohrCoulomb.h"


#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz material tests

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#endif

/**
 * @brief Create the stiffness matrix of a cube from -1 and 1 in cartesian coordinates, with elastic material
 * @param
 * @note The stiffness is created with the following basis function: (x,0,0), (0,x,0), (0,0,x), (y,0,0), (0,y,0), (0,0,y), (z,0,0), (0,z,0), (0,0,z) 
 * 
 */
TPZFMatrix<STATE> computeStressStrain() {
    TPZFMatrix<STATE> stressStrain(19, 2, 0);

    TPZManVector<STATE, 3> epsPnext(3), // plastic strain
            epsT(3), // total strain
            deleps(3), // Increments of strain
            epssol(3), // elastic strain
            deltaepsP(3),
            sigproj(3),
            sigtrial(3), // eigenvalues of trial stress tensor
            deltasigma(3);
    
    TPZSandlerDimaggioPV materialmodel(0.25, 0.67, 0.18, 0.67, 66.67, 40., 0.066, 2.5);

    deleps[0] = -0.005;

    for (int k = 0; k < 3; k++) {
        deltaepsP[k] = 0.;
        epsT[k] = 0.;
        epsPnext[k] = 0.;
        epssol[k] = 0.;
    }

    STATE kproj, kprev;
    TPZFMatrix<STATE> GradSigma(3, 3);

    // The point where we change from f1 (Sandler-DiMaggio) to f2 (strain-hardening cap) in the yield function
    kproj = 0.;
    kprev = 0.13304;
    for (int i = 0; i < 19; i++) {
        for (int k = 0; k < 3; k++) {
            epssol[k] = epsT[k] - epsPnext[k];
        }

        materialmodel.ApplyStrainComputeElasticStress(epssol, sigtrial); // Computes a "trial" stress tensor sigtrial
        materialmodel.ProjectSigmaDep(sigtrial, kprev, sigproj, kproj, GradSigma); // Determines the regime (elastic or plastic) and projects the trial stress tensor (if needed)

        stressStrain(i, 0) = -epsT[0];
        stressStrain(i, 1) = -sigproj[0];

        if (i == 12) {
            deleps[0] = -0.002;
            deleps[0] *= -1;
        }

        for (int k = 0; k < 3; k++) {
            deltasigma[k] = sigtrial[k] - sigproj[k];
        }

        materialmodel.ApplyStressComputeElasticStrain(deltasigma, deltaepsP);

        for (int k = 0; k < 3; k++) {
            epsPnext[k] += deltaepsP[k];
            epsT[k] += deleps[k];
        }
        kprev = kproj;
    }
    return stressStrain;
}

TPZFMatrix<STATE> readStressStrain(std::string &FileName) {
    std::ifstream in(FileName.c_str());
    TPZFMatrix<STATE> stressStrain(19, 2, 0.);
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 2; j++) {
            in >> stressStrain(i, j);
        }
    }
    return stressStrain;
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
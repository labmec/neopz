/**
 * @file
 * @brief Contains Unit Tests for methods of the matrices classes.
 */

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzblockdiag.h"
#include "pzbndmat.h"
#include "pzespmat.h"
#include "pzsbndmat.h"
#include "pzsfulmat.h"
#include "pzskylnsymmat.h"
#include "pzskylmat.h"
#include "pzysmp.h"
#include "pzsysmp.h"
#include "TPZTensor.h"

#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz tensor tests

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#endif

#ifdef USING_BOOST

/// Testing Eigenvalues of a tensor

/**
 * @brief Tests the Eigenvalues/eigenvectors of a symmetric tensor A. 
 * The eigenproblem Av=wv is solved. 
 * This is the case when all eigenvalues are distinct. */
template <class TTensor, class TNumber>
void TestingEigenDecompositionThreeDistinct() {
    TTensor tensor;

    tensor.Zero();
    tensor.fData[_XX_] = 2;
    tensor.fData[_XY_] = 1;
    tensor.fData[_YY_] = 2;
    tensor.fData[_ZZ_] = 2;

    bool check = true;

    typename TTensor::TPZDecomposed decomposedTensor;
    tensor.EigenSystem(decomposedTensor);

    const TPZManVector < TNumber > w(decomposedTensor.fEigenvalues);
    const TPZManVector<TTensor, 3> eigenTensors(decomposedTensor.fEigentensors);

    double mult = 10.;
    if (sizeof (TNumber) == 4) {
        mult *= 10.;
    }
    for (unsigned int i = 0; i < 3; i++) { // for each eigenvalue
        TPZFMatrix<TNumber> eigenSpace;
        eigenTensors[i].CopyToTensor(eigenSpace);

        //std::cout << "Eigenvalue " << i+1 << ": " << w[i] << " EigenVectors: " << std::endl;
        //eigenSpace.Print(std::cout);

        // pick a vector in the eigenspace
        TPZVec< TNumber > x(3, 0.);
        for (unsigned int j = 0; j < 3; j++) {
            x[j] = TNumber(0.);
            for (unsigned int k = 0; k < 3; ++k) {
                x[j] += eigenSpace(j, k);
            }
        }
        // residual vector = -wv + Av
        TPZVec< TNumber > res(3, 0.);
        res[0] = -w[i] * x[0] + 2 * x[0] + x[1];
        res[1] = -w[i] * x[1] + x[0] + 2 * x[1];
        res[2] = -w[i] * x[2] + 2 * x[2];

        for (int j = 0; j < 3; j++) {
            if (!IsZero(TNumber(res[j] / mult))) {
                std::cout << "diff = " << res[j] << " eigenvector = " << i << " component = " << j << std::endl;
                check = false;
            }
        }
    }

    BOOST_CHECK(check);
}

/**
 * @brief Tests the Eigenvalues/eigenvectors of a symmetric tensor A. 
 * The eigenproblem Av=wv is solved. 
 * It uses the AutoFill method to create a square 3x3 matrix, which will be the representation of the tensor. */
template <class TTensor, class TNumber>
void TestingEigenDecompositionAutoFill() {
    TPZFMatrix<TNumber> ma;
    ma.AutoFill(3, 3, true);
    
    TTensor tensor(ma);

    bool check = true;

    typename TTensor::TPZDecomposed decomposedTensor;
    tensor.EigenSystem(decomposedTensor);

    const TPZManVector < TNumber > w(decomposedTensor.fEigenvalues);
    const TPZManVector<TTensor, 3> eigenTensors(decomposedTensor.fEigentensors);

    double mult = 10.;
    if (sizeof (TNumber) == 4) {
        mult *= 10.;
    }
    for (unsigned int i = 0; i < 3; i++) { // for each eigenvalue
        TPZFMatrix<TNumber> eigenSpace;
        eigenTensors[i].CopyToTensor(eigenSpace);

        //std::cout << "Eigenvalue " << i+1 << ": " << w[i] << " EigenVectors: " << std::endl;
        //eigenSpace.Print(std::cout);
        
        // pick a vector in the eigenspace
        TPZVec< TNumber > x(3, 0.);
        for (unsigned int j = 0; j < 3; j++) {
            x[j] = TNumber(0.);
            for (unsigned int k = 0; k < 3; ++k) {
                x[j] += eigenSpace(j, k);
            }
        }
        // residual vector = -wv + Av
        TPZVec< TNumber > res(3, 0.);
        for (int k = 0; k < 3; k++) {
            res[k] = -w[i] * x[k];
            for (int l = 0; l < 3; l++) {
                res[k] += ma(k, l) * x[l];
            }
        }
        for (int j = 0; j < 3; j++) {
            if (!IsZero(TNumber(res[j] / mult))) {
                std::cout << "diff = " << res[j] << " eigenvector = " << i << " component = " << j << std::endl;
                check = false;
            }
        }
    }

    if (!check) {
        ma.Print("Matrix = ", std::cout, EMathematicaInput);
    }
    BOOST_CHECK(check);
}

/**
 * @brief Tests the Eigenvalues/eigenvectors of a null tensor A.
 * The eigenproblem Av=wv is solved. */
template <class TTensor, class TNumber>
void TestingEigenDecompositionTensorZero() {
    TTensor tensor;
    
    tensor.Zero();
    
    bool check = true;
    
    typename TTensor::TPZDecomposed decomposedTensor;
    tensor.EigenSystem(decomposedTensor);
    
    const TPZManVector < TNumber > w(decomposedTensor.fEigenvalues);
    const TPZManVector<TTensor, 3> eigenTensors(decomposedTensor.fEigentensors);
    
    double mult = 10.;
    if (sizeof (TNumber) == 4) {
        mult *= 10.;
    }
    for (unsigned int i = 0; i < 3; i++) { // for each eigenvalue
        TPZFMatrix<TNumber> eigenSpace;
        eigenTensors[i].CopyToTensor(eigenSpace);
        
        //std::cout << "Eigenvalue " << i+1 << ": " << w[i] << " EigenVectors: " << std::endl;
        //eigenSpace.Print(std::cout);
        
        // pick a vector in the eigenspace
        TPZVec< TNumber > x(3, 0.);
        for (unsigned int j = 0; j < 3; j++) {
            x[j] = TNumber(0.);
            for (unsigned int k = 0; k < 3; ++k) {
                x[j] += eigenSpace(j, k);
            }
        }
        // residual vector = -wv + v
        TPZVec< TNumber > res(3, 0.);
        for (int k = 0; k < 3; k++) {
            res[k] = -w[i] * x[k] + 0.0 * x[k];
        }
        for (int j = 0; j < 3; j++) {
            if (!IsZero(TNumber(res[j] / mult))) {
                std::cout << "diff = " << res[j] << " eigenvector = " << i << " component = " << j << std::endl;
                check = false;
            }
        }
    }
    
    BOOST_CHECK(check);
}

/**
 * @brief Tests the Eigenvalues/eigenvectors of a symmetric tensor A. 
 * The eigenproblem Av=wv is solved. */
template <class TTensor, class TNumber>
void TestingEigenDecompositionHydrostatic() {
    TTensor tensor;

    tensor.Identity();

    bool check = true;

    typename TTensor::TPZDecomposed decomposedTensor;
    tensor.EigenSystem(decomposedTensor);

    const TPZManVector < TNumber > w(decomposedTensor.fEigenvalues);
    const TPZManVector<TTensor, 3> eigenTensors(decomposedTensor.fEigentensors);

    double mult = 10.;
    if (sizeof (TNumber) == 4) {
        mult *= 10.;
    }
    for (unsigned int i = 0; i < 3; i++) { // for each eigenvalue
        TPZFMatrix<TNumber> eigenSpace;
        eigenTensors[i].CopyToTensor(eigenSpace);

        //std::cout << "Eigenvalue " << i+1 << ": " << w[i] << " EigenVectors: " << std::endl;
        //eigenSpace.Print(std::cout);
        
        // pick a vector in the eigenspace
        TPZVec< TNumber > x(3, 0.);
        for (unsigned int j = 0; j < 3; j++) {
            x[j] = TNumber(0.);
            for (unsigned int k = 0; k < 3; ++k) {
                x[j] += eigenSpace(j, k);
            }
        }
        // residual vector = -wv + v
        TPZVec< TNumber > res(3, 0.);
        for (int k = 0; k < 3; k++) {
            res[k] = -w[i] * x[k] + x[k];
        }
        for (int j = 0; j < 3; j++) {
            if (!IsZero(TNumber(res[j] / mult))) {
                std::cout << "diff = " << res[j] << " eigenvector = " << i << " component = " << j << std::endl;
                check = false;
            }
        }
    }

    BOOST_CHECK(check);
}

/**
 * @brief Tests the Eigenvalues/eigenvectors of a symmetric tensor A. 
 * The eigenproblem Av=wv is solved. 
 * This is the case when 2 of the eigenvalues are the same. */
template <class TTensor, class TNumber>
void TestingEigenDecompositionTwoEigenValues() {
    TTensor tensor;

    tensor.Zero();
    tensor.fData[_XX_] = 2;
    tensor.fData[_XY_] = 1;
    tensor.fData[_YY_] = 2;
    tensor.fData[_ZZ_] = 1;

    bool check = true;

    typename TTensor::TPZDecomposed decomposedTensor;
    tensor.EigenSystem(decomposedTensor);

    const TPZManVector < TNumber > w(decomposedTensor.fEigenvalues);
    const TPZManVector<TTensor, 3> eigenTensors(decomposedTensor.fEigentensors);

    double mult = 10.;
    if (sizeof (TNumber) == 4) {
        mult *= 10.;
    }
    for (unsigned int i = 0; i < 3; i++) { // for each eigenvalue
        TPZFMatrix<TNumber> eigenSpace;
        eigenTensors[i].CopyToTensor(eigenSpace);

        //std::cout << "Eigenvalue " << i+1 << ": " << w[i] << " EigenVectors: " << std::endl;
        //eigenSpace.Print(std::cout);

        // pick a vector in the eigenspace
        TPZVec< TNumber > x(3, 0.);
        for (unsigned int j = 0; j < 3; j++) {
            x[j] = TNumber(0.);
            for (unsigned int k = 0; k < 3; ++k) {
                x[j] += eigenSpace(j, k);
            }
        }
        // residual vector = -wv + Av
        TPZVec< TNumber > res(3, 0.);
        res[0] = -w[i] * x[0] + 2 * x[0] + x[1];
        res[1] = -w[i] * x[1] + x[0] + 2 * x[1];
        res[2] = -w[i] * x[2] + x[2];

        for (int j = 0; j < 3; j++) {
            if (!IsZero(TNumber(res[j] / mult))) {
                std::cout << "diff = " << res[j] << " eigenvector = " << i << " component = " << j << std::endl;
                check = false;
            }
        }
    }

    BOOST_CHECK(check);
}

/**
 * @brief Tests without shear stresses Eigenvalues/eigenvectors of a symmetric tensor A.
 * The eigenproblem Av=wv is solved.
 * This is the case when 2 of the eigenvalues are the same. */
template <class TTensor, class TNumber>
void TestingEigenDecompositionTwoEigenValuesNoShearStress() {
    TTensor tensor;
    
    tensor.Zero();
    tensor.fData[_XX_] = -0.21322009032156708;
    tensor.fData[_YY_] = -0.71941695995021415;
    tensor.fData[_ZZ_] = -0.21322009032156708;
    
    bool check = true;
    
    typename TTensor::TPZDecomposed decomposedTensor;
    tensor.EigenSystem(decomposedTensor);
    
    const TPZManVector < TNumber > w(decomposedTensor.fEigenvalues);
    const TPZManVector<TTensor, 3> eigenTensors(decomposedTensor.fEigentensors);
    
    double mult = 10.;
    if (sizeof (TNumber) == 4) {
        mult *= 10.;
    }
    for (unsigned int i = 0; i < 3; i++) { // for each eigenvalue
        TPZFMatrix<TNumber> eigenSpace;
        eigenTensors[i].CopyToTensor(eigenSpace);
        
        //std::cout << "Eigenvalue " << i+1 << ": " << w[i] << " EigenVectors: " << std::endl;
        //eigenSpace.Print(std::cout);
        
        // pick a vector in the eigenspace
        TPZVec< TNumber > x(3, 0.);
        for (unsigned int j = 0; j < 3; j++) {
            x[j] = TNumber(0.);
            for (unsigned int k = 0; k < 3; ++k) {
                x[j] += eigenSpace(j, k);
            }
        }
        // residual vector = -wv + Av
        TPZVec< TNumber > res(3, 0.);
        res[0] = -w[i] * x[0] + tensor.fData[_XX_] * x[0];
        res[1] = -w[i] * x[1] + tensor.fData[_YY_] * x[1];
        res[2] = -w[i] * x[2] + tensor.fData[_ZZ_] * x[2];
        
        for (int j = 0; j < 3; j++) {
            if (!IsZero(TNumber(res[j] / mult))) {
                std::cout << "diff = " << res[j] << " eigenvector = " << i << " component = " << j << std::endl;
                check = false;
            }
        }
    }
    
    BOOST_CHECK(check);
}

BOOST_AUTO_TEST_SUITE(tensor_tests)


BOOST_AUTO_TEST_CASE(eigenvalue_tests) {
    
    TestingEigenDecompositionThreeDistinct<TPZTensor<double >, double >();
    TestingEigenDecompositionThreeDistinct<TPZTensor<float>, float>();
    TestingEigenDecompositionAutoFill<TPZTensor<double >, double >();
    TestingEigenDecompositionAutoFill<TPZTensor<float>, float>();
    TestingEigenDecompositionTensorZero<TPZTensor<double >, double >();
    TestingEigenDecompositionTensorZero<TPZTensor<float >, float >();
    TestingEigenDecompositionHydrostatic<TPZTensor<double >, double >();
    TestingEigenDecompositionHydrostatic<TPZTensor<float>, float>();
    TestingEigenDecompositionTwoEigenValues<TPZTensor<double >, double >();
    TestingEigenDecompositionTwoEigenValues<TPZTensor<float>, float>();
    TestingEigenDecompositionTwoEigenValuesNoShearStress<TPZTensor<double >, double >();
    TestingEigenDecompositionTwoEigenValuesNoShearStress<TPZTensor<float>, float>();
}

BOOST_AUTO_TEST_SUITE_END()

#endif

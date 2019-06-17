//
// Created by natalia on 17/05/19.
//
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElastoPlasticMem.h"
#include "TPZMatElastoPlastic2D.h"

#ifdef USING_CUDA
#include "TPZVecGPU.h"
#include "TPZCudaCalls.h"
#endif

#ifndef INTPOINTSFEM_TPZMYLAMBDAEXPRESSION_H
#define INTPOINTSFEM_TPZMYLAMBDAEXPRESSION_H

class TPZMyLambdaExpression {

public:
    TPZMyLambdaExpression();

    TPZMyLambdaExpression(int npts, TPZVec<REAL> weight, TPZMaterial *material);

    ~TPZMyLambdaExpression();

    TPZMyLambdaExpression(const TPZMyLambdaExpression &copy);

    TPZMyLambdaExpression &operator=(const TPZMyLambdaExpression &copy);

    void SetIntPoints(int64_t npts);

    void SetWeightVector(TPZVec<REAL> weight);

    void SetMaterial(TPZMaterial *material);

    void ElasticStrain(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &elastic_strain);

    void ComputeStress(TPZFMatrix<REAL> &elastic_strain, TPZFMatrix<REAL> &sigma);

    void SpectralDecomposition(TPZFMatrix<REAL> &sigma_trial, TPZFMatrix<REAL> &eigenvalues, TPZFMatrix<REAL> &eigenvectors);

    void ProjectSigma(TPZFMatrix<REAL> &eigenvalues, TPZFMatrix<REAL> &sigma_projected);

    void StressCompleteTensor(TPZFMatrix<REAL> &sigma_projected, TPZFMatrix<REAL> &eigenvectors, TPZFMatrix<REAL> &sigma);

    void ComputeStrain(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &elastic_strain);

    void PlasticStrain(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &elastic_strainS);

    void ComputeSigma(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &sigma);

#ifdef USING_CUDA
    void ComputeSigma(TPZVecGPU<REAL> &delta_strain, TPZVecGPU<REAL> &sigma);
#endif

private:
    int64_t fNpts;

    TPZVec<REAL> fWeight;

    TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse>, TPZElastoPlasticMem> *fMaterial;

    TPZFMatrix<REAL> fPlasticStrain;

    TPZFMatrix<REAL> fMType;

    TPZFMatrix<REAL> fAlpha;

#ifdef USING_CUDA
    TPZCudaCalls *fCudaCalls;
    TPZVecGPU<REAL> dWeight;
#endif

};


#endif //INTPOINTSFEM_TPZMYLAMBDAEXPRESSION_H

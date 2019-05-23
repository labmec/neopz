//
// Created by natalia on 17/05/19.
//

#ifndef INTPOINTSFEM_TPZMYLAMBDAEXPRESSION_H
#define INTPOINTSFEM_TPZMYLAMBDAEXPRESSION_H

#include "TPZIntPointsStructMatrix.h"

class TPZMyLambdaExpression {

public:
    TPZMyLambdaExpression();

    TPZMyLambdaExpression(TPZIntPointsStructMatrix *intPointsStructMatrix, int imat);

    ~TPZMyLambdaExpression();

    TPZMyLambdaExpression(const TPZMyLambdaExpression &copy);

    TPZMyLambdaExpression &operator=(const TPZMyLambdaExpression &copy);

    void SetIntPoints(TPZIntPointsStructMatrix *IntPoints) {
        fIntPointsStructMatrix = IntPoints;
    }

    void SetMaterialId(int imat) {
        fIMat = imat;
        TPZMaterial *material = fIntPointsStructMatrix->Mesh()->FindMaterial(imat + 1);
        fMaterial = dynamic_cast<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> *>(material);
    }

    void ElasticStrain(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &elastic_strain);

    void ComputeStress(TPZFMatrix<REAL> &elastic_strain, TPZFMatrix<REAL> &sigma);

    void SpectralDecomposition(TPZFMatrix<REAL> &sigma_trial, TPZFMatrix<REAL> &eigenvalues, TPZFMatrix<REAL> &eigenvectors);

    void ProjectSigma(TPZFMatrix<REAL> &eigenvalues, TPZFMatrix<REAL> &sigma_projected);

    void StressCompleteTensor(TPZFMatrix<REAL> &sigma_projected, TPZFMatrix<REAL> &eigenvectors, TPZFMatrix<REAL> &sigma);

    void ComputeStrain(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &elastic_strain);

    void PlasticStrain(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &elastic_strainS);

    void ComputeSigma(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &sigma);

protected:
    int fIMat;

    int fNpts;

    TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse>, TPZElastoPlasticMem> *fMaterial;

    TPZFMatrix<REAL> fPlasticStrain;

    TPZFMatrix<REAL> fMType;

    TPZFMatrix<REAL> fAlpha;

    TPZIntPointsStructMatrix *fIntPointsStructMatrix;

};


#endif //INTPOINTSFEM_TPZMYLAMBDAEXPRESSION_H

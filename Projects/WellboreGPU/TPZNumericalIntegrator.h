//
// Created by natalia on 24/05/19.
//

#include "TPZIrregularBlocksMatrix.h"
#include "TPZConstitutiveLawProcessor.h"

#ifdef USING_CUDA
#include "TPZVecGPU.h"
#include "TPZCudaCalls.h"
#endif


#ifdef USING_MKL
#include <mkl.h>
#endif

#ifndef INTPOINTSFEM_TPZNUMERICALINTEGRATOR_H
#define INTPOINTSFEM_TPZNUMERICALINTEGRATOR_H


class TPZNumericalIntegrator {

public:

    TPZNumericalIntegrator();

    TPZNumericalIntegrator(TPZIrregularBlocksMatrix &irregularBlocksMatrix);

    ~TPZNumericalIntegrator();

    void Multiply(TPZFMatrix<REAL> &coef, TPZFMatrix<REAL> &delta_strain);

    void MultiplyTranspose(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &res);

    void ResidualIntegration(TPZFMatrix<REAL> & solution ,TPZFMatrix<REAL> &rhs);

    void ResidualIntegration(TPZFMatrix<REAL> & solution ,TPZFMatrix<REAL> &rhs, TPZFMatrix<REAL> &sigma);

    void ComputeConstitutiveMatrix(int64_t point_index, TPZFMatrix<STATE> &De);

    void ComputeTangentMatrix(int64_t iel, TPZFMatrix<REAL> &K);

    void ComputeTangentMatrix(int64_t iel,  TPZFMatrix<REAL> &Dep, TPZFMatrix<REAL> &K);

#ifdef USING_CUDA
    void Multiply(TPZVecGPU<REAL> &coef, TPZVecGPU<REAL> &delta_strain);
    
    void MultiplyTranspose(TPZVecGPU<REAL> &sigma, TPZVecGPU<REAL> &res); 

    void ResidualIntegration(TPZFMatrix<REAL> & solution ,TPZVecGPU<REAL> &rhs);
#endif

    void SetIrregularBlocksMatrix(TPZIrregularBlocksMatrix & irregularBlocksMatrix) {
        fBlockMatrix = irregularBlocksMatrix;
    }

    TPZIrregularBlocksMatrix & IrregularBlocksMatrix() {
        return fBlockMatrix;
    }

    void SetConstitutiveLawProcessor(TPZConstitutiveLawProcessor & processor){
        fConstitutiveLawProcessor = processor;
    }

    TPZConstitutiveLawProcessor & ConstitutiveLawProcessor(){
        return fConstitutiveLawProcessor;
    }

    void SetDoFIndexes(TPZVec<int> dof_indexes) {
        fDoFIndexes = dof_indexes;
    }

    TPZVec<int> & DoFIndexes() {
        return fDoFIndexes;
    }

    void SetColorIndexes(TPZVec<int> color_indexes) {
        fColorIndexes = color_indexes;
    }
    
    TPZVec<int> & ColorIndexes() {
     return fColorIndexes;
    }

    void SetNColors(int ncolor) {
        fNColor = ncolor;
    }

    int NColors() {
        return fNColor;
    }

#ifdef USING_CUDA
    TPZVecGPU<int> & DoFIndexesDev() {
        return dDoFIndexes;
    }

    TPZVecGPU<int> & ColorIndexesDev() {
        return dColorIndexes;
    }

    void TransferDataToGPU();

#endif

private:

    /// Irregular block matrix containing spatial gradients for scalar basis functions of order k
    TPZIrregularBlocksMatrix fBlockMatrix;

    /// Number of colors grouping no adjacent elements
    int64_t fNColor; //needed to do the assembly

    /// Degree of Freedom indexes organized element by element with stride ndof
    TPZVec<int> fDoFIndexes; // needed to do the gather operation

    /// Color indexes organized element by element with stride ndof
    TPZVec<int> fColorIndexes; //nedeed to scatter operation

    TPZConstitutiveLawProcessor fConstitutiveLawProcessor;
    
#ifdef USING_CUDA
    TPZVecGPU<int> dDoFIndexes;
    TPZVecGPU<int> dColorIndexes;
    TPZCudaCalls fCudaCalls;
#endif


};


#endif //INTPOINTSFEM_TPZNUMERICALINTEGRATOR_H

//
// Created by natalia on 24/05/19.
//

#include "TPZCoefToGradSol.h"
#include "pzcmesh.h"

TPZCoefToGradSol::TPZCoefToGradSol() : fBlockMatrix(0,0), fNColor(-1), fIndexes(0), fIndexesColor(0) {
#ifdef USING_CUDA
    dIndexes.resize(0);
    dIndexesColor.resize(0);
#endif
}

TPZCoefToGradSol::TPZCoefToGradSol(TPZIrregularBlocksMatrix &irregularBlocksMatrix) : fBlockMatrix(0,0), fNColor(-1), fIndexes(0), fIndexesColor(0) {
    SetIrregularBlocksMatrix(irregularBlocksMatrix);
}

TPZCoefToGradSol::~TPZCoefToGradSol() {

}

void TPZCoefToGradSol::SetIrregularBlocksMatrix(TPZIrregularBlocksMatrix & irregularBlocksMatrix) {
    fBlockMatrix = irregularBlocksMatrix;
}

#ifdef USING_CUDA
void TPZCoefToGradSol::Multiply(TPZVecGPU<REAL> &coef, TPZVecGPU<REAL> &grad_u) {
    int dim = 2;
    int64_t rows = fBlockMatrix.Rows();
    int64_t cols = fBlockMatrix.Cols();

    TPZVecGPU<REAL> gather_solution(dim * cols);
    fCudaCalls.GatherOperation(dim * cols, coef.getData(), gather_solution.getData(), dIndexes.getData());

    grad_u.resize(dim * rows);

    fBlockMatrix.MultiplyVector(&gather_solution.getData()[0], &grad_u.getData()[0], false);
    fBlockMatrix.MultiplyVector(&gather_solution.getData()[cols], &grad_u.getData()[rows], false);   
}
#endif 

void TPZCoefToGradSol::Multiply(TPZFMatrix<REAL> &coef, TPZFMatrix<REAL> &grad_u) {
    int dim = 2;
    int64_t rows = fBlockMatrix.Rows();
    int64_t cols = fBlockMatrix.Cols();

    TPZFMatrix<REAL> gather_solution(dim * cols, 1);
    cblas_dgthr(dim * cols, coef, &gather_solution(0, 0), &fIndexes[0]);

    grad_u.Resize(dim * rows, 1);

    fBlockMatrix.MultiplyVector(&gather_solution(0,0), &grad_u(0, 0), false);
    fBlockMatrix.MultiplyVector(&gather_solution(cols,0), &grad_u(rows, 0), false);

}

#ifdef USING_CUDA
void TPZCoefToGradSol::MultiplyTranspose(TPZVecGPU<REAL> &sigma, TPZVecGPU<REAL> &res) {
    int dim = 2;
    int64_t rows = fBlockMatrix.Rows();
    int64_t cols = fBlockMatrix.Cols();

    int64_t ncolor = fNColor;
    int64_t neq = res.getSize();    

    TPZVecGPU<REAL> forces(dim * cols);
    res.resize(ncolor * neq);
    res.Zero();

    fBlockMatrix.MultiplyVector(&sigma.getData()[0], &forces.getData()[0], true);
    fBlockMatrix.MultiplyVector(&sigma.getData()[rows], &forces.getData()[cols], true); 

    // Assemble forces
    fCudaCalls.ScatterOperation(dim * cols, forces.getData(), res.getData(), dIndexesColor.getData());

    int64_t colorassemb = ncolor / 2.;
    while (colorassemb > 0) {

        int64_t firsteq = (ncolor - colorassemb) * neq;
        fCudaCalls.DaxpyOperation(colorassemb * neq, 1., &res.getData()[firsteq], &res.getData()[0]); 

        ncolor -= colorassemb;
        colorassemb = ncolor/2;
    }
    // res.resize(neq);
}
#endif


void TPZCoefToGradSol::MultiplyTranspose(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &res) {
    int dim = 2;
    int64_t rows = fBlockMatrix.Rows();
    int64_t cols = fBlockMatrix.Cols();

    int64_t ncolor = fNColor;
    int64_t neq = res.Rows();

    TPZFMatrix<REAL> forces(dim * cols, 1);
    res.Resize(ncolor * neq, 1);
    res.Zero();

    fBlockMatrix.MultiplyVector(&sigma(0, 0), &forces(0, 0), true);
    fBlockMatrix.MultiplyVector(&sigma(rows, 0), &forces(cols, 0), true);

    // Assemble forces
    cblas_dsctr(dim * cols, forces, &fIndexesColor[0], &res(0,0));

    int64_t colorassemb = ncolor / 2.;
    while (colorassemb > 0) {

        int64_t firsteq = (ncolor - colorassemb) * neq;
        cblas_daxpy(colorassemb * neq, 1., &res(firsteq, 0), 1., &res(0, 0), 1.);

        ncolor -= colorassemb;
        colorassemb = ncolor/2;
    }
    res.Resize(neq, 1);
}

void TPZCoefToGradSol::TransferDataToGPU() {
#ifdef USING_CUDA
    fBlockMatrix.TransferDataToGPU();

    dIndexes.resize(fIndexes.size());
    dIndexes.set(&fIndexes[0], fIndexes.size());

    dIndexesColor.resize(fIndexesColor.size());
    dIndexesColor.set(&fIndexesColor[0], fIndexesColor.size());
#endif
}

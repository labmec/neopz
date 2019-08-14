//
// Created by natalia on 24/05/19.
//

#include "TPZNumericalIntegrator.h"
#include "pzcmesh.h"
#include "Timer.h"

TPZNumericalIntegrator::TPZNumericalIntegrator() : fBlockMatrix(0,0), fNColor(-1), fDoFIndexes(0), fColorIndexes(0), fConstitutiveLawProcessor(){
#ifdef USING_CUDA
    dDoFIndexes.resize(0);
    dColorIndexes.resize(0);
#endif
}

TPZNumericalIntegrator::TPZNumericalIntegrator(TPZIrregularBlocksMatrix &irregularBlocksMatrix) : fBlockMatrix(0,0), fNColor(-1), fDoFIndexes(0), fColorIndexes(0), fConstitutiveLawProcessor() {
    SetIrregularBlocksMatrix(irregularBlocksMatrix);
}

TPZNumericalIntegrator::~TPZNumericalIntegrator() {

}

#ifdef USING_CUDA
void TPZNumericalIntegrator::Multiply(TPZVecGPU<REAL> &coef, TPZVecGPU<REAL> &delta_strain) {

    int64_t rows = fBlockMatrix.Rows();
    int64_t cols = fBlockMatrix.Cols();

    TPZVecGPU<REAL> gather_solution(cols);
    gather_solution.Zero();
    fCudaCalls.GatherOperation(cols, coef.getData(), gather_solution.getData(), dDoFIndexes.getData());

    delta_strain.resize(rows);
    delta_strain.Zero();
    fBlockMatrix.MultiplyVector(&gather_solution.getData()[0], &delta_strain.getData()[0], false);
}
#endif 

void TPZNumericalIntegrator::Multiply(TPZFMatrix<REAL> &coef, TPZFMatrix<REAL> &delta_strain) {

    int64_t rows = fBlockMatrix.Rows();
    int64_t cols = fBlockMatrix.Cols();

    TPZFMatrix<REAL> gather_solution(cols, 1);
    gather_solution.Zero();
    cblas_dgthr(cols, coef, &gather_solution(0, 0), &fDoFIndexes[0]);

    delta_strain.Resize(rows, 1);
    fBlockMatrix.MultiplyVector(&gather_solution(0,0), &delta_strain(0, 0), false);
}

#ifdef USING_CUDA
void TPZNumericalIntegrator::MultiplyTranspose(TPZVecGPU<REAL> &sigma, TPZVecGPU<REAL> &res) {
    int64_t cols = fBlockMatrix.Cols();

    int64_t ncolor = fNColor;
    int64_t neq = res.getSize();    

    TPZVecGPU<REAL> forces(cols);
    forces.Zero();
    res.resize(ncolor * neq);
    res.Zero();

    fBlockMatrix.MultiplyVector(&sigma.getData()[0], &forces.getData()[0], true);

    fCudaCalls.ScatterOperation(cols, forces.getData(), res.getData(), dColorIndexes.getData());

    int64_t colorassemb = ncolor / 2.;
    while (colorassemb > 0) {

        int64_t firsteq = (ncolor - colorassemb) * neq;
        fCudaCalls.DaxpyOperation(colorassemb * neq, 1., &res.getData()[firsteq], &res.getData()[0]); 

        ncolor -= colorassemb;
        colorassemb = ncolor/2;
    }
}
#endif

void TPZNumericalIntegrator::MultiplyTranspose(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &res) {
    int64_t cols = fBlockMatrix.Cols();

    int64_t ncolor = fNColor;
    int64_t neq = res.Rows();

    TPZFMatrix<REAL> forces(cols, 1);
    res.Resize(ncolor * neq, 1);
    res.Zero();

    fBlockMatrix.MultiplyVector(&sigma(0, 0), &forces(0, 0), true);

    cblas_dsctr(cols, forces, &fColorIndexes[0], &res(0,0));

    int64_t colorassemb = ncolor / 2.;
    while (colorassemb > 0) {

        int64_t firsteq = (ncolor - colorassemb) * neq;
        cblas_daxpy(colorassemb * neq, 1., &res(firsteq, 0), 1., &res(0, 0), 1.);

        ncolor -= colorassemb;
        colorassemb = ncolor/2;
    }
    res.Resize(neq, 1);
}

void TPZNumericalIntegrator::ResidualIntegration(TPZFMatrix<REAL> & solution ,TPZFMatrix<REAL> &rhs) {
    TPZFMatrix<REAL> delta_strain;
    TPZFMatrix<REAL> sigma;

    Multiply(solution, delta_strain);
    fConstitutiveLawProcessor.ComputeSigma(delta_strain, sigma);
    MultiplyTranspose(sigma, rhs); // Perform Residual integration using a global linear application B
}

void TPZNumericalIntegrator::ResidualIntegration(TPZFMatrix<REAL> & solution ,TPZFMatrix<REAL> &rhs, TPZFMatrix<REAL> &dep) {
    TPZFMatrix<REAL> delta_strain;
    TPZFMatrix<REAL> sigma;

    Multiply(solution, delta_strain);
    fConstitutiveLawProcessor.ComputeSigmaDep(delta_strain, sigma, dep);
    MultiplyTranspose(sigma, rhs); // Perform Residual integration using a global linear application B
}

#ifdef USING_CUDA
void TPZNumericalIntegrator::ResidualIntegration(TPZFMatrix<REAL> & solution ,TPZVecGPU<REAL> &rhs) {

    TPZVecGPU<REAL> d_solution(solution.Rows());
    d_solution.set(&solution(0,0), solution.Rows());

    TPZVecGPU<REAL> d_delta_strain;
    TPZVecGPU<REAL> d_sigma;

    Multiply(d_solution, d_delta_strain);
    fConstitutiveLawProcessor.ComputeSigma(d_delta_strain, d_sigma);
    MultiplyTranspose(d_sigma, rhs);
}
#endif

#ifdef USING_CUDA
void TPZNumericalIntegrator::TransferDataToGPU() {
    fBlockMatrix.TransferDataToGPU();
    fConstitutiveLawProcessor.TransferDataToGPU();

    dDoFIndexes.resize(fDoFIndexes.size());
    dDoFIndexes.set(&fDoFIndexes[0], fDoFIndexes.size());

    dColorIndexes.resize(fColorIndexes.size());
    dColorIndexes.set(&fColorIndexes[0], fColorIndexes.size());
}
#endif


void TPZNumericalIntegrator::ComputeConstitutiveMatrix(int64_t point_index, TPZFMatrix<STATE> &De){
    
    De.Zero();
    REAL lambda = 555.555555555556;
    REAL mu = 833.333333333333;
    
    De(0,0) = lambda + 2.0*mu;
    De(1,1) = mu;
    De(2,2) = lambda + 2.0*mu;
    De(0,2) = lambda;
    De(2,0) = lambda;
}

void TPZNumericalIntegrator::ComputeTangentMatrix(int64_t iel, TPZFMatrix<REAL> &K){
    
    int n_sigma_comps = 3;
    int el_npts = fBlockMatrix.Blocks().fRowSizes[iel]/n_sigma_comps;
    int el_dofs = fBlockMatrix.Blocks().fColSizes[iel];
    int first_el_ip = fBlockMatrix.Blocks().fRowFirstIndex[iel]/n_sigma_comps;
    
    K.Resize(el_dofs, el_dofs);
    K.Zero();

    int pos = fBlockMatrix.Blocks().fMatrixPosition[iel];
    TPZFMatrix<STATE> De(3,3);
    TPZFMatrix<STATE> Bip(n_sigma_comps,el_dofs,0.0);
    TPZFMatrix<STATE> DeBip;
    int c = 0;
    for (int ip = 0; ip < el_npts; ip++) {
        for (int i = 0; i < n_sigma_comps; i++) {
            for (int j = 0; j < el_dofs; j++) {
                Bip(i,j) = fBlockMatrix.Blocks().fStorage[pos + c];
                c++;
            }
        }
        
        REAL omega = fConstitutiveLawProcessor.fWeight[first_el_ip + ip];
        ComputeConstitutiveMatrix(ip,De);
        De.Multiply(Bip, DeBip);
        Bip.MultAdd(DeBip, K, K, omega, 1.0, 1);
    }
}

void TPZNumericalIntegrator::ComputeTangentMatrix(int64_t iel, TPZFMatrix<REAL> &Dep, TPZFMatrix<REAL> &K){

    int n_sigma_comps = 3;
    int el_npts = fBlockMatrix.Blocks().fRowSizes[iel]/n_sigma_comps;
    int el_dofs = fBlockMatrix.Blocks().fColSizes[iel];
    int first_el_ip = fBlockMatrix.Blocks().fRowFirstIndex[iel]/n_sigma_comps;

    K.Resize(el_dofs, el_dofs);
    K.Zero();

    int pos = fBlockMatrix.Blocks().fMatrixPosition[iel];
    TPZFMatrix<STATE> dep(3,3);
    TPZFMatrix<STATE> Bip(n_sigma_comps,el_dofs,0.0);
    TPZFMatrix<STATE> DeBip;
    int c1 = 0;
    int c2 = 0;
    for (int ip = 0; ip < el_npts; ip++) {
        for (int i = 0; i < n_sigma_comps; i++) {
            for (int j = 0; j < el_dofs; j++) {
                Bip(i,j) = fBlockMatrix.Blocks().fStorage[pos + c1];
                c1++;
            }
        }

        for (int i = 0; i < n_sigma_comps; i++) {
            for (int j = 0; j < n_sigma_comps; j++) {
                dep(i,j) = Dep(first_el_ip * n_sigma_comps * n_sigma_comps + c2, 0);
                c2++;
            }
        }

//        dep.Print(std::cout);
        REAL omega = fConstitutiveLawProcessor.fWeight[first_el_ip + ip];
        dep.Multiply(Bip, DeBip);
        Bip.MultAdd(DeBip, K, K, omega, 1.0, 1);
    }
}
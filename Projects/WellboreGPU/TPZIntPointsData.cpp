//
// Created by natalia on 22/05/19.
//

#include "TPZIntPointsData.h"
#ifdef USING_MKL
#include <mkl.h>
#endif

TPZIntPointsData::TPZIntPointsData() : fNColor(-1), fIndexes(0), fIndexesColor(0), fWeight(0) {

}

TPZIntPointsData::~TPZIntPointsData() {

}

TPZIntPointsData::TPZIntPointsData(const TPZIntPointsData &copy) {
    fNColor = copy.fNColor;
    fIndexes = copy.fIndexes;
    fIndexesColor = copy.fIndexesColor;
    fWeight = copy.fWeight;
}

TPZIntPointsData &TPZIntPointsData::operator=(const TPZIntPointsData &copy) {
    if(&copy == this){
        return *this;
    }

    fNColor = copy.fNColor;
    fIndexes = copy.fIndexes;
    fIndexesColor = copy.fIndexesColor;
    fWeight = copy.fWeight;

    return *this;
}

void TPZIntPointsData::GatherSolution(TPZFMatrix<REAL> &solution, TPZFMatrix<REAL> &gather_solution) {
    int64_t nsize = gather_solution.Rows();
    cblas_dgthr(nsize, solution, &gather_solution(0, 0), &fIndexes[0]);
}

void TPZIntPointsData::ColoredAssemble(TPZFMatrix<REAL>  &nodal_forces, TPZFMatrix<REAL>  &residual) {
    int64_t ncolor = fNColor;
    int64_t nsize = nodal_forces.Rows();
    int64_t neq = residual.Rows();

    residual.Resize(neq*ncolor,1);
    residual.Zero();
    cblas_dsctr(nsize, nodal_forces, &fIndexesColor[0], &residual(0,0));

    int64_t colorassemb = ncolor / 2.;
    while (colorassemb > 0) {

        int64_t firsteq = (ncolor - colorassemb) * neq;
        cblas_daxpy(colorassemb * neq, 1., &residual(firsteq, 0), 1., &residual(0, 0), 1.);

        ncolor -= colorassemb;
        colorassemb = ncolor/2;
    }
    residual.Resize(neq, 1);
}
#ifndef TPZSHAPEHDIVCOLLAPSED_H
#define TPZSHAPEHDIVCOLLAPSED_H

#include "TPZShapeHDiv.h"

/// Class that implements the computation of H1 shape functions with variable connect order
template <class TSHAPE>
struct TPZShapeHDivCollapsed : public TPZShapeHDiv<TSHAPE>
{
    
    static void Initialize(const TPZVec<int64_t> &ids,
                           const TPZVec<int> &connectorders,
                           const TPZVec<int> &sideorient,
                           TPZShapeData &data);
    
    static int NShapeF(TPZShapeData &data);
    
    static void Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi);

};

#endif

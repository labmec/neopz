#ifndef TPZSHAPEHDIVBOUND_H
#define TPZSHAPEHDIVBOUND_H

#include "TPZShapeH1.h"
/// Class that implements the computation of H1 shape functions with variable connect order
template <class TSHAPE>
struct TPZShapeHDivBound
{
    
    static void Initialize(const TPZVec<int64_t> &ids, int connectorder,
                           int sideorient, TPZShapeData &data);
    
    static int NShape(const TPZShapeData &data);
    
    static void Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi);

};

#endif

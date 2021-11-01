#ifndef TPZSHAPEH1_H
#define TPZSHAPEH1_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
#include "TPZShapeData.h"

/// Class that implements the computation of H1 shape functions with variable connect order
template <class TSHAPE>
struct TPZShapeH1
{
    
    static void Initialize(const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZShapeData &data);
    
    static void Shape(const TPZVec<REAL> &pt, TPZShapeData &data);

};

#endif

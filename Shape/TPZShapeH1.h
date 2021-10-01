#ifndef TPZSHAPEH1_H
#define TPZSHAPEH1_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
#include "TPZShapeData.h"

template <class TSHAPE>
struct TPZShapeH1
{
//    TPZManVector<int,27> fConnectOrders;
//    TPZManVector<int,27> fNSideShape;
//    TPZVec<TPZTransform<REAL> > fTransformationVector;
    
    static void Initialize(const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZShapeData &data);
    
    static void Shape(TPZVec<REAL> &pt, TPZShapeData &data);

//    typedef typename TSHAPE::GraphElType GraphElType;
};

#endif

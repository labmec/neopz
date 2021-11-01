#ifndef TPZSHAPEHDIV_H
#define TPZSHAPEHDIV_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
#include "pzeltype.h"

template<class T>
class TPZFMatrix;

class TPZShapeData;

/// Traditional HDiv spaces, data structures that do not depend on the geometric map
template <class TSHAPE>
struct TPZShapeHDiv
{
    
    TPZShapeHDiv();
    
    static void Initialize(const TPZVec<int64_t> &ids,
                    const TPZVec<int> &connectorders,
                    const TPZVec<int>& sideorient, TPZShapeData &data);
    
    static void ComputeMasterDirections(TPZShapeData &data);
    
    static void ComputeVecandShape(TPZShapeData &data);
    
    static void IndexShapeToVec(TPZShapeData &data);
    
    static void FirstShapeIndex(TPZVec<int64_t> &Index, int &scalarorders);
    
    static void FillOrderScalarShapeFunctions(const TPZVec<int> &connectorders, TPZVec<int> &scalarOrder);
    
    static void HDivPermutation(int side, TPZShapeData &data, TPZVec<int> &permutegather);
    
    static void HDivPermutation(MElementType eltype, const TPZVec<int64_t> &ids, TPZVec<int> &permutegather);
    
    static void Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi);

    static int ComputeNConnectShapeF(int connect, int order);
    
    static int NConnectShapeF(int connect, const TPZShapeData &data);
    
    static int NShapeF(const TPZShapeData &shapedata);
};

#endif

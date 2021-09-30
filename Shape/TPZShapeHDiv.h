#ifndef TPZSHAPEHDIV_H
#define TPZSHAPEHDIV_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
template<class T>
class TPZFMatrix;

class TPZShapeData;

/// Traditional HDiv spaces, data structures that do not depend on the geometric map
template <class TSHAPE>
struct TPZShapeHDiv
{
    
    TPZShapeHDiv();
    
    static void Initialize(TPZVec<int64_t> &ids,
                    TPZVec<int> &connectorders,                    
                    TPZVec<int>& sideorient, TPZShapeData &data);
    
    static void ComputeMasterDirections(TPZShapeData &data);
    
    static void ComputeVecandShape(TPZShapeData &data);
    
    static void IndexShapeToVec(TPZShapeData &data);
    
    static void FirstShapeIndex(TPZVec<int64_t> &Index, TPZVec<int> &scalarorders);
    
    static void FillOrderScalarShapeFunctions(const TPZVec<int> &connectorders, TPZVec<int> &scalarOrder);
    
    static void HDivPermutation(int side, TPZShapeData &data, TPZVec<int> &permutegather);
    
    static void Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

};

#endif

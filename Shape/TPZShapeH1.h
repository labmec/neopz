#ifndef TPZSHAPEH1_H
#define TPZSHAPEH1_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
#include "TPZShapeData.h"
#include "pzshtmat.h"

/// Class that implements the computation of H1 shape functions with variable connect order
template <class TSHAPE>
struct TPZShapeH1
{
    
    static void Initialize(const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZShapeData &data);
    
    template<class T>
    static void Shape(const TPZVec<T> &pt, TPZShapeData &data, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi);

    static void Shape(const TPZVec<REAL> &pt, TPZShapeData &data)
    {
        TPZShapeH1<TSHAPE>::Shape(pt,data,data.fH1.fPhi,data.fH1.fDPhi);
    }

    static void ShapeOrders(TPZGenMatrix<int> &shapeorders, TPZShapeData &data);
    
};

template <class TSHAPE>
struct TPZSideShapeH1
{
    
    int fSide = -1;
    
    TPZSideShapeH1(int side) {
        fSide = side;
    }
    void Initialize(const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZShapeData &data);
    
    template<class T>
    static void Shape(const TPZVec<T> &pt, TPZShapeData &data, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi);

    void Shape(const TPZVec<REAL> &pt, TPZShapeData &data);

};


#endif

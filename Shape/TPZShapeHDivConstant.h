#ifndef TPZSHAPEHDIVCONSTANT_H
#define TPZSHAPEHDIVCONSTANT_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
template<class T>
class TPZFMatrix;

#include "TPZShapeData.h"
#include "TPZShapeH1.h"

/// Traditional HDiv spaces, data structures that do not depend on the geometric map
template <class TSHAPE>
struct TPZShapeHDivConstant : public TPZShapeH1<TSHAPE>
{
    
    static int NConnectShapeF(int connect, TPZShapeData &data);
    
    static int NHDivShapeF(TPZShapeData &data);

    static void Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi);
    
    static int ComputeNConnectShapeF(int connect, int order);
};

#endif

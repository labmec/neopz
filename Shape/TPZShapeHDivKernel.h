#ifndef TPZSHAPEHDIVKERNEL_H
#define TPZSHAPEHDIVKERNEL_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
template<class T>
class TPZFMatrix;

#include "TPZShapeData.h"
#include "TPZShapeHCurl.h"

/// Traditional HDiv spaces, data structures that do not depend on the geometric map
template <class TSHAPE>
struct TPZShapeHDivKernel : public TPZShapeHCurl<TSHAPE>
{
    
    
    static void ComputeVecandShape(TPZShapeData &data);
    
    static int NConnectShapeF(int connect, TPZShapeData &data);
    
    static int NHDivShapeF(TPZShapeData &data);

    static void Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi);

};

#endif

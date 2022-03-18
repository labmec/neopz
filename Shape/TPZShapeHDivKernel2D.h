#ifndef TPZSHAPEHDIVKERNEL2D_H
#define TPZSHAPEHDIVKERNEL2D_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
template<class T>
class TPZFMatrix;

#include "TPZShapeData.h"
#include "TPZShapeH1.h"

/// Traditional HDiv spaces, data structures that do not depend on the geometric map
template <class TSHAPE>
struct TPZShapeHDivKernel2D : public TPZShapeH1<TSHAPE>
{
    
    static int NHDivShapeF(TPZShapeData &data);

    static void Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi);

    
};

#endif

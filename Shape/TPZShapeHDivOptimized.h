/**
 * @file
 * @brief Contains TPZShapeHDivOptimized class which implements HDiv shape using HCurlNoGrads functions on the facet and HDiv standard interior functions
 * @author Giovane Avancini
 */

#pragma once

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
#include "pzeltype.h"
#include "TPZShapeHCurlNoGrads.h"
#include "TPZShapeHDiv.h"

template <class T>
class TPZFMatrix;

class TPZShapeData;

template <class TSHAPE>
struct TPZShapeHDivOptimized : public TPZShapeHDiv<TSHAPE>, TPZShapeHCurlNoGrads<TSHAPE>
{
    TPZShapeHDivOptimized();

    static void Initialize(const TPZVec<int64_t> &ids,
                           const TPZVec<int> &connectorders,
                           const TPZVec<int> &sideorient, TPZShapeData &data);

    static void Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi);

    static void Shape(const TPZVec<Fad<REAL>> &pt, TPZShapeData &data, TPZFMatrix<Fad<REAL>> &phi, TPZFMatrix<Fad<REAL>> &divphi);

    static int ComputeNConnectShapeF(int connect, int order);

    static int NConnectShapeF(int connect, const TPZShapeData &data);

    static int NShapeF(const TPZShapeData &shapedata);

    static void CheckH1ConnectOrder(const TPZVec<int> &connectorders, TPZVec<int> &H1Orders);
};

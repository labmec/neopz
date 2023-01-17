#ifndef TPZSHAPEH1_H
#define TPZSHAPEH1_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
#include "TPZShapeData.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzshapeprism.h"


/// Class that implements the computation of H1 shape functions with variable connect order
template <class TSHAPE>
struct TPZShapeH1
{
    
    static void Initialize(const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZShapeData &data);
    
    template<class T>
    static void Shape(const TPZVec<T> &pt, TPZShapeData &data, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi);

    template<class T>
    static void SideShape(int side, const TPZVec<T> &pt, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi);

    // compute the parameter transforms from the interior to the sides for all sides excluding the vertices
    static void ComputeTransforms(const TPZVec<int64_t> &id,  TPZVec<TPZTransform<REAL> > &transvec) ;

    /// compute the parameter transformation from the interior to a given side (of side dimension > 0)
    static TPZTransform<REAL> GetSideTransform(const int side, int trans_id);

    template<class T>
    static void ShapeInternal(int side, TPZVec<T> &x, int order, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
    {
        MElementType sidetype = TSHAPE::Type(side);
        switch (sidetype) {
            case EPoint:
                pzshape::TPZShapePoint::ShapeInternal(x, order, phi, dphi);
                break;
            case EOned:
                pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
                break;
            case EQuadrilateral:
                pzshape::TPZShapeQuad::ShapeInternal(x, order, phi, dphi);
                break;
            case ETriangle:
                pzshape::TPZShapeTriang::ShapeInternal(x, order, phi, dphi);
                break;
            case ECube:
                pzshape::TPZShapeCube::ShapeInternal(x, order, phi, dphi);
                break;
            case ETetraedro:
                pzshape::TPZShapeTetra::ShapeInternal(x, order, phi, dphi);
                break;
            case EPrisma:
                pzshape::TPZShapePrism::ShapeInternal(x, order, phi, dphi);
                break;
            case EPiramide:
                pzshape::TPZShapePiram::ShapeInternal(x, order, phi, dphi);
                break;
            default:
                DebugStop();
                break;
        }
    }
};

#endif

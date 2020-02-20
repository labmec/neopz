

#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapecube.h"
#include "pzshapepiram.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"

template <class TSHAPE>
TPZTransform<REAL> GetSideTransform(int side, int trans_id);
template <class TSHAPE>
void ComputeTransforms(TPZVec<int64_t> &id,  TPZVec<TPZTransform<REAL> > &transvec) ;
template <class TSHAPE>
void Shape(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) ;



template <class TSHAPE>
void ComputeTransforms(TPZVec<int64_t> &id,  TPZVec<TPZTransform<REAL> > &transvec) {
    int NSides = TSHAPE::NSides;
    int NCorners = TSHAPE::NCornerNodes;
    transvec.resize(NSides - NCorners);
    for (int iside = NCorners; iside< NSides ; iside++) {
        int pos = iside - NCorners;
        int trans_id = TSHAPE::GetTransformId(iside, id); // Foi criado
        transvec[pos] = GetSideTransform<TSHAPE>(iside, trans_id); // Foi criado
    }
}

template <class TSHAPE>
void Shape(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {

    int dim = TSHAPE::Dimension;
    int NSides = TSHAPE::NSides;
    int NCorners = TSHAPE::NCornerNodes;
    
    TPZFNMatrix<100> phiblend(NSides,1),dphiblend(dim,NSides);
    TSHAPE::ShapeCorner(pt,phi,dphi);
    
    phiblend = phi;
    dphiblend = dphi;
    
    TSHAPE::ShapeGenerating(pt,phiblend, dphiblend);
    int shape = NCorners;
    for (int side = NCorners; side<NSides ; side++)
    {
        TPZTransform<REAL> transform = transvec[side - NCorners];
        int sidedim = TSHAPE::SideDimension(side);
        int numshape =TSHAPE::NConnectShapeF( side, orders[side - NCorners]);
        TPZFNMatrix<100,REAL> phin(numshape,1), dphin(sidedim,numshape), dphiaux(sidedim,numshape);
        
        TPZManVector<REAL,1> outvec(sidedim);
        transform.Apply(pt, outvec);
        dphin.Zero();
        TSHAPE::ShapeInternal(side, outvec,orders[side - NCorners], phin, dphin);
        REAL alpha = 1.0, beta =0.0;
        int transpose = 1;
        transform.Mult().MultAdd(dphin, dphin, dphiaux, alpha, beta, transpose);
        for (int i = 0; i < numshape; i++) {
            phi(shape,0) = phiblend(side,0)*phin(i,0);
            for(int xj=0;xj<TSHAPE::Dimension;xj++) {
                dphi(xj,shape) = dphiblend(xj,side)*phin(i,0)+phiblend(side,0)*dphiaux(xj,i);
            }
            shape++;
        }
        
    }
    
}


template <class TSHAPE>
TPZTransform<REAL> GetSideTransform(int side, int trans_id) {
    
    MElementType type_side = TSHAPE::Type(side);
    TPZTransform<REAL> TransElToSide = TSHAPE::TransformElementToSide(side);
    
    TPZTransform<REAL> TransParametric(1,1);
    switch (type_side) {
        case EOned:
        {
            TransParametric = pzshape::TPZShapeLinear::ParametricTransform(trans_id);
        }
            break;
        case EQuadrilateral:
        {
            TransParametric = pzshape::TPZShapeQuad::ParametricTransform(trans_id);
        }
            break;
        case ETriangle:
        {
            TransParametric = pzshape::TPZShapeTriang::ParametricTransform(trans_id);
        }
            break;
        default:
            break;
    }
    
    if ((side == TSHAPE::NSides - 1) && TSHAPE::Dimension>2) {
        
        return TransElToSide;
    }
    else{
        TPZFMatrix<REAL> resul_mult;
        TransParametric.Mult().Multiply(TransElToSide.Mult(), resul_mult);
        TransElToSide.Mult() = resul_mult;
        TPZFMatrix<REAL> res1;
        TransParametric.Mult().Multiply( TransElToSide.Sum(), res1);
        TransElToSide.Sum() =res1;
        TransElToSide.Sum() += TransParametric.Sum();
        
    }
    
    return TransElToSide;
    
}


template void Shape<pzshape::TPZShapeLinear>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
template void Shape<pzshape::TPZShapeTriang>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
template void Shape<pzshape::TPZShapeQuad>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
template void Shape<pzshape::TPZShapePrism>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
template void Shape<pzshape::TPZShapePiram>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
template void Shape<pzshape::TPZShapeTetra>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
template void Shape<pzshape::TPZShapeCube>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);



//void Shape(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> phi,TPZFMatrix<REAL> dphi);

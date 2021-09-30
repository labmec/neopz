#include "pzgenericshape.h"

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"

template <class TSHAPE>
TPZTransform<REAL> GetSideTransform(const int side, int trans_id) {
    
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
    //why TSHAPE::Dimension>2 here??
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

template <class TSHAPE>
void ComputeTransforms(const TPZVec<int64_t> &id,  TPZVec<TPZTransform<REAL> > &transvec) {
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
void inline Shape(TPZVec<REAL> &pt, TParDefs &par, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi
) {

    
//    TPZVec<REAL> &pt = par.pt;
//    TPZFMatrix<REAL> &phi = par.phi;
//    TPZFMatrix<REAL> &dphi = par.dphi;
    
    TSHAPE::ShapeCorner(pt,phi,dphi);
    
    if(phi.Rows() == TSHAPE::NCornerNodes) return;
    
    const int dim = TSHAPE::Dimension;
    const int NSides = TSHAPE::NSides;
    const int NCorners = TSHAPE::NCornerNodes;

    TPZFNMatrix<27*3,REAL> phiblend(NSides,1),dphiblend(dim,NSides);
    for(int nod=0; nod<NCorners; nod++)
    {
        phiblend(nod,0) = phi(nod,0);
        for(int d=0; d< dim; d++)
        {
            dphiblend(d,nod) = dphi(d,nod);
        }
    }
    TSHAPE::ShapeGenerating(pt, par.nshape, phiblend, dphiblend);
    int shape = NCorners;
    TPZVec<int> &nshape = par.nshape;
    for (int side = NCorners; side<NSides ; side++)
    {
        int numshape =nshape[side - NCorners];
        if(numshape == 0) continue;
        
        phi(shape,0) = phiblend(side,0);
        for(int d=0; d<TSHAPE::Dimension; d++) dphi(d,shape) = dphiblend(d,side);
        shape++;
        
        if(numshape == 1) continue;
        
        TPZTransform<REAL> &transform = par.transvec[side - NCorners];
        int sidedim = TSHAPE::SideDimension(side);
        TPZFNMatrix<100,REAL> phin(numshape,1), dphin(sidedim,numshape), dphiaux(TSHAPE::Dimension,numshape),
            dphiaux2(TSHAPE::Dimension,numshape);
        TPZManVector<REAL,3> outvec(sidedim);
        transform.Apply(pt, outvec);
//        dphin.Zero();
        TSHAPE::ShapeInternal(side, outvec,par.orders[side - NCorners], phin, dphin);
        if(sidedim < 3)
        {
            constexpr REAL alpha = 1.;
            constexpr REAL beta = 0.;
            constexpr int opt = 1;
            TPZFNMatrix<1,REAL> auxmat(1,1,0);    
            TPZFMatrix<REAL> &mult = transform.Mult();
            mult.MultAdd(dphin, auxmat, dphiaux,alpha,beta,opt);
            
            for (int i = 1; i < numshape; i++) {
                phi(shape,0) = phiblend(side,0)*phin(i,0);
                for(int xj=0;xj<TSHAPE::Dimension;xj++) {
                    dphi(xj,shape) = dphiblend(xj,side)*phin(i,0)+phiblend(side,0)*dphiaux(xj,i);
                }
                shape++;
            }
        } else
        {
            for (int i = 1; i < numshape; i++) {
                phi(shape,0) = phiblend(side,0)*phin(i,0);
                for(int xj=0;xj<TSHAPE::Dimension;xj++) {
                    dphi(xj,shape) = dphiblend(xj,side)*phin(i,0)+phiblend(side,0)*dphin(xj,i);
                }
                shape++;
            }

        }
    }
    
}


#define IMPLEMENTSHAPE(TSHAPE)\
template \
void Shape<TSHAPE>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);\
template \
TPZTransform<REAL> GetSideTransform<TSHAPE>(const int side, int trans_id);\
template \
void ComputeTransforms<TSHAPE>(const TPZVec<int64_t> &id,  TPZVec<TPZTransform<REAL> > &transvec);

IMPLEMENTSHAPE(pzshape::TPZShapePoint)
IMPLEMENTSHAPE(pzshape::TPZShapeLinear)
IMPLEMENTSHAPE(pzshape::TPZShapeTriang)
IMPLEMENTSHAPE(pzshape::TPZShapeQuad)
IMPLEMENTSHAPE(pzshape::TPZShapeTetra)
IMPLEMENTSHAPE(pzshape::TPZShapeCube)
IMPLEMENTSHAPE(pzshape::TPZShapePrism)
IMPLEMENTSHAPE(pzshape::TPZShapePiram)

#undef IMPLEMENTSHAPE

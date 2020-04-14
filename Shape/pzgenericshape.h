

#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapecube.h"
#include "pzshapetriang.h"
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

struct TParDefs
{
    TPZManVector<int,27> orders;
    TPZManVector<int,27> nshape;
    TPZVec<TPZTransform<REAL> > transvec;
};

#include "TPZLapack.h"

template<int dim, int sidedim> void inline ExtendDerivative(TPZFMatrix<REAL> &dphin, TPZFMatrix<REAL> &dphivol, TPZFMatrix<REAL> &mult)
{
    REAL alpha = 1.;
    REAL beta = 0.;
    bool transpose = true;
    mult.MultAdd(dphin, dphin, dphivol, alpha, beta, transpose);
}

template<> void inline ExtendDerivative<1,1>(TPZFMatrix<REAL> &dphin, TPZFMatrix<REAL> &dphivol, TPZFMatrix<REAL> &mult)
{
    REAL val = mult(0,0);
    int64_t cols = dphin.Cols();
    for(int i = 0; i<cols; i++) dphivol(0,i) = val*dphin(0,i);
}

template<> void inline ExtendDerivative<2,1>(TPZFMatrix<REAL> &dphin, TPZFMatrix<REAL> &dphivol, TPZFMatrix<REAL> &mult)
{
    REAL val[2] = {mult(0,0),mult(0,1)};
    int64_t cols = dphin.Cols();
    for(int i = 0; i<cols; i++)
    {
        dphivol(0,i) = val[0]*dphin(0,i);
        dphivol(1,i) = val[1]*dphin(0,i);
    }
}

template<> void inline ExtendDerivative<2,2>(TPZFMatrix<REAL> &dphin, TPZFMatrix<REAL> &dphivol, TPZFMatrix<REAL> &mult)
{
    REAL val1[3] = {mult(0,0),mult(0,1)};
    REAL val2[3] = {mult(1,0),mult(1,1)};
    int64_t cols = dphin.Cols();
    for(int i = 0; i<cols; i++)
    {
        dphivol(0,i) = val1[0]*dphin(0,i)+val2[0]*dphin(1,i);
        dphivol(1,i) = val1[1]*dphin(0,i)+val2[1]*dphin(1,i);
    }
}


template<> void inline ExtendDerivative<3,1>(TPZFMatrix<REAL> &dphin, TPZFMatrix<REAL> &dphivol, TPZFMatrix<REAL> &mult)
{
    REAL val[3] = {mult(0,0),mult(0,1),mult(0,2)};
    int64_t cols = dphin.Cols();
    for(int i = 0; i<cols; i++)
    {
        dphivol(0,i) = val[0]*dphin(0,i);
        dphivol(1,i) = val[1]*dphin(0,i);
        dphivol(2,i) = val[2]*dphin(0,i);
    }
}

template<> void inline ExtendDerivative<3,2>(TPZFMatrix<REAL> &dphin, TPZFMatrix<REAL> &dphivol, TPZFMatrix<REAL> &mult)
{
    REAL val1[3] = {mult(0,0),mult(0,1),mult(0,2)};
    REAL val2[3] = {mult(1,0),mult(1,1),mult(1,2)};
    int64_t cols = dphin.Cols();
    for(int i = 0; i<cols; i++)
    {
        dphivol(0,i) = val1[0]*dphin(0,i)+val2[0]*dphin(1,i);
        dphivol(1,i) = val1[1]*dphin(0,i)+val2[1]*dphin(1,i);
        dphivol(2,i) = val1[2]*dphin(0,i)+val2[2]*dphin(1,i);
    }
}

template<int dim, int sidedim> void inline Apply(TPZVec<REAL> &pt, TPZVec<REAL> &project, TPZFMatrix<REAL> &mult)
{
    for(int is = 0; is<sidedim; is++) for(int d = 0; d<dim; d++)
        project[is] += mult(is,d)*pt[d];
}

template<> void inline Apply<1,1>(TPZVec<REAL> &pt, TPZVec<REAL> &project, TPZFMatrix<REAL> &mult)
{
    project[0] += mult(0,0)*pt[0];
}

template<> void inline Apply<2,1>(TPZVec<REAL> &pt, TPZVec<REAL> &project, TPZFMatrix<REAL> &mult)
{
    project[0] += mult(0,0)*pt[0]+mult(0,1)*pt[1];
}

template<> void inline Apply<2,2>(TPZVec<REAL> &pt, TPZVec<REAL> &project, TPZFMatrix<REAL> &mult)
{
    project[0] += mult(0,0)*pt[0]+mult(0,1)*pt[1];
    project[1] += mult(1,0)*pt[0]+mult(1,1)*pt[1];
}

template<> void inline Apply<3,1>(TPZVec<REAL> &pt, TPZVec<REAL> &project, TPZFMatrix<REAL> &mult)
{
    project[0] += mult(0,0)*pt[0]+mult(0,1)*pt[1]+mult(0,2)*pt[2];
}
template<> void inline Apply<3,2>(TPZVec<REAL> &pt, TPZVec<REAL> &project, TPZFMatrix<REAL> &mult)
{
    project[0] += mult(0,0)*pt[0]+mult(0,1)*pt[1]+mult(0,2)*pt[2];
    project[1] += mult(1,0)*pt[0]+mult(1,1)*pt[1]+mult(1,2)*pt[2];
}

template <class TSHAPE>
void inline Shape(TPZVec<REAL> &pt, TParDefs &par, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi
) {

    
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
        switch(sidedim)
        {
            case 1:
                outvec[0] = transform.Sum()(0,0);
                Apply<TSHAPE::Dimension,1>(pt,outvec,transform.Mult());
                break;
            case 2:
                outvec[0] = transform.Sum()(0,0);
                outvec[1] = transform.Sum()(1,0);
                Apply<TSHAPE::Dimension,2>(pt,outvec,transform.Mult());
                break;
            case 3:
                outvec = pt;
                break;
            default:
                break;
        }
//        dphin.Zero();
        TSHAPE::ShapeInternal(side, outvec,par.orders[side - NCorners], phin, dphin);
        if(sidedim < 3)
        {
            switch(sidedim)
            {
                case 1:
                    ExtendDerivative<TSHAPE::Dimension,1>(dphin,dphiaux,transform.Mult());
                    break;
                case 2:
                    ExtendDerivative<TSHAPE::Dimension,2>(dphin,dphiaux,transform.Mult());
                    break;
                default:
                    DebugStop();
            }
            
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

template void Shape<pzshape::TPZShapeCube>(TPZVec<REAL> &pt, TParDefs &par, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);


//void Shape(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> phi,TPZFMatrix<REAL> dphi);

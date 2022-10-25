#include "TPZShapeH1.h"
//#include "pzgenericshape.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"


template <class TSHAPE>
void TPZShapeH1<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
                                           const TPZVec<int> &connectorders,
                                           TPZShapeData &data) {
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    if((ids.size() != ncorner || connectorders.size() != nsides-ncorner) && TSHAPE::Type() != EPoint)
    {
        DebugStop();
    }
    data.fH1ConnectOrders = connectorders;
    ComputeTransforms(ids, data.fSideTransforms);
    data.fCornerNodeIds = ids;
    data.fH1NumConnectShape.resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
    int64_t nshape = TSHAPE::NCornerNodes;
    for(int i = TSHAPE::NCornerNodes; i < TSHAPE::NSides; i++)
    {
        int nshapeconnect = TSHAPE::NConnectShapeF(i,connectorders[i-TSHAPE::NCornerNodes]);
        data.fH1NumConnectShape[i-TSHAPE::NCornerNodes] = nshapeconnect;
        nshape += nshapeconnect;
    }
    data.fPhi.Resize(nshape,1);
    data.fDPhi.Resize(TSHAPE::Dimension, nshape);
}


/*
 /// this is the "old" non-templated version. It is here for the case the templated version may have a big impact on performance
template <class TSHAPE>
void TPZShapeH1<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {

        
    TSHAPE::ShapeCorner(pt,data.fPhi,data.fDPhi);
    
    if(data.fPhi.Rows() == TSHAPE::NCornerNodes) return;
    
    const int dim = TSHAPE::Dimension;
    const int NSides = TSHAPE::NSides;
    const int NCorners = TSHAPE::NCornerNodes;

    TPZFNMatrix<NSides*dim,REAL> phiblend(NSides,1),dphiblend(dim,NSides);
    for(int nod=0; nod<NCorners; nod++)
    {
        phiblend(nod,0) = data.fPhi(nod,0);
        for(int d=0; d< dim; d++)
        {
            dphiblend(d,nod) = data.fDPhi(d,nod);
        }
    }
    TSHAPE::ShapeGenerating(pt, phiblend, dphiblend);
    int shape = NCorners;
    for (int side = NCorners; side<NSides ; side++)
    {
        int numshape = TSHAPE::NConnectShapeF(side, data.fH1ConnectOrders[side-NCorners]);
        if(numshape == 0) continue;
        
        data.fPhi(shape,0) = phiblend(side,0);
        for(int d=0; d<dim; d++) data.fDPhi(d,shape) = dphiblend(d,side);
        shape++;
        
        if(numshape == 1) continue;
        
        TPZTransform<REAL> &transform = data.fSideTransforms[side - NCorners];
        int sidedim = TSHAPE::SideDimension(side);
        TPZFNMatrix<100,REAL> phin(numshape,1), dphin(sidedim,numshape), dphiaux(TSHAPE::Dimension,numshape),
            dphiaux2(TSHAPE::Dimension,numshape);
        TPZManVector<REAL,3> outvec(sidedim);
        transform.Apply(pt, outvec);
//        dphin.Zero();
        TSHAPE::ShapeInternal(side, outvec,data.fH1ConnectOrders[side - NCorners], phin, dphin);
        if(sidedim < 3)
        {
            constexpr REAL alpha = 1.;
            constexpr REAL beta = 0.;
            constexpr int opt = 1;
            TPZFNMatrix<1,REAL> auxmat(1,1,0);    
            TPZFMatrix<REAL> &mult = transform.Mult();
            mult.MultAdd(dphin, auxmat, dphiaux,alpha,beta,opt);
            
            for (int i = 1; i < numshape; i++) {
                data.fPhi(shape,0) = phiblend(side,0)*phin(i,0);
                for(int xj=0;xj<TSHAPE::Dimension;xj++) {
                    data.fDPhi(xj,shape) = dphiblend(xj,side)*phin(i,0)+phiblend(side,0)*dphiaux(xj,i);
                }
                shape++;
            }
        } else
        {
            for (int i = 1; i < numshape; i++) {
                data.fPhi(shape,0) = phiblend(side,0)*phin(i,0);
                for(int xj=0;xj<TSHAPE::Dimension;xj++) {
                    data.fDPhi(xj,shape) = dphiblend(xj,side)*phin(i,0)+phiblend(side,0)*dphin(xj,i);
                }
                shape++;
            }

        }
    }
    phi = data.fPhi;
    dphi = data.fDPhi;
}
*/

template <class TSHAPE>
template <class T>
void TPZShapeH1<TSHAPE>::Shape(const TPZVec<T> &pt, TPZShapeData &data, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
{
    TSHAPE::ShapeCorner(pt,phi,dphi);
    
    if(data.fPhi.Rows() == TSHAPE::NCornerNodes) return;
    
    const int dim = TSHAPE::Dimension;
    
    if constexpr (std::is_same_v<FADREAL, T>)
    {
        if(pt[0].size() != dim) DebugStop();
    }
    const int NSides = TSHAPE::NSides;
    const int NCorners = TSHAPE::NCornerNodes;

    TPZFNMatrix<NSides*dim,T> phiblend(NSides,1),dphiblend(dim,NSides);
    for(int nod=0; nod<NCorners; nod++)
    {
        phiblend(nod,0) = phi(nod,0);
        for(int d=0; d< dim; d++)
        {
            dphiblend(d,nod) = dphi(d,nod);
        }
    }
    TSHAPE::ShapeGenerating(pt, phiblend, dphiblend);
    int shape = NCorners;
    for (int side = NCorners; side<NSides ; side++)
    {
        int numshape = TSHAPE::NConnectShapeF(side, data.fH1ConnectOrders[side-NCorners]);
        if(numshape == 0) continue;
        
        phi(shape,0) = phiblend(side,0);
        for(int d=0; d<dim; d++) dphi(d,shape) = dphiblend(d,side);
        shape++;
        
        if(numshape == 1) continue;
        
        TPZTransform<REAL> &transformREAL = data.fSideTransforms[side - NCorners];
        TPZTransform<T> transform(transformREAL.Rows(),transformREAL.Cols());
        transform.CopyFrom(transformREAL);
//        for (int i=0; i<transformREAL.Rows(); i++) {
//            transform.Sum()(i,0) = FADREAL(dim,transformREAL.Sum()(i,0));
//            for (int j= 0; j<transformREAL.Cols(); j++) {
//                transform.Mult()(i,j) = FADREAL(dim,transformREAL.Mult()(i,j));
//            }
//        }
        int sidedim = TSHAPE::SideDimension(side);
        TPZFNMatrix<100,T> phin(numshape,1), dphin(sidedim,numshape), dphiaux(TSHAPE::Dimension,numshape),
            dphiaux2(TSHAPE::Dimension,numshape);
        TPZManVector<T,3> outvec(sidedim);
        transform.Apply(pt, outvec);
//        dphin.Zero();
        ShapeInternal(side, outvec,data.fH1ConnectOrders[side - NCorners], phin, dphin);
        if(sidedim < 3)
        {
            T alpha = 1.;
            T beta = 0.;
            constexpr int opt = 1;
            TPZFNMatrix<3,T> auxmat(1,1,0);
            TPZFMatrix<T> &mult = transform.Mult();
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
    if constexpr (std::is_same_v<REAL, T>)
    {
        data.fPhi = phi;
        data.fDPhi = dphi;
    }
}

template <class TSHAPE>
TPZTransform<REAL> TPZShapeH1<TSHAPE>::GetSideTransform(const int side, int trans_id) {
    
    MElementType type_side = TSHAPE::Type(side);
    /// this is the generic transform mapping the interior point to the side
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
        case ECube:
        case ETetraedro:
        case EPrisma:
        case EPiramide:
            break;
        default:
#ifdef PZDEBUG
            DebugStop();
#endif
            break;
    }
    //when dimension == 3 there is no permutation, because 3D elements are not supposed to overlap
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
void TPZShapeH1<TSHAPE>::ComputeTransforms(const TPZVec<int64_t> &id,  TPZVec<TPZTransform<REAL> > &transvec) {
    int NSides = TSHAPE::NSides;
    int NCorners = TSHAPE::NCornerNodes;
    transvec.resize(NSides - NCorners);
    for (int iside = NCorners; iside< NSides ; iside++) {
        int pos = iside - NCorners;
        /// permutation id as a function of the ids of the corner nodes
        int trans_id = TSHAPE::GetTransformId(iside, id); // Foi criado
        transvec[pos] = GetSideTransform(iside, trans_id); // Foi criado
    }
}

template <class TSHAPE>
template <class T>
void TPZShapeH1<TSHAPE>::SideShape(int side, const TPZVec<T> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
{
    MElementType Type = TSHAPE::Type(side);
    TPZShapeData sidedata;
    switch (Type) {
        case EPoint:
            TPZShapeH1<pzshape::TPZShapePoint>::Initialize(ids, connectorders, sidedata);
            TPZShapeH1<pzshape::TPZShapePoint>::Shape(point,sidedata,phi,dphi);
            break;
        case EOned:
            TPZShapeH1<pzshape::TPZShapeLinear>::Initialize(ids, connectorders, sidedata);
            TPZShapeH1<pzshape::TPZShapeLinear>::Shape(point,sidedata,phi,dphi);
            break;
        case EQuadrilateral:
            TPZShapeH1<pzshape::TPZShapeQuad>::Initialize(ids, connectorders, sidedata);
            TPZShapeH1<pzshape::TPZShapeQuad>::Shape(point,sidedata,phi,dphi);
            break;
        case ETriangle:
            TPZShapeH1<pzshape::TPZShapeTriang>::Initialize(ids, connectorders, sidedata);
            TPZShapeH1<pzshape::TPZShapeTriang>::Shape(point,sidedata,phi,dphi);
            break;
        case ECube:
            TPZShapeH1<pzshape::TPZShapeCube>::Initialize(ids, connectorders, sidedata);
            TPZShapeH1<pzshape::TPZShapeCube>::Shape(point,sidedata,phi,dphi);
            break;
        case ETetraedro:
            TPZShapeH1<pzshape::TPZShapeTetra>::Initialize(ids, connectorders, sidedata);
            TPZShapeH1<pzshape::TPZShapeTetra>::Shape(point,sidedata,phi,dphi);
            break;
        case EPrisma:
            TPZShapeH1<pzshape::TPZShapePrism>::Initialize(ids, connectorders, sidedata);
            TPZShapeH1<pzshape::TPZShapePrism>::Shape(point,sidedata,phi,dphi);
            break;
        case EPiramide:
            TPZShapeH1<pzshape::TPZShapePiram>::Initialize(ids, connectorders, sidedata);
            TPZShapeH1<pzshape::TPZShapePiram>::Shape(point,sidedata,phi,dphi);
            break;
        default:
            break;
    }
}


template
struct TPZShapeH1<pzshape::TPZShapePoint>;

template
struct TPZShapeH1<pzshape::TPZShapeLinear>;

template
struct TPZShapeH1<pzshape::TPZShapeTriang>;

template
struct TPZShapeH1<pzshape::TPZShapeQuad>;

template
struct TPZShapeH1<pzshape::TPZShapeTetra>;

template
struct TPZShapeH1<pzshape::TPZShapeCube>;

template
struct TPZShapeH1<pzshape::TPZShapePrism>;

template
struct TPZShapeH1<pzshape::TPZShapePiram>;

template
void TPZShapeH1<pzshape::TPZShapePoint>::Shape(const TPZVec<FADREAL> &pt, TPZShapeData &data, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeLinear>::Shape(const TPZVec<FADREAL> &pt, TPZShapeData &data, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeTriang>::Shape(const TPZVec<FADREAL> &pt, TPZShapeData &data, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeQuad>::Shape(const TPZVec<FADREAL> &pt, TPZShapeData &data, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeCube>::Shape(const TPZVec<FADREAL> &pt, TPZShapeData &data, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapePrism>::Shape(const TPZVec<FADREAL> &pt, TPZShapeData &data, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeTetra>::Shape(const TPZVec<FADREAL> &pt, TPZShapeData &data, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapePiram>::Shape(const TPZVec<FADREAL> &pt, TPZShapeData &data, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);

template
void TPZShapeH1<pzshape::TPZShapePoint>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeLinear>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeTriang>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeQuad>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeCube>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapePrism>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeTetra>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapePiram>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

template
void TPZShapeH1<pzshape::TPZShapePoint>::SideShape(int side, const TPZVec<REAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeLinear>::SideShape(int side, const TPZVec<REAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeQuad>::SideShape(int side, const TPZVec<REAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeTriang>::SideShape(int side, const TPZVec<REAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeCube>::SideShape(int side, const TPZVec<REAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeTetra>::SideShape(int side, const TPZVec<REAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapePrism>::SideShape(int side, const TPZVec<REAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapePiram>::SideShape(int side, const TPZVec<REAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

template
void TPZShapeH1<pzshape::TPZShapePoint>::SideShape(int side, const TPZVec<FADREAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeLinear>::SideShape(int side, const TPZVec<FADREAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeQuad>::SideShape(int side, const TPZVec<FADREAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeTriang>::SideShape(int side, const TPZVec<FADREAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeCube>::SideShape(int side, const TPZVec<FADREAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapeTetra>::SideShape(int side, const TPZVec<FADREAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapePrism>::SideShape(int side, const TPZVec<FADREAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);
template
void TPZShapeH1<pzshape::TPZShapePiram>::SideShape(int side, const TPZVec<FADREAL> &point, const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);

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
}

using namespace pzshape;

template <class TSHAPE>
void TPZShapeH1<TSHAPE>::ShapeOrders(TPZGenMatrix<int> &shapeorders, TPZShapeData &data)
{
    int nshape = data.fPhi.Rows();
    if(shapeorders.Rows() != nshape || shapeorders.Cols() != 3) DebugStop();
    for(int i = 0; i < TSHAPE::NCornerNodes; i++) {
        shapeorders(i,0) = 1;
        shapeorders(i,1) = 0;
        shapeorders(i,2) = 0;
    }
    if(TSHAPE::Dimension < 1) {
        if(TSHAPE::NCornerNodes != nshape) DebugStop();
        return;
    }
    int count = TSHAPE::NCornerNodes;
    int numoned = TSHAPE::NumSides(1);
    for(int side = TSHAPE::NCornerNodes; side < TSHAPE::NCornerNodes+numoned; side++) {
        int nshape = data.fH1NumConnectShape[side-TSHAPE::NCornerNodes];
        TPZGenMatrix<int> locorder(nshape,3);
        const int nsidecorners = TSHAPE::NContainedSides(side);
        TPZManVector<int64_t,8> ids(nsidecorners,0);
        for(int i=0; i<nsidecorners; i++) {
            int locid = TSHAPE::ContainedSideLocId(side,i);
        }
        int order = data.fH1ConnectOrders[side-TSHAPE::NCornerNodes];
        TPZShapeLinear::InternalShapeOrder(ids, order, locorder);
        for(int ish=0; ish<nshape; ish++) {
            for(int i=0; i<3; i++) {
                shapeorders(count,i) = locorder(ish,i);
            }
            count++;
        }
    }
    if(TSHAPE::Dimension < 2) {
        if(count != nshape) DebugStop();
        return;
    }
    int numtwod = TSHAPE::NumSides(2);
    for(int side = TSHAPE::NCornerNodes+numoned; side < TSHAPE::NCornerNodes+numoned+numtwod; side++) {
        int nshape = data.fH1NumConnectShape[side-TSHAPE::NCornerNodes];
        TPZGenMatrix<int> locorder(nshape,3);
        const int nsidecorners = TSHAPE::NContainedSides(side);
        TPZManVector<int64_t,8> ids(nsidecorners,0);
        for(int i=0; i<nsidecorners; i++) {
            int locid = TSHAPE::ContainedSideLocId(side,i);
        }
        int order = data.fH1ConnectOrders[side-TSHAPE::NCornerNodes];
        if(TSHAPE::Type(side) == EQuadrilateral) {
            TPZShapeQuad::InternalShapeOrder(ids, order, locorder);
        } else if(TSHAPE::Type(side) == ETriangle) {
            TPZShapeTriang::InternalShapeOrder(ids, order, locorder);
        } else {
            DebugStop();
        }
        for(int ish=0; ish<nshape; ish++) {
            for(int i=0; i<3; i++) {
                shapeorders(count,i) = locorder(ish,i);
            }
            count++;
        }
    }
    if(TSHAPE::Dimension == 3) {
        int side = TSHAPE::NSides - 1;
        int nshape = data.fH1NumConnectShape[side-TSHAPE::NCornerNodes];
        TPZGenMatrix<int> locorder(nshape,3);
        const int nsidecorners = TSHAPE::NContainedSides(side);
        TPZManVector<int64_t,8> ids(nsidecorners,0);
        for(int i=0; i<nsidecorners; i++) {
            int locid = TSHAPE::ContainedSideLocId(side,i);
        }
        int order = data.fH1ConnectOrders[side-TSHAPE::NCornerNodes];
        TSHAPE::InternalShapeOrder(ids, order, locorder);
        for(int ish=0; ish<nshape; ish++) {
            for(int i=0; i<3; i++) {
                shapeorders(count,i) = locorder(ish,i);
            }
            count++;
        }
    }
    if(count != nshape) DebugStop();
}


using namespace pzshape;

template<class TSHAPE>
static void Init(const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZShapeData &data)
{
    TPZShapeH1<TSHAPE>::Initialize(ids, connectorders, data);
}

template <class TSHAPE>
void TPZSideShapeH1<TSHAPE>::Initialize(const TPZVec<int64_t> &ids, const TPZVec<int> &connectorders, TPZShapeData &data) {
    if(fSide == TSHAPE::NSides-1) {
        TPZShapeH1<TSHAPE>::Initialize(ids, connectorders, data);
    } else {
        switch (TSHAPE::Type(fSide)) {
            case EPoint:
                Init<TPZShapePoint>(ids, connectorders, data);
                break;
            case EOned:
                Init<TPZShapeLinear>(ids, connectorders, data);
                break;
            case ETriangle:
                Init<TPZShapeTriang>(ids, connectorders, data);
                break;
            case EQuadrilateral:
                Init<TPZShapeQuad>(ids, connectorders, data);
                break;
            default:
                DebugStop();
                break;
        }
    }
    
        
    
}
    
template <class TSHAPE>
void TPZSideShapeH1<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data) {
    if(fSide == TSHAPE::NSides-1) {
        TPZShapeH1<TSHAPE>::Shape(pt, data);
    } else {
        switch (TSHAPE::Type(fSide)) {
            case EPoint:
                TPZShapeH1<TPZShapePoint>::Shape(pt, data);
                break;
            case EOned:
                TPZShapeH1<TPZShapeLinear>::Shape(pt, data);
                break;
            case ETriangle:
                TPZShapeH1<TPZShapeTriang>::Shape(pt, data);
                break;
            case EQuadrilateral:
                TPZShapeH1<TPZShapeQuad>::Shape(pt, data);
                break;
            default:
                DebugStop();
                break;
        }
    }
    
        
    

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
struct TPZSideShapeH1<pzshape::TPZShapePoint>;

template
struct TPZSideShapeH1<pzshape::TPZShapeLinear>;

template
struct TPZSideShapeH1<pzshape::TPZShapeTriang>;

template
struct TPZSideShapeH1<pzshape::TPZShapeQuad>;

template
struct TPZSideShapeH1<pzshape::TPZShapeTetra>;

template
struct TPZSideShapeH1<pzshape::TPZShapeCube>;

template
struct TPZSideShapeH1<pzshape::TPZShapePrism>;

template
struct TPZSideShapeH1<pzshape::TPZShapePiram>;

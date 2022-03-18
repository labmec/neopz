#include "TPZShapeH1.h"
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
    ComputeTransforms<TSHAPE>(ids, data.fSideTransforms);
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
//    fNodeIds = ids; Should we make it an attribute of the class?
}



template <class TSHAPE>
void TPZShapeH1<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data) {

    
//    TPZVec<REAL> &pt = par.pt;
//    TPZFMatrix<REAL> &phi = par.phi;
//    TPZFMatrix<REAL> &dphi = par.dphi;
    
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

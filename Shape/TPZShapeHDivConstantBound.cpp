#include "TPZShapeHDivConstantBound.h"
#include "TPZShapeHDiv.h"
#include "pzgenericshape.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"


template <class TSHAPE>
void TPZShapeHDivConstantBound<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
                                           int connectorder, int sideorient,
                                           TPZShapeData &data) {
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
#ifdef PZDEBUG
    if(ids.size() != ncorner || connectorder < 0 || std::abs(sideorient) != 1)
    {
        DebugStop();
    }
#endif
    data.fSideOrient.Resize(1);
    data.fSideOrient[0] = sideorient;
    data.fH1ConnectOrders.Resize(nsides-ncorner, connectorder);
    data.fH1ConnectOrders.Fill(connectorder);
    if(connectorder > 0)
    {
        TPZShapeH1<TSHAPE>::Initialize(ids,data.fH1ConnectOrders,data);
    }
}

template <class TSHAPE>
int TPZShapeHDivConstantBound<TSHAPE>::NShape(const TPZShapeData &data)
{
    if(data.fH1ConnectOrders[0] == 0) return 1;
    return data.fPhi.Rows();
}

template <>
void TPZShapeHDivConstantBound<class pzshape::TPZShapePoint>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi) {
    phi(0,0) = data.fSideOrient[0];
}


template <class TSHAPE>
void TPZShapeHDivConstantBound<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi) {

    
//    TPZVec<REAL> &pt = par.pt;
//    TPZFMatrix<REAL> &phi = par.phi;
//    TPZFMatrix<REAL> &dphi = par.dphi;
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const REAL volume = TSHAPE::RefElVolume();
    if(data.fH1ConnectOrders[0] == 0)
    {
        phi(0,0) = 1./volume;
        return;
    }
    const REAL scale = data.fSideOrient[0]*ncorner/volume;
    TPZShapeH1<TSHAPE>::Shape(pt,data);
    TPZManVector<int, 9> permutegather(nsides);
    TPZShapeHDiv<TSHAPE>::HDivPermutation(TSHAPE::Type(), data.fCornerNodeIds, permutegather);
    // for(int c=0; c < ncorner; c++) phi(c,0) = data.fPhi(permutegather[c],0)*scale;
    int count = 0;
    TPZManVector<int,9> firstshape(nsides-ncorner+1,ncorner);
    // phi(0,0) = data.fDPhi*data.axes;
    for(int side = ncorner; side < nsides; side++)
    {
        int nshape = data.fH1NumConnectShape[side-ncorner];
        firstshape[side+1-ncorner] = firstshape[side-ncorner]+ nshape;
    }
    for(int side = ncorner; side < nsides; side++)
    {
        int sidebound = permutegather[side];
        int nshape = data.fH1NumConnectShape[sidebound-ncorner];
        int fsh = firstshape[sidebound-ncorner];
        for(int sh = 0; sh<nshape; sh++)
        {
            phi(count,0) = data.fPhi(fsh+sh,0)*scale;
            count++;
        }
    }
}

template<class TSHAPE>
int TPZShapeHDivConstantBound<TSHAPE>::ComputeNConnectShapeF(int connect, int order)
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif
    // int order = data.fHDivConnectOrders[connect];
    MElementType thistype = TSHAPE::Type();

    if(thistype == EOned)
    {
        if(connect < 2) return order;
        else return order;
        // DebugStop();
    }
    else if(thistype == ETriangle)
    {
        if(connect < TSHAPE::NFacets) return 1;//(order+1);
        else return 1;//(order+1)*(order+1)-1;
    }
    else if(thistype == EQuadrilateral)
    {
        if(connect < TSHAPE::NFacets) return 0;//(order+1);
        else return 2*order*(order+1);
    }

    DebugStop();
    unreachable();
 }


template
struct TPZShapeHDivConstantBound<pzshape::TPZShapePoint>;

template
struct TPZShapeHDivConstantBound<pzshape::TPZShapeLinear>;

template
struct TPZShapeHDivConstantBound<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivConstantBound<pzshape::TPZShapeQuad>;


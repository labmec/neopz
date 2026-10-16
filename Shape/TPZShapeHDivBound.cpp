#include "TPZShapeHDivBound.h"
#include "TPZShapeHDiv.h"
#include "pzgenericshape.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"


template <class TSHAPE>
void TPZShapeHDivBound<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
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
    data.fHDiv.fSideOrient.Resize(1);
    data.fHDiv.fSideOrient[0] = sideorient;
    int connectOrdersSize = nsides-ncorner;
    if (TSHAPE::Type() == EPoint){
        connectOrdersSize = 1;
    }
    data.fH1.fConnectOrders.Resize(connectOrdersSize, connectorder);
    data.fH1.fConnectOrders.Fill(connectorder);
    if(connectorder > 0)
    {
        TPZShapeH1<TSHAPE>::Initialize(ids,data.fH1.fConnectOrders,data);
    }
}

template <class TSHAPE>
int TPZShapeHDivBound<TSHAPE>::NShape(const TPZShapeData &data)
{
    if(data.fH1.fConnectOrders[0] == 0) return 1;
    return data.fH1.fPhi.Rows();
}

template <>
void TPZShapeHDivBound<class pzshape::TPZShapePoint>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi) {
    phi(0,0) = data.fHDiv.fSideOrient[0];
}


template <class TSHAPE>
void TPZShapeHDivBound<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi) {

    
//    TPZVec<REAL> &pt = par.pt;
//    TPZFMatrix<REAL> &phi = par.phi;
//    TPZFMatrix<REAL> &dphi = par.dphi;
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const REAL volume = TSHAPE::RefElVolume();
    if(data.fH1.fConnectOrders[0] == 0)
    {
        phi(0,0) = 1./volume;
        return;
    }
    const REAL scale = data.fHDiv.fSideOrient[0]*ncorner/volume;
    TPZShapeH1<TSHAPE>::Shape(pt,data,data.fH1.fPhi,data.fH1.fDPhi);
    TPZManVector<int, 9> permutegather(nsides);
    TPZShapeHDiv<TSHAPE>::HDivPermutation(TSHAPE::Type(), data.fCornerNodeIds, permutegather);
    for(int c=0; c < ncorner; c++) phi(c,0) = data.fH1.fPhi(permutegather[c],0)*scale;
    int count = ncorner;
    TPZManVector<int,9> firstshape(nsides-ncorner+1,ncorner);
    for(int side = ncorner; side < nsides; side++)
    {
        int nshape = data.fH1.fNumConnectShape[side-ncorner];
        firstshape[side+1-ncorner] = firstshape[side-ncorner]+ nshape;
    }
    for(int side = ncorner; side < nsides; side++)
    {
        int sidebound = permutegather[side];
        int nshape = data.fH1.fNumConnectShape[sidebound-ncorner];
        int fsh = firstshape[sidebound-ncorner];
        for(int sh = 0; sh<nshape; sh++)
        {
            phi(count,0) = data.fH1.fPhi(fsh+sh,0)*scale;
            count++;
        }
    }
}

template
struct TPZShapeHDivBound<pzshape::TPZShapePoint>;

template
struct TPZShapeHDivBound<pzshape::TPZShapeLinear>;

template
struct TPZShapeHDivBound<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivBound<pzshape::TPZShapeQuad>;


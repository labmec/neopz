#include "TPZShapeHDivCollapsed.h"

#include "TPZShapeData.h"
#include "TPZShapeH1.h"
#include "TPZShapeHDivBound.h"
#include "pzgenericshape.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"


template <class TSHAPE>
void TPZShapeHDivCollapsed<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
                                    const TPZVec<int> &connectorders,
                                    const TPZVec<int> &sideorient,
                                        TPZShapeData &data) {
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    if(ids.size() != ncorner || connectorders.size() != TSHAPE::NFacets+3)
    {
        DebugStop();
    }
    TPZManVector<int,20> hdivconnectorders(connectorders);
    hdivconnectorders.Resize(TSHAPE::NFacets+1, 0);
    TPZManVector<int,20> hdivsideorient(sideorient);
    hdivsideorient.Resize(TSHAPE::NFacets+1, 0);
    TPZShapeHDiv<TSHAPE>::Initialize(ids, hdivconnectorders, hdivsideorient, data);
    int nphi = data.fH1.fPhi.Rows();
    // extend the dimension of the master element directions
    const int ndirections = data.fHDiv.fMasterDirections.Cols();
    data.fHDiv.fMasterDirections.Resize(TSHAPE::Dimension+1, ndirections+2);
    for (int idir = 0; idir < ndirections+2; idir++) {
        data.fHDiv.fMasterDirections(dim,idir) = 0.;
    }
    for (int d = 0; d<dim+1; d++) {
        data.fHDiv.fMasterDirections(d,ndirections) = 0.;
        data.fHDiv.fMasterDirections(d,ndirections+1) = 0.;
    }
    data.fHDiv.fMasterDirections(dim,ndirections) = 1.;
    data.fHDiv.fMasterDirections(dim,ndirections+1) = -1.;
    // extend the dimensions of the connect orders
    data.fHDiv.fConnectOrders = connectorders;
    data.fHDiv.fSideOrient = sideorient;
    data.fHDiv.fNumConnectShape.Resize(TSHAPE::NFacets+3);
    // compute the number of shape functions of the top and bottom fluxes
    {
        TPZShapeData locdata;
        const int connect = TSHAPE::NFacets+1;
        TPZShapeHDivBound<TSHAPE> bound;
        bound.Initialize(ids,connectorders[connect],sideorient[TSHAPE::NFacets],locdata);
        const int nphibound = bound.NShape(locdata);
        data.fHDiv.fNumConnectShape[connect] = nphibound;
        nphi += nphibound;
        data.fH1.fPhi.Resize(nphi,data.fH1.fPhi.Cols());
        data.fH1.fDPhi.Resize(data.fH1.fDPhi.Rows(),nphi);
    }
    {
        TPZShapeData locdata;
        const int connect = TSHAPE::NFacets+2;
        TPZShapeHDivBound<TSHAPE> bound;
        bound.Initialize(ids,connectorders[connect],sideorient[TSHAPE::NFacets+1],locdata);
        const int nphibound = bound.NShape(locdata);
        data.fHDiv.fNumConnectShape[connect] = nphibound;
        nphi += nphibound;
        data.fH1.fPhi.Resize(nphi,data.fH1.fPhi.Cols());
        data.fH1.fDPhi.Resize(data.fH1.fDPhi.Rows(),nphi);
    }
}

template <class TSHAPE>
int TPZShapeHDivCollapsed<TSHAPE>::NShapeF(TPZShapeData &data)
{
    int nshape = 0;
    for(int i = 0; i< TSHAPE::NFacets+3; i++)
    {
        nshape += data.fHDiv.fNumConnectShape[i];
    }
    return nshape;
}

template <class TSHAPE>
void TPZShapeHDivCollapsed<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{
    const int nshape = NShapeF(data);
    divphi.Resize(nshape,1);
    divphi.Zero();
    phi.Resize(TSHAPE::Dimension+1,nshape);
    phi.Zero();

    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    TPZShapeH1<TSHAPE>::Shape(pt,data,data.fH1.fPhi,data.fH1.fDPhi);
    for(int i = 0; i< data.fHDiv.fSDVecShapeIndex.size(); i++)
    {
        auto it = data.fHDiv.fSDVecShapeIndex[i];
        int vecindex = it.first;
        int scalindex = it.second;
        divphi(i,0) = 0.;
        for(int d = 0; d<= TSHAPE::Dimension; d++)
        {
            phi(d,i) = data.fH1.fPhi(scalindex,0)*data.fHDiv.fMasterDirections(d,vecindex);
        }
        for(int d = 0; d<TSHAPE::Dimension; d++)
        {
            divphi(i,0) += data.fH1.fDPhi(d,scalindex)*data.fHDiv.fMasterDirections(d,vecindex);
        }
    }
    // compute the number of shape functions of the top and bottom fluxes
    int firstshape = data.fHDiv.fSDVecShapeIndex.size();
    {
        TPZShapeData locdata;
        const int connect = TSHAPE::NFacets+1;
        TPZShapeHDivBound<TSHAPE>::Initialize(data.fCornerNodeIds,data.fHDiv.fConnectOrders[connect],data.fHDiv.fSideOrient[TSHAPE::NFacets],locdata);
        TPZFNMatrix<25> locphi(data.fHDiv.fNumConnectShape[connect],1);
        TPZShapeHDivBound<TSHAPE>::Shape(pt, locdata, locphi);
        for (int ish = 0; ish < locphi.Rows(); ish++) {
            phi(dim,firstshape+ish) = locphi(ish,0);
            divphi(firstshape+ish) = locphi(ish,0);
        }
        firstshape += locphi.Rows();
    }
    {
        TPZShapeData locdata;
        const int connect = TSHAPE::NFacets+2;
        TPZShapeHDivBound<TSHAPE>::Initialize(data.fCornerNodeIds,data.fHDiv.fConnectOrders[connect],data.fHDiv.fSideOrient[TSHAPE::NFacets+1],locdata);
        TPZFNMatrix<25> locphi(data.fHDiv.fNumConnectShape[connect],1);
        TPZShapeHDivBound<TSHAPE>::Shape(pt, locdata, locphi);
        for (int ish = 0; ish < locphi.Rows(); ish++) {
            phi(dim,firstshape+ish) = locphi(ish,0);
            divphi(firstshape+ish) = locphi(ish,0);
        }
    }

}

template
struct TPZShapeHDivCollapsed<pzshape::TPZShapeLinear>;

template
struct TPZShapeHDivCollapsed<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivCollapsed<pzshape::TPZShapeQuad>;


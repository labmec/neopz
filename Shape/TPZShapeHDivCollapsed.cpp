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
    int nphi = data.fPhi.Rows();
    // extend the dimension of the master element directions
    const int ndirections = data.fMasterDirections.Cols();
    data.fMasterDirections.Resize(TSHAPE::Dimension+1, ndirections+2);
    for (int idir = 0; idir < ndirections+2; idir++) {
        data.fMasterDirections(dim,idir) = 0.;
    }
    for (int d = 0; d<dim+1; d++) {
        data.fMasterDirections(d,ndirections) = 0.;
        data.fMasterDirections(d,ndirections+1) = 0.;
    }
    data.fMasterDirections(dim,ndirections) = 1.;
    data.fMasterDirections(dim,ndirections+1) = -1.;
    // extend the dimensions of the connect orders
    data.fHDivConnectOrders = connectorders;
    data.fSideOrient = sideorient;
    data.fHDivNumConnectShape.Resize(TSHAPE::NFacets+3);
    // compute the number of shape functions of the top and bottom fluxes
    {
        TPZShapeData locdata;
        const int connect = TSHAPE::NFacets+1;
        TPZShapeHDivBound<TSHAPE> bound;
        bound.Initialize(ids,connectorders[connect],sideorient[TSHAPE::NFacets],locdata);
        const int nphibound = bound.NShape(locdata);
        data.fHDivNumConnectShape[connect] = nphibound;
        nphi += nphibound;
        data.fPhi.Resize(nphi,data.fPhi.Cols());
        data.fDPhi.Resize(data.fDPhi.Rows(),nphi);
    }
    {
        TPZShapeData locdata;
        const int connect = TSHAPE::NFacets+2;
        TPZShapeHDivBound<TSHAPE> bound;
        bound.Initialize(ids,connectorders[connect],sideorient[TSHAPE::NFacets+1],locdata);
        const int nphibound = bound.NShape(locdata);
        data.fHDivNumConnectShape[connect] = nphibound;
        nphi += nphibound;
        data.fPhi.Resize(nphi,data.fPhi.Cols());
        data.fDPhi.Resize(data.fDPhi.Rows(),nphi);
    }
}

template <class TSHAPE>
int TPZShapeHDivCollapsed<TSHAPE>::NShapeF(TPZShapeData &data)
{
    int nshape = 0;
    for(int i = 0; i< TSHAPE::NFacets+3; i++)
    {
        nshape += data.fHDivNumConnectShape[i];
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
    TPZShapeH1<TSHAPE>::Shape(pt,data);
    for(int i = 0; i< data.fSDVecShapeIndex.size(); i++)
    {
        auto it = data.fSDVecShapeIndex[i];
        int vecindex = it.first;
        int scalindex = it.second;
        divphi(i,0) = 0.;
        for(int d = 0; d<= TSHAPE::Dimension; d++)
        {
            phi(d,i) = data.fPhi(scalindex,0)*data.fMasterDirections(d,vecindex);
        }
        for(int d = 0; d<TSHAPE::Dimension; d++)
        {
            divphi(i,0) += data.fDPhi(d,scalindex)*data.fMasterDirections(d,vecindex);
        }
    }
    // compute the number of shape functions of the top and bottom fluxes
    int firstshape = data.fSDVecShapeIndex.size();
    {
        TPZShapeData locdata;
        const int connect = TSHAPE::NFacets+1;
        TPZShapeHDivBound<TSHAPE>::Initialize(data.fCornerNodeIds,data.fHDivConnectOrders[connect],data.fSideOrient[TSHAPE::NFacets],locdata);
        TPZFNMatrix<25> locphi(data.fHDivNumConnectShape[connect],1);
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
        TPZShapeHDivBound<TSHAPE>::Initialize(data.fCornerNodeIds,data.fHDivConnectOrders[connect],data.fSideOrient[TSHAPE::NFacets+1],locdata);
        TPZFNMatrix<25> locphi(data.fHDivNumConnectShape[connect],1);
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


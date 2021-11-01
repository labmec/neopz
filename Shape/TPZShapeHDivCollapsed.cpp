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
    if(ids.size() != ncorner || connectorders.size() != nsides-ncorner+2)
    {
        DebugStop();
    }
    TPZManVector<int,20> hdivconnectorders(connectorders);
    hdivconnectorders.Resize(nsides-ncorner, 0);
    TPZManVector<int,20> hdivsideorient(sideorient);
    hdivsideorient.Resize(nsides-ncorner, 0);
    TPZShapeHDiv<TSHAPE>::Initialize(ids, hdivconnectorders, hdivsideorient, data);
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
    // compute the number of shape functions of the top and bottom fluxes
    {
        TPZShapeData locdata;
        const int connect = nsides-ncorner;
        TPZShapeHDivBound<TSHAPE> bound;
        bound.Initialize(ids,connectorders[connect],sideorient[connect],locdata);
        data.fHDivNumConnectShape[connect] = locdata.fHDivNumConnectShape[0];
    }
    {
        TPZShapeData locdata;
        const int connect = nsides-ncorner+1;
        TPZShapeHDivBound<TSHAPE> bound;
        bound.Initialize(ids,connectorders[connect],sideorient[connect],locdata);
        data.fHDivNumConnectShape[connect] = locdata.fHDivNumConnectShape[0];
    }
}



template <class TSHAPE>
void TPZShapeHDivCollapsed<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{
    divphi.Resize(data.fSDVecShapeIndex.size(),1);
    divphi.Zero();
    phi.Resize(TSHAPE::Dimension,data.fSDVecShapeIndex.size());
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
        const int connect = nsides-ncorner;
        TPZShapeHDivBound<TSHAPE> bound;
        bound.Initialize(data.fCornerNodeIds,data.fHDivConnectOrders[connect],data.fSideOrient[connect],locdata);
        TPZFNMatrix<25> locphi(data.fHDivNumConnectShape[connect],1);
        bound.Shape(pt, locdata, locphi);
        int firstshape = data.fHDivNumConnectShape[connect];
        for (int ish = 0; ish < locphi.Rows(); ish++) {
            phi(dim,firstshape+ish) = locphi(ish,0);
            divphi(firstshape+ish) = locphi(ish,0);
        }
        firstshape += locphi.Rows();
    }
    {
        TPZShapeData locdata;
        const int connect = nsides-ncorner+1;
        TPZShapeHDivBound<TSHAPE> bound;
        bound.Initialize(data.fCornerNodeIds,data.fHDivConnectOrders[connect],data.fSideOrient[connect],locdata);
        TPZFNMatrix<25> locphi(data.fHDivNumConnectShape[connect],1);
        bound.Shape(pt, locdata, locphi);
        for (int ish = 0; ish < locphi.Rows(); ish++) {
            phi(dim,firstshape+ish) = -locphi(ish,0);
            divphi(firstshape+ish) = -locphi(ish,0);
        }
    }

}

template
struct TPZShapeHDivCollapsed<pzshape::TPZShapeLinear>;

template
struct TPZShapeHDivCollapsed<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivCollapsed<pzshape::TPZShapeQuad>;


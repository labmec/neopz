#include "TPZShapeNewHDiv.h"
#include "TPZShapeH1.h"
#include "TPZShapeHCurl.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "TPZShapeData.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.shapehdiv");
#endif

template <class TSHAPE>
TPZShapeNewHDiv<TSHAPE>::TPZShapeNewHDiv() {}
template <class TSHAPE>
int TPZShapeNewHDiv<TSHAPE>::NConnectShapeF(int connect, const TPZShapeData &shapedata)
{
    return shapedata.fHDiv.fNumConnectShape[connect];
}

template <class TSHAPE>
int TPZShapeNewHDiv<TSHAPE>::NShapeF(const TPZShapeData &shapedata)
{
    const int nconnect = shapedata.fHDiv.fNumConnectShape.size();
    int nshape = 0;
    for (int ic = 0; ic < nconnect; ic++)
        nshape += shapedata.fHDiv.fNumConnectShape[ic];
    return nshape;
}

template <class TSHAPE>
int TPZShapeNewHDiv<TSHAPE>::ComputeNConnectShapeF(int connect, int order)
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets)
    {
        DebugStop();
    }
#endif
    MElementType thistype = TSHAPE::Type();

    // For the facet connects, it should return the same number of functions as HDivConstant
    // For the internal connect, it should return the same number of functions as HDiv Standard

    if (thistype == EOned)
    {
        if (connect < 2)
            return 0;
        else
            return order;
    }
    else if (thistype == ETriangle)
    {
        if (connect < TSHAPE::NFacets)
            return (order + 1);
        else
            return (order + 1) * (order + 1) - 1;
    }
    else if (thistype == EQuadrilateral)
    {
        if (connect < TSHAPE::NFacets)
            return (order + 1);
        else
            return 2 * order * (order + 1);
    }
    else if (thistype == ETetraedro)
    {
        if (connect < TSHAPE::NFacets)
            return (order + 1) * (order + 2) / 2;
        else
            return order * (order + 2) * (order + 3) / 2;
    }
    else if (thistype == EPrisma)
    {
        if (connect == 0 || connect == 4)
            return (order + 1) * (order + 2) / 2;
        else if (connect < TSHAPE::NFacets)
            return (order + 1) * (order + 1);
        else
            return order * order * (3 * order + 5) / 2 + 7 * order - 2;
    }
    else if (thistype == ECube)
    {
        if (connect < TSHAPE::NFacets)
            return (order + 1) * (order + 1);
        else
            return 3 * order * (order + 1) * (order + 1);
    }
    DebugStop();
    unreachable();
}

template struct TPZShapeNewHDiv<pzshape::TPZShapeLinear>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeTriang>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeQuad>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeTetra>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeCube>;

template struct TPZShapeNewHDiv<pzshape::TPZShapePrism>;
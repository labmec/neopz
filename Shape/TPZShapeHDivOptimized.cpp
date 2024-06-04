#include "TPZShapeHDivOptimized.h"
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
TPZShapeHDivOptimized<TSHAPE>::TPZShapeHDivOptimized() {}
template <class TSHAPE>
int TPZShapeHDivOptimized<TSHAPE>::NConnectShapeF(int connect, const TPZShapeData &shapedata)
{
    return shapedata.fHDiv.fNumConnectShape[connect];
}

template <class TSHAPE>
int TPZShapeHDivOptimized<TSHAPE>::NShapeF(const TPZShapeData &shapedata)
{
    const int nconnect = shapedata.fHDiv.fNumConnectShape.size();
    int nshape = 0;
    for (int ic = 0; ic < nconnect; ic++)
        nshape += shapedata.fHDiv.fNumConnectShape[ic];
    return nshape;
}

template <class TSHAPE>
int TPZShapeHDivOptimized<TSHAPE>::ComputeNConnectShapeF(int connect, int order)
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

template <class TSHAPE>
void TPZShapeHDivOptimized<TSHAPE>::CheckH1ConnectOrder(const TPZVec<int> &connectorders, TPZVec<int> &H1Orders)
{
    constexpr int nHDivcon = TSHAPE::NFacets + 1;
    constexpr int nHCurlcon = TSHAPE::NSides - TSHAPE::NCornerNodes;
    constexpr int dim = TSHAPE::Dimension;

    // H1 order required by HDiv Shape
    const int maxorder = connectorders[TSHAPE::NFacets] + 1;
    TPZManVector<int, 27> H1HDivOrders(nHCurlcon, maxorder);

    if constexpr (dim == 2)
    {
        H1Orders.resize(nHDivcon);
        for (int i = 0; i < nHDivcon; i++)
        {
            H1Orders[i] = H1HDivOrders[i];
        }
    }
    else if constexpr (dim == 3)
    {
        // H1 order required by HCurl Shape
        constexpr int nedges = TSHAPE::NSides - TSHAPE::NFacets - TSHAPE::NCornerNodes - 1;
        TPZManVector<int, 27> HCurlOrders(nHCurlcon, 1);
        for (int ic = nedges; ic < nHCurlcon; ic++)
        {
            HCurlOrders[ic] = connectorders[ic - nedges];
        }
        if (TSHAPE::Type() == ETetraedro)
        {
            for (int ic = nedges; ic < nHCurlcon; ic++)
            {
                HCurlOrders[ic]++;
            }
            HCurlOrders[nHCurlcon - 1]++;
        }
        TPZManVector<int, 27> H1HCurlOrders;
        TPZShapeHCurl<TSHAPE>::CalcH1ShapeOrders(HCurlOrders, H1HCurlOrders);

        // Maximum between H1HDiv and H1HCurl orders
        H1Orders.resize(nHCurlcon);
        for (int i = 0; i < nHCurlcon; i++)
        {
            H1Orders[i] = std::max(H1HDivOrders[i], H1HCurlOrders[i]);
        }
    }
}

template struct TPZShapeHDivOptimized<pzshape::TPZShapeLinear>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapeTriang>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapeQuad>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapeTetra>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapeCube>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapePrism>;
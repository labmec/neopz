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
void TPZShapeHDivOptimized<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
                                         const TPZVec<int> &connectorders,
                                         const TPZVec<int> &sideorient,
                                         TPZShapeData &data)
{
    constexpr int nHDivcon = TSHAPE::NFacets + 1;
    constexpr int nHCurlcon = TSHAPE::NSides - TSHAPE::NCornerNodes;
    if (connectorders.size() != nHDivcon)
        DebugStop();

    data.fCornerNodeIds = ids;
    data.fHDiv.fSideOrient = sideorient;

    // The H1 connects order is used for both HDiv and HCurl shapes.
    // Thus, we compute the required order for both shapes and take the maximum between them to construct VecAndShape
    TPZManVector<int, 27> H1Orders;
    CheckH1ConnectOrder(connectorders, H1Orders);

    TPZShapeH1<TSHAPE>::Initialize(data.fCornerNodeIds, H1Orders, data);

    // Initialize the HDiv structure
    data.fHDiv.fConnectOrders = connectorders;

    data.fHDiv.fNumConnectShape.Resize(nHDivcon);
    int nShape = 0;
    for (int i = 0; i < nHDivcon; i++)
    {
        const int order = data.fHDiv.fConnectOrders[i];
        data.fHDiv.fNumConnectShape[i] = ComputeNConnectShapeF(i, order);
        nShape += data.fHDiv.fNumConnectShape[i];
    }

    data.fHDiv.fSDVecShapeIndex.Resize(nShape);

    TPZShapeHDiv<TSHAPE>::ComputeMasterDirections(data);
    TPZShapeHDiv<TSHAPE>::ComputeVecandShape(data);

    // Checks if the last connect order is >= then the other connects
    const int maxOrder = data.fHDiv.fConnectOrders[nHDivcon - 1];
    for (int i = 0; i < nHDivcon - 1; i++)
    {
        if (data.fHDiv.fConnectOrders[i] > maxOrder)
        {
            DebugStop();
        }
    }

    // Initialize the HCurl structure
    constexpr int nedges = TSHAPE::NSides - TSHAPE::NFacets - TSHAPE::NCornerNodes - 1;
    data.fHCurl.fConnectOrders.resize(nHCurlcon);
    data.fHCurl.fConnectOrders.Fill(1);
    for (int ic = nedges; ic < nHCurlcon; ic++)
    {
        data.fHCurl.fConnectOrders[ic] = connectorders[ic - nedges];
    }
    if (TSHAPE::Type() == ETetraedro)
    {
        for (int ic = nedges; ic < nHCurlcon; ic++)
        {
            data.fHCurl.fConnectOrders[ic]++;
        }
        data.fHCurl.fConnectOrders[nHCurlcon - 1]++;
    }

    data.fH1.fSideTransformationId.Resize(nHCurlcon, 0);
    for (int iside = TSHAPE::NCornerNodes; iside < TSHAPE::NSides; iside++)
    {
        int pos = iside - TSHAPE::NCornerNodes;
        int trans_id = TSHAPE::GetTransformId(iside, ids); // Foi criado
        data.fH1.fSideTransformationId[iside - TSHAPE::NCornerNodes] = trans_id;
    }

    data.fHCurl.fNumConnectShape.Resize(nHCurlcon);
    nShape = 0;
    for (int i = 0; i < nHCurlcon; i++)
    {
        data.fHCurl.fNumConnectShape[i] = TPZShapeHCurl<TSHAPE>::ComputeNConnectShapeF(i, data.fHCurl.fConnectOrders[i]);
        nShape += data.fHCurl.fNumConnectShape[i];
    }

    data.fHCurl.fSDVecShapeIndex.Resize(nShape);
    TPZFNMatrix<9, REAL> gradX(TSHAPE::Dimension, TSHAPE::Dimension, 0);
    gradX.Identity();

    data.fHCurl.fMasterDirections.Redim(TSHAPE::Dimension, 3 * TSHAPE::NSides);
    TSHAPE::ComputeHCurlDirections(gradX, data.fHCurl.fMasterDirections, data.fH1.fSideTransformationId);

    TPZShapeHCurl<TSHAPE>::ComputeVecandShape(data);

    data.fHCurl.fNumConnectShape.Resize(nHCurlcon);
    // we need to update the number of filtered hcurl functions
    for (int i = 0; i < nHCurlcon; i++)
    {
        data.fHCurl.fNumConnectShape[i] = TPZShapeHCurlNoGrads<TSHAPE>::ComputeNConnectShapeF(i, data.fHCurl.fConnectOrders[i]);
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        data.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}
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
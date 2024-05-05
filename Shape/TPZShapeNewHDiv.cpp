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
void TPZShapeNewHDiv<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
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
void TPZShapeNewHDiv<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{

    constexpr int ncorner = TSHAPE::NCornerNodes;
    constexpr int nsides = TSHAPE::NSides;
    constexpr int dim = TSHAPE::Dimension;
    constexpr int nfacets = TSHAPE::NFacets;
    const int nedges = TSHAPE::NumSides(1);

    if (phi.Rows() != dim || phi.Cols() != data.fHDiv.fSDVecShapeIndex.size())
    {
        phi.Resize(dim, data.fHDiv.fSDVecShapeIndex.size());
        phi.Zero();
    }
    if (divphi.Rows() != data.fHDiv.fSDVecShapeIndex.size())
    {
        divphi.Resize(data.fHDiv.fSDVecShapeIndex.size(), 1);
        divphi.Zero();
    }

    // Compute constant Hdiv functions
    TPZFNMatrix<dim * nfacets, REAL> RT0phi(dim, nfacets);
    TPZManVector<REAL, nfacets> RT0div(nfacets);
    RT0phi.Zero();
    RT0div.Fill(0.);
    TSHAPE::ComputeConstantHDiv(pt, RT0phi, RT0div);

    // For dim = 2, we use the gradient of H1 functions to compute the facet HDiv functions, while the internal functions come from the standard HDiv
    if constexpr (dim == 2)
    {
        TPZShapeH1<TSHAPE>::Shape(pt, data, data.fH1.fPhi, data.fH1.fDPhi);

        int count = 0;
        int countKernel = ncorner;
        // Edge functions
        for (int i = 0; i < nedges; i++)
        {
            // RT0 Function
            phi(0, count) = RT0phi(0, i) * data.fHDiv.fSideOrient[i];
            phi(1, count) = RT0phi(1, i) * data.fHDiv.fSideOrient[i];
            divphi(count, 0) = RT0div[i] * data.fHDiv.fSideOrient[i];
            count++;

            // Kernel Hdiv
            for (int j = 1; j < data.fHDiv.fNumConnectShape[i]; j++)
            {
                phi(0, count) = -data.fH1.fDPhi(1, countKernel);
                phi(1, count) = data.fH1.fDPhi(0, countKernel);
                count++;
                countKernel++;
            }
        }

        // Internal functions
        const int nfacetfunc = count;
        for (int i = nfacetfunc; i < data.fHDiv.fSDVecShapeIndex.size(); i++)
        {
            auto it = data.fHDiv.fSDVecShapeIndex[i];
            int vecindex = it.first;
            int scalindex = it.second;
            divphi(i, 0) = 0.;
            for (int d = 0; d < dim; d++)
            {
                phi(d, i) = data.fH1.fPhi(scalindex, 0) * data.fHDiv.fMasterDirections(d, vecindex);
                divphi(i, 0) += data.fH1.fDPhi(d, scalindex) * data.fHDiv.fMasterDirections(d, vecindex);
            }
            count++;
        }
    }
    // For dim = 3, the facet functions come from HCurlNoGrads, while the internal functions are computed according to standard HDiv
    else if constexpr (dim == 3)
    {
        divphi.Zero();
        int nshapehcurl = TPZShapeHCurlNoGrads<TSHAPE>::NHCurlShapeF(data);
        int nshape = NShapeF(data);

        TPZFNMatrix<200, REAL> phiAux(dim, nshapehcurl), curlPhiAux(3, nshapehcurl);
        phiAux.Zero();
        curlPhiAux.Zero();

        TPZShapeHCurlNoGrads<TSHAPE>::Shape(pt, data, phiAux, curlPhiAux);

        int count = 0;
        int countKernel = nedges;

        // Face functions
        for (int i = 0; i < nfacets; i++)
        {
            // RT0 Function
            for (auto d = 0; d < dim; d++)
            {
                phi(d, count) = RT0phi(d, i) * data.fHDiv.fSideOrient[i];
            }
            divphi(count, 0) = RT0div[i] * data.fHDiv.fSideOrient[i];
            count++;

            // Kernel HDiv functions
            for (int k = 0; k < data.fHCurl.fNumConnectShape[nedges + i]; k++)
            {
                for (auto d = 0; d < dim; d++)
                {
                    phi(d, count) = curlPhiAux(d, countKernel);
                }
                countKernel++;
                count++;
            }
        }
        // Internal functions
        const int nfacetfunc = count;
        TPZShapeH1<TSHAPE>::Shape(pt, data, data.fH1.fPhi, data.fH1.fDPhi);
        for (int i = nfacetfunc; i < data.fHDiv.fSDVecShapeIndex.size(); i++)
        {
            auto it = data.fHDiv.fSDVecShapeIndex[i];
            int vecindex = it.first;
            int scalindex = it.second;
            divphi(i, 0) = 0.;
            for (int d = 0; d < dim; d++)
            {
                phi(d, i) = data.fH1.fPhi(scalindex, 0) * data.fHDiv.fMasterDirections(d, vecindex);
                divphi(i, 0) += data.fH1.fDPhi(d, scalindex) * data.fHDiv.fMasterDirections(d, vecindex);
            }
            count++;
        }
        if (count != nshape)
            DebugStop();
    }
    else
    {
        DebugStop();
    }
}

template <class TSHAPE>
void TPZShapeNewHDiv<TSHAPE>::Shape(const TPZVec<Fad<REAL>> &pt, TPZShapeData &data, TPZFMatrix<Fad<REAL>> &phi, TPZFMatrix<Fad<REAL>> &divphi)
{

    constexpr int ncorner = TSHAPE::NCornerNodes;
    constexpr int nsides = TSHAPE::NSides;
    constexpr int dim = TSHAPE::Dimension;
    constexpr int nfacets = TSHAPE::NFacets;
    const int nedges = TSHAPE::NumSides(1);
    int fadsize = pt[0].size();

    if (phi.Rows() != dim || phi.Cols() != data.fHDiv.fSDVecShapeIndex.size())
    {
        phi.Resize(dim, data.fHDiv.fSDVecShapeIndex.size());
        phi.Zero();
    }
    if (divphi.Rows() != data.fHDiv.fSDVecShapeIndex.size())
    {
        divphi.Resize(data.fHDiv.fSDVecShapeIndex.size(), 1);
        divphi.Zero();
    }

    // Compute constant Hdiv functions
    TPZFNMatrix<dim * nfacets, Fad<REAL>> RT0phi(dim, nfacets);
    TPZManVector<Fad<REAL>, nfacets> RT0div(nfacets);
    RT0phi.Zero();
    RT0div.Fill(0.);
    TSHAPE::ComputeConstantHDiv(pt, RT0phi, RT0div);

    // For dim = 2, we use the gradient of H1 functions to compute the facet HDiv functions, while the internal functions come from the standard HDiv
    if constexpr (dim == 2)
    {
        TPZFNMatrix<9, Fad<REAL>> locphi(data.fH1.fPhi.Rows(), 1), locdphi(dim, data.fH1.fPhi.Rows());
        TPZShapeH1<TSHAPE>::Shape(pt, data, locphi, locdphi);

        int count = 0;
        int countKernel = ncorner;
        // Edge functions
        for (int i = 0; i < nedges; i++)
        {
            // RT0 Function
            phi(0, count) = RT0phi(0, i) * data.fHDiv.fSideOrient[i];
            phi(1, count) = RT0phi(1, i) * data.fHDiv.fSideOrient[i];
            divphi(count, 0) = RT0div[i] * data.fHDiv.fSideOrient[i];
            count++;

            // Kernel Hdiv
            for (int j = 1; j < data.fHDiv.fNumConnectShape[i]; j++)
            {
                phi(0, count) = -locdphi(1, countKernel);
                phi(1, count) = locdphi(0, countKernel);
                count++;
                countKernel++;
            }
        }

        // Internal functions
        const int nfacetfunc = count;
        for (int i = nfacetfunc; i < data.fHDiv.fSDVecShapeIndex.size(); i++)
        {
            auto it = data.fHDiv.fSDVecShapeIndex[i];
            int vecindex = it.first;
            int scalindex = it.second;
            divphi(i, 0) = 0.;
            for (int d = 0; d < dim; d++)
            {
                phi(d, i) = locphi(scalindex, 0) * data.fHDiv.fMasterDirections(d, vecindex);
                divphi(i, 0) += locdphi(d, scalindex) * data.fHDiv.fMasterDirections(d, vecindex);
            }
            count++;
        }
    }
    // For dim = 3, the facet functions come from HCurlNoGrads, while the internal functions are computed according to standard HDiv
    else if constexpr (dim == 3)
    {
        divphi.Zero();
        int nshapehcurl = TPZShapeHCurlNoGrads<TSHAPE>::NHCurlShapeF(data);
        int nshape = NShapeF(data);

        TPZFNMatrix<200, Fad<REAL>> phiAux(dim, nshapehcurl), curlPhiAux(3, nshapehcurl);
        phiAux.Zero();
        curlPhiAux.Zero();

        TPZShapeHCurlNoGrads<TSHAPE>::Shape(pt, data, phiAux, curlPhiAux);

        int count = 0;
        int countKernel = nedges;

        // Face functions
        for (int i = 0; i < nfacets; i++)
        {
            // RT0 Function
            for (auto d = 0; d < dim; d++)
            {
                phi(d, count) = RT0phi(d, i) * data.fHDiv.fSideOrient[i];
            }
            divphi(count, 0) = RT0div[i] * data.fHDiv.fSideOrient[i];
            count++;

            // Kernel HDiv functions
            for (int k = 0; k < data.fHCurl.fNumConnectShape[nedges + i]; k++)
            {
                for (auto d = 0; d < dim; d++)
                {
                    phi(d, count) = curlPhiAux(d, countKernel);
                }
                countKernel++;
                count++;
            }
        }
        // Internal functions
        const int nfacetfunc = count;
        TPZFNMatrix<9, Fad<REAL>> locphi(data.fH1.fPhi.Rows(), 1), locdphi(dim, data.fH1.fPhi.Rows());
        TPZShapeH1<TSHAPE>::Shape(pt, data, locphi, locdphi);
        for (int i = nfacetfunc; i < data.fHDiv.fSDVecShapeIndex.size(); i++)
        {
            auto it = data.fHDiv.fSDVecShapeIndex[i];
            int vecindex = it.first;
            int scalindex = it.second;
            divphi(i, 0) = 0.;
            for (int d = 0; d < dim; d++)
            {
                phi(d, i) = locphi(scalindex, 0) * data.fHDiv.fMasterDirections(d, vecindex);
                divphi(i, 0) += locdphi(d, scalindex) * data.fHDiv.fMasterDirections(d, vecindex);
            }
            count++;
        }
        if (count != nshape)
            DebugStop();
        if (countKernel != nshapehcurl)
            DebugStop();
    }
    else
    {
        DebugStop();
    }
}

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

template <class TSHAPE>
void TPZShapeNewHDiv<TSHAPE>::CheckH1ConnectOrder(const TPZVec<int> &connectorders, TPZVec<int> &H1Orders)
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

template struct TPZShapeNewHDiv<pzshape::TPZShapeLinear>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeTriang>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeQuad>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeTetra>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeCube>;

template struct TPZShapeNewHDiv<pzshape::TPZShapePrism>;
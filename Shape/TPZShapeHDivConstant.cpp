#include "TPZShapeHDivConstant.h"
#include "TPZShapeHCurlNoGrads.h"

#include "TPZShapeH1.h"
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
#include "TPZCompElHCurl.h"

//! Should be called once per element. Initializes the data structure
template <class TSHAPE>
void TPZShapeHDivConstant<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
                                              const TPZVec<int> &connectorders,
                                              const TPZVec<int> &sideorient,
                                              TPZShapeData &data)
{
    data.fHDiv.fSideOrient = sideorient;
    constexpr int conSize = TSHAPE::NFacets+1;
    data.fHDiv.fConnectOrders.Resize(conSize);
    data.fHDiv.fNumConnectShape.Resize(conSize);
    for (int ic = 0; ic < conSize; ic++)
    {
        data.fHDiv.fConnectOrders[ic] = connectorders[ic];
    }
    
    if constexpr(TSHAPE::Dimension == 2)
    {
        data.fH1.fConnectOrders.resize(conSize);

        // Compatibilize the polynomial order
        for (int i = 0; i < conSize; i++)
        {
            data.fH1.fConnectOrders[i] = connectorders[i]+1;
        }
        // if (TSHAPE::Type() == ETriangle){
        //     data.fH1.fConnectOrders[i][conSize-1]++;
        // }

        // Initialize H1 data structures
        TPZShapeH1<TSHAPE>::Initialize(ids, data.fH1.fConnectOrders, data);

        for (int ic = 0; ic < TSHAPE::NFacets; ic++)
        {
            data.fHDiv.fNumConnectShape[ic] = data.fH1.fNumConnectShape[ic] + 1;
        }
        int ic = TSHAPE::NFacets;
        data.fHDiv.fNumConnectShape[ic] = data.fH1.fNumConnectShape[ic];
    }
    else if constexpr(TSHAPE::Dimension == 3)
    {
        // Compatibilize the polynomial order
        constexpr int nHcurlCon = TSHAPE::NSides - TSHAPE::NCornerNodes;
        constexpr int nedges = TSHAPE::NSides - TSHAPE::NFacets - TSHAPE::NCornerNodes - 1;
        data.fHCurl.fConnectOrders.resize(nHcurlCon);
        data.fHCurl.fConnectOrders.Fill(1);
        for (int ic = nedges; ic < nHcurlCon; ic++)
        {
            data.fHCurl.fConnectOrders[ic] = connectorders[ic - nedges];
        }
        if constexpr(TSHAPE::Type() == ETetraedro)
        {
            for (int ic = nedges; ic < nHcurlCon; ic++)
            {
                data.fHCurl.fConnectOrders[ic]++;
            }
            data.fHCurl.fConnectOrders[nHcurlCon-1]++;
        }

        TPZShapeHCurlNoGrads<TSHAPE>::Initialize(ids, data.fHCurl.fConnectOrders, data);

        for (int ic = nedges; ic < TSHAPE::NSides - TSHAPE::NCornerNodes - 1; ic++)
        {
            const int numshape = data.fHCurl.fNumConnectShape[ic] + 1;
            data.fHDiv.fNumConnectShape[ic - nedges] = numshape;
        }
        const int ic = TSHAPE::NSides - TSHAPE::NCornerNodes - 1;
        const int numshape = data.fHCurl.fNumConnectShape[ic];
        data.fHDiv.fNumConnectShape[ic - nedges] = numshape;
    }
}

template <class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NHDivShapeF(TPZShapeData &data)
{
    int nshape = 0;
    int nc = data.fHDiv.fNumConnectShape.size();
    for (int ic = 0; ic < nc; ic++)
        nshape += NConnectShapeF(ic, data);
    return nshape;
}

template <class TSHAPE>
void TPZShapeHDivConstant<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{

    constexpr int ncorner = TSHAPE::NCornerNodes;
    constexpr int nsides = TSHAPE::NSides;
    constexpr int dim = TSHAPE::Dimension;
    constexpr int nfacets = TSHAPE::NFacets;
    const int nedges = TSHAPE::NumSides(1);

    // Compute constant Hdiv functions
    TPZFNMatrix<dim*nfacets, REAL> vecDiv(dim, nfacets);
    TPZManVector<REAL,nfacets> div(nfacets);
    vecDiv.Zero();
    div.Fill(0.);
    TSHAPE::ComputeConstantHDiv(pt, vecDiv, div);

    int nshape = data.fH1.fPhi.Rows();

    if constexpr(dim == 2)
    {
        TPZShapeH1<TSHAPE>::Shape(pt, data, data.fH1.fPhi, data.fH1.fDPhi);
        divphi.Zero();

        int count = 0;
        int countKernel = ncorner;
        // Edge functions
        for (int i = 0; i < nedges; i++)
        {
            // RT0 Function
            phi(0, count) = vecDiv(0, i) * data.fHDiv.fSideOrient[i];
            phi(1, count) = vecDiv(1, i) * data.fHDiv.fSideOrient[i];
            divphi(count, 0) = div[i] * data.fHDiv.fSideOrient[i];
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
        for (int i = countKernel; i < nshape; i++)
        {
            phi(0, count) = -data.fH1.fDPhi(1, countKernel);
            phi(1, count) = data.fH1.fDPhi(0, countKernel);
            count++;
            countKernel++;
        }
    }
    else if constexpr(dim == 3)
    {
        divphi.Zero();
        int nshapehcurl = TPZShapeHCurlNoGrads<TSHAPE>::NHCurlShapeF(data);
        int nshape = NHDivShapeF(data);

        TPZFNMatrix<200,REAL> phiAux(dim, nshapehcurl), curlPhiAux(3, nshapehcurl);
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
                phi(d, count) = vecDiv(d, i) * data.fHDiv.fSideOrient[i];
            }
            divphi(count, 0) = div[i] * data.fHDiv.fSideOrient[i];
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
        // Internal Functions - HDivKernel
        for (int i = 0; i < data.fHCurl.fNumConnectShape[TSHAPE::NSides - TSHAPE::NCornerNodes - 1]; i++)
        {
            for (auto d = 0; d < dim; d++)
            {
                phi(d, count) = curlPhiAux(d, countKernel);
            }
            countKernel++;
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
void TPZShapeHDivConstant<TSHAPE>::Shape(const TPZVec<Fad<REAL>> &pt, TPZShapeData &data, TPZFMatrix<Fad<REAL>> &phi, TPZFMatrix<Fad<REAL>> &divphi)
{

    constexpr int ncorner = TSHAPE::NCornerNodes;
    constexpr int nsides = TSHAPE::NSides;
    constexpr int dim = TSHAPE::Dimension;
    constexpr int nfacets = TSHAPE::NFacets;
    const int nedges = TSHAPE::NumSides(1);

    // Compute constant Hdiv functions
    TPZFNMatrix<dim*nfacets, Fad<REAL>> vecDiv(dim, nfacets);
    TPZManVector<Fad<REAL>,nfacets> div(nfacets);
    vecDiv.Zero();
    div.Fill(0.);
    TSHAPE::ComputeConstantHDiv(pt, vecDiv, div);

    int nshape = data.fH1.fPhi.Rows();

    if constexpr(dim == 2)
    {
        TPZFNMatrix<9, Fad<REAL>> locphi(data.fH1.fPhi.Rows(), data.fH1.fPhi.Cols()), dphi(data.fH1.fDPhi.Rows(), data.fH1.fDPhi.Cols());
        TPZShapeH1<TSHAPE>::Shape(pt, data, locphi, dphi);
        divphi.Zero();

        int count = 0;
        int countKernel = ncorner;
        // Edge functions
        for (int i = 0; i < nedges; i++)
        {
            // RT0 Function
            phi(0, count) = vecDiv(0, i) * data.fHDiv.fSideOrient[i];
            phi(1, count) = vecDiv(1, i) * data.fHDiv.fSideOrient[i];
            divphi(count, 0) = div[i] * data.fHDiv.fSideOrient[i];
            count++;

            // Kernel Hdiv
            for (int j = 1; j < data.fHDiv.fNumConnectShape[i]; j++)
            {
                phi(0, count) = -dphi(1, countKernel);
                phi(1, count) = dphi(0, countKernel);
                count++;
                countKernel++;
            }
        }

        // Internal functions
        for (int i = countKernel; i < nshape; i++)
        {
            phi(0, count) = -dphi(1, countKernel);
            phi(1, count) = dphi(0, countKernel);
            count++;
            countKernel++;
        }
    }
    else if constexpr(dim == 3)
    {
        divphi.Zero();
        int nshapehcurl = TPZShapeHCurlNoGrads<TSHAPE>::NHCurlShapeF(data);
        int nshape = NHDivShapeF(data);

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
                phi(d, count) = vecDiv(d, i) * data.fHDiv.fSideOrient[i];
            }
            divphi(count, 0) = div[i] * data.fHDiv.fSideOrient[i];
            count++;

            // Kernel HDiv functions
            for (int k = 0; k < data.fHCurl.fNumConnectShape[nedges]; k++)
            {
                for (auto d = 0; d < dim; d++)
                {
                    phi(d, count) = curlPhiAux(d, countKernel);
                }
                countKernel++;
                count++;
            }
        }
        // Internal Functions - HDivKernel
        for (int i = 0; i < data.fHCurl.fNumConnectShape[TSHAPE::NSides - TSHAPE::NCornerNodes - 1]; i++)
        {
            for (auto d = 0; d < dim; d++)
            {
                phi(d, count) = curlPhiAux(d, countKernel);
            }
            countKernel++;
            count++;
        }
    }
    else
    {
        DebugStop();
    }
}

// icon is the connect index of the hdiv element
template <class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NConnectShapeF(int icon, TPZShapeData &data)
{
    // const int firstcon = TSHAPE::NSides-TSHAPE::NFacets-TSHAPE::NCornerNodes-1;
    // int faceconnect = icon+firstcon;
    int nshape = data.fHDiv.fNumConnectShape[icon];
    // if(icon < TSHAPE::NFacets) nshape++;
    return nshape;
}

template <class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(int connect, int order)
{

#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets)
    {
        DebugStop();
    }
#endif
    // int order = data.fHDiv.fConnectOrders[connect];
    MElementType thistype = TSHAPE::Type();

    if (thistype == EOned)
    {
        order++;
        if (connect < 2)
            return 0;
        else
            return order;
        // DebugStop();
    }
    else if (thistype == ETriangle)
    {
        order++;
        if (connect < TSHAPE::NFacets)
            return (order);
        else
        {
            // order++;
            return (order - 1) * (order - 2) / 2;
        }
    }
    else if (thistype == EQuadrilateral)
    {
        order++;
        if (connect < TSHAPE::NFacets)
            return (order);
        else
            return (order - 1) * (order - 1);
    }
    else if (thistype == ETetraedro)
    {
        order++;
        if (connect < TSHAPE::NFacets)
            return 1 + (order - 1) * (2 * order + 4) / 4;
        else
        {
            order++;
            return (order - 1) * (order - 2) * (2 * order + 3) / 6;
        }
    }
    else if (thistype == EPrisma)
    {
        DebugStop();
        // if(connect == 0 || connect == 4) return (order+1)*(order+2)/2;
        // else if(connect < TSHAPE::NFacets) return (order+1)*(order+1);
        // else return order*order*(3*order+5)/2+7*order-2;
    }
    else if (thistype == ECube)
    {
        if (connect < TSHAPE::NFacets)
            return 1 + order * (order + 2);
        else
            return order * order * (3 + 2 * order);
    }
    DebugStop();
    unreachable();
}

template struct TPZShapeHDivConstant<pzshape::TPZShapeLinear>;

template struct TPZShapeHDivConstant<pzshape::TPZShapeTriang>;

template struct TPZShapeHDivConstant<pzshape::TPZShapeQuad>;

template struct TPZShapeHDivConstant<pzshape::TPZShapeTetra>;

template struct TPZShapeHDivConstant<pzshape::TPZShapeCube>;

template struct TPZShapeHDivConstant<pzshape::TPZShapePrism>;

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
    
    if (TSHAPE::Dimension == 2)
    {
        //Creating a copy of the connect orders to avoid modifying the original data when calling SHAPEH1::Initialize()
        TPZManVector<int, 19> locconnectorders(connectorders);

        // Compatibilize the polynomial order
        int conSize = connectorders.size();
        for (int i = 0; i < conSize; i++)
        {
            locconnectorders[i]++;
        }
        // if (TSHAPE::Type() == ETriangle){
        //     locconnectorders[conSize-1]++;
        // }

        // Initialize data structures
        TPZShapeH1<TSHAPE>::Initialize(ids, locconnectorders, data);

        data.fHDiv.fConnectOrders.Resize(TSHAPE::NFacets + 1);
        data.fHDiv.fNumConnectShape.Resize(TSHAPE::NFacets + 1);
        for (int ic = 0; ic < TSHAPE::NFacets; ic++)
        {
            data.fHDiv.fConnectOrders[ic] = connectorders[ic];
            data.fHDiv.fNumConnectShape[ic] = data.fH1.fNumConnectShape[ic] + 1;
        }
        int ic = TSHAPE::NFacets;
        data.fHDiv.fConnectOrders[ic] = connectorders[ic];
        data.fHDiv.fNumConnectShape[ic] = data.fH1.fNumConnectShape[ic];
    }
    else
    {
        // Compatibilize the polynomial order
        TPZManVector<int, 19> adjustorder(connectorders);
        if (TSHAPE::Type() == ETetraedro)
        {
            int conSize = adjustorder.size();
            for (int i = 0; i < conSize; i++)
            {
                adjustorder[i]++;
            }
            adjustorder[conSize - 1]++;
        }

        constexpr int nedges = TSHAPE::NSides - TSHAPE::NFacets - TSHAPE::NCornerNodes - 1;
        TPZManVector<int, 19> locconnectorders(TSHAPE::NSides - TSHAPE::NCornerNodes, 1);
        for (int ic = nedges; ic < TSHAPE::NSides - TSHAPE::NCornerNodes; ic++)
        {
            locconnectorders[ic] = adjustorder[ic - nedges];
        }

        // ShapeHCurlNoGrads modifies ShapeData according to the HCurl space data structure.
        // after its initialization, we need to adjust the data structure to the HDiv space.
        TPZShapeHCurlNoGrads<TSHAPE>::Initialize(ids, locconnectorders, data);

        TPZManVector<int,60> HCurlConnectOrders = data.fHDiv.fConnectOrders;
        TPZManVector<int,60> HCurlNumConnectShapeF = data.fHDiv.fNumConnectShape;
        data.fHDiv.fConnectOrders.Resize(TSHAPE::NFacets + 1);
        data.fHDiv.fNumConnectShape.Resize(TSHAPE::NFacets + 1);
        for (int ic = 0; ic < TSHAPE::NFacets+1; ic++)
        {
            data.fHDiv.fConnectOrders[ic] = connectorders[ic];
        }
        for (int ic = nedges; ic < TSHAPE::NSides - TSHAPE::NCornerNodes - 1; ic++)
        {
            int numshape = HCurlNumConnectShapeF[ic] + 1;
            data.fHDiv.fNumConnectShape[ic - nedges] = numshape;
        }
        int ic = TSHAPE::NSides - TSHAPE::NCornerNodes - 1;
        int numshape = HCurlNumConnectShapeF[ic];
        data.fHDiv.fNumConnectShape[ic - nedges] = numshape;
    }
}

template <class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NHDivShapeF(TPZShapeData &data)
{
    int nshape = 0;
    // constexpr int firstConnect = TSHAPE::NSides-TSHAPE::NCornerNodes-TSHAPE::NFacets-1;
    int nc = data.fHDiv.fNumConnectShape.size();
    for (int ic = 0; ic < nc; ic++)
        nshape += NConnectShapeF(ic, data);
    // nshape += TSHAPE::NFacets;
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
    // std::cout << "fHDiv.fSide trans ID = " << data.fHDiv.fSideTransformationId << std::endl;
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
        //Adjusting the data structure to HCurl pattern
        TPZShapeData dataHCurl = data;
        constexpr int nHCurlCon = nsides - ncorner;
        dataHCurl.fHDiv.fConnectOrders.resize(nHCurlCon);
        dataHCurl.fHDiv.fConnectOrders.Fill(1);
        for (int ic = nedges; ic < nHCurlCon; ic++)
        {
            dataHCurl.fHDiv.fConnectOrders[ic] = data.fHDiv.fConnectOrders[ic - nedges];
        }
        if constexpr(TSHAPE::Type() == ETetraedro)
        {
            for (int ic = nedges; ic < nHCurlCon; ic++)
            {
                dataHCurl.fHDiv.fConnectOrders[ic]++;
            }
            int ic = nHCurlCon - 1;
            dataHCurl.fHDiv.fConnectOrders[ic]++;
        }

        dataHCurl.fHDiv.fNumConnectShape.resize(nHCurlCon);
        dataHCurl.fHDiv.fNumConnectShape.Fill(1);
        //For the facets, we subtract the constant function
        for (int ic = 0; ic < nfacets; ic++)
        {
            int numshape = data.fHDiv.fNumConnectShape[ic] - 1;
            dataHCurl.fHDiv.fNumConnectShape[nedges + ic] = numshape;
        }
        int ic = nHCurlCon - 1;
        int numshape = data.fHDiv.fNumConnectShape[nfacets];
        dataHCurl.fHDiv.fNumConnectShape[ic] = numshape;

        divphi.Zero();
        int nshapehcurl = TPZShapeHCurlNoGrads<TSHAPE>::NHCurlShapeF(dataHCurl);
        int nshape = NHDivShapeF(data);

        TPZFNMatrix<200,REAL> phiAux(dim, nshapehcurl), curlPhiAux(3, nshapehcurl);
        phiAux.Zero();
        curlPhiAux.Zero();

        TPZShapeHCurlNoGrads<TSHAPE>::Shape(pt, dataHCurl, phiAux, curlPhiAux);

        int count = 0;
        int countKernel = nedges;

        // Face functions
        for (int i = 0; i < nfacets; i++)
        {
            // std::cout << "Side orient - " << i << " " << data.fHDiv.fSideOrient[i] << std::endl;
            // RT0 Function
            for (auto d = 0; d < dim; d++)
            {
                phi(d, count) = vecDiv(d, i) * data.fHDiv.fSideOrient[i];
            }
            divphi(count, 0) = div[i] * data.fHDiv.fSideOrient[i];
            count++;

            // Kernel HDiv functions
            for (int k = 0; k < dataHCurl.fHDiv.fNumConnectShape[nedges + i]; k++)
            {
                for (auto d = 0; d < dim; d++)
                {
                    phi(d, count) = curlPhiAux(d, countKernel);
                }
                countKernel++;
                count++;
            }
            // if(i==0) {
            //     for(int k=0; k<count; k++) {
            //         std::cout << "phi ";
            //         for(int d=0; d<dim; d++) {
            //             std::cout << phi(d,k) << " ";
            //         }
            //         std::cout << std::endl;
            //     }
            // }
        }
        // Internal Functions - HDivKernel
        for (int i = 0; i < dataHCurl.fHDiv.fNumConnectShape[TSHAPE::NSides - TSHAPE::NCornerNodes - 1]; i++)
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
        // std::cout << "VecDiv = " << vecDiv << std::endl;
        // std::cout << "divphi = " << divphi << std::endl;
        // std::cout << "phi = " << phi << std::endl;
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
    // std::cout << "fHDiv.fSide trans ID = " << data.fHDiv.fSideTransformationId << std::endl;
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
        //Adjusting the data structure to HCurl pattern
        TPZShapeData dataHCurl = data;
        constexpr int nHCurlCon = nsides - ncorner;
        dataHCurl.fHDiv.fConnectOrders.resize(nHCurlCon);
        dataHCurl.fHDiv.fConnectOrders.Fill(1);
        for (int ic = nedges; ic < nHCurlCon; ic++)
        {
            dataHCurl.fHDiv.fConnectOrders[ic] = data.fHDiv.fConnectOrders[ic - nedges];
        }
        if constexpr(TSHAPE::Type() == ETetraedro)
        {
            for (int ic = nedges; ic < nHCurlCon; ic++)
            {
                dataHCurl.fHDiv.fConnectOrders[ic]++;
            }
            int ic = nHCurlCon - 1;
            dataHCurl.fHDiv.fConnectOrders[ic]++;
        }

        dataHCurl.fHDiv.fNumConnectShape.resize(nHCurlCon);
        dataHCurl.fHDiv.fNumConnectShape.Fill(1);
        //For the facets, we subtract the constant function
        for (int ic = 0; ic < nfacets; ic++)
        {
            int numshape = data.fHDiv.fNumConnectShape[ic] - 1;
            dataHCurl.fHDiv.fNumConnectShape[nedges + ic] = numshape;
        }
        int ic = nHCurlCon - 1;
        int numshape = data.fHDiv.fNumConnectShape[nfacets];
        dataHCurl.fHDiv.fNumConnectShape[ic] = numshape;

        divphi.Zero();
        int nshapehcurl = TPZShapeHCurlNoGrads<TSHAPE>::NHCurlShapeF(dataHCurl);
        int nshape = NHDivShapeF(data);

        TPZFNMatrix<200, Fad<REAL>> phiAux(dim, nshapehcurl), curlPhiAux(3, nshapehcurl);
        phiAux.Zero();
        curlPhiAux.Zero();

        TPZShapeHCurlNoGrads<TSHAPE>::Shape(pt, dataHCurl, phiAux, curlPhiAux);

        int count = 0;
        int countKernel = nedges;

        // Face functions
        for (int i = 0; i < nfacets; i++)
        {
            // std::cout << "Side orient - " << i << " " << data.fHDiv.fSideOrient[i] << std::endl;
            // RT0 Function
            for (auto d = 0; d < dim; d++)
            {
                phi(d, count) = vecDiv(d, i) * data.fHDiv.fSideOrient[i];
            }
            divphi(count, 0) = div[i] * data.fHDiv.fSideOrient[i];
            count++;

            // Kernel HDiv functions
            for (int k = 0; k < dataHCurl.fHDiv.fNumConnectShape[nedges]; k++)
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
        for (int i = 0; i < dataHCurl.fHDiv.fNumConnectShape[TSHAPE::NSides - TSHAPE::NCornerNodes - 1]; i++)
        {
            for (auto d = 0; d < dim; d++)
            {
                phi(d, count) = curlPhiAux(d, countKernel);
            }
            countKernel++;
            count++;
        }
        // std::cout << "VecDiv = " << vecDiv << std::endl;
        // std::cout << "divphi = " << divphi << std::endl;
        // std::cout << "phi = " << phi << std::endl;
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

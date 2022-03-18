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
template<class TSHAPE>
void TPZShapeHDivConstant<TSHAPE>::Initialize(TPZVec<int64_t> &ids,
                TPZVec<int> &connectorders,
                const TPZVec<int>& sideorient, 
                TPZShapeData &data)
{
    data.fSideOrient = sideorient;
    data.fHDivConnectOrders.Resize(TSHAPE::NFacets+1);
    data.fHDivNumConnectShape.Resize(TSHAPE::NFacets+1);
    if(TSHAPE::Dimension == 2)
    {
        //Compatibilize the polynomial order
        int conSize = connectorders.size();
        for (int i = 0; i < conSize; i++){
            connectorders[i]++;
        }
        if (TSHAPE::Type() == ETriangle){
            connectorders[conSize-1]++;
        }
        //Initialize data structures
        TPZShapeH1<TSHAPE>::Initialize(ids,connectorders,data);

        for(int ic = 0; ic<TSHAPE::NFacets;ic++)
        {
            data.fHDivConnectOrders[ic] =  connectorders[ic];
            data.fHDivNumConnectShape[ic] = data.fH1NumConnectShape[ic]+1;
        }
        int ic = TSHAPE::NFacets;
        data.fHDivConnectOrders[ic] =  connectorders[ic];
        data.fHDivNumConnectShape[ic] = data.fH1NumConnectShape[ic];
    }
    else
    {
        //Compatibilize the polynomial order
        if (TSHAPE::Type() == ETetraedro) {
            int conSize = connectorders.size();
            for (int i = 0; i < conSize; i++){
                connectorders[i]++;
            }
            connectorders[conSize-1]++;
        }

        const int nedges = TSHAPE::NSides-TSHAPE::NFacets-TSHAPE::NCornerNodes-1;
        TPZManVector<int,15> locconorders(TSHAPE::NSides-TSHAPE::NCornerNodes,1);
        for(int ic = nedges; ic<TSHAPE::NSides-TSHAPE::NCornerNodes; ic++)
        {
            locconorders[ic] = connectorders[ic-nedges];
        }
        TPZShapeHCurlNoGrads<TSHAPE>::Initialize(ids, locconorders, data);
//        TPZManVector<int> HCurlNumConnectShapeF = data.fHDivNumConnectShape;
//        data.fHDivNumConnectShape.Resize(TSHAPE::NFacets+1);
//        for(int ic=nedges; ic<TSHAPE::NSides-TSHAPE::NCornerNodes-1; ic++)
//        {
//            int numshape = HCurlNumConnectShapeF[ic]+1;
//            data.fHDivNumConnectShape[ic-nedges] = numshape;
//        }
//        int ic = TSHAPE::NSides-TSHAPE::NCornerNodes-1;
//        int numshape = HCurlNumConnectShapeF[ic];
//        data.fHDivNumConnectShape[ic-nedges] = numshape;
    }
}


template<class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NHDivShapeF(TPZShapeData &data)
{
    // int nshape = TPZShapeH1<TSHAPE>::NShape(data);
    int nshape = 0;
    constexpr int firstConnect = TSHAPE::NSides-TSHAPE::NCornerNodes-TSHAPE::NFacets-1;
    int nc = data.fHDivNumConnectShape.size();
    for(int ic = firstConnect; ic<nc; ic++) nshape += data.fHDivNumConnectShape[ic];
    nshape += TSHAPE::NFacets;
    return nshape;
}

    

template<class TSHAPE>
void TPZShapeHDivConstant<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{

    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    const int nfacets = TSHAPE::NFacets;

    // Compute constant Hdiv functions
    TPZFNMatrix<12,REAL> vecDiv(dim,nfacets);
    TPZVec<REAL> div(nfacets);
    vecDiv.Zero();
    div.Fill(0.);
    // std::cout << "FSide trans ID = " << data.fSideTransformationId << std::endl;
    TSHAPE::ComputeConstantHDiv(pt, vecDiv, div);

    int nshape = data.fPhi.Rows();
    
    if (dim == 2){
        TPZShapeH1<TSHAPE>::Shape(pt,data);    
        divphi.Zero();
        const auto nEdges = TSHAPE::NumSides(1);

        int count = 0;
        int countKernel=ncorner;
        //Edge functions
        for (int i = 0; i < nEdges; i++)
        {
            //RT0 Function
            phi(0,count) = vecDiv(0,i) * data.fSideOrient[i];
            phi(1,count) = vecDiv(1,i) * data.fSideOrient[i];
            divphi(count,0) = div[i] * data.fSideOrient[i];
            count++;

            //Kernel Hdiv
            for (int j = 1; j < data.fHDivConnectOrders[i]; j++)
            {
                phi(0,count) = -data.fDPhi(1,countKernel);
                phi(1,count) =  data.fDPhi(0,countKernel);
                count++;
                countKernel++;
            }
        }
        
        //Internal functions
        for (int i = countKernel; i < nshape; i++){
            phi(0,count) = -data.fDPhi(1,countKernel);
            phi(1,count) =  data.fDPhi(0,countKernel);
            count++;
            countKernel++;
        }
    } else if (dim == 3){
        
        divphi.Zero();
        const auto nEdges = TSHAPE::NumSides(1);
        int nshapehcurl = TPZShapeHCurlNoGrads<TSHAPE>::NHCurlShapeF(data);
        int nshape = NHDivShapeF(data);
        
        TPZFMatrix<REAL> phiAux(dim,nshapehcurl),curlPhiAux(3,nshapehcurl);
        phiAux.Zero(); curlPhiAux.Zero();

        TPZShapeHCurlNoGrads<TSHAPE>::Shape(pt,data,phiAux,curlPhiAux);
        
        int count = 0;
        int countKernel=nEdges;

        //Face functions
        for (int i = 0; i < nfacets; i++)
        {
            // std::cout << "Side orient - " << i << " " << data.fSideOrient[i] << std::endl;
            //RT0 Function
            for(auto d = 0; d < dim; d++) {
                phi(d,count) = vecDiv(d,i) * data.fSideOrient[i];
            }
            divphi(count,0) = div[i] * data.fSideOrient[i];
            count++;
        
            //Kernel HDiv functions
            for (int k = 0; k < data.fHDivNumConnectShape[nEdges]; k++){
                for(auto d = 0; d < dim; d++) {
                    phi(d,count) = curlPhiAux(d,countKernel);               
                }
                countKernel++;
                count++;
            }
        }
        //Internal Functions - HDivKernel
        for (int i = 0; i < data.fHDivNumConnectShape[TSHAPE::NSides-TSHAPE::NCornerNodes-1]; i++)
        {
            for (auto d = 0; d < dim; d++)
            {
                phi(d,count) = curlPhiAux(d,countKernel);  
            }
            countKernel++;
            count++;
        }
        // std::cout << "VecDiv = " << vecDiv << std::endl;
        // std::cout << "divphi = " << divphi << std::endl;
        // std::cout << "phi = " << phi << std::endl;
    } else {
        DebugStop();
    }
    

}

// icon is the connect index of the hdiv element
template<class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NConnectShapeF(int icon, TPZShapeData &data)
{
    const int firstcon = TSHAPE::NSides-TSHAPE::NFacets-TSHAPE::NCornerNodes-1;
    int faceconnect = icon+firstcon;
    int nshape = data.fHDivNumConnectShape[faceconnect] + 1;
    return nshape;
}

template<class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(int connect, int order)
{

#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif
    // int order = data.fHDivConnectOrders[connect];
    MElementType thistype = TSHAPE::Type();

    if(thistype == EOned)
    {   
        order++;
        if(connect < 2) return 0;
        else return order;
        // DebugStop();
    }
    else if(thistype == ETriangle)
    {   
        order++;
        if(connect < TSHAPE::NFacets) return (order);
        else {
            order++;
            return (order-1)*(order-2)/2;
        }
    }
    else if(thistype == EQuadrilateral)
    {
        order++;
        if(connect < TSHAPE::NFacets) return (order);
        else return (order-1)*(order-1);
    }
    else if(thistype == ETetraedro)
    {
        order++;
        if(connect < TSHAPE::NFacets) return 1 + (order-1)*(2*order+4)/4;
        else {
            order++;
            return (order-1)*(order-2)*(2*order+3)/6;
        }
    }
    else if(thistype == EPrisma)
    {
        DebugStop();
        // if(connect == 0 || connect == 4) return (order+1)*(order+2)/2;
        // else if(connect < TSHAPE::NFacets) return (order+1)*(order+1);
        // else return order*order*(3*order+5)/2+7*order-2;
    }
    else if(thistype == ECube)
    {
        if(connect < TSHAPE::NFacets) return 1 + order*(order+2);
        else return order*order*(3+2*order);
    }
    DebugStop();
    unreachable();
 }

template
struct TPZShapeHDivConstant<pzshape::TPZShapeLinear>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapeQuad>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapeTetra>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapeCube>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapePrism>;

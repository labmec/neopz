#include "TPZShapeHDivConstant.h"
#include "TPZShapeHDivKernel.h" 

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


template<class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NHDivShapeF(TPZShapeData &data)
{
    // int nshape = TPZShapeH1<TSHAPE>::NShape(data);
    int nshape = 0;
    int nc = data.fHDivNumConnectShape.size();
    for(int ic = 0; ic<nc; ic++) nshape += data.fHDivNumConnectShape[ic];
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
    TSHAPE::ComputeConstantHDiv(pt, vecDiv, div, data.fSideTransformationId);

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
            for (int j = 1; j < data.fHDivConnectOrders[0]; j++)
            {
                phi(0,count) =  data.fDPhi(1,countKernel);
                phi(1,count) = -data.fDPhi(0,countKernel);
                count++;
                countKernel++;
            }
        }
        
        //Internal functions
        for (int i = countKernel; i < nshape; i++){
            phi(0,count) =  data.fDPhi(1,countKernel);
            phi(1,count) = -data.fDPhi(0,countKernel);
            count++;
            countKernel++;
        }
    } else if (dim == 3){
        
        divphi.Zero();
        const auto nEdges = TSHAPE::NumSides(1);
        int nshape = TPZShapeHDivKernel<TSHAPE>::NHDivShapeF(data);
        
        TPZFMatrix<REAL> phiAux(dim,nshape),divphiAux(nshape,1);
        phiAux.Zero(); divphiAux.Zero();

        TPZShapeHDivKernel<TSHAPE>::Shape(pt,data,phiAux,divphiAux);
        
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
            for (int k = 0; k < data.fHDivNumConnectShape[TSHAPE::NumSides(1)]; k++){
                for(auto d = 0; d < dim; d++) {
                    phi(d,count) = phiAux(d,countKernel);               
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
                phi(d,count) = phiAux(d,countKernel);  
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


template<class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NConnectShapeF(int icon, TPZShapeData &data)
{
    int nshape = data.fHDivNumConnectShape[icon] + 1;
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
        if(connect < 2) return 0;
        else return order;
        // DebugStop();
    }
    else if(thistype == ETriangle)
    {
        if(connect < TSHAPE::NFacets) return (order);
        else return (order-1)*(order-2)/2;
    }
    else if(thistype == EQuadrilateral)
    {
        if(connect < TSHAPE::NFacets) return (order);
        else return (order-1)*(order-1);
    }
    else if(thistype == ETetraedro)
    {
        if(connect < TSHAPE::NFacets) return 1 + (order-1)*(2*order+4)/4;
        else return (order-1)*(order-2)*(2*order+3)/6;
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

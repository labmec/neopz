#include "TPZShapeHDivConstantBound.h"
#include "TPZShapeHDiv.h"
#include "pzgenericshape.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "TPZShapeHDivKernel.h"

template <class TSHAPE>
void TPZShapeHDivConstantBound<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
                                           int connectorder, int sideorient,
                                           TPZShapeData &data) {
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
#ifdef PZDEBUG
    if(ids.size() != ncorner || connectorder < 0 || std::abs(sideorient) != 1)
    {
        DebugStop();
    }
#endif
    data.fSideOrient.Resize(1);
    data.fSideOrient[0] = sideorient;
    data.fH1ConnectOrders.Resize(nsides-ncorner, connectorder);
    data.fH1ConnectOrders.Fill(connectorder);
    if(connectorder > 0)
    {
        TPZShapeH1<TSHAPE>::Initialize(ids,data.fH1ConnectOrders,data);
    }
}

template <class TSHAPE>
int TPZShapeHDivConstantBound<TSHAPE>::NShape(const TPZShapeData &data)
{
    if(data.fH1ConnectOrders[0] == 0) return 1;
    return data.fPhi.Rows();
}

template <>
void TPZShapeHDivConstantBound<class pzshape::TPZShapePoint>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi) {
    phi(0,0) = data.fSideOrient[0];
}


template <class TSHAPE>
void TPZShapeHDivConstantBound<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi) {

    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    const auto nEdges = TSHAPE::NumSides(1);
    TPZShapeH1<TSHAPE>::Shape(pt,data);

    if (dim == 1){
        phi(0,0) = -0.5;

        int nshape = data.fPhi.Rows();

        for (int i = 0; i < nshape-ncorner; i++){
            phi(i+1,0) = -data.fDPhi(0,i+ncorner);
        }
    } else if (dim == 2) {
        
        REAL sign = 1.;
        
        if (TSHAPE::Type() == EQuadrilateral){
            sign = -1.;
            phi(0,0) = -0.25;
        } else if (TSHAPE::Type() == ETriangle){
            phi(0,0) = -2.;
        }

        auto nshape = TPZShapeHDivKernel<TSHAPE>::NHCurlShapeF(data);
        
        TPZFMatrix<REAL> phiAux(2,nshape),divphiAux(nshape,1);
        auto qsi = pt;
        TPZShapeHDivKernel<TSHAPE>::Shape(qsi,data,phiAux,divphiAux);

        auto nEdges = TSHAPE::NumSides(1); 
        // std::cout << "phiaux = " << phiAux << std::endl;
       
        for (int i = nEdges; i < nshape; i++)
        {
            phi(i-nEdges+1,0) = phiAux(0,i) * sign;
            phi(i-nEdges+1,1) = phiAux(1,i) * sign;
        }
        
    } else {
        DebugStop();
    }


}

template<class TSHAPE>
int TPZShapeHDivConstantBound<TSHAPE>::ComputeNConnectShapeF(int connect, int order)
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
        return order;
    }
    else if(thistype == ETriangle)
    {
        if(connect < TSHAPE::NFacets) return 1 + (order-1)*(2*order+4)/4;
        else return 0;
    }
    else if(thistype == EQuadrilateral)
    {
        if(connect < TSHAPE::NFacets) return 1 + order*(order+2);
        else return 0;
    }

    DebugStop();
    unreachable();
 }


template
struct TPZShapeHDivConstantBound<pzshape::TPZShapePoint>;

template
struct TPZShapeHDivConstantBound<pzshape::TPZShapeLinear>;

template
struct TPZShapeHDivConstantBound<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivConstantBound<pzshape::TPZShapeQuad>;


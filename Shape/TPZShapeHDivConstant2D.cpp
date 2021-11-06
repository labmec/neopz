#include "TPZShapeHDivConstant2D.h"

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


template<class TSHAPE>
int TPZShapeHDivConstant2D<TSHAPE>::NHDivShapeF(TPZShapeData &data)
{
    // int nshape = TPZShapeH1<TSHAPE>::NShape(data);
    int nshape = 0;
    int nc = data.fHDivNumConnectShape.size();
    for(int ic = 0; ic<nc; ic++) nshape += data.fHDivNumConnectShape[ic];
    return nshape;
}

    

template<class TSHAPE>
void TPZShapeHDivConstant2D<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{

    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    const int nfacets = TSHAPE::NFacets;

    TPZShapeH1<TSHAPE>::Shape(pt,data);
    divphi.Zero();
    const auto nEdges = TSHAPE::NumSides(1);

    if (dim != 2)
    {
        DebugStop();
    }

    // Compute constant Hdiv functions
    TPZFNMatrix<12,REAL> vecDiv(dim,nfacets);
    TPZVec<REAL> div(nfacets);
    vecDiv.Zero();
    div.Fill(0.);
    TSHAPE::ComputeConstantHDiv(pt, vecDiv, div, data.fSideTransformationId);

    int nshape = data.fPhi.Rows();
    phi.Resize(2,nshape);
    //Div free functions
    for (int i = ncorner; i < nshape; i++){
		phi(0,i) =  data.fDPhi(1,i);
        phi(1,i) = -data.fDPhi(0,i);
	}

    //Div-Constant Functions
    for (int i = 0; i < ncorner; i++){
		phi(0,i) = vecDiv(0,i);
        phi(1,i) = vecDiv(1,i);
	}

    divphi.Zero();
}


template<class TSHAPE>
int TPZShapeHDivConstant2D<TSHAPE>::NConnectShapeF(int icon, TPZShapeData &data)
{
    int order = data.fH1ConnectOrders[icon];
    const int side = icon + TSHAPE::NCornerNodes;
    #ifdef PZDEBUG
    if (side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) {
        DebugStop();
    }
    #endif
    if(order == 0) {
        PZError<<__PRETTY_FUNCTION__
            <<"\nERROR: polynomial order not compatible.\nAborting..."
            <<std::endl;
        DebugStop();
        return 0;
    }
    const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    const int nShapeF = [&](){
        if (side < TSHAPE::NCornerNodes + nEdges) {//edge connect
        return 1;
        }
        else if(side < TSHAPE::NCornerNodes + nEdges + nFaces){//face connect
        switch(TSHAPE::Type(side)){
            case ETriangle://triangular face
                /**
                 we remove one internal function for each h1 face function of order k+1
                since there are (k-1)(k-2)/2 functions per face in a face with order k,
                we remove k(k-1)/2.
                so:
                (k-1)*(k+1)-k*(k-1)/2
                */
                return (order - 1) * (order+2) / 2;
            default:
                PZError<<__PRETTY_FUNCTION__<<" error. Not yet implemented"<<std::endl;
                DebugStop();
                return 0;
            }
        }
        else{//internal connect (3D element only)
        if constexpr (TSHAPE::Type() == ETetraedro){
            DebugStop();
        }
        return 0;
        }
    }();
    return nShapeF;
    // return data.fHDivNumConnectShape[icon];
}



template
struct TPZShapeHDivConstant2D<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivConstant2D<pzshape::TPZShapeQuad>;


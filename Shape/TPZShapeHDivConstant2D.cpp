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
    // phi.Resize(2,nshape);

    for (int i = 0; i < nshape; i++){
		phi(0,i) =  data.fDPhi(1,i);
        phi(1,i) = -data.fDPhi(0,i);
        divphi(i,0) = div[i];
	}
    

    

    //Div-Constant Functions
    // for (int i = 0; i < ncorner; i++){
	// 	phi(0,i) = vecDiv(0,i);
    //     phi(1,i) = vecDiv(1,i);
	// }

    // divphi.Zero();
}


template<class TSHAPE>
int TPZShapeHDivConstant2D<TSHAPE>::NConnectShapeF(int icon, TPZShapeData &data)
{
    int nshape = data.fHDivNumConnectShape[icon] + 1;
    return nshape;
}

template<class TSHAPE>
int TPZShapeHDivConstant2D<TSHAPE>::ComputeNConnectShapeF(int connect, int order)
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
        if(connect < TSHAPE::NFacets) return (order+1);
        else return (order+1)*(order+1)-1;
    }
    else if(thistype == EQuadrilateral)
    {
        if(connect < TSHAPE::NFacets) return (order);
        else return (order-1)*(order-1);// 2*order*(order+1);
    }
    else if(thistype == ETetraedro)
    {
        if(connect < TSHAPE::NFacets) return (order+1)*(order+2)/2;
        else return order*(order+2)*(order+3)/2;
    }
    else if(thistype == EPrisma)
    {
        if(connect == 0 || connect == 4) return (order+1)*(order+2)/2;
        else if(connect < TSHAPE::NFacets) return (order+1)*(order+1);
        else return order*order*(3*order+5)/2+7*order-2;
    }
    else if(thistype == ECube)
    {
        if(connect < TSHAPE::NFacets) return (order+1)*(order+1);
        else return 3*order*(order+1)*(order+1);
    }
    DebugStop();
    unreachable();
 }


template
struct TPZShapeHDivConstant2D<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivConstant2D<pzshape::TPZShapeQuad>;


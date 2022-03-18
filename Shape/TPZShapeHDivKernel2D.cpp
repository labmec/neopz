#include "TPZShapeHDivKernel2D.h"

#include "TPZShapeH1.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "TPZShapeData.h"


template<class TSHAPE>
int TPZShapeHDivKernel2D<TSHAPE>::NHDivShapeF(TPZShapeData &data)
{
    int nshape = 0;
    int nc = data.fHDivNumConnectShape.size();
    for(int ic = 0; ic<nc; ic++) nshape += data.fHDivNumConnectShape[ic];
    return nshape;
}

    

template<class TSHAPE>
void TPZShapeHDivKernel2D<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{

    const int dim = TSHAPE::Dimension;
    TPZShapeH1<TSHAPE>::Shape(pt,data);
    divphi.Zero();

    int nshape = data.fPhi.Rows();
    phi.Resize(2,nshape);

    switch (dim)
    {
    case 1:
        for (int i = 0; i < nshape; i++){
            phi(0,i) = data.fDPhi(0,i);
        }
        break;
    
    case 2:
        for (int i = 0; i < nshape; i++){
            phi(0,i) =  data.fDPhi(1,i);
            phi(1,i) = -data.fDPhi(0,i);
        }
        break;
    
    default:
        DebugStop();
        break;
    }


}

template
struct TPZShapeHDivKernel2D<pzshape::TPZShapeLinear>;

template
struct TPZShapeHDivKernel2D<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivKernel2D<pzshape::TPZShapeQuad>;


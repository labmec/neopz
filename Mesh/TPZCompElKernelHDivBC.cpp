/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHDivBC.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"
#include "TPZShapeHDivKernel2DBound.h"

template<class TSHAPE>
TPZCompElKernelHDivBC<TSHAPE>::TPZCompElKernelHDivBC(TPZCompMesh &mesh, TPZGeoEl *gel) :
TPZRegisterClassId(&TPZCompElKernelHDivBC::ClassId), fSideOrient(1), TPZCompElH1<TSHAPE>(mesh,gel) {

}

template<class TSHAPE>
TPZCompElKernelHDivBC<TSHAPE>::~TPZCompElKernelHDivBC(){
    this->~TPZCompElH1<TSHAPE>();
}

template<class TSHAPE>
void TPZCompElKernelHDivBC<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZCompElH1<TSHAPE>::InitMaterialData(data);

    int nshape = this->NShapeF();
    data.fVecShapeIndex.Resize(nshape);
    for (int i=0; i<nshape; i++) {
		data.fVecShapeIndex[i] = std::make_pair(i,1);
    }

}

template<class TSHAPE>
void TPZCompElKernelHDivBC<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    data.fNeedsSol = needsol;
    TPZCompElH1<TSHAPE>::ComputeShape(qsi,data);

    TPZShapeData &shapedata = data;
    int nshape = this->NShapeF();
    // data.phi.Resize(nshape,1);
    // data.dphix.Resize(1,nshape);

    TPZFMatrix<REAL> auxPhi(1,nshape);
    auxPhi.Zero();
    TPZShapeHDivKernel2DBound<TSHAPE>::Shape(qsi, shapedata, auxPhi, data.divphi);
    if (data.phi.Rows()>1){
      for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = -auxPhi(0,i)/data.detjac;
	  }
    }
    // std::cout << "AUXphi = " << data.phi << std::endl;


}//void


// NAO TESTADO
template<class TSHAPE>
void TPZCompElKernelHDivBC<TSHAPE>::SetSideOrient(int side, int sideorient)
{
    if (side != TSHAPE::NSides - 1) {
        DebugStop();
    }
    fSideOrient = sideorient;
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElKernelHDivBC<TSHAPE>::GetSideOrient(int side)
{
    if (side != TSHAPE::NSides - 1) {
        DebugStop();
    }
    return fSideOrient;
}

#include "pzshapelinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelt2dmapped.h"

using namespace pztopology;

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;

// template class TPZCompElKernelHDivBC<TPZShapePoint>;
template class TPZCompElKernelHDivBC<TPZShapeLinear>;
// template class TPZCompElKernelHDivBC<TPZShapeTriang>;
// template class TPZCompElKernelHDivBC<TPZShapeQuad>;

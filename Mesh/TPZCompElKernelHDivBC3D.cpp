/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHDivBC3D.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"
#include "TPZShapeHDivKernel.h"
#include "TPZShapeHDivConstant.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

template<class TSHAPE>
TPZCompElKernelHDivBC3D<TSHAPE>::TPZCompElKernelHDivBC3D(TPZCompMesh &mesh, TPZGeoEl *gel, int shapetype) :
TPZRegisterClassId(&TPZCompElKernelHDivBC3D::ClassId), TPZCompElHCurlNoGrads<TSHAPE>(mesh,gel), fShapeType(shapetype)  {

}

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
    if (fShapeType == ECurlNoGrads){
	    TPZCompElHCurlNoGrads<TSHAPE>::InitMaterialData(data);
    } else if (fShapeType == EHDivKernel) {
        //Init the material data of Hcurl
        TPZCompElHCurl<TSHAPE>::InitMaterialData(data);
        
        TPZShapeData dataaux = data;
        data.fVecShapeIndex=dataaux.fSDVecShapeIndex;
        data.divphi.Resize(data.fVecShapeIndex.size(),1);
        TPZShapeHDivKernel<TSHAPE>::ComputeVecandShape(data);
        
        //setting the type of shape functions as vector shape functions
        data.fShapeType = TPZMaterialData::EVecShape;
    } else {
        DebugStop();
    }


}

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::ComputeRequiredDataT(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi){

    if (fShapeType == ECurlNoGrads){
        bool needsol = data.fNeedsSol;
        data.fNeedsSol = true;
        TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredData(data,qsi);
        data.fNeedsSol = needsol;
    } else {
        //Compute the element geometric data
        TPZGeoEl * ref = this->Reference();
        if (!ref){
            PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
            return;
        }

        ref->Jacobian(qsi, data.jacobian, data.axes, data.detjac , data.jacinv);
        ref->X(qsi, data.x);
        data.xParametric = qsi;
        
        constexpr auto dim{TSHAPE::Dimension};
        auto nshape = TPZShapeHDivKernel<TSHAPE>::NHCurlShapeF(data);
        // auto nshape = TPZShapeHDivConstant<TSHAPE>::NHDivShapeF(data);
        
        TPZFMatrix<REAL> phiAux(dim,nshape),divphiAux(nshape,1);

        TPZShapeHDivKernel<TSHAPE>::Shape(qsi,data,phiAux,divphiAux);

        TPZCompElHCurl<TSHAPE>::TransformCurl(phiAux, data.detjac, data.jacobian, data.curlphi);

        const int ncorner = TSHAPE::NCornerNodes;
        int nEdges = TSHAPE::NumSides(1);
        const int nsides = TSHAPE::NSides;

        // std::cout << "dphi" << data.fDPhi << std::endl;
        // std::cout << "master" << data.fMasterDirections << std::endl;
    
        // std::cout << "data.curlphi = " << data.curlphi << std::endl;
        // std::cout << "data.phi = " << data.phi << std::endl;

        data.divphi = divphiAux;
        // if (fSideOrient != 1) DebugStop();
        // if (data.fSideOrient[0] != 1) DebugStop();
        // data.phi = phiHCurl;
        if (data.fNeedsSol) {
            this->ReallyComputeSolution(data);
        }
    }
    
    data.phi.Resize(data.curlphi.Cols(),3);
    data.phi.Zero();
    if (data.phi.Rows()>1){
      for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = -data.curlphi(0,i);
	  }
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        //	this->Print(sout);
        // sout << "\nVecshape = " << data.fVecShapeIndex << std::endl;
        // sout << "MASTER = " << data.fMasterDirections << std::endl;
        // sout << "TransformIDs = " << data.fSideTransformationId << std::endl;
        sout << "Phi = " << data.phi << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
        
    }
#endif

}//void

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::SetSideOrient(int orient){
    fSideOrient = orient;
}

template<class TSHAPE>
int TPZCompElKernelHDivBC3D<TSHAPE>::GetSideOrient(){
    return fSideOrient;
}

// template<class TSHAPE>
// int TPZCompElKernelHDivBC3D<TSHAPE>::NConnectShapeF(int icon, int order) const{
  
//         const int side = icon;// + TSHAPE::NCornerNodes;
//     // #ifdef PZDEBUG
//     // if (side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) {
//     //     DebugStop();
//     // }
//     // #endif
//     if(order == 0) {
//         PZError<<__PRETTY_FUNCTION__
//             <<"\nERROR: polynomial order not compatible.\nAborting..."
//             <<std::endl;
//         DebugStop();
//         return 0;
//     }
//     const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
//     const auto nEdges = TSHAPE::NumSides(1);
//     const int nShapeF = [&](){
//         if (side < TSHAPE::NCornerNodes + nEdges) {//edge connect
//         return 1;
//         }
//         else if(side < TSHAPE::NCornerNodes + nEdges + nFaces){//face connect
//         switch(TSHAPE::Type(side)){
//         case ETriangle://triangular face
//             /**
//              we remove one internal function for each h1 face function of order k+1
//             since there are (k-1)(k-2)/2 functions per face in a face with order k,
//             we remove k(k-1)/2.
//             so:
//             (k-1)*(k+1)-k*(k-1)/2
//             */
//             return (order - 1) * (order+2) / 2;
//         default:
//             PZError<<__PRETTY_FUNCTION__<<" error. Not yet implemented"<<std::endl;
//             DebugStop();
//             return 0;
//         }
//         }
//         else{//internal connect (3D element only)
//         if constexpr (TSHAPE::Type() == ETetraedro){
//             return (order-1)*(order-2)*(2*order+3)/6;
//         }
//         return 0;
//         }
//     }();
//     return nShapeF;


// }



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

// #include "tpzpoint.h"
// #include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"

// #include "TPZCompElHCurl.h"

using namespace pzgeom;
using namespace pzshape;

template class TPZCompElKernelHDivBC3D<TPZShapeTriang>;
template class TPZCompElKernelHDivBC3D<TPZShapeQuad>;

/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHDivBC3D.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"
#include "TPZShapeHDivKernel.h"
#include "TPZShapeHDivConstant.h"
#include "TPZShapeHCurlNoGrads.h"
#include <pzconnect.h>
#include <pzcmesh.h>

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

template<class TSHAPE>
TPZCompElKernelHDivBC3D<TSHAPE>::TPZCompElKernelHDivBC3D(TPZCompMesh &mesh, TPZGeoEl *gel, 
                                                         const HDivFamily hdivfam, const HCurlFamily hcurlfam) :
TPZRegisterClassId(&TPZCompElKernelHDivBC3D::ClassId), 
TPZCompElHCurl<TSHAPE>(mesh,gel), fhdivfam(hdivfam), fhcurlfam(hcurlfam)  {
    if (fhdivfam != HDivFamily::EHDivKernel && fhcurlfam != HCurlFamily::EHCurlNoGrads){
        std::cout << "You need to chose EHivKernel or EHCurlNoGrads approximation space to use TPZCompElKernelHDivBC3D" << std::endl;
        DebugStop();
    }
    this->AdjustConnects();
}

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
    if (fhcurlfam == HCurlFamily::EHCurlNoGrads){
	    //Init the material data of Hcurl
        TPZCompElHCurl<TSHAPE>::InitMaterialData(data);
        /*
            we first compute ALL the hcurl traditional functions(scalar+vector)
            then, at each integration point we will filter them.
            
            some of this code is copied from the TPZCompElHCurlFull since we want to avoid
            calls for the virtual methods.
        */
        
        //computes the index that will associate each scalar function to a constant vector field
        constexpr auto nConnects = TSHAPE::NSides - TSHAPE::NCornerNodes;
        TPZManVector<int,nConnects> connectOrders(nConnects,-1);
        int unfiltnshape = 0;
        for(auto i = 0; i < nConnects; i++){
            const auto conorder = this->EffectiveSideOrder(i + TSHAPE::NCornerNodes);
            connectOrders[i] = conorder;
            unfiltnshape += TPZCompElHCurl<TSHAPE>::NConnectShapeF(i, conorder);
        }

        const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
        const auto nEdges = TSHAPE::NumSides(1);
        constexpr auto nNodes = TSHAPE::NCornerNodes;

        TPZManVector<int64_t,nNodes> nodes(nNodes, 0);

        for (auto iNode = 0; iNode < nNodes; iNode++){
            nodes[iNode] = this->Reference()->NodeIndex(iNode);
        }
    
        TPZManVector<int64_t, TSHAPE::NSides - nNodes>
            firstH1ShapeFunc(TSHAPE::NSides - nNodes,0);
        
        //calculates the first shape function associated with each side of dim > 0
        TPZManVector<int,TSHAPE::NSides-nNodes> sidesH1Ord(TSHAPE::NSides - nNodes,-1);
        // TPZShapeHDivKernel<TSHAPE>::CalcH1ShapeOrders(connectOrders,sidesH1Ord);
        TPZShapeHCurl<TSHAPE>::CalcH1ShapeOrders(connectOrders,sidesH1Ord);
        firstH1ShapeFunc[0] = nNodes;
        
        for (int iSide = nNodes + 1; iSide < TSHAPE::NSides; iSide++) {
            const int iCon = iSide - nNodes;
            firstH1ShapeFunc[iCon] =
            firstH1ShapeFunc[iCon - 1] +
            TSHAPE::NConnectShapeF(iSide - 1, sidesH1Ord[iCon-1]);

        }
        
        auto &indexVecShape = data.fVecShapeIndex;
        indexVecShape.Resize(unfiltnshape);

        TPZVec<unsigned int> shapeCountVec(TSHAPE::NSides - nNodes, 0);
        TPZShapeHCurl<TSHAPE>::StaticIndexShapeToVec(connectOrders,
                                                    firstH1ShapeFunc,sidesH1Ord, nodes, shapeCountVec,indexVecShape );

        //setting the type of shape functions as vector shape functions
        data.fShapeType = TPZMaterialData::EVecShape;

    } else if (fhdivfam == HDivFamily::EHDivKernel) {
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

    if (fhcurlfam == HCurlFamily::EHCurlNoGrads){
        bool needsol = data.fNeedsSol;
        data.fNeedsSol = true;
        // We can't use ComputeRequired data from TPZCompElHCurl, because we are dealing with a number of shape functions
        // uncompatible with the number of HCurl shape functions
    
        /******************************************************************************************************************
         * at this point, we already have the basis functions on the deformed element, since we have data.phi,
         * data.fVecShapeIndex and data.fDeformedDirections. Now it is time to compute the curl, which will be stored in
         * data.curlphi.
         *******************************************************************************************************************/
        
        //Compute the element geometric data
        TPZGeoEl * ref = this->Reference();
        if (!ref){
            PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
            return;
        }
        ref->Jacobian(qsi, data.jacobian, data.axes, data.detjac , data.jacinv);

        ref->X(qsi, data.x);
        data.xParametric = qsi;
        
        const int rows = data.phi.Rows();
        const int cols = data.phi.Cols();
        TPZFMatrix<REAL> phiHCurl(cols,rows,0.);
        TPZFMatrix<REAL> curlHCurl(cols,rows,0.);

        TPZShapeHCurlNoGrads<TSHAPE>::Shape(qsi,data,phiHCurl,curlHCurl);

        TPZCompElHCurl<TSHAPE>::TransformCurl(curlHCurl, data.detjac, data.jacobian, data.curlphi);

        //data.phi.Transpose(&phiHCurl);

        if (data.fNeedsSol) {
            this->ReallyComputeSolution(data);
        };
        
    } else if (fhdivfam == HDivFamily::EHDivKernel) {
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
    } else {
        DebugStop();
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
        // sout << "Phi = " << data.phi << std::endl;
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

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::AdjustConnects()
{
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto ncon = TSHAPE::NSides - nNodes;
    for(int icon = 0; icon < ncon; icon++){
        const int connect = this->MidSideConnectLocId(icon+nNodes);
        TPZConnect &c = this->Connect(connect);
        const int nshape =this->NConnectShapeF(connect,c.Order());
        c.SetNShape(nshape);
        const auto seqnum = c.SequenceNumber();
        const int nStateVars = [&](){
            TPZMaterial * mat =this-> Material();
            if(mat) return mat->NStateVariables();
            else {
                return 1;
            }
        }();
        this-> Mesh()->Block().Set(seqnum,nshape*nStateVars);
    }
}


template<class TSHAPE>
int TPZCompElKernelHDivBC3D<TSHAPE>::NConnectShapeF(int icon, int order) const
{
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
        case EQuadrilateral://quadrilateral face
            //Following the same logic:
            /**
             we remove one internal function for each h1 face function of order k+1
            since there are (k-1)^2 functions per face in a face with order k,
            we remove k^2.
            so:
            2k(k+1) - k^2 = k(k+2)
            
            */
            // if (order == 2) return 0;
            return order * (order + 2);
            // return 2*order * (order + 1);
        default:
            PZError<<__PRETTY_FUNCTION__<<" error. Not yet implemented"<<std::endl;
            DebugStop();
            return 0;
        }
        }
        else{//internal connect (3D element only)
        if constexpr (TSHAPE::Type() == ETetraedro){
        /**
             we remove one internal function for each h1 internal function of order k+1
            since there are (k-1)(k-2)(k-3)/6 functions in a h1 element with order k,
            we remove k(k-1)(k-2)/6.
            so:
            (k-1)(k-2)(k+1)/2 - k(k-1)(k-2)/6 = (k-1)(k-2)(2k+3)/6.

            since we will remove k(k-1)(k-2)/6, for each     we remove (k-1)(k-2)/2 funcs.

            we have two kinds of internal functions. phi_kf and phi_ki.
            func        k                   k-1                 new funcs
            phi_kf      2(k-1)(k-2)         2(k-2)(k-3)         4(k-2)
            phi_ki      (k-1)(k-2)(k-3)/2   (k-2)(k-3)(k-4)/2   3(k-2)(k-3)/2
            
            that means that if we remove, for each k, (k-2) phi_kf 
            (for instance, all phi_kf associated with a given face),
            we need to remove (k-1)(k-2)/2 - (k-2) = (k-2)(k-3)/2
            which is exactly one third of the phi_ki.
            
        */
            return (order-1)*(order-2)*(2*order+3)/6;
        } else 
        if constexpr (TSHAPE::Type() == ECube){
        /**
             we remove one internal function for each h1 internal function of order k+1
            since there are (k-1)^3 functions in a h1 element with order k,
            we remove k^3.
            so:
            3k^2(k+1) - k^3 = k^2 (3 + 2 k).

            we have two kinds of internal functions. phi_kf and phi_ki.
            func        k                   k-1                 new funcs
            phi_kf      6k^2                6(k-1)^2            k(8k-k^2-1)
            phi_ki      3k^2*(k-1)          3(k-1)(k-1)(k-2)/2  3(2-5k+2k^2+k^3)/2
           
        */
            return order*order*(3+2*order);
        }
        return 0;
        }
    }();
    return nShapeF;
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

// #include "tpzpoint.h"
// #include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"

// #include "TPZCompElHCurl.h"

using namespace pzgeom;
using namespace pzshape;

template class TPZCompElKernelHDivBC3D<TPZShapeTriang>;
template class TPZCompElKernelHDivBC3D<TPZShapeQuad>;

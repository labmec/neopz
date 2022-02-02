/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHDiv3D.h"
#include "pzcmesh.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "TPZMatSingleSpace.h"
#include "pzlog.h"
#include "pzgeoquad.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZMaterialDataT.h"
#include "pzshapepiram.h"
// #include "tpzline.h"
#include "tpztriangle.h"
#include "TPZShapeHDivKernel.h"
#include "TPZShapeHCurlNoGrads.h"
#include <pzconnect.h>
#include <pzcmesh.h>

#include "pzshtmat.h"


#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

using namespace std;

template<class TSHAPE>
TPZCompElKernelHDiv3D<TSHAPE>::TPZCompElKernelHDiv3D(TPZCompMesh &mesh, TPZGeoEl *gel, 
                                                     const HDivFamily hdivfam, const HCurlFamily hcurlfam) :
TPZRegisterClassId(&TPZCompElKernelHDiv3D::ClassId), TPZCompElHCurl<TSHAPE>(mesh,gel),
                  fSideOrient(TSHAPE::NFacets,1), fhdivfam(hdivfam), fhcurlfam(hcurlfam) {
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    for(int side = firstside ; side < TSHAPE::NSides-1; side++ )
    {
        fSideOrient[side-firstside] = this->Reference()->NormalOrientation(side);
    }
    
    if (fhdivfam != HDivFamily::EHDivKernel && fhcurlfam != HCurlFamily::EHCurlNoGrads){
        std::cout << "You need to chose EHivKernel or EHCurlNoGrads approximation space to use TPZCompElKernelHDivBC3D" << std::endl;
        DebugStop();
    }
    this->AdjustConnects();
}

template<class TSHAPE>
TPZCompElKernelHDiv3D<TSHAPE>::~TPZCompElKernelHDiv3D(){
    this->~TPZCompElHCurl<TSHAPE>();
}
 

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv3D<TSHAPE>::ComputeRequiredDataT(TPZMaterialDataT<TVar> &data,
                                                TPZVec<REAL> &qsi){
                                                
    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    if (fhcurlfam == HCurlFamily::EHCurlNoGrads) {
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
        int nshape = 0;
        nshape = TPZShapeHDivKernel<TSHAPE>::NHDivShapeF(data);
        
        TPZFMatrix<REAL> phiAux(dim,nshape),divphiAux(nshape,1);
        phiAux.Zero(); divphiAux.Zero();

        TPZShapeHDivKernel<TSHAPE>::Shape(qsi,data,phiAux,divphiAux);

        TPZCompElHCurl<TSHAPE>::TransformCurl(phiAux, data.detjac, data.jacobian, data.curlphi);

        const int ncorner = TSHAPE::NCornerNodes;
        int nEdges = TSHAPE::NumSides(1);
        const int nsides = TSHAPE::NSides;

        data.divphi = divphiAux;
    
        if (data.fNeedsSol) {
            this->ReallyComputeSolution(data);
        }
    } else {
        DebugStop();
    }
    
    data.fNeedsSol = needsol;

    int nshape = this->NShapeF();
    data.fDeformedDirections.Resize(3,nshape);
    data.fVecShapeIndex.Resize(nshape);
    TPZShapeData &shapedata = data;
    int size = data.curlphi.Cols();

    if (size != nshape) DebugStop();
    int ncorner = TSHAPE::NCornerNodes;
    for (int j = 0; j < nshape; j++){
        data.fVecShapeIndex[j].first = j;
        data.fVecShapeIndex[j].second = j;
        for (int i = 0; i < 3; i++){
            data.fDeformedDirections(i,j)=data.curlphi(i,j);
        }
    }
    
    // data.fDeformedDirections=data.curlphi;

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        //	this->Print(sout);
        // sout << "\nVecshape = " << data.fVecShapeIndex << std::endl;
        // sout << "Phi = " << data.fDeformedDirections << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
        
    }
#endif

    data.phi.Resize(nshape,3);
    REAL val = 1.;

    for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = val;
        data.phi(i,1) = val;
        data.phi(i,2) = val;
	}

	for (int i = 0; i < data.dphix.Rows(); i++)
        for (int j = 0; j < data.dphix.Cols(); j++)
    	    	data.dphix(i,j) = 1.;
    
    // for (int i=0; i<data.fVecShapeIndex.size(); i++) {
	// 	data.fVecShapeIndex[i] = std::make_pair(i,1);
    // }
    if (data.fNeedsSol) {
        this->ReallyComputeSolution(data);
    }

}//void

template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	data.fNeedsSol = true;
	if (fhcurlfam == HCurlFamily::EHCurlNoGrads) {
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

    data.fShapeType = data.EVecandShape;
    
}

/**
 * @brief It returns the normal orientation of the reference element by the side.
 * Only side that has dimension larger than zero and smaller than me.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
int TPZCompElKernelHDiv3D<TSHAPE>::GetSideOrient(int side){

    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    return fSideOrient[side-firstside];
}

/**
 * @brief It set the normal orientation of the element by the side.
 * Only side that has dimension equal to my dimension minus one.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::SetSideOrient(int side, int sideorient){

    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    fSideOrient[side-firstside] = sideorient;
}

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv3D<TSHAPE>::ComputeSolutionKernelHdivT(TPZMaterialDataT<TVar> &data)
{
    this->ComputeSolutionHCurlT(data.phi, data.curlphi,
                                data.sol, data.curlsol);
    // data.fDeformedDirections=data.curlphi;

    const int dim = 3; // Hdiv vectors are always in R3
    const int nstate = this->Material()->NStateVariables();
    const int ncon = this->NConnects();

    TPZFMatrix<TVar> &MeshSol = this->Mesh()->Solution();

    int64_t numbersol = MeshSol.Cols();

    if(numbersol != 1)
    {
        DebugStop();
    }
    data.sol.Resize(numbersol);
    data.dsol.Resize(numbersol);
    data.divsol.Resize(numbersol);

    for (int64_t is=0; is<numbersol; is++)
    {
        data.sol[is].Resize(dim*nstate);
        data.sol[is].Fill(0);
        data.dsol[is].Redim(dim*nstate, dim);
        data.divsol[is].Resize(nstate);
        data.divsol[is].Fill(0.);
    }
    data.sol = data.curlsol;

}

template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::AdjustConnects()
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
int TPZCompElKernelHDiv3D<TSHAPE>::NConnectShapeF(int icon, int order) const
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
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"

using namespace pztopology;

#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"

using namespace pzgeom;
using namespace pzshape;


template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeTetra>>;
template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapePrism>>;

template class TPZCompElKernelHDiv3D<TPZShapeTetra>;
template class TPZCompElKernelHDiv3D<TPZShapeCube>;
template class TPZCompElKernelHDiv3D<TPZShapePrism>;

#include "TPZCompElKernelHDivBC.h"
#include "TPZCompElKernelHDivBC3D.h"
#include "TPZCompElKernelHDiv.h"


//BC
TPZCompEl * CreateHDivKernelBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElH1<TPZShapePoint>(mesh,gel);
}
TPZCompEl * CreateHDivKernelBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDivBC< TPZShapeLinear>(mesh,gel);
}
TPZCompEl * CreateHDivKernelBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam, const HCurlFamily hcurlfam) {
    return new TPZCompElKernelHDivBC3D< TPZShapeQuad>(mesh,gel,hdivfam,hcurlfam);
}
TPZCompEl * CreateHDivKernelBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam, const HCurlFamily hcurlfam) {
    return new TPZCompElKernelHDivBC3D< TPZShapeTriang >(mesh,gel,hdivfam,hcurlfam);
}

//Domain
TPZCompEl * CreateHDivKernelQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv< TPZShapeQuad>(mesh,gel);
}
TPZCompEl * CreateHDivKernelTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv< TPZShapeTriang >(mesh,gel);
}
TPZCompEl * CreateHDivKernelCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam, const HCurlFamily hcurlfam) {
	return new TPZCompElKernelHDiv3D< TPZShapeCube >(mesh,gel,hdivfam,hcurlfam);
}
TPZCompEl * CreateHDivKernelPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam, const HCurlFamily hcurlfam) {
	return new TPZCompElKernelHDiv3D< TPZShapePrism>(mesh,gel,hdivfam,hcurlfam);
}
TPZCompEl * CreateHDivKernelTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam, const HCurlFamily hcurlfam) {
	return new TPZCompElKernelHDiv3D< TPZShapeTetra >(mesh,gel,hdivfam,hcurlfam);
}

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
#include "TPZShapeHCurlNoGrads.h"
#include <pzconnect.h>
#include <pzcmesh.h>

#include "pzshtmat.h"


#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

using namespace std;

template<class TSHAPE>
TPZCompElKernelHDiv3D<TSHAPE>::TPZCompElKernelHDiv3D(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElKernelHDiv3D::ClassId), TPZCompElHCurl<TSHAPE>(mesh,gel,HCurlFamily::EHCurlNoGrads),
                  fhdivfam(hdivfam) {
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    
    //Updates the number of shape functions and also the integration rule
    if (TSHAPE::Type() == ETetraedro || TSHAPE::Type() == ETriangle){
        for (int icon = 0; icon < this->NConnects(); icon++)
        {
            TPZConnect &c = this->Connect(icon);
            int nShapeF = NConnectShapeF(icon,c.Order());
            c.SetNShape(nShapeF);
            int64_t seqnum = c.SequenceNumber();
            int nvar = 1;
            TPZMaterial * mat = this->Material();
            if (mat) nvar = mat->NStateVariables();
            this->Mesh()->Block().Set(seqnum, nvar * nShapeF);
            this->AdjustIntegrationRule();
        }
    }
}
 

template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::InitMaterialData(TPZMaterialData &data){

	TPZIntelGen<TSHAPE>::InitMaterialData(data);
#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
        LOGPZ_DEBUG(logger,"Initializing MaterialData of TPZCompElHCurl")
    }
#endif
    data.fShapeType = TPZMaterialData::MShapeFunctionType::EVecShape;
    TPZShapeData & shapedata = data;

    TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes,0);
    TPZGeoEl *ref = this->Reference();
    for(auto i=0; i<TSHAPE::NCornerNodes; i++) {
        ids[i] = ref->NodePtr(i)->Id();
    }
    
    auto &conOrders = shapedata.fHDivConnectOrders;
    constexpr auto nConnects = TSHAPE::NSides - TSHAPE::NCornerNodes;
    conOrders.Resize(nConnects,-1);
    // For Tetrahedra we increase the polynomial order by one. The internal connect is increased by two
    int ordTop = 0;
    if (TSHAPE::Type() == ETetraedro || TSHAPE::Type() == ETriangle) ordTop = 1;
    for(auto i = 0; i < nConnects; i++){
        conOrders[i] = this->EffectiveSideOrder(i + TSHAPE::NCornerNodes) + ordTop;
    }
    if (TSHAPE::Type() == ETetraedro){
        conOrders[nConnects-1]++;
    }

    TPZShapeHCurlNoGrads<TSHAPE>::Initialize(ids, conOrders, shapedata);

    //resizing of TPZMaterialData structures
    constexpr int dim = TSHAPE::Dimension;
    constexpr int curldim = [dim](){
        if constexpr (dim == 1) return 1;
        else{
            return 2*dim - 3;//1 for 2D 3 for 3D
        }
    }();
    const int nshape = this->NShapeF();
    
    auto &phi = data.phi;
    auto &curlphi = data.curlphi;
    
    phi.Redim(nshape,3);
    curlphi.Redim(curldim,nshape);
    data.divphi.Redim(nshape,1);

    data.axes.Redim(dim,3);
    data.jacobian.Redim(dim,dim);
    data.jacinv.Redim(dim,dim);
    data.x.Resize(3);
    data.fDeformedDirections.Resize(3,nshape);
    data.fVecShapeIndex.Resize(nshape);

}


template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data){
                                                
    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;

    auto nshape = this->NShapeF();
    TPZFMatrix<REAL> phiHCurl(3,nshape,0.);
    TPZFMatrix<REAL> curlHCurl(3,nshape,0.);

    TPZShapeHCurlNoGrads<TSHAPE>::Shape(qsi,data,phiHCurl,curlHCurl);
    TPZCompElHCurl<TSHAPE>::TransformCurl(curlHCurl, data.detjac, data.jacobian, data.curlphi);    
    
    data.fNeedsSol = needsol;
    if (TSHAPE::Dimension == 3){        
        TPZShapeData &shapedata = data;
        int size = data.curlphi.Cols();
        if (size != nshape) DebugStop();
        int ncorner = TSHAPE::NCornerNodes;
        data.fDeformedDirections = data.curlphi;
        data.phi = 1.;
        data.dphix = 1.;
        
    } else if (TSHAPE::Dimension == 2) {
        if (data.phi.Rows()>1){
            for (int i = 0; i < data.phi.Rows(); i++){
                data.phi(i,0) = data.curlphi(0,i);
            }
        }
    } else {
        DebugStop();
    }

    for (int j = 0; j < nshape; j++){
        data.fVecShapeIndex[j].first = j;
        data.fVecShapeIndex[j].second = j;
    }

}//void


template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv3D<TSHAPE>::ComputeSolutionKernelHdivT(TPZMaterialDataT<TVar> &data)
{
    this->ComputeSolutionHCurlT(data.phi, data.curlphi,
                                data.sol, data.curlsol);

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
int TPZCompElKernelHDiv3D<TSHAPE>::NConnectShapeF(int icon, int order) const
{
    if (TSHAPE::Type() == ETetraedro || TSHAPE::Type() == ETriangle){
        order++;
        if (icon == 10){//the internal connect
            order++;
        }
    }
    return TPZShapeHCurlNoGrads<TSHAPE>::ComputeNConnectShapeF(icon,order);
}

// The MaxOrder of the boundary elements is increased by one to be compatible with the approximation space;
// For Tetrahedron, it needs to be increased by two as the internal connect has higher polynomial order
template<class TSHAPE>
int TPZCompElKernelHDiv3D<TSHAPE>::MaxOrder(){
    int maxorder = TPZInterpolationSpace::MaxOrder();
    
    return maxorder +1;
}


#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"

using namespace pzgeom;
using namespace pzshape;

template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeQuad>>;
template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeTetra>>;
template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapePrism>>;

template class TPZCompElKernelHDiv3D<TPZShapeTriang>;
template class TPZCompElKernelHDiv3D<TPZShapeQuad>;
template class TPZCompElKernelHDiv3D<TPZShapeTetra>;
template class TPZCompElKernelHDiv3D<TPZShapeCube>;
template class TPZCompElKernelHDiv3D<TPZShapePrism>;

#include "TPZCompElKernelHDiv.h"

//BC
TPZCompEl * CreateHDivKernelBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElKernelHDiv3D< TPZShapeQuad>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivKernelBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElKernelHDiv3D< TPZShapeTriang >(mesh,gel,hdivfam);
}

//Domain
TPZCompEl * CreateHDivKernelCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv3D< TPZShapeCube >(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivKernelPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv3D< TPZShapePrism>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivKernelTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv3D< TPZShapeTetra >(mesh,gel,hdivfam);
}

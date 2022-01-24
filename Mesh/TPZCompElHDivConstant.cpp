/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElHDivConstant.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"
#include "TPZShapeHDivConstant.h"
#include "TPZShapeHDivKernel2D.h"
#include "pzcmesh.h"
#include "TPZCompElH1.h"
#include "TPZCompElKernelHDiv.h"
#include "TPZShapeHCurl.h"
#include "TPZShapeHDivKernel.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

template<class TSHAPE>
TPZCompElHDivConstant<TSHAPE>::TPZCompElHDivConstant(TPZCompMesh &mesh, TPZGeoEl *gel) :
TPZRegisterClassId(&TPZCompElHDivConstant::ClassId), TPZCompElHDiv<TSHAPE>(mesh,gel) {
    this->AdjustConnects();
}

template<class TSHAPE>
void TPZCompElHDivConstant<TSHAPE>::AdjustConnects()
{
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    auto ncon = this->NConnects();
    if (TSHAPE::Dimension == 2) ncon =TSHAPE::NSides - nNodes;
    for(int icon = 0; icon < ncon; icon++){
        const int connect = icon;//this->MidSideConnectLocId(icon+1);
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
void TPZCompElHDivConstant<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZCompElHDiv<TSHAPE>::InitMaterialData(data);

    data.fH1ConnectOrders = data.fHDivConnectOrders;
    if (TSHAPE::Dimension == 2) TPZShapeH1<TSHAPE>::Initialize(data.fCornerNodeIds, data.fH1ConnectOrders, data);
    
    const int nSides = TSHAPE::NSides;
    const int nCorner = TSHAPE::NCornerNodes;
    TPZManVector<int64_t,nCorner> ids(nCorner,0);
    for(auto i=0; i<nCorner; i++) ids[i] = i;
    
    data.fSideTransformationId.Resize(nSides-nCorner, 0);
    for (int iside =nCorner; iside< nSides ; iside++) {
        int pos = iside - nCorner;
        int trans_id = TSHAPE::GetTransformId(iside, ids); // Foi criado
        data.fSideTransformationId[iside-nCorner] = trans_id;
    }

    int nshape = this->NShapeF();
    // int64_t numvec = TSHAPE::Dimension*TSHAPE::NSides;
    data.fMasterDirections.Resize(3, nshape);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < nshape; j++)
			data.fMasterDirections(i,j) = 1;

    data.divphi.Zero();
    
    data.fVecShapeIndex.Resize(nshape);
    for (int i=0; i<nshape; i++) {
		data.fVecShapeIndex[i] = std::make_pair(i,1);
    }
    data.fDeformedDirections.Resize(3,nshape);
}


template<class TSHAPE>
void TPZCompElHDivConstant<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    data.fNeedsSol = needsol;
    // TPZCompElHDiv<TSHAPE>::ComputeShape(qsi,data);

    TPZShapeData &shapedata = data;
    int nshape = this->NShapeF();
    TPZFMatrix<REAL> auxPhi(TSHAPE::Dimension,nshape);
    data.divphi.Resize(nshape,1);
    auxPhi.Zero();

    //Initialize data for HDiv Kernel
    auto &dataKernel = data;
    if(TSHAPE::Dimension == 3){
        TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes,0);
        TPZGeoEl *ref = this->Reference();
        for(auto i=0; i<TSHAPE::NCornerNodes; i++) {
            ids[i] = ref->NodePtr(i)->Id();
        }
        
        auto &conOrders = dataKernel.fHDivConnectOrders;
        constexpr auto nConnects = TSHAPE::NSides - TSHAPE::NCornerNodes;
        conOrders.Resize(nConnects,-1);
        for(auto i = 0; i < nConnects; i++){
            conOrders[i] = this->EffectiveSideOrder(i + TSHAPE::NCornerNodes);
        }
        for (int i = 0; i < TSHAPE::NumSides(1); i++)
        {
            conOrders[i] = conOrders[TSHAPE::NumSides(1)];
        }
        
        TPZShapeHCurl<TSHAPE>::Initialize(ids, conOrders, dataKernel);
        dataKernel.divphi.Resize(dataKernel.fVecShapeIndex.size(),1);
        TPZShapeHDivKernel<TSHAPE>::ComputeVecandShape(dataKernel);
    }

    // std::cout << "Trans1 = " << data.fTra

    if (TSHAPE::Dimension == 2) {
        TPZShapeHDivConstant<TSHAPE>::Shape(qsi, shapedata, auxPhi, data.divphi);
    } else if (TSHAPE::Dimension == 3){
        TPZShapeHDivConstant<TSHAPE>::Shape(qsi, dataKernel, auxPhi, data.divphi);
    } else {
        DebugStop();
    }
    // std::cout << "AuxPhi = " << auxPhi << std::endl;

    TPZFMatrix<REAL> gradx;
    this->Reference()->GradX(qsi,gradx);

    TPZFMatrix<REAL> phiSHdiv;
    gradx.Multiply(auxPhi,phiSHdiv);
    phiSHdiv *= 1./data.detjac;
    data.divphi *= 1./data.detjac;
    // std::cout << "gradx = " << gradx << std::endl;

    // std::cout << "div = " << data.divphi << std::endl;
    data.phi.Resize(data.fVecShapeIndex.size(),1);
    data.phi = 1.;
    data.fDeformedDirections = phiSHdiv;

}//void

template<class TSHAPE>
void TPZCompElHDivConstant<TSHAPE>::SetSideOrient(int orient){
    fSideOrient = orient;
}

template<class TSHAPE>
int TPZCompElHDivConstant<TSHAPE>::GetSideOrient(){
    return fSideOrient;
}

template<class TSHAPE>
void TPZCompElHDivConstant<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
{
    
    TPZMaterialDataT<STATE> data;
    constexpr bool hasPhi{false};
    this->ComputeSolution(qsi,data,hasPhi);

    sol.Resize(3);
    
    // REAL Sol = data.sol[0];
    // data.sol.resize(3);
    // data.sol[0] = Sol;
    sol[0] = data.sol[0][0];
    data.sol[0].Resize(3);
    // sol = std::move(data.sol[0]);
}


template<class TSHAPE>
int TPZCompElHDivConstant<TSHAPE>::NConnectShapeF(int icon, int order) const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif
    return TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(icon,order);
}


#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"

using namespace pztopology;

// #include "tpzpoint.h"
// #include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

using namespace pzgeom;
using namespace pzshape;

template class TPZCompElHDivConstant<TPZShapeTriang>;
template class TPZCompElHDivConstant<TPZShapeQuad>;
template class TPZCompElHDivConstant<TPZShapeTetra>;
template class TPZCompElHDivConstant<TPZShapeCube>;
template class TPZCompElHDivConstant<TPZShapePrism>;

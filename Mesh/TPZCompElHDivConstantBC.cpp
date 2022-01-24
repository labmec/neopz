/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElHDivConstantBC.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"
#include "TPZShapeHDivKernel2DBound.h"
#include "TPZShapeHDivConstantBound.h"
#include "TPZCompElHCurl.h"
#include "TPZShapeHDivKernel.h"

#include "pzcmesh.h"

template<class TSHAPE>
void TPZCompElHDivConstantBC<TSHAPE>::AdjustConnects()
{
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    auto ncon = 1;//this->NConnects();
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
TPZCompElHDivConstantBC<TSHAPE>::TPZCompElHDivConstantBC(TPZCompMesh &mesh, TPZGeoEl *gel, int shapetype) :
TPZRegisterClassId(&TPZCompElHDivConstantBC::ClassId), TPZCompElHDivBound2<TSHAPE>(mesh,gel), fShapeType(shapetype)  {
    AdjustConnects();
}

template<class TSHAPE>
TPZCompElHDivConstantBC<TSHAPE>::~TPZCompElHDivConstantBC(){
    this->~TPZCompElHDivBound2<TSHAPE>();
}

template<class TSHAPE>
void TPZCompElHDivConstantBC<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZCompElHDivBound2<TSHAPE>::InitMaterialData(data);
    // this->AdjustConnects();
    int nshape = this->NShapeF();
    data.fVecShapeIndex.Resize(nshape);
    for (int i=0; i<nshape; i++) {
		data.fVecShapeIndex[i] = std::make_pair(i,1);
    }

}

template<class TSHAPE>
void TPZCompElHDivConstantBC<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    data.fNeedsSol = needsol;

    TPZShapeData &shapedata = data;
    int nshape = this->NShapeF();
    data.phi.Resize(nshape,1);

    TPZFMatrix<REAL> auxPhi(nshape,TSHAPE::Dimension);
    auxPhi.Zero();

    auto &dataKernel = data;
    if (TSHAPE::Dimension == 2){
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

    if (TSHAPE::Dimension == 1){
        TPZShapeHDivConstantBound<TSHAPE>::Shape(qsi, shapedata, auxPhi);
    } else if (TSHAPE::Dimension == 2){
        TPZShapeHDivConstantBound<TSHAPE>::Shape(qsi, dataKernel, auxPhi);
    } else {
        DebugStop();
    }

    // if (this->GetSideOrient() != 1) DebugStop();
    if (data.fSideOrient[0] != 1) DebugStop();
    
    // TPZFMatrix<REAL> gradx;
    // this->Reference()->GradX(qsi,gradx);

    // if (TSHAPE::Dimension == 1) DebugStop();//Verificar se nesse caso também multiplica pelo gradiente.

    // std::cout << "gradx " << gradx << std::endl;
    // std::cout << "phi " << auxPhi << std::endl;
    // TPZFMatrix<REAL> phiSHdiv, phit;
    // auxPhi.Transpose(&phit);
    // gradx.Multiply(phit,phiSHdiv);
    // phiSHdiv *= 1./data.detjac;

    for (int i = 0; i < data.phi.Rows(); i++){
        data.phi(i,0) = auxPhi(i,0) / data.detjac;
    }
 
}//void

template<class TSHAPE>
int TPZCompElHDivConstantBC<TSHAPE>::NConnectShapeF(int connect, int connectorder) const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif
    return TPZShapeHDivConstantBound<TSHAPE>::ComputeNConnectShapeF(connect,connectorder);
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
template class TPZCompElHDivConstantBC<TPZShapeLinear>;
template class TPZCompElHDivConstantBC<TPZShapeTriang>;
template class TPZCompElHDivConstantBC<TPZShapeQuad>;

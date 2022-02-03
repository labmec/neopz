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
TPZCompElKernelHDiv3D<TSHAPE>::TPZCompElKernelHDiv3D(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElKernelHDiv3D::ClassId), TPZCompElHCurl<TSHAPE>(mesh,gel,HCurlFamily::EHCurlNoGrads),
                  fSideOrient(TSHAPE::NFacets,1), fhdivfam(hdivfam) {
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    for(int side = firstside ; side < TSHAPE::NSides-1; side++ )
    {
        fSideOrient[side-firstside] = this->Reference()->NormalOrientation(side);
    }
    
    if (fhdivfam != HDivFamily::EHDivKernel && fhdivfam != HDivFamily::EHCurlNoGrads){
        std::cout << "You need to chose EHivKernel approximation space to use TPZCompElKernelHDiv3D" << std::endl;
        DebugStop();
    }

}
 

template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data){
                                                
    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;

    int rows = data.phi.Rows();
    int cols = data.phi.Cols();
    TPZFMatrix<REAL> phiHCurl(cols,rows,0.);
    TPZFMatrix<REAL> curlHCurl(cols,rows,0.);
    TPZFMatrix<REAL> DivPhi(rows,1,0.);

    switch (fhdivfam)
    {
    case HDivFamily::EHCurlNoGrads:
        TPZShapeHCurlNoGrads<TSHAPE>::Shape(qsi,data,phiHCurl,curlHCurl);
        break;
    case HDivFamily::EHDivKernel:
        TPZShapeHDivKernel<TSHAPE>::Shape(qsi,data,curlHCurl,DivPhi);
        break;

    default:
        DebugStop();
        break;
    }
    
    data.divphi.Zero();
    TPZCompElHCurl<TSHAPE>::TransformCurl(curlHCurl, data.detjac, data.jacobian, data.curlphi);    
    
    data.fNeedsSol = needsol;
    if (TSHAPE::Dimension == 3){
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
        data.phi.Resize(nshape,3);
        data.phi = 1.;
        data.dphix = 1.;
        
    } else if (TSHAPE::Dimension == 2) {
        data.phi.Resize(data.curlphi.Cols(),3);
        data.phi.Zero();
        if (data.phi.Rows()>1){
            for (int i = 0; i < data.phi.Rows(); i++){
                data.phi(i,0) = -data.curlphi(0,i);
            }
        }
    } else {
        DebugStop();
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        //	this->Print(sout);
        // sout << "\nVecshape = " << data.fVecShapeIndex << std::endl;
        // sout << "fDeformedDirections = " << data.fDeformedDirections << std::endl;
        // sout << "phi = " << data.phi << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
        
    }
#endif

}//void

template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	data.fNeedsSol = true;
    
    // We create the derived element as a HCurlFamily::EHCurlNoGrads, so it will have the right connects.
    // However, at this point we need to initialize it as a HCurlFamily::EHCurlStandard in order to get
    // the hardcoded filtered functions a.k.a. TPZShapeHDivKernel.
    if (fhdivfam == HDivFamily::EHDivKernel) this->fhcurlfam = HCurlFamily::EHCurlStandard;

    //Init the material data of Hcurl
    TPZCompElHCurl<TSHAPE>::InitMaterialData(data);
	
    if (fhdivfam == HDivFamily::EHDivKernel) {   
        data.fVecShapeIndex = data.fSDVecShapeIndex;
        data.divphi.Resize(data.fVecShapeIndex.size(),1);
        TPZShapeHDivKernel<TSHAPE>::ComputeVecandShape(data);
            
        //setting the type of shape functions as vector shape functions
        data.fShapeType = TPZMaterialData::EVecShape;
    }
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
int TPZCompElKernelHDiv3D<TSHAPE>::NConnectShapeF(int icon, int order) const
{
    return TPZShapeHCurlNoGrads<TSHAPE>::ComputeNConnectShapeF(icon,order);
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

#include "TPZCompElKernelHDivBC.h"
#include "TPZCompElKernelHDiv.h"


//BC
TPZCompEl * CreateHDivKernelBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElH1<TPZShapePoint>(mesh,gel);
}
TPZCompEl * CreateHDivKernelBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDivBC< TPZShapeLinear>(mesh,gel);
}
TPZCompEl * CreateHDivKernelBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElKernelHDiv3D< TPZShapeQuad>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivKernelBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElKernelHDiv3D< TPZShapeTriang >(mesh,gel,hdivfam);
}

//Domain
TPZCompEl * CreateHDivKernelQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv< TPZShapeQuad>(mesh,gel);
}
TPZCompEl * CreateHDivKernelTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv< TPZShapeTriang >(mesh,gel);
}
TPZCompEl * CreateHDivKernelCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv3D< TPZShapeCube >(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivKernelPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv3D< TPZShapePrism>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivKernelTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElKernelHDiv3D< TPZShapeTetra >(mesh,gel,hdivfam);
}

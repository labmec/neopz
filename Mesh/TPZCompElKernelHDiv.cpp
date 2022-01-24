/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHDiv.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"
#include "TPZShapeHDivKernel2D.h"
#include <TPZCompElHCurl.h>
#include "pzcmesh.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::TPZCompElKernelHDiv(TPZCompMesh &mesh, TPZGeoEl *gel) :
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId), TPZCompElH1<TSHAPE>(mesh,gel) {

}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZCompElH1<TSHAPE>::InitMaterialData(data);

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
void TPZCompElKernelHDiv<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    data.fNeedsSol = needsol;
    TPZCompElH1<TSHAPE>::ComputeShape(qsi,data);

    TPZShapeData &shapedata = data;
    int nshape = this->NShapeF();
    TPZFMatrix<REAL> auxPhi(2,nshape);
    auxPhi.Zero();

    TPZShapeHDivKernel2D<TSHAPE>::Shape(qsi, shapedata, auxPhi, data.divphi);

    TPZFMatrix<REAL> gradx, aux2(3,auxPhi.Cols());
    this->Reference()->GradX(qsi,gradx);

    const auto nCorner = TSHAPE::NCornerNodes;

    aux2.Zero();
    gradx.MultAdd(auxPhi,data.phi,aux2,1/data.detjac,0);
    // std::cout << "AUX2 " << aux2 << std::endl;
    data.fDeformedDirections = aux2;

    data.phi.Resize(auxPhi.Cols(),3);
    data.divphi.Resize(auxPhi.Cols(),1);
    data.divphi.Zero();
    data.phi = 1.;
    // data.phi.Transpose(&aux2);
    // std::cout << "PHI = " << data.phi << std::endl;
// #ifdef PZ_LOG
//     if (logger.isDebugEnabled())
//     {
//         std::stringstream sout;
//         //	this->Print(sout);
//         // sout << "\nVecshape = " << data.fVecShapeIndex << std::endl;
//         // sout << "MASTER = " << data.fMasterDirections << std::endl;
//         sout << "\nPhi = " << data.phi << std::endl;
//         LOGPZ_DEBUG(logger,sout.str())
        
//     }
// #endif

}//void

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetSideOrient(int orient){
    fSideOrient = orient;
}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::GetSideOrient(){
    return fSideOrient;
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
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
template<class TVar>
void TPZCompElKernelHDiv<TSHAPE>::ComputeSolutionKernelHdivT(TPZMaterialDataT<TVar> &data)
{
    
    const int dim = 3; 
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
    TPZFNMatrix<220,REAL> dphix(3,data.dphix.Cols());
    TPZFMatrix<REAL> &dphi = data.dphix;;

    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, data.axes);

    TPZBlock &block =this->Mesh()->Block();
    int ishape=0,ivec=0,counter=0;



    // TPZFMatrix<TVar> dsolX(3,1);
    // // const auto &sol = data.sol[0];
    // const auto &dsol = data.dsol[0];
    // TPZAxesTools<TVar>::Axes2XYZ(dsol,dsolX,data.axes);

//     int nshapeV = data.fVecShapeIndex.NElements();

    for(int in=0; in<ncon; in++)
    {
        TPZConnect *df = &this->Connect(in);
        int64_t dfseq = df->SequenceNumber();
        int dfvar = block.Size(dfseq);
//         // pos : position of the block in the solution matrix
        int64_t pos = block.Position(dfseq);

//         /// ish loops of the number of shape functions associated with the block
        for(int ish=0; ish<dfvar/nstate; ish++)
        {
            ishape  = data.fVecShapeIndex[counter].first;
            for(int idf=0; idf<nstate; idf++)
            {
                TVar meshsol = MeshSol(pos+ish*nstate+idf,0);
                REAL phival = data.phi(ishape,0);
                //Computes sol and dsol
                // data.sol[0][dim*idf] += phival*meshsol;
                // data.dsol[0](dim*idf,0)+= meshsol * dphix(0,ishape);
                // data.dsol[0](dim*idf,1)+= meshsol * dphix(1,ishape);

                //Compute rotated flux
                data.sol[0][dim*idf+0] -= meshsol * dphix(1,ishape);
                data.sol[0][dim*idf+1] += meshsol * dphix(0,ishape);
            }
            counter++;
        }
    }
    // data.sol[1][0] = 0.;
    // data.sol[0][0] = -data.dsol[0](0,1);
    // data.sol[0][1] = data.dsol[0](0,0);
    // data.sol[0][2] = 0.;
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

template class TPZCompElKernelHDiv<TPZShapeTriang>;
template class TPZCompElKernelHDiv<TPZShapeQuad>;

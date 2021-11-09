#include "TPZShapeHDivConstant.h"

#include "TPZShapeH1.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "TPZShapeData.h"


template<class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NHDivShapeF(TPZShapeData &data)
{
    int nshape = 0;
    int nc = data.fHDivNumConnectShape.size();
    for(int ic = 0; ic<nc; ic++) nshape += data.fHDivNumConnectShape[ic];
    return nshape;
}

    

template<class TSHAPE>
void TPZShapeHDivConstant<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{
    
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    const int nfacets = TSHAPE::NFacets;
    TPZShapeH1<TSHAPE>::Shape(pt,data);
    divphi.Zero();
    phi.Zero();
    const auto nEdges = TSHAPE::NumSides(1);

    // Compute constant Hdiv functions
    TPZFNMatrix<12,REAL> vecDiv(dim,nfacets);
    TPZVec<REAL> div(nfacets);
    vecDiv.Zero();
    div.Fill(0.);
    TSHAPE::ComputeConstantHDiv(pt, vecDiv, div, data.fSideTransformationId);

    int countHDiv = 0;
    int countHCurl = 0;
    //Percorre as faces, atribui as funções Hdiv e depois Hcurl. Não tem função de aresta.
    //Atribui função HDiv
    //Atribui função HCurl
    int nEdgeF = 0;
    for (int i = 0; i < nEdges; i++)
    {
        nEdgeF += data.fHDivNumConnectShape[i];
    }

    int size = data.fSDVecShapeIndex.size();
    for(int i = 0; i < nfacets; i++){

        if(dim == 3){
            //Function HDiv
            phi(0,i) = div[i];
        } else if (dim == 2){
            phi(0,i) = div[i];
        } else {
            DebugStop();
        }   
    }

    int nconnects = data.fHDivNumConnectShape.size();
    int nInternal = data.fHDivNumConnectShape[nconnects-1];
    //Internal Functions    
    for(int i = size-nInternal; i<size; i++)
    {
        
        auto it = data.fSDVecShapeIndex[i];
        int vecindex = it.first;
        int scalindex = it.second;
        
        if(dim == 3)
        {
            for(auto d = 0; d < dim; d++) {
                const auto di = (d+1)%dim;
                const auto dj = (d+2)%dim;
                phi(d,i+nfacets) = data.fDPhi.GetVal(di,scalindex) * data.fMasterDirections.GetVal(dj,vecindex)-
                                   data.fDPhi.GetVal(dj,scalindex) * data.fMasterDirections.GetVal(di,vecindex);
            }
        }  else if (dim == 2)
        {
            phi(0,i+nfacets) = data.fDPhi.GetVal(0,scalindex) * data.fMasterDirections.GetVal(1,vecindex) -
                               data.fDPhi.GetVal(1,scalindex) * data.fMasterDirections.GetVal(0,vecindex);
        } else {
            DebugStop();
        }
    }

    //Atribui Função interna


    divphi.Zero();
}


template<class TSHAPE>
int TPZShapeHDivConstant<TSHAPE>::NConnectShapeF(int icon, TPZShapeData &data)
{
    int nshape = data.fHDivNumConnectShape[icon] + 1;
    return nshape;
}

template
struct TPZShapeHDivConstant<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapeQuad>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapeTetra>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapeCube>;

template
struct TPZShapeHDivConstant<pzshape::TPZShapePrism>;


#include "TPZShapeHDivKernel.h"

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
void TPZShapeHDivKernel<TSHAPE>::ComputeVecandShape(TPZShapeData &data) {
    
    const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto nConnects = TSHAPE::NSides - nNodes;

    for (int iedge = 0; iedge < nEdges; iedge++) {
        data.fHDivNumConnectShape[iedge] = 2;
        data.fHDivConnectOrders[iedge] = 1;
    }
    TPZShapeHCurl<TSHAPE>::ComputeVecandShape(data);
    
    typedef std::pair<MElementType,int> orderpair;
    std::map<orderpair ,std::set<int>> ShapeRemove;
    ShapeRemove[orderpair(ETriangle,2)] = {0};
    ShapeRemove[orderpair(ETriangle,3)] = {0,1,7};
    ShapeRemove[orderpair(ETetraedro,3)] ={0};

    for (int iedge = 0; iedge < nEdges; iedge++) {
        data.fHDivNumConnectShape[iedge] = 1;
        data.fHDivConnectOrders[iedge] = 0;
    }
    for (int connect = nEdges; connect < nConnects; connect++) {
        const int side = connect+nNodes;
        MElementType sidetype = TSHAPE::Type(side);
        int connectorder = data.fHDivConnectOrders[connect];
        if(ShapeRemove.find(orderpair(sidetype,connectorder)) == ShapeRemove.end()) DebugStop();
        auto remove = ShapeRemove[orderpair(sidetype,connectorder)];
        data.fHDivNumConnectShape[connect] = TPZShapeHCurl<TSHAPE>::ComputeNConnectShapeF(connect, connectorder)-remove.size();
        
        
    }
    TPZManVector<std::pair<int,int64_t>, 100> VecShapeIndex(NHDivShapeF(data)-nEdges);
    int countHDiv = 0;
    int countHCurl = 2 * nEdges;
    for (int connect = nEdges; connect < nConnects; connect++) {
        const int side = connect+nNodes;
        MElementType sidetype = TSHAPE::Type(side);
        int connectorder = data.fHDivConnectOrders[connect];
        auto remove = ShapeRemove[orderpair(sidetype,connectorder)];
        int nHCurlShape = TPZShapeHCurl<TSHAPE>::ComputeNConnectShapeF(connect, connectorder);
        for(int ishape = 0; ishape < nHCurlShape; ishape++)
        {
            // skip the removed shape functions
            if(remove.find(ishape) != remove.end())
            {
                countHCurl++;
                continue;
            }
            VecShapeIndex[countHDiv] = data.fVecShapeIndex[countHCurl];
            countHDiv++;
            countHCurl++;
        }
    }
    data.fVecShapeIndex = VecShapeIndex;
}

template<class TSHAPE>
int TPZShapeHDivKernel<TSHAPE>::NHDivShapeF(TPZShapeData &data)
{
    int nshape = 0;
    int nc = data.fHDivNumConnectShape.size();
    for(int ic = 0; ic<nc; ic++) nshape += data.fHDivNumConnectShape[ic];
    return nshape;
}

    

template<class TSHAPE>
void TPZShapeHDivKernel<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{
    


    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    TPZShapeH1<TSHAPE>::Shape(pt,data);
    divphi.Zero();
    const auto nEdges = TSHAPE::NumSides(1);

    // AQUI TEM QUE CALCULAR AS FUNCOES HCURL CONSTANTES
    // Jeferson deve ter colocado em outro branch
    DebugStop();
    
    for(int i = 0; i< data.fVecShapeIndex.size(); i++)
    {
        auto it = data.fVecShapeIndex[i];
        int vecindex = it.first;
        int scalindex = it.second;
        
        if(dim == 3)
        {
            for(int a=0; a<3; a++)
            {
                phi(a,i+nEdges) = data.fDPhi((a+2)%3,scalindex)*data.fMasterDirections((a+1)%3,vecindex) -
                            data.fDPhi((a+1)%3,scalindex)*data.fMasterDirections((a+2)%3,vecindex);
            }
        } else
        {
            DebugStop();
        }
    }
}


template<class TSHAPE>
int TPZShapeHDivKernel<TSHAPE>::NConnectShapeF(int icon, TPZShapeData &data)
{
    return data.fHDivNumConnectShape[icon];
}


template
struct TPZShapeHDivKernel<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDivKernel<pzshape::TPZShapeQuad>;

template
struct TPZShapeHDivKernel<pzshape::TPZShapeTetra>;

template
struct TPZShapeHDivKernel<pzshape::TPZShapeCube>;

template
struct TPZShapeHDivKernel<pzshape::TPZShapePrism>;


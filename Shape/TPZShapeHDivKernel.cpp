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
    ShapeRemove[orderpair(ETriangle,4)] = {6,9,10,11,13,14};
    ShapeRemove[orderpair(ETriangle,5)] = {8,12,13,14,16,17,20,21,22,23};
    ShapeRemove[orderpair(ETetraedro,3)] = {0};
    ShapeRemove[orderpair(ETetraedro,4)] = {9,12,13,14};
    ShapeRemove[orderpair(ETetraedro,5)] = {18,24,25,26,29,31,32,33,34,35};
    ShapeRemove[orderpair(EQuadrilateral,1)] = {3};
    ShapeRemove[orderpair(EQuadrilateral,2)] = {6,8,9,11};
    ShapeRemove[orderpair(EQuadrilateral,3)] = {9,12,13,15,17,18,20,22,23};
    ShapeRemove[orderpair(EQuadrilateral,4)] = {12,16,17,19,21,23,24,26,28,30,31,33,35,37,38,39};
    ShapeRemove[orderpair(EQuadrilateral,5)] = {15,20,21,23,25,27,29,30,32,34,36,38,39,41,43,45,47,48,50,52,54,56,57,58,59};
    ShapeRemove[orderpair(ECube,1)] = {5};
    ShapeRemove[orderpair(ECube,2)] = {11,17,19,20,21,24,25,26,27,28,29,30,31,32,33,34,35};
    ShapeRemove[orderpair(ECube,3)] = {22,23,24,25,26,37,38,40,41,42,43,44,45,46,47,50,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,
                                       71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107};
    ShapeRemove[orderpair(ECube,4)] = {19,23,27,31,37,38,40,41,42,43,46,47,51,55,59,63,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,
                                       86,87,90,91,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,
                                       120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,
                                       148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
                                       176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,
                                       204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239};
    ShapeRemove[orderpair(EPrisma,1)] = {};
    ShapeRemove[orderpair(EPrisma,2)] = {7,8};
    ShapeRemove[orderpair(EPrisma,3)] = {21,24,25,26,29,32,33,34,35};
    ShapeRemove[orderpair(EPrisma,4)] = {5,38,39,40,41,42,44,45,46,47,48,49,50,52,53,55,56,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,
                                         75,76,77,78,79,80,81,82,83,84,85,86,87,88,89};
    ShapeRemove[orderpair(EPrisma,5)] = {5,6,8,9,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,
                                         86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,
                                         116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
                                         144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179};
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
            VecShapeIndex[countHDiv] = data.fSDVecShapeIndex[countHCurl];
            countHDiv++;
            countHCurl++;
        }
    }
    data.fSDVecShapeIndex = VecShapeIndex;
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

    int curldim = dim < 3 ? 1 : dim;
    // AQUI TEM QUE CALCULAR AS FUNCOES HCURL CONSTANTES
    TPZFNMatrix<12,REAL> vecDiv(dim,nEdges), curl(curldim,nEdges);
    TSHAPE::ComputeConstantHCurl(pt, vecDiv, curl);
    phi.Zero();
    for (int connect = 0; connect < nEdges; connect++) {
        if(dim == 3)
        {
            for(int a=0; a<3; a++)
            {
                phi(a,connect) += curl(a,connect);
            }
        } else
        {
            DebugStop();
        }

    }

    int size = data.fSDVecShapeIndex.size();
    for(int i = nEdges; i< size; i++)
    {
        auto it = data.fSDVecShapeIndex[i];
        int vecindex = it.first;
        int scalindex = it.second;
        
        if(dim == 3)
        {
            for(int a=0; a<3; a++)
            {
                phi(a,i) += data.fDPhi((a+2)%3,scalindex)*data.fMasterDirections((a+1)%3,vecindex) -
                            data.fDPhi((a+1)%3,scalindex)*data.fMasterDirections((a+2)%3,vecindex);
            }
        } else
        {
            DebugStop();
        }
    }
    divphi.Zero();
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


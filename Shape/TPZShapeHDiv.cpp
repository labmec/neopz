#include "TPZShapeHDiv.h"

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
TPZShapeHDiv<TSHAPE>::TPZShapeHDiv()
{
    
}

template<class TSHAPE>
void TPZShapeHDiv<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
                                      const TPZVec<int> &connectorders,
                                      const TPZVec<int>& sideorient,
                                      TPZShapeData &data)
{
    
        
    data.fCornerNodeIds = ids;
    data.fSideOrient = sideorient;
    const int ncon = TSHAPE::NFacets+1;
    //data.fHDivConnectOrders = connectorders;TODOPHIL
    int scalarorder = connectorders[TSHAPE::NFacets]+1;
    TPZManVector<int,27> scalarOrders(TSHAPE::NSides-TSHAPE::NCornerNodes,scalarorder);
    

    TPZShapeH1<TSHAPE>::Initialize(data.fCornerNodeIds, scalarOrders, data);
    if(connectorders.size() != TSHAPE::NFacets+1) DebugStop();
    
    data.fHDivConnectOrders = connectorders;

    data.fHDivNumConnectShape.Resize(TSHAPE::NFacets+1);
    int nShape = 0;
    for (int i = 0; i < TSHAPE::NFacets+1; i++)
    {
        int order = data.fHDivConnectOrders[i];
        data.fHDivNumConnectShape[i] = ComputeNConnectShapeF(i,order);
        nShape += data.fHDivNumConnectShape[i];
    }
    
    data.fSDVecShapeIndex.Resize(nShape);

    ComputeMasterDirections(data);
    ComputeVecandShape(data);
    
    //Checks if the last connect order is >= then the other connects
    int size = data.fHDivConnectOrders.size();
    int maxOrder = data.fHDivConnectOrders[size-1];
    for (int i = 0; i < size-1; i++)
    {
        if (data.fHDivConnectOrders[i] > maxOrder){
            DebugStop();
        }
    }    
}

template<class TSHAPE>
void TPZShapeHDiv<TSHAPE>::ComputeMasterDirections(TPZShapeData &data)
{
    int64_t numvec = TSHAPE::Dimension*TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    data.fMasterDirections.Resize(3,numvec);
    TPZFNMatrix<9,REAL> gradx(3,TSHAPE::Dimension,0.);
    for (int i = 0; i < TSHAPE::Dimension; i++) {
        gradx(i,i) = 1.;
    }
    TSHAPE::ComputeHDivDirections(gradx, data.fMasterDirections);
    
    int firstface = TSHAPE::NSides - TSHAPE::NFacets - 1;
    int lastface = TSHAPE::NSides - 1;
    int cont = 0;
    for(int side = firstface; side < lastface; side++)
    {
        int nvec = TSHAPE::NContainedSides(side);
        for (int ivet = 0; ivet<nvec; ivet++)
        {
            for (int il = 0; il<dim; il++)
            {
              data.fMasterDirections(il,ivet+cont) *= data.fSideOrient[side-firstface];
            }
        }
        cont += nvec;
    }
    
}



template<class TSHAPE>
void TPZShapeHDiv<TSHAPE>::ComputeVecandShape(TPZShapeData &data) {
    const int64_t numvec = TSHAPE::Dimension*TSHAPE::NSides;
    TPZManVector<int,numvec> VectorSides(numvec),bilinear(numvec),directions(numvec);
    TPZManVector<int,numvec> normalsides(numvec);

    TSHAPE::GetSideHDivDirections(VectorSides,directions,bilinear,normalsides);

    if (data.fSDVecShapeIndex.size() == 0) {
        DebugStop();
    }
    
    // const int pressureorder = data.fHDivConnectOrders[TSHAPE::NFacets];//TODOPHIL
    const int pressureorder = data.fHDivConnectOrders[TSHAPE::NFacets];

    int nshape = TSHAPE::NShapeF(data.fH1ConnectOrders);

    int nexternalvectors = 0;
    TPZManVector<int> facevector(VectorSides.size(),TSHAPE::NSides-1);
    // compute the permutation which needs to be applied to the vectors to enforce compatibility between neighbours
    TPZManVector<int,81> vecpermute(TSHAPE::NSides*TSHAPE::Dimension);
    if (TSHAPE::Type() == EPiramide) {
        DebugStop();
    }
    int count = 0;
    for (int side = 0; side < TSHAPE::NSides; side++) {
        if (TSHAPE::SideDimension(side) != TSHAPE::Dimension -1) {
            continue;
        }
        TPZStack<int> smallsides;
        TSHAPE::LowerDimensionSides(side,smallsides);
//        TPZGeoElSide gelside(gel,side);
//        int nlowdim = gelside.NSides();
        const int nlowdim = smallsides.NElements()+1;
        TPZManVector<int,27> permgather(nlowdim);
        HDivPermutation(side,data, permgather);
        int counthold = count;
        for (int i=0; i<nlowdim; i++) {
            vecpermute[counthold+i] = permgather[i]+counthold;
            facevector[count] = side;
            count++;
        }
    }

    nexternalvectors = count;
    for (; count < vecpermute.size(); count++) {
        vecpermute[count] = count;
    }
    TPZGenMatrix<int> shapeorders(nshape,3);
    TSHAPE::ShapeOrder(data.fCornerNodeIds, data.fH1ConnectOrders, shapeorders);

//    int nshapeflux = NFluxShapeF();
//    IndexVecShape.Resize(nshapeflux);

    int scalarorder = data.fHDivConnectOrders[TSHAPE::NFacets]+1;
    // VectorSide indicates the side associated with each vector entry
    TPZManVector<int64_t,27> FirstIndex(TSHAPE::NSides+1);
    // the first index of the shape functions
    FirstShapeIndex(FirstIndex,scalarorder);

    int64_t nvec = VectorSides.NElements();
    count = 0;

    const int firstface = TSHAPE::NSides-TSHAPE::NFacets-1;
    for (int locvec = 0; locvec<nexternalvectors; locvec++) {
        int ivec = vecpermute[locvec];
        int side = VectorSides[ivec];
        int face = facevector[locvec];
//        int connectindex = SideConnectLocId(0, face);
        int connectindex = face-firstface;
//        int order = this->Connect(connectindex).Order();
        //int order = data.fHDivConnectOrders[connectindex];
        int order = data.fHDivConnectOrders[connectindex];//TODOPHIL

        int firstshape = FirstIndex[side];
        int lastshape = FirstIndex[side+1];

        for (int ish = firstshape; ish<lastshape; ish++) {
            int sidedimension = TSHAPE::SideDimension(side);
            int include=true;
            for (int d=0; d<sidedimension; d++)
            {
                if (shapeorders(ish,d) > order) {
                    include = false;
                }
            }
            if (include)
            {
                data.fSDVecShapeIndex[count] = std::make_pair(ivec, ish);
                count++;
            }
        }
    }


    for (int locvec = nexternalvectors; locvec<nvec; locvec++) {
        int ivec = vecpermute[locvec];
        int side = VectorSides[ivec];
        int bil = bilinear[ivec];
        int dir = directions[ivec];
        MElementType tipo = TSHAPE::Type(side);

        int firstshape = FirstIndex[side];
        int lastshape = FirstIndex[side+1];

        for (int ish = firstshape; ish<lastshape; ish++) {
            int sidedimension = TSHAPE::SideDimension(side);
            int maxorder[3] = {pressureorder,pressureorder,pressureorder};
            if (bil) {
                maxorder[dir]++;
            }
            int shord[3] = {0};
            int include=true;
            for (int d=0; d<sidedimension; d++)
            {
                if (tipo==ETriangle||tipo==ETetraedro)//
                {
                    shord[d] = shapeorders(ish,d);
                    int maxd = maxorder[d]+1;
                    if (shord[d] > maxd) {
                        include = false;
                    }
                }
                else if(tipo==EQuadrilateral)
                {
                    shord[d] = shapeorders(ish,d);
                    int maxd = maxorder[d];
                    if (shord[d] > maxd) {
                        include = false;
                    }
                }
                else if(tipo==ECube)
                {
                    shord[d] = shapeorders(ish,d);
                    if (shord[d] > maxorder[d]) {
                        include = false;
                    }
                }
                else if (tipo==EPrisma)
                {
                    //DebugStop();
                    if (shapeorders(ish,d) > maxorder[d]) {
                        include = false;
                    }
                }
                else if (tipo == EOned)
                {
                    shord[0] = shapeorders(ish,0);
                    if (shord[0] > maxorder[d]) {
                        include = false;
                    }
                }
                else if (tipo == EPiramide)
                {
                    shord[d] = shapeorders(ish,d);
                    if (shord[d] > maxorder[d]) {
                        include = false;
                    }
                }
                else
                {
                    DebugStop();
                }
            }
            if (include)
            {
                data.fSDVecShapeIndex[count] = std::make_pair(ivec, ish);
                count++;
            }
        }
    }

    int ivs =  data.fSDVecShapeIndex.size();
    if (count != ivs) {
        std::cout<<"count "<<count
                 <<"\nivs "<<ivs<<std::endl;
        DebugStop();
    }
    
}

template<class TSHAPE>
void TPZShapeHDiv<TSHAPE>::FillOrderScalarShapeFunctions(const TPZVec<int> &connectorders, TPZVec<int> &scalarorder)
{
    scalarorder.resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
    const int ncon = TSHAPE::NFacets;
    int internalorder = connectorders[ncon];
    scalarorder.Fill(internalorder+1);
#ifdef PZDEBUG
    for(int ic=0; ic<ncon-1; ic++)
    {
        if(connectorders[ic] > internalorder) DebugStop();
    }
#endif
    return;
}
    

template<class TSHAPE>
void TPZShapeHDiv<TSHAPE>::Shape(const TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{
    
    divphi.Resize(data.fSDVecShapeIndex.size(),1);
    divphi.Zero();
    phi.Resize(TSHAPE::Dimension,data.fSDVecShapeIndex.size());
    phi.Zero();

    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    TPZShapeH1<TSHAPE>::Shape(pt,data);
    for(int i = 0; i< data.fSDVecShapeIndex.size(); i++)
    {
        auto it = data.fSDVecShapeIndex[i];
        int vecindex = it.first;
        int scalindex = it.second;
        divphi(i,0) = 0.;
        for(int d = 0; d<TSHAPE::Dimension; d++)
        {
            phi(d,i) = data.fPhi(scalindex,0)*data.fMasterDirections(d,vecindex);
            divphi(i,0) += data.fDPhi(d,scalindex)*data.fMasterDirections(d,vecindex);
        }
    }
}

template<class TSHAPE>
void TPZShapeHDiv<TSHAPE>::FirstShapeIndex(TPZVec<int64_t> &Index, int &scalarorders) {
    Index[0] = 0;

    for(int iside=0;iside<TSHAPE::NSides;iside++)
    {
        int sideorder = 1;
        if (iside >= TSHAPE::NCornerNodes) {
            sideorder = scalarorders;
        }
        int temp = Index[iside] + TSHAPE::NConnectShapeF(iside,sideorder);
        Index[iside+1] = temp;
    }
}


template<class TSHAPE>
void TPZShapeHDiv<TSHAPE>::HDivPermutation(int side, TPZShapeData &data, TPZVec<int> &permutegather) {
    int dimension = TSHAPE::Dimension;
    int sidedimension = TSHAPE::SideDimension(side);

    if(dimension != sidedimension+1)
    {
        std::cout << "HDivPermutation called with wrong side parameter " << side << std::endl;
        DebugStop();
    }
    

    const int nsidenodes = TSHAPE::NSideNodes(side);
    TPZManVector<int64_t,4> id(nsidenodes);
    for(int inode=0; inode<nsidenodes; inode++)
    {

        int64_t nodeindex = TSHAPE::SideNodeLocId(side, inode);
//        id[inode] = NodePtr(nodeindex)->Id();
        id[inode] = data.fCornerNodeIds[nodeindex];
    }
    
    MElementType sidetype = TSHAPE::Type(side);
    HDivPermutation(sidetype, id, permutegather);

}

template<class TSHAPE>
void TPZShapeHDiv<TSHAPE>::HDivPermutation(MElementType eltype, const TPZVec<int64_t> &ids, TPZVec<int> &permutegather)
{
    int transformid;
    switch (eltype) {
        case EOned:
            transformid = pztopology::TPZLine::GetTransformId(ids);
            pztopology::TPZLine::GetSideHDivPermutation(transformid, permutegather);
            break;
        case EQuadrilateral:
            transformid = pztopology::TPZQuadrilateral::GetTransformId(ids);
            pztopology::TPZQuadrilateral::GetSideHDivPermutation(transformid, permutegather);
            break;
        case ETriangle:
            transformid = pztopology::TPZTriangle::GetTransformId(ids);
            pztopology::TPZTriangle::GetSideHDivPermutation(transformid, permutegather);
            break;
        case EPoint:
            transformid = 0;
            permutegather[0] = 0;
            break;
        default:
            DebugStop();
            break;
    }

}

template<class TSHAPE>
int TPZShapeHDiv<TSHAPE>::NConnectShapeF(int connect,const TPZShapeData &shapedata)
{
    return shapedata.fHDivNumConnectShape[connect];
}

template<class TSHAPE>
int TPZShapeHDiv<TSHAPE>::NShapeF(const TPZShapeData &shapedata)
{
    int nconnect = shapedata.fHDivNumConnectShape.size();
    int nshape = 0;
    for(int ic = 0; ic<nconnect; ic++) nshape += shapedata.fHDivNumConnectShape[ic];
    return nshape;
}


template<class TSHAPE>
int TPZShapeHDiv<TSHAPE>::ComputeNConnectShapeF(int connect, int order)
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif
    // int order = data.fHDivConnectOrders[connect];
    MElementType thistype = TSHAPE::Type();

    if(thistype == EOned)
    {
        if(connect < 2) return 1;
        else return order;
    }
    else if(thistype == ETriangle)
    {
        if(connect < TSHAPE::NFacets) return (order+1);
        else return (order+1)*(order+1)-1;
    }
    else if(thistype == EQuadrilateral)
    {
        if(connect < TSHAPE::NFacets) return (order+1);
        else return 2*order*(order+1);
    }
    else if(thistype == ETetraedro)
    {
        if(connect < TSHAPE::NFacets) return (order+1)*(order+2)/2;
        else return order*(order+2)*(order+3)/2;
    }
    else if(thistype == EPrisma)
    {
        if(connect == 0 || connect == 4) return (order+1)*(order+2)/2;
        else if(connect < TSHAPE::NFacets) return (order+1)*(order+1);
        else return order*order*(3*order+5)/2+7*order-2;
    }
    else if(thistype == ECube)
    {
        if(connect < TSHAPE::NFacets) return (order+1)*(order+1);
        else return 3*order*(order+1)*(order+1);
    }
    DebugStop();
    unreachable();
 }


template
struct TPZShapeHDiv<pzshape::TPZShapeLinear>;

template
struct TPZShapeHDiv<pzshape::TPZShapeTriang>;

template
struct TPZShapeHDiv<pzshape::TPZShapeQuad>;

template
struct TPZShapeHDiv<pzshape::TPZShapeTetra>;

template
struct TPZShapeHDiv<pzshape::TPZShapeCube>;

template
struct TPZShapeHDiv<pzshape::TPZShapePrism>;

template
struct TPZShapeHDiv<pzshape::TPZShapePiram>;

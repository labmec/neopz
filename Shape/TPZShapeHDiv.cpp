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
void TPZShapeHDiv<TSHAPE>::Initialize(TPZVec<int64_t> &ids,
                                             TPZVec<int> &connectorders,
                                             TPZVec<int>& sideorient, TPZShapeData &data)
{
    
        
    data.fCornerNodeIds = ids;
    data.fSideOrient = sideorient;
    const int ncon = TSHAPE::NFacets+1;
    data.fConnectOrders = connectorders;
    //data.fHDivConnectOrders = connectorders;TODOPHIL
    int scalarorder = connectorders[TSHAPE::NFacets];
    data.fConnectOrders.resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
    data.fConnectOrders.Fill(scalarorder);
    TPZShapeH1<TSHAPE>::Initialize(data.fCornerNodeIds, data.fConnectOrders, data.fSideOrient, data);

    ComputeMasterDirections(data);
    ComputeVecandShape(data);
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

    if (data.fVecShapeIndex.size() == 0) {
        DebugStop();
    }
    
    //const int pressureorder = data.fHDivConnectOrders[TSHAPE::NFacets];TODOPHIL
    const int pressureorder = data.fConnectOrders[TSHAPE::NFacets];

    TPZManVector<int,TSHAPE::NFacets+1> scalarorder(TSHAPE::NFacets+1,0); 
    FillOrderScalarShapeFunctions(data.fConnectOrders, scalarorder);
    int nshape = TSHAPE::NShapeF(scalarorder);

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
    TSHAPE::ShapeOrder(data.fCornerNodeIds, scalarorder, shapeorders);

//    int nshapeflux = NFluxShapeF();
//    IndexVecShape.Resize(nshapeflux);

    // VectorSide indicates the side associated with each vector entry
    TPZManVector<int64_t,27> FirstIndex(TSHAPE::NSides+1);
    // the first index of the shape functions
    FirstShapeIndex(FirstIndex,scalarorder);

    int64_t nvec = VectorSides.NElements();
    count = 0;

    for (int locvec = 0; locvec<nexternalvectors; locvec++) {
        int ivec = vecpermute[locvec];
        int side = VectorSides[ivec];
        int face = facevector[locvec];
//        int connectindex = SideConnectLocId(0, face);
        int connectindex = face-(TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1);
//        int order = this->Connect(connectindex).Order();
        //int order = data.fHDivConnectOrders[connectindex];
        int order = data.fConnectOrders[connectindex];//TODOPHIL

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
                data.fVecShapeIndex[count] = std::make_pair(ivec, ish);
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
                data.fVecShapeIndex[count] = std::make_pair(ivec, ish);
                count++;
            }
        }
    }

    int ivs =  data.fVecShapeIndex.size();
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
void TPZShapeHDiv<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi)
{
    
//    int64_t numvec = TSHAPE::Dimension*TSHAPE::NSides;
//    TPZManVector<int64_t,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
//    TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
//    int i;
//    TPZGeoEl *ref = this->Reference();
//    for(i=0; i<TSHAPE::NCornerNodes; i++) {
//        id[i] = ref->NodePtr(i)->Id();
//    }


//    FillOrder(ord);
//    int nshape= this->NShapeContinuous(ord);

    TPZShapeH1<TSHAPE>::Shape(pt,data);
    for(int i = 0; i< data.fVecShapeIndex.size(); i++)
    {
        auto it = data.fVecShapeIndex[i];
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
void TPZShapeHDiv<TSHAPE>::FirstShapeIndex(TPZVec<int64_t> &Index, TPZVec<int> &scalarorders) {
    Index[0] = 0;

    for(int iside=0;iside<TSHAPE::NSides;iside++)
    {
        int sideorder = 1;
        if (iside >= TSHAPE::NCornerNodes) {
            sideorder = scalarorders[iside-TSHAPE::NCornerNodes];
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
    int transformid;
    switch (sidetype) {
        case EOned:
            transformid = pztopology::TPZLine::GetTransformId(id);
            pztopology::TPZLine::GetSideHDivPermutation(transformid, permutegather);
            break;
        case EQuadrilateral:
            transformid = pztopology::TPZQuadrilateral::GetTransformId(id);
            pztopology::TPZQuadrilateral::GetSideHDivPermutation(transformid, permutegather);
            break;
        case ETriangle:
            transformid = pztopology::TPZTriangle::GetTransformId(id);
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

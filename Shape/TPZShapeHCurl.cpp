#include "TPZShapeHCurl.h"

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
TPZShapeHCurl<TSHAPE>::TPZShapeHCurl()
{
    
}

template<class TSHAPE>
void TPZShapeHCurl<TSHAPE>::Initialize(TPZVec<int64_t> &ids,
                                             TPZVec<int> &connectorders,
                                             TPZShapeData &data)
{
    
        
    constexpr int ncon = TSHAPE::NSides-TSHAPE::NCornerNodes;
    constexpr int NCorners = TSHAPE::NCornerNodes;
    constexpr int NSides = TSHAPE::NSides;
    data.fCornerNodeIds = ids;
    if(connectorders.size() != ncon) DebugStop();
    //data.fHDivConnectOrders = connectorders;TODOPHIL
    StaticCalcH1ShapeOrders(connectorders, data.fH1ConnectOrders);
    TPZShapeH1<TSHAPE>::Initialize(data.fCornerNodeIds, data.fH1ConnectOrders, data);
    data.fSideTransformationId.Resize(NSides-NCorners, 0);
    for (int iside = NCorners; iside< NSides ; iside++) {
        int pos = iside - NCorners;
        int trans_id = TSHAPE::GetTransformId(iside, ids); // Foi criado
        data.fSideTransformationId[iside-NCorners] = trans_id;
    }

    data.fHDivConnectOrders = connectorders;

    data.fHDivNumConnectShape.Resize(ncon);
    int nShape = 0;
    for (int i = 0; i < ncon; i++)
    {
        data.fHDivNumConnectShape[i] = NConnectShapeF(i,data);
        nShape += data.fHDivNumConnectShape[i];
    }
    
    data.fVecShapeIndex.Resize(nShape);
    TPZFNMatrix<9,REAL> gradX(TSHAPE::Dimension, TSHAPE::Dimension, 0);
    gradX.Identity();

    data.fMasterDirections.Redim(TSHAPE::Dimension, 3*NSides);
    TSHAPE::ComputeHCurlDirections(gradX,data.fMasterDirections,data.fSideTransformationId);

    ComputeVecandShape(data);
}




template<class TSHAPE>
void TPZShapeHCurl<TSHAPE>::ComputeVecandShape(TPZShapeData &data) {
    /******************************************************************************************************************
    * The directions are calculated based on the LOCAL side ids (and SideNodeLocId), such as the H1 shape functions.
    * For instance, for the triangle, the vectors are:
    * vea30 vea31 vea41 vea42 vea52 vea51 vet3 vet4 vet5 vfe63 vfe64 vfe65 vft1 vft2
    * and the shapes will be organized as follows:
    * phi0 phi1 phi2  phi31 phi32 ... phi3i   phi41 phi42 ... phi4j   phi51 phi52 ... phi5k   phi61 phi62 ... phi6n
    *^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^
    *  corner funcs       edge 3 funcs            edge 4 funcs            edge 5 funcs            internal funcs
    *
    * In order to ensure that the functions will coincide in two neighbouring elements, they will be then sorted by
    * their sides' GLOBAL ids instead of their LOCAL ids
    ******************************************************************************************************************/

#ifdef PZ_LOG2
    if (logger.isDebugEnabled()) {
        std::ostringstream sout;
        sout << std::endl;
        sout << "ELEMENT ID: " << this->Reference()->Id() << std::endl;
        sout << "ELEMENT TYPE: " << this->Reference()->TypeName() << std::endl;
        sout << "CONNECT ORDERS: " << std::endl;
        for (auto &iCon : connectOrder) {
            sout << "\t" << iCon << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto nConnects = TSHAPE::NSides - nNodes;
    TPZManVector<int64_t,nNodes> nodes(nNodes, 0);
    nodes = data.fCornerNodeIds;


#ifdef PZ_LOG2
    if (logger.isDebugEnabled()) {
        std::ostringstream sout;
        sout << "first H1 shape function:" << std::endl;
        for (auto &iShape : firstH1ShapeFunc) {
            sout << "\t" << iShape << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int nshape = NHCurlShapeF(data);
    data.fVecShapeIndex.Resize(nshape);

    TPZVec<unsigned int> shapeCountVec(TSHAPE::NSides - nNodes, 0);
    TPZVec<std::pair<int,int64_t>> & indexVecShape = data.fVecShapeIndex;
    TPZVec<int> &connOrder = data.fHDivConnectOrders;
    TPZManVector<int64_t, TSHAPE::NSides - nNodes> firstH1ShapeFunc(TSHAPE::NSides - nNodes,
                                                                                  0);
    firstH1ShapeFunc[0] = nNodes;
    for (int iSide = nNodes + 1; iSide < TSHAPE::NSides; iSide++) {
        const int iCon = iSide - nNodes;
        firstH1ShapeFunc[iCon] = firstH1ShapeFunc[iCon - 1] + data.fH1NumConnectShape[iCon-1];
    }
    TPZVec<int> &sidesH1Ord = data.fH1ConnectOrders;
    auto &nodeIds = data.fCornerNodeIds;
    StaticIndexShapeToVec(data.fVecShapeIndex, connOrder, firstH1ShapeFunc, sidesH1Ord, shapeCountVec, nodeIds);
}

template<class TSHAPE>
int TPZShapeHCurl<TSHAPE>::NHCurlShapeF(TPZShapeData &data)
{
    int nshape = 0;
    int nc = data.fHDivNumConnectShape.size();
    for(int ic = 0; ic<nc; ic++) nshape += data.fHDivNumConnectShape[ic];
    return nshape;
}

    

template<class TSHAPE>
void TPZShapeHCurl<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &curlphi)
{
    


    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    const int dim = TSHAPE::Dimension;
    TPZShapeH1<TSHAPE>::Shape(pt,data);
    for(int i = 0; i< data.fVecShapeIndex.size(); i++)
    {
        auto it = data.fVecShapeIndex[i];
        int vecindex = it.first;
        int scalindex = it.second;
        
        for(int d = 0; d<TSHAPE::Dimension; d++)
        {
            phi(d,i) = data.fPhi(scalindex,0)*data.fMasterDirections(d,vecindex);
        }
        if(dim == 3)
        {
            for(int a=0; a<3; a++)
            {
                curlphi(a,i) = data.fDPhi((a+2)%3,scalindex)*data.fMasterDirections((a+1)%3,vecindex) -
                            data.fDPhi((a+1)%3,scalindex)*data.fMasterDirections((a+2)%3,vecindex);
            }
        } else
        {
            int a = 2;
            curlphi(0,i) = data.fDPhi((a+2)%3,scalindex)*data.fMasterDirections((a+1)%3,vecindex) -
            data.fDPhi((a+1)%3,scalindex)*data.fMasterDirections((a+2)%3,vecindex);
        }
    }
}


template<class TSHAPE>
int TPZShapeHCurl<TSHAPE>::NConnectShapeF(int icon, TPZShapeData &data)
{
    const int order = data.fHDivConnectOrders[icon];
    const int side = icon + TSHAPE::NCornerNodes;
    #ifdef PZDEBUG
        if (side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) {
            DebugStop();
        }
    #endif
        //given the choice of implementation, there are no shape functions for k=0
        if(order == 0) return 0;
        const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
        const auto nEdges = TSHAPE::NumSides(1);
        const int nShapeF = [&](){
            if (side < TSHAPE::NCornerNodes + nEdges) {//edge connect
                return 1 + order;
            }
            else if(side < TSHAPE::NCornerNodes + nEdges + nFaces){//face connect
                switch(TSHAPE::Type(side)){
                case ETriangle://triangular face
                    return (order - 1) * (order+1);
                case EQuadrilateral://quadrilateral face
                    return 2 * order * (order+1);
                default:
                    PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
                    DebugStop();
                    return 0;
                }
            }
            else{//internal connect (3D element only)
                int count = 0;
                //first we count the face-based interior functions \phi^{K,F}
                for(int iFace = 0; iFace < nFaces; iFace++){
                    switch(TSHAPE::Type(TSHAPE::NCornerNodes+nEdges + iFace)){
                        case ETriangle://triangular face
                            count +=  (order - 1) * ( order - 2)/2;
                            break;
                        case EQuadrilateral://quadrilateral face
                            //we need the functions of k+1
                            count +=  order * order;
                            break;
                        default:
                            PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
                            DebugStop();
                    }
                }

                const int nVkf = count;
                //now we count the interior bubble functions
                
                //number of H1 funcs
                const auto orderH1 = TSHAPE::Type() == ECube ?
                    order +1 : order;
                
                const auto nFuncsH1 = TSHAPE::NConnectShapeF(side, orderH1);
                if constexpr (TSHAPE::Type() == ECube) {
                  //shapeorders[k] contains the polynomial orders
                  //in xi, eta and zeta for each of the internal funcs
                  TPZGenMatrix<int> shapeorders(nFuncsH1,3);
                  //not really used since we are interested in internal funcs
                  TPZManVector<int64_t,0> idvec(0);
                  TSHAPE::SideShapeOrder(side,idvec,orderH1,shapeorders);
                  /*for compatibility, we need the spaces
                    Q^{k,k+1,k+1} \times Q^{k+1,k,k+1} \times Q^{k+1,k+1,k}
                    I am supposing that #funcs in xi-dir is equal to
                    #funcs in eta-dir and in zeta-dir. So I will just
                    see how many we will have in xi-dir after skipping
                    internals that were included in the k-1 set of funcs*/
                  int countInt = 0;
                  for(int ifunc = 0; ifunc < nFuncsH1; ifunc++){
                    if(shapeorders(ifunc,0) <= order) {countInt++;}
                  }
                  count += 3* countInt;

                  const auto nVki = count - nVkf;
                  // std::stringstream sout;
                  // sout << __PRETTY_FUNCTION__<<'\n'
                  //      <<"\tside "<<side<<"\tcon "<<icon<<"\torder "<<order<<'\n'
                  //      <<"\tnfuncs "<<count<<"\tnvkf "<<nVkf<<"\tnvki "<<nVki<<'\n';
                  // std::cout<<sout.str()<<std::flush;
                }
                else{
                  count += 3 * nFuncsH1;
                }
                return count;
            }
        }();

    #ifdef PZ_LOG2
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__<<std::endl;
            sout<<"\tside "<<side<<"\tcon "<<icon<<std::endl;
            sout<<"\torder "<<order<<std::endl;
            sout<<"\tn shape funcs "<<nShapeF;
            LOGPZ_DEBUG(logger,sout.str())
        }
    #endif
        return nShapeF;
    
 }

template<class TSHAPE>
template<class TSIDESHAPE>
void TPZShapeHCurl<TSHAPE>::StaticCalcH1ShapeOrders(
    const TPZVec<int> &ordHCurl,
    TPZVec<int> &ordH1)
{
    constexpr auto nConnects = TSIDESHAPE::NSides-TSIDESHAPE::NCornerNodes;
    if(ordH1.size() != nConnects)  ordH1.Resize(nConnects,-1);
    for(auto iCon = 0; iCon < nConnects; iCon++){
        const auto iSide = iCon + TSIDESHAPE::NCornerNodes;
        const auto sideDim = TSIDESHAPE::SideDimension(iSide);
        const bool quadSide = TSIDESHAPE::Type(iSide) == EQuadrilateral ||
            TSIDESHAPE::Type(iSide) == ECube;
        /*some H1 functions associated with the side iSide of dimension dim
          might be needed for computing the shape functions of a side with
          dimension dim+1 that contains the side iSide.
          It is also worth noting that quadrilateral sides require functions
          of ordH1er k+1*/
        TPZStack<int> highDimSides;
        TSIDESHAPE::HigherDimensionSides(iSide, highDimSides);
        const auto sideOrder = ordHCurl[iCon];
        auto maxOrder = quadSide ? sideOrder + 1: sideOrder;
        for(auto &iHighSide : highDimSides){
            if(TSIDESHAPE::SideDimension(iHighSide) != sideDim+1) break;
            else {
                const auto hSideOrder = ordHCurl[iHighSide-TSIDESHAPE::NCornerNodes];
                const auto hQuadSide = TSIDESHAPE::Type(iHighSide) == EQuadrilateral ||
                    TSIDESHAPE::Type(iSide) == ECube;
                const auto hMaxOrder = hQuadSide ? hSideOrder + 1 : hSideOrder;
                maxOrder = std::max(maxOrder, hMaxOrder);
            }
        }
        ordH1[iCon] = maxOrder;
    }
}

template<class TSHAPE>
template<class TSIDESHAPE>
void TPZShapeHCurl<TSHAPE>::StaticIndexShapeToVec(TPZVec<std::pair<int,int64_t>> & indexVecShape,
                                                  const TPZVec<int>& connectOrder,
                                                  const TPZVec<int64_t>& firstH1ShapeFunc,
                                                  const TPZVec<int>& sidesH1Ord,
                                                  TPZVec<unsigned int>& shapeCountVec,
                                                  const TPZVec<int64_t>& nodeIds) {



//                                                  TPZVec<std::pair<int,int64_t>> & indexVecShape,
//                                                       const TPZVec<int>& connectOrder,
//                                                       const TPZVec<int64_t>& firstH1ShapeFunc,
//                                                       const TPZVec<int>& sidesH1Ord,
//                                                       TPZVec<unsigned int>& shapeCountVec,
//                                                       const TPZVec<int64_t>& nodeIds) {
    /******************************************************************************************************************
    * The directions are calculated based on the LOCAL side ids (and SideNodeLocId), such as the H1 shape functions.
    * For instance, for the triangle, the vectors are:
    * vea30 vea31 vea41 vea42 vea52 vea51 vet3 vet4 vet5 vfe63 vfe64 vfe65 vft1 vft2
    * and the shapes will be organized as follows:
    * phi0 phi1 phi2  phi31 phi32 ... phi3i   phi41 phi42 ... phi4j   phi51 phi52 ... phi5k   phi61 phi62 ... phi6n
    *^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^
    *  corner funcs       edge 3 funcs            edge 4 funcs            edge 5 funcs            internal funcs
    *
    * In order to ensure that the functions will coincide in two neighbouring elements, they will be then sorted by
    * their sides' GLOBAL ids instead of their LOCAL ids
    ******************************************************************************************************************/

    
    const auto nFaces = TSIDESHAPE::Dimension < 2 ? 0 : TSIDESHAPE::NumSides(2);
    const auto nEdges = TSIDESHAPE::NumSides(1);
    constexpr auto nNodes = TSIDESHAPE::NCornerNodes;

    const auto nConnects = connectOrder.size();
    TPZManVector<int, TSIDESHAPE::NSides - nNodes>
        transformationIds(TSHAPE::NSides - nNodes, -1);
    //computing transformation id for sides.
    for (auto iCon = 0; iCon < nConnects; iCon++) {
        transformationIds[iCon] = TSIDESHAPE::GetTransformId(nNodes + iCon, nodeIds);
    }
    
    unsigned int shapeCount = 0;


    //calculates edge functions
    for (auto iCon = 0; iCon < nEdges; iCon++) {
        const auto pOrder = connectOrder[iCon];
        //there will be 2 + pOrder - 1 = pOrder + 1 functions for each edge

        for (auto iNode = 0; iNode < 2; iNode++) {
            auto whichNode = transformationIds[iCon] == 0 ? iNode : (iNode + 1) % 2;
            const int vecIndex = iCon * 2 + whichNode;
            const int64_t shapeIndex = TSIDESHAPE::SideNodeLocId(iCon + nNodes, whichNode);
            indexVecShape[shapeCount] = std::make_pair(vecIndex, shapeIndex);
            shapeCount++;
            shapeCountVec[iCon]++;
        }//phi^ea funcs
        const int vecIndex = nEdges * 2 + iCon;
        for (int iEdgeInternal = 0; iEdgeInternal < pOrder - 1; iEdgeInternal++) {
            const int shapeIndex = firstH1ShapeFunc[iCon] + iEdgeInternal;
            indexVecShape[shapeCount] = std::make_pair(vecIndex, shapeIndex);
            shapeCount++;
            shapeCountVec[iCon]++;
        }//phi^et funcs
    }

    const int firstFaceShape = shapeCount;

#ifdef PZ_LOG2
    if (logger.isDebugEnabled()) {
        std::ostringstream sout;
        sout << __PRETTY_FUNCTION__ << '\n'
             << " n shape funcs (edge connects): "<<shapeCount<<'\n';
        //thats way too much info, uncomment if needed
        
        // sout << "vec shape index (edge connects):" << std::endl;
        // for (int iShape = 0; iShape < firstFaceShape; iShape++) {
        //     auto pair = indexVecShape[iShape];
        //     sout << "\tvec: " << pair.first << "\tshape: " << pair.second << std::endl;
        // }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    if(TSIDESHAPE::Dimension < 2) return;
    /**
     * In order to ease the calculation of the indexes, these structures will store, respectively:
     * a) firstVfeVec the first Vfe vector for each face
     * b) faceEdges the edges contained in a certain face
     * c) firstVftVec the first Vft vector. Then, it is firstVftVec + 2*iFace + iVec
     * d) firstVfOrthVec the first VfOrthVec. Then, it is firstVfOrthVec + iFace
     */
    TPZManVector<int> firstVfeVec(nFaces,-1);
    TPZManVector<TPZStack<int>> faceEdges(nFaces,TPZStack<int>(0,0));
    {
        //we skip v^{e,a} and v^{e,T} vectors
        const int nEdgeVectors = nEdges * 3;
        firstVfeVec[0] = nEdgeVectors;
        for(auto iFace = 0; iFace < nFaces; iFace++){
            TSIDESHAPE::LowerDimensionSides(iFace + nEdges + nNodes, faceEdges[iFace], 1);
            const int nFaceEdges = faceEdges[iFace].size();
            firstVfeVec[iFace] = iFace == 0 ?
                firstVfeVec[iFace] : firstVfeVec[iFace - 1] + nFaceEdges;
        }
    }
    const int firstVftVec = firstVfeVec[nFaces-1] + faceEdges[nFaces-1].size();
    const int firstVfOrthVec = firstVftVec + 2 * nFaces;

    for(auto iCon = nEdges; iCon < nEdges + nFaces; iCon++){
        const auto iSide = iCon + nNodes;
        const auto iFace = iCon - nEdges;
        TPZManVector<int,4> permutedSideSides = [&](){
            TPZVec<int> perm;
            switch(TSIDESHAPE::Type(iSide)){
                case ETriangle://triangular face
                    pztopology::GetPermutation<pztopology::TPZTriangle>(transformationIds[iCon],perm);
                    break;
                case EQuadrilateral://quadrilateral face
                    pztopology::GetPermutation<pztopology::TPZQuadrilateral>(transformationIds[iCon],perm);
                    break;
                default:
                    PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
                    DebugStop();
            }
            return perm;
        }();
#ifdef PZ_LOG2
        if (logger.isDebugEnabled()) {
            std::ostringstream sout;
            sout << "face :"<< iSide <<" permutation:"<< std::endl;
            for (auto i = 0; i < permutedSideSides.size(); i++) sout << permutedSideSides[i]<<"\t";
            sout<<std::endl;
            sout<<"transformation id:"<<transformationIds[iCon]<<std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif


        const int nFaceNodes = TSIDESHAPE::NSideNodes(iSide);
        //this is not a mistake, since for faces nEdges = nNodes
        const int &nFaceEdges = nFaceNodes;

        //first the phi Fe functions
        for(auto iEdge = 0; iEdge < nFaceEdges; iEdge++ ){
            const auto currentLocalEdge = permutedSideSides[iEdge+nFaceNodes];
            const auto currentEdge = TSIDESHAPE::ContainedSideLocId(iSide, currentLocalEdge);
            const auto vecIndex = firstVfeVec[iFace] + currentLocalEdge - nFaceNodes;
            const auto pOrder = sidesH1Ord[currentEdge-nNodes];
            for(auto iEdgeInternal = 0; iEdgeInternal < pOrder - 1; iEdgeInternal++){
                const int shapeIndex = firstH1ShapeFunc[currentEdge - nNodes] + iEdgeInternal;
                indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
                shapeCount++;
                shapeCountVec[iCon]++;
            }
        }

        const auto nVfeFuncs = shapeCountVec[iCon];
        
        /**now the phi Fi functions
         They are calculated differently for quad and triangular faces.
         For triangular faces, all the face functions of order k are taken,
         and each is multiplied by each of the tangent vectors associated with
         the face.
         For quadrilateral faces, the face functions of order k+1 are taken:
         since we are interested in the Q_{k,k+1}\times Q_{k+1,k} space,
         we check their orders.*/

        //max order of h1 functions associated with face
        const auto h1FaceOrder = sidesH1Ord[nEdges+iFace];
        //number of h1 face funcs
        const auto nH1FaceFuncs =
            TSIDESHAPE::NConnectShapeF(iSide,h1FaceOrder);
        
        const auto quadFace = TSIDESHAPE::Type(iSide) == EQuadrilateral;
        
        if(quadFace){
            TPZGenMatrix<int> shapeorders(nH1FaceFuncs,3);
            TSIDESHAPE::SideShapeOrder(iSide,nodeIds,h1FaceOrder,shapeorders);
            const auto hCurlFaceOrder = h1FaceOrder-1;
            /**now we assume that the first vft vec is in the x direction
               and that the next one is in the y direction.*/
            const auto xVecIndex = firstVftVec + 2*iFace;
            const auto yVecIndex = xVecIndex + 1;
            for(auto iFunc = 0; iFunc < nH1FaceFuncs; iFunc++ ){
                const auto shapeIndex = firstH1ShapeFunc[iCon] + iFunc;
                if(shapeorders(iFunc,0) <= hCurlFaceOrder){
                    indexVecShape[shapeCount] = std::make_pair(xVecIndex,shapeIndex);
                    shapeCount++;
                    shapeCountVec[iCon]++;
                }
                if(shapeorders(iFunc,1) <= hCurlFaceOrder){
                    indexVecShape[shapeCount] = std::make_pair(yVecIndex,shapeIndex);
                    shapeCount++;
                    shapeCountVec[iCon]++;
                }
            }
        }
        else{
            //ok that one is easy to guess
            const auto nFaceInternalFuncs =
                2 * nH1FaceFuncs;
            for(auto iFunc = 0; iFunc < nFaceInternalFuncs; iFunc++ ){
                //it should alternate between them
                const auto vecIndex = firstVftVec + 2*iFace + iFunc % 2;
                //they should repeat
                const auto shapeIndex = firstH1ShapeFunc[iCon] + iFunc / 2;
                indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
                shapeCount++;
                shapeCountVec[iCon]++;
            }
        }

        const auto nVfiFuncs = shapeCountVec[iCon] - nVfeFuncs;
#ifdef PZ_LOG2
    if (logger.isDebugEnabled()) {
        std::ostringstream sout;
        sout << "iFace: "<<iFace<<' '
             << "nVfeFuncs: "<< nVfeFuncs <<' '
             << "nVfiFuncs: "<< nVfiFuncs << '\n';
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    }
    const int firstInternalShape = shapeCount;

#ifdef PZ_LOG2
    if (logger.isDebugEnabled()) {
        std::ostringstream sout;
        sout << "n shape funcs (face connects): "
             << firstInternalShape - firstFaceShape << '\n';
        //way too much info, uncomment if needed
        // sout << "vec shape index (face connects):" << std::endl;
        // for (int iShape = firstFaceShape; iShape < firstInternalShape; iShape++) {
        //     auto pair = indexVecShape[iShape];
        //     sout << "\tvec: " << pair.first << "\tshape: " << pair.second << std::endl;
        // }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    if(TSIDESHAPE::Dimension < 3) return;
    else {
        const int iCon = nEdges + nFaces;
        //first, the phi KF functions
        for(int iFace = 0; iFace < nFaces; iFace++){
            const auto faceOrder = sidesH1Ord[nEdges+iFace];
            const auto nFaceInternalHCurlFuncs =
                TSIDESHAPE::NConnectShapeF(iFace + nEdges + nNodes,faceOrder);
            const auto vecIndex = firstVfOrthVec + iFace;
            for(auto iFunc = 0; iFunc < nFaceInternalHCurlFuncs; iFunc++ ){
                const auto shapeIndex = firstH1ShapeFunc[nEdges + iFace] + iFunc;
                indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
                shapeCount++;
                shapeCountVec[iCon]++;
            }
        }

        const auto nKfFuncs = shapeCount - firstInternalShape;
        //for the hexahedral element we need some internal functions of k+1
        const auto h1SideOrder = sidesH1Ord[nEdges+nFaces];
        //hcurl connect order
        const auto sideOrder = connectOrder[iCon];
        //now the phi Ki funcs
        const int firstInternalVec = firstVfOrthVec + nFaces;
        //ALL H1 internal functions
        const auto nH1Internal =
            TSIDESHAPE::NConnectShapeF(TSIDESHAPE::NSides - 1, h1SideOrder);
        //only for hexahedron
        TPZGenMatrix<int> shapeorders(nH1Internal,3);
        if constexpr(TSIDESHAPE::Type() == ECube){
            //not really used since we are interested in internal funcs
            TSHAPE::SideShapeOrder(TSIDESHAPE::NSides-1,nodeIds,
                                   h1SideOrder,shapeorders);
        }
        const auto xVecIndex = firstInternalVec + 0;
        const auto yVecIndex = firstInternalVec + 1;
        const auto zVecIndex = firstInternalVec + 2;

        auto addFunc = [&indexVecShape,&shapeCount,&shapeCountVec,iCon](
            int vIndex, int sIndex){
            indexVecShape[shapeCount] = std::make_pair(vIndex, sIndex);
            shapeCount++;
            shapeCountVec[iCon]++;
        };
        
        for(auto iFunc = 0; iFunc < nH1Internal; iFunc++ ){
          //it should alternate between them
          const auto shapeIndex = firstH1ShapeFunc[iCon] + iFunc;
          if constexpr(TSIDESHAPE::Type() == ECube)
            {
              if(shapeorders(iFunc,0) <= sideOrder){
                  addFunc(xVecIndex,shapeIndex);
              }
              if(shapeorders(iFunc,1) <= sideOrder){
                  addFunc(yVecIndex,shapeIndex);
              }
              if(shapeorders(iFunc,2) <= sideOrder){
                  addFunc(zVecIndex,shapeIndex);
              }
            }
          else{
              addFunc(xVecIndex,shapeIndex);
              addFunc(yVecIndex,shapeIndex);
              addFunc(zVecIndex,shapeIndex);
          }
        }

#ifdef PZ_LOG2
        if (logger.isDebugEnabled()) {
            const auto nInternalFuncs = shapeCount - firstInternalShape;
            const auto nKiFuncs = nInternalFuncs - nKfFuncs;
            std::ostringstream sout;
            sout << "n shape funcs (internal connect): "
                 << shapeCount - firstInternalShape << '\n'
                 << "\t n kf funcs : " << nKfFuncs
                 << "\t n ki funcs : " << nKiFuncs << '\n';
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
    }
}


template
struct TPZShapeHCurl<pzshape::TPZShapeTriang>;

template
struct TPZShapeHCurl<pzshape::TPZShapeQuad>;

template
struct TPZShapeHCurl<pzshape::TPZShapeTetra>;

template
struct TPZShapeHCurl<pzshape::TPZShapeCube>;

template
struct TPZShapeHCurl<pzshape::TPZShapePrism>;

//#define IMPLEMENTHCURLFULL(TSHAPE) \
//\
//template void TPZShapeHCurl<TSHAPE>::StaticIndexShapeToVec<pzshape::TPZShapeLinear>(TPZVec<std::pair<int,int64_t>> & indexVecShape,\
//const TPZVec<int>& connectOrder,\
//const TPZVec<int64_t>& firstH1ShapeFunc,const TPZVec<int>& sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,\
//const TPZVec<int64_t>& nodeIds);\
//template void TPZShapeHCurl<TSHAPE>::StaticIndexShapeToVec<pzshape::TPZShapeTriang>(TPZVec<std::pair<int,int64_t>> & indexVecShape,\
//const TPZVec<int>& connectOrder,\
//const TPZVec<int64_t>& firstH1ShapeFunc,const TPZVec<int>& sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,\
//const TPZVec<int64_t>& nodeIds);\
//template void TPZShapeHCurl<TSHAPE>::StaticIndexShapeToVec<pzshape::TPZShapeQuad>(TPZVec<std::pair<int,int64_t>> & indexVecShape,\
//const TPZVec<int>& connectOrder,\
//const TPZVec<int64_t>& firstH1ShapeFunc,const TPZVec<int>& sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,\
//const TPZVec<int64_t>& nodeIds);\
//
//
//IMPLEMENTHCURLFULL(pzshape::TPZShapeLinear)
//IMPLEMENTHCURLFULL(pzshape::TPZShapeTriang)
//IMPLEMENTHCURLFULL(pzshape::TPZShapeQuad)
//IMPLEMENTHCURLFULL(pzshape::TPZShapeCube)
//IMPLEMENTHCURLFULL(pzshape::TPZShapeTetra)
//IMPLEMENTHCURLFULL(pzshape::TPZShapePrism)
//
//#undef IMPLEMENTHCURLFULL

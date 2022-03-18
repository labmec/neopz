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
void TPZShapeHCurl<TSHAPE>::Initialize(const TPZVec<int64_t> &ids,
                                       TPZVec<int> &connectorders,
                                       TPZShapeData &data)
{
    
        
    constexpr int ncon = TSHAPE::NSides-TSHAPE::NCornerNodes;
    constexpr int NCorners = TSHAPE::NCornerNodes;
    constexpr int NSides = TSHAPE::NSides;
    if(connectorders.size() != ncon) DebugStop();
    CalcH1ShapeOrders(connectorders, data.fH1ConnectOrders);
    TPZShapeH1<TSHAPE>::Initialize(ids, data.fH1ConnectOrders, data);
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
        data.fHDivNumConnectShape[i] = ComputeNConnectShapeF(i,connectorders[i]);
        nShape += data.fHDivNumConnectShape[i];
    }
    
    data.fSDVecShapeIndex.Resize(nShape);
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
    data.fSDVecShapeIndex.Resize(nshape);

    TPZVec<unsigned int> shapeCountVec(TSHAPE::NSides - nNodes, 0);
    TPZVec<std::pair<int,int64_t>> & indexVecShape = data.fSDVecShapeIndex;
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
    StaticIndexShapeToVec(connOrder, firstH1ShapeFunc, sidesH1Ord, nodeIds, shapeCountVec, data.fSDVecShapeIndex);
}

template<class TSHAPE>
int TPZShapeHCurl<TSHAPE>::NHCurlShapeF(const TPZShapeData &data)
{
    int nshape = 0;
    int nc = data.fHDivNumConnectShape.size();
    for(int ic = 0; ic<nc; ic++) nshape += data.fHDivNumConnectShape[ic];
    return nshape;
}

    

template<class TSHAPE>
void TPZShapeHCurl<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &curlphi)
{

    constexpr int ncorner = TSHAPE::NCornerNodes;
    constexpr int nsides = TSHAPE::NSides;
    constexpr int dim = TSHAPE::Dimension;
    constexpr int curldim = [dim](){
        if constexpr (dim == 1) return 1;
        else{
            return 2*dim - 3;//1 for 2D 3 for 3D
        }
    }();
    
    TPZShapeH1<TSHAPE>::Shape(pt,data);
    for(int i = 0; i< data.fSDVecShapeIndex.size(); i++)
    {
        const auto &it = data.fSDVecShapeIndex[i];
        const int vecindex = it.first;
        const int scalindex = it.second;
        
        for(int d = 0; d<TSHAPE::Dimension; d++)
        {
            phi(d,i) = data.fPhi(scalindex,0)*data.fMasterDirections(d,vecindex);
        }

        if constexpr (dim==1){
            curlphi(0,i) =
                data.fDPhi.GetVal( 0,vecindex) *
                data.fMasterDirections.GetVal(0,vecindex);
        }else if constexpr (dim==2){
            curlphi(0,i) =
                data.fDPhi.GetVal(0,scalindex) *
                data.fMasterDirections.GetVal(1,vecindex) -
                data.fDPhi.GetVal(1,scalindex) *
                data.fMasterDirections.GetVal(0,vecindex);
            }
        else if constexpr(dim==3){
            for(auto d = 0; d < dim; d++) {
                const auto di = (d+1)%dim;
                const auto dj = (d+2)%dim;
                curlphi(d,i) =
                    data.fDPhi.GetVal(di,scalindex) *
                    data.fMasterDirections.GetVal(dj,vecindex)-
                    data.fDPhi.GetVal(dj,scalindex) *
                    data.fMasterDirections.GetVal(di,vecindex);
            }
        }else{
            if constexpr (std::is_same_v<TSHAPE,TSHAPE>){
                static_assert(!sizeof(TSHAPE),"Invalid curl dimension");
            }
        }        
    }
}


template<class TSHAPE>
int TPZShapeHCurl<TSHAPE>::ComputeNConnectShapeF(const int icon, const int order)

{
    constexpr int nNodes = TSHAPE::NCornerNodes;
    const int side = icon + nNodes;
#ifdef PZDEBUG
    if (side < nNodes || side >= TSHAPE::NSides) {
        DebugStop();
    }
#endif
    //given the choice of implementation, there are no shape functions for k=0
    if(order == 0) return 0;
    const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    const int nShapeF = [&](){
        if (side < nNodes + nEdges) {//edge connect
            return 1 + order;
        }
        else if(side < nNodes + nEdges + nFaces){//face connect
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
            if constexpr (TSHAPE::Type() == ETetraedro){
                return (order-1)*(order-2)*(order+1)/2;
            }
            else if constexpr (TSHAPE::Type() == ECube){
                return 3*order*order*(order+1);
            }
            else if constexpr (TSHAPE::Type() == EPrisma){
                return 3*order*(order-1)*(order+1)/2;
            }
            else{
                PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
                DebugStop();
                return 0;
            }
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
void TPZShapeHCurl<TSHAPE>::CalcH1ShapeOrders(
    const TPZVec<int> &ordHCurl,
    TPZVec<int> &ordH1)
{
    constexpr auto nConnects = TSHAPE::NSides-TSHAPE::NCornerNodes;
    if(ordH1.size() != nConnects)  ordH1.Resize(nConnects,-1);
    for(auto iCon = 0; iCon < nConnects; iCon++){
        const auto iSide = iCon + TSHAPE::NCornerNodes;
        const auto sideDim = TSHAPE::SideDimension(iSide);
        const bool quadSide =
            TSHAPE::Type(iSide) == EQuadrilateral ||
            TSHAPE::Type(iSide) == ECube ||
            TSHAPE::Type(iSide) == EPrisma;
        /*some H1 functions associated with the side iSide of dimension dim
          might be needed for computing the shape functions of a side with
          dimension dim+1 that contains the side iSide.
          It is also worth noting that quadrilateral sides require functions
          of ordH1er k+1*/
        TPZStack<int> highDimSides;
        TSHAPE::HigherDimensionSides(iSide, highDimSides);
        const auto sideOrder = ordHCurl[iCon];
        auto maxOrder = quadSide ? sideOrder + 1: sideOrder;
        for(auto &ihSide : highDimSides){
            if(TSHAPE::SideDimension(ihSide) != sideDim+1) break;
            else {
                const auto hSideOrder = ordHCurl[ihSide-TSHAPE::NCornerNodes];
                const auto hQuadSide =
                    TSHAPE::Type(ihSide) == EQuadrilateral ||
                    TSHAPE::Type(ihSide) == ECube ||
                    TSHAPE::Type(ihSide) == EPrisma;
                const auto hMaxOrder = hQuadSide ? hSideOrder + 1 : hSideOrder;
                maxOrder = std::max(maxOrder, hMaxOrder);
            }
        }
        ordH1[iCon] = maxOrder;
    }
}

template<class TSHAPE>
void TPZShapeHCurl<TSHAPE>::StaticIndexShapeToVec(const TPZVec<int>& connectOrder,
                                                  const TPZVec<int64_t>& firstH1ShapeFunc,
                                                  const TPZVec<int>& sidesH1Ord,
                                                  const TPZVec<int64_t>& nodeIds,
                                                  TPZVec<unsigned int>& shapeCountVec,
                                                  TPZVec<std::pair<int,int64_t>> & indexVecShape
                                                  ) {



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

    
    const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    constexpr auto nNodes = TSHAPE::NCornerNodes;

    const auto nConnects = connectOrder.size();
    TPZManVector<int, TSHAPE::NSides - nNodes>
        transformationIds(TSHAPE::NSides - nNodes, -1);
    //computing transformation id for sides.
    for (auto iCon = 0; iCon < nConnects; iCon++) {
        transformationIds[iCon] = TSHAPE::GetTransformId(nNodes + iCon, nodeIds);
    }
    
    unsigned int shapeCount = 0;


    //calculates edge functions
    for (auto iCon = 0; iCon < nEdges; iCon++) {
        const auto pOrder = connectOrder[iCon];
        //there will be 2 + pOrder - 1 = pOrder + 1 functions for each edge

        for (auto iNode = 0; iNode < 2; iNode++) {
            auto whichNode = transformationIds[iCon] == 0 ? iNode : (iNode + 1) % 2;
            const int vecIndex = iCon * 2 + whichNode;
            const int64_t shapeIndex = TSHAPE::SideNodeLocId(iCon + nNodes, whichNode);
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
    if(TSHAPE::Dimension < 2) return;
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
            TSHAPE::LowerDimensionSides(iFace + nEdges + nNodes, faceEdges[iFace], 1);
            const int nFaceEdges = iFace == 0 ? 0 : faceEdges[iFace-1].size();
            firstVfeVec[iFace] = iFace == 0 ?
                firstVfeVec[iFace] : firstVfeVec[iFace - 1] + nFaceEdges;
        }
    }
    const int firstVftVec = firstVfeVec[nFaces-1] + faceEdges[nFaces-1].size();
    const int firstVfOrthVec = firstVftVec + 2 * nFaces;

    const int nH1Funcs = TSHAPE::NShapeF(sidesH1Ord);
    TPZGenMatrix<int> shapeorders(nH1Funcs,3);
    TSHAPE::ShapeOrder(nodeIds,sidesH1Ord, shapeorders);
    
    for(auto iCon = nEdges; iCon < nEdges + nFaces; iCon++){
        const auto iSide = iCon + nNodes;
        const auto iFace = iCon - nEdges;

        int h1FaceOrder = -1;
        TPZManVector<int,4> permutedSideSides(4,-1);
        
        switch(TSHAPE::Type(iSide)){
        case ETriangle://triangular face
            pztopology::GetPermutation<pztopology::TPZTriangle>(transformationIds[iCon],
                                                                permutedSideSides);
            h1FaceOrder = connectOrder[iCon];
            break;
        case EQuadrilateral://quadrilateral face
            pztopology::GetPermutation<pztopology::TPZQuadrilateral>(transformationIds[iCon],
                                                                     permutedSideSides);
            h1FaceOrder = connectOrder[iCon]+1;
            break;
        default:
            PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
            DebugStop();
        }
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


        const int nFaceNodes = TSHAPE::NSideNodes(iSide);
        //this is not a mistake, since for faces nEdges = nNodes
        const int &nFaceEdges = nFaceNodes;

        //first the phi Fe functions
        for(auto iEdge = 0; iEdge < nFaceEdges; iEdge++ ){
            const auto currentLocalEdge = permutedSideSides[iEdge+nFaceNodes];
            const auto currentEdge = TSHAPE::ContainedSideLocId(iSide, currentLocalEdge);
            const auto vecIndex = firstVfeVec[iFace] + currentLocalEdge - nFaceNodes;
            
            for(auto iEdgeInternal = 0; iEdgeInternal < h1FaceOrder - 1; iEdgeInternal++){
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
        
        //number of h1 face funcs
        const auto nH1FaceFuncs =
            TSHAPE::NConnectShapeF(iSide,h1FaceOrder);
        
        const auto quadFace = TSHAPE::Type(iSide) == EQuadrilateral;
        
        if(quadFace){
            const auto transid = transformationIds[iCon];
            /*
              transformation id is 0 or 1 if starting at node 0
              transformation id is 2 or 3 if starting at node 1
              transformation id is 4 or 5 if starting at node 2
              transformation id is 6 or 7 if starting at node 3
              even numbers are counterclockwise dir,
              odd numbers are clockwise direction
              so ids 0, 3, 4, 7 have x as first dir
              and 1, 2, 5, 6 have y as first dir.

              we know that the vft vectors are created
              such that the first one is in x and the second one
              is in y (local face coordinates).
              therefore, we need the variables in the h1 functions
              to match.
              warning: sketchy integer arithmetic will follow:
             */
            const auto xdir = ((transid+1)/2)%2;//i told you so
            const auto ydir = 1-xdir;

            const auto hCurlFaceOrder = h1FaceOrder-1;
            const auto nfuncsk = 2 * (hCurlFaceOrder - 1) * (hCurlFaceOrder - 1);
            const auto nfuncsk1 = hCurlFaceOrder - 1;
            TPZVec<std::pair<int,int>> funcXY(nfuncsk);
            TPZVec<std::pair<int,int>> funcX(nfuncsk1);
            TPZVec<std::pair<int,int>> funcY(nfuncsk1);
            int countxy{0}, countx{0},county{0};
            /**now we assume that the first vft vec is in the x direction
               and that the next one is in the y direction.*/
            const int vecindex[] = {firstVftVec + 2*iFace,firstVftVec + 2*iFace+1};
            for(auto iFunc = 0; iFunc < nH1FaceFuncs; iFunc++ ){
                const auto shapeIndex = firstH1ShapeFunc[iCon] + iFunc;

                //functions of degree k
                if((shapeorders(shapeIndex,xdir) <= hCurlFaceOrder) &&
                   (shapeorders(shapeIndex,ydir) <= hCurlFaceOrder)){
                    funcXY[countxy++] = {vecindex[0],shapeIndex};
                    funcXY[countxy++] = {vecindex[1],shapeIndex};
                }else if(shapeorders(shapeIndex,xdir) <= hCurlFaceOrder){
                    funcX[countx++] = {vecindex[0],shapeIndex};
                }else if(shapeorders(shapeIndex,ydir) <= hCurlFaceOrder){
                    funcY[county++] = {vecindex[1],shapeIndex};
                }
            }

            auto AddFromVec = [&indexVecShape,&shapeCountVec, &shapeCount,iCon]
                (TPZVec<std::pair<int,int>> myvec){
                for(auto [vi,si] : myvec){
                    indexVecShape[shapeCount] = std::make_pair(vi,si);
                    shapeCount++;
                    shapeCountVec[iCon]++;
                }
            };
            AddFromVec(funcXY);
            AddFromVec(funcX);
            AddFromVec(funcY);
            
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

    if(TSHAPE::Dimension < 3) return;

    const int firstInternalShape = shapeCount;
    
    const int iCon = nEdges + nFaces;
    //hcurl connect order
    const auto sideOrder = connectOrder[iCon];
    //first, the phi KF functions
    for(int iFace = 0; iFace < nFaces; iFace++){
        const auto faceSide = iFace + nEdges + nNodes;
        const auto faceType = TSHAPE::Type(faceSide);
        const auto faceDim = TSHAPE::SideDimension(faceSide);
        const auto faceOrderH1 =
            TSHAPE::Type() == EPrisma || TSHAPE::Type() == ECube ?
            sideOrder + 1 : sideOrder;
        const auto nH1FaceFuncs =
            TSHAPE::NConnectShapeF(faceSide,faceOrderH1);
        const auto vecIndex = firstVfOrthVec + iFace;

        for(auto iFunc = 0; iFunc < nH1FaceFuncs; iFunc++ ){
            const auto shapeIndex = firstH1ShapeFunc[nEdges + iFace] + iFunc;

            bool skip = false;
            if constexpr(TSHAPE::Type() == EPrisma){
                const int ordvec[] =
                    {shapeorders(shapeIndex,0),shapeorders(shapeIndex,1)};
                skip = true;
                if((TSHAPE::Type(faceSide) == EQuadrilateral) &&
                   (ordvec[0] <= sideOrder) &&
                   (ordvec[1] <= sideOrder+1)){
                    //for quad faces ordvec[0] = xord and ordvec[1] = zord
                    skip = false;
                    
                }
                else if(TSHAPE::Type(faceSide) == ETriangle){
                    skip = false;
                }
            }
            if(!skip){
                indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
                shapeCount++;
                shapeCountVec[iCon]++;
            }
        }   
    }

    const auto nKfFuncs = shapeCount - firstInternalShape;

    const auto h1InternalOrd = sidesH1Ord[nEdges+nFaces];
    //now the phi Ki funcs
    const int firstInternalVec = firstVfOrthVec + nFaces;
    //ALL H1 internal functions
    const auto nH1Internal =
        TSHAPE::NConnectShapeF(TSHAPE::NSides - 1, h1InternalOrd);

    const auto xVecIndex = firstInternalVec + 0;
    const auto yVecIndex = firstInternalVec + 1;
    const auto zVecIndex = firstInternalVec + 2;

    auto addFunc = [&indexVecShape,&shapeCount,&shapeCountVec,iCon](
      int vIndex, int sIndex){
        indexVecShape[shapeCount] = std::make_pair(vIndex, sIndex);
        shapeCount++;
        shapeCountVec[iCon]++;
    };
    /**in order to filter the gradients of hcurl functions,
       we need to sort the shape functions for the hexahedral el*/
    auto AddToVec = [](TPZVec<std::pair<int,int>> &v,std::pair<int,int> vs){
        const auto vi = v.size();
        v.Resize(vi+1);
        v[vi] = vs;
    };

    TPZVec<std::pair<int,int>> funcXYZ, funcX, funcY, funcZ;
    
    
    for(auto iFunc = 0; iFunc < nH1Internal; iFunc++ ){
        const auto shapeIndex = firstH1ShapeFunc[iCon] + iFunc;
        const int xord = shapeorders(shapeIndex,0);
        const int yord = shapeorders(shapeIndex,1);
        const int zord = shapeorders(shapeIndex,2);
        if constexpr(TSHAPE::Type() == ECube){
            if(xord <= sideOrder && yord <= sideOrder && zord <= sideOrder){
                AddToVec(funcXYZ, std::make_pair(xVecIndex,shapeIndex));
                AddToVec(funcXYZ, std::make_pair(yVecIndex,shapeIndex));
                AddToVec(funcXYZ, std::make_pair(zVecIndex,shapeIndex));
            }
            else{
                if(xord <= sideOrder){
                    AddToVec(funcX, std::make_pair(xVecIndex,shapeIndex));
                }
                if(yord <= sideOrder){
                    AddToVec(funcY, std::make_pair(yVecIndex,shapeIndex));
                }
                if(zord <= sideOrder){
                    AddToVec(funcZ, std::make_pair(zVecIndex,shapeIndex));
                }
            }
        }
        else if constexpr(TSHAPE::Type() == EPrisma){
            if((xord <= sideOrder) &&
               (yord <= sideOrder) &&
               (zord <= sideOrder+1)){
                addFunc(xVecIndex,shapeIndex);
                addFunc(yVecIndex,shapeIndex);
            }
            if((xord <= sideOrder+1) &&
               (yord <= sideOrder+1) &&
               (zord <= sideOrder)){
                addFunc(zVecIndex,shapeIndex);
            }
        }
        else{
            addFunc(xVecIndex,shapeIndex);
            addFunc(yVecIndex,shapeIndex);
            addFunc(zVecIndex,shapeIndex);
        }
    }
    //now we actually add the functions
    if constexpr(TSHAPE::Type() == ECube){
        for(auto [v,s] : funcXYZ){
            addFunc(v,s);
        }
        for(auto [v,s] : funcX){
            addFunc(v,s);
        }
        for(auto [v,s] : funcY){
            addFunc(v,s);
        }
        for(auto [v,s] : funcZ){
            addFunc(v,s);
        }
    }
    
    if(shapeCount != indexVecShape.size()){
        DebugStop();
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


template<class TSHAPE>
int TPZShapeHCurl<TSHAPE>::MaxOrder(const int ordh1){

    if constexpr (std::is_same_v<TSHAPE,pzshape::TPZShapeCube> ||
                  std::is_same_v<TSHAPE,pzshape::TPZShapeQuad>){
        return ordh1+1;
    }else{
        return ordh1;
    }
}

template
struct TPZShapeHCurl<pzshape::TPZShapeLinear>;

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

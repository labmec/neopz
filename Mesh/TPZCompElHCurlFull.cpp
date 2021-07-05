/**
 * @file
 * @brief Contains the implementation of the TPZCompElHCurlmethods.
 */
#include <TPZCompElHCurlFull.h>
#include <pzshapetetra.h>
#include <pzshapequad.h>
#include <pzshapetriang.h>
#include <pzshapelinear.h>

#include "pzcmesh.h"
#include "TPZTopologyUtils.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "pzgenericshape.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHCurl");
#endif

template<class TSHAPE>
TPZCompElHCurlFull<TSHAPE>::TPZCompElHCurlFull(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElHCurlFull::ClassId),
TPZCompElHCurl<TSHAPE>(mesh,gel,index){
    TPZCompElHCurl<TSHAPE>::CreateHCurlConnects(mesh);
}

template<class TSHAPE>
TPZCompElHCurlFull<TSHAPE>::TPZCompElHCurlFull(TPZCompMesh &mesh, const TPZCompElHCurlFull<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHCurlFull::ClassId),
TPZCompElHCurl<TSHAPE>(mesh,copy){
}

template<class TSHAPE>
TPZCompElHCurlFull<TSHAPE>::TPZCompElHCurlFull(TPZCompMesh &mesh,
									 const TPZCompElHCurlFull<TSHAPE> &copy,
									 std::map<int64_t,int64_t> & gl2lcConMap,
									 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHCurlFull::ClassId),
TPZCompElHCurl<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap){
}

template<class TSHAPE>
TPZCompElHCurlFull<TSHAPE>::TPZCompElHCurlFull() :
TPZRegisterClassId(&TPZCompElHCurlFull::ClassId),
TPZCompElHCurl<TSHAPE>(){
}

template<class TSHAPE>
int TPZCompElHCurlFull<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHCurlFull") ^ TPZCompElHCurl<TSHAPE>::ClassId() << 1;
}

template<class TSHAPE>
int TPZCompElHCurlFull<TSHAPE>::NConnectShapeF(int icon, int order) const {
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

#ifdef PZ_LOG
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
void TPZCompElHCurlFull<TSHAPE>::IndexShapeToVec(TPZVec<std::pair<int,int64_t>> & indexVecShape, const TPZVec<int>& connectOrder) const {
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

#ifdef PZ_LOG
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
    {
        for (auto iNode = 0; iNode < nNodes; iNode++) nodes[iNode] = this->Reference()->NodeIndex(iNode);
    }

    TPZManVector<int64_t, TSHAPE::NSides - nNodes> firstH1ShapeFunc(TSHAPE::NSides - nNodes,
                                                                                  0);
    //calculates the first shape function associated with each side of dim > 0
    TPZManVector<int,TSHAPE::NSides-nNodes> sidesH1Ord(TSHAPE::NSides - nNodes,-1);
    this->CalcH1ShapeOrders(sidesH1Ord);
    firstH1ShapeFunc[0] = nNodes;
    for (int iSide = nNodes + 1; iSide < TSHAPE::NSides; iSide++) {
        const int iCon = iSide - nNodes;
        firstH1ShapeFunc[iCon] = firstH1ShapeFunc[iCon - 1] + TSHAPE::NConnectShapeF(iSide - 1, sidesH1Ord[iCon-1]);
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::ostringstream sout;
        sout << "first H1 shape function:" << std::endl;
        for (auto &iShape : firstH1ShapeFunc) {
            sout << "\t" << iShape << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif


    indexVecShape.Resize(this->NShapeF());
    TPZVec<unsigned int> shapeCountVec(TSHAPE::NSides - nNodes, 0);

    StaticIndexShapeToVec(indexVecShape, connectOrder, firstH1ShapeFunc,sidesH1Ord, shapeCountVec, nodes);



#ifdef PZDEBUG
    int nFuncs = 0;
    for (int iEdge = 0; iEdge < nEdges; iEdge++){
        nFuncs += shapeCountVec[iEdge];
        if (shapeCountVec[iEdge] != this->NConnectShapeF(iEdge, connectOrder[iEdge])) {
            std::ostringstream soutAbort;
            soutAbort << __PRETTY_FUNCTION__ << std::endl;
            soutAbort << "\tError with the number of shape functions of edge " << iEdge + nNodes << std::endl;
            soutAbort << "\tCalculated " << shapeCountVec[iEdge] << " instead of "
                      << NConnectShapeF(iEdge, connectOrder[iEdge]) << std::endl;
#ifdef PZ_LOG
            LOGPZ_ERROR(logger, soutAbort.str())
#endif
            PZError << soutAbort.str() << std::endl;
        }
    }

    for (int iFace = 0; iFace < nFaces; iFace++){
        nFuncs += shapeCountVec[iFace + nEdges];
        if (shapeCountVec[iFace + nEdges] != this->NConnectShapeF(iFace + nEdges, connectOrder[iFace + nEdges])) {
            std::ostringstream soutAbort;
            soutAbort << __PRETTY_FUNCTION__ << std::endl;
            soutAbort << "\tError with the number of shape functions of face " << iFace + nEdges + nNodes << std::endl;
            soutAbort << "\tCalculated " << shapeCountVec[iFace + nEdges] << " instead of "
                      << NConnectShapeF(iFace + nEdges, connectOrder[iFace + nEdges]) << std::endl;
#ifdef PZ_LOG
            LOGPZ_ERROR(logger, soutAbort.str())
#endif
            PZError << soutAbort.str() << std::endl;
        }
    }
    if(TSHAPE::Dimension < 3) return;
    const auto lastSide = TSHAPE::NSides - 1;
    if(indexVecShape.size() - nFuncs !=
    this->NConnectShapeF(lastSide - nNodes, connectOrder[lastSide - nNodes]) ){
        std::ostringstream soutAbort;
        soutAbort << __PRETTY_FUNCTION__ << std::endl;
        soutAbort << "\tError with the number of internal shape functions"<< std::endl;
        soutAbort << "\tCalculated " << indexVecShape.size() - nFuncs << " instead of "
                  << this->NConnectShapeF(lastSide - nNodes, connectOrder[lastSide - nNodes]) << std::endl;
#ifdef PZ_LOG
        LOGPZ_ERROR(logger, soutAbort.str())
#endif
        PZError << soutAbort.str() << std::endl;
    }
#endif

}

template<class TSHAPE>
template<class TSIDESHAPE>
void TPZCompElHCurlFull<TSHAPE>::StaticIndexShapeToVec(TPZVec<std::pair<int,int64_t>> & indexVecShape,
                                                       const TPZVec<int>& connectOrder,
                                                       const TPZVec<int64_t>& firstH1ShapeFunc,
                                                       const TPZVec<int>& sidesH1Ord,
                                                       TPZVec<unsigned int>& shapeCountVec,
                                                       const TPZVec<int64_t>& nodeIds) {
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

#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::ostringstream sout;
        sout << "vec shape index (edge connects):" << std::endl;
        for (int iShape = 0; iShape < firstFaceShape; iShape++) {
            auto pair = indexVecShape[iShape];
            sout << "\tvec: " << pair.first << "\tshape: " << pair.second << std::endl;

        }
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
#ifdef PZ_LOG
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
#ifdef PZ_LOG
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

#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::ostringstream sout;
        sout << "vec shape index (face connects):" << std::endl;
        for (int iShape = firstFaceShape; iShape < firstInternalShape; iShape++) {
            auto pair = indexVecShape[iShape];
            sout << "\tvec: " << pair.first << "\tshape: " << pair.second << std::endl;

        }
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
    }
}


template<class TSHAPE>
void TPZCompElHCurlFull<TSHAPE>::CalcH1ShapeOrders(TPZVec<int> &ord) const{
    const auto nConnects = this->NConnects();
    TPZVec<int> hcurlOrders(nConnects);
    for(auto i = 0; i < nConnects; i++){
        hcurlOrders[i] = this->EffectiveSideOrder(TSHAPE::NCornerNodes + i);
    }
    StaticCalcH1ShapeOrders<TSHAPE>(hcurlOrders, ord);
}

template<class TSHAPE>
template<class TSIDESHAPE>
void TPZCompElHCurlFull<TSHAPE>::StaticCalcH1ShapeOrders(
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
void TPZCompElHCurlFull<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
    const int dim = TSHAPE::Dimension;
    const int sideDim = TSHAPE::SideDimension(side);
    const MElementType sidetype = TSHAPE::Type(side);
#ifdef PZDEBUG
    if(side >= TSHAPE::NSides || side < TSHAPE::NCornerNodes){
        PZError<<__PRETTY_FUNCTION__<<"\n";
        PZError<<"There is no side shape associated to this side"<<"\n";
        DebugStop();
    }
	if( sideDim != point.size() ){
        PZError<<__PRETTY_FUNCTION__<<"\n";
        PZError<<"Wrong dimension of point! point dim: "<<point.size()<<"\tside dim: "<<sideDim<<"\n";
		return ;
	}
#endif

    const int connectLocalId = this->MidSideConnectLocId(side);
    const int connectOrder = this->Connect(connectLocalId).Order();
    const int nContainedSides = TSHAPE::NContainedSides(side);
    const int nSideNodes = TSHAPE::NSideNodes(side);
    const int nSideConnects = nContainedSides - nSideNodes;
    //vector with the connect order for each side of dim >=1 contained in side
    TPZVec<int> ordHCurl(nSideConnects);
    for (auto is=nSideNodes; is< nContainedSides; is++) {
        const int subSide = TSHAPE::ContainedSideLocId(side,is);
        ordHCurl[is-nSideNodes] = this->EffectiveSideOrder(subSide);
    }
    //get h1 connect orders
    TPZVec<int> ordH1(nSideConnects);
    switch(sidetype){
    case EOned:
        StaticCalcH1ShapeOrders<pztopology::TPZLine>(ordHCurl,
                                                            ordH1);
        break;
    case ETriangle:
        StaticCalcH1ShapeOrders<pztopology::TPZTriangle>(ordHCurl,
                                                                ordH1);
        break;
    case EQuadrilateral:
        StaticCalcH1ShapeOrders<pztopology::TPZQuadrilateral>(ordHCurl,
                                                                     ordH1);
        break;
    default:
        PZError<<__PRETTY_FUNCTION__
               <<"\n invalid side type.Aborting...\n";
        DebugStop();
    }
    //number of HCURL shape functions
    const int nSideShapes = [&]{

        int nShapes = 0;
        for (auto is=nSideNodes; is< nContainedSides; is++) {
            const int subSide = TSHAPE::ContainedSideLocId(side,is);
            const int subSideDim = TSHAPE::SideDimension(subSide);
            if(subSideDim < 1) continue;
            const int subConnectLocalId = this->MidSideConnectLocId(subSide);
            const int subConnectOrder = ordHCurl[is-nSideNodes];
            nShapes += this->NConnectShapeF(subConnectLocalId,subConnectOrder);
        }
        return nShapes;
    }();
#ifdef PZDEBUG
    if (nSideShapes != this->NSideShapeF(side)) {
        DebugStop();
    }
#endif
    //get the ids of the side nodes
    TPZGeoEl *gel = this->Reference();
    TPZManVector<int64_t,8> sideNodesId(nSideNodes);
    TPZManVector<int, 5> transformationIds(nSideConnects, -1);
    for (int ic=0; ic< nSideNodes; ic++) {
        const int localId = TSHAPE::SideNodeLocId(side,ic);
        sideNodesId[ic] = gel->Node(localId).Id();
    }
    for (auto iSubSide = nSideNodes; iSubSide < nContainedSides; iSubSide++) {
        transformationIds[iSubSide-nSideNodes] = [&](){
            switch(sidetype){
            case EOned:
                return pztopology::TPZLine::GetTransformId(iSubSide,
                                                           sideNodesId);
                break;
            case ETriangle:
                return pztopology::TPZTriangle::GetTransformId(iSubSide,
                                                               sideNodesId);
                break;
            case EQuadrilateral:
                return pztopology::TPZQuadrilateral::GetTransformId(iSubSide,
                                                                    sideNodesId);
                break;
            default:
                DebugStop();
                return -1;
            }
        }();
//        const int localId = TSHAPE::ContainedSideLocId(side,iSubSide);
//        transformationIds[iSubSide-nSideNodes] = TSHAPE::GetTransformId(localId, elNodes);
    }


    //calculates the directions on the master element associated with the side
    //and the indexes associating the side with the respective h1 scalar function
    TPZFMatrix<REAL> sideMasterDirections(sideDim,sideDim * nContainedSides,0);
    TPZManVector<std::pair<int,int64_t>> indexVecShape(nSideShapes);
    {
        MElementType sidetype = TSHAPE::Type(side);
        TPZManVector<unsigned int,5> shapeCountVec(nSideConnects,-1);
        TPZManVector<int64_t,5> firstH1ShapeFunc(nSideConnects,-1);

        //calculates the first SCALAR shape function associated with each
        // side of dim > 0
        firstH1ShapeFunc[0] = nSideNodes;
        for (int iSide = nSideNodes + 1; iSide < nContainedSides; iSide++) {
            const int prevLocalId = TSHAPE::ContainedSideLocId(side,iSide-1);
            const int &lastFirstH1 = firstH1ShapeFunc[iSide - nSideNodes - 1];
            const int nShapeF = TSHAPE::NConnectShapeF(prevLocalId, ordH1[iSide-nSideNodes-1]);
            firstH1ShapeFunc[iSide - nSideNodes] = lastFirstH1 + nShapeF;
//            firstH1ShapeFunc[iSide - nSideNodes] = firstH1ShapeFunc[iSide - nSideNodes - 1]
//                                                + TSHAPE::NConnectShapeF(prevCon, sidesH1Ord[iSide-nSideNodes-1]);
        }

#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::ostringstream sout;
            sout << __PRETTY_FUNCTION__ << std::endl;
            sout << "side :"<< side <<std::endl;
            sout << "transformation id :"<< transformationIds[nContainedSides-nSideNodes - 1] <<std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        //calculates hcurl dirs and vec shape index
        TPZFNMatrix<9,REAL> gradxSide(sideDim,sideDim,0);
        for(auto ix = 0; ix < sideDim; ix++) gradxSide(ix,ix) = 1;
        switch (sidetype) {
            case EOned://these wont be really used, just for signs purposes
                pztopology::TPZLine::ComputeHCurlDirections(gradxSide,sideMasterDirections,transformationIds);
                StaticIndexShapeToVec<pzshape::TPZShapeLinear>(indexVecShape,ordHCurl,firstH1ShapeFunc,ordH1,shapeCountVec,sideNodesId);
                break;
            case EQuadrilateral:
                pztopology::TPZQuadrilateral::ComputeHCurlDirections(gradxSide,sideMasterDirections,transformationIds);
                StaticIndexShapeToVec<pzshape::TPZShapeQuad>(indexVecShape,ordHCurl,firstH1ShapeFunc,ordH1,shapeCountVec,sideNodesId);
                break;
            case ETriangle:
                pztopology::TPZTriangle::ComputeHCurlDirections(gradxSide,sideMasterDirections,transformationIds);
                StaticIndexShapeToVec<pzshape::TPZShapeTriang>(indexVecShape,ordHCurl,firstH1ShapeFunc,ordH1,shapeCountVec,sideNodesId);
                break;
            default:
                DebugStop();
                break;
        }
    }

    //get the jacobian of the side transformation
    TPZGeoElSide gelside = TPZGeoElSide(this->Reference(),side);
    TPZFNMatrix<9,REAL> jac(sideDim,sideDim),jacinv(sideDim,sideDim),axes(sideDim,3);
    REAL detjac = 0;
    gelside.Jacobian(point, jac, axes, detjac, jacinv);

    //compute side scalar shape functions
    TPZFNMatrix<50,REAL> phiSide(nSideShapes,1),dPhiSide(sideDim,nSideShapes);
    TSHAPE::SideShape(side, point, sideNodesId, ordH1, phiSide, dPhiSide);
    
    const int phiDim = [&](){
        switch(sideDim){
            case 1: return 1;
            case 2: return 3;
            default:
                DebugStop();
                return -1;
        }
    }();
    //gets sideDeformedDirections
    TPZFMatrix<REAL> sideDeformedDirections(3,sideDim * nContainedSides,0);
    const int nVec = sideMasterDirections.Cols();
    for (auto iVec = 0; iVec < nVec; iVec++) {
        TPZManVector<REAL, 3> tempDirection(sideDim, 0);
        for (auto i = 0; i < sideDim; i++) {
            //covariant piola transform: J^{-T}
            tempDirection[i] = 0;
            for (auto j = 0; j < sideDim; j++) tempDirection[i] += jacinv(j, i) * sideMasterDirections(j, iVec);
        }
        for (auto i = 0; i < 3; i++) {
            sideDeformedDirections(i, iVec) = 0;
            for (auto j = 0; j < sideDim; j++) sideDeformedDirections(i, iVec) += axes(j, i) * tempDirection[j];
        }
    }

    //calculates phi
    phi.Redim(nSideShapes,3);
    TPZHCurlAuxClass::ComputeShape(indexVecShape, phiSide, sideDeformedDirections,phi);

    /*
     * The curl of a function express two different things in 2D/3D. In 2D, $curl\phi: \mathbb{R}^2 \rightarrow \mathbb{R}$,
     * while in 3D $curl\phi: \mathbb{R}^3 \rightarrow \mathbb{R}^3$. Therefore, they are computed in different ways
     * in this function. For 3D elements, the curl follows basically the same procedure as the shape functions:
     * It's 3D expression is computed and then its trace is evaluated (in the XYZ space as well).
     * For 2D elements, however, the curl is but a scalar value.
     * */
    const auto transfSide = TSHAPE::TransformElementToSide(side).Mult();
    TPZFMatrix<REAL> sideMasterDirectionsOnElement(dim,sideDim * nContainedSides,0);
    for (auto iVec = 0; iVec < nVec; iVec++) {
        TPZManVector<REAL, 3> tempDirection(sideDim, 0);
        for (auto i = 0; i < dim; i++) {
            sideMasterDirectionsOnElement(i, iVec) = 0;
            for (auto j = 0; j < sideDim; j++) sideMasterDirectionsOnElement(i, iVec) += transfSide.GetVal(j, i) * sideMasterDirections(j,iVec);
        }
    }
    TPZFMatrix<REAL>dPhiSideOnElement(dim,nSideShapes);
    {
        TPZManVector<int64_t,8> elNodesIds(TSHAPE::NCornerNodes);
        for (int ic=0; ic< TSHAPE::NCornerNodes; ic++) {
            elNodesIds[ic] = gel->Node(ic).Id();
        }
        const REAL alpha = 1.0, beta =0.0;
        const int transpose = 1;
        GetSideTransform<TSHAPE>(side,TSHAPE::GetTransformId(side,elNodesIds))
                .Mult().MultAdd(dPhiSide,dPhiSide,dPhiSideOnElement,1,0,1);
    }
    const int curlDim = sideDim ==1 ? 1 : 2* sideDim - 3;
    dphi.Resize(curlDim,nSideShapes);
    switch(sideDim){

    case 1:
      TPZHCurlAuxClass::ComputeCurl<1>(indexVecShape, dPhiSideOnElement,
                                       sideMasterDirectionsOnElement, jac,
                                       detjac, axes, dphi);
      break;
    case 2:
      TPZHCurlAuxClass::ComputeCurl<2>(indexVecShape, dPhiSideOnElement,
                                       sideMasterDirectionsOnElement, jac,
                                       detjac, axes, dphi);
      break;
    case 3:
      TPZHCurlAuxClass::ComputeCurl<3>(indexVecShape, dPhiSideOnElement,
                                       sideMasterDirectionsOnElement, jac,
                                       detjac, axes, dphi);
      break;
    default:
      DebugStop();
    }
}

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"

#define IMPLEMENTHCURLFULL(TSHAPE) \
\
template class \
TPZRestoreClass< TPZCompElHCurlFull<TSHAPE> >; \
template class TPZCompElHCurlFull<TSHAPE>; \
template void TPZCompElHCurlFull<TSHAPE>::StaticIndexShapeToVec<pzshape::TPZShapeLinear>(TPZVec<std::pair<int,int64_t>> & indexVecShape,\
const TPZVec<int>& connectOrder,\
const TPZVec<int64_t>& firstH1ShapeFunc,const TPZVec<int>& sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,\
const TPZVec<int64_t>& nodeIds);\
template void TPZCompElHCurlFull<TSHAPE>::StaticIndexShapeToVec<pzshape::TPZShapeTriang>(TPZVec<std::pair<int,int64_t>> & indexVecShape,\
const TPZVec<int>& connectOrder,\
const TPZVec<int64_t>& firstH1ShapeFunc,const TPZVec<int>& sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,\
const TPZVec<int64_t>& nodeIds);\
template void TPZCompElHCurlFull<TSHAPE>::StaticIndexShapeToVec<pzshape::TPZShapeQuad>(TPZVec<std::pair<int,int64_t>> & indexVecShape,\
const TPZVec<int>& connectOrder,\
const TPZVec<int64_t>& firstH1ShapeFunc,const TPZVec<int>& sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,\
const TPZVec<int64_t>& nodeIds);


IMPLEMENTHCURLFULL(pzshape::TPZShapeLinear)
IMPLEMENTHCURLFULL(pzshape::TPZShapeTriang)
IMPLEMENTHCURLFULL(pzshape::TPZShapeQuad)
IMPLEMENTHCURLFULL(pzshape::TPZShapeCube)
IMPLEMENTHCURLFULL(pzshape::TPZShapeTetra)
IMPLEMENTHCURLFULL(pzshape::TPZShapePrism)

#undef IMPLEMENTHCURLFULL

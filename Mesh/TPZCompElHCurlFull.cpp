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
    if(order == 0) return 0;//given the choice of implementation, there are no shape functions for k=0
    const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    const auto sideOrder = this->EffectiveSideOrder(side);
    const int nShapeF = [&](){
        if (side < TSHAPE::NCornerNodes + nEdges) {//edge connect
            return 1 + sideOrder;
        }
        else if(side < TSHAPE::NCornerNodes + nEdges + nFaces){//face connect
            const int factor = [&](){

                switch(TSHAPE::Type(side)){
                    case ETriangle://triangular face
                        return 1;
                    case EQuadrilateral://quadrilateral face
                        return 2;
                    default:
                        PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
                        DebugStop();
                        return 0;
                }
            }();
            return factor * (sideOrder - 1) * (sideOrder + 1);
        }
        else{//internal connect (3D element only)
            int count = 0;
            for(int iFace = 0; iFace < nFaces; iFace++){
                switch(TSHAPE::Type(TSHAPE::NCornerNodes+nEdges + iFace)){
                    case ETriangle://triangular face
                        count +=  (sideOrder - 1) * ( sideOrder - 2)/2;
                        break;
                    case EQuadrilateral://quadrilateral face
                        count +=  (sideOrder - 1) * ( sideOrder - 1);
                        break;
                    default:
                        PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
                        DebugStop();
                }
            }
            //number of H1 funcs
            const auto nFuncsK = TSHAPE::NConnectShapeF(side, sideOrder);
            if constexpr (TSHAPE::Type() == ECube){
              //make sure there are no negative vals
              const int nFuncsKminusOne =
                TSHAPE::NConnectShapeF(side, sideOrder-1) >= 0 ? 
                TSHAPE::NConnectShapeF(side, sideOrder-1) : 
                0;
              //include all internal functions of k-1 (3*(k-2)^2 funcs)
              count += 3 * (nFuncsKminusOne);
              //shapeorders[k] contains the polynomial orders
              //in xi, eta and zeta for each of the internal funcs
              TPZGenMatrix<int> shapeorders(nFuncsK,3);
              //not really used since we are interested in internal funcs
              TPZManVector<int64_t,0> idvec(0);
              TSHAPE::SideShapeOrder(side,idvec,sideOrder,shapeorders);
              /*for compatibility, we need the spaces
                Q^{k-1,k,k} \times Q^{k,k-1,k} \times Q^{k,k,k-1}
               */
              int countIntK = 0;
              /*i am supposing that #funcs in xi-dir is equal to 
                #funcs in eta-dir and in zeta-dir. So I will just
                see how many we will have in xi-dir afte rskipping
                internals that were included in the k-1 set of funcs*/
              for(int ifunc = nFuncsKminusOne; ifunc < nFuncsK; ifunc++){
                if(shapeorders(ifunc,0) < sideOrder) {countIntK++;}
              }
              count += 3* countIntK;
                
            }
            else{
              count += 3 * nFuncsK;
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
        sout<<"\torder "<<order<<"\teffective order "<<sideOrder<<std::endl;
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
    TPZManVector<int, TSHAPE::NSides - nNodes> transformationIds(TSHAPE::NSides - nNodes, -1);
    {
        TPZManVector<int64_t,nNodes> nodes(nNodes, 0);
        for (auto iNode = 0; iNode < nNodes; iNode++) nodes[iNode] = this->Reference()->NodeIndex(iNode);
        //computing transformation id for sides.
        for (auto iCon = 0; iCon < nConnects; iCon++) {
            transformationIds[iCon] = TSHAPE::GetTransformId(nNodes + iCon, nodes);
        }
    }

    TPZManVector<int64_t, TSHAPE::NSides - nNodes> firstH1ShapeFunc(TSHAPE::NSides - nNodes,
                                                                                  0);
    //calculates the first shape function associated with each side of dim > 0
    TPZManVector<int,TSHAPE::NSides> sidesH1Ord(TSHAPE::NSides);
    this->CalculateSideShapeOrders(sidesH1Ord);
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

    StaticIndexShapeToVec(indexVecShape, connectOrder, firstH1ShapeFunc,sidesH1Ord, shapeCountVec, transformationIds);


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
                                                       const TPZVec<int>& transformationIds) {
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
    const auto nNodes = TSIDESHAPE::NCornerNodes;
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
     * e) nFaceInternalFunctions number of internal SCALAR functions associated with each face
     */
    TPZManVector<int> firstVfeVec(nFaces,-1);
    TPZManVector<TPZStack<int>> faceEdges(nFaces,TPZStack<int>(0,0));
    TPZManVector<int> nFaceInternalFunctions(nFaces,-1);
    {
        const int nEdgeVectors = nEdges * 3;
        firstVfeVec[0] = nEdgeVectors;
        nFaceInternalFunctions[0] = TSIDESHAPE::NConnectShapeF(0 + nEdges + nNodes,connectOrder[nEdges + 0]);
        for(auto iFace = 1; iFace < nFaces; iFace++){
            TSIDESHAPE::LowerDimensionSides(iFace - 1 + nEdges + nNodes, faceEdges[iFace-1], 1);
            const int nFaceEdges = faceEdges[iFace-1].size();
            firstVfeVec[iFace] = firstVfeVec[iFace - 1] + nFaceEdges;
            nFaceInternalFunctions[iFace] = TSIDESHAPE::NConnectShapeF(iFace + nEdges + nNodes,connectOrder[nEdges + iFace]);
        }
        TSIDESHAPE::LowerDimensionSides(nFaces - 1 + nEdges + nNodes, faceEdges[nFaces-1], 1);
    }
    const int firstVftVec = firstVfeVec[nFaces-1] + faceEdges[nFaces-1].size();
    const int firstVfOrthVec = firstVftVec + 2 * nFaces;

    //first the phi Fe functions
    for(auto iCon = nEdges; iCon < nEdges + nFaces; iCon++){
        const auto iSide = iCon + nNodes;
        const auto iFace = iCon - nEdges;
        const auto pOrder = connectOrder[iCon];
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
        const int &nFaceEdges = nFaceNodes;//this is not a mistake, since for faces nEdges = nNodes
        for(auto iEdge = 0; iEdge < nFaceEdges; iEdge++ ){
            const auto currentLocalEdge = permutedSideSides[iEdge+nFaceNodes];
            const auto currentEdge = TSIDESHAPE::ContainedSideLocId(iSide, currentLocalEdge);
            const auto vecIndex = firstVfeVec[iFace] + currentLocalEdge - nFaceNodes;
            for(auto iEdgeInternal = 0; iEdgeInternal < pOrder - 1; iEdgeInternal++){
                const int shapeIndex = firstH1ShapeFunc[currentEdge - nNodes] + iEdgeInternal;
                indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
                shapeCount++;
                shapeCountVec[iCon]++;
            }
        }
        //now the phi Fi functions
        const auto nFaceInternalFuncs = 2 * nFaceInternalFunctions[iFace];
        for(auto iFunc = 0; iFunc < nFaceInternalFuncs; iFunc++ ){
            const auto vecIndex = firstVftVec + 2*iFace + iFunc % 2;//it should alternate between them
            const auto shapeIndex = firstH1ShapeFunc[iCon] + iFunc / 2;
            indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
            shapeCount++;
            shapeCountVec[iCon]++;
        }
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
            const auto nFaceInternalHCurlFuncs = TSIDESHAPE::NConnectShapeF(iFace + nEdges + nNodes,sidesH1Ord[nEdges + iFace]);
            const auto vecIndex = firstVfOrthVec + iFace;
            for(auto iFunc = 0; iFunc < nFaceInternalHCurlFuncs; iFunc++ ){
                const auto shapeIndex = firstH1ShapeFunc[nEdges + iFace] + iFunc;
                indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
                shapeCount++;
                shapeCountVec[iCon]++;
            }
        }
        const auto sideOrder = connectOrder[iCon];
        //now the phi Ki funcs
        const int firstInternalVec = firstVfOrthVec + nFaces;
        //ALL H1 internal functions
        const auto nH1Internal = TSIDESHAPE::NConnectShapeF(TSIDESHAPE::NSides - 1, sideOrder);
        const auto nInternalFuncs = 3 * nH1Internal;
        //only for hexahedron
        TPZGenMatrix<int> shapeorders(nH1Internal,3);
        if constexpr(TSIDESHAPE::Type() == ECube)
          {
            //not really used since we are interested in internal funcs
            TPZManVector<int64_t,0> idvec(0);
            TSHAPE::SideShapeOrder(TSIDESHAPE::NSides-1,idvec,
                                   sideOrder,shapeorders);
          }
        for(auto iFunc = 0; iFunc < nInternalFuncs; iFunc++ ){
          //it should alternate between them
          const auto dir = iFunc % 3;
          const auto vecIndex = firstInternalVec + dir;
          const auto intShapeIndex = iFunc/3;
          const auto shapeIndex = firstH1ShapeFunc[iCon] + intShapeIndex;
          if constexpr(TSIDESHAPE::Type() == ECube)
            {
              if(shapeorders(intShapeIndex,dir) >= sideOrder)
                {continue;}
            }
          indexVecShape[shapeCount] = std::make_pair(vecIndex, shapeIndex);
          shapeCount++;
          shapeCountVec[iCon]++;
        }
    }
}

template<class TSHAPE>
void TPZCompElHCurlFull<TSHAPE>::CalculateSideShapeOrders(TPZVec<int> &ord) const{
    constexpr auto nConnects = TSHAPE::NSides-TSHAPE::NCornerNodes;
    if(ord.size() != nConnects)  ord.Resize(nConnects,-1);
    for(auto iCon = 0; iCon < nConnects; iCon++){
        const auto iSide = iCon + TSHAPE::NCornerNodes;
        const auto sideDim = TSHAPE::SideDimension(iSide);
        /*some H1 functions associated with the side iSide of dimension dim might be needed for computing
        the shape functions of a side with dimension dim+1 that contains the side iSide*/
        TPZStack<int> highDimSides;
        TSHAPE::HigherDimensionSides(iSide, highDimSides);
        auto maxOrder = this->EffectiveSideOrder(iSide);
        for(auto &iHighSide : highDimSides){
            if(TSHAPE::SideDimension(iHighSide) != sideDim+1) break;
            else if(maxOrder < this->EffectiveSideOrder(iHighSide)){
                maxOrder = this->EffectiveSideOrder(iHighSide);
            }
        }
        ord[iCon] = maxOrder;
    }
}

template<class TSHAPE>
void TPZCompElHCurlFull<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
    const int sideDim = TSHAPE::SideDimension(side);
    const int dim = TSHAPE::Dimension;
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
    const int nSideShapes = [&]{

        int nShapes = 0;
        for (auto is=0; is< nContainedSides; is++) {
            const int subSide = TSHAPE::ContainedSideLocId(side,is);
            const int subSideDim = TSHAPE::SideDimension(subSide);
            if(subSideDim < 1) continue;
            const int subConnectLocalId = this->MidSideConnectLocId(subSide);
            const int subConnectOrder = this->Connect(subConnectLocalId).Order();
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
    const int  nSideNodes = TSHAPE::NSideNodes(side);
    const int nSideConnects = nContainedSides - nSideNodes;
    TPZGeoEl *gel = this->Reference();
    TPZManVector<int64_t,8> sideNodesId(nSideNodes);
    TPZManVector<int, 5> transformationIds(nSideConnects, -1);
    for (int ic=0; ic< nSideNodes; ic++) {
        const int localId = TSHAPE::SideNodeLocId(side,ic);
        sideNodesId[ic] = gel->Node(localId).Id();
    }
    for (auto iSide = nSideNodes; iSide < nContainedSides; iSide++) {
        MElementType sidetype = TSHAPE::Type(side);
        transformationIds[iSide-nSideNodes] = [&](){
            switch(sidetype){
                case EOned:
                    return pztopology::TPZLine::GetTransformId(iSide, sideNodesId);
                    break;
                case ETriangle:
                    return pztopology::TPZTriangle::GetTransformId(iSide, sideNodesId);
                    break;
                case EQuadrilateral:
                    return pztopology::TPZQuadrilateral::GetTransformId(iSide, sideNodesId);
                    break;
                default:
                    DebugStop();
                    return -1;
            }
        }();
//        const int localId = TSHAPE::ContainedSideLocId(side,iSide);
//        transformationIds[iSide-nSideNodes] = TSHAPE::GetTransformId(localId, elNodes);
    }


    //calculates the directions on the master element associated with the side and the indexes associating
    //the side with the respective h1 scalar function
    TPZFMatrix<REAL> sideMasterDirections(sideDim,sideDim * nContainedSides,0);
    TPZManVector<std::pair<int,int64_t>> indexVecShape(nSideShapes);
    {
        MElementType sidetype = TSHAPE::Type(side);
        TPZManVector<unsigned int,5> shapeCountVec(nSideConnects,-1);
        TPZManVector<int64_t,5> firstH1ShapeFunc(nSideConnects,-1);
        //calculates the first SCALAR shape function associated with each side of dim > 0
        TPZVec<int> sidesH1Ord(nSideConnects,-1);
        TPZVec<int> connectOrders(nSideConnects,-1);
        for(auto iCon = 0; iCon < nSideConnects; iCon++){
            const auto iSide = iCon + TSHAPE::NCornerNodes;
            const auto sideDim = TSHAPE::SideDimension(iSide);
            TPZStack<int> highDimSides;
            TSHAPE::HigherDimensionSides(iSide, highDimSides);
            const auto ord = this->EffectiveSideOrder(iSide);
            sidesH1Ord[iCon] = ord;
            connectOrders[iCon] = ord;

        }

        firstH1ShapeFunc[0] = nSideNodes;
        for (int iSide = nSideNodes + 1; iSide < nContainedSides; iSide++) {
            const int prevLocalId = TSHAPE::ContainedSideLocId(side,iSide-1);
            const int &lastFirstH1 = firstH1ShapeFunc[iSide - nSideNodes - 1];
            const int nShapeF = TSHAPE::NConnectShapeF(prevLocalId, sidesH1Ord[iSide-nSideNodes-1]);
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

        TPZFNMatrix<9,REAL> gradxSide(sideDim,sideDim,0);
        for(auto ix = 0; ix < sideDim; ix++) gradxSide(ix,ix) = 1;
        switch (sidetype) {
            case EOned://these wont be really used, just for signs purposes
                pztopology::TPZLine::ComputeHCurlDirections(gradxSide,sideMasterDirections,transformationIds);
                StaticIndexShapeToVec<pzshape::TPZShapeLinear>(indexVecShape,connectOrders,firstH1ShapeFunc,sidesH1Ord,shapeCountVec,transformationIds);
                break;
            case EQuadrilateral:
                pztopology::TPZQuadrilateral::ComputeHCurlDirections(gradxSide,sideMasterDirections,transformationIds);
                StaticIndexShapeToVec<pzshape::TPZShapeQuad>(indexVecShape,connectOrders,firstH1ShapeFunc,sidesH1Ord,shapeCountVec,transformationIds);
                break;
            case ETriangle:
                pztopology::TPZTriangle::ComputeHCurlDirections(gradxSide,sideMasterDirections,transformationIds);
                StaticIndexShapeToVec<pzshape::TPZShapeTriang>(indexVecShape,connectOrders,firstH1ShapeFunc,sidesH1Ord,shapeCountVec,transformationIds);
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
    TPZManVector<int,TSHAPE::NSides> ord(TSHAPE::NSides,connectOrder);
    TPZFNMatrix<50,REAL> phiSide(nSideShapes,1),dPhiSide(sideDim,nSideShapes);
    TSHAPE::SideShape(side, point, sideNodesId, ord, phiSide, dPhiSide);
    const int phiDim = [&](){
        switch(sideDim){
            case 1: return 1;
            case 2: return 3;
            default:
                DebugStop();
                return -1;
        }
    }();

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
const TPZVec<int>& transformationIds);\
template void TPZCompElHCurlFull<TSHAPE>::StaticIndexShapeToVec<pzshape::TPZShapeTriang>(TPZVec<std::pair<int,int64_t>> & indexVecShape,\
const TPZVec<int>& connectOrder,\
const TPZVec<int64_t>& firstH1ShapeFunc,const TPZVec<int>& sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,\
const TPZVec<int>& transformationIds);\
template void TPZCompElHCurlFull<TSHAPE>::StaticIndexShapeToVec<pzshape::TPZShapeQuad>(TPZVec<std::pair<int,int64_t>> & indexVecShape,\
const TPZVec<int>& connectOrder,\
const TPZVec<int64_t>& firstH1ShapeFunc,const TPZVec<int>& sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,\
const TPZVec<int>& transformationIds);


IMPLEMENTHCURLFULL(pzshape::TPZShapeLinear)
IMPLEMENTHCURLFULL(pzshape::TPZShapeTriang)
IMPLEMENTHCURLFULL(pzshape::TPZShapeQuad)
IMPLEMENTHCURLFULL(pzshape::TPZShapeCube)
IMPLEMENTHCURLFULL(pzshape::TPZShapeTetra)
IMPLEMENTHCURLFULL(pzshape::TPZShapePrism)

#undef IMPLEMENTHCURLFULL

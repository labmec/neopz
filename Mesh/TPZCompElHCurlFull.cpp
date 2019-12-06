/**
 * @file
 * @brief Contains the implementation of the TPZCompElHCurlmethods.
 */
#include <TPZCompElHCurlFull.h>

#include "pzcmesh.h"
#include "TPZTopologyUtils.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHCurl"));
static LoggerPtr loggercurl(Logger::getLogger("pz.mesh.tpzinterpolatedelement.divide"));
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
    const auto nFaces = TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NSides - 2 + TSHAPE::Dimension - nFaces - TSHAPE::NCornerNodes;
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
        else{//internal connect
            int count = 0;
            for(int iFace = 0; iFace < nFaces; iFace++){
                const int faceOrder = this->EffectiveSideOrder(TSHAPE::NCornerNodes+nEdges + iFace);
                switch(TSHAPE::Type(TSHAPE::NCornerNodes+nEdges + iFace)){
                    case ETriangle://triangular face
                        count += 0.5 * (faceOrder - 1) * ( faceOrder - 2);
                        break;
                    case EQuadrilateral://quadrilateral face
                        count +=  (faceOrder - 1) * ( faceOrder - 1);
                        break;
                    default:
                        PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
                        DebugStop();
                }
            }
            count += 3 * (TSHAPE::NConnectShapeF(side, sideOrder) );
            return count;
        }
    }();

#ifdef LOG4CXX
    if (logger->isDebugEnabled())
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
    const auto nFaces = TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NSides - 2 + TSHAPE::Dimension - nFaces - TSHAPE::NCornerNodes;
    const auto nNodes = TSHAPE::NCornerNodes;


#ifdef LOG4CXX
    std::stringstream sout;
    if (logger->isDebugEnabled()) {
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


    TPZVec<int64_t> nodes(nNodes, 0);
    for (auto iNode = 0; iNode < nNodes; iNode++) nodes[iNode] = this->Reference()->NodeIndex(iNode);
    //computing transformation id for sides
    TPZManVector<int, 21> transformationIds(TSHAPE::NSides - nNodes, -1);
    for (auto iSide = 0; iSide < nEdges + nFaces + TSHAPE::Dimension - 2; iSide++) {
        transformationIds[iSide] = TSHAPE::GetTransformId(nNodes + iSide, nodes);
    }

    TPZManVector<int64_t, TSHAPE::NSides - TSHAPE::NCornerNodes> firstH1ShapeFunc(TSHAPE::NSides - TSHAPE::NCornerNodes,
                                                                                  0);
    //calculates the first shape function associated with each side of dim > 0
    TPZVec<int> sidesH1Ord;
    this->CalculateSideShapeOrders(sidesH1Ord);
    firstH1ShapeFunc[0] = TSHAPE::NCornerNodes;
    for (int iSide = TSHAPE::NCornerNodes + 1; iSide < TSHAPE::NSides; iSide++) {
        const int iCon = iSide - TSHAPE::NCornerNodes;
        firstH1ShapeFunc[iCon] = firstH1ShapeFunc[iCon - 1] + TSHAPE::NConnectShapeF(iSide - 1, sidesH1Ord[iCon]);
    }

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        sout << "first H1 shape function:" << std::endl;
        for (auto &iShape : firstH1ShapeFunc) {
            sout << "\t" << iShape << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif


    indexVecShape.Resize(this->NShapeF());
    TPZVec<uint> shapeCountVec(TSHAPE::NSides - TSHAPE::NCornerNodes, 0);
    uint shapeCount = 0;


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
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        sout << "vec shape index (edge connects):" << std::endl;
        for (int iShape = 0; iShape < shapeCount; iShape++) {
            auto pair = indexVecShape[iShape];
            sout << "\tvec: " << pair.first << "\tshape: " << pair.second << std::endl;

        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

#ifdef PZDEBUG
    for (int iEdge = 0; iEdge < nEdges; iEdge++){
        if (shapeCountVec[iEdge] != this->NConnectShapeF(iEdge, connectOrder[iEdge])) {
            std::ostringstream soutAbort;
            soutAbort << __PRETTY_FUNCTION__ << std::endl;
            soutAbort << "\tError with the number of shape functions of edge " << iEdge + nNodes << std::endl;
            soutAbort << "\tCalculated " << shapeCountVec[iEdge] << " instead of "
                      << NConnectShapeF(iEdge, connectOrder[iEdge]) << std::endl;
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                LOGPZ_DEBUG(logger, sout.str() + soutAbort.str())
            }
#endif
            PZError << soutAbort.str() << std::endl;
        }
    }
#endif
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
        for(auto iFace = 1; iFace < nFaces; iFace++){
            TSHAPE::LowerDimensionSides(iFace - 1 + nEdges + nNodes, faceEdges[iFace-1], 1);
            const int nFaceEdges = faceEdges[iFace-1].size();
            firstVfeVec[iFace] = firstVfeVec[iFace - 1] + nFaceEdges;
            nFaceInternalFunctions[iFace] = TSHAPE::NConnectShapeF(iFace + nEdges + nNodes,connectOrder[nEdges + iFace]);
        }
        TSHAPE::LowerDimensionSides(nFaces - 1 + nEdges + nNodes, faceEdges[nFaces-1], 1);
        nFaceInternalFunctions[nFaces - 1] = TSHAPE::NConnectShapeF(nFaces - 1 + nEdges + nNodes,connectOrder[nEdges + nFaces - 1]);
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
            switch(TSHAPE::Type(iSide)){
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
        const int nFaceEdges = TSHAPE::NSideNodes(iSide);//this is not a mistake, since for faces nEdges = nNodes
        for(auto iEdge = 0; iEdge < nFaceEdges; iEdge++ ){
            const auto currentEdge = permutedSideSides[iEdge + nFaceEdges];
            const auto vecIndex = firstVfeVec[iFace] + iEdge;
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

    #ifdef PZDEBUG
    for (int iFace = 0; iFace < nFaces; iFace++){
        if (shapeCountVec[iFace + nEdges] != this->NConnectShapeF(iFace + nEdges, connectOrder[iFace + nEdges])) {
            std::ostringstream soutAbort;
            soutAbort << __PRETTY_FUNCTION__ << std::endl;
            soutAbort << "\tError with the number of shape functions of face " << iFace + nEdges + nNodes << std::endl;
            soutAbort << "\tCalculated " << shapeCountVec[iFace + nEdges] << " instead of "
                      << NConnectShapeF(iFace + nEdges, connectOrder[iFace + nEdges]) << std::endl;
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                LOGPZ_DEBUG(logger, sout.str() + soutAbort.str())
            }
#endif
            PZError << soutAbort.str() << std::endl;
        }
    }
    #endif

    if(TSHAPE::Dimension < 3) return;
    else {
        const int iCon = nEdges + nFaces;
        //first, the phi KF functions
        for(int iFace = 0; iFace < nFaces; iFace++){
            const auto nFaceInternalFuncs = nFaceInternalFunctions[iFace];
            const auto vecIndex = firstVfOrthVec + iFace;
            for(auto iFunc = 0; iFunc < nFaceInternalFuncs; iFunc++ ){
                const auto shapeIndex = firstH1ShapeFunc[nEdges + iFace] + iFunc;
                indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
                shapeCount++;
                shapeCountVec[iCon]++;
            }
        }
        //now the phi Ki funcs
        const int firstInternalVec = firstVfOrthVec + nFaces;
        const auto nInternalFuncs = 3 * TSHAPE::NConnectShapeF(TSHAPE::NSides - 1, connectOrder[iCon]);
        for(auto iFunc = 0; iFunc < nInternalFuncs; iFunc++ ){
            const auto vecIndex = firstInternalVec + iFunc % 3;//it should alternate between them
            const auto shapeIndex = firstH1ShapeFunc[iCon] + iFunc / 3;
            indexVecShape[shapeCount] = std::make_pair(vecIndex,shapeIndex);
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
        ord[iCon] = this->EffectiveSideOrder(iCon + TSHAPE::NCornerNodes);
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
template class TPZCompElHCurlFull<TSHAPE>;

IMPLEMENTHCURLFULL(pzshape::TPZShapeLinear)
IMPLEMENTHCURLFULL(pzshape::TPZShapeTriang)
IMPLEMENTHCURLFULL(pzshape::TPZShapeQuad)
IMPLEMENTHCURLFULL(pzshape::TPZShapeCube)
IMPLEMENTHCURLFULL(pzshape::TPZShapeTetra)
IMPLEMENTHCURLFULL(pzshape::TPZShapePrism)

#undef IMPLEMENTHCURLFULL
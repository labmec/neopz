/**
 * @file
 * @brief Contains the implementation of the TPZCompElHCurlmethods.
 */
#include <TPZCompElHCurlFull.h>

#include "pzcmesh.h"

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
void TPZCompElHCurlFull<TSHAPE>::IndexShapeToVec(TPZVec<std::pair<int,int64_t> > & indexVecShape, const TPZVec<int>& connectOrder) const{
    DebugStop();
}

template<class TSHAPE>
void TPZCompElHCurlFull<TSHAPE>::CalculateSideShapeOrders(TPZVec<int> &ord) const{
    DebugStop();
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
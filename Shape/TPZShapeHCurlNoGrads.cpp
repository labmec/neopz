#include "TPZShapeHCurlNoGrads.h"
#include "TPZShapeHCurl.h"
#include "pzeltype.h"

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"

template<class TSHAPE>
void TPZShapeHCurlNoGrads<TSHAPE>::Initialize(TPZVec<int64_t> &ids,
                                       TPZVec<int> &connectorders,
                                       TPZShapeData &data)
{
  TPZShapeHCurl<TSHAPE>::Initialize(ids,connectorders,data);
  constexpr int ncon = TSHAPE::NSides-TSHAPE::NCornerNodes;
  data.fHDivNumConnectShape.Resize(ncon);
  //we need to update the number of filtered hcurl functions
  for (int i = 0; i < ncon; i++){
    data.fHDivNumConnectShape[i] = ComputeNConnectShapeF(i,connectorders[i]);
  }
}



template<class TSHAPE>
int TPZShapeHCurlNoGrads<TSHAPE>::NHCurlShapeF(const TPZShapeData &data)
{
  return TPZShapeHCurl<TSHAPE>::NHCurlShapeF(data);
}

    

template<class TSHAPE>
void TPZShapeHCurlNoGrads<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &curlphi)
{

    constexpr int ncorner = TSHAPE::NCornerNodes;
    constexpr int nsides = TSHAPE::NSides;
    constexpr int ncon = nsides - ncorner;
    constexpr int dim = TSHAPE::Dimension;
    constexpr int curldim = [dim](){
        if constexpr (dim == 1) return 1;
        else{
            return 2*dim - 3;//1 for 2D 3 for 3D
        }
    }();
    const int nedges = TSHAPE::NumSides(1);
    
    //calculates # of unfiltered hcurl functions
    const auto &connectorders = data.fHDivConnectOrders;
    //first_hcurl_side[i] is the index of the first shape function associated with side i
    TPZManVector<int,ncon> first_hcurl_side(ncon,0);
    //total number of unfiltered funcs
    int n_unfilt = 0;

    {
      n_unfilt += TPZShapeHCurl<TSHAPE>::ComputeNConnectShapeF(0,connectorders[0]);
      for (int i = 1; i < ncon; i++){
        first_hcurl_side[i] = n_unfilt;
        n_unfilt += TPZShapeHCurl<TSHAPE>::ComputeNConnectShapeF(i,connectorders[i]);
      }
    }
    
    //computes  unfiltered hcurl funcs
    TPZFMatrix<REAL> phi_unfilt(n_unfilt,dim), curlphi_unfilt(curldim,n_unfilt);
    TPZShapeHCurl<TSHAPE>::Shape(pt, data, phi_unfilt, curlphi_unfilt);

    //now we filter the functions
    int fcount = 0;
    //edges: we sum the lowest order functions on the edge
    for(auto ie = 0; ie < nedges; ie++){
      //fss = first side shape
      const auto fss = first_hcurl_side[ie];
      for(auto x = 0; x < dim; x++){
        phi(fcount,x) += phi_unfilt(fss,x) + phi_unfilt(fss+1,x);
      }
      for(auto x = 0; x < curldim; x++){
        curlphi(x,fcount) += curlphi_unfilt(x,fss) + curlphi_unfilt(x,fss+1);
      }
      fcount++;
    }
    if constexpr (dim < 2) return;

    TPZVec<int> filtVecShape;
    HighOrderFunctionsFilter(first_hcurl_side, connectorders,filtVecShape);
    
    const auto newfuncs = filtVecShape.size();
    for(int ifunc = 0; ifunc < newfuncs; ifunc++){
      const auto fi = filtVecShape[ifunc];
      for(auto x = 0; x < curldim; x++){
        curlphi(x,fcount) += curlphi_unfilt(x,fi);
      }
      for(auto x = 0; x < dim; x++){
        phi(fcount,x) += phi(fi,x);
      }
      fcount++;
    } 
    
    
}


template<class TSHAPE>
int TPZShapeHCurlNoGrads<TSHAPE>::ComputeNConnectShapeF(const int icon, const int order)

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
          return 1;
        }
        else if(side < nNodes + nEdges + nFaces){//face connect
            switch(TSHAPE::Type(side)){
            case ETriangle://triangular face
              return (order - 1) * (order+2)/2;
            default:
              PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
              DebugStop();
              return 0;
            }
        }
        else{//internal connect (3D element only)
            if constexpr (TSHAPE::Type() == ETetraedro){
                return (order-1)*(order-2)*(2*order+3)/6;
            }
            else{
                PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
                DebugStop();
                return 0;
            }
        }
    }();
    return nShapeF;
 }

template<class TSHAPE>
int TPZShapeHCurlNoGrads<TSHAPE>::MaxOrder(const int ordh1){
  return TPZShapeHCurl<TSHAPE>::MaxOrder(ordh1);
}

template<class TSHAPE>
void TPZShapeHCurlNoGrads<TSHAPE>::HighOrderFunctionsFilter(
    const TPZVec<int> &firstHCurlFunc,
    const TPZVec<int> &conOrders,
    TPZVec<int> &filteredFuncs)
{
    constexpr auto dim = TSHAPE::Dimension;
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto nConnects = TSHAPE::NSides - nNodes;
    const auto nFaces = TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);

    filteredFuncs.Resize(0);
    /*
        face connects: 
    */
    for(auto iface = 0; iface < nFaces; iface++){
        const auto icon = nEdges + iface;
        const auto side = icon + nNodes;
        const auto firstSideShape = firstHCurlFunc[icon];
        const auto order = conOrders[icon];

        switch(TSHAPE::Type(side)){
        case ETriangle:{
        //there are no face functions for k < 2
        if(order < 2) break;
        /**
             we remove one internal function for each h1 face function of order k+1
            since there are (k-1)(k-2)/2 functions per face in a face with order k,
            we remove k(k-1)/2.
            so:
            (k-1)*(k+1)-k*(k-1)/2 = (k-1)(k+2)/2.
        */
        const auto nfacefuncs =  (order - 1) * (order+2) / 2;

        auto fcount = filteredFuncs.size();
        filteredFuncs.Resize(fcount+nfacefuncs);

        /**
             we will iterate over the phi_fe hcurl functions. there are 3(k-1) vfe funcs.
            that means that, at each k, there are 3 new vfe functions. we can remove 
            one of them for each polynomial order.
        */
        auto firstVfeK = firstSideShape;
        for(auto ik = 2; ik <= order; ik++){
            for(auto ifunc = 0; ifunc < 2; ifunc++){
            filteredFuncs[fcount] = firstVfeK+ifunc;
            fcount++;
            }
            //we skip to the higher order ones
            firstVfeK += 3;
        }
        /**
             there are already 2(k-1) filtered functions. we discarded (k-1).
            so that means we need more (k-1)(k-2)/2 functions, because
            2(k-1) + (k-1)(k-2)/2 = (k-1)(k+2)/2. 
            That means we can take exactly half of the phi_fi functions. 
            for each k, there are 2(k-2) phi_fi func. so...
        **/
        auto firstVfiK = firstSideShape + 3*(order-1);
        for(auto ik = 3; ik < order; ik++){
            const auto nPhiFi = 2*(ik-1);
            for(auto ifunc = 0; ifunc < nPhiFi; ifunc+=2){
            //we always skip one of them
            filteredFuncs[fcount] = firstVfiK+ifunc;
            fcount++;
            }
            //we skip to the higher order ones
            firstVfeK += nPhiFi;
        }
        break;
        }
        default:
        PZError<<__PRETTY_FUNCTION__<<" error. Not yet implemented"<<std::endl;
        DebugStop();
        }
    }
    
    if constexpr (dim < 3) return;

    if constexpr (TSHAPE::Type() == ETetraedro){
        const auto icon = nEdges + nFaces;
        const auto order = conOrders[icon];
        const auto firstSideShape = firstHCurlFunc[icon];
        /**
             we remove one internal function for each h1 internal function of order k+1
            since there are (k-1)(k-2)(k-3)/6 functions in a h1 element with order k,
            we remove k(k-1)(k-2)/6.
            so:
            (k-1)(k-2)(k+1)/2 - k(k-1)(k-2)/6 = (k-1)(k-2)(2k+3)/6.

            since we will remove k(k-1)(k-2)/6, for each     we remove (k-1)(k-2)/2 funcs.

            we have two kinds of internal functions. phi_kf and phi_ki.
            func        k                   k-1                 new funcs
            phi_kf      2(k-1)(k-2)         2(k-2)(k-3)         4(k-2)
            phi_ki      (k-1)(k-2)(k-3)/2   (k-2)(k-3)(k-4)/2   3(k-2)(k-3)/2
            
            that means that if we remove, for each k, (k-2) phi_kf 
            (for instance, all phi_kf associated with a given face),
            we need to remove (k-1)(k-2)/2 - (k-2) = (k-2)(k-3)/2
            which is exactly one third of the phi_ki.
            
        */

        const auto nintfuncs =  (order - 1) * (order-2) * (2*order+ 3) / 6;

        auto fcount = filteredFuncs.size();
        filteredFuncs.Resize(fcount+nintfuncs);

        /**
         we will iterate over the phi_kf hcurl functions.
        */
        auto firstVkf = firstSideShape;
        for(auto ik = 3; ik <= order; ik++){
        //we chose to remove all the functions for a given face
        const auto newvkf = 4*(ik-2);
        for(auto ifunc = ik-2; ifunc < newvkf; ifunc++){
            filteredFuncs[fcount] = firstVkf+ifunc;
            fcount++;
        }
        //we skip to the higher order ones
        firstVkf += newvkf;
        }

        /**
         we now iterate over the phi_ki hcurl functions
        */
        const auto nvkf = 2*(order-1)*(order-2);
        auto firstVki = firstSideShape + nvkf;
        for(auto ik = 4; ik <= order; ik++){
        //we chose to remove all the functions for a given direction
        const auto newvki = 3*(ik-2)*(ik-3)/2;
        for(auto ifunc = 0; ifunc < newvki; ifunc++){
            if(ifunc%3 == 0) continue;
            filteredFuncs[fcount] = firstVki+ifunc;
            fcount++;
        }
        //we skip to the higher order ones
        firstVkf += newvki;
        }
    }
}


template
struct TPZShapeHCurlNoGrads<pzshape::TPZShapeLinear>;

template
struct TPZShapeHCurlNoGrads<pzshape::TPZShapeTriang>;

template
struct TPZShapeHCurlNoGrads<pzshape::TPZShapeQuad>;

template
struct TPZShapeHCurlNoGrads<pzshape::TPZShapeTetra>;

template
struct TPZShapeHCurlNoGrads<pzshape::TPZShapeCube>;

template
struct TPZShapeHCurlNoGrads<pzshape::TPZShapePrism>;
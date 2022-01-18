#include "TPZCompElHCurlNoGrads.h"

#include <TPZMaterial.h>
#include <pzconnect.h>
#include <pzcmesh.h>
#include "TPZShapeHCurl.h"
#include "TPZShapeHDivKernel.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif


template<class TSHAPE>
TPZCompElHCurlNoGrads<TSHAPE>::TPZCompElHCurlNoGrads() : TPZCompElHCurl<TSHAPE>()
{
    this->AdjustConnects();
}

template<class TSHAPE>
TPZCompElHCurlNoGrads<TSHAPE>::TPZCompElHCurlNoGrads(
  TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
  TPZCompElHCurl<TSHAPE>(mesh,gel,index)
{
    this->AdjustConnects();
}

template<class TSHAPE>
void TPZCompElHCurlNoGrads<TSHAPE>::AdjustConnects()
{
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto ncon = TSHAPE::NSides - nNodes;
    for(int icon = 0; icon < ncon; icon++){
        const int connect = this->MidSideConnectLocId(icon+nNodes);
        TPZConnect &c = this->Connect(connect);
        const int nshape =this->NConnectShapeF(connect,c.Order());
        c.SetNShape(nshape);
        const auto seqnum = c.SequenceNumber();
        const int nStateVars = [&](){
            TPZMaterial * mat =this-> Material();
            if(mat) return mat->NStateVariables();
            else {
                return 1;
            }
        }();
        this-> Mesh()->Block().Set(seqnum,nshape*nStateVars);
    }
}
  
template<class TSHAPE>
void TPZCompElHCurlNoGrads<TSHAPE>::InitMaterialData(TPZMaterialData &data){

    //Init the material data of Hcurl
    TPZCompElHCurl<TSHAPE>::InitMaterialData(data);

    /*
        we first compute ALL the hcurl traditional functions(scalar+vector)
        then, at each integration point we will filter them.
        
        some of this code is copied from the TPZCompElHCurlFull since we want to avoid
        calls for the virtual methods.
    */
    
    //computes the index that will associate each scalar function to a constant vector field
    constexpr auto nConnects = TSHAPE::NSides - TSHAPE::NCornerNodes;
    TPZManVector<int,nConnects> connectOrders(nConnects,-1);
    int unfiltnshape = 0;
    for(auto i = 0; i < nConnects; i++){
        const auto conorder = this->EffectiveSideOrder(i + TSHAPE::NCornerNodes);
        connectOrders[i] = conorder;
        unfiltnshape += TPZCompElHCurl<TSHAPE>::NConnectShapeF(i, conorder);
    }

    const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    constexpr auto nNodes = TSHAPE::NCornerNodes;

    TPZManVector<int64_t,nNodes> nodes(nNodes, 0);

    for (auto iNode = 0; iNode < nNodes; iNode++){
        nodes[iNode] = this->Reference()->NodeIndex(iNode);
    }
   
    TPZManVector<int64_t, TSHAPE::NSides - nNodes>
        firstH1ShapeFunc(TSHAPE::NSides - nNodes,0);
    
    //calculates the first shape function associated with each side of dim > 0
    TPZManVector<int,TSHAPE::NSides-nNodes> sidesH1Ord(TSHAPE::NSides - nNodes,-1);
    // TPZShapeHDivKernel<TSHAPE>::CalcH1ShapeOrders(connectOrders,sidesH1Ord);
    TPZShapeHCurl<TSHAPE>::CalcH1ShapeOrders(connectOrders,sidesH1Ord);
    firstH1ShapeFunc[0] = nNodes;
    
    for (int iSide = nNodes + 1; iSide < TSHAPE::NSides; iSide++) {
        const int iCon = iSide - nNodes;
        firstH1ShapeFunc[iCon] =
        firstH1ShapeFunc[iCon - 1] +
        TSHAPE::NConnectShapeF(iSide - 1, sidesH1Ord[iCon-1]);

    }
    
    auto &indexVecShape = data.fVecShapeIndex;
    indexVecShape.Resize(unfiltnshape);

    TPZVec<unsigned int> shapeCountVec(TSHAPE::NSides - nNodes, 0);
    TPZShapeHCurl<TSHAPE>::StaticIndexShapeToVec(connectOrders,
                                                 firstH1ShapeFunc,sidesH1Ord, nodes, shapeCountVec,indexVecShape );

    //setting the type of shape functions as vector shape functions
    data.fShapeType = TPZMaterialData::EVecShape;
}

template<class TSHAPE>
template<class TVar>
void TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredDataT(
  TPZMaterialDataT<TVar> &data, TPZVec<REAL> &qsi)
{
    // We can't use ComputeRequired data from TPZCompElHCurl, because we are dealing with a number of shape functions
    // uncompatible with the number of HCurl shape functions
  
    /******************************************************************************************************************
     * at this point, we already have the basis functions on the deformed element, since we have data.phi,
     * data.fVecShapeIndex and data.fDeformedDirections. Now it is time to compute the curl, which will be stored in
     * data.curlphi.
     *******************************************************************************************************************/
    
    //Compute the element geometric data
    TPZGeoEl * ref = this->Reference();
    if (!ref){
        PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
        return;
    }
    ref->Jacobian(qsi, data.jacobian, data.axes, data.detjac , data.jacinv);

    ref->X(qsi, data.x);
    data.xParametric = qsi;
    

    TPZFMatrix<REAL> phiHCurl;
    TPZFMatrix<REAL> unfiltCurl;
    constexpr auto dim{TSHAPE::Dimension};
    data.curlphi.Redim(2*dim - 3 > 0 ? 2*dim - 3 : 1, this->NShapeF());
    //number of FILTERED hcurl functions
    const auto nshape = this->NShapeF();
    const auto &vecShapeIndex = data.fVecShapeIndex;
    const auto &deformedDirs = data.fDeformedDirections;
    const auto &phiH1 = data.phi;
    phiHCurl.Resize(dim,vecShapeIndex.size());
    unfiltCurl.Resize(dim,vecShapeIndex.size());

    TPZShapeHCurl<TSHAPE>::Shape(qsi,data,phiHCurl,unfiltCurl);
    
    TPZCompElHCurl<TSHAPE>::TransformShape(phiHCurl, data.detjac, data.jacinv, data.axes, data.phi);
    TPZCompElHCurl<TSHAPE>::TransformCurl(unfiltCurl, data.detjac, data.jacobian, data.curlphi);

    phiHCurl = data.phi;
    unfiltCurl = data.curlphi;

// #ifdef PZ_LOG
//     if (logger.isDebugEnabled())
//     {
//         std::stringstream sout;
//         //	this->Print(sout);
//         sout << "\nUnfiltCurl = " << unfiltCurl << std::endl;
//         LOGPZ_DEBUG(logger,sout.str())
        
//     }
// #endif

    // constexpr auto dim = TSHAPE::Dimension;
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto nConnects = TSHAPE::NSides - nNodes;
    const auto nFaces = TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);

    //number of FILTERED hcurl functions
    // const int nshape = this->NShapeF();
    //unfiltered shape functions count for each connect
    TPZManVector<int,nConnects> firstHCurlFunc(nConnects,0);
    TPZManVector<int,nConnects> conOrders(nConnects,0);
    for(auto icon = 1; icon < nConnects; icon++){
        const auto conorder = conOrders[icon] = this->ConnectOrder(icon);
        firstHCurlFunc[icon] = firstHCurlFunc[icon-1] +
        TPZCompElHCurl<TSHAPE>::NConnectShapeF(icon,conorder);
    }
    
    constexpr auto curlDim= 2*dim - 3 > 0 ? 2*dim -3 : 1;

    auto &curlPhi = data.curlphi;
    curlPhi.Redim(curlDim,nshape);

    int fcount = 0;
    //edges
    for(auto ie = 0; ie < nEdges; ie++){
        //fss = first side shape
        const auto fss = firstHCurlFunc[ie];
        for(auto x = 0; x < curlDim; x++){
            curlPhi(x,fcount) += unfiltCurl(x,fss) + unfiltCurl(x,fss+1);
        }
        fcount++;
    }
    if constexpr (dim < 2) return;

    TPZVec<int> filtVecShape;
    this->HighOrderFunctionsFilter(firstHCurlFunc, conOrders,
                                    filtVecShape);
    
    const auto newfuncs = filtVecShape.size();
    data.phi.Resize(nshape,dim);
    data.phi.Zero();
    for(int ifunc = 0; ifunc < newfuncs; ifunc++){
        const auto fi = filtVecShape[ifunc];
        for(auto x = 0; x < curlDim; x++){
            curlPhi(x,fcount) += unfiltCurl(x,fi);
        }
        for(auto x = 0; x < dim; x++){
            data.phi(fcount,x) += phiHCurl(fi,x);
        }
        fcount++;
    } 

    // data.phi = phiHCurl;
    if (data.fNeedsSol) {
        this->ReallyComputeSolution(data);
    }
}

template<class TSHAPE>
template<class TVar>
void TPZCompElHCurlNoGrads<TSHAPE>::ReallyComputeSolutionT(TPZMaterialDataT<TVar> &data)
{
    this->ComputeSolutionHCurlT(data.phi, data.curlphi,
                         data.sol, data.curlsol);
}


template<class TSHAPE>
int TPZCompElHCurlNoGrads<TSHAPE>::NConnectShapeF(int icon, int order) const
{
    const int side = icon + TSHAPE::NCornerNodes;
    #ifdef PZDEBUG
    if (side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) {
        DebugStop();
    }
    #endif
    if(order == 0) {
        PZError<<__PRETTY_FUNCTION__
            <<"\nERROR: polynomial order not compatible.\nAborting..."
            <<std::endl;
        DebugStop();
        return 0;
    }
    const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    const int nShapeF = [&](){
        if (side < TSHAPE::NCornerNodes + nEdges) {//edge connect
        return 1;
        }
        else if(side < TSHAPE::NCornerNodes + nEdges + nFaces){//face connect
        switch(TSHAPE::Type(side)){
        case ETriangle://triangular face
            /**
             we remove one internal function for each h1 face function of order k+1
            since there are (k-1)(k-2)/2 functions per face in a face with order k,
            we remove k(k-1)/2.
            so:
            (k-1)*(k+1)-k*(k-1)/2
            */
            return (order - 1) * (order+2) / 2;
        case EQuadrilateral://quadrilateral face
            //Following the same logic:
            /**
             we remove one internal function for each h1 face function of order k+1
            since there are (k-1)^2 functions per face in a face with order k,
            we remove k^2.
            so:
            2k(k+1) - k^2 = k(k+2)
            
            */

            return order * (order + 2);
        default:
            PZError<<__PRETTY_FUNCTION__<<" error. Not yet implemented"<<std::endl;
            DebugStop();
            return 0;
        }
        }
        else{//internal connect (3D element only)
        if constexpr (TSHAPE::Type() == ETetraedro){
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
            return (order-1)*(order-2)*(2*order+3)/6;
        } else 
        if constexpr (TSHAPE::Type() == ECube){
        /**
             we remove one internal function for each h1 internal function of order k+1
            since there are (k-1)^3 functions in a h1 element with order k,
            we remove k^3.
            so:
            3k^2(k+1) - k^3 = k^2 (3 + 2 k).

            we have two kinds of internal functions. phi_kf and phi_ki.
            func        k                   k-1                 new funcs
            phi_kf      6k^2                6(k-1)^2            k(8k-k^2-1)
            phi_ki      3k^2*(k-1)          3(k-1)(k-1)(k-2)/2  3(2-5k+2k^2+k^3)/2
           
        */
            switch (order)
            {
            case 1:
                return 5;
            case 2:
                return 19;
            case 3: 
                return 37;
            case 4:
                return 56;
            default:
                break;
            }
            return order*order*(3+2*order);
        }
        return 0;
        // int count = 0;
        // //first we count the face-based interior functions \phi^{K,F}
        // for(int iFace = 0; iFace < nFaces; iFace++){
        //   switch(TSHAPE::Type(TSHAPE::NCornerNodes+nEdges + iFace)){
        //   case ETriangle://triangular face
        //     count +=  (order - 1) * ( order - 2)/2;
        //     break;
        //   case EQuadrilateral://quadrilateral face
        //     //we need the functions of k+1
        //     count +=  order * order;
        //     break;
        //   default:
        //     PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
        //     DebugStop();
        //   }
        // }

        // const int nVkf = count;
        // //now we count the interior bubble functions
                
        // //number of H1 funcs
        // const auto orderH1 = TSHAPE::Type() == ECube ?
        //   order +1 : order;
                
        // const auto nFuncsH1 = TSHAPE::NConnectShapeF(side, orderH1);
        // if constexpr (TSHAPE::Type() == ECube) {
        //   //shapeorders[k] contains the polynomial orders
        //   //in xi, eta and zeta for each of the internal funcs
        //   TPZGenMatrix<int> shapeorders(nFuncsH1,3);
        //   //not really used since we are interested in internal funcs
        //   TPZManVector<int64_t,0> idvec(0);
        //   TSHAPE::SideShapeOrder(side,idvec,orderH1,shapeorders);              
        //   /*for compatibility, we need the spaces
        //     Q^{k,k+1,k+1} \times Q^{k+1,k,k+1} \times Q^{k+1,k+1,k}
        //     I am supposing that #funcs in xi-dir is equal to 
        //     #funcs in eta-dir and in zeta-dir. So I will just
        //     see how many we will have in xi-dir after skipping
        //     internals that were included in the k-1 set of funcs*/
        //   int countInt = 0;
        //   for(int ifunc = 0; ifunc < nFuncsH1; ifunc++){
        //     if(shapeorders(ifunc,0) <= order) {countInt++;}
        //   }
        //   count += 3* countInt;

        //   const auto nVki = count - nVkf;
        //   // std::stringstream sout;
        //   // sout << __PRETTY_FUNCTION__<<'\n'
        //   //      <<"\tside "<<side<<"\tcon "<<icon<<"\torder "<<order<<'\n'
        //   //      <<"\tnfuncs "<<count<<"\tnvkf "<<nVkf<<"\tnvki "<<nVki<<'\n';
        //   // std::cout<<sout.str()<<std::flush;
        // }
        // else{
        //   count += 3 * nFuncsH1;
        // }
        // return count;
        }
    }();
    return nShapeF;
}

template<class TSHAPE>
void TPZCompElHCurlNoGrads<TSHAPE>::ComputeShape(TPZMaterialData &data,
                                                 TPZFMatrix<REAL> &phiHCurl)
{
    constexpr auto dim = TSHAPE::Dimension;
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto nConnects = TSHAPE::NSides - nNodes;
    const auto nFaces = TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);
    
    //unfiltered shape functions count for each connect
    TPZManVector<int,nConnects> firstHCurlFunc(nConnects,0);
    TPZManVector<int,nConnects> conOrders(nConnects,0);
    for(auto icon = 1; icon < nConnects; icon++){
        const auto conorder = conOrders[icon] = this->ConnectOrder(icon);
        firstHCurlFunc[icon] = firstHCurlFunc[icon-1] +
        TPZCompElHCurl<TSHAPE>::NConnectShapeF(icon,conorder);
    }

    //number of FILTERED hcurl functions
    const auto nshape = this->NShapeF();
    
    const auto &vecShapeIndex = data.fVecShapeIndex;
    const auto &deformedDirs = data.fDeformedDirections;
    const auto &phiH1 = data.phi;

            
    phiHCurl.Resize(nshape, dim);
    int fcount = 0;
            
    /*
        edge connects: we combine the lowest order edge functions in order to
        create a function with constant tangential trace
    */

    for(auto ie = 0; ie < nEdges; ie++){
        const auto firstSideShape = firstHCurlFunc[ie];
        const auto vIndex1 = vecShapeIndex[firstSideShape].first;
        const auto sIndex1 = vecShapeIndex[firstSideShape].second;
        const auto vIndex2 = vecShapeIndex[firstSideShape+1].first;
        const auto sIndex2 = vecShapeIndex[firstSideShape+1].second;
        for(auto x = 0; x < dim; x++){
        phiHCurl(fcount,x) = phiH1.GetVal(sIndex1,0) *
            deformedDirs.GetVal(x,vIndex1) +
            phiH1.GetVal(sIndex2,0) *
            deformedDirs.GetVal(x,vIndex2);
        }
        fcount++;
    }
    if constexpr (dim < 2) return;

    // TPZCompElHCurl<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &curlphi);

    TPZVec<int> filtVecShape;
    this->HighOrderFunctionsFilter(firstHCurlFunc, conOrders,
                                    filtVecShape);
    
    const auto newfuncs = filtVecShape.size();

    for(int ifunc = 0; ifunc < newfuncs; ifunc++){
        const auto vs = filtVecShape[ifunc];
        const auto vi = vecShapeIndex[vs].first;
        const auto si = vecShapeIndex[vs].second;
        for(auto x = 0; x < dim; x++){
        phiHCurl(fcount,x) = phiH1.GetVal(si,0) *
            deformedDirs.GetVal(x,vi);
        }
        fcount++;
    } 
}

template<class TSHAPE>
template<int DIM>
void TPZCompElHCurlNoGrads<TSHAPE>::ComputeCurl(TPZMaterialData &data)
{

    
    TPZFMatrix<REAL> unfiltCurl;
    //compute curl of ALL HCurl functions
    // TPZCompElHCurl<TSHAPE>::ComputeCurl(
    //     data.fVecShapeIndex,data.fDPhi,data.fMasterDirections,
    //     data.jacobian,data.detjac,data.axes,unfiltCurl);//Jeferson
    

    constexpr auto dim = TSHAPE::Dimension;
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto nConnects = TSHAPE::NSides - nNodes;
    const auto nFaces = TSHAPE::NumSides(2);
    const auto nEdges = TSHAPE::NumSides(1);

    //number of FILTERED hcurl functions
    const int nshape = this->NShapeF();
    //unfiltered shape functions count for each connect
    TPZManVector<int,nConnects> firstHCurlFunc(nConnects,0);
    TPZManVector<int,nConnects> conOrders(nConnects,0);
    for(auto icon = 1; icon < nConnects; icon++){
        const auto conorder = conOrders[icon] = this->ConnectOrder(icon);
        firstHCurlFunc[icon] = firstHCurlFunc[icon-1] +
        TPZCompElHCurl<TSHAPE>::NConnectShapeF(icon,conorder);
    }
    
    constexpr auto curlDim= 2*dim - 3 > 0 ? 2*dim -3 : 1;

    auto &curlPhi = data.curlphi;
    curlPhi.Redim(curlDim,nshape);

    int fcount = 0;
    //edges
    for(auto ie = 0; ie < nEdges; ie++){
        //fss = first side shape
        const auto fss = firstHCurlFunc[ie];
        for(auto x = 0; x < curlDim; x++){
        curlPhi(x,fcount) += unfiltCurl(x,fss) + unfiltCurl(x,fss+1);
        }
        fcount++;
    }
    if constexpr (dim < 2) return;

    TPZVec<int> filtVecShape;
    this->HighOrderFunctionsFilter(firstHCurlFunc, conOrders,
                                    filtVecShape);
    
    const auto newfuncs = filtVecShape.size();

    for(int ifunc = 0; ifunc < newfuncs; ifunc++){
        const auto fi = filtVecShape[ifunc];
        for(auto x = 0; x < curlDim; x++){
        curlPhi(x,fcount) += unfiltCurl(x,fi);
        }
        fcount++;
    } 
}

template<class TSHAPE>
void TPZCompElHCurlNoGrads<TSHAPE>::HighOrderFunctionsFilter(
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

#include <pzshapetriang.h>
#include <pzshapetetra.h>
#include <pzshapecube.h>
#include <pzshapequad.h>
#include <pzshapeprism.h>


#define IMPLEMENTHCURLNOGRADS(TSHAPE)                           \
                                                                \
    template class                                                \
    TPZRestoreClass< TPZCompElHCurlNoGrads<TSHAPE> >;             \
    template class TPZCompElHCurlNoGrads<TSHAPE>;                 \
                                                                    \
    template void                                                 \
    TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredDataT<STATE>(   \
        TPZMaterialDataT<STATE> &data,TPZVec<REAL> &qsi);           \
    template void                                                 \
    TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredDataT<CSTATE>(  \
        TPZMaterialDataT<CSTATE> &data, TPZVec<REAL> &qsi);

IMPLEMENTHCURLNOGRADS(pzshape::TPZShapeTriang)
IMPLEMENTHCURLNOGRADS(pzshape::TPZShapeTetra)
IMPLEMENTHCURLNOGRADS(pzshape::TPZShapeQuad)
IMPLEMENTHCURLNOGRADS(pzshape::TPZShapeCube)
IMPLEMENTHCURLNOGRADS(pzshape::TPZShapePrism)

#undef IMPLEMENTHCURLNOGRADS
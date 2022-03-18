/**
 * @file
 * @brief Contains the implementation of the TPZCompElHCurlmethods.
 */
#include "TPZCompElHCurl.h"

#include "TPZShapeHCurl.h"
#include "TPZShapeHCurlNoGrads.h"
#include "TPZMaterial.h"
#include "pzcmesh.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzstepsolver.h"
#include "pzcheckrestraint.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHCurl");
static TPZLogger loggercurl("pz.mesh.tpzinterpolatedelement.divide");
#else
static int logger;
#endif

/*********************************************************************************************************
                                       TPZCompElHCurl methods
 *********************************************************************************************************/

template<class TSHAPE>
TPZCompElHCurl<TSHAPE>::TPZCompElHCurl(TPZCompMesh &mesh, TPZGeoEl *gel,
                                       const HCurlFamily hcurlfam) :
TPZRegisterClassId(&TPZCompElHCurl::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel,1), fhcurlfam(hcurlfam)
{
    gel->SetReference(this);
    this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
    this->CreateHCurlConnects(mesh);
}

template<class TSHAPE>
TPZCompElHCurl<TSHAPE>::TPZCompElHCurl(TPZCompMesh &mesh, const TPZCompElHCurl<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHCurl::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy), fhcurlfam(copy.fhcurlfam)
{

}

template<class TSHAPE>
TPZCompElHCurl<TSHAPE>::TPZCompElHCurl(TPZCompMesh &mesh,
									 const TPZCompElHCurl<TSHAPE> &copy,
									 std::map<int64_t,int64_t> & gl2lcConMap,
									 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHCurl::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fhcurlfam(copy.fhcurlfam)
{
	int i;
	for(i=0;i<NConnects();i++)
	{
		int lcIdx = -1;
		int glIdx = copy.fConnectIndexes[i];
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			this-> fConnectIndexes[i] = -1;
			return;
		}
		this-> fConnectIndexes[i] = lcIdx;
	}
}

template<class TSHAPE>
TPZCompElHCurl<TSHAPE>::TPZCompElHCurl() :
TPZRegisterClassId(&TPZCompElHCurl::ClassId),
TPZIntelGen<TSHAPE>()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides-TSHAPE::NCornerNodes;i++) {
		this-> fConnectIndexes[i] = -1;
	}

}

template<class TSHAPE>
TPZCompElHCurl<TSHAPE>::~TPZCompElHCurl(){
    TPZGeoEl *gel = this->Reference();
    if (gel && gel->Reference() != this) {
        return;
    }
    if (gel) {
        TPZCompEl *cel = gel->Reference();
        if (cel == this) {
            this->RemoveSideRestraintsII(TPZInterpolatedElement::MInsertMode::EDelete);
        }
        this->Reference()->ResetReference();
    }
    TPZStack<int64_t > connectlist;
    this->BuildConnectList(connectlist);
    int64_t nconnects = connectlist.size();
    for (int ic = 0; ic < nconnects; ic++) {
        if (connectlist[ic] != -1){
            this->fMesh->ConnectVec()[connectlist[ic]].DecrementElConnected();
        }
    }
    if (gel){
        gel->ResetReference();
    }
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHCurl") ^ TPZIntelGen<TSHAPE>::ClassId() << 1;
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsHCurl();
}


template<class TSHAPE>
MElementType TPZCompElHCurl<TSHAPE>::Type() {
    return TSHAPE::Type();
}


template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::NConnects() const {
    return TSHAPE::NSides - TSHAPE::NCornerNodes;
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::SideConnectLocId(int con,int side) const {
#ifdef PZDEBUG
    if(TSHAPE::SideDimension(side)< TSHAPE::Dimension - 2 || con >= NSideConnects(side)) {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " there is no connect with local (side) id "<<con<<" on side "<<side << std::endl;
        PZError<<sout.str();
#ifdef PZ_LOG
        LOGPZ_ERROR(logger,sout.str())
#endif
        DebugStop();
        return -1;
    }
#endif
    int conSide = -1;
    TPZStack<int> sideClosure;
    TSHAPE::LowerDimensionSides(side,sideClosure);
    sideClosure.Push(side);
    int iCon = -1;
    for(auto &subSide :sideClosure){
        if(TSHAPE::SideDimension(subSide)) iCon++;
        if(iCon == con) {
            conSide = subSide;
            break;
        }
    }
    if(conSide<0){
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " ERROR: could not find subside associated with connect "<<con<<" on side "<<side << std::endl;
        PZError<<sout.str();
#ifdef PZ_LOG
        LOGPZ_ERROR(logger,sout.str())
#endif
        DebugStop();
        return -1;
    }
    return conSide-TSHAPE::NCornerNodes;
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::NSideConnects(int side) const{
#ifdef PZDEBUG
    if(side <0 || side >= TSHAPE::NSides){
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << "Side: " << side <<"unhandled case\n";
        PZError<<sout.str();
#ifdef PZ_LOG
        LOGPZ_ERROR(logger,sout.str())
#endif
    }
#endif
    int nCons = 0;
    TPZStack<int> sideClosure;
    TSHAPE::LowerDimensionSides(side,sideClosure);
    sideClosure.Push(side);
    for(auto &subSide :sideClosure){
        if(TSHAPE::SideDimension(subSide)) nCons++;
    }
    return nCons;
}

template<class TSHAPE>
int64_t TPZCompElHCurl<TSHAPE>::ConnectIndex(int con) const{
#ifndef PZNODEBUG
    if(con <0 || con >= this->NConnects()) {
        std::cout <<__PRETTY_FUNCTION__ <<" wrong parameter connect " << con <<
                  " NConnects " << this-> NConnects() << std::endl;
        DebugStop();
        return -1;
    }

#endif
    return this->fConnectIndexes[con];
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::SetConnectIndex(int i, int64_t connectindex){
#ifndef PZNODEBUG
    if(i<0 || i>= this->NConnects()) {
        std::cout << " TPZCompElHCurl<TSHAPE>::SetConnectIndex index " << i <<
                  " out of range\n";
        DebugStop();
        return;
    }
#endif
    this-> fConnectIndexes[i] = connectindex;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << std::endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::NConnectShapeF(int connect, int order) const
{
    switch (fhcurlfam)
    {
    case HCurlFamily::EHCurlStandard:
        return TPZShapeHCurl<TSHAPE>::ComputeNConnectShapeF(connect, order);
        break;
    case HCurlFamily::EHCurlNoGrads:
        return TPZShapeHCurlNoGrads<TSHAPE>::ComputeNConnectShapeF(connect, order);
        break;
    }/**there is no default case on purpose, because now the compiler
      will warn us if a new hcurl family is added*/
    return -1;
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::ConnectOrder(int connect) const {
    if (connect < 0 || connect >= this->NConnects()) {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__<<std::endl;
        sout << "Connect index out of range connect " << connect << " nconnects " << NConnects();
        PZError<<sout.str()<<std::endl;
#ifdef PZ_LOG
        LOGPZ_ERROR(logger, sout.str())
#endif
        DebugStop();
        return -1;
    }
    TPZConnect &c = this-> Connect(connect);
    return c.Order();
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::EffectiveSideOrder(int side) const{
	if(!NSideConnects(side)) return -1;
	const auto connect = this->MidSideConnectLocId( side);
	if(connect >= 0 || connect < NConnects()){
        return ConnectOrder(connect);
	}
    else{
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__<<std::endl;
        sout << "Connect index out of range connect " << connect << " nconnects " << NConnects();
        PZError<<sout.str()<<std::endl;
#ifdef PZ_LOG
        LOGPZ_ERROR(logger, sout.str())
#endif
        DebugStop();
    }
	return -1;
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
	ord.Resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
	for(auto i=0; i<TSHAPE::NSides-TSHAPE::NCornerNodes; i++) {
		ord[i] = this->Connect(i).Order();
	}
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::SetSideOrder(int side, int order){
    const int connect= this->MidSideConnectLocId(side);
    if(connect<0 || connect > this-> NConnects()) {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__<<"Bad parameter side " << side << " order " << order;
        PZError << sout.str()<<std::endl;
#ifdef PZ_LOG
        LOGPZ_ERROR(logger,sout.str())
#endif
        DebugStop();
        return;
    }
    TPZConnect &c = this->Connect(connect);
    c.SetOrder(order,this->fConnectIndexes[connect]);
    int64_t seqnum = c.SequenceNumber();
    const int nStateVars = [&](){
        TPZMaterial * mat =this-> Material();
        if(mat) return mat->NStateVariables();
        else {
#ifdef PZ_LOG
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__<<"\tAssuming only one state variable since no material has been set";
            LOGPZ_DEBUG(logger,sout.str())
#endif
            return 1;
        }
    }();
    c.SetNState(nStateVars);
    const int nshape =this->NConnectShapeF(connect,order);
    c.SetNShape(nshape);
    this-> Mesh()->Block().Set(seqnum,nshape*nStateVars);
    this->AdjustIntegrationRule();
    //for the hcurl and hdiv spaces to be compatible, the approximation order of a face must be max(k,ke), where
    //k is the (attempted) order of the face, and ke the maximum order of the edges contained in it.
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::InitMaterialData(TPZMaterialData &data){
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
        LOGPZ_DEBUG(logger,"Initializing MaterialData of TPZCompElHCurl")
    }
#endif
    data.fShapeType = TPZMaterialData::MShapeFunctionType::EVecShape;
    TPZShapeData & shapedata = data;

    TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes,0);
    TPZGeoEl *ref = this->Reference();
    for(auto i=0; i<TSHAPE::NCornerNodes; i++) {
        ids[i] = ref->NodePtr(i)->Id();
    }
    
    auto &conOrders = shapedata.fHDivConnectOrders;
    constexpr auto nConnects = TSHAPE::NSides - TSHAPE::NCornerNodes;
    conOrders.Resize(nConnects,-1);
    for(auto i = 0; i < nConnects; i++){
        conOrders[i] = this->EffectiveSideOrder(i + TSHAPE::NCornerNodes);
    }

    switch (fhcurlfam)
    {
    case HCurlFamily::EHCurlStandard:
        TPZShapeHCurl<TSHAPE>::Initialize(ids, conOrders, shapedata);
        break;
    case HCurlFamily::EHCurlNoGrads:
        TPZShapeHCurlNoGrads<TSHAPE>::Initialize(ids, conOrders, shapedata);
        break;
    }/**there is no default case on purpose, because now the compiler
      will warn us if a new hcurl family is added*/
    

    //resizing of TPZMaterialData structures

    constexpr int dim = TSHAPE::Dimension;
    constexpr int curldim = [dim](){
        if constexpr (dim == 1) return 1;
        else{
            return 2*dim - 3;//1 for 2D 3 for 3D
        }
    }();
    const int nshape = this->NShapeF();
    
    auto &phi = data.phi;
    auto &curlphi = data.curlphi;
    
    phi.Redim(nshape,3);
    curlphi.Redim(curldim,nshape);
    
    data.axes.Redim(dim,3);
    data.jacobian.Redim(dim,dim);
    data.jacinv.Redim(dim,dim);
    data.x.Resize(3);
// #ifdef PZ_LOG
//     if(logger.isDebugEnabled()){
// 		std::stringstream sout;
// 		sout << "Vector/Shape indexes \n";
//         for (int i = 0; i < shapedata.fVecShapeIndex.size(); i++) {
//             sout << i << '|' << shapedata.fVecShapeIndex[i] << " ";
//         }
//         sout << std::endl;
// 		LOGPZ_DEBUG(logger,sout.str())
// 	}
// #endif

}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data) {

    constexpr int dim = TSHAPE::Dimension;
    constexpr int curldim = [dim](){
        if constexpr (dim == 1) return 1;
        else{
            return 2*dim - 3;//1 for 2D 3 for 3D
        }
    }();

    const int nshape = this->NShapeF();
    TPZFNMatrix<dim*80,REAL> phiref(dim,nshape);
    TPZFNMatrix<curldim*80,REAL> curlphiref(curldim,nshape);

    TPZShapeData &shapedata = data;
    switch (fhcurlfam)
    {
    case HCurlFamily::EHCurlStandard:
        TPZShapeHCurl<TSHAPE>::Shape(qsi, shapedata, phiref, curlphiref);
        break;
    case HCurlFamily::EHCurlNoGrads:
        TPZShapeHCurlNoGrads<TSHAPE>::Shape(qsi,shapedata,phiref,curlphiref);
        break;
    }
    /**there is no default case on purpose, because now the compiler
      will warn us if a new hcurl family is added*/
    
    //these are resized in InitMaterialData
    auto &phi = data.phi;
    auto &curlphi = data.curlphi;
    
    TransformShape(phiref, data.detjac, data.jacinv, data.axes, phi);
    switch(dim){
    case 1:
        TransformCurl<1>(curlphiref, data.detjac, data.jacobian, curlphi);
        break;
    case 2:
        TransformCurl<2>(curlphiref, data.detjac, data.jacobian, curlphi);
        break;
    case 3:
        TransformCurl<3>(curlphiref, data.detjac, data.jacobian, curlphi);
        break;
    }
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
    //this method is not really useful right now
    TPZShapeData data;

    TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes,0);
    TPZGeoEl *ref = this->Reference();
    for(auto i=0; i<TSHAPE::NCornerNodes; i++) {
        ids[i] = ref->NodePtr(i)->Id();
    }
    
    constexpr auto nConnects{TSHAPE::NSides-TSHAPE::NCornerNodes};
    TPZManVector<int, nConnects> conorders(nConnects,0);

    for(int ic = 0; ic < nConnects; ic++){
        conorders[ic] = this->ConnectOrder(ic);
    }

    TPZShapeHCurl<TSHAPE>::Initialize(ids,conorders,data);
    const auto nShape = this->NShapeF();

    constexpr int dim = TSHAPE::Dimension;

    constexpr int curldim = [dim](){
        if constexpr (dim == 1) return 1;
        else{
            return 2*dim - 3;//1 for 2D 3 for 3D
        }
    }();
    
    phi.Redim(dim,nShape);
    dphi.Redim(curldim, nShape);
    
    TPZShapeHCurl<TSHAPE>::Shape(pt, data, phi, dphi);
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &curlphi)
{
    constexpr int dim = TSHAPE::Dimension;
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
    TPZManVector<int,TSHAPE::NSides> ordHCurl(nSideConnects,-1);
    for (auto is=nSideNodes; is< nContainedSides; is++) {
        const int subSide = TSHAPE::ContainedSideLocId(side,is);
        ordHCurl[is-nSideNodes] = this->EffectiveSideOrder(subSide);
    }
    
    TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(nSideNodes,0);
    TPZGeoEl *ref = this->Reference();
    for (auto is=0; is< nSideNodes; is++) {
        const int subSide = TSHAPE::ContainedSideLocId(side,is);
        ids[is] = ref->NodePtr(subSide)->Id();
    }

    TPZShapeData shapedata;

    
    int sidedim{-1};
    int nshape{-1};
    switch(sidetype){
    case EOned:
        TPZShapeHCurl<pzshape::TPZShapeLinear>::Initialize(ids,ordHCurl,shapedata);
        nshape = TPZShapeHCurl<pzshape::TPZShapeLinear>::NHCurlShapeF(shapedata);
        sidedim = 1;
        break;
    case ETriangle:
        TPZShapeHCurl<pzshape::TPZShapeTriang>::Initialize(ids,ordHCurl,shapedata);
        nshape = TPZShapeHCurl<pzshape::TPZShapeTriang>::NHCurlShapeF(shapedata);
        sidedim = 2;
        break;
    case EQuadrilateral:
        TPZShapeHCurl<pzshape::TPZShapeQuad>::Initialize(ids,ordHCurl,shapedata);
        nshape = TPZShapeHCurl<pzshape::TPZShapeQuad>::NHCurlShapeF(shapedata);
        sidedim = 2;
        break;
    default:
        PZError<<__PRETTY_FUNCTION__
               <<"\n invalid side type.Aborting...\n";
        DebugStop();
    }

    const int curldim = [sidedim](){
        if (sidedim == 1) return 1;
        else{
            return 2*sidedim - 3;//1 for 2D 3 for 3D
        }
    }();
    
    TPZFMatrix<REAL> phiref(sidedim,nshape);
    TPZFMatrix<REAL> curlphiref(curldim,nshape);
    
    switch(sidetype){
    case EOned:
        TPZShapeHCurl<pzshape::TPZShapeLinear>::Shape(point, shapedata, phiref, curlphiref);
        break;
    case ETriangle:
        TPZShapeHCurl<pzshape::TPZShapeTriang>::Shape(point, shapedata, phiref, curlphiref);
        break;
    case EQuadrilateral:
        TPZShapeHCurl<pzshape::TPZShapeQuad>::Shape(point, shapedata, phiref, curlphiref);
        break;
    default:
        PZError<<__PRETTY_FUNCTION__
               <<"\n invalid side type.Aborting...\n";
        DebugStop();
    }

    //get the jacobian of the side transformation
    TPZGeoElSide gelside = TPZGeoElSide(this->Reference(),side);
    TPZFNMatrix<9,REAL> jac(sideDim,sideDim),jacinv(sideDim,sideDim),axes(sideDim,3);
    REAL detjac = 0;
    gelside.Jacobian(point, jac, axes, detjac, jacinv);


    phi.Redim(nshape,3);
    curlphi.Redim(curldim,nshape);
    
    TransformShape(phiref, detjac, jacinv, axes, phi);
    switch(sideDim){
    case 1:
        TransformCurl<1>(curlphiref, detjac, jac, curlphi);
        break;
    case 2:
        TransformCurl<2>(curlphiref, detjac, jac, curlphi);
        break;
    case 3:
        TransformCurl<3>(curlphiref, detjac, jac, curlphi);
        break;
    }
    
    
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::TransformShape(const TPZFMatrix<REAL> &phiref,
                                            const REAL detjac,
                                            const TPZFMatrix<REAL> &jacinv,
                                            const TPZFMatrix<REAL> &axes,
                                            TPZFMatrix<REAL> &phi)
{

    //applies covariant piola transform and compute the deformed vectors
    TPZFMatrix<REAL> axest, jacinvt;
    jacinv.Transpose(&jacinvt);
    axes.Transpose(&axest);

    (axest * (jacinvt * phiref)).Transpose(&phi);
}

template<class TSHAPE>
template<int TDIM>
void TPZCompElHCurl<TSHAPE>::TransformCurl(const TPZFMatrix<REAL> &curlphiref,
                                           const REAL detjac,
                                           const TPZFMatrix<REAL> &jacobian,
                                           TPZFMatrix<REAL> &curlphi)
{
    if constexpr(TDIM==3){
        curlphi = jacobian * curlphiref;
        curlphi *= 1./detjac;
    }else {
        curlphi = curlphiref;
        curlphi *= 1./detjac;
    }
}


template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::CreateHCurlConnects(TPZCompMesh &mesh){
    constexpr int nNodes = TSHAPE::NCornerNodes;
    constexpr int nConnects = TSHAPE::NSides - nNodes;
    this->fConnectIndexes.Resize(nConnects);
    for(auto i = 0; i < nConnects; i++){
        const int sideId = nNodes + i;
        this->fConnectIndexes[i] = this->CreateMidSideConnect(sideId);
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "After creating last HCurl connect " << i << std::endl;
            //	this->Print(sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        mesh.ConnectVec()[this->fConnectIndexes[i]].IncrementElConnected();
        this->IdentifySideOrder(sideId);
    }
    this->AdjustIntegrationRule();
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::RestrainSide(int side, TPZInterpolatedElement *large, int neighbourside) {
    const TPZCompElSide thisCompSide(this, side);
    const TPZCompElSide largeCompSide(large, neighbourside);
    TPZGeoElSide thisGeoSide(this->Reference(), side);
    const TPZGeoElSide largeGeoSide = largeCompSide.Reference();
    const auto thisSideDimension = thisGeoSide.Dimension();
    const auto largeSideDimension = largeGeoSide.Dimension();


    TPZInterpolatedElement *largeCel = nullptr;
    if (largeCompSide.Exists()) largeCel = dynamic_cast<TPZInterpolatedElement *> (largeCompSide.Element());
    if (!largeCel) {
        LOGPZ_ERROR(logger, "Exiting RestrainSide - null computational element.");
        return;
    }
    const auto connectLocId = this->MidSideConnectLocId(side);
    if (connectLocId < 0) {
        DebugStop();
    }
    const TPZConnect &myConnect = this->Connect(connectLocId);
    if (myConnect.NShape() == 0) {
        /// no shape functions to restrain
        return;
    }
    if (myConnect.HasDependency() && largeSideDimension > 0) {
        LOGPZ_ERROR(logger, "RestrainSide - unnecessary call to restrainside");
        DebugStop();
    }
    if (largeCel->ConnectIndex(largeCel->MidSideConnectLocId(largeCompSide.Side())) == -1) {
        LOGPZ_ERROR(logger, "Exiting RestrainSide - Side of large element not initialized");
        DebugStop();
        return;
    }
    if (largeSideDimension == 0) {
        LOGPZ_ERROR(logger, "Exiting RestrainSide - dimension of large element is 0");
        DebugStop();
        return;
    }
    TPZTransform<> t(thisSideDimension);
    thisGeoSide.SideTransform3(largeGeoSide, t);
    const auto nSideConnects = NSideConnects(side);
    int maxord = 1;
    for (int sidecon = 0; sidecon < nSideConnects; sidecon++) {
        TPZConnect &c = this->SideConnect(sidecon, side);
        int sideord = c.Order();
        maxord = maxord < sideord ? sideord : maxord;
    }
    const auto largeOrder = large->EffectiveSideOrder(neighbourside);
    const auto sideOrder = this->MidSideConnect(side).Order();
    if (sideOrder < largeOrder && thisSideDimension && largeSideDimension) {
        DebugStop();
    }
    TPZIntPoints *intrule = this->Reference()->CreateSideIntegrationRule(side, maxord * 2);
    if (!intrule) {
        LOGPZ_ERROR(logger, "Exiting RestrainSide - cannot create side integration rule");
        return;
    }
    const auto numint = intrule->NPoints();
    const auto numshape = this->NSideShapeF(side);
    const auto numshapel = large->NSideShapeF(neighbourside);
    TPZFNMatrix<1000, REAL> MSL(numshape, numshapel, 0.);
    TPZFNMatrix<1000, REAL> *M = new TPZFNMatrix<1000, REAL>(numshape, numshape, 0.);
    TPZManVector<REAL, 3> par(thisSideDimension), pointl(largeSideDimension);

    {
        REAL detjac;
        TPZVec<REAL> centerPoint(thisSideDimension,0);
        thisGeoSide.CenterPoint(centerPoint);
        TPZFNMatrix<9,REAL> jac(thisSideDimension,thisSideDimension),jacinv(thisSideDimension,thisSideDimension),axes(thisSideDimension,3);
        thisGeoSide.Jacobian(centerPoint, jac, axes, detjac, jacinv);
        REAL weight;
        TPZFNMatrix<100, REAL> phis(numshape, 3), dphis(3, numshape), phil(numshapel, 3), dphil(3, numshapel);
        TPZFNMatrix<100, REAL> thisTrace(numshape, 3), largeTrace(numshapel, 3);
        for (int it = 0; it < numint; it++) {
            intrule->Point(it, par, weight);
            SideShapeFunction(side, par, thisTrace, dphis);
            t.Apply(par, pointl);
            large->SideShapeFunction(neighbourside, pointl, largeTrace, dphil);
            for (auto in = 0; in < numshape; in++) {
                for (auto jn = 0; jn < numshape; jn++) {
                    REAL dotProduct = 0;
                    for(auto iaxes = 0; iaxes < 3; iaxes ++){
                        dotProduct += thisTrace(in,iaxes) * thisTrace(jn,iaxes);
                    }
                    (*M)(in, jn) += dotProduct * weight;
                }
                for (auto jn = 0; jn < numshapel; jn++) {
                    REAL dotProduct = 0;
                    for(auto iaxes = 0; iaxes < 3; iaxes ++){
                        dotProduct += thisTrace(in,iaxes) * largeTrace(jn,iaxes);
                    }
                    MSL(in, jn) += dotProduct * weight;
                }
            }
        }
    }

#ifdef PZ_LOG_keep
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        M->Print("MSS = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZStepSolver<REAL> MSolve(M);
    MSolve.SetDirect(ELU);

    MSolve.Solve(MSL, MSL);

    const auto thisNumSideNodes = NSideConnects(side);
    const auto largeNumSideNodes = large->NSideConnects(neighbourside);
    TPZBlock MBlocksmall(0, thisNumSideNodes), MBlocklarge(0, largeNumSideNodes);
    for (auto in = 0; in < thisNumSideNodes; in++) {
        int locid = SideConnectLocId(in, side);
        TPZConnect &c = this->Connect(locid);
#ifdef PZDEBUG
        if (NConnectShapeF(locid, c.Order()) != c.NShape()) {
            DebugStop();
        }
#endif
        unsigned int nshape = c.NShape();
        MBlocksmall.Set(in, nshape);
    }
    for (auto in = 0; in < largeNumSideNodes; in++) {
        int locid = large->SideConnectLocId(in, neighbourside);
        TPZConnect &c = large->Connect(locid);
        unsigned int nshape = c.NShape();
#ifdef PZDEBUG
        if (large->NConnectShapeF(locid, c.Order()) != nshape) {
            DebugStop();
        }
#endif
        MBlocklarge.Set(in, nshape);
    }

    MBlocksmall.Resequence();
    MBlocklarge.Resequence();
    TPZFNMatrix<20, REAL> blocknorm(thisNumSideNodes, largeNumSideNodes, 0.);
    for (auto in = 0; in < thisNumSideNodes; in++) {
        int ibl = MBlocksmall.Size(in);
        if (!ibl) continue;
        for (auto jn = 0; jn < largeNumSideNodes; jn++) {
            int jbl = MBlocklarge.Size(jn);
            if (!jbl) continue;
            int i, j;
            int64_t ipos = MBlocksmall.Position(in);
            int64_t jpos = MBlocklarge.Position(jn);
            for (i = 0; i < ibl; i++) for (j = 0; j < jbl; j++) blocknorm(in, jn) += fabs(MSL(ipos + i, jpos + j)) * fabs(MSL(ipos + i, jpos + j));
            blocknorm(in, jn) /= (ibl * jbl);
            blocknorm(in, jn) = sqrt(blocknorm(in, jn));
        }
    }
#ifdef PZDEBUG
    this->CheckConstraintConsistency(side);
#endif
    TPZConnect &inod = this->Connect(this->MidSideConnectLocId(side));
    const auto inodindex = ConnectIndex(this->MidSideConnectLocId(side));
    int64_t ndepend = 0;
    const auto in = thisNumSideNodes - 1;
#ifdef PZDEBUG
    if (MBlocksmall.Size(in) == 0) {
        DebugStop();
    }
#endif
    for (auto jn = 0; jn < largeNumSideNodes; jn++) {
        if (MBlocksmall.Size(in) == 0 || MBlocklarge.Size(jn) == 0) {
            continue;
        }
        int64_t jnodindex = large->SideConnectIndex(jn, neighbourside);
        TPZConnect::TPZDepend *depend = inod.AddDependency(inodindex, jnodindex, MSL, MBlocksmall.Position(in), MBlocklarge.Position(jn),
                                                           MBlocksmall.Size(in), MBlocklarge.Size(jn));
        if (blocknorm(in, jn) < 1.e-8) {
            depend->fDepMatrix.Zero();
        }
        ndepend++;
    }

    if (!ndepend) {
        for (auto jn = 0; jn < largeNumSideNodes; jn++) {
            int64_t jnodindex = large->SideConnectIndex(jn, neighbourside);
            if (MBlocklarge.Size(jn)) {
                inod.AddDependency(inodindex, jnodindex, MSL, MBlocksmall.Position(in), MBlocklarge.Position(jn),
                                   MBlocksmall.Size(in), MBlocklarge.Size(jn));
            }
            ndepend++;
        }
    }
    delete intrule;

#ifdef HUGE_DEBUG
    // restraint matrix should be equal to MSL
    {
        TPZCheckRestraint test(thisCompSide, largeCompSide);
        int64_t imsl, jmsl;
        const int64_t rmsl = MSL.Rows();
        const int64_t cmsl = MSL.Cols();
        const int64_t rtest = test.RestraintMatrix().Rows();
        const int64_t ctest = test.RestraintMatrix().Cols();

        if (rtest != rmsl || ctest != cmsl) {
            std::stringstream sout;
            sout << "Exiting - Restraint matrix side incompatibility: MSL (rows,cols): ( " << rmsl
                 << " , " << cmsl << " )" << " RestraintMatrix (rows,cols): (" << rtest << " , " << ctest << " )\n";
            LOGPZ_ERROR(logger, sout.str())
            return;
        }

        TPZFMatrix<REAL> mslc(MSL);
        mslc -= test.RestraintMatrix();

        REAL normmsl = 0.;
        for (imsl = 0; imsl < rmsl; imsl++) {
            for (jmsl = 0; jmsl < cmsl; jmsl++) {
                normmsl += sqrt(mslc(imsl, jmsl) * mslc(imsl, jmsl));
            }
        }
        if (normmsl > 1.E-6) {
            std::stringstream sout;
            sout << "TPZInterpolatedElement::Error::MSL matrix has non zero norm " << normmsl << "\n";
            mslc.Print("Difference Matrix ", sout);
            for (imsl = 0; imsl < rmsl; imsl++) {
                for (jmsl = 0; jmsl < cmsl; jmsl++) {
                    if (fabs(MSL(imsl, jmsl) - test.RestraintMatrix()(imsl, jmsl)) > 1.E-6) {
                        sout << "msl[ " << imsl << " , " << jmsl << " ] = " << MSL(imsl, jmsl) << "\t "
                             << test.RestraintMatrix()(imsl, jmsl) << std::endl;
                    }
                }
            }
            LOGPZ_ERROR(logger, sout.str())
        }

        if (test.CheckRestraint()) {
            std::stringstream sout;
            sout << "TPZInterpolatedElement::Error::Bad restraints detected\n";
            test.Print(sout);
            test.Diagnose();
            LOGPZ_ERROR(logger, sout.str())
        }
    }
#endif
}

template<class TSHAPE>
template<class TVar>
void TPZCompElHCurl<TSHAPE>::ComputeSolutionHCurlT(
  const TPZFMatrix<REAL> &phiHCurl, const TPZFMatrix<REAL> &curlPhi,
    TPZSolVec<TVar> &sol, TPZSolVec<TVar> &curlSol)
{
    constexpr int dim = TSHAPE::Dimension;
    constexpr int curlDim = [dim](){
        if constexpr (dim == 1) return 1;
        else{
            return 2*dim - 3;//1 for 2D 3 for 3D
        }
    }();
    const int nVar = this->Material()->NStateVariables();
    const int nConnects = this->NConnects();

    TPZFMatrix<TVar> &meshSol = this->Mesh()->Solution();

    long numberSol = meshSol.Cols();
#ifdef PZDEBUG
    if (numberSol != 1 || nVar != 1) {
        DebugStop();
    }
#endif

    sol.Resize(numberSol);
    curlSol.Resize(numberSol);

    for (long iSol = 0; iSol < numberSol; iSol++) {
        sol[iSol].Resize(dim);
        sol[iSol].Fill(0);
        curlSol[iSol].Resize(curlDim);
        curlSol[iSol].Fill(0);
    }

    TPZBlock &block = this->Mesh()->Block();
    int ishape = 0;
    for (int iCon = 0; iCon < nConnects; iCon++) {
        TPZConnect *con = &this->Connect(iCon);
        const auto conSeqN = con->SequenceNumber();
        const auto nShapeCon = block.Size(conSeqN);
        const auto pos = block.Position(conSeqN);

        for (int jShape = 0; jShape < nShapeCon; jShape++) {

            for (long iSol = 0; iSol < numberSol; iSol++) {
                for (int coord = 0; coord < dim; coord++) {
                    sol[iSol][coord] +=
                            (TVar)meshSol(pos + jShape, iSol) * phiHCurl(ishape, coord);
                }
                for (int coord = 0; coord < curlDim; coord++) {
                    curlSol[iSol][coord] +=
                            (TVar)meshSol(pos + jShape, iSol) * curlPhi.Get(coord, ishape);
                }
            }
            ishape++;
        }
    }
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::MaxOrder() {
    const int ordh1 = TPZInterpolationSpace::MaxOrder();
    return TPZShapeHCurl<TSHAPE>::MaxOrder(ordh1);
}




#define IMPLEMENTHCURL(TSHAPE)                                          \
    template class TPZCompElHCurl<TSHAPE>;                              \
    template void                                                       \
    TPZCompElHCurl<TSHAPE>::TransformCurl<TSHAPE::Dimension>(const TPZFMatrix<REAL> &curlphiref, \
                                                             const REAL detjac, \
                                                             const TPZFMatrix<REAL> &jacobian, \
                                                             TPZFMatrix<REAL> &curlphi);

IMPLEMENTHCURL(pzshape::TPZShapeLinear)
IMPLEMENTHCURL(pzshape::TPZShapeTriang)
IMPLEMENTHCURL(pzshape::TPZShapeQuad)
IMPLEMENTHCURL(pzshape::TPZShapeCube)
IMPLEMENTHCURL(pzshape::TPZShapeTetra)
IMPLEMENTHCURL(pzshape::TPZShapePrism)

#undef IMPLEMENTHCURL


#define HCURL_EL_NOT_AVAILABLE \
    PZError<<__PRETTY_FUNCTION__;\
    PZError<<"Element not available.\n";\
    PZError<<"Aborting...\n";\
    DebugStop();\
    return nullptr;

TPZCompEl *CreateHCurlBoundPointEl(TPZGeoEl *gel, TPZCompMesh &mesh){HCURL_EL_NOT_AVAILABLE}

TPZCompEl *CreateHCurlBoundLinearEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HCurlFamily hcurlfam)
{
    return new TPZCompElHCurl<pzshape::TPZShapeLinear>(mesh, gel, hcurlfam);
}

TPZCompEl *CreateHCurlBoundQuadEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HCurlFamily hcurlfam)
{
    return new TPZCompElHCurl<pzshape::TPZShapeQuad>(mesh, gel, hcurlfam);
}

TPZCompEl *CreateHCurlLinearEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HCurlFamily hcurlfam)
{
    return new TPZCompElHCurl<pzshape::TPZShapeLinear>(mesh, gel, hcurlfam);
}

TPZCompEl *CreateHCurlTriangleEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HCurlFamily hcurlfam)
{
    return new TPZCompElHCurl<pzshape::TPZShapeTriang>(mesh, gel, hcurlfam);
}

TPZCompEl *CreateHCurlQuadEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HCurlFamily hcurlfam)
{
    return new TPZCompElHCurl<pzshape::TPZShapeQuad>(mesh, gel, hcurlfam);
}

TPZCompEl *CreateHCurlTetraEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HCurlFamily hcurlfam)
{
    return new TPZCompElHCurl<pzshape::TPZShapeTetra>(mesh, gel, hcurlfam);
}

TPZCompEl *CreateHCurlCubeEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HCurlFamily hcurlfam)
{
    return new TPZCompElHCurl<pzshape::TPZShapeCube>(mesh, gel, hcurlfam);
}

TPZCompEl *CreateHCurlPrismEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HCurlFamily hcurlfam)
{
    return new TPZCompElHCurl<pzshape::TPZShapePrism>(mesh, gel, hcurlfam);
}

TPZCompEl *CreateHCurlPyramEl(TPZGeoEl *gel, TPZCompMesh &mesh) {
  HCURL_EL_NOT_AVAILABLE
}

#undef HCURL_EL_NOT_AVAILABLE

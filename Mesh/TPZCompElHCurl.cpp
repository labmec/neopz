/**
 * @file
 * @brief Contains the implementation of the TPZCompElHCurlmethods.
 */
#include "TPZCompElHCurl.h"

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
#endif
/*********************************************************************************************************
                                       TPZHCurlAuxClass methods
 *********************************************************************************************************/

void TPZHCurlAuxClass::ComputeShape(const TPZVec<std::pair<int, int64_t>> &vecShapeIndex, const TPZFMatrix<REAL> &phi,
                         const TPZMatrix<REAL> &deformedDirections, TPZMatrix<REAL> &phiHCurl){
    const auto nFuncs = vecShapeIndex.NElements();
    constexpr auto dim = 3;//always 3D in the deformed element
    phiHCurl.Redim(nFuncs,3);
    for(auto iFunc = 0; iFunc < nFuncs; iFunc++){
        const auto vIndex = vecShapeIndex[iFunc].first;
        const auto sIndex = vecShapeIndex[iFunc].second;
        for(auto iDim = 0; iDim < dim; iDim++){
            phiHCurl(iFunc,iDim) = phi.GetVal(sIndex,0) * deformedDirections.GetVal(iDim,vIndex);
        }
    }
}
TPZHCurlAuxClass::EHCurlFamily TPZHCurlAuxClass::hCurlFamily = TPZHCurlAuxClass::EHCurlFamily::EFullOrder;

/*********************************************************************************************************
                                       TPZCompElHCurl methods
 *********************************************************************************************************/

template<class TSHAPE>
TPZCompElHCurl<TSHAPE>::TPZCompElHCurl(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElHCurl::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel,index,1),
        fSidePermutation(TSHAPE::NSides - TSHAPE::NCornerNodes,1),
        fMasterDirections(TSHAPE::Dimension,TSHAPE::Dimension * TSHAPE::NSides,0){
    constexpr int nNodes = TSHAPE::NCornerNodes;
    gel->SetReference(this);
    this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
    /*************************************************************************************
     THE CONNECTS SHOULD BE CREATED IN THE DERIVED CLASS'S CONSTRUCTOR CALLING THE METHOD
                    TPZCompElHCurl<TSHAPE>::CreateHCurlConnects
     ************************************************************************************/
    //compute transform ids for all sides
    TPZVec<int64_t> nodes(nNodes, 0);
    for (auto i = 0; i < nNodes; i++) nodes[i] = gel->NodeIndex(i);
    //computing transformation id for sides
    for(auto iSide = 0 ; iSide < TSHAPE::NSides - TSHAPE::NCornerNodes; iSide++){
        fSidePermutation[iSide] = TSHAPE::GetTransformId(nNodes + iSide, nodes);
    }
    TPZFMatrix<REAL> gradX(TSHAPE::Dimension, TSHAPE::Dimension, 0);
    for (auto x = 0; x < TSHAPE::Dimension; x++) gradX(x, x) = 1;
    TSHAPE::ComputeHCurlDirections(gradX,fMasterDirections,fSidePermutation);
}

template<class TSHAPE>
TPZCompElHCurl<TSHAPE>::TPZCompElHCurl(TPZCompMesh &mesh, const TPZCompElHCurl<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHCurl::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy), fSidePermutation(copy.fSidePermutation), fMasterDirections(copy.fMasterDirections)
{

}

template<class TSHAPE>
TPZCompElHCurl<TSHAPE>::TPZCompElHCurl(TPZCompMesh &mesh,
									 const TPZCompElHCurl<TSHAPE> &copy,
									 std::map<int64_t,int64_t> & gl2lcConMap,
									 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHCurl::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fSidePermutation(copy.fSidePermutation),
fMasterDirections(copy.fMasterDirections)
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
TPZIntelGen<TSHAPE>(), fSidePermutation(TSHAPE::NSides - TSHAPE::NCornerNodes,-1),
fMasterDirections(TSHAPE::Dimension,TSHAPE::Dimension * TSHAPE::NSides,0)
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
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
//    for (int side=TSHAPE::NCornerNodes; side < TSHAPE::NSides; side++) {
//        TPZGeoElSide gelside(this->Reference(),side);
//        TPZStack<TPZCompElSide> celstack;
//        TPZCompElSide largecel = gelside.LowerLevelCompElementList2(0);
//        if (largecel) {
//            int cindex = this->MidSideConnectLocId(side);
//            TPZConnect &c = this->Connect(cindex);
//            c.RemoveDepend();
//        }
//        if (gelside.Element()){
//            gelside.HigherLevelCompElementList3(celstack, 0, 1);
//        }
//        int64_t ncel = celstack.size();
//        for (int64_t el=0; el<ncel; el++) {
//            TPZCompElSide celside = celstack[el];
//            TPZCompEl *celsmall = celside.Element();
//            TPZGeoEl *gelsmall = celsmall->Reference();
//            if (gelsmall->SideDimension(celside.Side()) != gel->Dimension()-1) {
//                continue;
//            }
//            TPZInterpolatedElement *intelsmall = dynamic_cast<TPZInterpolatedElement *>(celsmall);
//            if (!intelsmall) {
//                DebugStop();
//            }
//            int cindex = intelsmall->MidSideConnectLocId(celside.Side());
//            TPZConnect &c = intelsmall->Connect(cindex);
//            c.RemoveDepend();
//        }
//    }
    if (gel){
        gel->ResetReference();
    }
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::ClassId() const{
    return TPZCompElHCurl<TSHAPE>::StaticClassId();
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::StaticClassId(){
    return Hash("TPZCompElHCurl") ^ TPZIntelGen<TSHAPE>().ClassId() << 1;
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
    constexpr int dim = 0;
    return TSHAPE::NSides - TSHAPE::NumSides(dim);
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
#ifndef NODEBUG
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
#ifndef NODEBUG
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
    //setting the type of shape functions as scalar functions + constant vector fields
    data.fShapeType = TPZMaterialData::EVecandShape;

    data.fMasterDirections = this->fMasterDirections;

    //computes the index that will associate each scalar function to a constant vector field
    constexpr auto nConnects = TSHAPE::NSides - TSHAPE::NCornerNodes;
    TPZManVector<int,nConnects> connectOrders(nConnects,-1);
    for(auto i = 0; i < nConnects; i++){
        connectOrders[i] = this->EffectiveSideOrder(i + TSHAPE::NCornerNodes);
    }
    IndexShapeToVec(data.fVecShapeIndex, connectOrders);


#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
		std::stringstream sout;
		sout << "Vector/Shape indexes \n";
        for (int i = 0; i < data.fVecShapeIndex.size(); i++) {
            sout << i << '|' << data.fVecShapeIndex[i] << " ";
        }
        sout << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif

}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi){

    {
        const bool needsSol = data.fNeedsSol;
        data.fNeedsSol = false;
        TPZIntelGen<TSHAPE>::ComputeRequiredData(data,qsi);//in this method, Shape will be called
        data.fNeedsSol = needsSol;
    }
    constexpr auto nVec{TSHAPE::Dimension*TSHAPE::NSides};
    data.fDeformedDirections.Resize(3,nVec);
    constexpr auto dim{TSHAPE::Dimension};
    //applies covariant piola transform and compute the deformed vectors
    for (auto iVec = 0; iVec < nVec; iVec++) {
        TPZManVector<REAL, 3> tempDirection(dim, 0);
        for (auto i = 0; i < dim; i++) {
            //covariant piola transform: J^{-T}
            tempDirection[i] = 0;
            for (auto j = 0; j < dim; j++) tempDirection[i] += data.jacinv(j, i) * data.fMasterDirections(j, iVec);
        }
        for (auto i = 0; i < 3; i++) {
            data.fDeformedDirections(i, iVec) = 0;
            for (auto j = 0; j < dim; j++) data.fDeformedDirections(i, iVec) += data.axes(j, i) * tempDirection[j];
        }
    }
    /******************************************************************************************************************
    * at this point, we already have the basis functions on the deformed element, since we have data.phi,
    * data.fVecShapeIndex and data.fDeformedDirections. Now it is time to compute the curl, which will be stored in
    * data.curlphi.
    *******************************************************************************************************************/
    data.curlphi.Redim(2*dim - 3 > 0 ? 2*dim - 3 : 1, this->NShapeF());
    ComputeCurl(data.fVecShapeIndex,data.dphi,this->fMasterDirections,data.jacobian,data.detjac,data.axes,data.curlphi);
    if (data.fNeedsSol) {
        ComputeSolution(qsi, data);
    }
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {

	TPZManVector<int64_t,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
    TPZGeoEl *ref = this->Reference();
    for(auto i=0; i<TSHAPE::NCornerNodes; i++) {
        id[i] = ref->NodePtr(i)->Id();
    }
    constexpr auto nConnects{TSHAPE::NSides-TSHAPE::NCornerNodes};
    TPZManVector<int, nConnects> ord(nConnects,0);
    CalculateSideShapeOrders(ord);
    const int nShape = TSHAPE::NShapeF(ord);

    phi.Redim(nShape, 1);
    dphi.Redim(TSHAPE::Dimension, nShape);
    TSHAPE::Shape(pt,id,ord,phi,dphi);
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data){
    TPZCompElHCurl<TSHAPE>::ComputeRequiredData(data,qsi);
    ComputeSolutionHCurl(qsi, data.fVecShapeIndex, data.fDeformedDirections, data.phi, data.curlphi,
            data.sol, data.curlsol);
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
#ifdef PZ_LOG
        LOGPZ_WARN(logger, "RestrainSide - unnecessary call to restrainside");
#endif
        DebugStop();
    }
    if (largeCel->ConnectIndex(largeCel->MidSideConnectLocId(largeCompSide.Side())) == -1) {
#ifdef PZ_LOG
        LOGPZ_ERROR(logger, "Exiting RestrainSide - Side of large element not initialized");
#endif
        DebugStop();
        return;
    }
    if (largeSideDimension == 0) {
#ifdef PZ_LOG
        LOGPZ_ERROR(logger, "Exiting RestrainSide - dimension of large element is 0");
#endif
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
void TPZCompElHCurl<TSHAPE>::ComputeSolutionHCurl(const TPZVec<REAL> &qsi, const TPZVec<std::pair<int,int64_t> > &vecShapeIndex,
                          const TPZFMatrix<REAL> &deformedDirections, TPZFMatrix<REAL> &phi,
                          TPZFMatrix<REAL> &curlPhi, TPZSolVec &sol, TPZSolVec &curlSol){
    TPZFMatrix<REAL> phiHCurl;
    TPZHCurlAuxClass::ComputeShape(vecShapeIndex,phi,deformedDirections,phiHCurl);
    constexpr int dim = TSHAPE::Dimension;
    constexpr int curlDim = 2*dim - 3;//1 for 2D 3 for 3D
    const int nVar = this->Material()->NStateVariables();
    const int nConnects = this->NConnects();

    TPZFMatrix<STATE> &meshSol = this->Mesh()->Solution();

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
                            (STATE)meshSol(pos + jShape, iSol) * phiHCurl(ishape, coord);
                }
                for (int coord = 0; coord < curlDim; coord++) {
                    curlSol[iSol][coord] +=
                            (STATE)meshSol(pos + jShape, iSol) * curlPhi(coord, ishape);
                }
            }
            ishape++;
        }
    }
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::CalcShapeSideTraces(const int side, const TPZFMatrix<REAL> &phi,
//        const TPZFMatrix<REAL> &curlPhi,
        TPZFMatrix<REAL> &phiTrace
//        , TPZFMatrix<REAL> &curlTrace
        ) const {
    const TPZGeoElSide geoSide(this->Reference(), side);
    const auto sideDim = geoSide.Dimension();

    //calculate the axes associated with the side side
    TPZFNMatrix<9,REAL> axes(sideDim,3);
    {
        REAL detjac;
        TPZVec<REAL> centerPoint(sideDim,0);
        geoSide.CenterPoint(centerPoint);
        TPZFNMatrix<9,REAL> jac(sideDim,sideDim),jacinv(sideDim,sideDim);
        geoSide.Jacobian(centerPoint, jac, axes, detjac, jacinv);
    }

    const auto numshape = phi.Rows();
    TPZFMatrix<REAL> phiTraceAxes(numshape,sideDim,0);
    for(auto iphi = 0; iphi < numshape; iphi++){
        for(auto iaxes = 0; iaxes < sideDim; iaxes ++){
            phiTrace(iphi,iaxes) = 0;
            for(auto jaxes = 0; jaxes < 3; jaxes ++){
                phiTrace(iphi,iaxes) += axes(iaxes,jaxes) * phi.GetVal(iphi,jaxes);
            }
        }
    }
    phiTrace.Resize(numshape,3);
    //converting back to XYZ
    TPZAxesTools<REAL>::Axes2XYZ(phiTraceAxes, phiTrace, axes, false);
//    //if the functions are being calculated over a face, the curl of the shape funcs are 3d as well. otherwise, no transformation
//    //should be performed
//    if(sideDim == 2){
//        TPZFMatrix<REAL> curlTraceAxes(numshape,sideDim,0);
//        for(auto iphi = 0; iphi < numshape; iphi++){
//            for(auto iaxes = 0; iaxes < sideDim; iaxes ++){
//                curlTraceAxes(iphi,iaxes) = 0;
//                for(auto jaxes = 0; jaxes < 3; jaxes ++){
//                    curlTraceAxes(iphi,iaxes) += axes(iaxes,jaxes) * curlPhi.GetVal(jaxes,iphi);
//                }
//            }
//        }
//        TPZAxesTools<REAL>::Axes2XYZ(curlTraceAxes, curlTrace, axes, true);
//    }
//    else{
//        curlTrace = curlPhi;
//    }
}

template <>
void TPZHCurlAuxClass::ComputeCurl<3>(const TPZVec<std::pair<int, int64_t>> &vecShapeIndex, const TPZFMatrix<REAL> &dphi,
                                         const TPZFMatrix<REAL> &masterDirections, const TPZFMatrix<REAL> &jacobian,
                                         REAL detJac, const TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &curlPhi) {
    constexpr auto dim = 3;
    const auto nShapeFuncs = vecShapeIndex.size();
    const REAL jacInv = 1/detJac;
    const TPZFMatrix<REAL> * dphiPtr = &dphi;
    const TPZFMatrix<REAL> * masterDirectionsPtr = &masterDirections;
    const TPZFMatrix<REAL> * jacobianPtr = &jacobian;
    const TPZFMatrix<REAL> * axesPtr = &axes;
    for(auto iShapeFunc = 0; iShapeFunc < nShapeFuncs; iShapeFunc++) {
        const auto iVec = vecShapeIndex[iShapeFunc].first;
        const auto iShape = vecShapeIndex[iShapeFunc].second;
        TPZManVector<REAL, dim> gradPhiCrossDirections(dim, 0);
        for(auto ix = 0; ix < dim; ix++) {
            gradPhiCrossDirections[ix] =
                    GETVAL(dphiPtr,dim,(ix+1)%dim,iShape) * GETVAL(masterDirectionsPtr,dim, (ix + 2) % dim,iVec) -
                    GETVAL(masterDirectionsPtr,dim,(ix + 1) % dim,iVec) * GETVAL(dphiPtr,dim,(ix + 2) % dim,iShape);
        }
        TPZManVector<REAL, dim> tempCurl(dim, 0);
        for (auto i = 0; i < dim; i++) {
            tempCurl[i] = 0;
            for (auto j = 0; j < dim; j++) tempCurl[i] += GETVAL(jacobianPtr,dim,i,j) * gradPhiCrossDirections[j];
        }
        for (auto i = 0; i < 3; i++) {
            curlPhi(i, iShapeFunc) = 0;
            for (auto j = 0; j < dim; j++) curlPhi(i, iShapeFunc) += jacInv * GETVAL(axesPtr,dim,j,i) * tempCurl[j];
        }
    }
}

template <>
void TPZHCurlAuxClass::ComputeCurl<2>(const TPZVec<std::pair<int, int64_t>> &vecShapeIndex, const TPZFMatrix<REAL> &dphi,
                                      const TPZFMatrix<REAL> &masterDirections, const TPZFMatrix<REAL> &jacobian,
                                      REAL detJac, const TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &curlPhi) {
    const auto nShapeFuncs = vecShapeIndex.size();
    const REAL jacInv = 1/detJac;
    for(auto iShapeFunc = 0; iShapeFunc < nShapeFuncs; iShapeFunc++) {
        const auto iVec = vecShapeIndex[iShapeFunc].first;
        const auto iShape = vecShapeIndex[iShapeFunc].second;
        const REAL gradPhiCrossDirections = dphi.GetVal( 1,iShape) * masterDirections.GetVal(0,iVec) -
                                            masterDirections.GetVal(1,iVec) * dphi.GetVal( 0,iShape);
        curlPhi(0, iShapeFunc) += jacInv * gradPhiCrossDirections;
    }
}

template <>
void TPZHCurlAuxClass::ComputeCurl<1>(const TPZVec<std::pair<int, int64_t>> &vecShapeIndex, const TPZFMatrix<REAL> &dphi,
                                      const TPZFMatrix<REAL> &masterDirections, const TPZFMatrix<REAL> &jacobian,
                                      REAL detJac, const TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &curlPhi) {
    const auto nShapeFuncs = vecShapeIndex.size();
    for(auto iShapeFunc = 0; iShapeFunc < nShapeFuncs; iShapeFunc++) {
        const auto iVec = vecShapeIndex[iShapeFunc].first;
        const auto iShape = vecShapeIndex[iShapeFunc].second;
        curlPhi(0, iShapeFunc) = dphi.GetVal( 0,iShape) * masterDirections.GetVal(0,iVec);
    }
}

template <int TDIM>
void TPZHCurlAuxClass::ComputeCurl(const TPZVec<std::pair<int, int64_t>> &vecShapeIndex, const TPZFMatrix<REAL> &dphi,
                                      const TPZFMatrix<REAL> &masterDirections, const TPZFMatrix<REAL> &jacobian,
                                      REAL detJac, const TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &curlPhi) {
    DebugStop();
}


#define IMPLEMENTHCURL(TSHAPE) \
\
template class TPZCompElHCurl<TSHAPE>;

IMPLEMENTHCURL(pzshape::TPZShapeLinear)
IMPLEMENTHCURL(pzshape::TPZShapeTriang)
IMPLEMENTHCURL(pzshape::TPZShapeQuad)
IMPLEMENTHCURL(pzshape::TPZShapeCube)
IMPLEMENTHCURL(pzshape::TPZShapeTetra)
IMPLEMENTHCURL(pzshape::TPZShapePrism)

#undef IMPLEMENTHCURL

#include <TPZCompElHCurlFull.h>


#define HCURL_EL_NOT_AVAILABLE \
    PZError<<__PRETTY_FUNCTION__;\
    PZError<<"Element not available.\n";\
    PZError<<"Aborting...\n";\
    DebugStop();\
    return nullptr;

TPZCompEl *CreateHCurlBoundPointEl(TPZGeoEl *gel, TPZCompMesh &mesh,
                                   int64_t &index){HCURL_EL_NOT_AVAILABLE}

TPZCompEl *CreateHCurlBoundLinearEl(TPZGeoEl *gel, TPZCompMesh &mesh,
                                    int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapeLinear>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlBoundTriangleEl(TPZGeoEl *gel, TPZCompMesh &mesh,
                                      int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapeTriang>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlBoundQuadEl(TPZGeoEl *gel, TPZCompMesh &mesh,
                                  int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapeQuad>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlLinearEl(TPZGeoEl *gel, TPZCompMesh &mesh,
                               int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapeLinear>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlTriangleEl(TPZGeoEl *gel, TPZCompMesh &mesh,
                                 int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapeTriang>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlQuadEl(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapeQuad>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlTetraEl(TPZGeoEl *gel, TPZCompMesh &mesh,int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapeTetra>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlCubeEl(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapeCube>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlPrismEl(TPZGeoEl *gel, TPZCompMesh &mesh,
                              int64_t &index) {
  switch (TPZHCurlAuxClass::GetHCurlFamily()) {
  case TPZHCurlAuxClass::EHCurlFamily::EFullOrder:
    return new TPZCompElHCurlFull<pzshape::TPZShapePrism>(mesh, gel, index);
    break;
  default:
    HCURL_EL_NOT_AVAILABLE
  }
}

TPZCompEl *CreateHCurlPyramEl(TPZGeoEl *gel, TPZCompMesh &mesh,
                              int64_t &index) {
  HCURL_EL_NOT_AVAILABLE
}

#undef HCURL_EL_NOT_AVAILABLE
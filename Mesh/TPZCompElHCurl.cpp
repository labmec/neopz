/**
 * @file
 * @brief Contains the implementation of the TPZCompElHCurlmethods.
 */
#include <TPZCompElHCurl.h>

#include <TPZMaterial.h>
#include <pzcmesh.h>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHCurl"));
static LoggerPtr loggercurl(Logger::getLogger("pz.mesh.tpzinterpolatedelement.divide"));
#endif


TPZHCurlSettings::EHCurlFamily TPZHCurlSettings::hCurlFamily = TPZHCurlSettings::EHCurlFamily::EFullOrder;

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
void TPZCompElHCurl<TSHAPE>::CreateHCurlConnects(TPZCompMesh &mesh){
    constexpr int nNodes = TSHAPE::NCornerNodes;
    constexpr int nConnects = TSHAPE::NSides - nNodes;
    this->fConnectIndexes.Resize(nConnects);
    for(auto i = 0; i < nConnects; i++){
        const int sideId = nNodes + i;
        this->fConnectIndexes[i] = this->CreateMidSideConnect(sideId);
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
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
    for (int side=TSHAPE::NCornerNodes; side < TSHAPE::NSides; side++) {
        if (TSHAPE::SideDimension(side) != TSHAPE::Dimension-1) {
            continue;
        }
        TPZGeoElSide gelside(this->Reference(),side);
        TPZStack<TPZCompElSide> celstack;
        TPZCompElSide largecel = gelside.LowerLevelCompElementList2(0);
        if (largecel) {
            int cindex = SideConnectLocId(0, side);
            TPZConnect &c = this->Connect(cindex);
            c.RemoveDepend();
        }
        if (gelside.Element()){
            gelside.HigherLevelCompElementList3(celstack, 0, 1);
        }
        int64_t ncel = celstack.size();
        for (int64_t el=0; el<ncel; el++) {
            TPZCompElSide celside = celstack[el];
            TPZCompEl *celsmall = celside.Element();
            TPZGeoEl *gelsmall = celsmall->Reference();
            if (gelsmall->SideDimension(celside.Side()) != gel->Dimension()-1) {
                continue;
            }
            TPZInterpolatedElement *intelsmall = dynamic_cast<TPZInterpolatedElement *>(celsmall);
            if (!intelsmall) {
                DebugStop();
            }
            int cindex = intelsmall->SideConnectLocId(0, celside.Side());
            TPZConnect &c = intelsmall->Connect(cindex);
            c.RemoveDepend();
        }
    }
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
int TPZCompElHCurl<TSHAPE>::SideConnectLocId(int node,int side) const {
#ifdef PZDEBUG
    if(TSHAPE::SideDimension(side)< TSHAPE::Dimension - 2 || node >= NSideConnects(side)) {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " no connects associated with node "<<node<<" and side "<<side << std::endl;
        PZError<<sout.str();
#ifdef LOG4CXX
        LOGPZ_ERROR(logger,sout.str())
#endif
        DebugStop();
        return -1;
    }
#endif

    return node+side-(TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1);
}

template<class TSHAPE>
int TPZCompElHCurl<TSHAPE>::NSideConnects(int side) const{
#ifdef PZDEBUG
    if(side <0 || side >= TSHAPE::NSides){
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << "Side: " << side <<"unhandled case\n";
        PZError<<sout.str();
#ifdef LOG4CXX
        LOGPZ_ERROR(logger,sout.str())
#endif
    }
#endif
    if(TSHAPE::SideDimension(side) == 0) return 0;
    else return 1;
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
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
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
#ifdef LOG4CXX
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
	const auto connect = SideConnectLocId(0, side);
	if(connect >= 0 || connect < NConnects()){
        return ConnectOrder(connect);
	}
    else{
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__<<std::endl;
        sout << "Connect index out of range connect " << connect << " nconnects " << NConnects();
        PZError<<sout.str()<<std::endl;
#ifdef LOG4CXX
        LOGPZ_ERROR(logger, sout.str())
#endif
        DebugStop();
    }
	return -1;
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::SetSideOrder(int side, int order){
    const int connect= SideConnectLocId(0,side);
    if(connect<0 || connect > this-> NConnects()) {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__<<"Bad parameter side " << side << " order " << order;
        PZError << sout.str()<<std::endl;
#ifdef LOG4CXX
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
#ifdef LOG4CXX
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
}

template<class TSHAPE>
void TPZCompElHCurl<TSHAPE>::InitMaterialData(TPZMaterialData &data){
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
#ifdef LOG4CXX
    if(logger->isDebugEnabled()){
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


#ifdef LOG4CXX
    if(logger->isDebugEnabled()){
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
     * data.dphi.
     ******************************************************************************************************************/
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
    //@TODOFran::IMPLEMENT ME PLEASE
}

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"

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

TPZCompEl * CreateHCurlBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	return nullptr;//return new TPZCompElHCurlBound2<TPZShapePoint>(mesh,gel,index);
}

TPZCompEl * CreateHCurlBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	return nullptr;//return new TPZCompElHCurlBound2< TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateHCurlBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return nullptr;//return new TPZCompElHCurlBound2< TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateHCurlBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return nullptr;//return new TPZCompElHCurlBound2< TPZShapeQuad>(mesh,gel,index);
}

TPZCompEl * CreateHCurlLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    switch(TPZHCurlSettings::GetHCurlFamily()){
        case TPZHCurlSettings::EHCurlFamily::EFullOrder:
            return new TPZCompElHCurlFull< pzshape::TPZShapeLinear>(mesh,gel,index);
            break;
        default:
            DebugStop();
    }
}

TPZCompEl * CreateHCurlTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    switch(TPZHCurlSettings::GetHCurlFamily()){
        case TPZHCurlSettings::EHCurlFamily::EFullOrder:
            return new TPZCompElHCurlFull< pzshape::TPZShapeTriang >(mesh,gel,index);
            break;
        default:
            DebugStop();
    }
}

TPZCompEl * CreateHCurlQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	switch(TPZHCurlSettings::GetHCurlFamily()){
	    case TPZHCurlSettings::EHCurlFamily::EFullOrder:
	        return new TPZCompElHCurlFull< pzshape::TPZShapeQuad>(mesh,gel,index);
	        break;
	    default:
	        DebugStop();
	}
}

TPZCompEl * CreateHCurlTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    switch(TPZHCurlSettings::GetHCurlFamily()){
        case TPZHCurlSettings::EHCurlFamily::EFullOrder:
            return new TPZCompElHCurlFull< pzshape::TPZShapeTetra >(mesh,gel,index);
            break;
        default:
            DebugStop();
    }
}

TPZCompEl * CreateHCurlCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	switch(TPZHCurlSettings::GetHCurlFamily()){
	    case TPZHCurlSettings::EHCurlFamily::EFullOrder:
	        return new TPZCompElHCurlFull< pzshape::TPZShapeCube >(mesh,gel,index);
	        break;
	    default:
	        DebugStop();
	}
}

TPZCompEl * CreateHCurlPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	switch(TPZHCurlSettings::GetHCurlFamily()){
	    case TPZHCurlSettings::EHCurlFamily::EFullOrder:
	        return new TPZCompElHCurlFull< pzshape::TPZShapePrism>(mesh,gel,index);
	        break;
	    default:
	        DebugStop();
	}
}

TPZCompEl * CreateHCurlPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	DebugStop();
    return nullptr;
}

//
//template<class TSHAPE>
//MElementType TPZCompElHCurl<TSHAPE>::Type() {
//	return TSHAPE::Type();
//}
//
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::NConnects() const {
//	int dimension = Dimension()-1;
//
//	return TSHAPE::NumSides(dimension) + 1;
//}
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::SetConnectIndex(int i, int64_t connectindex){
//#ifndef NODEBUG
//	if(i<0 || i>= this->NConnects()) {
//		std::cout << " TPZCompElHCurl<TSHAPE>::SetConnectIndex index " << i <<
//		" out of range\n";
//		DebugStop();
//		return;
//	}
//#endif
//	this-> fConnectIndexes[i] = connectindex;
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled())
//	{
//		std::stringstream sout;
//		sout << endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
//		LOGPZ_DEBUG(logger,sout.str())
//	}
//#endif
//}
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::NConnectShapeF(int connect, int order)const
//{
//#ifdef DEBUG
//    if (connect < 0 || connect >= NConnects()) {
//        DebugStop();
//    }
//#endif
////    if (connect < TPZCompElHCurl<TSHAPE>::NConnects()-1) {
////         int64_t connectindex = ConnectIndex(connect);
//////         int order = 0;
//////         if (connectindex >= 0) {
//////             order = this->Connect(connect).Order();
//////         }
//////         else
//////         {
//////             order = this->fPreferredOrder;
//////         }
////         const int nfaces = TSHAPE::NumSides(TSHAPE::Dimension-1);
////         int face = TSHAPE::NSides-nfaces+connect-1;
////         TPZStack<int> lowerdimensionsides;
////         TSHAPE::LowerDimensionSides(face,lowerdimensionsides);
////         int nshape = TSHAPE::NConnectShapeF(face,order);
////         for (int is=0; is<lowerdimensionsides.size(); is++) {
////             nshape += TSHAPE::NConnectShapeF(lowerdimensionsides[is],order);
////         }
////         return nshape;
////     }
////
////     const int nfaces = TSHAPE::NumSides(TSHAPE::Dimension-1);
////     int face = TSHAPE::NSides-nfaces-1;
////     int nvecignore = 0;
////     for (int side = face; side < TSHAPE::NSides-1; side++) {
////         nvecignore += TSHAPE::NContainedSides(side);
////     }
////
////
////     TPZManVector<int,TSHAPE::Dimension*TSHAPE::NSides+1> vecside(TSHAPE::Dimension*TSHAPE::NSides),bilinear(TSHAPE::Dimension*TSHAPE::NSides),directions(TSHAPE::Dimension*TSHAPE::NSides);
////     TSHAPE::GetSideHDivDirections(vecside,directions,bilinear);
//////     if (TSHAPE::Type()==ETriangle||TSHAPE::Type()==ETetraedro) {
////////         pressureorder=this->fPreferredOrder-1;
//////     }
//////     else if (TSHAPE::Type()==ECube||TSHAPE::Type()==EQuadrilateral||TSHAPE::Type()==EPrisma) {
//////         pressureorder=this->fPreferredOrder;
//////     }
//////     else
//////     {
//////         // Tipo nao implementado
//////         DebugStop();
//////     }
////     TPZManVector<int,27> orders(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
////     FillOrder(orders);
////     int nshape = TSHAPE::NShapeF(orders);
////
////     TPZManVector<int64_t, TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes);
////     for (int i=0; i<id.size(); i++) {
////         id[i] = this->Reference()->NodePtr(i)->Id();
////     }
////
////
////     int nexternalvectors = 0;
////     nexternalvectors = nvecignore;
////
////     TPZGenMatrix<int> shapeorders(nshape,3);
////     TSHAPE::ShapeOrder(id, orders, shapeorders);
////    {
////        static int first = 0;
////        if (first==0) {
////            shapeorders.Print("ShapeOrders");
////            first++;
////        }
////    }
////     // VectorSide indicates the side associated with each vector entry
////     TPZManVector<int64_t,30> FirstIndex(TSHAPE::NSides+1);
////     // the first index of the shape functions
////     FirstShapeIndex(FirstIndex);
////     //FirstIndex.Print();
////
////#ifdef LOG4CXX
////     if (logger->isDebugEnabled())
////     {
////         std::stringstream sout;
////         sout << "FirstIndex "<<FirstIndex << std::endl;
////         LOGPZ_DEBUG(logger,sout.str())
////     }
////#endif
////
////     int count = 0;
////
////    int internalorder = this->Connect(connect).Order();
////
////     int64_t nvec = vecside.NElements();
////     for (int locvec = nexternalvectors; locvec<nvec; locvec++)
////     {
////         int side = vecside[locvec];
////         int bil = bilinear[locvec];
////         int dir = directions[locvec];
////
////         MElementType tipo = TSHAPE::Type(side);
////
////         int firstshape = FirstIndex[side];
////         int lastshape = FirstIndex[side+1];
////
////         for (int ish = firstshape; ish<lastshape; ish++)
////         {
////             int sidedimension = TSHAPE::SideDimension(side);
////             int maxorder[3] = {internalorder,internalorder,internalorder};
////             if (bil) {
////                 maxorder[dir]++;
////             }
////             int include=true;
////             for (int d=0; d<sidedimension; d++)
////             {
////                 if (tipo==ETriangle||tipo==ETetraedro)//
////                 {
////                     if (shapeorders(ish,d) > maxorder[d]+1) {
////                         include = false;
////                     }
////                 }
////                 else if(tipo==EQuadrilateral)
////                 {
////                     if (shapeorders(ish,d) > maxorder[d]) {
////                         include = false;
////                     }
////                 }
////                 else if(tipo==ECube)
////                 {
////                     if (shapeorders(ish,d) > maxorder[d]) {
////                         include = false;
////                     }
////                 }
////                 else if (tipo==EPrisma)
////                 {
////                     //DebugStop();
////                     if (shapeorders(ish,d) > maxorder[d]) {
////                         include = false;
////                     }
////                 }
////                 else if (tipo == EOned)
////                 {
////                     if (shapeorders(ish,0) > maxorder[d]) {
////                         include = false;
////                     }
////                 }
////                 else if (tipo == EPoint)
////                 {
////                     DebugStop();
////                 }
////                 else if (tipo == EPiramide)
////                 {
////                     if (shapeorders(ish,d) > maxorder[d]) {
////                         include = false;
////                     }
////                 }
////                 else
////                 {
////                     DebugStop();
////                 }
////             }
////             if (include)
////             {
////                 count++;
////             }
////         }
////     }
////    return count;
////
//
// #ifdef LOG4CXX
//     {
//         std::stringstream sout;
//         sout <<__PRETTY_FUNCTION__<< "unhandled case ";
//         LOGPZ_ERROR(logger,sout.str())
//     }
// #endif
// return -1;
//
// }
//
//
//
//
//////
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::SetIntegrationRule(int ord) {
//	TPZManVector<int,3> order(TSHAPE::Dimension,ord);
//	this->fIntRule.SetOrder(order);
//}
//
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::NSideConnects(int side) const{
//	if(TSHAPE::SideDimension(side)<= Dimension()-2) return 0;
//	if(TSHAPE::SideDimension(side)==Dimension()-1) return 1;
//	if(TSHAPE::SideDimension(side)== Dimension()) {
//        int ncon = 1;
//        return ncon;
//    }
//#ifdef LOG4CXX
//	{
//		std::stringstream sout;
//		sout << __PRETTY_FUNCTION__ << "Side: " << side <<"unhandled case ";
//		LOGPZ_ERROR(logger,sout.str())
//	}
//#endif
//	return -1;
//
//}
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::SideConnectLocId(int node,int side) const {
//#ifdef PZDEBUG
//	if(TSHAPE::SideDimension(side)<= TSHAPE::Dimension - 2 || node >= NSideConnects(side)) {
//		PZError << "TPZCompElHCurl<TSHAPE>::SideConnectLocId no connect associate " <<  endl;
//		return -1;
//	}
//#endif
//
//    return node+side-(TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1);
//}
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::ConnectSideLocId(int connect) const{
//
//    int side = connect+TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1 ;
//    return side;
//}
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
//	ord.Resize(NConnects());
//	int i;
//	for(i=0; i<NConnects(); i++) {
//		ord[i] = ConnectOrder(i);
//	}
//}
//
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::PreferredSideOrder(int side) {
//	if(TSHAPE::SideDimension(side) < Dimension()-1)
//	{
//		PZError << __PRETTY_FUNCTION__ << " side " << side << std::endl;
//	}
//	int connect= SideConnectLocId(0,side);
//	if(connect<0 || connect > NConnects()) {
//		PZError << "TPZCompElHCurl<TSHAPE>::PreferredSideOrder no polynomial associate " <<  endl;
//		return -1;
//	}
//	if(connect<NConnects()) {
//			int order =this->fPreferredOrder;
//			return order;//this->AdjustPreferredSideOrder(side,order);
//	}
//	PZError << "TPZCompElHCurl<TSHAPE>::PreferredSideOrder called for connect = " << connect << "\n";
//	return 0;
//
//}
//
//template<class TSHAPE>
//int64_t TPZCompElHCurl<TSHAPE>::ConnectIndex(int con) const{
//#ifndef NODEBUG
//	if(con<0 || con>= this->NConnects()) {
//		std::cout << "TPZCompElHCurl::ConnectIndex wrong parameter connect " << con <<
//		" NConnects " << this-> NConnects() << std::endl;
//		DebugStop();
//		return -1;
//	}
//
//#endif
//
////    #ifndef NODEBUG
////    	if(con<0) {
////    		std::cout << "TPZCompElHCurl::ConnectIndex wrong parameter connect " << con <<
////    		" NConnects " << this-> NConnects() << std::endl;
////    		DebugStop();
////    		return -1;
////    	}
////
////    #endif
////
////
////    if(con>= this->NConnects())
////    {
////        int con2= con-TSHAPE::NCornerNodes;
////        return this->fConnectIndexes[con2];
////    }
//
//	return this->fConnectIndexes[con];
//}
//
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::SetPreferredOrder(int order)
//{
//		TPZIntelGen<TSHAPE>:: SetPreferredOrder(order);
//	//this->fPreferredOrder = order;
//}
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::SetSideOrder(int side, int order){
//	int connectaux= SideConnectLocId(0,side);
//	if(connectaux<0 || connectaux > this-> NConnects()) {
//		PZError << "TPZCompElHCurl::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
//#ifdef LOG4CXX
//		std::stringstream sout;
//		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
//		LOGPZ_DEBUG(logger,sout.str())
//#endif
//		return;
//	}
//	TPZConnect &c = this->Connect(connectaux);
//    c.SetOrder(order,this->fConnectIndexes[connectaux]);
//    int64_t seqnum = c.SequenceNumber();
//    int nvar = 1;
//    TPZMaterial * mat =this-> Material();
//    if(mat) nvar = mat->NStateVariables();
//    c.SetNState(nvar);
//    int nshape =this->NConnectShapeF(connectaux,order);
//    c.SetNShape(nshape);
//	this-> Mesh()->Block().Set(seqnum,nshape*nvar);
//}
//
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::ConnectOrder(int connect) const{
//	if (connect < 0 || connect >= this->NConnects()){
//#ifdef LOG4CXX
//		{
//			std::stringstream sout;
//			sout << "Connect index out of range connect " << connect <<
//			" nconnects " << NConnects();
//			LOGPZ_DEBUG(logger,sout.str())
//		}
//#endif
//		return -1;
//	}
//
//	if (this->fConnectIndexes[connect] == -1) {
//		std::stringstream sout;
//		sout << __PRETTY_FUNCTION__ << " connect " << connect
//		<< " is not initialized" << std::endl;
//#ifdef LOG4CXX
//		LOGPZ_ERROR(logger,sout.str());
//#else
//		std::cout << sout.str() << std::endl;
//#endif
//		return 0;
//	}
//
//    TPZConnect &c = this-> Connect(connect);
//    return c.Order();
//}
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::EffectiveSideOrder(int side) const
//{
//	if(!NSideConnects(side)) return -1;
//	int corder =SideConnectLocId(0, side);
//	int maxorder = 0;
//	int conectaux;
//	if(corder>=0 || corder <= NConnects()) return ConnectOrder(corder);
//
//    DebugStop();
//	TPZStack< int > high;
//	TSHAPE::HigherDimensionSides(side, high);
//	int highside= high.NElements();
//
//
//	for(int j=0;j<highside;j++)
//	{
//		conectaux =SideConnectLocId(0, high[j]);
//		maxorder = (ConnectOrder(conectaux) > maxorder) ? ConnectOrder(conectaux) : maxorder;
//	}
//
//	return maxorder;
//}
//
///**return the first shape associate to each side*/
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::FirstShapeIndex(TPZVec<int64_t> &Index) const {
//
//    TPZManVector<int> orders(TSHAPE::NSides-TSHAPE::NCornerNodes);
//    FillOrder(orders);
//    Index[0] = 0;
//	for(int iside=0;iside<TSHAPE::NSides;iside++)
//	{
//        int sideorder = 1;
//        if (iside >= TSHAPE::NCornerNodes) {
//            sideorder = orders[iside-TSHAPE::NCornerNodes];
//        }
//        int temp = Index[iside] + TSHAPE::NConnectShapeF(iside,sideorder);
//        Index[iside+1] = temp;
//	}
//
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled()) {
//        std::stringstream sout;
//        sout << "First  Index " << Index;
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
//}
//
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::NFluxShapeF() const{
//    int in,result=0;
//    int nn=TPZCompElHCurl::NConnects();
//    for(in=0;in<nn;in++){
////#ifdef LOG4CXX
////				std::stringstream sout;
////				sout << "conect " << in<< " seq number "<<seqnum<<" num func "<<TPZCompElHCurl::NConnectShapeF(in);
////				LOGPZ_DEBUG(logger,sout.str())
////#endif
//        int order = this->Connect(in).Order();
//        result += TPZCompElHCurl::NConnectShapeF(in,order);
//    }
//
//
////#ifdef LOG4CXX
////    std::stringstream sout;
////    sout << "Num funcoes associada ao fluxo " << result;
////    LOGPZ_DEBUG(logger,sout.str())
////#endif
//    return result;
//
//
//}
//
///**
// * @brief Returns a matrix index of the shape and vector  associate to element
// * @param[in] VectorSide Indicates the side associated with each vector
// * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
// * @param[in] pressureorder Order of the pressure (to select shape functions?)
// */
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide, TPZVec<int> &bilinear, TPZVec<int> &direction, TPZVec<std::pair<int,int64_t> > & IndexVecShape, int pressureorder)
//{
//
//        DebugStop();
//
//}
//
//
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int64_t> > & ShapeAndVec, int pressureorder)
//{
//    DebugStop();
//}
//
//
///**
// * @brief It returns the normal orientation of the reference element by the side.
// * Only side that has dimension larger than zero and smaller than me.
// * @param side: side of the reference elemen
// */
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::GetSideOrient(int side){
//
//    int firstside = TSHAPE::NSides-TSHAPE::NFaces-1;
//    if (side < firstside || side >= TSHAPE::NSides - 1) {
//        DebugStop();
//    }
//    return fSideOrient[side-firstside];
//}
//
///**
// * @brief It set the normal orientation of the element by the side.
// * Only side that has dimension equal to my dimension minus one.
// * @param side: side of the reference elemen
// */
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::SetSideOrient(int side, int sideorient){
//
//    int firstside = TSHAPE::NSides-TSHAPE::NFaces-1;
//    if (side < firstside || side >= TSHAPE::NSides - 1) {
//        DebugStop();
//    }
//    fSideOrient[side-firstside] = sideorient;
//}
//
////compute the values of the shape function of the side
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
//
//    if(side==TSHAPE::NSides || point.size() != TSHAPE::Dimension-1){
//        std::cout<<"Don't have side shape associated to this side";
//        DebugStop();
//    }
//	if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension -1 ){
//		return ;
//	}
//    int ncontained = TSHAPE::NContainedSides(side);
//    int nsideshape = 0;
//    int connectlocid = SideConnectLocId(0, side);
//    int order = this->Connect(connectlocid).Order();
//    int is;
//    for (is=0; is<ncontained; is++) {
//        int ic = TSHAPE::ContainedSideLocId(side,is);
//        nsideshape += TSHAPE::NConnectShapeF(ic,order);
//    }
//#ifdef PZDEBUG
//    if (nsideshape != this->NSideShapeF(side)) {
//        DebugStop();
//    }
//#endif
//
//    TPZGeoEl *gel = this->Reference();
//    //int nc = gel->NCornerNodes();
//    int nsn = TSHAPE::NSideNodes(side);
//    TPZManVector<int64_t,8> id(nsn);
//    for (int ic=0; ic<nsn; ic++) {
//        int locid = TSHAPE::SideNodeLocId(side,ic);
//        id[ic] = gel->Node(locid).Id();
//    }
//
//    //int idsize = id.size();
//    TPZManVector<int,9> permutegather(ncontained);
//    int transformid;
//
//
//    MElementType sidetype = TSHAPE::Type(side);
//    switch (sidetype) {
//        case EOned:
//            transformid = pztopology::TPZLine::GetTransformId(id);
//            pztopology::TPZLine::GetSideHDivPermutation(transformid, permutegather);
//            break;
//        case EQuadrilateral:
//            transformid = pztopology::TPZQuadrilateral::GetTransformId(id);
//            pztopology::TPZQuadrilateral::GetSideHDivPermutation(transformid, permutegather);
//            break;
//        case ETriangle:
//            transformid = pztopology::TPZTriangle::GetTransformId(id);
//            pztopology::TPZTriangle::GetSideHDivPermutation(transformid, permutegather);
//            break;
//        default:
//            DebugStop();
//            break;
//    }
//
//    TPZManVector<int,TSHAPE::NSides> ord(TSHAPE::NSides,order);
//
//    int sidedimension = TSHAPE::SideDimension(side);
//    TPZFNMatrix<50,REAL> philoc(nsideshape,1),dphiloc(sidedimension,nsideshape);
//
//    TSHAPE::SideShape(side,point,id,ord,philoc,dphiloc);
//
//    int ncs = TSHAPE::NContainedSides(side);
//    TPZManVector<int64_t,28> FirstIndex(ncs+1,0);
//    for (int ls=0; ls<ncs; ls++) {
//        int localside = TSHAPE::ContainedSideLocId(side,ls);
//        FirstIndex[ls+1] = FirstIndex[ls]+TSHAPE::NConnectShapeF(localside,order);
//    }
//
//    REAL detjac = 1.;
//    {
//        TPZGeoElSide gelside = TPZGeoElSide(this->Reference(),side);
//        int dim = gel->SideDimension(side);
//        TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim),axes(dim,3);
//        gelside.Jacobian(point, jac, axes, detjac, jacinv);
//    }
//
//    for (int side=0; side < ncs; side++) {
//        int ifirst = FirstIndex[side];
//        int kfirst = FirstIndex[permutegather[side]];
//        int nshape = FirstIndex[side+1]-FirstIndex[side];
//        for (int i=0; i<nshape; i++) {
//            phi(ifirst+i,0) = philoc(kfirst+i,0)/detjac;
//            for (int d=0; d< sidedimension; d++) {
//                dphi(d,ifirst+i) = dphiloc(d,kfirst+i)/detjac;
//            }
//        }
//    }
//
//}
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
//{
//    if (var == 99) {
//        return TPZIntelGen<TSHAPE>::Solution(qsi,var,sol);
//    }
//    TPZMaterialData data;
//	InitMaterialData(data);
//	//this->ComputeSolutionHDiv(data);
//    this->ComputeRequiredData(data,qsi);
//    this->ComputeSolutionHDiv(qsi,data);
//	this->Material()->Solution(data,var,sol);
//}
//

//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes) {
//
//	TPZGeoEl * ref = this->Reference();
//	const int nshape = this->NShapeF();
//	const int dim = ref->Dimension();
//
//    TPZMaterialData data;
//    InitMaterialData(data);
//    data.fNeedsSol = true;
//    ComputeRequiredData(data,qsi);
//    sol = data.sol;
//    dsol = data.dsol;
//    axes = data.axes;
//}
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::ComputeSolutionHDiv(TPZMaterialData &data)
//{
//    const int dim = 3; // Hdiv vectors are always in R3
//    const int nstate = this->Material()->NStateVariables();
//    const int ncon = this->NConnects();
//
//    TPZFMatrix<STATE> &MeshSol = this->Mesh()->Solution();
//
//    int64_t numbersol = MeshSol.Cols();
//    if(numbersol != 1)
//    {
//        DebugStop();
//    }
//    data.sol.Resize(numbersol);
//    data.dsol.Resize(numbersol);
//    data.divsol.Resize(numbersol);
//
//    for (int64_t is=0; is<numbersol; is++)
//    {
//        data.sol[is].Resize(dim*nstate);
//        data.sol[is].Fill(0);
//        data.dsol[is].Redim(dim*nstate, dim);
//        data.divsol[is].Resize(nstate);
//        data.divsol[is].Fill(0.);
//    }
//    TPZFNMatrix<220,REAL> dphix(3,data.dphix.Cols());
//    TPZFMatrix<REAL> &dphi = data.dphix;;
//
//    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, data.axes);
//
//    TPZFMatrix<STATE> GradOfPhiHdiv(dim,dim);
//    GradOfPhiHdiv.Zero();
//
//
//    int normvecRows = data.fDeformedDirections.Rows();
//    int normvecCols = data.fDeformedDirections.Cols();
//    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
//    TPZManVector<TPZFNMatrix<9,REAL>,18> GradNormalvec(normvecCols);
//    for (int i=0; i<GradNormalvec.size(); i++) {
//        GradNormalvec[i].Redim(dim,dim);
//    }
//
//    if (data.fNeedsDeformedDirectionsFad) {
//#ifdef _AUTODIFF
//        for (int e = 0; e < normvecRows; e++) {
//            for (int s = 0; s < normvecCols; s++) {
//                Normalvec(e,s)=data.fDeformedDirectionsFad(e,s).val();
//            }
//        }
//
//        TPZFNMatrix<4,REAL> Grad0(3,3,0.);
//
//        for (int s = 0; s < normvecCols; s++) {
//
//            if (data.fDeformedDirectionsFad(0,s)>0||data.fDeformedDirectionsFad(1,s)>0) {
//                Grad0(0,0)=data.fDeformedDirectionsFad(0,s).fastAccessDx(0);
//                Grad0(0,1)=data.fDeformedDirectionsFad(0,s).fastAccessDx(1);
//                Grad0(1,0)=data.fDeformedDirectionsFad(1,s).fastAccessDx(0);
//                Grad0(1,1)=data.fDeformedDirectionsFad(1,s).fastAccessDx(1);
//            }
//
//            GradNormalvec[s] = Grad0;
//        }
//
//#else
//        DebugStop();
//#endif
//    }else{
//        Normalvec=data.fDeformedDirections;
//    }
//
//    TPZBlock<STATE> &block =this->Mesh()->Block();
//    int ishape=0,ivec=0,counter=0;
//
//    int nshapeV = data.fVecShapeIndex.NElements();
//
//    for(int in=0; in<ncon; in++)
//    {
//        TPZConnect *df = &this->Connect(in);
//        int64_t dfseq = df->SequenceNumber();
//        int dfvar = block.Size(dfseq);
//        // pos : position of the block in the solution matrix
//        int64_t pos = block.Position(dfseq);
//
//        /// ish loops of the number of shape functions associated with the block
//        for(int ish=0; ish<dfvar/nstate; ish++)
//        {
//            ivec    = data.fVecShapeIndex[counter].first;
//            ishape  = data.fVecShapeIndex[counter].second;
//
//
//            // portion of the gradient coming from the gradient of the scalar function
//            for (int e = 0; e < dim; e++) {
//                for (int f = 0; f< dim; f++) {
//                    GradOfPhiHdiv(e,f) = Normalvec(e,ivec)*dphix(f,ishape);
//                }
//            }
//
//            for (int64_t is=0; is<numbersol; is++)
//            {
//                for(int idf=0; idf<nstate; idf++)
//                {
//                    STATE meshsol = MeshSol(pos+ish*nstate+idf,is);
//                    REAL phival = data.phi(ishape,0);
//                    TPZManVector<REAL,3> normal(3);
//
//                    for (int i=0; i<3; i++)
//                    {
//                        if (data.fNeedsDeformedDirectionsFad) {
//#ifdef _AUTODIFF
//                            normal[i] = data.fDeformedDirectionsFad(i,ivec).val();
//#else
//                            DebugStop();
//#endif
//                        }else{
//                            normal[i] = data.fDeformedDirections(i,ivec);
//                        }
//                    }
//
//#ifdef LOG4CXX
//                    if(logger->isDebugEnabled() && abs(meshsol) > 1.e-6)
//                    {
//                        std::stringstream sout;
//                        sout << "meshsol = " << meshsol << " ivec " << ivec << " ishape " << ishape << " x " << data.x << std::endl;
//                        sout << " phi = " << data.phi(ishape,0) << " dphix " << dphix(0,ishape) << " " << dphix(1,ishape) << std::endl;
//                        sout << "normal = " << normal << std::endl;
//                        sout << "GradOfPhiHdiv " << GradOfPhiHdiv << std::endl;
//                        sout << "GradNormalVec " << GradNormalvec[ivec] << std::endl;
//                        LOGPZ_DEBUG(logger,sout.str())
//                    }
//#endif
//                    data.divsol[is][idf] += data.divphi(counter,0)*meshsol;
//                    for (int ilinha=0; ilinha<dim; ilinha++) {
//                        data.sol[is][ilinha+dim*idf] += normal[ilinha]*phival*meshsol;
//                        for (int kdim = 0 ; kdim < dim; kdim++) {
//                            data.dsol[is](ilinha+dim*idf,kdim)+= meshsol * GradOfPhiHdiv(ilinha,kdim);
//                            if(data.fNeedsDeformedDirectionsFad){
//                                data.dsol[is](ilinha+dim*idf,kdim)+=meshsol *GradNormalvec[ivec](ilinha,kdim)*data.phi(ishape,0);
//                            }
//                        }
//                    }
//
//                }
//            }
//            counter++;
//        }
//    }
//
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout << "x " << data.x << " sol " << data.sol[0] << std::endl;
//        data.dsol[0].Print("dsol",sout);
//        sout << "divsol" << data.divsol[0] << std::endl;
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
//
//}
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12)
//{
//
//	bool Is_u1PHI = (u1.Cols() == 1) ? true : false;
//	bool Is_u2PHI = (u2.Cols() == 1) ? true : false;
//
//	if(Is_u1PHI && Is_u2PHI)
//	{
//		int64_t nu1 = u1.Rows(),nu2 = u2.Rows();
//		u12.Redim(nu1+nu2,1);
//		int64_t i;
//		for(i=0; i<nu1; i++) u12(i,0) = u1(i,0);
//		for(i=0; i<nu2; i++) u12(i+nu1,0) = u2(i,0);
//
//
//	}
//	else if(!Is_u1PHI || !Is_u2PHI)
//	{
//		int64_t ru1 = u1.Rows(), cu1 = u1.Cols(), ru2 = u2.Rows(), cu2 = u2.Cols();
//		int64_t ru12 = ru1 < ru2 ? ru2 : ru1;
//		int64_t cu12 = cu1+cu2;
//		u12.Redim(ru12,cu12);
//		int64_t i,j;
//		for(i=0; i<ru1; i++) for(j=0; j<cu1; j++) u12(i,j) = u1(i,j);
//		for(i=0; i<ru2; i++) for(j=0; j<cu2; j++) u12(i,j+cu1) = u2(i,j);
//	}
//	else
//	{
//		PZError << "TPZCompElHCurl::Append. Bad input parameters " << std::endl;
//
//	}
//
//}
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::FillOrder(TPZVec<int> &order) const
//{
//    int nvecs = TSHAPE::Dimension*TSHAPE::NSides;
//    TPZManVector<int,3*27> associated_side(nvecs),bilinear(nvecs),direction(nvecs);
//    TSHAPE::GetSideHDivDirections(associated_side,direction,bilinear);
//    TPZManVector<int,27> sideinc(TSHAPE::NSides,0);
//    for (int iv=0; iv<nvecs; iv++) {
//        int side = associated_side[iv];
//        int bil = bilinear[iv];
//        if (bil) {
//            sideinc[side] = 1;
//        }
//    }
//    order.resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
//    int ncon = NConnects();
//
//    TPZConnect &c = this->Connect(ncon-1);
//    int internalorder = c.Order();
//    int nsides = TSHAPE::NSides;
//    for (int is=0; is<nsides; is++) {
//        if (TSHAPE::SideDimension(is) ==0) {
//            continue;
//        }
//        else if(TSHAPE::SideDimension(is) == TSHAPE::Dimension -1)
//        {
//            int intorder = internalorder;
//            if (sideinc[is]) {
//                intorder++;
//            }
//            int connectindex = SideConnectLocId(0, is);
//            if (connectindex < 0) {
//                DebugStop();
//            }
//            TPZConnect &c = this->Connect(connectindex);
//            if (c.Order() > intorder) {
//                intorder = c.Order();
//            }
//            order[is-TSHAPE::NCornerNodes] = intorder;
//        }
//        else
//        {
//            int intorder = internalorder;
//            if (sideinc[is]) {
//                intorder++;
//            }
//            order[is-TSHAPE::NCornerNodes] = intorder;
//        }
//    }
//}
//
//
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::NShapeContinuous(TPZVec<int> &order ){
//
//    return TSHAPE::NShapeF(order);
////#ifdef LOG4CXX
////		{
////				std::stringstream sout;
////				sout << "ordem max "<<maxorder<< " vec order " << order<<" num func cont "<< nshape<<std::endl;
////				LOGPZ_DEBUG(logger,sout.str())
////		}
////#endif
//
//
//}
//
//
//template<class TSHAPE>
//TPZTransform<> TPZCompElHCurl<TSHAPE>::TransformSideToElement(int side){
//	return TSHAPE::TransformSideToElement(side);
//}
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex){
//
//    DebugStop();
//	TPZManVector<int64_t> firstshapeindex;
//	FirstShapeIndex(firstshapeindex);
//	int nshape = TPZIntelGen<TSHAPE>::NShapeF();
//	shapeindex.Resize(nshape);
//	int64_t nsides = sides.NElements();
//	int64_t is, count=0;
//	for(is=0 ; is<nsides; is++)
//	{
//		int side = sides[is];
//        int conind = SideConnectLocId(0, side);
//		int sideorder = this->Connect(conind).Order();
//		int NShapeFace = TSHAPE::NConnectShapeF(side,sideorder);
//		int ishapeface;
//		for(ishapeface=0; ishapeface<NShapeFace; ishapeface++)
//		{
//			shapeindex[count++] = is;
//		}
//	}
//	shapeindex.Resize(count);
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled())
//	{
//		std::stringstream sout;
//		sout << "count = " << count << " nshape " << nshape;
//		sout << endl<<"sides associated with the normals "<< sides <<
//		"\nnormal associated with each shape function : shape function indexes " << shapeindex;
//		LOGPZ_DEBUG(logger,sout.str())
//	}
//#endif
//}
//

//
///** Initialize a material data and its attributes based on element dimension, number
// * of state variables and material definitions
// */
//
//
//
//// Save the element data to a stream
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::Write(TPZStream &buf, int withclassid) const
//{
//	TPZInterpolatedElement::Write(buf,withclassid);
//	TPZManVector<int,3> order(3,0);
//	this->fIntRule.GetOrder(order);
//	buf.Write(order);
//    buf.Write(fSideOrient);
//
//	buf.Write(this->fConnectIndexes.begin(),TSHAPE::NSides);
//	buf.Write(&this->fPreferredOrder,1);
//    buf.Write(fSideOrient);
//    int sz = fRestraints.size();
//    buf.Write(&sz);
//    for (std::list<TPZOneShapeRestraint>::const_iterator it = fRestraints.begin(); it != fRestraints.end(); it++) {
//        it->Write(buf);
//    }
//	int classid = this->ClassId();
//	buf.Write ( &classid, 1 );
//}
//
//
//// Read the element data from a stream
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::Read(TPZStream &buf, void *context)
//{
//	TPZInterpolatedElement::Read(buf,context);
//	TPZManVector<int,3> order;
//	buf.Read(order);
//	this-> fIntRule.SetOrder(order);
//    TPZManVector<int, TSHAPE::NFaces> SideOrient;
//    buf.Read(SideOrient);
//    fSideOrient = SideOrient;
//	buf.Read(this->fConnectIndexes.begin(),TSHAPE::NSides);
//	buf.Read(&this->fPreferredOrder,1);
//    buf.Read(fSideOrient);
//    int sz;
//    buf.Read(&sz);
//    for (int i=0; i<sz; i++) {
//        TPZOneShapeRestraint one;
//        one.Read(buf);
//        fRestraints.push_back(one);
//    }
//	int classid = -1;
//	buf.Read( &classid, 1 );
//	if ( classid != this->ClassId())
//	{
//		std::stringstream sout;
//		sout << "ERROR - " << __PRETTY_FUNCTION__
//        << " trying to restore an object id " << this->ClassId() << " and classid read = " << classid;
//		LOGPZ_ERROR ( logger, sout.str().c_str() );
//	}
//}
////refinamento
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::PRefine(int order)
//{
//    this->SetPreferredOrder(order);
//    int side;
//    int icon;
//    int ncon=NConnects();
//    TPZCompElHCurlPressure<TSHAPE> *hdivpressure = dynamic_cast<TPZCompElHCurlPressure<TSHAPE> *>(this);
//
//    if (hdivpressure) {
//        ncon--;
//    }
//    int nnodes = this->Reference()->NNodes();
//    for(icon=0; icon<ncon; icon++)
//    {//somente para os conects de fluxo
////        TPZConnect &con = this->Connect(icon);
////        con.SetOrder(order);
//        side= ConnectSideLocId(icon);
//
//#ifdef LOG4CXX
//        if (logger->isDebugEnabled())
//        {
//                std::stringstream sout;
//                sout << "side " << side << " order " << this->PreferredSideOrder(side)<<std::endl;
//                LOGPZ_DEBUG(logger,sout.str())
//        }
//#endif
//
//        this->IdentifySideOrder(side);
//    }
//    #ifdef LOG4CXX
//    if (loggerdiv->isDebugEnabled()) {
//        std::stringstream sout;
//        sout << (void*) this->Mesh() << "PRefine elindex " << this->Index() << " gel index " << this->Reference()->Index() << " " << order;
//        sout << "\nPRefine connect orders ";
//        int nc = this->NConnects();
//        for(int ic=0; ic<nc; ic++) sout << (int)this->Connect(ic).Order() << " ";
//        LOGPZ_DEBUG(loggerdiv, sout.str())
//    }
//#endif
//
//		// conect da pressao
//
//    if(ncon>nnodes+1)
//    {
//		TPZCompElHCurlPressure<TSHAPE> *hdivpressure = dynamic_cast<TPZCompElHCurlPressure<TSHAPE> *>(this);
//		TPZConnect &con = this->Connect(ncon-1);
//
//		if (TSHAPE::Type()==EQuadrilateral) {
//				hdivpressure->SetPressureOrder(order);
//				con.SetOrder(order,this->fConnectIndexes[ncon-1]);
//
//		}
//		else {
//				hdivpressure->SetPressureOrder(order-1);
//				con.SetOrder(order-1,this->fConnectIndexes[ncon-1]);
//
//		}
//		int nshape = hdivpressure-> NConnectShapeF(ncon-1,con.Order());
//		con.SetNShape(nshape);
//		int64_t seqnum = con.SequenceNumber();
//		this->Mesh()->Block().Set(seqnum,nshape);
//    }
//
//
//}
//
///** @brief Prints the relevant data of the element to the output stream */
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::Print(std::ostream &out) const
//{
//    out << __PRETTY_FUNCTION__ << std::endl;
//    TPZIntelGen<TSHAPE>::Print(out);
//    out << "Side orientation " << fSideOrient << std::endl;
//    if (fRestraints.size()) {
//        out << "One shape restraints associated with the element\n";
//        for (std::list<TPZOneShapeRestraint>::const_iterator it = fRestraints.begin(); it != fRestraints.end(); it++)
//        {
//            it->Print(out);
//        }
//    }
//
//
//
//}
//
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::MaxOrder(){
//
//    int maxorder = TPZInterpolationSpace::MaxOrder();
//    return maxorder+1;
//}
//
//#include "pzshapecube.h"
//#include "TPZRefCube.h"
//#include "pzshapelinear.h"
//#include "TPZRefLinear.h"
//#include "pzrefquad.h"
//#include "pzshapequad.h"
//#include "pzgeoquad.h"
//#include "pzshapetriang.h"
//#include "pzreftriangle.h"
//#include "pzgeotriangle.h"
//#include "pzshapeprism.h"
//#include "pzrefprism.h"
//#include "pzgeoprism.h"
//#include "pzshapetetra.h"
//#include "pzreftetrahedra.h"
//#include "pzgeotetrahedra.h"
//#include "pzshapepiram.h"
//#include "pzrefpyram.h"
//#include "pzgeopyramid.h"
//#include "pzrefpoint.h"
//#include "pzgeopoint.h"
//#include "pzshapepoint.h"
//#include "pzgraphelq2dd.h"
//#include "tpzgraphelt3d.h"
//#include "pzgraphel1dd.h"
//#include "pztrigraphd.h"
//#include "pzgraphelq3dd.h"
//#include "tpzgraphelprismmapped.h"
//#include "tpzgraphelpyramidmapped.h"
//#include "tpzgraphelt2dmapped.h"
//
//using namespace pztopology;
//
//#include "tpzpoint.h"
//#include "tpzline.h"
//#include "tpzquadrilateral.h"
//#include "tpztriangle.h"
//#include "tpzcube.h"
//#include "tpztetrahedron.h"
//#include "tpzprism.h"
//
//#include "pzelchdivbound2.h"
//
//using namespace pzgeom;
//using namespace pzshape;
//
//
//template<class TSHAPE>
//void TPZCompElHCurl<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
//	if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
//		new typename TSHAPE::GraphElType(this,&grafgrid);
//	}
//}
//
///// return the first one dof restraint
//template<class TSHAPE>
//int TPZCompElHCurl<TSHAPE>::RestrainedFace()
//{
//    return -1;
//}
//
//template<>
//int TPZCompElHCurl<TPZShapePiram>::RestrainedFace()
//{
//    if (fRestraints.size() == 0) {
//        return -1;
//        DebugStop(); //AQUIPHIL
//    }
//    std::list<TPZOneShapeRestraint>::iterator it = fRestraints.begin();
//    int foundis = -1;
//    bool found = false;
//    while (found == false && it != fRestraints.end()) {
//        int64_t connectindex = it->fFaces[3].first;
//        int64_t cindex = -1;
//        for (int is = 14; is<18; is++) {
//            cindex = ConnectIndex(is-13);
//            if (connectindex == cindex) {
//                found = true;
//                foundis = is+1;
//                if (foundis == 18) {
//                    foundis = 14;
//                }
//            }
//        }
//        it++;
//    }
//    if (found == false) {
//        DebugStop();
//    }
//    return foundis;
//}
//

//
//#include "pzreferredcompel.h"
//
//template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHCurl<TPZShapeLinear>>>;
//template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHCurl<TPZShapeTriang>>>;
//template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHCurl<TPZShapeQuad>>>;
//template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHCurl<TPZShapeCube>>>;
//template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHCurl<TPZShapeTetra>>>;
//template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHCurl<TPZShapePrism>>>;
//template class TPZRestoreClass< TPZReferredCompEl<TPZCompElHCurl<TPZShapePiram>>>;
//
//
//template class TPZReferredCompEl<TPZCompElHCurl<TPZShapeLinear>>;
//template class TPZReferredCompEl<TPZCompElHCurl<TPZShapeTriang>>;
//template class TPZReferredCompEl<TPZCompElHCurl<TPZShapeQuad>>;
//template class TPZReferredCompEl<TPZCompElHCurl<TPZShapeTetra>>;
//template class TPZReferredCompEl<TPZCompElHCurl<TPZShapePrism>>;
//template class TPZReferredCompEl<TPZCompElHCurl<TPZShapePiram>>;
//template class TPZReferredCompEl<TPZCompElHCurl<TPZShapeCube>>;
//
//TPZCompEl * CreateRefHCurlLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//    return new TPZReferredCompEl<TPZCompElHCurl< TPZShapeLinear>>(mesh,gel,index);
//}
//
//TPZCompEl * CreateRefHCurlQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//    return new TPZReferredCompEl<TPZCompElHCurl< TPZShapeQuad>>(mesh,gel,index);
//}
//
//TPZCompEl * CreateRefHCurlTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//    return new TPZReferredCompEl<TPZCompElHCurl< TPZShapeTriang >>(mesh,gel,index);
//}
//
//TPZCompEl * CreateRefHCurlCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//    return new TPZReferredCompEl<TPZCompElHCurl< TPZShapeCube >>(mesh,gel,index);
//}
//
//TPZCompEl * CreateRefHCurlPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//    return new TPZReferredCompEl<TPZCompElHCurl< TPZShapePrism>>(mesh,gel,index);
//}
//
//TPZCompEl * CreateRefHCurlPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//    return new TPZReferredCompEl<TPZCompElHCurl< TPZShapePiram >>(mesh,gel,index);
//}
//
//TPZCompEl * CreateRefHCurlTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//    return new TPZReferredCompEl<TPZCompElHCurl< TPZShapeTetra >>(mesh,gel,index);
//}

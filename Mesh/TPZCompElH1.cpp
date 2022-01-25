#include "TPZCompElH1.h"
#include "pzlog.h"
#include "pzconnect.h"
#include "TPZMaterial.h"
#include "TPZMatSingleSpace.h"
#include "TPZMaterialData.h"
#include "pzcmesh.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzintelgen");
#else
static int logger;
#endif


template<class TSHAPE>
TPZCompElH1<TSHAPE>::TPZCompElH1(TPZCompMesh &mesh, TPZGeoEl *gel, const H1Family h1fam) : TPZRegisterClassId(&TPZCompElH1::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel), fh1fam(h1fam){

	for(int i=0;i<TSHAPE::NSides;i++) {
		fConnectIndexes[i] = this->CreateMidSideConnect(i);
		mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
	}
}

template<class TSHAPE>
TPZCompElH1<TSHAPE>::TPZCompElH1(TPZCompMesh &mesh, TPZGeoEl *gel,
                                 int nocreate, const H1Family h1fam) :
TPZRegisterClassId(&TPZCompElH1::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel), fh1fam(h1fam)
{
	
}

template<class TSHAPE>
TPZCompElH1<TSHAPE>::TPZCompElH1(TPZCompMesh &mesh, const TPZCompElH1<TSHAPE> &copy) : TPZRegisterClassId(&TPZCompElH1::ClassId),
                                                                                       TPZIntelGen<TSHAPE>(mesh,copy), fh1fam(copy.fh1fam){
	this->fConnectIndexes = copy.fConnectIndexes;
}


template<class TSHAPE>
TPZCompElH1<TSHAPE>::TPZCompElH1(TPZCompMesh &mesh,
								 const TPZCompElH1<TSHAPE> &copy,
								 std::map<int64_t,int64_t> & gl2lcConMap,
								 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElH1::ClassId), 
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fh1fam(copy.fh1fam)
{
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		int64_t lcIdx = -1;
		int64_t glIdx = copy.fConnectIndexes[i];
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			fConnectIndexes[i] = -1;
			return;
		}
		fConnectIndexes[i] = lcIdx;
	}
  //TODO:NATHANFRAN
}

template<class TSHAPE>
TPZCompElH1<TSHAPE>::~TPZCompElH1() {
    TPZGeoEl *gel = this->Reference();
    if (gel) {
        TPZCompEl *cel = gel->Reference();
        if (cel == this) {
            this->RemoveSideRestraintsII(this->EDelete);
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
}

#include "pzgenericshape.h"
#include "TPZShapeData.h"
#include "TPZShapeH1.h"

template<class TSHAPE>
void TPZCompElH1<TSHAPE>::InitMaterialData(TPZMaterialData &data){
  data.gelElId = this->Reference()->Id();
    auto *mat =
        dynamic_cast<TPZMatSingleSpace*>(this->Material());
#ifdef PZDEBUG
    if(!mat)
    {
        DebugStop();
    }
#endif
    TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes);
    TPZManVector<int,TSHAPE::NSides> orders(TSHAPE::NSides-TSHAPE::NCornerNodes);
    TPZManVector<int,TSHAPE::NFacets> sideorient(TSHAPE::NFacets,0);
    TPZGeoEl *gel = this->Reference();
    for(int i=0; i<TSHAPE::NCornerNodes; i++) ids[i] = gel->NodeIndex(i);
    for(int i=TSHAPE::NCornerNodes; i<TSHAPE::NSides; i++) orders[i-TSHAPE::NCornerNodes] = this->Connect(i).Order();
    TPZShapeData &shapedata = data;
    TPZShapeH1<TSHAPE>::Initialize(ids, orders, shapedata);
    mat->FillDataRequirements(data);
    const int dim = this->Dimension();
    const int nshape = data.fPhi.Rows();
    const int nstate = this->Material()->NStateVariables();
    data.fShapeType = TPZMaterialData::EScalarShape;
    data.phi.Redim(nshape,1);
//    data.dphi.Redim(dim,nshape);
    data.dphix.Redim(dim,nshape);
    data.axes.Redim(dim,3);
    data.jacobian.Redim(dim,dim);
    data.jacinv.Redim(dim,dim);
    data.x.Resize(3);
    if (data.fNeedsSol){
        uint64_t ulen,durow,ducol;
        mat->GetSolDimensions(ulen,durow,ducol);
        data.SetSolSizes(nstate, ulen, durow, ducol);
    }
    if(data.fNeedsNeighborCenter)
    {
        TPZGeoElSide gelside(gel);
        gelside.CenterX(data.XCenter);
    }
    
}//void

template<class TSHAPE>
void TPZCompElH1<TSHAPE>::ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data){
    
    TPZShapeData &shapedata = data;
    TPZShapeH1<TSHAPE>::Shape(intpoint,shapedata);
    int tranpose = 1;
    REAL alpha = 1.;
    REAL beta = 0.;
    data.jacinv.MultAdd(shapedata.fDPhi, data.dphix, data.dphix,alpha,beta,tranpose);
    data.phi = shapedata.fPhi;
}

template<class TSHAPE>
void TPZCompElH1<TSHAPE>::SetConnectIndex(int i, int64_t connectindex){
#ifndef PZNODEBUG
	if(i<0 || i>= TSHAPE::NSides) {
		std::cout << " TPZIntelGen<TSHAPE>::SetConnectIndex index " << i <<
		" out of range\n";
		return;
	}
#endif
	fConnectIndexes[i] = connectindex;
}

/**sets the interpolation order of side to order*/
template<class TSHAPE>
void TPZCompElH1<TSHAPE>::SetSideOrder(int side, int order) {
	if(side<0 || side >= TSHAPE::NSides || (side >= TSHAPE::NCornerNodes && order <1)) {
		PZError << "TPZIntelGen::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
		DebugStop();
#ifdef PZ_LOG
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_ERROR(logger,sout.str())
#endif
		return;
	}
	if(side>= TSHAPE::NCornerNodes) {
		if(fConnectIndexes[side] == -1) return;
        int prevorder = this->fPreferredOrder;
        if (this->ConnectIndex(TSHAPE::NSides-1) != -1) {
            prevorder = EffectiveSideOrder(TSHAPE::NSides-1);
        }
		TPZConnect &c = this->Connect(side);
        int previousconnectorder = c.Order();
        if (order != previousconnectorder)
        {
            c.SetOrder(order,this->ConnectIndex(side));
            int64_t seqnum = c.SequenceNumber();
            int nvar = 1;
            TPZMaterial * mat = this->Material();
            if(mat) nvar = mat->NStateVariables();
            int nshape = TSHAPE::NConnectShapeF(side, order);
            c.SetNShape(nshape);
            c.SetNState(nvar);
            this->Mesh()->Block().Set(seqnum,nshape*nvar);
            TPZGeoElSide gelside(this->Reference(),side);
            TPZStack<TPZCompElSide> equal;
            gelside.EqualLevelCompElementList(equal, 1, 0);
            int64_t neq = equal.size();
            for (int64_t eq=0; eq<neq; eq++) {
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(equal[eq].Element());
                intel->AdjustIntegrationRule();
            }
        }
        
        int neworder = prevorder;
        if (this->ConnectIndex(TSHAPE::NSides-1) != -1) {
            neworder = EffectiveSideOrder(TSHAPE::NSides-1);
        }
		if(neworder != prevorder) {
            this->AdjustIntegrationRule();
		}
	}
}

/**returns the actual interpolation order of the polynomial along the side*/
template<class TSHAPE>
int TPZCompElH1<TSHAPE>::EffectiveSideOrder(int side) const {

	if(side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) return 0;
	if(fConnectIndexes[side] == -1)
	{
		std::stringstream sout ;
		sout << __PRETTY_FUNCTION__ << " side " << side << std::endl;
		//Print(sout);
#ifdef PZ_LOG
		LOGPZ_ERROR(logger,sout.str());
#else
		std::cout << sout.str() << std::endl;
#endif
		DebugStop();
		return -1;
	}
    TPZStack<int> lowdim;
    this->Reference()->LowerDimensionSides(side,lowdim);
	TPZConnect &c = this->Connect(side);
	int order = c.Order();
    for (int is=0; is<lowdim.size(); is++) {
        TPZConnect &c = this->MidSideConnect(lowdim[is]);
        if(c.Order() > order) order = c.Order();
    }
    return order;
}

template<class TSHAPE>
int TPZCompElH1<TSHAPE>::NConnectShapeF(int connect, int order) const{
    
    if(connect < TSHAPE::NCornerNodes) return TSHAPE::NConnectShapeF(connect,0);
	if(order < 0) return 0;
    int nshape = TSHAPE::NConnectShapeF(connect, order);
#ifdef PZDEBUG
    if(nshape < 0 )
    {
        nshape = TSHAPE::NConnectShapeF(connect, order);
        DebugStop();
    }
#endif
	return nshape;
}

template<class TSHAPE>
int TPZCompElH1<TSHAPE>::NSideConnects(int side) const {
	return TSHAPE::NContainedSides(side);
}

template<class TSHAPE>
int TPZCompElH1<TSHAPE>::SideConnectLocId(int node, int side) const {
	return TSHAPE::ContainedSideLocId(side,node);
}

/**Identifies the interpolation order on the interior of the element*/
template<class TSHAPE>
void TPZCompElH1<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
	ord.Resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
	int i;
	for(i=0; i<TSHAPE::NSides-TSHAPE::NCornerNodes; i++) {
		ord[i] = this->Connect(i+TSHAPE::NCornerNodes).Order();
	}
}


/**compute the values of the shape function of the side*/
template<class TSHAPE>
void TPZCompElH1<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
	int nc = TSHAPE::NContainedSides(side);
	int nn = TSHAPE::NSideNodes(side);
	TPZManVector<int64_t,27> id(nn);
	TPZManVector<int,27> order(nc-nn);
	int n,c;
	TPZGeoEl *ref = this->Reference();
	for (n=0;n<nn;n++){
		int nodloc = TSHAPE::SideNodeLocId(side,n);
		id [n] = ref->NodePtr(nodloc)->Id();
	}
	for (c=nn;c<nc;c++){
		int conloc = TSHAPE::ContainedSideLocId(side,c);
        order[c-nn] = this->Connect(conloc).Order();
	}
	TSHAPE::SideShape(side, point, id, order, phi, dphi);
}

template<class TSHAPE>
void TPZCompElH1<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
	TPZManVector<int64_t,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
	TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
	int i;
	TPZGeoEl *ref = this->Reference();
	for(i=0; i<TSHAPE::NCornerNodes; i++) {
		id[i] = ref->NodePtr(i)->Id();
	}
	for(i=0; i<TSHAPE::NSides-TSHAPE::NCornerNodes; i++) {
		ord[i] = this->Connect(i+TSHAPE::NCornerNodes).Order();
	}
	TSHAPE::Shape(pt,id,ord,phi,dphi);
}

#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"
#include "pzshapecube.h"
#include "pzshapepiram.h"
#include "pzshapepiramHdiv.h"



using namespace pzshape;

template class TPZCompElH1<TPZShapeTriang>;
template class TPZCompElH1<TPZShapePoint>;
template class TPZCompElH1<TPZShapeLinear>;
template class TPZCompElH1<TPZShapeQuad>;
template class TPZCompElH1<TPZShapeTetra>;
template class TPZCompElH1<TPZShapePrism>;
template class TPZCompElH1<TPZShapePiram>;
template class TPZCompElH1<TPZShapeCube>;

template class TPZRestoreClass< TPZCompElH1<TPZShapePoint>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeQuad>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeTetra>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapePrism>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapePiram>>;

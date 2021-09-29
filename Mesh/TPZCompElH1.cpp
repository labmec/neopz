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
TPZCompElH1<TSHAPE>::TPZCompElH1(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) : TPZRegisterClassId(&TPZCompElH1::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel,index){

	for(int i=0;i<TSHAPE::NSides;i++) {
		fConnectIndexes[i] = this->CreateMidSideConnect(i);
		mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
	}
}

template<class TSHAPE>
TPZCompElH1<TSHAPE>::TPZCompElH1(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index, int nocreate) : TPZRegisterClassId(&TPZCompElH1::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel,index)
{
	
}

template<class TSHAPE>
TPZCompElH1<TSHAPE>::TPZCompElH1(TPZCompMesh &mesh, const TPZCompElH1<TSHAPE> &copy) : TPZRegisterClassId(&TPZCompElH1::ClassId),
                                                                                       TPZIntelGen<TSHAPE>(mesh,copy){
	this->fConnectIndexes = copy.fConnectIndexes;
}


template<class TSHAPE>
TPZCompElH1<TSHAPE>::TPZCompElH1(TPZCompMesh &mesh,
								 const TPZCompElH1<TSHAPE> &copy,
								 std::map<int64_t,int64_t> & gl2lcConMap,
								 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElH1::ClassId), 
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap)
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
    mat->FillDataRequirements(data);
    const int dim = this->Dimension();
    const int nshape = this->NShapeF();
    const int nstate = this->Material()->NStateVariables();
    data.fShapeType = TPZMaterialData::EScalarShape;
    data.phi.Redim(nshape,1);
    data.dphi.Redim(dim,nshape);
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
    //Completing for three dimensional elements
    TPZManVector<REAL,3> x_center(3,0.0);
    TPZVec<REAL> center_qsi(dim,0.0);
    if (dim == 2) {
        if (this->Reference()->Type() == ETriangle) {
            center_qsi[0] = 0.25;
            center_qsi[1] = 0.25;
        }
    }
    else if (dim == 3) {
        if (this->Reference()->Type() == EPrisma) {
            center_qsi[0] = 1./3.;
            center_qsi[1] = 1./3.;
            center_qsi[2] = 0.0;
        }
        else if (this->Reference()->Type() == ETetraedro) {
            center_qsi[0] = 0.25;
            center_qsi[1] = 0.25;
            center_qsi[2] = 0.25;
        }
        else if (this->Reference()->Type() == EPiramide) {
            center_qsi[0] = 0.0;
            center_qsi[1] = 0.0;
            center_qsi[2] = 1./5.;
        }
    }
    this->Reference()->X(center_qsi, x_center);
    data.XCenter = x_center;
    
}//void


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
template class TPZCompElH1<TPZShapePiramHdiv>;
template class TPZCompElH1<TPZShapeCube>;

#ifndef BORLAND
template class TPZRestoreClass< TPZCompElH1<TPZShapePoint>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeQuad>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapeTetra>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapePrism>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapePiram>>;
template class TPZRestoreClass< TPZCompElH1<TPZShapePiramHdiv>>;
#endif
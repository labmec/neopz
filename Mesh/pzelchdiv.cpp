/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "pzcmesh.h"
#include "pzelchdiv.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "TPZMatSingleSpace.h"
#include "pzlog.h"
#include "pzgeoquad.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZMaterialDataT.h"
#include "pzshapepiram.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeH1.h"
#include "TPZShapeHDivConstant.h"
#include "TPZShapeHCurlNoGrads.h"

#include "pzshtmat.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDiv");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#endif

using namespace std;


template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDiv::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel,1), fSideOrient(TSHAPE::NFacets,1), fhdivfam(hdivfam) {
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	int nconflux= TPZCompElHDiv::NConnects();
    this->fConnectIndexes.Resize(nconflux);
	gel->SetReference(this);

//    int nfaces = TSHAPE::NumSides(TSHAPE::Dimension-1);
    TPZStack<int> facesides;
    TSHAPE::LowerDimensionSides(TSHAPE::NSides-1,facesides,TSHAPE::Dimension-1);
    facesides.Push(TSHAPE::NSides-1);
	for(int i=0;i< facesides.size(); i++)
	{
        int sideaux= facesides[i];
		this->fConnectIndexes[i] = this->CreateMidSideConnect(sideaux);
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "After creating last flux connect " << i << std::endl;
            //	this->Print(sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif

		mesh.ConnectVec()[this->fConnectIndexes[i]].IncrementElConnected();
		this->IdentifySideOrder(sideaux);
    }


    int sideorder = EffectiveSideOrder(TSHAPE::NSides-1);
//    if(TSHAPE::Type()==EQuadrilateral)
//    {
//        sideorder++;
//    }

    sideorder++;

	sideorder = 2*sideorder;
	if (sideorder > this->fIntRule.GetMaxOrder()) sideorder = this->fIntRule.GetMaxOrder();
	TPZManVector<int,3> order(3,sideorder);
	this->fIntRule.SetOrder(order);
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    for(int side = firstside ; side < TSHAPE::NSides-1; side++ )
    {
        fSideOrient[side-firstside] = this->Reference()->NormalOrientation(side);
    }
    auto *mat =
        dynamic_cast<TPZMatSingleSpace *>(this->Material());
    if (mat)
    {
        int order = mat->IntegrationRuleOrder(MaxOrder());
        TPZManVector<int,3> ord(gel->Dimension(),order);
        this->fIntRule.SetOrder(ord);
    }

    if (fhdivfam == HDivFamily::EHDivConstant) {
        // For HDiv constant, polynomial order was compatibilized in connectorders, 
        // see TPZShapeHDivConstant<TSHAPE>::Initialize. So now we need to update
        // the number of shape functions and also the integration rule
        if (TSHAPE::Type() == ETriangle || TSHAPE::Type() == EQuadrilateral || TSHAPE::Type() == ETetraedro){
            for (int icon = 0; icon < this->NConnects(); icon++)
            {
                TPZConnect &c = this->Connect(icon);
                int nShapeF = NConnectShapeF(icon,c.Order());
                if (c.NShape() != nShapeF){
                    DebugStop();
                }
            }
        }
    }
}

template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh, const TPZCompElHDiv<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHDiv::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy), fSideOrient(copy.fSideOrient), fhdivfam(copy.fhdivfam)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
    this->fConnectIndexes = copy.fConnectIndexes;
}

template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh,
									 const TPZCompElHDiv<TSHAPE> &copy,
									 std::map<int64_t,int64_t> & gl2lcConMap,
									 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHDiv::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fSideOrient(copy.fSideOrient), fhdivfam(copy.fhdivfam)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
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
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv() :
TPZRegisterClassId(&TPZCompElHDiv::ClassId),
TPZIntelGen<TSHAPE>()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
		this-> fConnectIndexes[i] = -1;
	}

}

template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::~TPZCompElHDiv(){
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
MElementType TPZCompElHDiv<TSHAPE>::Type() {
	return TSHAPE::Type();
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NConnects() const {
	return TSHAPE::NFacets + 1;
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetConnectIndex(int i, int64_t connectindex){
#ifndef PZNODEBUG
	if(i<0 || i>= this->NConnects()) {
		std::cout << " TPZCompElHDiv<TSHAPE>::SetConnectIndex index " << i <<
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
		sout << endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NConnectShapeF(int connect, int order)const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif
    switch (fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        return TPZShapeHDiv<TSHAPE>::ComputeNConnectShapeF(connect,order);    
        break;
    case HDivFamily::EHDivConstant:
        return TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(connect,order);
        break;
    
    default:
        return -1;
        break;
    }
    return -1;
 }

////
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetIntegrationRule(int ord) {
	TPZManVector<int,3> order(TSHAPE::Dimension,ord);
	this->fIntRule.SetOrder(order);
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NSideConnects(int side) const{
	if(TSHAPE::SideDimension(side)<= Dimension()-2) return 0;
	if(TSHAPE::SideDimension(side)==Dimension()-1) return 1;
	if(TSHAPE::SideDimension(side)== Dimension()) {
        int ncon = 1;
        return ncon;
    }
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << "Side: " << side <<"unhandled case ";
		LOGPZ_ERROR(logger,sout.str())
	}
#endif
	return -1;

}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::SideConnectLocId(int node,int side) const {
#ifdef PZDEBUG
	if(TSHAPE::SideDimension(side)<= TSHAPE::Dimension - 2 || node >= NSideConnects(side)) {
		PZError << "TPZCompElHDiv<TSHAPE>::SideConnectLocId no connect associate " <<  endl;
		return -1;
	}
#endif

    return node+side-(TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1);
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::ConnectSideLocId(int connect) const{

    int side = connect+TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1 ;
    return side;
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
	ord.Resize(NConnects());
	int i;
	for(i=0; i<NConnects(); i++) {
		ord[i] = ConnectOrder(i);
	}
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::PreferredSideOrder(int side) {
	if(TSHAPE::SideDimension(side) < Dimension()-1)
	{
		PZError << __PRETTY_FUNCTION__ << " side " << side << std::endl;
	}
	int connect= SideConnectLocId(0,side);
	if(connect<0 || connect > NConnects()) {
		PZError << "TPZCompElHDiv<TSHAPE>::PreferredSideOrder no polynomial associate " <<  endl;
		return -1;
	}
	if(connect<NConnects()) {
			int order =this->fPreferredOrder;
			return order;//this->AdjustPreferredSideOrder(side,order);
	}
	PZError << "TPZCompElHDiv<TSHAPE>::PreferredSideOrder called for connect = " << connect << "\n";
	return 0;

}

template<class TSHAPE>
int64_t TPZCompElHDiv<TSHAPE>::ConnectIndex(int con) const{
#ifndef PZNODEBUG
	if(con<0 || con > TSHAPE::NFacets) {
		std::cout << "TPZCompElHDiv::ConnectIndex wrong parameter connect " << con <<
		" NConnects " << TSHAPE::NFacets << std::endl;
		DebugStop();
		return -1;
	}

#endif

	return this->fConnectIndexes[con];
}


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetPreferredOrder(int order)
{
		TPZIntelGen<TSHAPE>:: SetPreferredOrder(order);
	//this->fPreferredOrder = order;
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetSideOrder(int side, int order){
	int connectaux= SideConnectLocId(0,side);
	if(connectaux<0 || connectaux > this-> NConnects()) {
		PZError << "TPZCompElHDiv::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef PZ_LOG
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
	}
	TPZConnect &c = this->Connect(connectaux);
    c.SetOrder(order,this->fConnectIndexes[connectaux]);
    int64_t seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial * mat =this-> Material();
    if(mat) nvar = mat->NStateVariables();
    c.SetNState(nvar);
    int nshape =this->NConnectShapeF(connectaux,order);
    c.SetNShape(nshape);
	this-> Mesh()->Block().Set(seqnum,nshape*nvar);
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::ConnectOrder(int connect) const{
	if (connect < 0 || connect >= this->NConnects()){
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "Connect index out of range connect " << connect <<
			" nconnects " << NConnects();
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return -1;
	}

	if (this->fConnectIndexes[connect] == -1) {
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " connect " << connect
		<< " is not initialized" << std::endl;
#ifdef PZ_LOG
		LOGPZ_ERROR(logger,sout.str());
#else
		std::cout << sout.str() << std::endl;
#endif
		return 0;
	}

    TPZConnect &c = this-> Connect(connect);
    return c.Order();
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::EffectiveSideOrder(int side) const
{
	if(!NSideConnects(side)) return -1;
	int cindex =SideConnectLocId(0, side);
#ifdef PZDEBUG
    if(cindex<0 || cindex >= NConnects()) DebugStop();
#endif
    return ConnectOrder(cindex);

}

/**
 * @brief It returns the normal orientation of the reference element by the side.
 * Only side that has dimension larger than zero and smaller than me.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::GetSideOrient(int side){

    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    return fSideOrient[side-firstside];
}

/**
 * @brief It set the normal orientation of the element by the side.
 * Only side that has dimension equal to my dimension minus one.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetSideOrient(int side, int sideorient){

    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    fSideOrient[side-firstside] = sideorient;
}

//compute the values of the shape function of the side
// this method is used by the method RestrainSide. It does not consider the side orientation
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {

    if(side==TSHAPE::NSides || point.size() != TSHAPE::Dimension-1){
        std::cout<<"Don't have side shape associated to this side";
        DebugStop();
    }
	if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension -1 ){
		return ;
	}
    int ncontained = TSHAPE::NContainedSides(side);
    int nsideshape = 0;
    int connectlocid = SideConnectLocId(0, side);
    int order = this->Connect(connectlocid).Order();
    int is;
    for (is=0; is<ncontained; is++) {
        int ic = TSHAPE::ContainedSideLocId(side,is);
        nsideshape += TSHAPE::NConnectShapeF(ic,order);
    }
#ifdef PZDEBUG
    if (nsideshape != this->NSideShapeF(side)) {
        DebugStop();
    }
#endif

    TPZGeoEl *gel = this->Reference();
    //int nc = gel->NCornerNodes();
    int nsn = TSHAPE::NSideNodes(side);
    TPZManVector<int64_t,8> id(nsn);
    for (int ic=0; ic<nsn; ic++) {
        int locid = TSHAPE::SideNodeLocId(side,ic);
        id[ic] = gel->Node(locid).Id();
    }

    //int idsize = id.size();
    TPZManVector<int,9> permutegather(ncontained);
    int transformid;


    MElementType sidetype = TSHAPE::Type(side);
    switch (sidetype) {
        case EOned:
            transformid = pztopology::TPZLine::GetTransformId(id);
            pztopology::TPZLine::GetSideHDivPermutation(transformid, permutegather);
            break;
        case EQuadrilateral:
            transformid = pztopology::TPZQuadrilateral::GetTransformId(id);
            pztopology::TPZQuadrilateral::GetSideHDivPermutation(transformid, permutegather);
            break;
        case ETriangle:
            transformid = pztopology::TPZTriangle::GetTransformId(id);
            pztopology::TPZTriangle::GetSideHDivPermutation(transformid, permutegather);
            break;
        default:
            DebugStop();
            break;
    }

    TPZManVector<int,TSHAPE::NSides> ord(TSHAPE::NSides,order);

    int sidedimension = TSHAPE::SideDimension(side);
    TPZFNMatrix<50,REAL> philoc(nsideshape,1),dphiloc(sidedimension,nsideshape);

    TSHAPE::SideShape(side,point,id,ord,philoc,dphiloc);

    int ncs = TSHAPE::NContainedSides(side);
    TPZManVector<int64_t,28> FirstIndex(ncs+1,0);
    for (int ls=0; ls<ncs; ls++) {
        int localside = TSHAPE::ContainedSideLocId(side,ls);
        FirstIndex[ls+1] = FirstIndex[ls]+TSHAPE::NConnectShapeF(localside,order);
    }

    REAL detjac = 1.;
    {
        TPZGeoElSide gelside = TPZGeoElSide(this->Reference(),side);
        int dim = gel->SideDimension(side);
        TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim),axes(dim,3);
        gelside.Jacobian(point, jac, axes, detjac, jacinv);
    }
    if(sidetype == ETriangle) detjac /= 6.;
    for (int side=0; side < ncs; side++) {
        int ifirst = FirstIndex[side];
        int kfirst = FirstIndex[permutegather[side]];
        int nshape = FirstIndex[side+1]-FirstIndex[side];
        for (int i=0; i<nshape; i++) {
            phi(ifirst+i,0) = philoc(kfirst+i,0)/detjac;
            for (int d=0; d< sidedimension; d++) {
                dphi(d,ifirst+i) = dphiloc(d,kfirst+i)/detjac;
            }
        }
    }

}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
    //this method is not really useful right now
    // TPZMaterialData data;

    // this->InitMaterialData(data);
    
    // TPZShapeHDiv<TSHAPE>::Shape(pt, data, data.phi, data.dphi);
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
{
    //TODOCOMPLEX
    if (var == 99) {
        return TPZIntelGen<TSHAPE>::Solution(qsi,var,sol);
    }
    TPZMaterialDataT<STATE> data;
    constexpr bool hasPhi{false};
    this->ComputeSolution(qsi,data,hasPhi);
    sol = std::move(data.sol[0]);
}

template<class TSHAPE>
template<class TVar>
void TPZCompElHDiv<TSHAPE>::ComputeSolutionHDivT(TPZMaterialDataT<TVar> &data)
{
    
    const int dim = 3; // Hdiv vectors are always in R3
    const int nstate = this->Material()->NStateVariables();
    const int ncon = this->NConnects();

    TPZFMatrix<TVar> &MeshSol = this->Mesh()->Solution();

    int64_t numbersol = MeshSol.Cols();

    if(numbersol != 1)
    {
        DebugStop();
    }
    data.sol.Resize(numbersol);
    data.dsol.Resize(numbersol);
    data.divsol.Resize(numbersol);

    for (int64_t is=0; is<numbersol; is++)
    {
        data.sol[is].Resize(dim*nstate);
        data.sol[is].Fill(0);
        data.dsol[is].Redim(dim*nstate, dim);
        data.divsol[is].Resize(nstate);
        data.divsol[is].Fill(0.);
    }


    TPZFMatrix<TVar> GradOfPhiHdiv(dim,dim);
    GradOfPhiHdiv.Zero();


    int normvecRows = data.fDeformedDirections.Rows();
    int normvecCols = data.fDeformedDirections.Cols();
    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
    TPZManVector<TPZFNMatrix<9,REAL>,18> GradNormalvec(normvecCols);
    for (int i=0; i<GradNormalvec.size(); i++) {
        GradNormalvec[i].Redim(dim,dim);
    }

    if (data.fNeedsDeformedDirectionsFad) {
        // Needs to be rethought
        DebugStop();
        for (int e = 0; e < normvecRows; e++) {
            for (int s = 0; s < normvecCols; s++) {
                Normalvec(e,s)=data.fDeformedDirectionsFad(e,s).val();
            }
        }

        TPZFNMatrix<4,REAL> Grad0(3,3,0.);
        TPZGeoEl *ref = this->Reference();
        const int gel_dim = ref->Dimension();

        for (int s = 0; s < normvecCols; s++) {
            for (int i = 0; i < gel_dim; i++) {
                for (int j = 0; j < gel_dim; j++) {
                    Grad0(i,j)=data.fDeformedDirectionsFad(i,s).fastAccessDx(j);
                }
            }
            GradNormalvec[s] = Grad0;
        }

    }else{
        Normalvec=data.fDeformedDirections;
    }

    TPZBlock &block =this->Mesh()->Block();
    int ishape=0,ivec=0,counter=0;

    int nshapeV = data.fVecShapeIndex.NElements();

    for(int in=0; in<ncon; in++)
    {
        TPZConnect *df = &this->Connect(in);
        int64_t dfseq = df->SequenceNumber();
        int dfvar = block.Size(dfseq);
        // pos : position of the block in the solution matrix
        int64_t pos = block.Position(dfseq);

        /// ish loops of the number of shape functions associated with the block
        for(int ish=0; ish<dfvar/nstate; ish++)
        {
            for (int64_t is=0; is<numbersol; is++)
            {
                for(int idf=0; idf<nstate; idf++)
                {
                    TVar meshsol = MeshSol(pos+ish*nstate+idf,is);
                    REAL phival = data.phi(ishape,0);
                    TPZManVector<REAL,3> normal(3);

                    for (int i=0; i<3; i++)
                    {
                        if (data.fNeedsDeformedDirectionsFad) {
                            normal[i] = data.fDeformedDirectionsFad(i,ivec).val();
                        }else{
                            normal[i] = data.fDeformedDirections(i,ivec);
                        }
                    }

#ifdef PZ_LOG
                    if(logger.isDebugEnabled() && abs(meshsol) > 1.e-6)
                    {
                        std::stringstream sout;
                        sout << "meshsol = " << meshsol << " ivec " << ivec << " ishape " << ishape << " x " << data.x << std::endl;
                        sout << " phi = " << data.phi(ishape,0)  << std::endl;
                        sout << "normal = " << normal << std::endl;
//                        sout << "GradOfPhiHdiv " << GradOfPhiHdiv << std::endl;
                        sout << "GradNormalVec " << GradNormalvec[ivec] << std::endl;
                        LOGPZ_DEBUG(logger,sout.str())
                    }
#endif

                    data.divsol[is][idf] += data.divphi(counter,0)*meshsol;
                    for (int ilinha=0; ilinha<dim; ilinha++) {
                        data.sol[is][ilinha+dim*idf] += normal[ilinha]*phival*meshsol;
                        for (int kdim = 0 ; kdim < dim; kdim++) {
                            data.dsol[is](ilinha+dim*idf,kdim)+= meshsol * GradOfPhiHdiv(ilinha,kdim);
                            if(data.fNeedsDeformedDirectionsFad){
                                data.dsol[is](ilinha+dim*idf,kdim)+=meshsol *GradNormalvec[ivec](ilinha,kdim)*data.phi(ishape,0);
                            }
                        }

                    }

                }
            }
            counter++;
            ivec++;
            ishape++;
        }
    }






#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "x " << data.x << " sol " << data.sol[0] << std::endl;
        data.dsol[0].Print("dsol",sout);
        sout << "divsol" << data.divsol[0] << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif

}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12)
{

	bool Is_u1PHI = (u1.Cols() == 1) ? true : false;
	bool Is_u2PHI = (u2.Cols() == 1) ? true : false;

	if(Is_u1PHI && Is_u2PHI)
	{
		int64_t nu1 = u1.Rows(),nu2 = u2.Rows();
		u12.Redim(nu1+nu2,1);
		int64_t i;
		for(i=0; i<nu1; i++) u12(i,0) = u1(i,0);
		for(i=0; i<nu2; i++) u12(i+nu1,0) = u2(i,0);


	}
	else if(!Is_u1PHI || !Is_u2PHI)
	{
		int64_t ru1 = u1.Rows(), cu1 = u1.Cols(), ru2 = u2.Rows(), cu2 = u2.Cols();
		int64_t ru12 = ru1 < ru2 ? ru2 : ru1;
		int64_t cu12 = cu1+cu2;
		u12.Redim(ru12,cu12);
		int64_t i,j;
		for(i=0; i<ru1; i++) for(j=0; j<cu1; j++) u12(i,j) = u1(i,j);
		for(i=0; i<ru2; i++) for(j=0; j<cu2; j++) u12(i,j+cu1) = u2(i,j);
	}
	else
	{
		PZError << "TPZCompElHDiv::Append. Bad input parameters " << std::endl;

	}

}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NShapeContinuous(TPZVec<int> &order ){

    return TSHAPE::NShapeF(order);
}


template<class TSHAPE>
TPZTransform<> TPZCompElHDiv<TSHAPE>::TransformSideToElement(int side){
	return TSHAPE::TransformSideToElement(side);
}

/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{

    TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes);
    TPZManVector<int,TSHAPE::NSides> orders(TSHAPE::NFacets+1,0);
    TPZManVector<int,TSHAPE::NFacets> sideorient(TSHAPE::NFacets,0);
    TPZGeoEl *gel = this->Reference();
    for(int i=0; i<TSHAPE::NCornerNodes; i++) ids[i] = gel->NodePtr(i)->Id();
    for(int i=0; i<TSHAPE::NFacets+1; i++) orders[i] = this->Connect(i).Order();
    for(int i=0; i<TSHAPE::NFacets; i++) sideorient[i] = this->SideOrient(i);
    TPZShapeData &shapedata = data;
    int nvec_shape = 0;

    switch (fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        TPZShapeHDiv<TSHAPE>::Initialize(ids, orders, sideorient, data);
        nvec_shape = TPZShapeHDiv<TSHAPE>::NShapeF(shapedata);
        break;
    case HDivFamily::EHDivConstant:
        TPZShapeHDivConstant<TSHAPE>::Initialize(ids, orders, sideorient, data);
        nvec_shape = this->NShapeF();
        break;
    
    default:
        break;
    }

//    int nshapescalar = shapedata.fPhi.Rows();
//    data.dphi.Resize(TSHAPE::Dimension, nshapescalar);
//    data.dphix.Resize(TSHAPE::Dimension, nshapescalar);
    // Trick to make actual hdiv materials work.
    // phi are all = 1. VecShapeIndex is 1 to 1 with its size the number of vec shapes
    
    data.phi.Resize(nvec_shape,1);
    data.fVecShapeIndex.Resize(nvec_shape);
    for (int ish = 0; ish<nvec_shape; ish++) {
        data.phi(ish,0) = 1.;
        data.fVecShapeIndex[ish] = {ish,ish};
    }
    
#ifdef PZ_LOG
        if(logger.isDebugEnabled())
		{
				LOGPZ_DEBUG(logger,"Initializing MaterialData of TPZCompElHDiv")
		}
#endif

    data.fShapeType = TPZMaterialData::EVecShape;

#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		data.fDeformedDirections.Print("Normal vector ", sout,EMathematicaInput);
        for (int i=0; i<TSHAPE::NCornerNodes; i++) {
            sout << "Id[" << i << "] = " << this->Reference()->NodePtr(i)->Id() << " ";
        }

        sout << std::endl;
		sout << "NormalVector/Shape indexes \n";
        for (int i=0; i<data.fVecShapeIndex.size(); i++) {
            sout << i << '|' << data.fVecShapeIndex[i] << " ";
        }
        sout << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif

}


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data) {

    constexpr int dim = TSHAPE::Dimension;
    const int nshape = this->NShapeF();
    TPZShapeData &shapedata = data;
    TPZFMatrix<REAL> phiMaster;

    switch (fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        TPZShapeHDiv<TSHAPE>::Shape(qsi, shapedata, phiMaster, data.divphi);
        break;

    case HDivFamily::EHDivConstant:
        phiMaster.Resize(TSHAPE::Dimension,nshape);
        data.divphi.Resize(nshape,1);
        TPZShapeHDivConstant<TSHAPE>::Shape(qsi, shapedata, phiMaster, data.divphi);
        break;

    default:
        DebugStop();
        break;
    }

    //This is for legacy material classes - will be deleted in the future
    int shapeSize = data.divphi.Rows();

    for (int i = 0; i < shapeSize; i++){
        data.fVecShapeIndex[i] = make_pair(i,i);
    }
    //
    
    TPZFNMatrix<9,REAL> gradx(3,TSHAPE::Dimension,0.);
    this->Reference()->GradX(qsi, gradx);
    gradx.MultAdd(phiMaster,data.fDeformedDirections,data.fDeformedDirections,1./fabs(data.detjac));
    data.divphi *= 1/fabs(data.detjac);
    data.phi = 1.;

}

// Save the element data to a stream
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	TPZInterpolatedElement::Write(buf,withclassid);
  buf.Write(fConnectIndexes.begin(),TSHAPE::NSides);
	TPZManVector<int,3> order(3,0);
	this->fIntRule.GetOrder(order);
	buf.Write(order);
    buf.Write(fSideOrient);

	buf.Write(this->fConnectIndexes.begin(),TSHAPE::NSides);
	buf.Write(&this->fPreferredOrder,1);
    buf.Write(fSideOrient);
    int sz = fRestraints.size();
    buf.Write(&sz);
    for (std::list<TPZOneShapeRestraint>::const_iterator it = fRestraints.begin(); it != fRestraints.end(); it++) {
        it->Write(buf);
    }
	int classid = this->ClassId();
	buf.Write ( &classid, 1 );
}


// Read the element data from a stream
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZInterpolatedElement::Read(buf,context);
  buf.Read(fConnectIndexes.begin(),TSHAPE::NSides);
	TPZManVector<int,3> order;
	buf.Read(order);
	this-> fIntRule.SetOrder(order);
    TPZManVector<int, TSHAPE::NFacets> SideOrient;
    buf.Read(SideOrient);
    fSideOrient = SideOrient;
	buf.Read(this->fConnectIndexes.begin(),TSHAPE::NSides);
	buf.Read(&this->fPreferredOrder,1);
    buf.Read(fSideOrient);
    int sz;
    buf.Read(&sz);
    for (int i=0; i<sz; i++) {
        TPZOneShapeRestraint one;
        one.Read(buf);
        fRestraints.push_back(one);
    }
	int classid = -1;
	buf.Read( &classid, 1 );
	if ( classid != this->ClassId())
	{
		std::stringstream sout;
		sout << "ERROR - " << __PRETTY_FUNCTION__
        << " trying to restore an object id " << this->ClassId() << " and classid read = " << classid;
		LOGPZ_ERROR ( logger, sout.str().c_str() );
	}
}
//refinamento
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::PRefine(int order)
{
    // This function has been deprecated since TPZCompElHDivPressure has been deprecated
    // TODO: Delete or update?
    /*
    this->SetPreferredOrder(order);
    int side;
    int icon;
    int ncon=NConnects();
    TPZCompElHDivPressure<TSHAPE> *hdivpressure = dynamic_cast<TPZCompElHDivPressure<TSHAPE> *>(this);

    if (hdivpressure) {
        ncon--;
    }
    int nnodes = this->Reference()->NNodes();
    for(icon=0; icon<ncon; icon++)
    {//somente para os conects de fluxo
//        TPZConnect &con = this->Connect(icon);
//        con.SetOrder(order);
        side= ConnectSideLocId(icon);

#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
                std::stringstream sout;
                sout << "side " << side << " order " << this->PreferredSideOrder(side)<<std::endl;
                LOGPZ_DEBUG(logger,sout.str())
        }
#endif

        this->IdentifySideOrder(side);
    }
    #ifdef PZ_LOG
    if (loggerdiv.isDebugEnabled()) {
        std::stringstream sout;
        sout << (void*) this->Mesh() << "PRefine elindex " << this->Index() << " gel index " << this->Reference()->Index() << " " << order;
        sout << "\nPRefine connect orders ";
        int nc = this->NConnects();
        for(int ic=0; ic<nc; ic++) sout << (int)this->Connect(ic).Order() << " ";
        LOGPZ_DEBUG(loggerdiv, sout.str())
    }
#endif

		// conect da pressao

    if(ncon>nnodes+1)
    {
		TPZCompElHDivPressure<TSHAPE> *hdivpressure = dynamic_cast<TPZCompElHDivPressure<TSHAPE> *>(this);
		TPZConnect &con = this->Connect(ncon-1);

		if (TSHAPE::Type()==EQuadrilateral) {
				hdivpressure->SetPressureOrder(order);
				con.SetOrder(order,this->fConnectIndexes[ncon-1]);

		}
		else {
				hdivpressure->SetPressureOrder(order-1);
				con.SetOrder(order-1,this->fConnectIndexes[ncon-1]);

		}
		int nshape = hdivpressure-> NConnectShapeF(ncon-1,con.Order());
		con.SetNShape(nshape);
		int64_t seqnum = con.SequenceNumber();
		this->Mesh()->Block().Set(seqnum,nshape);
    }

     */
}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZIntelGen<TSHAPE>::Print(out);
    out << "Side orientation " << fSideOrient << std::endl;
    if (fRestraints.size()) {
        out << "One shape restraints associated with the element\n";
        for (std::list<TPZOneShapeRestraint>::const_iterator it = fRestraints.begin(); it != fRestraints.end(); it++)
        {
            it->Print(out);
        }
    }



}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::MaxOrder(){

    int maxorder = TPZInterpolationSpace::MaxOrder();
    // if (fhdivfam == HDivFamily::EHDivConstant){
    //     if (TSHAPE::Type() == ETetraedro || TSHAPE::Type() == ETriangle){
    //         maxorder++;
    //     }
    // }
    return maxorder+1;
}

#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt2dmapped.h"

using namespace pztopology;

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}

/// return the first one dof restraint
template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::RestrainedFace()
{
    return -1;
}

// template<>
// int TPZCompElHDiv<TPZShapePiram>::RestrainedFace()
// {
//     if (fRestraints.size() == 0) {
//         return -1;
//         DebugStop(); //AQUIPHIL
//     }
//     std::list<TPZOneShapeRestraint>::iterator it = fRestraints.begin();
//     int foundis = -1;
//     bool found = false;
//     while (found == false && it != fRestraints.end()) {
//         int64_t connectindex = it->fFaces[3].first;
//         int64_t cindex = -1;
//         for (int is = 14; is<18; is++) {
//             cindex = ConnectIndex(is-13);
//             if (connectindex == cindex) {
//                 found = true;
//                 foundis = is+1;
//                 if (foundis == 18) {
//                     foundis = 14;
//                 }
//             }
//         }
//         it++;
//     }
//     if (found == false) {
//         DebugStop();
//     }
//     return foundis;
// }


template class TPZRestoreClass< TPZCompElHDiv<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDiv<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDiv<TPZShapeQuad>>;
template class TPZRestoreClass< TPZCompElHDiv<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElHDiv<TPZShapeTetra>>;
template class TPZRestoreClass< TPZCompElHDiv<TPZShapePrism>>;
// template class TPZRestoreClass< TPZCompElHDiv<TPZShapePiram>>;


template class TPZCompElHDiv<TPZShapeLinear>;
template class TPZCompElHDiv<TPZShapeTriang>;
template class TPZCompElHDiv<TPZShapeQuad>;
template class TPZCompElHDiv<TPZShapeTetra>;
template class TPZCompElHDiv<TPZShapePrism>;
// template class TPZCompElHDiv<TPZShapePiram>;
template class TPZCompElHDiv<TPZShapeCube>;


TPZCompEl * CreateHDivBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivBound2<TPZShapePoint>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivBound2< TPZShapeLinear>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivBound2< TPZShapeQuad>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivBound2< TPZShapeTriang >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDiv< TPZShapeLinear>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDiv< TPZShapeQuad>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDiv< TPZShapeTriang >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDiv< TPZShapeCube >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDiv< TPZShapePrism>(mesh,gel,hdivfam);
}

// TPZCompEl * CreateHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
// 	return new TPZCompElHDiv< TPZShapePiram >(mesh,gel,hdivfam);
// }

TPZCompEl * CreateHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDiv< TPZShapeTetra >(mesh,gel,hdivfam);
}


#include "TPZCompElHDivDuplConnects.h"
#include "TPZCompElHDivDuplConnectsBound.h"
#include "TPZMaterial.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeHDivConstant.h"
#include "pzlog.h"
#include "pzconnect.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDiv");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#endif

template<class TSHAPE>
TPZCompElHDivDuplConnects<TSHAPE>::TPZCompElHDivDuplConnects(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivDuplConnects::ClassId), TPZCompElHDiv<TSHAPE>(mesh,gel,hdivfam) {
    
    
}


template<class TSHAPE>
int TPZCompElHDivDuplConnects<TSHAPE>::NConnects() const {
	if(fDuplicationActive){
        return this->fConnectIndexes.size();
    } else {
        return TPZCompElHDiv<TSHAPE>::NConnects();
    }
}

template<class TSHAPE>
int TPZCompElHDivDuplConnects<TSHAPE>::NSideConnects(int side) const{
    
	if(TSHAPE::SideDimension(side)<= this->Dimension()-2) return 0;
	if(TSHAPE::SideDimension(side)== this->Dimension()-1) {
        if (fDuplicationActive){
            return 2;
        } else {
            return TPZCompElHDiv<TSHAPE>::NSideConnects(side);
        };
    }
	if(TSHAPE::SideDimension(side)== this->Dimension()) {
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
int TPZCompElHDivDuplConnects<TSHAPE>::NConnectShapeF(int connect, int order)const
{
    
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets*2) {
        DebugStop();
    }
#endif
    if (connect >= 2*TSHAPE::NFacets+1) return 0;

    switch (this->fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        {
            int conCorrect = connect/2;
            int res = connect % 2;
            int nshape;
            if (!fDuplicationActive){
                return TPZShapeHDiv<TSHAPE>::ComputeNConnectShapeF(connect,order);
            } else {
                nshape = TPZShapeHDiv<TSHAPE>::ComputeNConnectShapeF(conCorrect,order);
            }
            if (res == 1){ 
                nshape -= 1;
            } else {
                if (connect != 2*TSHAPE::NFacets){
                    nshape = 1;
                }
            }
            return nshape;   
        }
        break;
    case HDivFamily::EHDivConstant:
        {
            int conCorrect = connect/2;
            int res = connect % 2;
            int nshape;
            if (!fDuplicationActive){
                return TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(connect,order);
            } else {
                nshape = TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(conCorrect,order);
            }
             
            if (res == 1){ 
                nshape -= 1;
            } else {
                if (connect != 2*TSHAPE::NFacets){
                    nshape = 1;
                }
            }
            return nshape;
        }
        break;
    
    default:
        return -1;
        break;
    }
    return -1;
 }

template<class TSHAPE>
void TPZCompElHDivDuplConnects<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{

    TPZManVector<int64_t,TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes);
    TPZManVector<int,TSHAPE::NSides> orders(TSHAPE::NFacets+1,0);
    TPZManVector<int,TSHAPE::NFacets> sideorient(TSHAPE::NFacets,0);
    TPZGeoEl *gel = this->Reference();
    for(int i=0; i<TSHAPE::NCornerNodes; i++) ids[i] = gel->NodePtr(i)->Id();
    if (fDuplicationActive){
        for(int i=0; i<TSHAPE::NFacets; i++) orders[i] = this->Connect(2*i).Order();
        orders[TSHAPE::NFacets] = this->Connect(2*TSHAPE::NFacets).Order();
    } else {
        for(int i=0; i<TSHAPE::NFacets+1; i++) orders[i] = this->Connect(i).Order();
    }
    for(int i=0; i<TSHAPE::NFacets; i++) sideorient[i] = this->SideOrient(i);
    TPZShapeData &shapedata = data;
    int nvec_shape = 0;

    switch (this->fhdivfam)
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
int64_t TPZCompElHDivDuplConnects<TSHAPE>::ConnectIndex(int con) const{
    
#ifndef PZNODEBUG
	if(con<0 || con > TSHAPE::NFacets*2) {
		std::cout << "TPZCompElHDivDuplConnects::ConnectIndex wrong parameter connect " << con <<
		" NConnects " << TSHAPE::NFacets << std::endl;
		DebugStop();
		return -1;
	}

#endif

    if(fDuplicationActive){
        return this->fConnectIndexes[con];
    } else {
        return TPZCompElHDiv<TSHAPE>::ConnectIndex(con);
    }	
}


template<class TSHAPE>
int TPZCompElHDivDuplConnects<TSHAPE>::SideConnectLocId(int node,int side) const {
    if (fDuplicationActive){
        if (TSHAPE::Dimension == 2){
            return 2*(side-TSHAPE::NCornerNodes);
        } else if (TSHAPE::Dimension == 3){
            return 2*(side-(TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1));
        } else {
            DebugStop();
        }
    } else {
        return TPZCompElHDiv<TSHAPE>::SideConnectLocId(node,side);
    }
    return -1;
}


template<class TSHAPE>
void TPZCompElHDivDuplConnects<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
    if (fDuplicationActive){
        this->fConnectIndexes[i] = connectindex;
    }else{
        TPZCompElHDiv<TSHAPE>::SetConnectIndex(i,connectindex);
    }
}

template<class TSHAPE>
void TPZCompElHDivDuplConnects<TSHAPE>::ActiveDuplConnects(std::map<int64_t,int64_t> &fConnDuplicated){

    TPZManVector<int> conOrders(TSHAPE::NFacets,0);
    for (int i = 0; i < TSHAPE::NFacets; i++){
        conOrders[i] = this->ConnectOrder(i);
    }

    fDuplicationActive = true;

    this->fConnectIndexes.Resize(TSHAPE::NFacets*2+1);
    
    // std::cout << "Connects before = " << this->fConnectIndexes << std::endl;

    // Reorder the connects
    auto prevCon = this->fConnectIndexes;
    this->fConnectIndexes.Fill(-1);
    for (int i = 0; i < TSHAPE::NFacets; i++)
    {
        this->fConnectIndexes[2*i  ] = prevCon[i];    
        this->fConnectIndexes[2*i+1] = -1;//prevCon[i+TSHAPE::NFacets+1];

        
    }
    this->fConnectIndexes[TSHAPE::NFacets*2] = prevCon[TSHAPE::NFacets];
    // std::cout << "Connects after = " << this->fConnectIndexes << std::endl;

    auto nFacets = this->Reference()->NSides(this->Dimension()-1);
    //Loop over the element facets - which are the connects the be duplicated (edges in 2D and faces in 3D)
    for (int iFacet = 0; iFacet < nFacets; iFacet++)
    {
        
        // Algorithm description: for each element facet, checks if the corresponding original connect is in the map fConnDuplicated.
        // If Yes, just sets the returning value from fConnDuplicated to the duplicated connect in the current element;
        // If No, allocate a new connect and inserts its index to fConnDuplicated using the original connect as key

        int64_t conn = ConnectIndex(2*iFacet);           

        if (fConnDuplicated.find(conn) == fConnDuplicated.end()){
            //not found, so allocate a new connect
            auto pOrder = conOrders[iFacet];//this->Mesh()->GetDefaultOrder();
            int nshape = 0;//It is updated in the next loop
            int nstate = 1;//It can possibly change
            int64_t newConnect = this->Mesh()->AllocateNewConnect(nshape,nstate,pOrder);
            fConnDuplicated[conn] = newConnect;
            SetConnectIndex(2*iFacet+1,newConnect);
        } else {
            //found, so just set the proper index of the duplicated connect
            SetConnectIndex(2*iFacet+1,fConnDuplicated[conn]);
        }
    }

    //Updates the number of shape functions and also the integration rule. 
    //We need different casts because the element can be volumetric or boundary
    TPZInterpolatedElement *celHybrid = dynamic_cast<TPZInterpolatedElement *> (this); 
    if (celHybrid){
        int nConnects = celHybrid->NConnects();
        for (int icon = 0; icon < nConnects; icon++)
        {
            TPZConnect &c = celHybrid->Connect(icon);
            int nShapeF = celHybrid->NConnectShapeF(icon,c.Order());
            c.SetNShape(nShapeF);
            int64_t seqnum = c.SequenceNumber();
            int nvar = 1;
            TPZMaterial * mat = celHybrid->Material();
            if (mat) nvar = mat->NStateVariables();
            c.SetNState(nvar);
            celHybrid->Mesh()->Block().Set(seqnum, nvar * nShapeF);
        }
    }
}

template<class TSHAPE>
void TPZCompElHDivDuplConnects<TSHAPE>::InactiveDuplConnects(){

    auto auxCon = this->fConnectIndexes;
    
    for (int i = 0; i < TSHAPE::NFacets; i++)
    {
        auxCon[i] = this->fConnectIndexes[2*i  ];
        TPZConnect &c = this->Connect(2*i  );
        TPZConnect &c2 = this->Connect(2*i+1);
        int nshape = c2.NShape();
        nshape++;
        c.SetNShape(nshape);
        c2.SetNShape(0);
        // int64_t seqnum = c.SequenceNumber();
        // int nvar = 1;
        // TPZMaterial * mat = this->Material();
        // if (mat) nvar = mat->NStateVariables();
        // c.SetNState(nvar);
        // this->Mesh()->Block().Set(seqnum, nvar * nshape);
    }
    auxCon[TSHAPE::NFacets] = this->fConnectIndexes[TSHAPE::NFacets*2];

    // std::cout << "this->fConnectIndexes" << this->fConnectIndexes << std::endl;
    this->fConnectIndexes = auxCon;
    // std::cout << "this->fConnectIndexes" << this->fConnectIndexes << std::endl;
    this->fConnectIndexes.Resize(TSHAPE::NFacets+1);
    // std::cout << "this->fConnectIndexes" << this->fConnectIndexes << std::endl;
    fDuplicationActive = false;

    

}

#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
using namespace pzshape;

// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeLinear>>;
// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeTriang>>;
// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeQuad>>;
// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeCube>>;
// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeTetra>>;

template class TPZCompElHDivDuplConnects<TPZShapeLinear>;
template class TPZCompElHDivDuplConnects<TPZShapeTriang>;
template class TPZCompElHDivDuplConnects<TPZShapeQuad>;
template class TPZCompElHDivDuplConnects<TPZShapeTetra>;
template class TPZCompElHDivDuplConnects<TPZShapeCube>;

//BC elements
TPZCompEl * CreateHDivDuplConnectsBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivDuplConnectsBound< TPZShapePoint>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivDuplConnectsBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivDuplConnectsBound< TPZShapeLinear>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivDuplConnectsBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivDuplConnectsBound< TPZShapeQuad>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivDuplConnectsBoundTriangEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivDuplConnectsBound< TPZShapeTriang>(mesh,gel,hdivfam);
}



//Volumetric elements
TPZCompEl * CreateHDivDuplConnectsLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeLinear>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivDuplConnectsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeQuad>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivDuplConnectsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeTriang >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivDuplConnectsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeCube >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivDuplConnectsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeTetra >(mesh,gel,hdivfam);
}



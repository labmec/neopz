//
//  pzhdivpressurebound.cpp
//  PZ
//
//  Created by Agnaldo Farias on 25/06/13.
//
//

#include "pzhdivpressurebound.h"
#include "pzelchdivbound2.h"
#include "pzlog.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivPressureBound"));
#endif


template <class TSHAPE>
TPZCompElHDivPressureBound<TSHAPE>::TPZCompElHDivPressureBound(TPZCompMesh &mesh, TPZGeoEl *gel, int &index) :
TPZCompElHDivBound2<TSHAPE>(mesh, gel, index){
    
    
    //Creating connect of the pressure's variable
	this->fConnectIndexes.Resize(NConnects());
    this->fConnectIndexes[this->NConnects()-1] = this->CreateMidSideConnect(2);
    
    int nshape = 0;
    if(this->Dimension!= mesh.Dimension()-1) DebugStop();
    nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc::ETensorial);
    
    if(this->Dimension==2 && TSHAPE::Type()==EQuadrilateral){
        nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc::  ETensorial);
    }
    if (this->Dimension==2 && TSHAPE::Type()==ETriangle){
        nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc::  EOrdemTotal);
    }
    
    int nstate = 1;
    
    int newnodeindex = mesh.AllocateNewConnect(nshape,nstate,this->fPressureOrder);
    TPZConnect &newnod = mesh.ConnectVec()[newnodeindex];
    newnod.SetLagrangeMultiplier(1);
    this->fConnectIndexes[this->NConnects()-1] = this->CreateMidSideConnect(2);
    int seqnum = newnod.SequenceNumber();
    newnod.SetLagrangeMultiplier(1);
    mesh.Block().Set(seqnum,nshape);
    mesh.ConnectVec()[this->fConnectIndexes[this->NConnects()-1]].IncrementElConnected();
}

template<class TSHAPE>
int TPZCompElHDivPressureBound<TSHAPE>::NConnects() const {
	
	return TPZCompElHDivBound2<TSHAPE>::NConnects()+1;
}

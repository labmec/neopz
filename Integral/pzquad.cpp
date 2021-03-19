/**
 * @file
 * @brief Contains the implementation of the TPZInt1d methods. 
 */

#include "pzquad.h"
#include "tpzintrulet.h"
#include "tpzintrulet3d.h"
#include "tpzintrulep3d.h"
#include "tpzintrulelist.h"
#include "pzmanvector.h"
#include "pzerror.h"

#include <math.h>
#include <stdlib.h>

#ifdef VC
#include <algorithm>
#endif

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.integral.pzquad");
#endif

using namespace std;

/** TPZIntPoints method that returns a coherent maxime order to work */
/** It is necessary because now can to be computed rule with integration points up to one thousand. */
int TPZIntPoints::GetMaxOrder() const {
#ifdef VC
	return Max<int>(TPZIntRuleT::NRULESTRIANGLE_ORDER, Max<int>(TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER,TPZIntRuleP3D::NRULESPYRAMID_ORDER));
#else
	return fmaxl(TPZIntRuleT::NRULESTRIANGLE_ORDER, fmaxl(TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER,TPZIntRuleP3D::NRULESPYRAMID_ORDER));
#endif
}

void TPZIntPoints::Print(std::ostream &out) const {
	int np = NPoints();
	std::string namerule;
	Name(namerule);
	TPZVec<int> order(3,0);
    GetOrder(order);
	out << "Cubature rule (" << namerule << ") " << np << " : Order ( ";
	for(int i=0;i<Dimension();i++) out << order[i] << " ";
	out << ") \nNumber of points " << NPoints() << std::endl;
	int ip;
	TPZVec<REAL> pos(Dimension());
	REAL w;
	for(ip=0; ip<np; ip++) {
		Point(ip,pos,w);
		out << "ip " << ip << " pos " << pos << " w " << w << std::endl;
	}
}

//***** number of integration rules in PZINTVEC

int TPZInt1d::GetRealMaxOrder() const {
	if(fIntP->Type() == 1)
		return (TPZGaussRule::NRULESLOBATTO_ORDER);
	return TPZGaussRule::NRULESLEGENDRE_ORDER;
}

int TPZIntQuad::GetRealMaxOrder() const {
	if(fIntKsi->Type() == 1 || fIntEta->Type() == 1)
		return (TPZGaussRule::NRULESLOBATTO_ORDER - 1);
	return TPZGaussRule::NRULESLEGENDRE_ORDER - 1;
}

int TPZIntCube3D::GetRealMaxOrder() const {
	if(fIntKsi->Type() == 1 || fIntEta->Type() == 1 || fIntZeta->Type() == 1)
		return (TPZGaussRule::NRULESLOBATTO_ORDER - 2);
	return TPZGaussRule::NRULESLEGENDRE_ORDER - 2;
}

int TPZIntTriang::GetMaxOrder() const {	
	return TPZIntRuleT::NRULESTRIANGLE_ORDER;
}
int TPZIntTetra3D::GetMaxOrder() const {
	return TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER;
}
int TPZIntPyram3D::GetMaxOrder() const{
	return TPZIntRuleP3D::NRULESPYRAMID_ORDER;
}
int TPZIntPrism3D::GetMaxOrder() const {
	return TPZIntRuleT::NRULESTRIANGLE_ORDER;
}

//**************************************
TPZInt1d::TPZInt1d(int OrdX,int type) {
	fIntP = TPZIntRuleList::gIntRuleList.GetRule(OrdX,type);
	fOrdKsi 	= fIntP->Order();
}

int TPZInt1d::NPoints() const {
	if(fIntP) return fIntP->NInt();
	PZError << "Null Pointer passed to method TPZInt1d::TPZInt1d(TPZGaussRule *)\n";
	return 0;
}

void TPZInt1d::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
	if((fIntP) && ((ip >= 0) && (ip < NPoints()))){
		pos[0]	= fIntP->Loc(ip);
		w       = fIntP->W(ip);
		return;
	}
	if(!fIntP)
		PZError 	<< "Null Pointer passed to method " << "TPZInt1d::TPZInt1d(TPZGaussRule *)\n";
	if((ip < 0) || (ip >= NPoints()))
		PZError 	<< "ip = " << ip << ", Out of Range: 0 -> " << NPoints() << std::endl;
}

void TPZInt1d::SetOrder(TPZVec<int> &ord,int type){
	if(ord.NElements() < 1 || ord[0] < 0) {
		std::cout << "TPZINt1d::SetOrder: NULL number of integration points specified\n";
		return;
	}
	if(ord[0] > GetRealMaxOrder())
		ord[0] = GetRealMaxOrder();
	fOrdKsi = ord[0];
	fIntP   = TPZIntRuleList::gIntRuleList.GetRule(fOrdKsi,type);
    fOrdKsi = fIntP->Order();
}

void TPZInt1d::GetOrder(TPZVec<int> &ord) const{
	ord[0] =  fOrdKsi;
}

//**************************************
TPZIntQuad::TPZIntQuad(int OrdK, int OrdE){
	fIntKsi = TPZIntRuleList::gIntRuleList.GetRule(OrdK);
	fIntEta = TPZIntRuleList::gIntRuleList.GetRule(OrdE);
	fOrdKsi = fIntKsi->Order();
	fOrdEta = fIntEta->Order();
}

//**************************************
TPZIntQuad::TPZIntQuad(int OrdK){
    fIntKsi = TPZIntRuleList::gIntRuleList.GetRule(OrdK);
    fIntEta = TPZIntRuleList::gIntRuleList.GetRule(OrdK);
    fOrdKsi = fIntKsi->Order();
    fOrdEta = fIntEta->Order();
}

int TPZIntQuad::NPoints() const {
	if (!fIntKsi || !fIntEta){
		PZError << "Null Pointer passed to method TPZInt1d::TPZInt1d(TPZGaussRule *)\n";
		return 0;
	}
	return (fIntKsi->NInt() * fIntEta->NInt());
}

void TPZIntQuad::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
    
#ifdef PZDEBUG
    if(pos.size() == 1 || pos.size() > 3 )
    {
        DebugStop();
    }
#endif
	if((fIntEta) && (fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
		int ik, ie;
		ik = ip/fIntEta->NInt();
		ie = ip - (ip/fIntEta->NInt())*(fIntEta->NInt());
		pos[0] = fIntKsi->Loc(ik);
		pos[1] = fIntEta->Loc(ie);
		w      = fIntKsi->W(ik)*fIntEta->W(ie);
		return;
	}
	if(!fIntKsi || !fIntEta)
		PZError 	<< "Null Pointer passed to method " << "TPZInt1d::TPZInt1d(TPZGaussRule *)\n";
	if((ip < 0) || (ip >= NPoints()))
		PZError 	<< "ip = " << ip << ", Out of Range: 0 -> " << NPoints() << std::endl;
}

void TPZIntQuad::SetOrder(TPZVec<int> &ord,int type) {
    int prevtype = fIntKsi->Type();
    if (fOrdKsi != ord[0] || type != prevtype) {
        fOrdKsi = ord[0];
        fIntKsi = TPZIntRuleList::gIntRuleList.GetRule(fOrdKsi,type);
//        fOrdKsi = fIntKsi->Order();
    }
    prevtype = fIntEta->Type();
    if (fOrdEta != ord[1] || prevtype != type) {
        fOrdEta = ord[1];
        fIntEta = TPZIntRuleList::gIntRuleList.GetRule(fOrdEta,type);
//        fOrdEta = fIntEta->Order();
    }
}

void TPZIntQuad::GetOrder(TPZVec<int> &ord) const {
	ord[0] = fOrdKsi;
	ord[1] = fOrdEta;
}

//**************************************
TPZIntTriang::TPZIntTriang(int OrdK) {
	fIntKsi = TPZIntRuleList::gIntRuleList.GetRuleT(OrdK);
	fOrdKsi = fIntKsi->Order();
}

int TPZIntTriang::NPoints() const {
#ifdef PZDEBUG
	if (!fIntKsi){
		PZError << "Null Pointer passed to method TPZIntTriang::NPoints()\n";
        DebugStop();
		return 0;
	}
#endif
	return fIntKsi->NInt();
}

void TPZIntTriang::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
	if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
		fIntKsi->Loc(ip, pos);
		w = fIntKsi->W(ip);
		return;
	}
	if(!fIntKsi)
		PZError 	<< "Null Pointer passed to method " << "TPZIntTriang::Point(..)\n";
	if((ip < 0) || (ip >= NPoints()))
		PZError 	<< "ip = " << ip << ", Out of Range: 0 -> " << NPoints() << std::endl;
}

void TPZIntTriang::SetOrder(TPZVec<int> &ord,int type) {
#ifdef PZDEBUG
    if (type != 0) {
        DebugStop();
    }
#endif
	fOrdKsi = ord[0];
	if(ord[1] > ord[0]) fOrdKsi = ord[1];
	if(ord[0] < 0 || ord[0] > TPZIntRuleT::NRULESTRIANGLE_ORDER || ord[1] < 0 || ord[1] > TPZIntRuleT::NRULESTRIANGLE_ORDER) 
	{
#ifdef PZDEBUG
#ifdef PZ_LOG
		LOGPZ_WARN(logger,"Integration rule for triangle - Order is bigger than NRULESTRIANGLE_ORDER (Max)");
#endif
#endif
		
		fOrdKsi = TPZIntRuleT::NRULESTRIANGLE_ORDER;//havendo erro assume a maxima ordem
	}
	fIntKsi = TPZIntRuleList::gIntRuleList.GetRuleT(fOrdKsi);
    fOrdKsi = fIntKsi->Order();
}

void TPZIntTriang::GetOrder(TPZVec<int> &ord) const {
	ord[0] = fOrdKsi;
	ord[1] = ord[0];
}

//##############################################################################
TPZIntCube3D::TPZIntCube3D(int OrdK, int OrdE, int OrdZ) {
	fIntKsi  = TPZIntRuleList::gIntRuleList.GetRule(OrdK);
	fIntEta  = TPZIntRuleList::gIntRuleList.GetRule(OrdE);
	fIntZeta = TPZIntRuleList::gIntRuleList.GetRule(OrdZ);
	fOrdKsi  = fIntKsi->Order();
	fOrdEta  = fIntEta->Order();
	fOrdZeta = fIntZeta->Order();
}

int TPZIntCube3D::NPoints() const {
	if (!fIntKsi || !fIntEta|| !fIntZeta){
		PZError << "Null Pointer passed to method TPZIntCube3D::NPoints()\n";
		return 0;
	}
	return (fIntKsi->NInt() * fIntEta->NInt() * fIntZeta->NInt());
}

void TPZIntCube3D::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
	if((fIntZeta) && (fIntEta) && (fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
		int ik, ie , iz;
		ik = ip % fIntKsi->NInt();
		ie = (ip % (fIntKsi->NInt()*fIntEta->NInt()))/fIntKsi->NInt();
		iz = ip/(fIntKsi->NInt()*fIntEta->NInt());
		
		pos[0] 	= fIntKsi->Loc(ik);
		pos[1]	= fIntEta->Loc(ie);
		pos[2]	= fIntZeta->Loc(iz);
		w       = fIntKsi->W(ik)*fIntEta->W(ie)*fIntZeta->W(iz);
		return;
	}
	if(!fIntKsi || !fIntEta || !fIntZeta)
		PZError << "Null Pointer passed to method " << "TPZIntCube3D::Point(..)\n";
	if((ip < 0) || (ip >= NPoints()))
		PZError 	<< "ip = " << ip << ", Out of Range: 0 -> " << NPoints() << std::endl;
}

void TPZIntCube3D::SetOrder(TPZVec<int> &ord,int type) {
	fOrdKsi  = ord[0];
	fOrdEta  = ord[1];
	fOrdZeta = ord[2];
	fIntKsi  = TPZIntRuleList::gIntRuleList.GetRule(fOrdKsi,type);
	fIntEta  = TPZIntRuleList::gIntRuleList.GetRule(fOrdEta,type);
	fIntZeta = TPZIntRuleList::gIntRuleList.GetRule(fOrdZeta,type);
    fOrdKsi = fIntKsi->Order();
    fOrdEta = fIntEta->Order();
    fOrdZeta = fIntZeta->Order();
}

void TPZIntCube3D::GetOrder(TPZVec<int> &ord) const {
	ord[0] = fOrdKsi;
	ord[1] = fOrdEta;
	ord[2] = fOrdZeta;
}

//##############################################################################
TPZIntTetra3D::TPZIntTetra3D(int OrdK) {
	fIntKsi = TPZIntRuleList::gIntRuleList.GetRuleT3D(OrdK);
	fOrdKsi = fIntKsi->Order();
}

int TPZIntTetra3D::NPoints() const {
	if (!fIntKsi){
		PZError << "Null Pointer passed to method TPZIntTetra3D::NPoints()\n";
		return 0;
	}
	return fIntKsi->NInt();
}

void TPZIntTetra3D::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
	if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))) {
		fIntKsi->Loc(ip, pos);
		w = fIntKsi->W(ip);
		return;
	}
	if(!fIntKsi)
		PZError 	<< "Null Pointer passed to method " << "TPZIntTetra3D::Point(..)\n";
	if((ip < 0) || (ip >= NPoints()))
		PZError 	<< "ip = " << ip << ", Out of Range: 0 -> " << NPoints() << std::endl;
}

void TPZIntTetra3D::SetOrder(TPZVec<int> &ord,int type) {
	fOrdKsi = (ord[1] > ord[0]) ? ord[1] : ord[0];
	fOrdKsi = (fOrdKsi > ord[2]) ? fOrdKsi : ord[2];
	if(fOrdKsi > TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER) {
#ifdef PZDEBUG
#ifdef PZ_LOG
            if(logger.isWarnEnabled())
            {
                std::stringstream sout;
                sout << "Integration rule for tetrahedra - Order is bigger than NRULESTETRAHEDRA_ORDER (Max)";
                sout << " fOrdKsi " << fOrdKsi << " TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER " << TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER;
                LOGPZ_WARN(logger, sout.str())
            }
#endif
#endif
		fOrdKsi = TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER;
	}
	fIntKsi = TPZIntRuleList::gIntRuleList.GetRuleT3D(fOrdKsi);
    fOrdKsi = fIntKsi->Order();
}

void TPZIntTetra3D::GetOrder(TPZVec<int> &ord) const {
	ord[0] = fOrdKsi;
	ord[1] = ord[0];
	ord[2] = ord[0];
}

//##############################################################################
TPZIntPyram3D::TPZIntPyram3D(int OrdK) {
	fIntKsi = TPZIntRuleList::gIntRuleList.GetRuleP3D(OrdK);
	fOrdKsi = fIntKsi->Order();
}

int TPZIntPyram3D::NPoints() const {
	if(!fIntKsi) {
		PZError << "Null Pointer passed to method TPZIntPyram3D::NPoints()\n";
		return 0;
	}
	return fIntKsi->NInt();
}

void TPZIntPyram3D::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
	if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))) {
		fIntKsi->Loc(ip, pos);
		w = fIntKsi->W(ip);
		return;
	}
	if(!fIntKsi)
		PZError << "Null Pointer passed to method " << "TPZIntPyram3D::Point(..)\n";
	if((ip < 0) || (ip >= NPoints()))
		PZError << "ip = " << ip << ", Out of Range: 0 -> " << NPoints() << std::endl;
}

void TPZIntPyram3D::SetOrder(TPZVec<int> &ord,int type) {
	fOrdKsi = (ord[1] > ord[0]) ? ord[1] : ord[0];
	fOrdKsi = (fOrdKsi > ord[2]) ? fOrdKsi : ord[2];
	if(fOrdKsi > TPZIntRuleP3D::NRULESPYRAMID_ORDER) {
#ifdef PZDEBUG
#ifdef PZ_LOG
		LOGPZ_WARN(logger,"Integration rule for pyramid - Order is bigger than NRULESPYRAMID_ORDER (Max)");
#endif
#endif
		fOrdKsi = TPZIntRuleP3D::NRULESPYRAMID_ORDER;
	}
	fIntKsi = TPZIntRuleList::gIntRuleList.GetRuleP3D(fOrdKsi);
    fOrdKsi = fIntKsi->Order();
}

void TPZIntPyram3D::GetOrder(TPZVec<int> &ord) const {
	ord[0] = fOrdKsi;
	ord[1] = ord[0];
	ord[2] = ord[0];
}

//##############################################################################
TPZIntPrism3D::TPZIntPrism3D(int OrdK,int OrdL) : fIntRule1D(OrdK), fIntTriang(OrdL) {
    TPZManVector<int,2> ord1d(1),ordt(2);
    fIntRule1D.GetOrder(ord1d);
    fIntTriang.GetOrder(ordt);
	fOrdKsi = ord1d[0];
	fOrdKti = ordt[0];
}

TPZIntPrism3D::~TPZIntPrism3D() {
}

int TPZIntPrism3D::NPoints() const {
	return (fIntRule1D.NPoints()*fIntTriang.NPoints());
}

void TPZIntPrism3D::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
	if((ip >= 0) && (ip < NPoints())) {		
		REAL v;
		TPZVec<REAL> ps(2);
		fIntTriang.Point(ip % fIntTriang.NPoints(), ps, v);
		pos[0] = ps[0]; pos[1] = ps[1];
		fIntRule1D.Point(ip / fIntTriang.NPoints(), ps, w);
		pos[2] = ps[0];
		w *= v;
		return;
	}
	if((ip < 0) || (ip >= NPoints()))
		PZError 	<< "ip = " << ip << ", Out of Range: 0 -> " << NPoints() << std::endl;
}

void TPZIntPrism3D::SetOrder(TPZVec<int> &ord,int type) {	
	fOrdKsi = ord[0];   //ordem na reta : zeta
	fOrdKti = (ord[1] > ord[2]) ? ord[1] : ord[2];   //ordem no plano XY
	TPZManVector<int,2> prc1(1),prc2(2);
	prc1[0] = ord[0];
	prc2[0] = ord[1];
	prc2[1] = ord[2];
	fIntRule1D.SetOrder(prc1,type);
	fIntTriang.SetOrder(prc2);
    fIntRule1D.GetOrder(prc1);
    fIntTriang.GetOrder(prc2);
    fOrdKsi = prc1[0];
    fOrdKti = prc2[0];
}

void TPZIntPrism3D::GetOrder(TPZVec<int> &ord) const {
	ord[0] = fOrdKsi;
	ord[1] = fOrdKti;
	ord[2] = ord[1];
}

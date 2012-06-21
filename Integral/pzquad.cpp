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

using namespace std;

/** TPZIntPoints method that returns a coherent maxime order to work */
/** It is necessary because now can to be computed rule with integration points up to one thousand. */
int TPZIntPoints::GetMaxOrder() const {
#ifdef VC
	return min(TPZIntRuleT::NRULESTRIANGLE_ORDER, \
			   min(TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER,TPZIntRuleP3D::NRULESPYRAMID_ORDER));
#else
	return fminl(TPZIntRuleT::NRULESTRIANGLE_ORDER, \
			fminl(TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER,TPZIntRuleP3D::NRULESPYRAMID_ORDER));
#endif
}

void TPZIntPoints::Print(std::ostream &out) {
	int np = Dimension();
	std::string namerule;
	Name(namerule);
	TPZVec<int> order(3,0);
	out << "Cubature rule (" << namerule << ") " << np << "d : Order ( ";
	for(int i=0;i<np;i++) out << order[i] << " ";
	out << ") \t Number of points " << NPoints() << std::endl;
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
	fIntP = gIntRuleList.GetRule(OrdX,type);
	fOrdKsi 	= OrdX;
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
	fIntP   = gIntRuleList.GetRule(fOrdKsi,type);
}

void TPZInt1d::GetOrder(TPZVec<int> &ord) const{
	ord[0] =  fOrdKsi;
}

//**************************************
TPZIntQuad::TPZIntQuad(int OrdK, int OrdE){
	fIntKsi = gIntRuleList.GetRule(OrdK);
	fIntEta = gIntRuleList.GetRule(OrdE);
	fOrdKsi = OrdK;
	fOrdEta = OrdE;
}

int TPZIntQuad::NPoints() const {
	if (!fIntKsi || !fIntEta){
		PZError << "Null Pointer passed to method TPZInt1d::TPZInt1d(TPZGaussRule *)\n";
		return 0;
	}
	return (fIntKsi->NInt() * fIntEta->NInt());
}

void TPZIntQuad::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
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
	fOrdKsi = ord[0];
	fOrdEta = ord[1];
	fIntKsi = gIntRuleList.GetRule(fOrdKsi,type);
	fIntEta = gIntRuleList.GetRule(fOrdEta,type);
}

void TPZIntQuad::GetOrder(TPZVec<int> &ord) const {
	ord[0] = (short) fOrdKsi;
	ord[1] = (short) fOrdEta;
}

//**************************************
TPZIntTriang::TPZIntTriang(int OrdK) {
	fIntKsi = gIntRuleList.GetRuleT(OrdK);
	fOrdKsi = OrdK;
}

int TPZIntTriang::NPoints() const {
	if (!fIntKsi){
		PZError << "Null Pointer passed to method TPZIntTriang::NPoints()\n";
		return 0;
	}
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
	fOrdKsi = ord[0];
	if(ord[1] > ord[0]) fOrdKsi = ord[1];
	if(ord[0] < 0 || ord[0] > TPZIntRuleT::NRULESTRIANGLE_ORDER || ord[1] < 0 || ord[1] > TPZIntRuleT::NRULESTRIANGLE_ORDER) 
		fOrdKsi = TPZIntRuleT::NRULESTRIANGLE_ORDER;//havendo erro assume a maxima ordem
	fIntKsi = gIntRuleList.GetRuleT(fOrdKsi);
}

void TPZIntTriang::GetOrder(TPZVec<int> &ord) const {
	ord[0] = (short) fOrdKsi;
	ord[1] = ord[0];
}

//##############################################################################
TPZIntCube3D::TPZIntCube3D(int OrdK, int OrdE, int OrdZ) {
	fIntKsi  = gIntRuleList.GetRule(OrdK);
	fIntEta  = gIntRuleList.GetRule(OrdE);
	fIntZeta = gIntRuleList.GetRule(OrdZ);
	fOrdKsi  = OrdK;
	fOrdEta  = OrdE;
	fOrdZeta = OrdZ;
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
	fIntKsi  = gIntRuleList.GetRule(fOrdKsi,type);
	fIntEta  = gIntRuleList.GetRule(fOrdEta,type);
	fIntZeta = gIntRuleList.GetRule(fOrdZeta,type);
}

void TPZIntCube3D::GetOrder(TPZVec<int> &ord) const {
	ord[0] = (short) fOrdKsi;
	ord[1] = (short) fOrdEta;
	ord[2] = (short) fOrdZeta;
}

//##############################################################################
TPZIntTetra3D::TPZIntTetra3D(int OrdK) {
	fIntKsi = gIntRuleList.GetRuleT3D(OrdK);
	fOrdKsi = OrdK;
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
	if(fOrdKsi > TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER)
		fOrdKsi = TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER;
	fIntKsi = gIntRuleList.GetRuleT3D(fOrdKsi);
}

void TPZIntTetra3D::GetOrder(TPZVec<int> &ord) const {
	ord[0] = (short) fOrdKsi;
	ord[1] = ord[0];
	ord[2] = ord[0];
}

//##############################################################################
TPZIntPyram3D::TPZIntPyram3D(int OrdK) {
	fIntKsi = gIntRuleList.GetRuleP3D(OrdK);
	fOrdKsi = OrdK;
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
	if(fOrdKsi > TPZIntRuleP3D::NRULESPYRAMID_ORDER)
		fOrdKsi = TPZIntRuleP3D::NRULESPYRAMID_ORDER;
	fIntKsi = gIntRuleList.GetRuleP3D(fOrdKsi);
}

void TPZIntPyram3D::GetOrder(TPZVec<int> &ord) const {
	ord[0] = (short) fOrdKsi;
	ord[1] = ord[0];
	ord[2] = ord[0];
}

//##############################################################################
TPZIntPrism3D::TPZIntPrism3D(int OrdK,int OrdL) : fIntRule1D(OrdK), fIntTriang(OrdL) {
	fOrdKsi = OrdK;
	fOrdKti = OrdL;
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
	TPZVec<int> prc1(1),prc2(2);
	prc1[0] = ord[0];
	prc2[0] = ord[1];
	prc2[1] = ord[2];
	fIntRule1D.SetOrder(prc1,type);
	fIntTriang.SetOrder(prc2);
}

void TPZIntPrism3D::GetOrder(TPZVec<int> &ord) const {
	ord[0] = fOrdKsi;
	ord[1] = fOrdKti;
	ord[2] = ord[1];
}

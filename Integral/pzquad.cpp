
#include "pzquad.h"
#include "tpzintrule.h"
#include "tpzintrulet.h"
#include "tpzintrulet3d.h"
#include "tpzintrulep3d.h"
#include "tpzintrulelist.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include <math.h>
#include <stdlib.h>

//#pragma segment UTIL
//***** number of integration rules in PZINTVEC



void TPZIntPyram3D::GetOrder(TPZVec<int> &ord) const {
  ord[0] = (short) fOrdKsi;
  ord[1] = ord[0];
  ord[2] = ord[0];
}


int TPZInt1d::GetMaxOrder() const{return TPZIntRule::NUMINT_RULES - 1;}
int TPZIntTriang::GetMaxOrder() const{return TPZIntRuleT::NUMINT_RULEST - 1 ;}
int TPZIntQuad::GetMaxOrder() const{return TPZIntRule::NUMINT_RULES - 2 ;}
int TPZIntCube3D::GetMaxOrder() const {return TPZIntRule::NUMINT_RULES - 3;}
int TPZIntTetra3D::GetMaxOrder() const {return TPZIntRuleT3D::NUMINT_RULEST3D - 1;}
int TPZIntPyram3D::GetMaxOrder() const{return TPZIntRuleT3D::NUMINT_RULEST3D - 1;}
int TPZIntPrism3D::GetMaxOrder() const {return TPZIntRuleP3D::NUMINT_RULESP3D - 1;}





//**************************************
//**************************************
TPZInt1d::TPZInt1d(int OrdX){
  fOrdKsi 	= OrdX;
  fIntP 	= gIntRuleList.GetRule(OrdX);
}

void TPZInt1d::SetOrder(TPZVec<int> &ord){
  if(ord.NElements() < 1) {
    std::cout << "TPZINt1d::SetOrder: NULL number of integration points specified\n";
    return;
  }
  fOrdKsi = ord[0];
  fIntP   = gIntRuleList.GetRule(ord[0]);
}

void TPZInt1d::GetOrder(TPZVec<int> &ord) const{
  ord[0] =  fOrdKsi;
}

int TPZInt1d::NPoints() const {
  if(fIntP) return fIntP->NInt();
  PZError << "Null Pointer passed to method TPZInt1d::TPZInt1d(TPZIntRule *)\n";
  return 0;
}

void TPZInt1d::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
  if((fIntP) && ((ip >= 0) && (ip < NPoints()))){
    pos[0] 	= fIntP->Loc(ip);
    w        = fIntP->W(ip);
    return;
  }
  if(!fIntP)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZInt1d::TPZInt1d(TPZIntRule *)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << std::endl;
}

//**************************************
//**************************************
TPZIntQuad::TPZIntQuad(int OrdK, int OrdE){
  fOrdKsi = OrdK;
  fOrdEta = OrdE;
  fIntKsi = gIntRuleList.GetRule(OrdK);
  fIntEta = gIntRuleList.GetRule(OrdE);
}


int TPZIntQuad::NPoints() const {
  if (!fIntKsi || !fIntEta){
    PZError << "Null Pointer passed to method TPZInt1d::TPZInt1d(TPZIntRule *)\n";
    return 0;
  }
  return fIntKsi->NInt() * fIntEta->NInt();
}

void TPZIntQuad::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
  if((fIntEta) && (fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    int ik, ie;
    ik = ip/fIntEta->NInt();
    ie = ip - (ip/fIntEta->NInt())*(fIntEta->NInt());
    pos[0] 	= fIntKsi->Loc(ik);
    pos[1]	= fIntEta->Loc(ie);
    w        = fIntKsi->W(ik)*fIntEta->W(ie);
    return;
  }
  if(!fIntKsi || !fIntEta)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZInt1d::TPZInt1d(TPZIntRule *)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << std::endl;
}

void TPZIntQuad::SetOrder(TPZVec<int> &ord){
  fOrdKsi = ord[0];
  fOrdEta = ord[1];
  fIntKsi = gIntRuleList.GetRule(ord[0]);
  fIntEta = gIntRuleList.GetRule(ord[1]);
}

void TPZIntQuad::GetOrder(TPZVec<int> &ord) const {
  ord[0] = (short) fOrdKsi;
  ord[1] = (short) fOrdEta;
}

//**************************************
//**************************************
TPZIntTriang::TPZIntTriang(	int OrdK){
  fOrdKsi = OrdK;
  fIntKsi = gIntRuleList.GetRuleT(OrdK);
}

int  TPZIntTriang::NPoints() const{
  if (!fIntKsi){
    PZError << "Null Pointer passed to method TPZIntTriang::NPoints()\n";
    return 0;
  }
  return fIntKsi->NInt();
}

void TPZIntTriang::Point(int ip, TPZVec<REAL> &pos, REAL &w) const{
  if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    fIntKsi->Loc(ip, pos);
    w = fIntKsi->W(ip);
    return;
  }
  if(!fIntKsi)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntTriang::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << std::endl;
}

void TPZIntTriang::SetOrder(TPZVec<int> &ord){
  fOrdKsi = ord[0];
  if(ord[1]>ord[0]) fOrdKsi = ord[1];
  fIntKsi = gIntRuleList.GetRuleT(fOrdKsi);
  if(ord[0] <0 || ord[0]>11 || ord[1] <0 || ord[1]>11) fOrdKsi = 10;//havendo erro assume a maxima ordem
}

void TPZIntTriang::GetOrder(TPZVec<int> &ord) const {
  ord[0] = (short) fOrdKsi;
  ord[1] = ord[0];
}
//##############################################################################
//##############################################################################
TPZIntCube3D::TPZIntCube3D(int OrdK, int OrdE, int OrdZ){
  fOrdKsi = OrdK;
  fOrdEta = OrdE;
  fOrdZeta = OrdZ;
  fIntKsi = gIntRuleList.GetRule(OrdK);
  fIntEta = gIntRuleList.GetRule(OrdE);
  fIntZeta = gIntRuleList.GetRule(OrdZ);
}
//------------------------------------------------------------------------------
void TPZIntCube3D::SetOrder(TPZVec<int> &ord){
  fOrdKsi = ord[0];
  fOrdEta = ord[1];
  fOrdZeta = ord[2];
  fIntKsi = gIntRuleList.GetRule(ord[0]);
  fIntEta = gIntRuleList.GetRule(ord[1]);
  fIntZeta = gIntRuleList.GetRule(ord[2]);
}
//------------------------------------------------------------------------------
int TPZIntCube3D::NPoints() const {
  if (!fIntKsi || !fIntEta|| !fIntZeta){
    PZError << "Null Pointer passed to method TPZIntCube3D::NPoints()\n";
    return 0;
  }
  return fIntKsi->NInt() * fIntEta->NInt() * fIntZeta->NInt();
}
//------------------------------------------------------------------------------
void TPZIntCube3D::Point(int ip, TPZVec<REAL> &pos, REAL &w) const{

  if((fIntZeta) && (fIntEta) && (fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    int ik, ie , iz;
    
    ik = ip % fIntKsi->NInt();
    ie = (ip % (fIntKsi->NInt()*fIntEta->NInt()))/fIntKsi->NInt();
    iz = ip/(fIntKsi->NInt()*fIntEta->NInt());
    
    pos[0] 	= fIntKsi->Loc(ik);
    pos[1]	= fIntEta->Loc(ie);
    pos[2]	= fIntZeta->Loc(iz);
    w        = fIntKsi->W(ik)*fIntEta->W(ie)*fIntZeta->W(iz);
    return;
  }
  if(!fIntKsi || !fIntEta || !fIntZeta)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntCube3D::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << std::endl;
}
//------------------------------------------------------------------------------
void TPZIntCube3D::GetOrder(TPZVec<int> &ord) const {
  ord[0] = (short) fOrdKsi;
  ord[1] = (short) fOrdEta;
  ord[2] = (short) fOrdZeta;
}
//##############################################################################
//##############################################################################
TPZIntTetra3D::TPZIntTetra3D(int OrdK){
  fOrdKsi = OrdK;
  fIntKsi = gIntRuleList.GetRuleT3D(OrdK);
}
//------------------------------------------------------------------------------
void TPZIntTetra3D::SetOrder(TPZVec<int> &ord){

  fOrdKsi = (ord[1] > ord[0]) ? ord[1] : ord[0];
  fOrdKsi = (fOrdKsi > ord[2]) ? fOrdKsi : ord[2];
  fIntKsi = gIntRuleList.GetRuleT3D(fOrdKsi);
}
//------------------------------------------------------------------------------
int TPZIntTetra3D::NPoints() const {
  if (!fIntKsi){
    PZError << "Null Pointer passed to method TPZIntTetra3D::NPoints()\n";
    return 0;
  }
  return fIntKsi->NInt();
}
//------------------------------------------------------------------------------
void TPZIntTetra3D::Point(int ip, TPZVec<REAL> &pos, REAL &w) const{

  if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    fIntKsi->Loc(ip, pos);
    w = fIntKsi->W(ip);
    return;
  }
  if(!fIntKsi)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntTetra3D::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << std::endl;
}
//------------------------------------------------------------------------------
void TPZIntTetra3D::GetOrder(TPZVec<int> &ord) const {
  ord[0] = (short) fOrdKsi;
  ord[1] = ord[0];
  ord[2] = ord[0];
}
//##############################################################################
//##############################################################################
TPZIntPrism3D::TPZIntPrism3D(int OrdK,int OrdL) : fIntRule1D(OrdK), fIntTriang(OrdL) {
  fOrdKsi = OrdK;
  fOrdKti = OrdL;
}
TPZIntPrism3D::~TPZIntPrism3D() {}
//------------------------------------------------------------------------------
void TPZIntPrism3D::SetOrder(TPZVec<int> &ord){

  fOrdKsi = ord[0];//ordem na reta : zeta
  fOrdKti = (ord[1] > ord[2]) ? ord[1] : ord[2];//ordem no plano XY
  TPZVec<int> prc1(1),prc2(2);
  prc1[0] = ord[0];
  prc2[0] = ord[1];
  prc2[1] = ord[2];
  fIntRule1D.SetOrder(prc1);
  fIntTriang.SetOrder(prc2);
}
//------------------------------------------------------------------------------
int TPZIntPrism3D::NPoints() const {
  if (!&fIntRule1D || !&fIntTriang){
    PZError << "Null Pointer passed to method TPZIntPrism3D::NPoints()\n";
    return 0;
  }
  return fIntRule1D.NPoints()*fIntTriang.NPoints();
}
//------------------------------------------------------------------------------
void TPZIntPrism3D::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {

  if((&fIntRule1D && &fIntTriang) && ((ip >= 0) && (ip < NPoints()))){

    REAL v;
    TPZVec<REAL> ps(2);
    fIntTriang.Point(ip % fIntTriang.NPoints(), ps, v);
    pos[0] = ps[0]; pos[1] = ps[1];
    fIntRule1D.Point(ip / fIntTriang.NPoints(), ps, w);
    pos[2] = ps[0];
    w *= v;
    return;
  }
  if(!&fIntRule1D || !&fIntTriang)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntPrism3D::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << std::endl;
}
//------------------------------------------------------------------------------
void TPZIntPrism3D::GetOrder(TPZVec<int> &ord) const{
  ord[0] = fOrdKsi;
  ord[1] = fOrdKti;
  ord[2] = ord[1];
}
//##############################################################################
//##############################################################################
TPZIntPyram3D::TPZIntPyram3D(int OrdK){
  fOrdKsi = OrdK;
  fIntKsi = gIntRuleList.GetRuleP3D(OrdK);
}
//**************************************
TPZIntRuleP3D* TPZIntRuleList::GetRuleP3D(int precision) {

  // <<<<<>>>>>
  if (precision < 1 || precision > intavailP3D) {
    PZError << "\nERROR(TPZIntRuleList::getrule)-> precision = " << precision;
    //		PZError.show();
    precision = intavailP3D;
  }

  return intlistP3D[precision-1];
}
//------------------------------------------------------------------------------
void TPZIntPyram3D::SetOrder(TPZVec<int> &ord){

  fOrdKsi = (ord[1] > ord[0]) ? ord[1] : ord[0];
  fOrdKsi = (fOrdKsi > ord[2]) ? fOrdKsi : ord[2];
  fIntKsi = gIntRuleList.GetRuleP3D(fOrdKsi);
}
//------------------------------------------------------------------------------
int TPZIntPyram3D::NPoints() const{
  if (!fIntKsi){
    PZError << "Null Pointer passed to method TPZIntPyram3D::NPoints()\n";
    return 0;
  }
  return fIntKsi->NInt();
}
//------------------------------------------------------------------------------
void TPZIntPyram3D::Point(int ip, TPZVec<REAL> &pos, REAL &w) const{

  if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    fIntKsi->Loc(ip, pos);
    w = fIntKsi->W(ip);
    return;
  }
  if(!fIntKsi)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntPyram3D::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << std::endl;
}
//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////



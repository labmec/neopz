/**
 * @file
 * @brief Contains the implementation of the TPZGeoElRefLess methods.
 */
//
// C++ Interface: pzgeoelrefless
//
// Description:
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef PZGEOELREFLESS_H_H
#define PZGEOELREFLESS_H_H

#include "pzgeoelrefless.h"

#include <sstream>

#include "pzlog.h" // test 
#ifdef LOG4CXX
static LoggerPtr loggerrefless(Logger::getLogger("pz.mesh.tpzgeoelrefless"));
#endif

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess():TPZGeoEl(){
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSide();
	//  fSubElement = -1;
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(const TPZGeoElRefLess<TGeo>  &gel):TPZGeoEl(gel), fGeo(gel.fGeo){
	int i;
	for(i=0;i<TGeo::NSides;i++){
		TPZGeoElSide thisside(this->fNeighbours[i], this->Mesh());
		TPZGeoElSide gelside(gel.fNeighbours[i], this->Mesh());
		thisside.SetConnectivity(gelside);
		//     fNeighbours[i].SetConnectivity(gel.fNeighbours[i]);
	}
	//  fSubElement = -1;
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::~TPZGeoElRefLess(){
	//RemoveConnectivities();
}
/** divides the element and puts the resulting elements in the vector */
//template<class TGeo>
//void TPZGeoElRefLess<TGeo>::Divide(TPZVec < TPZGeoEl * > & pv){
//  pv.Resize (1);
//  TPZGeoEl *subel;
//  subel = new TPZGeoElRefLess(*this);
//  fSubElement = subel->Index();
//  subel->SetFather(this);
//  subel->SetFather(fIndex);
//}

/** return 1 if the element has subelements along side */
//template<class TGeo>
//int TPZGeoElRefLess<TGeo>::HasSubElement(){
//  return 0;
//}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh) :
TPZGeoEl(matind,mesh), fGeo(nodeindices) {
	
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSide();
    fGeo.Initialize(this);
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(TGeo &geo,int matind,TPZGeoMesh &mesh) :
TPZGeoEl(matind,mesh), fGeo(geo) {
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSide();
    fGeo.Initialize(this);
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index) :
TPZGeoEl(matind,mesh,index) , fGeo(nodeindices) 
{
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSide();
    fGeo.Initialize(this);
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
TPZGeoEl(id,matind,mesh) , fGeo(nodeindexes) {
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSide();
    fGeo.Initialize(this);
}

/*
 template< class TGeo>
 void TPZGeoElRefLess<TGeo>::Initialize(TPZVec<int> &nodeindices) {
 fGeo.Initialize(nodeindices);
 int i;
 for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSide();
 }
 */

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NodeIndex(int node) const {
	if(node<0 || node>=fGeo.NNodes) return -1;
	return fGeo.fNodeIndexes[node];
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideNodeIndex(int side,int node) {
	if(side<0 || side>(TGeo::NSides - 1) || node<0) {
		PZError << "TPZGeoElRefLess::SideNodeIndex. Bad parameter side.\n";
		return -1;
	}
	return fGeo.fNodeIndexes[TGeo::SideNodeLocId(side,node)];
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideNodeLocIndex(int side,int node) {
	
	if(side<0 || side>(TGeo::NSides - 1) || node<0) {
		PZError << "TPZGeoElRefLess::SideNodeIndex. Bad parameter side.\n";
		return -1;
	}
	return TGeo::SideNodeLocId(side,node);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::SetSubElement(int id, TPZGeoEl *el){
	if (id!=0)
		PZError << "TPZGeoElRefLess<TGeo>::SetSubElement - Fodeu!\n";
	else PZError << "TPZGeoElRefLess<TGeo>::SetSubElement - Por enquanto eu no faco nada!\n";
	return;
}

template<class TGeo>
TPZIntPoints *TPZGeoElRefLess<TGeo>::CreateSideIntegrationRule(int side, int order){
	return TGeo::CreateSideIntegrationRule(side,order);
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NNodes() {
	return TGeo::NNodes;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NCornerNodes(){
	return TGeo::NCornerNodes;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSides(){
	return TGeo::NSides;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideNodeLocId(int side, int node){
	return TGeo::SideNodeLocId(side,node);
}

template<class TGeo>
REAL
TPZGeoElRefLess<TGeo>::RefElVolume(){
	return TGeo::RefElVolume();
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSideNodes(int side){
	return TGeo::NSideNodes(side);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::MidSideNodeIndex(int side,int &index){
	//TRef::MidSideNodeIndex(this,side,index);
	index = -1;
	if(side<0 || side>NSides()-1) {
		PZError << "TPZGeoElRefLess<TGeo>::MidSideNodeIndex. Bad parameter side = " << side << std::endl;
		return;
	}
	if(side<NNodes()) {//o n�medio do lado 0 �o 0 etc.
		index = NodeIndex(side);
		return;
	}
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideIsUndefined(int side){
	if (side < 0 || side > NSides()){
		PZError << "TPZGeoElRefLess<TGeo>::SideIsUndefined - bad side: " << side << std::endl;
	}
	return (fNeighbours[side].Side() == -1);
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSubElements(){
	//return TRef::NSubEl;
	return 0;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSideSubElements(int side){
	//return TRef::NSideSubElements(side);
	return 0;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSideSubElements2(int side){
	// return TRef::NSideSubElements(side);
	return 0;
}

template<class TGeo>
TPZGeoEl *
TPZGeoElRefLess<TGeo>::CreateBCGeoEl(int side, int bc){
	TPZGeoEl * result = fGeo.CreateBCGeoEl(this,side,bc);
	//result->BuildBlendConnectivity();
	return result;  
}

template<class TGeo>
TPZGeoEl * TPZGeoElRefLess<TGeo>::CreateGeoElement(MElementType type,
												   TPZVec<int>& nodeindexes,
												   int matid,
												   int& index)
{
	return fGeo.CreateGeoElement(*Mesh(),type,nodeindexes,matid,index);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::SetNodeIndex(int i,int nodeindex){
	if(i<0 || i>(TGeo::NNodes - 1)){
		std::cout << "TPZGeoElRefLess::SetNodeIndex index error i = " << i << std::endl;
		return;
	}
	fGeo.fNodeIndexes[i] = nodeindex;
}

template<class TGeo>
TPZTransform
TPZGeoElRefLess<TGeo>::SideToSideTransform(int sidefrom,int sideto){
	return TGeo::SideToSideTransform(sidefrom,sideto);
}

template<class TGeo>
TPZGeoEl *
TPZGeoElRefLess<TGeo>::SubElement(int is){
	if(is<0 || is>1){//(TRef::NSubEl - 1)){
		std::cout << "TPZGeoElRefLess::SubElement index error is= " << is << std::endl;;
	}
	//  return fSubEl[is];
	return 0;
}
/*
 template<class TGeo>
 TPZGeoElSide
 TPZGeoElRefLess<TGeo>::SideSubElement(int side,int position){
 TPZStack<TPZGeoElSide> subs;
 TRef::GetSubElements(this,side,subs);
 return subs[position];
 if (side < 0 || side > NSides()){
 PZError << "TPZGeoElRefLess<TGeo>::SideSubElement - bad side: " << side << std::endl;
 return TPZGeoElSide();
 }
 return (fNeighbours[side].Side() == -1);
 }  */

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideDimension(int side){
	return TGeo::SideDimension(side);
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::Dimension(){
	return TGeo::Dimension;
}

template<class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TGeo>::HigherDimensionSides(int side,int targetdimension){
	std::cout << "TPZGeoElRefLess::HigherDimensionSides nao deve ser usado\n";
	return TPZGeoElSide();
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
	TPZStack<int> highsides;
	TGeo::HigherDimensionSides(side,highsides);
	int i,size = highsides.NElements();
	for(i=0;i<size;i++) {
		if(SideDimension(highsides[i]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,highsides[i]));
		}
	}
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::LowerDimensionSides(int side,TPZStack<int> &smallsides){
	TGeo::LowerDimensionSides(side,smallsides);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::BuildTransform(int side, TPZGeoEl *father,TPZTransform &t){
	BuildTransform2(side,father,t);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){
	fGeo.Jacobian(*this,coordinate,jac,axes,detjac,jacinv);
	
#ifdef DEBUG
	if(IsZero(detjac)){
		std::stringstream sout;
		sout << "Jacobiano nulo\n";
		LOGPZ_ERROR(loggerrefless,sout.str())
		detjac = ZeroTolerance();
	}
#endif
	//   if(TGeo::NNodes == 2) {
	//     detjac = 1.;
	//     jacinv(0,0) = 1.;
	//     jac(0,0) = 1.;
	//   }
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result){
	fGeo.X(*this,coordinate,result);
}

template<class TGeo>
bool TPZGeoElRefLess<TGeo>::IsLinearMapping() const 
{ 
	return fGeo.IsLinearMapping();
}

template<class TGeo>
bool TPZGeoElRefLess<TGeo>::IsGeoBlendEl() const 
{ 
	return fGeo.IsGeoBlendEl();
}

template<class TGeo>
TPZTransform
TPZGeoElRefLess<TGeo>::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t)
{
	if(this == father) return t;
	//Augusto:09/01/01
	TPZGeoEl *myfather = Father();
	if(side<0 || side>(TGeo::NSides-1) || !myfather){
		PZError << "TPZGeoElRefLess::BuildTransform2 side out of range or father null\n";
		return TPZTransform(0,0);
	}
	TPZGeoElSide fathloc = Father2(side);
	int son = WhichSubel();
	TPZTransform trans=myfather->GetTransform(side,son);
	trans = trans.Multiply(t);
	if(fathloc.Element() == father) return trans;
	trans = myfather->BuildTransform2(fathloc.Side(),father,trans);
	return trans;
}


template<class TGeo>
TPZTransform
TPZGeoElRefLess<TGeo>::GetTransform(int /*side*/,int /*son*/){
	//  return TRef::GetTransform(side,son);
	//  if(side<0 || side>NSides()-1){
    PZError << "TPZGeoElRefLess<TGeo>::GetTransform::Never should be called\n";
    return TPZTransform(0,0);
	//  }
	/*  int smalldim = TGeo::SideDimension(side);
	 int fatherside = FatherSide(side,son);
	 int largedim = TGeo::SideDimension(fatherside);
	 TPZTransform trans(largedim,smalldim);
	 int i,j;
	 for(i=0; i<largedim; i++) {
	 for(j=0; j<smalldim; j++) {
	 trans.Mult()(i,j) = buildt[whichsubel][side][j][i];
	 }
	 trans.Sum() (i,0) = buildt[whichsubel][side][2][i];
	 }
	 return trans;
	 */
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::CenterPoint(int side, TPZVec<REAL> &cent){
	TGeo::CenterPoint(side,cent);
}

template<class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TGeo>::Father2(int side)
{
	//std::cout << " Father2 teste Cedric: 08/05/2003\n";
	TPZGeoEl *father = Father();
	if(!father) return TPZGeoElSide();
	int son = WhichSubel();
	if(son<0) return TPZGeoElSide();
	int fathsid = father->FatherSide(side,son);
	return TPZGeoElSide(father,fathsid);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){
	//  TRef::GetSubElements(this,side,subel);
	return;
}

template<class TGeo>
void TPZGeoElRefLess<TGeo>::Read(TPZStream &buf, void *context){
	TPZGeoEl::Read(buf,context);
    fGeo.Read(buf,context);
	int i, n = TGeo::NSides;
	for(i = 0; i < n; i++){
		this->fNeighbours[i].Read(buf);
	}
}//Read

template<class TGeo>
void TPZGeoElRefLess<TGeo>::Write(TPZStream &buf, int withclassid){
	TPZGeoEl::Write(buf,withclassid);
    fGeo.Write(buf);
	int i, n = TGeo::NSides;
	for(i = 0; i < n; i++){
		this->fNeighbours[i].Write(buf);
	}
}//Write

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(TPZGeoMesh &DestMesh, const TPZGeoElRefLess &cp):TPZGeoEl(DestMesh, cp), fGeo(cp.fGeo) {
	int i;
	const int n = TGeo::NSides;
	for(i = 0; i < n; i++){
		this->fNeighbours[i] = cp.fNeighbours[i];
	}
}


template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess( TPZGeoMesh &DestMesh,
									   const TPZGeoElRefLess &cp,
									   std::map<int,int> & gl2lcNdMap,
									   std::map<int,int> & gl2lcElMap ) :
TPZGeoEl(DestMesh, cp, gl2lcElMap), fGeo(cp.fGeo, gl2lcNdMap)
{
	int i;
	const int n = TGeo::NSides;
	// #ifdef DEBUG2
	//   std::stringstream sout;
	//   sout << __PRETTY_FUNCTION__ << " for element " << Index() << std::endl;
	// #endif
	
	for(i = 0; i < n; i++)
	{
		TPZGeoElSide neigh (cp.fNeighbours[i],cp.Mesh());
		int neighIdx = neigh.Element()->Index();
		int side = neigh.Side();
		/*#ifdef DEBUG2
		 sout << "neighbour data = " << neighIdx << " / " << side << std::endl;
		 #endif*/
		while (gl2lcElMap.find(neighIdx)==gl2lcElMap.end())
		{
			neigh = neigh.Neighbour();
			neighIdx = neigh.Element()->Index();
			side = neigh.Side();
			/*#ifdef DEBUG2
			 sout << "\t\tnot neighbour = " << neighIdx << " / " << side << std::endl;
			 #endif*/
		}
		this->fNeighbours[i] = TPZGeoElSideIndex ( gl2lcElMap [ neighIdx ] , side );
		/*#ifdef DEBUG2
		 sout << "defining neighbour = " << i << " = " << gl2lcElMap [ neighIdx ] << " / " << side << std::endl;
		 #endif*/
	}
	/*#ifdef DEBUG2
	 LOGPZ_DEBUG( logger, sout.str().c_str() );
	 #endif*/
}

#include "pzgeoquad.h"
/**
 * Compute the permutation for an HDiv side
 */
template<>
inline void TPZGeoElRefLess<pzgeom::TPZGeoQuad>::HDivPermutation(int side, TPZVec<int> &permutegather)
{
	if(side < 4 || side > 7)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " called with wrong side parameter " << side;
#ifdef LOG4CXX
		LOGPZ_ERROR(loggerrefless,sout.str())
#endif
        std::cout << sout.str() << std::endl;
	}
	permutegather.Resize(3);
	int id1 = NodePtr(SideNodeLocIndex(side,0))->Id();
	int id2 = NodePtr(SideNodeLocIndex(side,1))->Id();
	if(id1<id2)
	{
		permutegather[0] = 0;
		permutegather[1] = 1;
		permutegather[2] = 2;
	}
	else
	{
		permutegather[0] = 1;
		permutegather[1] = 0;
		permutegather[2] = 2;
	}
}

#include "pzgeotriangle.h"
/**
 * Compute the permutation for an HDiv side
 */
template<>
inline void TPZGeoElRefLess<pzgeom::TPZGeoTriangle>::HDivPermutation(int side, TPZVec<int> &permutegather)
{
	if(side < 3 || side > 5)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " called with wrong side parameter " << side;
#ifdef LOG4CXX
		LOGPZ_ERROR(loggerrefless,sout.str())
#endif
        std::cout << sout.str() << std::endl;
	}
	permutegather.Resize(3);
	int id1 = NodePtr(SideNodeLocIndex(side,0))->Id();
	int id2 = NodePtr(SideNodeLocIndex(side,1))->Id();
	if(id1<id2)
	{
		permutegather[0] = 0;
		permutegather[1] = 1;
		permutegather[2] = 2;
	}
	else
	{
		permutegather[0] = 1;
		permutegather[1] = 0;
		permutegather[2] = 2;
	}
}

/**
 * Compute the permutation for an HDiv side
 */
template<class TGeo>
inline void TPZGeoElRefLess<TGeo>::HDivPermutation(int side, TPZVec<int> &permutegather)
{
	int dimension = TGeo::Dimension;
	int sidedimension = TGeo::SideDimension(side);
	
	if(dimension != sidedimension+1)
	{
		std::stringstream sout;
		sout << "HDivPermutation called with wrong side parameter " << side;
#ifdef LOG4CXX
		LOGPZ_ERROR(loggerrefless,sout.str())
#endif
	}
	TPZManVector<int,TGeo::NCornerNodes> id(TGeo::NCornerNodes);
	for(int i=0; i<TGeo::NCornerNodes; i++) id[i] = fGeo.fNodeIndexes[i];
	TGeo::GetSideHDivPermutation(side, id, permutegather);
}

//HDiv
template<>
inline void TPZGeoElRefLess<pzgeom::TPZGeoQuad>::VecHdiv(TPZFMatrix<REAL> &normalvec,TPZVec<int> &sidevector )
{
    fGeo.VecHdiv(*this,normalvec,sidevector);
}

template<>
inline void TPZGeoElRefLess<pzgeom::TPZGeoTriangle>::VecHdiv(TPZFMatrix<REAL> &normalvec,TPZVec<int> &sidevector )
{
	fGeo.VecHdiv(*this,normalvec,sidevector);
}


template<class TGeo>
inline void TPZGeoElRefLess<TGeo>::VecHdiv(TPZFMatrix<REAL> &normalvec,TPZVec<int> &sidevector )
{
    PZError << __PRETTY_FUNCTION__ << " nao esta implementado\n";
}


#endif

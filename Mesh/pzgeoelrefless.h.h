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
//
#include "pzgeoelrefless.h"

#include <sstream>

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzgeoelrefless"));
#endif

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess():TPZGeoEl(){
  int i;
  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = -1;
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
//  fSubElement = -1;
}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(const TPZGeoElRefLess<TShape,TGeo>  &gel):TPZGeoEl(gel){
  int i;
  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = gel.fNodeIndexes[i];
  for(i=0;i<TShape::NSides;i++){
    TPZGeoElSide thisside(this->fNeighbours[i], this->Mesh());
    TPZGeoElSide gelside(gel.fNeighbours[i], this->Mesh());
    thisside.SetConnectivity(gelside);
//     fNeighbours[i].SetConnectivity(gel.fNeighbours[i]);
  }
//  fSubElement = -1;
}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::~TPZGeoElRefLess(){
  //RemoveConnectivities();
}
/** divides the element and puts the resulting elements in the vector */
//template<class TShape, class TGeo>
//void TPZGeoElRefLess<TShape,TGeo>::Divide(TPZVec < TPZGeoEl * > & pv){
//  pv.Resize (1);
//  TPZGeoEl *subel;
//  subel = new TPZGeoElRefLess(*this);
//  fSubElement = subel->Index();
//  subel->SetFather(this);
//  subel->SetFather(fIndex);
//}

/** return 1 if the element has subelements along side */
//template<class TShape, class TGeo>
//int TPZGeoElRefLess<TShape,TGeo>::HasSubElement(){
//  return 0;
//}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(matind,mesh) {

  int i,nnod = nodeindices.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElRefLess<TShape,TGeo>::Constuctor, number of nodes : " << nnod << std::endl;
    return;
  }

  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindices[i];
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index) :
  TPZGeoEl(matind,mesh,index) {

  int i,nnod = nodeindices.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElRefLess<TShape,TGeo>::Constuctor, number of nodes : " << nnod << std::endl;
    return;
  }

  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindices[i];
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
}

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoEl(id,matind,mesh) {
  int i,nnod = nodeindexes.NElements();
  if(nnod!=TGeo::NNodes) {
    PZError << "TPZGeoElRefLess<TShape,TGeo>::Constuctor, number of nodes : " << nnod << std::endl;
    return;
  }

  for(i=0;i<TGeo::NNodes;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
}

template< class TShape, class TGeo>
void TPZGeoElRefLess<TShape,TGeo>::Initialize(TPZVec<int> &nodeindices, int matind, TPZGeoMesh& mesh, int& index ) {
  int i;
  for(i = 0; i < TGeo::NNodes; i++ ){
     fNodeIndexes[ i ] = nodeindices[i];
  }
// for( int i = 0; i < TRef::NSubEl; i++ ){
//     fSubEl[ i ] = 0;
 // }
  for(i=0;i<TShape::NSides;i++)fNeighbours[i] = TPZGeoElSide();
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NodeIndex(int node) {
  if(node<0 || node>7) return -1;
  return fNodeIndexes[node];
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideNodeIndex(int side,int node) {
  if(side<0 || side>(TShape::NSides - 1) || node<0) {
    PZError << "TPZGeoElRefLess::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  return fNodeIndexes[TShape::SideNodeLocId(side,node)];
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideNodeLocIndex(int side,int node) {

  if(side<0 || side>(TShape::NSides - 1) || node<0) {
    PZError << "TPZGeoElRefLess::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  return TShape::SideNodeLocId(side,node);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::SetSubElement(int id, TPZGeoEl *el){
  if (id!=0)
    PZError << "TPZGeoElRefLess<TShape,TGeo>::SetSubElement - Fodeu!\n";
  else PZError << "TPZGeoElRefLess<TShape,TGeo>::SetSubElement - Por enquanto eu no faco nada!\n";
  return;
}

template<class TShape, class TGeo>
TPZIntPoints *TPZGeoElRefLess<TShape,TGeo>::CreateSideIntegrationRule(int side, int order){
  return TGeo::CreateSideIntegrationRule(side,order);
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NNodes() {
  return TGeo::NNodes;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NCornerNodes(){
  return TShape::NNodes;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSides(){
  return TGeo::NSides;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideNodeLocId(int side, int node){
  return TShape::SideNodeLocId(side,node);
}

template<class TShape, class TGeo>
REAL
TPZGeoElRefLess<TShape,TGeo>::RefElVolume(){
  return TShape::RefElVolume();
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSideNodes(int side){
  return TShape::NSideNodes(side);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::MidSideNodeIndex(int side,int &index){
  //TRef::MidSideNodeIndex(this,side,index);
  index = -1;
  if(side<0 || side>NSides()-1) {
    PZError << "TPZGeoElRefLess<TShape,TGeo>::MidSideNodeIndex. Bad parameter side = " << side << std::endl;
    return;
  }
  if(side<NNodes()) {//o n�medio do lado 0 �o 0 etc.
    index = NodeIndex(side);
    return;
  }
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideIsUndefined(int side){
  if (side < 0 || side > NSides()){
    PZError << "TPZGeoElRefLess<TShape,TGeo>::SideIsUndefined - bad side: " << side << std::endl;
  }
  return (fNeighbours[side].Side() == -1);
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSubElements(){
  //return TRef::NSubEl;
  return 0;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSideSubElements(int side){
  //return TRef::NSideSubElements(side);
  return 0;
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::NSideSubElements2(int side){
  // return TRef::NSideSubElements(side);
  return 0;
}

template<class TShape, class TGeo>
TPZGeoEl *
TPZGeoElRefLess<TShape,TGeo>::CreateBCGeoEl(int side, int bc){
  return TGeo::CreateBCGeoEl(this,side,bc);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::SetNodeIndex(int i,int nodeindex){
  if(i<0 || i>(TGeo::NNodes - 1)){
    std::cout << "TPZGeoElRefLess::SetNodeIndex index error i = " << i << std::endl;
    return;
  }
  fNodeIndexes[i] = nodeindex;
}

template<class TShape, class TGeo>
TPZTransform
TPZGeoElRefLess<TShape,TGeo>::SideToSideTransform(int sidefrom,int sideto){
  return TShape::SideToSideTransform(sidefrom,sideto);
}

template<class TShape, class TGeo>
TPZGeoEl *
TPZGeoElRefLess<TShape,TGeo>::SubElement(int is){
  if(is<0 || is>1){//(TRef::NSubEl - 1)){
    std::cout << "TPZGeoElRefLess::SubElement index error is= " << is << std::endl;;
  }
//  return fSubEl[is];
  return 0;
}
/*
template<class TShape, class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TShape,TGeo>::SideSubElement(int side,int position){
 TPZStack<TPZGeoElSide> subs;
  TRef::GetSubElements(this,side,subs);
  return subs[position];
  if (side < 0 || side > NSides()){
    PZError << "TPZGeoElRefLess<TShape,TGeo>::SideSubElement - bad side: " << side << std::endl;
    return TPZGeoElSide();
  }
  return (fNeighbours[side].Side() == -1);
}  */

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::SideDimension(int side){
  return TShape::SideDimension(side);
}

template<class TShape, class TGeo>
int
TPZGeoElRefLess<TShape,TGeo>::Dimension(){
  return TShape::Dimension;
}

template<class TShape, class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TShape,TGeo>::HigherDimensionSides(int side,int targetdimension){
  std::cout << "TPZGeoElRefLess::HigherDimensionSides nao deve ser usado\n";
  return TPZGeoElSide();
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
  TPZStack<int> highsides;
  TShape::HigherDimensionSides(side,highsides);
  int i,size = highsides.NElements();
  for(i=0;i<size;i++) {
    if(SideDimension(highsides[i]) == targetdimension) {
      elsides.Push(TPZGeoElSide(this,highsides[i]));
    }
  }
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::LowerDimensionSides(int side,TPZStack<int> &smallsides){
  int nsidecon = TShape::NSideConnects(side);
  int is;
  for(is=0; is<nsidecon-1; is++)
    smallsides.Push(TShape::SideConnectLocId(side,is));
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::BuildTransform(int side, TPZGeoEl *father,TPZTransform &t){
  BuildTransform2(side,father,t);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix &jac,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
  TPZFNMatrix<3*TGeo::NNodes> nodes(3,TGeo::NNodes);
  TPZGeoNode *np;
  TPZAdmChunkVector<TPZGeoNode> &nodevec = Mesh()->NodeVec();
  int i,j;
  for(i=0;i<TGeo::NNodes;i++) {
    np = &nodevec[fNodeIndexes[i]];
    for(j=0;j<3;j++) {
      nodes(j,i) = np->Coord(j);
    }
  }
  TGeo::Jacobian(nodes,coordinate,jac,axes,detjac,jacinv);
//   if(TGeo::NNodes == 2) {
//     detjac = 1.;
//     jacinv(0,0) = 1.;
//     jac(0,0) = 1.;
//   }
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result){
  TPZFNMatrix<3*TGeo::NNodes> nodes(3,TGeo::NNodes);
  TPZGeoNode *np;
  TPZAdmChunkVector<TPZGeoNode> &nodevec = Mesh()->NodeVec();
  int i,j;
  for(i=0;i<TGeo::NNodes;i++) {
    np = &nodevec[fNodeIndexes[i]];
    for(j=0;j<3;j++) {
      nodes(j,i) = np->Coord(j);
    }
  }
  TGeo::X(nodes,coordinate,result);
}

template<class TShape, class TGeo>
TPZTransform
TPZGeoElRefLess<TShape,TGeo>::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t)
{
  //Augusto:09/01/01
  TPZGeoEl *myfather = Father();
  if(side<0 || side>(TShape::NSides-1) || !myfather){
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


template<class TShape, class TGeo>
TPZTransform
TPZGeoElRefLess<TShape,TGeo>::GetTransform(int side,int son){
//  return TRef::GetTransform(side,son);
//  if(side<0 || side>NSides()-1){
    PZError << "TPZGeoElRefLess<TShape,TGeo>::GetTransform::Never should be called\n";
    return TPZTransform(0,0);
//  }
/*  int smalldim = TShape::SideDimension(side);
  int fatherside = FatherSide(side,son);
  int largedim = TShape::SideDimension(fatherside);
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

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::CenterPoint(int side, TPZVec<REAL> &cent){
  TShape::CenterPoint(side,cent);
}

template<class TShape, class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TShape,TGeo>::Father2(int side)
{
  //std::cout << " Father2 teste Cedric: 08/05/2003\n";
  TPZGeoEl *father = Father();
  if(!father) return TPZGeoElSide();
  int son = WhichSubel();
  if(son<0) return TPZGeoElSide();
  int fathsid = father->FatherSide(side,son);
  return TPZGeoElSide(father,fathsid);
}

template<class TShape, class TGeo>
void
TPZGeoElRefLess<TShape,TGeo>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){
//  TRef::GetSubElements(this,side,subel);
  return;
}

template<class TShape, class TGeo>
void TPZGeoElRefLess<TShape,TGeo>::Read(TPZStream &buf, void *context){
  TPZGeoEl::Read(buf,context);
  buf.Read(fNodeIndexes,TGeo::NNodes);
  int i, n = TShape::NSides;
  for(i = 0; i < n; i++){
    this->fNeighbours[i].Read(buf);
  }
}//Read

template<class TShape, class TGeo>
void TPZGeoElRefLess<TShape,TGeo>::Write(TPZStream &buf, int withclassid){
  TPZGeoEl::Write(buf,withclassid);
  buf.Write(fNodeIndexes,TGeo::NNodes);
  int i, n = TShape::NSides;
  for(i = 0; i < n; i++){
    this->fNeighbours[i].Write(buf);
  }
}//Write

template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(TPZGeoMesh &DestMesh, const TPZGeoElRefLess &cp):TPZGeoEl(DestMesh, cp){
  int i, n = TGeo::NNodes;
  for(i = 0; i < n; i++){
    this->fNodeIndexes[i] = cp.fNodeIndexes[i];
  }
  n = TShape::NSides;
  for(i = 0; i < n; i++){
    this->fNeighbours[i] = cp.fNeighbours[i];
  }
}


template<class TShape, class TGeo>
TPZGeoElRefLess<TShape,TGeo>::TPZGeoElRefLess(TPZGeoMesh &DestMesh,
                                              const TPZGeoElRefLess &cp,
                                              std::map<int,int> & gl2lcNdMap,
                                              std::map<int,int> & gl2lcElMap ) :
                                              TPZGeoEl(DestMesh, cp, gl2lcElMap)
{
  int i, n = TGeo::NNodes;
  for(i = 0; i < n; i++)
  {
    if (gl2lcNdMap.find(cp.fNodeIndexes[i]) == gl2lcNdMap.end())
    {
      std::stringstream sout;
      sout << "ERROR in - " << __PRETTY_FUNCTION__
           << " trying to clone a node " << i << " index " << cp.fNodeIndexes[i]
           << " wich is not mapped";
      LOGPZ_ERROR(logger,sout.str().c_str());
      exit(-1);
    }
    this->fNodeIndexes[i] = gl2lcNdMap [ cp.fNodeIndexes[i] ];
  }
  n = TShape::NSides;
  for(i = 0; i < n; i++)
  {
    TPZGeoElSide neigh (cp.fNeighbours[i],cp.Mesh());
    int neighIdx = neigh.Element()->Index();
    int side = neigh.Side();
    while (gl2lcElMap.find(neighIdx)==gl2lcElMap.end())
    {
      neigh = neigh.Neighbour();
      neighIdx = neigh.Element()->Index();
      side = neigh.Side();
    }
    this->fNeighbours[i] = TPZGeoElSideIndex ( gl2lcElMap [ neighIdx ] , side );
  }
}




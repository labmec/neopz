/***************************************************************************
                          tpzgeoelrefpattern.cc  -  description
                             -------------------
    begin                : Tue Dec 23 2003
    copyright            : (C) 2003 by LabMeC - DES - FEC - UNICAMP (Edimar Cesar Rylo) & EMBRAER
    email                : cesar@labmec.fec.unicamp.br
 ***************************************************************************/

#include "tpzgeoelrefpattern.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
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
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "TPZGeoElement.h"
#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern():TPZGeoElRefLess<TShape,TGeo>(){
  fSubEl = 0;
  fRefPattern = 0;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::~TPZGeoElRefPattern(){
  if (fRefPattern) delete fRefPattern;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(nodeindices,matind,mesh) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl = 0;
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(nodeindices,matind,mesh,index) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(0);
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(id,nodeindexes,matind,mesh) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(0);
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::Initialize : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::Initialize(TPZVec<int> &nodeindices, int matind, TPZGeoMesh& mesh, int& index,TPZRefPattern *refpat) {
  TPZGeoElRefLess<TShape,TGeo>::Initialize(nodeindices,matind,mesh,index);
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(0);
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::Initialize : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = 0;
}

template<class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::SetSubElement(int id, TPZGeoEl *el){
  int nsubel = fRefPattern->NSubElements();
  if (id<0 || id >nsubel){
    PZError << "TPZGeoElRefPattern::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id] = el;
  return;
}

template<class TShape, class TGeo>
REAL TPZGeoElRefPattern<TShape,TGeo>::RefElVolume(){
  return TShape::RefElVolume();
}

template<class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::MidSideNodeIndex(int side,int &index){
  index = -1;
  int i,j;
  if(side<0 || side>NSides()-1) {
    PZError << "TPZGeoElRefPattern::MidSideNodeIndex. Bad parameter side = " << side << endl;
    return;
  }
  //corner side
  if(side<NCornerNodes()) {//o nó medio do lado 0 é o 0 etc.
    index = NodeIndex(side);
    return;
  }
  //o nó medio da face é o centro e o nó medio do centro é o centro
  //como nó de algum filho se este existir
  //caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
  if(HasSubElement()) {
    //side -= NCornerNodes();
    //index =(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
    int nsubel = NSideSubElements(side);
    TPZStack <TPZGeoElSide> subels;
    GetSubElements2(side,subels);
    if (nsubel == 1) {
      subels[0].Element()->MidSideNodeIndex(subels[0].Side(),index);
      return;
    }
    TPZStack <int> msnindex;
    int subnodeindex = subels[0].SideNodeIndex(subels[0].Side());
    msnindex.Push(subnodeindex);
    for (i=1;i<nsubel;i++){
      subnodeindex = subels[i].SideNodeIndex(subels[i].Side());
      for (j=0;j<msnindex.NElements();j++){
        if (subnodeindex != msnindex[j]){
          cout << "TPZGeoElRefPattern<TShape,TGeo>::MidSideNodeIndex element with more than one midsidenodeindex..." << endl;
          msnindex.Push(subnodeindex);
        }
      }
    }
    if (msnindex.NElements() == 1) index = msnindex[0];
    else {
      TPZVec<REAL> centerpt(Dimension(),0.);
      TPZVec<REAL> auxpt(Dimension(),0.);      
      CenterPoint(side,auxpt);
      X(auxpt,centerpt);

      REAL dif = 1e6;
      REAL aux = 0;
      index = msnindex[0];
      for (i=0;i<msnindex.NElements();i++){
        for (j=0;j<Dimension();j++)
          aux += (Mesh()->NodeVec()[msnindex[i]].Coord(j) - centerpt[j]) *
                    (Mesh()->NodeVec()[msnindex[i]].Coord(j) - centerpt[j]);
        if (aux < dif){
          dif = aux;
          index = msnindex[i];
        }
      }
    }
  }
}

template<class TShape, class TGeo>
int TPZGeoElRefPattern<TShape,TGeo>::NSubElements(){
  return fRefPattern->NSubElements();
}

template<class TShape, class TGeo>
int TPZGeoElRefPattern<TShape,TGeo>::NSideSubElements2(int side){
  return fRefPattern->NSideSubElements(side);
}

template<class TShape, class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TShape,TGeo>::SubElement(int is){
  int nsubel = fRefPattern->NSubElements();
  if(is<0 || is>nsubel){
    cout << "TPZGeoElRefPattern::SubElement index error is= " << is << endl;;
  }
  return fSubEl[is];
}

template<class TShape, class TGeo>
TPZGeoElSide TPZGeoElRefPattern<TShape,TGeo>::SideSubElement(int side,int position){
//  TPZStack<TPZGeoElSide> subs;
//  
//  TRef::GetSubElements(this,side,subs);
//  return subs[position];
  int sub, sideout;
  fRefPattern->SideSubElement(side,position,sub,sideout);
  if (!fSubEl[sub]) {
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::SideSubElement : Error subelement not found for side "
            << side << " position " << position << endl;
    return TPZGeoElSide();
  }
  return TPZGeoElSide (fSubEl[sub],sideout);
}

template<class TShape, class TGeo>
TPZTransform TPZGeoElRefPattern<TShape,TGeo>::GetTransform(int side,int son){
//  return TRef::GetTransform(side,son);
  return fRefPattern->Transform(side,son);  
}

template<class TShape, class TGeo>
void
TPZGeoElRefPattern<TShape,TGeo>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){
  //TRef::GetSubElements(this,side,subel);
  int i,nsidesubel = fRefPattern->NSideSubElements(side);
  for (i=0;i<nsidesubel;i++){
    subel.Push(SideSubElement(side,i));  
  }
}


template<class TShape, class TGeo>
void
TPZGeoElRefPattern<TShape,TGeo>::Divide(TPZVec<TPZGeoEl *> &SubElVec){
  int i;
  int NSubEl = fRefPattern->NSubElements();
  if(HasSubElement()) {
    SubElVec.Resize(NSubEl);
    for(i=0;i<NSubEl;i++) SubElVec[i] = SubElement(i);
    return;//If exist fSubEl return this sons
  }
  int index,k,j,sub,matid=MaterialId();

  const int totalnodes = fRefPattern->Mesh()->NodeVec().NElements();  
  int np[totalnodes];

  for(j=0;j<NNodes();j++) np[j] = NodeIndex(j);
  //verifico no padro de refinamento os nos que compoem o lado
  //caso tenha nos que nao sejam nos de canto eu verificarei a necessidade
  //de cria-los.
  //Sera necessario nos seguintes casos: nenhum subelemento o criou ainda,
  //nenhum vizinho o possui ou ainda esses nao coincidam com os nos dos vizinhos
  //o que criaria uma inconsistencia nos elementos continuos ainda nao contemplada
  //no momento ??mas que pode vir a ser considerada??
  int sum = NCornerNodes(); 
  for (i=NNodes();i<NSides();i++){
    if (fRefPattern->NSideSubElements(i)<=1) continue;

    //equivalente as coordenadas no elemento mestre...
    TPZVec<REAL> refnodecoord(3,0.);
    for (k=0;k<3;k++) refnodecoord[k] = fRefPattern->Mesh()->NodeVec()[i].Coord(k);
    //passando para as coordenadas do elemento da malha real...
    TPZVec<REAL> newnodecoord(3,0.);
    X(refnodecoord,newnodecoord);
    
    TPZVec<int> sidenodes;
    fRefPattern->SideNodes(i,sidenodes);
    for (j=0;j<sidenodes.NElements();j++){
      index = -1;
      //De fato nao seria necessario pois criarei todos os nos antes de criar os subelementos
      MidSideNodeIndex(i,index);
      if(index < 0) {
        //verificar se um vizinho ja criou o no
        TPZGeoElSide gelside = Neighbour(i);
        if(gelside.Element() && gelside.Side()>-1) {
          bool isRefNode = false;
          while(gelside.Element() != this && !isRefNode) {
            gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
            gelside = gelside.Neighbour();
            if (index > -1){
              REAL dif = 0.;
              for (k=0;k<3;k++) {
                REAL indexcoord = Mesh()->NodeVec()[index].Coord(k);
                dif += (newnodecoord[k] - indexcoord) * (newnodecoord[k] - indexcoord);
              }
              if (dif < 1e-12) isRefNode = true;
              else index = -1;
            }
          }
        }
        if (index > -1) {
          np[sum] = index;
          sum ++;
          continue;
        }
        //Caso o no nao exista nos vizinhos sera necessario cria-lo...
        index = Mesh()->NodeVec().AllocateNewElement();
        Mesh()->NodeVec()[index].Initialize(newnodecoord,*Mesh());
        np[sum] = index;
        sum ++;
      }
    }
  }

  // creating new subelements
  int nsubel = fRefPattern->NSubElements();
  for(i=0;i<nsubel;i++) {
    int subcorner = fRefPattern->Element(i+1)->NCornerNodes();
    TPZManVector<int> cornerindexes(subcorner);
    for(j=0;j<subcorner;j++) {
      int cornerid = fRefPattern->Element(i+1)->NodeIndex(j);
      cornerindexes[j] = np[cornerid];
    }
    TPZGeoEl *subel = Mesh()->CreateGeoElement((MElementType)fRefPattern->Element(i+1)->Type(),cornerindexes,matid,index);
    SetSubElement(i,subel);
  }

  SubElVec.Resize(NSubEl);
  for(sub=0;sub<NSubEl;sub++) {
    SubElVec[sub] = SubElement(sub);
    SubElVec[sub]->SetFather(this);
  }

  //conectividades entre os filhos : viz interna  
  for(i=0;i<NSubEl;i++) {
    TPZGeoEl *refsubel = fRefPattern->Element(i+1);
    for (j=0;j<refsubel->NSides();j++){
      TPZGeoElSide refneigh = refsubel->Neighbour(j);
      if (refneigh.Element() && refneigh.Side() > -1){
        //como nao me deixam guardar o index do elemento geometrico...
        for (k=0;k<NSubEl;k++){
          if (fRefPattern->Element(k+1) == refneigh.Element()){
            TPZGeoElSide neighbour (SubElVec[k],refneigh.Side());
            SubElVec[i]->SetNeighbour(j,neighbour);
          }
        }
      }
    }
  }  
  SetSubElementConnectivities();
}


template class TPZGeoElRefPattern<TPZShapeCube,TPZGeoCube>;
template class TPZGeoElRefPattern<TPZShapeLinear,TPZGeoLinear>;
template class TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad>;
template class TPZGeoElRefPattern<TPZShapeTriang,TPZGeoTriangle>;
template class TPZGeoElRefPattern<TPZShapePrism,TPZGeoPrism>;
template class TPZGeoElRefPattern<TPZShapeTetra,TPZGeoTetrahedra>;
template class TPZGeoElRefPattern<TPZShapePiram,TPZGeoPyramid>;
template class TPZGeoElRefPattern<TPZShapePoint,TPZGeoPoint>;

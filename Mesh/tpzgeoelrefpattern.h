/***************************************************************************
                          tpzgeoelrefpattern.h  -  description
                             -------------------
    begin                : Tue Dec 23 2003
    copyright            : (C) 2003 by cesar
    email                : cesar@labmec.fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TPZGEOELREFPATTERN_H
#define TPZGEOELREFPATTERN_H
#include "pzgeoelrefless.h"
#include "TPZRefPattern.h"
#include "pzvec.h"


/**Implements a geometric element which division is given by a TPZRefPattern object.
  *@author Edimar Cesar Rylo
  *@since 2003-12-23
  */
class TPZGeoElSide;
class TPZCompMesh;
class TPZCompEl;
template<class T,int N>
class TPZStack;
class TPZRefPattern;

template <class TShape, class TGeo>  

class TPZGeoElRefPattern : public TPZGeoElRefLess<TShape,TGeo>  {

    TPZVec<int> fSubEl;
    TPZRefPattern *fRefPattern;
  
public: 
  TPZGeoElRefPattern();
  ~TPZGeoElRefPattern();
  TPZGeoElRefPattern(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat);
  TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat);
  TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index,TPZRefPattern *refpat);
  void Initialize(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index,TPZRefPattern *refpat);

  /** return 1 if the element has subelements along side*/
  int HasSubElement() {return fSubEl.NElements() && fSubEl[0]!=-1;}

  void SetSubElement(int id, TPZGeoEl *el);

  /**volume of the master element*/
  REAL RefElVolume();

  /**returns the midside node index along a side of the element*/
  void MidSideNodeIndex(int side,int &index);

  /**returns the midside node indices along a side of the element*/
  void MidSideNodeIndices(int side,TPZVec<int> &indices);

  /**
   * return the number of subelements of the element independent of the
   * fact hether the element has already been refined or not
   */
  int NSubElements();

  /**
  * return the number of subelements as returned by GetSubElements2(side)
  */
  int NSideSubElements2(int side);

  /**returns a pointer to the subelement is*/
  TPZGeoEl *SubElement(int is);

  /**return a pointer and a side of the subelement of the element at the side
     and the indicated position. position = 0 indicate first subelement, ...*/
  TPZGeoElSide SideSubElement(int side,int position);

  TPZTransform GetTransform(int side,int son);

  virtual int FatherSide(int side, int son);

  /**divides the element and puts the resulting elements in the vector*/
  virtual void Divide(TPZVec<TPZGeoEl *> &pv);

  virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel);
  
  /** Defines the refinement pattern. It's used only in TPZGeoElRefPattern objects. */
  virtual void SetRefPattern(TPZRefPattern *refpat );

};

//--| IMPLEMENTATION |----------------------------------------------------------

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern():TPZGeoElRefLess<TShape,TGeo>(){
  fSubEl.Resize(1);
  fSubEl[0] = -1;
  fRefPattern = 0;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::~TPZGeoElRefPattern(){
//#warning "Acredito que os objetos deste tipo devem ser administrados pela malha"
//  if (fRefPattern) delete fRefPattern;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(nodeindices,matind,mesh) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(1);
    fSubEl[0] = -1;
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern : NULL refinement pattern given" << endl;
    return;
  }
  this->Mesh()->InsertRefPattern(refpat);
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = -1;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh, int &index,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(nodeindices,matind,mesh,index) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(1);
    fSubEl [0] = -1;
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = -1;
}

template<class TShape, class TGeo>
TPZGeoElRefPattern<TShape,TGeo>::TPZGeoElRefPattern(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat) :
  TPZGeoElRefLess<TShape,TGeo>(id,nodeindexes,matind,mesh) {
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(1);
    fSubEl[0] = -1;
   // PZError << "TPZGeoElRefPattern<TShape,TGeo>::Initialize : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = -1;
}

template<class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::Initialize(TPZVec<int> &nodeindices, int matind, TPZGeoMesh& mesh, int& index,TPZRefPattern *refpat) {
  TPZGeoElRefLess<TShape,TGeo>::Initialize(nodeindices,matind,mesh,index);
  if (!refpat){
    fRefPattern = 0;
    fSubEl.Resize(1);
    fSubEl [0] = -1;
 //   PZError << "TPZGeoElRefPattern<TShape,TGeo>::Initialize : NULL refinement pattern given" << endl;
    return;
  }
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = -1;
}

template<class TShape, class TGeo>
void TPZGeoElRefPattern<TShape,TGeo>::SetSubElement(int id, TPZGeoEl *el){
  int nsubel = fRefPattern->NSubElements();
  if (id<0 || id >nsubel){
    PZError << "TPZGeoElRefPattern::Trying do define subelement :" << id << endl;
    return;
  }
  fSubEl[id] = el->Index();
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
  if(side<0 || side>this->NSides()-1) {
    PZError << "TPZGeoElRefPattern::MidSideNodeIndex. Bad parameter side = " << side << endl;
    return;
  }
  //corner side
  if(side<this->NCornerNodes()) {//o nó medio do lado 0 é o 0 etc.
    index = this->NodeIndex(side);
    return;
  }
  //o nó medio da face é o centro e o nó medio do centro é o centro
  //como nó de algum filho se este existir
  //caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
  if(HasSubElement()) {
    //side -= NCornerNodes();
    //index =(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
    int nsubel = this->NSideSubElements(side);
    TPZStack <TPZGeoElSide> subels;
    GetSubElements2(side,subels);
    // Nao sei porque deve se fazer esta execao
    if (nsubel == 1) {
      subels[0].Element()->MidSideNodeIndex(subels[0].Side(),index);
      return;
    }
    TPZStack <int> msnindex;
    // este sidenodeindex pode nao existir. Normalmente o numero de nos de um elemento e igual
    // NNodes. Quer dizer se o lado e maior igual NNodes, este metodo nao devolvera nada

    int subnodeindex;
    for (i=0;i<nsubel;i++){
      if(subels[i].Side() >= subels[i].Element()->NCornerNodes()) continue;
      subnodeindex = subels[i].SideNodeIndex(0);
      msnindex.Push(subnodeindex);
    }
    if(msnindex.NElements() > 1) {
      cout << "TPZGeoElRefPattern<TShape,TGeo>::MidSideNodeIndex element with more than one midsidenodeindex..." << endl;
    }
    if (msnindex.NElements() == 1) index = msnindex[0];
    else if(msnindex.NElements() == 0) {
      index = -1;
    } else {
      TPZVec<REAL> centerpt(this->Dimension(),0.);
      TPZVec<REAL> auxpt(this->Dimension(),0.);      
      this->CenterPoint(side,auxpt);
      this->X(auxpt,centerpt);

      REAL dif = 1e6;
      REAL aux = 0;
      index = msnindex[0];
      for (i=0;i<msnindex.NElements();i++){
        for (j=0;j<this->Dimension();j++)
          aux += (this->Mesh()->NodeVec()[msnindex[i]].Coord(j) - centerpt[j]) *
                    (this->Mesh()->NodeVec()[msnindex[i]].Coord(j) - centerpt[j]);
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
  if (!fRefPattern) return 0;
  return fRefPattern->NSubElements();
}

template<class TShape, class TGeo>
int TPZGeoElRefPattern<TShape,TGeo>::NSideSubElements2(int side){
  if (!fRefPattern) return 0;
  return fRefPattern->NSideSubElements(side);
}

template<class TShape, class TGeo>
TPZGeoEl * TPZGeoElRefPattern<TShape,TGeo>::SubElement(int is){
  if (!fRefPattern) return 0;
  int nsubel = fRefPattern->NSubElements();
  if(is<0 || is>nsubel){
    cout << "TPZGeoElRefPattern::SubElement index error is= " << is << endl;;
  }
  return (fSubEl[is] == -1) ? 0 : this->Mesh()->ElementVec()[fSubEl[is]];
}

template<class TShape, class TGeo>
TPZGeoElSide TPZGeoElRefPattern<TShape,TGeo>::SideSubElement(int side,int position){
//  TPZStack<TPZGeoElSide> subs;
//  
//  TRef::GetSubElements(this,side,subs);
//  return subs[position];
  int sub, sideout;
  fRefPattern->SideSubElement(side,position,sub,sideout);
  if (fSubEl[sub] == -1) {
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::SideSubElement : Error subelement not found for side "
            << side << " position " << position << endl;
    return TPZGeoElSide();
  }
  return TPZGeoElSide (SubElement(sub),sideout);
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
//  int i,nsidesubel;

  //A classe TPZRefPattern nao esta contemplando os nos
  TPZGeoEl * reffather = fRefPattern->Element(0);
  if (side < reffather->NCornerNodes()){
    TPZGeoElSide thisside (reffather,side);
    TPZGeoElSide neighbour = thisside.Neighbour();
    while(neighbour.Exists() && neighbour != thisside){
      TPZGeoEl *gel = neighbour.Element();
      TPZGeoEl *father = gel->Father();
      if (father == reffather) {
        int sonid = neighbour.Element()->Id()-1; //o id 0 e sempre do pai por definicao da classe
        int sonside = neighbour.Side();
        TPZGeoElSide sideson (SubElement(sonid),sonside);
        subel.Push(sideson);
      }
      neighbour = neighbour.Neighbour();
    }
    return;
  }
  fRefPattern->SidePartition(subel,side);
  int size = subel.NElements();
  int el;
  for(el=0; el<size; el++) {
    TPZGeoEl *gelref = subel[el].Element();
    int subelindex = gelref->Id()-1;
    TPZGeoEl *mysub = SubElement(subelindex);
    subel[el] = TPZGeoElSide(mysub,subel[el].Side());
  }
}


template<class TShape, class TGeo>
void
TPZGeoElRefPattern<TShape,TGeo>::Divide(TPZVec<TPZGeoEl *> &SubElVec){
  if (!fRefPattern) {
    PZError << "TPZGeoElRefPattern<TShape,TGeo>::Divide ERROR : Undefined Refinement Pattern!" << endl;
    SubElVec.Resize(0);
    return;
  }
  int i;
  int NSubEl = fRefPattern->NSubElements();
  if(HasSubElement()) {
    SubElVec.Resize(NSubEl);
    for(i=0;i<NSubEl;i++) SubElVec[i] = SubElement(i);
    return;//If exist fSubEl return this sons
  }
  int index,k,j,sub,matid=this->MaterialId();

  int totalnodes = fRefPattern->NNodes();  
  TPZManVector<int> np(totalnodes);
  int nnodes = this->NNodes();

  for(j=0;j<nnodes;j++) {
    np[j] = this->NodeIndex(j);
  }
  fRefPattern->CreateNewNodes (this, np);
  cout << "NewNodeIndexes : " << np << endl;


  // I dont think the previous part will work...
  // It should be better structured. One method which is not present within
  // the geometric element is returning the set of nodes present at the
  // interior of a side.
  // I suggest : MidsideNodeIndices(int side, TPZVec<int> &indices)
  // Its default behaviour can be to insert a single node.

  // creating new subelements
  for(i=0;i<NSubEl;i++) {
    int subcorner = fRefPattern->Element(i+1)->NCornerNodes();
    TPZManVector<int> cornerindexes(subcorner);
    for(j=0;j<subcorner;j++) {
      int cornerid = fRefPattern->Element(i+1)->NodeIndex(j);
      cornerindexes[j] = np[cornerid];
    }
    cout << "Subel corner " << cornerindexes << endl;
    TPZGeoEl *subel = this->Mesh()->CreateGeoElement((MElementType)fRefPattern->Element(i+1)->Type(),cornerindexes,matid,index,1);
    SetSubElement(i,subel);
  }

  SubElVec.Resize(NSubEl);
  for(sub=0;sub<NSubEl;sub++) {
    SubElVec[sub] = SubElement(sub);
    SubElVec[sub]->SetFather(this);
    SubElVec[sub]->SetFather(this->fIndex);
  }

  this->Mesh()->Print(cout);

//  fRefPattern->Mesh()->Print(cout);
  
  cout << "conectividades entre os filhos : viz interna  " << endl;
  for(i=0;i<NSubEl;i++) {
    cout << "SubElement " << i << endl;
    TPZGeoEl *refsubel = fRefPattern->Element(i+1);
    for (j=0;j<refsubel->NSides();j++){
      cout << "\t Side " << j << endl;
      TPZGeoElSide thisside (refsubel,j);
      TPZGeoElSide refneigh = refsubel->Neighbour(j);
      if(refneigh.Exists() && refneigh != thisside) {
        //como nao me deixam guardar o index do elemento geometrico...
        for (k=0;k<NSubEl+1;k++){
          if (fRefPattern->Element(k) == refneigh.Element() && k-1 != i){
            cout << "\t\t\tFound subel " << k << "  side of subel " << refneigh.Side()
                << "For subelement " << i << " side " << j << endl;
            TPZGeoElSide neighbour;
            if(k) {
              neighbour = TPZGeoElSide(SubElement(k-1),refneigh.Side());
            } else {
              neighbour = TPZGeoElSide(this, refneigh.Side());
            }
            TPZGeoElSide gs(SubElement(i),j);
            if(!gs.NeighbourExists(neighbour)) gs.SetConnectivity(neighbour);
            break;
          }
        }
        if(k == NSubEl+1) {
          cout << "TPZGeoElRefPattern<TShape,TGeo>::Divide inconsistent datastructure subelement "
            << i << " Side " << j << endl;
        }
      } else {
        SubElement(i)->SetNeighbour(j,TPZGeoElSide(SubElement(i),j));
      }
    }
  }

  this->Print(cout);
  for (i=0;i<NSubEl;i++){
    TPZGeoEl *subel = SubElement(i);
    if (subel) cout << "Subel " << i << " OK " << endl ; //subel->Print(cout);
    else  cout << "Erro no subelement0 " << i << endl;
  }

  cout << "Chegou aqui ..." << endl;
  this->SetSubElementConnectivities();
}

template<class TShape, class TGeo> int TPZGeoElRefPattern<TShape,TGeo>::FatherSide(int side, int son){
  return fRefPattern->FatherSide(side,son);

}

template<class TShape, class TGeo> 
void TPZGeoElRefPattern<TShape,TGeo>::MidSideNodeIndices(int side,TPZVec<int> &indices){
  if(!HasSubElement() || !fRefPattern || side < this->NCornerNodes()) {
    indices.Resize(0);
    return;
  }
  TPZManVector<TPZGeoElSide> gelsides;
  fRefPattern->SidePartition(gelsides,side);
  int nsub = gelsides.NElements();
  indices.Resize(nsub);
  int is;
  int counter = 0;
  for(is=0; is<nsub; is++) {
    if(gelsides[is].Side() >= gelsides[is].Element()->NCornerNodes()) continue;
    TPZGeoEl *subel = SubElement(gelsides[is].Element()->Id()-1);
    indices[counter] = subel->NodeIndex(gelsides[is].Side());
    counter++;
  }
  indices.Resize(counter);  
}

/** Defines the element refinement pattern  */
template<class TShape, class TGeo> 
void TPZGeoElRefPattern<TShape,TGeo>::SetRefPattern (TPZRefPattern *refpat){
#ifdef HUGE_DEBUG
  if (!refpat) {
    PZError << "Error trying to set a null refinement pattern objetct" << endl;
    return;
  }
#endif
  this->Mesh()->InsertRefPattern(refpat);
  //para nao manter copias, caso exista uso o existente;
  int eltype = refpat->Element(0)->Type();
  string refname = refpat->GetName();
  refpat = this->Mesh()->GetRefPattern(eltype,refname);
  fRefPattern = refpat;
  int i;
  int nsubel = refpat->NSubElements();
  fSubEl.Resize(nsubel);
  for(i=0;i<nsubel;i++) fSubEl[i] = -1;
}

#endif

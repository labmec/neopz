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

    TPZVec<TPZGeoEl *> fSubEl;
    TPZRefPattern *fRefPattern;
  
public: 
  TPZGeoElRefPattern();
  ~TPZGeoElRefPattern();
  TPZGeoElRefPattern(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat);
  TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,TPZRefPattern *refpat);
  TPZGeoElRefPattern(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index,TPZRefPattern *refpat);
  void Initialize(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index,TPZRefPattern *refpat);

  /** return 1 if the element has subelements along side*/
  int HasSubElement() {return fSubEl[0]!=0;}

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
  /** Defines the element refinement pattern  */
  void SetRefPattern (TPZRefPattern *refpat);

};

#endif

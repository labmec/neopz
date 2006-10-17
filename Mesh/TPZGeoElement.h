
// -*- c++ -*-

// $Id: TPZGeoElement.h,v 1.15 2006-10-17 01:53:26 phil Exp $

#ifndef TPZGEOELEMENTH
#define TPZGEOELEMENTH

//#include "pzgeoel.h"
#include "pzgeoelrefless.h"

class TPZGeoElSide;
class TPZCompMesh;
class TPZCompEl;
template<class T,int N>
class TPZStack;

/// Implements a generic geometric element with a uniform refinement pattern
template <class TShape, class TGeo, class TRef>
class TPZGeoElement : public TPZGeoElRefLess<TShape,TGeo> {

  int fSubEl[TRef::NSubEl];
//int fNodeIndexes[TGeo::NNodes];
//TPZGeoElSide fNeighbours[TShape::NSides];
public:
//static TPZCompEl *(*fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index);
  //  static int fTest;

public:

  TPZGeoElement();
  TPZGeoElement(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh);
  TPZGeoElement(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh);
  TPZGeoElement(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index);

  /** Copy constructor */
  TPZGeoElement(TPZGeoMesh &DestMesh, const TPZGeoElement &cp);
  
  void Initialize(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index);
  //  TPZGeoElement( int* nodeindices, int matind, TPZGeoMesh& mesh );
  //  TPZGeoElement( int* nodeindices, int matind, TPZGeoMesh& mesh, int& index );

  ~TPZGeoElement(){};

  virtual int ClassId() const;
  
  virtual void Read(TPZStream &str, void *context);

  virtual void Write(TPZStream &str, int withclassid);
  
  virtual TPZGeoEl * Clone(TPZGeoMesh &DestMesh) const;

  /** return 1 if the element has subelements along side*/
  int HasSubElement() {return fSubEl[0]!=-1;}


  void SetSubElement(int id, TPZGeoEl *el);


  /**volume of the master element*/
  REAL RefElVolume();


  /**returns the midside node index along a side of the element*/
  void MidSideNodeIndex(int side,int &index);


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


  virtual int FatherSide(int side, int son) {
    return TRef::FatherSide(side,son);
  }

  /**divides the element and puts the resulting elements in the vector*/
  virtual void Divide(TPZVec<TPZGeoEl *> &pv);

  virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel);

};


#endif 



// -*- c++ -*-

// $Id: TPZGeoElement.h,v 1.11 2004-04-22 13:13:52 phil Exp $

#ifndef TPZGEOELEMENTH
#define TPZGEOELEMENTH

//#include "pzgeoel.h"
#include "pzgeoelrefless.h"

class TPZGeoElSide;
class TPZCompMesh;
class TPZCompEl;
template<class T,int N>
class TPZStack;


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

  void Initialize(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index);
  //  TPZGeoElement( int* nodeindices, int matind, TPZGeoMesh& mesh );
  //  TPZGeoElement( int* nodeindices, int matind, TPZGeoMesh& mesh, int& index );

  ~TPZGeoElement(){};

//  TPZCompEl *CreateCompEl(TPZCompMesh &mesh,int &index);

//  static void SetCreateFunction(TPZCompEl *(*f)(TPZGeoEl *el,TPZCompMesh &mesh,int &index));


  /** return 1 if the element has subelements along side*/
  int HasSubElement() {return fSubEl[0]!=-1;}

  /** 
   * returns a pointer to the neighbour and the neighbourside 
   * along side of the current element
   */
//  TPZGeoElSide Neighbour(int side) { return fNeighbours[side]; }

//  int NodeIndex(int node);

  /**fill in the data structure for the neighbouring information*/
//  void SetNeighbour(int side,const TPZGeoElSide &neighbour)
//								{ fNeighbours[side]=neighbour; }

//  int SideNodeIndex(int side,int node);

//  int SideNodeLocIndex(int side,int node);

  /**flags the side as defined, this means no neighbouring element
   * was found*/
//  void SetSideDefined(int side) { fNeighbours[side] = TPZGeoElSide(this,side); }

  void SetSubElement(int id, TPZGeoEl *el);

  /**
   * Creates an integration rule for the topology of the corresponding side
   * and able to integrate a polynom of order exactly
   */
//  TPZIntPoints * CreateSideIntegrationRule(int side, int order);

  /**
   * returns the type of the element acording to the definition in pzeltype.h
   */
//  int Type() {
//    return TGeo::Type();
//  }

  /**return the number of nodes of the element*/
//  int NNodes();

  /**return the number of corner nodes of the element*/
//  int NCornerNodes();

  /**return the number of connectivities of the element*/
//  int NSides();

 /**
  * returns the local node number of the node "node" along side "side"
  */
//  int SideNodeLocId(int side, int node);

  /**volume of the master element*/
  REAL RefElVolume();

  /**return the number of nodes for a particular side*/
//  int NSideNodes(int side);

  /**returns the midside node index along a side of the element*/
  void MidSideNodeIndex(int side,int &index);

  /**returns 1 if the side has not been defined by buildconnectivity
     After construction the side is undefined. The buildconnectivity method
     loops over all elements and tries to identify neighbours along their
     uninitialized sides*/
//  int SideIsUndefined(int side);

  /**
   * return the number of subelements of the element independent of the
   * fact hether the element has already been refined or not 
   */
  int NSubElements();

  /**return the number of subelements of the same dimension of the element at the side*/
//  int NSideSubElements(int side);

  /**
  * return the number of subelements as returned by GetSubElements2(side)
  */
  int NSideSubElements2(int side);

  /**
   * method which creates a computational boundary condition element based
   * on the current geometric element, a side and a boundary condition number
   */
//  TPZGeoEl *CreateBCGeoEl(int side, int bc);

  /**initializes the node i of the element*/
//  void SetNodeIndex(int i,int nodeindex);

  /**
   * compute the transformation between the master element space of one side
   * of an element to the master element space of a higher dimension side
   */
//  TPZTransform SideToSideTransform(int sidefrom,int sideto);

  /**returns a pointer to the subelement is*/
  TPZGeoEl *SubElement(int is);

  /**return a pointer and a side of the subelement of the element at the side
     and the indicated position. position = 0 indicate first subelement, ...*/
  TPZGeoElSide SideSubElement(int side,int position);

  /**return the dimension of side*/
//  int SideDimension(int side);

  /**Returns the dimension of the element*/
//  virtual int Dimension();

//  TPZGeoElSide HigherDimensionSides(int side,int targetdimension);
//  void AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides);
//  void LowerDimensionSides(int side,TPZStack<int> &smallsides);

  /**accumulates the transformation of the jacobian which maps the current
     master element space into the space of the master element of the father*/
//  void BuildTransform(int side, TPZGeoEl *father,TPZTransform &t);
//  TPZTransform BuildTransform2(int side, TPZGeoEl *father,TPZTransform &t);

  /**return the Jacobian matrix at the point*/
//  void Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix &jac,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

  /**return the coordinate in real space of the point coordinate in the master element space*/
//  void X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result);

  /**return the normal vector at the position loc in the master element space of side*/
//  void NormalVector(int side,TPZVec<REAL> &loc,TPZVec<REAL> &normal,
//									TPZFMatrix &axes,TPZFMatrix &jac);

  TPZTransform GetTransform(int side,int son);

  /**
   * It returns the coordinates from the center of the side of the element
   */
//  virtual void CenterPoint(int side, TPZVec<REAL> &masscent);

//  virtual TPZGeoElSide Father2(int side);

  virtual int FatherSide(int side, int son) {
    return TRef::FatherSide(side,son);
  }

  /**divides the element and puts the resulting elements in the vector*/
  virtual void Divide(TPZVec<TPZGeoEl *> &pv);

  virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel);

};

//template<class TShape, class TGeo, class TRef>
//inline
//TPZCompEl *TPZGeoElement<TShape,TGeo,TRef>::CreateCompEl(TPZCompMesh &mesh,int &index){
//  return fp(this,mesh,index);
//}
//
//template<class TShape, class TGeo, class TRef>
//inline
//void TPZGeoElement<TShape,TGeo,TRef>::SetCreateFunction(TPZCompEl *(*f)(TPZGeoEl *el,TPZCompMesh &mesh,int &index)){
//  fp = f;
//}

#endif 


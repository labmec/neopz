//$Id: pzgeoelside.h,v 1.11 2004-03-03 23:15:46 cesar Exp $

#ifndef PZGEOELSIDEH
#define PZGEOELSIDEH

// -*- c++ -*-

/*******       TPZGeoElSide       *******/

//#include "pzcompel.h"	// Added by ClassView
class TPZGeoEl;
class TPZTransform;
class TPZCompElSide;

#include "pzvec.h"
#include "pzstack.h"

class TPZGeoElSide {

  TPZGeoEl *fGeoEl;
  int fSide;
 public:
	 int NSubElements2();
	 void GetSubElements2(TPZStack<TPZGeoElSide> &subelements);
	 TPZGeoElSide StrictFather();
	 TPZGeoElSide Father2();
	 TPZCompElSide LowerLevelCompElementList2(int onlyinterpolated);

	 /**
	  * Will return all elements of equal or higher level than than the current element
	  * All elements/sides have the same dimension
	  */
	 //void EqualorHigherCompElementList(TPZStack<TPZCompElSide> &celside, int onlyinterpolated, int removeduplicates);
	 /**
	  * Will return all elements of equal or higher level than than the current element
	  */
	 void EqualorHigherCompElementList2(TPZStack<TPZCompElSide> &celside, int onlyinterpolated, int removeduplicates);

  TPZGeoElSide(){ fGeoEl = 0; fSide  = -1;}
  //TPZGeoElSide(const TPZGeoElSide &gelside);
  TPZGeoElSide(TPZGeoEl *gel,int side){  fGeoEl = gel; fSide  = side;}
  TPZGeoEl *Element()const{return fGeoEl;}
  int Side() const {return fSide;}
  void SetSide(int side) { fSide = side; }
  int Exists() const {return (fGeoEl != 0 && fSide > -1);}
  TPZGeoElSide Neighbour() const;//return neighbour of the side fSide

  /**
   * returns the set of neighbours which can directly be accessed by the datastructure
   */
  void AllNeighbours(TPZStack<TPZGeoElSide> &allneigh);
  
  /**
   * returns the set of neighbours as computed by the intersection of neighbours along vertices
   */
  void ComputeNeighbours(TPZStack<TPZGeoElSide> &compneigh);
 
  int Id();
  int Dimension() const;

  int operator==(const TPZGeoElSide &other) const {
    return fGeoEl == other.fGeoEl && fSide == other.fSide;
  }
  int operator!=(const TPZGeoElSide &other) const {
    return fGeoEl != other.fGeoEl || fSide != other.fSide;
  }
  
  int operator<(const TPZGeoElSide &other) const {
	  return (fGeoEl < other.fGeoEl || (fGeoEl == other.fGeoEl && fSide < other.fSide));
  }
  
  int operator>(const TPZGeoElSide &other) const {
	  return (fGeoEl > other.fGeoEl || (fGeoEl == other.fGeoEl && fSide > other.fSide));
  }

  //void SideTransform(TPZGeoElSide neighbour,TPZTransform &t);
  /**set the element neighbour along side i equal to el
     the side of the neighbour which is connected to the current element
     is neighbourside*/
  //void SideTransform2(TPZGeoElSide neighbour,TPZTransform &t);

  /**
  Accumulates the transformations from the current element/side to the
  neighbour/side
  Third improved version
  */
  void SideTransform3(TPZGeoElSide neighbour,TPZTransform &t);

  void SetConnectivity(const TPZGeoElSide &neighbour) const;

  void RemoveConnectivity();

static void BuildConnectivities(TPZVec<TPZGeoElSide> &elvec, TPZVec<TPZGeoElSide> &neighvec);

/**fill in the data structure for the neighbouring information*/
  void SetNeighbour(const TPZGeoElSide &neighbour) const;

  /**returns a pointer to computational element referenced by a geometric
   element which is a son along side and has higher level than level
   if onlyinterpolated != 0 only elements which return IsInterpolated() != 0 will be put on the stack*/
//  void SmallConnect(int level,TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated);

  /**search elements of high dimension connected to the actual element*/
  //TPZCompElSide HigherDimensionConnected(int targetdimension,int onlyinterpolated);

  //TPZGeoElSide HigherDimensionSide(int targetdimension);

  /**return pointers to the subelements and their sides along side
   in the order determined by the referencenodes*/
  //void GetSubElement(TPZVec<int> &referencenodes,TPZVec<TPZGeoElSide> &subelements);

  /**return all subelements along the current side*/
  //void GetSubElements(TPZVec<TPZGeoElSide> &subelements);

  TPZTransform NeighbourSideTransform(TPZGeoElSide &neighbour);

/**compute the transformation between the master element space of one side
   of an element to the master element space of a higher dimension side
*/
  TPZTransform SideToSideTransform(TPZGeoElSide &higherdimensionside);
	
  ////////////////////////////////////////////////////////////////////////////////
  /**returns the number of fathers that can be followed*/
  //int Level();
  /**return a pointer to the elementside referenced by the geometric elementside*/
  TPZCompElSide Reference() const;
  /*return 1 if the element has subelements along side*/
  int HasSubElement();
	
  /*return the number of nodes for a particular side*/
  int NSideNodes();

  /**returns the index of the nodenum node of side*/
  int SideNodeIndex(int nodenum);

  /**returns 1 if neighbour is a neighbour of the element along side*/
  int NeighbourExists(const TPZGeoElSide &neighbour) const;

  /**Pushes all connected computational elements which have higher dimension than the current element/side
     if onlyinterpolated == 1 only elements which return IsInterpolated()== 1 will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
  */
void HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated);


  /**Returns all connected computational elements which have same dimension and level higher to the current element
     if onlyinterpolated == 1 only elements which return IsInterpolated()== 1 will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
  //void HigherLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
  /**Returns all connected computational elements which have level higher to the current element
     if onlyinterpolated == 1 only elements which return IsInterpolated()== 1 will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
  void HigherLevelCompElementList2(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
  /**Returns all connected computational elements to the current element
     if onlyinterpolated == 1 only elements which return IsInterpolated()== 1 will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
  void ConnectedCompElementList(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated, int removeduplicates);
  /**Returns all connected computational elements which have equal level to the current element
     This method will not put this on the stack
     if onlyinterpolated == 1 only elements which return IsInterpolated()== 1 will be put on the stack
     if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
  void EqualLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec,	int onlyinterpolated, int removeduplicates);
  //void Dim0EqualLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated, int removeduplicates);

  /**Returns all connected computational elements which have level lower to the current element
     if onlyinterpolated == 1 only elements which return IsInterpolated()== 1 will be put on the stack
     This method will only return element/sides of dimension > 0*/
  //TPZCompElSide LowerLevelCompElementList(int onlyinterpolated);
  /**Retorna o elem. computacional de nivel menor (elemento grande) ao qual estou restrito*/
  //TPZCompElSide TPZGeoElSide::Dim0LowerLevelCompElement(int onlyinterpolated);

  /** Jorge 13/01/2000 */
  /** Retorna todos os subelementos computacionais do atual elemento geometrico
   * ao longo do lado fSide
	 */
  //void GetCompSubElements(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated,int removeduplicates);

};

ostream  &operator << (ostream & out,const TPZGeoElSide &geoside);

#include "pzgeoel.h"

inline TPZGeoElSide TPZGeoElSide::Neighbour() const {
  return fGeoEl ? fGeoEl->Neighbour(fSide) : TPZGeoElSide();
}

inline void TPZGeoElSide::AllNeighbours(TPZStack<TPZGeoElSide> &allneigh) {
#ifndef NDEBUG
  if(! Exists() || ! this->Neighbour().Exists()) 
    {
      cout << "TPZGeoElSide AllNeighbours inconsistent\n";
      return;
    }
#endif
  TPZGeoElSide neigh = Neighbour();
  while(neigh != *this)
    {
      allneigh.Push(neigh);
      neigh = neigh.Neighbour();
    }
}

#endif

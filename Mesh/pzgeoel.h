//$Id: pzgeoel.h,v 1.11 2003-12-12 19:59:20 phil Exp $

// -*- c++ -*-

#ifndef GEOELEMHPP
#define GEOELEMHPP

#include <iostream>

using namespace std;

#include "pzerror.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pztrnsform.h"
#include "doxmesh.h"

class TPZGeoNode;
class TPZCompMesh;
class TPZCompEl;
class TPZFMatrix;
class TPZGeoMesh;
class TPZCompElSide;
class TPZGeoElSide;
class TPZIntPoints;

template<class T>
class TPZVec;
template<class T, int N>
class TPZStack;


/** 
 * TPZGeoEl is the common denominator for all geometric elements.
 * @ingroup geometry
 */
class TPZGeoEl { 	// header file for the element class

 private:

  /**pointer to the mesh to which the element belongs*/
  TPZGeoMesh *fMesh;
  /**traditional element number*/
  int		fId;
  /**material index*/
  int		fMatId;
  /**pointer to the element currently loaded*/
  TPZCompEl *fReference;

 protected:

  /**the element from which the element is a subelement*/
  TPZGeoEl *fFather;
  /**3x3 unit matrix to be copied to the axes if the geometric element
     does not have a particular orientation*/
  static TPZFMatrix gGlobalAxes;

public:

	/**
	* Creates an integration rule for the topology of the corresponding side
	* and able to integrate a polynom of order exactly
	*/
	virtual TPZIntPoints * CreateSideIntegrationRule(int side, int order) =0;

  /**
   * computes the values of the "num" one dimensional shapefunctions
   * and derivatives at point x, using lagrangian interpolation
   * @param x point at which the shapefunction is computed
   * @param num number of shapefunctions
   * @param phi values of the shapefunctions
   * @param dphi values of the derivatives of the shapefunctions
   */
  void Shape1d(double x,int num,TPZFMatrix &phi,TPZFMatrix &dphi);

  /**
   * computes the values of the "num" one dimensional shapefunctions
   * at point x, using lagrangian interpolation
   * @param x point at which the shapefunction is computed
   * @param num number of shapefunctions
   * @param phi values of the shapefunctions
   */
  void ShapePhi1d(double x,int num,TPZFMatrix &phi);

  /**constructor : 
   * @param Id is the number of the element 
   * @param materialindex is the material index
   * @param mesh is a pointer to the mesh to which the element belongs
   */
  TPZGeoEl(int id,int materialindex,TPZGeoMesh &mesh);
  /**This constructor generates a unique Id
   * @param materialindex is the material index
   * @param mesh is a pointer to the mesh to which the element belongs
   */
  TPZGeoEl(int materialindex,TPZGeoMesh &mesh);

  /**This constructor generates a unique Id
   * @param materialindex is the material index
   * @param mesh is a pointer to the mesh to which the element belongs
   * @param index index of the new element in the element vector
   */
  TPZGeoEl(int materialindex,TPZGeoMesh &mesh,int &index);

  /**
   * Copy constructor
   */
  TPZGeoEl(const TPZGeoEl &el) ;

  TPZGeoEl() {
    fId = -1;
    fMesh = 0;
    fMatId = 0;
    fReference = 0;
    fFather = 0;
  }

  virtual void Initialize(int materialindex, TPZGeoMesh &mesh, int &index);

  /**Destructor*/
  virtual ~TPZGeoEl() { }

  /**it removes the connectivities of the element*/
  void RemoveConnectivities();


  /** @name data access methods
   * methods which allow to access the internal data structure of the element
   */
  //@{
  /**return the mesh to which the element belongs*/
  TPZGeoMesh *Mesh() { return fMesh;}

  /**return the Id of the element*/
  int Id() { return fId; }

  /**return the number of nodes of the element*/
  virtual int NNodes() = 0;

  /**return the number of corner nodes of the element*/
  virtual int NCornerNodes() = 0;

  /**return a pointer to the ith node of the element*/
  TPZGeoNode* NodePtr(int i) {return &(fMesh->NodeVec()[NodeIndex(i)]); }

  /**return the index of the ith node
     the index is the location of the node in the nodevector of the mesh*/
  virtual int NodeIndex(int i) = 0;

    /**returns the face node ids */
 // virtual void NodeFaceIds(TPZVec<int> &ids,int face);

  /**return the material index of the element*/
  int MaterialId() { return fMatId; }

  /**return a pointer to the element referenced by the geometric element*/
  TPZCompEl *Reference() const { return fReference; }

  /**
   * returns the element type acording to pzeltype.h
   */
virtual int Type() =0;
// {
//    cout << "ElementType should never be called\n";
//    return -1;
//  }
  /**return the number of connectivities of the element*/
  virtual int NSides() = 0;

  /**return the number of nodes for a particular side*/
  virtual int NSideNodes(int side) = 0;

  /**returns the pointer to the nodenum node of side*/
  virtual TPZGeoNode *SideNodePtr(int side,int nodenum) {
    return &(fMesh->NodeVec()[SideNodeIndex(side,nodenum)]);
  }

  /**returns the midside node index along a side of the element*/
  virtual void MidSideNodeIndex(int side,int &index) = 0;

  /**returns the index of the nodenum node of side*/
  virtual int SideNodeIndex(int side,int nodenum) = 0;
  
  /**returns the local index of a node on a side*/
  virtual int SideNodeLocIndex(int side, int nodenum) = 0;

  /**returns 1 if the side has not been defined by buildconnectivity
     After construction the side is undefined. The buildconnectivity method
     loops over all elements and tries to identify neighbours along their
     uninitialized sides*/
  virtual int SideIsUndefined(int side) = 0;

  /**return the number of subelements of the element
  independent of the fact hether the element has already been refined or not */
  virtual int NSubElements() = 0;

  /**return the number of subelements of the same dimension of the element at the side*/
//  virtual int NSideSubElements(int side) =0;

  /**
  * return the number of subelements as returned by GetSubElements2(side)
  */
  virtual int NSideSubElements2(int side) = 0;

  /**return a pointer to the father*/
  TPZGeoEl *Father() { return fFather; }

  //@}

  /**method which creates a computational element based on the current
     geometric element*/
  virtual TPZCompEl *CreateCompEl(TPZCompMesh &cmesh,int &index) = 0;

  /**method which creates a computational boundary condition element
     based on the current geometric element, a side and a boundary condition number*/
  virtual TPZCompEl *CreateBCCompEl(int side, int bc, TPZCompMesh &cmesh);

  /** method which creates a geometric element on the side of an existing element */
  virtual TPZGeoEl *CreateBCGeoEl(int side, int bc) = 0;

  /**returns the side number which is connected to the SideNodes
     returns -1 if no side is found*/
  int WhichSide(TPZVec<int> &SideNodeIds);

  /**returns 1 if gel is a neighbour of the element along side*/
  int NeighbourExists(int side,const TPZGeoElSide &gel);

  /**Sets the material index of the element*/
  void SetMaterialId(int id) { fMatId = id;}

  /**initializes the node i of the element*/
  virtual void SetNodeIndex(int i,int nodeindex) = 0;

  /**flags the side as defined, this means no neighbouring element was found*/
  virtual void SetSideDefined(int side) = 0;

  /**returns a pointer to the neighbour and the neighbourside along side of the current element*/
  virtual TPZGeoElSide Neighbour(int side) = 0;

  /**fill in the data structure for the neighbouring information*/
  virtual void SetNeighbour(int side,const TPZGeoElSide &neighbour) = 0;

  /**Print all relevant data of the element to cout*/
virtual  void Print(ostream & out = cout);

  /**Make the current element reference to the computational element*/
  void SetReference(TPZCompEl *elp) { fReference = elp; }

  /**Set the subelement of index i*/
  virtual void SetSubElement(int id, TPZGeoEl* gel) = 0;

  /**Initializes the external connectivities of the subelements*/
  virtual void SetSubElementConnectivities();

  /**reset the element referenced by the geometric element to NULL*/
  void ResetReference() { fReference = NULL; }

  /**equivalent to Print*/
  friend ostream& operator<<(ostream &s,TPZGeoEl &el);

  /**divides the element and puts the resulting elements in the vector*/
  virtual void Divide(TPZVec<TPZGeoEl *> &pv);

  /**return the middlenode side if there is one*/
  //virtual TPZGeoNode *MiddleSideNode(int side);

  /**return pointers to the subelements and their sides along side
     in the order determined by the referencenodes*/
  //virtual void GetSubElement(int side,TPZVec<int> &referencenodeindex,
//			     TPZVec<TPZGeoElSide> &subelements);

  /**return 1 if the element has subelements along side*/
  virtual int HasSubElement() = 0;

  /**compute the transformation for a point on the master element to a point
     in the master element of the neighbour*/
  void SideTransform(int side,TPZGeoElSide neighbour,TPZTransform &t);

  /**compute the transformation between the master element space of one side
     of an element to the master element space of a higher dimension side*/
  virtual TPZTransform SideToSideTransform(int sidefrom,int sideto)= 0;

  /**get the transform id the face to face*/
	int GetTransformId2dQ(TPZVec<int> &idfrom,TPZVec<int> &idto);

  /**get the transform id the face to face*/
	int GetTransformId2dT(TPZVec<int> &idfrom,TPZVec<int> &idto);

virtual	TPZTransform GetTransform(int side,int son) = 0;

  /**returns a pointer to computational element referenced by a geometric
     element which is a son along side and has higher level than level
     if onlyinterpolated != 0 only elements which return IsInterpolated() != 0 will be put on the stack*/
  //void SmallConnect(int side,int level,TPZStack<TPZCompElSide> &elvec,
//		    int onlyinterpolated);

  /**Sets the father element*/
  void SetFather(TPZGeoEl *father);

  /**returns a pointer to the subelement is*/
  virtual TPZGeoEl *SubElement(int is) = 0;

  /**return a pointer and a side of the subelement of the element at the side
     and the indicated position. position = 0 indicate first subelement, ...*/
//  virtual TPZGeoElSide SideSubElement(int side,int position) = 0;

  /**returns the element from which the current element is a subelement
     climbing up along side
     This method may return NULL even if the element has a father*/
//  virtual TPZGeoElSide Father(int side);

  /**returns the number of fathers that can be followed*/
  int Level();

  /**return the dimension of side*/
  virtual int SideDimension(int side) = 0;

  /**Returns the dimension of the element*/
  virtual int Dimension() =0;

  /** */
  virtual TPZGeoElSide HigherDimensionSides(int side,int targetdimension);//SÓ PARA TESTAR CONTINUIDADE - APAGAR DEPOIS
  virtual void AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides) = 0;
  virtual void LowerDimensionSides(int side,TPZStack<int> &smallsides) = 0;

  /**accumulates the transformation of the jacobian which maps the current
     master element space into the space of the master element of the father*/
//  virtual void BuildTransform(int side, TPZGeoEl *father,TPZTransform &t) = 0;

  /**return the Jacobian matrix at the point*/
  virtual void Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix &jac,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv)=0;

  /**return the coordinate in real space of the point coordinate in the master element space*/
  virtual void X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result)=0;

  /**return the normal vector at the position loc in the master element space of side*/
//  virtual void NormalVector(int side,TPZVec<REAL> &loc,TPZVec<REAL> &normal,
//			    TPZFMatrix &axes,TPZFMatrix &jac) = 0;
  //para testar continuidade
  int ElementExists(TPZGeoEl *elem,int id);

  /** 
   * @name reftopology
   * Methods which will implement the declaration of a refinemnt topology
   */

  //@{

  /**
  * returns the father/side of the father which contains the side of the
  * sub element
  **/
  virtual TPZGeoElSide Father2(int side);

  virtual int FatherSide(int side, int son);

  /**
  * returns the transformation which maps the parameter side of the element/side
  * into the parameter space of the father element/side
  **/
virtual TPZTransform BuildTransform2(int side, TPZGeoEl *father, TPZTransform &t);

	/**
	* This method will return a partition of the side of the current element
	* as the union of sub elements/side which are put in the stack
	**/

	/**returns the side number which is connected to the point pt
     *returns -1 if no side is found
     */
     int WhichSide(TPZVec<REAL> &pt);

  /**
   * It returns the coordinates from the center of the side of the element
   */
  virtual void CenterPoint(int side, TPZVec<REAL> &masscent) = 0;


virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel);

void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel, int dimension);

/*return the son number of the sub element gel*/
int WhichSubel();

//checa a estrutura de dados de Father() e GetSubelement2()
void CheckSubelDataStructure();

     /*testa as transformacoes entre lados de pais e filhos*/
/* int main(TPZGeoEl *gel,int type); */
  //@}
	
  /**Jorge 17/7/99*/
  /** Return the measure of the geometrical element - Area */
//  virtual REAL Mesure(int dim) { return fMesure; }
	/** 
	 * Return into the center a especial point of the geometrical element
	 * If 1-d => middle point,
	 * If 2-d => barycenter, orthocenter, etc
	 */
//  virtual void Center(TPZVec<REAL> &center) { center[0] = 0.; }

void ComputeXInverse(TPZVec<REAL> &XD, TPZVec<REAL> &ksi);

TPZTransform ComputeParamTrans(TPZGeoEl *fat,int fatside, int sideson);


  /**return the volume of the element*/
  REAL Volume();

  /**volume of the master element*/
  virtual REAL RefElVolume() = 0;

  /**returns the area from the face*/
  virtual REAL SideArea(int side);

  /**return the área from a quadrilateral face*/
  static  REAL QuadArea(TPZVec<TPZGeoNode *> &nodes);

  /**return the área from the triangular face*/
  static REAL TriangleArea(TPZVec<TPZGeoNode *> &nodes);

  virtual REAL ElementRadius();//TPZGeoEl

  static REAL Distance(TPZVec<REAL> &centel,TPZVec<REAL> &centface);

 protected:
//  REAL fMesure;
 private:
  /**
   * To be used after the buid connectivity. If some neighbour isn't initialized\
   * It will be initialized as the TPZGeoElSide (this, thisside)
   */
  void InitializeNeighbours();
};

inline void TPZGeoEl::Divide(TPZVec<TPZGeoEl *> &) {
  PZError << "TPZGeoEl::Divide is called.\n";
}
//inline void TPZGeoEl::NodeFaceIds(TPZVec<int> &ids,int face) {
//	cout << "TPZGeoEl::NodeFaceIds is called." << endl;
//}




#include "pzgeoelside.h"
#include "pzgeoelbc.h"

inline TPZGeoElSide TPZGeoEl::HigherDimensionSides(int side,int targetdimension){//SÓ PARA TESTAR CONTINUIDADE - APAGAR DEPOIS
  cout << "TPZGeoEl::HigherDimensionSides is called." << endl;
  return TPZGeoElSide();
}
#endif

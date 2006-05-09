// -*- c++ -*- 

//$Id: TPZCompElDisc.h,v 1.44 2006-05-09 17:48:17 tiago Exp $

////////////////////////////////////////////////////////////////////////////////
// Discontinou Element
////////////////////////////////////////////////////////////////////////////////
#ifndef ELCOMPDISCHPP
#define ELCOMPDISCHPP

#include <iostream>
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzreal.h"
#include "TPZShapeDisc.h"

struct TPZElementMatrix;
class TPZFMatrix;
class TPZBndCond;
class TPZConnect;
class TPZMaterial;
class TPZGeoEl;
class TPZMatrix;
class TPZCompMesh;
class TPZGeoElC3d;
class TPZGeoElT3d;
class TPZGeoElPi3d;
class TPZGeoElPr3d;
class TPZGeoElT2d;
class TPZGeoElQ2d;
class TPZGeoEl1d;
class TPZGeoElPoint;

/// This class implements a discontinuous element (for use with discontinuous Galerkin)
/**
@ingroup CompElement
*/
class TPZCompElDisc : public TPZCompEl{

protected:

  /**
   * Shape function type used by the element
   */
  pzshape::TPZShapeDisc::MShapeType fShapefunctionType;

  /**
   * it preserves index of connect associated to the element 
   */
  int fConnectIndex;

  /**
   * Normalizing constant for shape functions 
   */
  REAL fConstC;

  protected:
  /**
   * Material object of this element
   */
  TPZMaterial *fMaterial;

  /**
   * it keeps the interior point coordinations of the element 
   */
  TPZManVector<REAL,3> fCenterPoint;

  /**
   * it creates new conect that it associates the degrees of freedom of the
   * element and returns its index 
   */
  virtual int CreateMidSideConnect();

 public:

  int GetMaterial( const TPZGeoElSide& gside );

  static TPZCompEl *CreateDisc(TPZGeoEl *geo, TPZCompMesh &mesh, int &index);

  /**return the geometric element to which this element references*/
//  TPZGeoEl *Reference() const { return fReference;}

  /**set the geometric element to which this element references*/
//  void SetReference(TPZGeoEl *ref) {fReference = ref;}

  /**
   * Sets the orthogonal function which will be used throughout the program.
   * @param orthogonal pointer to a function which will be used to generate the shape functions
   */
  static void SetOrthogonalFunction(void (*orthogonal)(REAL C, REAL x0, REAL x,int degree, TPZFMatrix & phi, TPZFMatrix & dphi, int n)){
    pzshape::TPZShapeDisc::fOrthogonal = orthogonal;
  }

  /**
   * Set the inner radius value.
   */
  virtual void SetInnerRadius(REAL InnerRadius) {PZError << "TPZCompElDisc::SetInnerRadius - This method should never be called because the inner" << std::endl 
							 << "radius is not stored in TPZCompElDisc. It is stored in TPZAgglomerateElement." << std::endl;}

  /**
   * Set element's number of interfaces.
   */
  virtual void SetNInterfaces(int nfaces) {PZError << "TPZCompElDisc::SetNFaces - This method should never be called because the number of interfaces" << std::endl 
				      << "is not stored in TPZCompElDisc. It is only stored by TPZAgglomerateElement." << std::endl;}

  /**
   * Retunrs the number of interfaces.
   */
  virtual int NInterfaces();

  /**
   * Returns the inner radius value.
   */
  virtual REAL InnerRadius() {return this->Reference()->ElementRadius();}

  TPZCompElDisc();

  TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int &index);//original
  TPZCompElDisc(TPZCompMesh &mesh,int &index);//construtor do aglomerado

  TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy);
  TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy,int &index);

virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    return new TPZCompElDisc(mesh,*this);
  }

virtual TPZCompEl *Clone(TPZCompMesh &mesh,int &index) const {
    return new TPZCompElDisc(mesh,*this,index);
  }

  ~TPZCompElDisc() {
    if(Reference()->Reference() == this) Reference()->ResetReference();
   }

  /**
   * Divide the computational element
   */
  void Divide(int index, TPZVec<int> &subindex, int interpolate = 0);

  /**
   * CalcStiff computes the element stiffness matrix and right hand side
   * @param ek element matrix
   * @param ef element right hand side
   */
  virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * value of the bases and derivatives of the element deformed in point X 
   */
  void Shape(TPZVec<REAL> &X, TPZFMatrix &phi, TPZFMatrix &dphi);

  /**
   * Type of the element 
   */
  virtual MElementType Type() {return EDiscontinuous;}

  /**
   * it returns the material object 
   */
  virtual TPZMaterial *Material() const {return fMaterial;}

  /**
   *
   */
  void InterpolateSolution(TPZCompElDisc &coarsel);

  /**
   * it returns the constant that normalizes the bases of the element 
   */
  REAL ConstC(){return fConstC;}

  void SetConstC(REAL c){fConstC = c;}

  /**
   */
  void InternalPoint(TPZVec<REAL> &point);

  /**set the material of the element*/
  virtual void SetMaterial(TPZMaterial *mat) {fMaterial = mat;}

  /**
   * it prints the features of the element 
   */
  virtual void Print(std::ostream & out = std::cout);

  /**
   * it returns the degree of interpolation of the element 
   */
  virtual int Degree(){
    if (fConnectIndex < 0) return -1;
    return this->Connect(0).Order() ;
  }

  /**
   * it assigns the degree of the element 
   */
  virtual void SetDegree(int degree);// {fDegree = degree;}

  int NConnects();

  /**
   * amount of vertices of the element 
   */
  int NCornerConnects() { return Reference()->NNodes();}

  /**
   * it returns dimension from the element
   */
  int Dimension() const { return Reference()->Dimension();}

  /**
   * it calculates the normalizing constant of the bases of the element 
   */
  virtual REAL NormalizeConst();

  /**
   * it returns the connect index from the element 
   */
  int ConnectIndex(int side = 0);
  void  SetConnectIndex(int /*inode*/, int index) {fConnectIndex = index;}

  /**
   * it returns the shapes number of the element 
   */
  int  NShapeF();

  REAL CenterPoint(int index) {return fCenterPoint[index];}

  virtual void CenterPoint(TPZVec<REAL> &center);
  
  void SetCenterPoint(int i,REAL x){fCenterPoint[i] = x;}

  REAL SizeOfElement();

  /** 
   * Returns the volume of the geometric element associated.
   */
  virtual  REAL VolumeOfEl() { return Reference()->Volume(); }

  /**
   * Creates corresponding graphical element(s) if the dimension matches
   * graphical elements are used to generate output files
   * @param graphmesh graphical mesh where the element will be created
   * @param dimension target dimension of the graphical element
   */
  void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension);

  /**
  * \brief Computes the solution in function of a point in cartesian space
  */
  void Solution(TPZVec<REAL> &x,TPZVec<REAL> &uh);

  /**
   * Calculates the solution - sol - for the variable var
   * at point qsi, where qsi is expressed in terms of the
   * master element coordinates
   * @param qsi master element coordinate
   * @param var variable name
   * @param sol vetor for the solution
   */
  virtual void Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol);

  virtual void AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight);

  /** accumulate the vertices of the agglomerated elements */
  virtual void AccumulateVertices(TPZStack<TPZGeoNode *> &nodes);

  int NSides();

  void CalcResidual(TPZElementMatrix &ef);

  void BuildTransferMatrix(TPZCompElDisc &coarsel, TPZTransfer &transfer);

  void EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
		     TPZVec< REAL> &errors, TPZBlock * /*flux */ );
  
  /**
  * returns the unique identifier for reading/writing objects to streams
  */
  virtual int ClassId() const;
  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);
  
  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);


};

inline TPZCompEl *TPZCompElDisc::CreateDisc(TPZGeoEl *geo, TPZCompMesh &mesh, int &index) {
  return new TPZCompElDisc(mesh,geo,index);
}
//Exemplo do quadril�ero: 
//acessar com -> TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
#endif

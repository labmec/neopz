// -*- c++ -*-

//$Id: TPZCompElDisc.h,v 1.57 2007-05-01 20:27:18 phil Exp $

////////////////////////////////////////////////////////////////////////////////
// Discontinous Elements
////////////////////////////////////////////////////////////////////////////////
#ifndef ELCOMPDISCHPP
#define ELCOMPDISCHPP

#include <iostream>
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "pzgeoel.h"
#include "pzreal.h"
#include "TPZShapeDisc.h"
#include "tpzautopointer.h"
#include "pzmaterial.h"
#include "pzquad.h"

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
class TPZCompElDisc : public TPZInterpolationSpace{

protected:

  TPZIntPoints * fIntRule;

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

  TPZCompElDisc();

  TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int &index);//original
  TPZCompElDisc(TPZCompMesh &mesh,int &index);//construtor do aglomerado

  TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy);


  /**
   * Creates a clone of the given element in a pathc mesh
   * @param mesh patch mesh
   * @param copy element to be copied
   * @param gl2lcConMap map between the connect indexes in original and patch mesh
   * @param gl2lcElMap map between the element indexes in original an patch mesh
   */
  TPZCompElDisc(TPZCompMesh &mesh,
                const TPZCompElDisc &copy,
                std::map<int,int> &gl2lcConMap,
                std::map<int,int> &gl2lcElMap);


  TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy,int &index);

  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    return new TPZCompElDisc(mesh,*this);
  }

  virtual TPZCompEl *Clone(TPZCompMesh &mesh,int &index) const {
    return new TPZCompElDisc(mesh,*this,index);
  }

  /**
   * @see class TPZCompEl
   */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                  std::map<int,int> & gl2lcConMap,
                                  std::map<int,int> & gl2lcElMap) const {
    return new TPZCompElDisc(mesh,*this,gl2lcConMap,gl2lcElMap);
  }

  ~TPZCompElDisc();

  /**
   * Divide the computational element
   */
  void Divide(int index, TPZVec<int> &subindex, int interpolate = 0);

  /**
   * value of the bases and derivatives of the element deformed in point X
   */
  virtual void ShapeX(TPZVec<REAL> &X, TPZFMatrix &phi, TPZFMatrix &dphi);

  /**computes the shape function set at the point x. This method uses the order of interpolation
   * of the element along the sides to compute the number of shapefunctions
   * @param qsi point in master element coordinates
   * @param phi vector of values of shapefunctions, dimension (numshape,1)
   * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
   */
  virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix &phi,TPZFMatrix &dphi){
    TPZManVector<REAL,4> x(3);
    Reference()->X(qsi,x);
    ShapeX(x,phi,dphi);
  }

  virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv,
                            TPZFMatrix &phi, TPZFMatrix &dphix);

  /**returns a reference to an integration rule suitable for integrating
     the interior of the element */
  virtual TPZIntPoints &GetIntegrationRule();

  /**
   * Type of the element
   */
  virtual MElementType Type() {return EDiscontinuous;}

  /**
   * it returns the constant that normalizes the bases of the element
   */
  REAL ConstC(){return fConstC;}

  void SetConstC(REAL c){fConstC = c;}

  /**
   */
  void InternalPoint(TPZVec<REAL> &point);

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
  virtual int NShapeF();

  /**returns the number of shapefunctions associated with a connect*/
  virtual int NConnectShapeF(int inod);

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
  void SolutionX(TPZVec<REAL> &x,TPZVec<REAL> &uh);

  /**
   * Computes solution and its derivatives in the local coordinate qsi.
   * @param qsi master element coordinate
   * @param sol finite element solution
   * @param dsol solution derivatives
   */
  virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix & axes);

 /**
  * Computes solution and its derivatives in local coordinate qsi
  * @param qsi master element coordinate
  * @param phi matrix containing shape functions compute in qsi point
  * @param dphix matrix containing the derivatives of shape functions with respect of global coordinates: D[phi,x], D[phi,y], D[phi,z]
  * @param sol finite element solution
  * @param dsol solution derivatives
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                               const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol);

 /**
   * Computes solution and its derivatives in the local coordinate qsi.
   * This method will function for both volumetric and interface elements
   * @param qsi master element coordinate of the interface element
   * @param sol finite element solution
   * @param dsol solution derivatives
   * @param axes axes associated with the derivative of the solution
   * @param leftsol finite element solution
   * @param dleftsol solution derivatives
   * @param leftaxes axes associated with the left solution
   * @param rightsol finite element solution
   * @param drightsol solution derivatives
   * @param rightaxes axes associated with the right solution
  */
virtual void ComputeSolution(TPZVec<REAL> &qsi,
                             TPZVec<REAL> &normal,
                             TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                             TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes);

  virtual void AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight);

  /** accumulate the vertices of the agglomerated elements */
  virtual void AccumulateVertices(TPZStack<TPZGeoNode *> &nodes);

  int NSides();

  void BuildTransferMatrix(TPZCompElDisc &coarsel, TPZTransfer &transfer);

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

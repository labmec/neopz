////////////////////////////////////////////////////////////////////////////////
// Discontinou Element
////////////////////////////////////////////////////////////////////////////////
#ifndef ELCC3DDISCHPP
#define ELCC3DDISCHPP

#include <iostream>
using namespace std;
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzreal.h"

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


class TPZCompElDisc : public TPZCompEl{	// header file for the computational element class

  /**
   * Geometric element to which this element refers
   */
  TPZGeoEl *fReference;

  /**
   * Interpolation order of the discontinous element 
   */
  int fDegree;

  /**
   * it preserves index of connect associated to the element 
   */
  int fConnectIndex;

  /**
   * it keeps the interior point coordinations of the element 
   */
  TPZVec<REAL> fCenterPoint;

  /**
   * Normalizing constant for shape functions 
   */
  REAL fConstC;

  /**
   * Material object of this element
   */
  TPZMaterial *fMaterial;

protected:

  /**
   * it creates new conect that it associates the degrees of freedom of the
   * element and returns its index 
   */
  virtual int CreateMidSideConnect();

  /**
   * it returns the normal to the face from the element 
   */
  //  virtual void NormalVector(int side,TPZVec<REAL> &int_point,
  //		    TPZVec<REAL> &normal,TPZFMatrix &axes, TPZFMatrix &norm_r3);

 public:

  static TPZCompEl *CreateC3Disc(TPZGeoElC3d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreateT3Disc(TPZGeoElT3d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreatePi3Disc(TPZGeoElPi3d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreatePr3Disc(TPZGeoElPr3d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreateT2Disc(TPZGeoElT2d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreateQ2Disc(TPZGeoElQ2d *geo, TPZCompMesh &mesh, int &index);

  /**return the geometric element to which this element references*/
  virtual TPZGeoEl *Reference() { return fReference;}

  /**
   * default degree of imterpolation
   */
  static int gDegree;

   TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int &index);
/*    TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElC3d *ref,int &index); */
/*    TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElT3d *ref,int &index); */
/*    TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElPi3d *ref,int &index); */
/*    TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElPr3d *ref,int &index); */
/*    TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElT2d *ref,int &index); */
/*    TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElQ2d *ref,int &index); */
  ~TPZCompElDisc() {
      if(Reference()) {
         Reference()->ResetReference();
      }
   }

  /**
   * Divide the computational element
   */
  void Divide(int index,TPZVec<int> &subs,int degree = 0);

  /**
   * CalcStiff computes the element stiffness matrix and right hand side
   * param ek element matrix
   * param ef element right hand side
   */
  void CalcStiffDisc(TPZFMatrix &ek, TPZFMatrix &ef);

  /**
   * value of the bases and derivatives of the element deformed in point X 
   */
  void Shape(TPZVec<REAL> X, TPZFMatrix &phi, TPZFMatrix &dphi);

  /**
   * Type of the element 
   */
  MElementType Type() { return EDiscontinous; }

  /**
   * it returns the material object 
   */
  virtual TPZMaterial *Material() {return fMaterial;}

  /**
   *
   */
  void InterpolateSolution(TPZCompElDisc &coarsel);

  /**
   * it returns the constant that normalizes the bases of the element 
   */
  double ConstC(){return fConstC;}

  /**
   */
  void InternalPoint(TPZVec<REAL> &point);

  /**set the material of the element*/
  virtual void SetMaterial(TPZMaterial *mat) {fMaterial = mat;}

  /**
   * it prints the features of the element 
   */
  virtual void Print(ostream & out = cout);

  /**
   * it returns the degree of interpolation of the element 
   */
  virtual int Degree() {return fDegree;}

  /**
   * it assigns the degree of the element 
   */
  virtual void SetDegree(int degree) {fDegree = degree;}

  int NConnects() { return 1;}

  /**
   * amount of vertices of the element 
   */
  int NCornerConnects() { return Reference()->NNodes();}

  /**
   * it returns dimension from the element 
   */
  int Dimension() { return Reference()->Dimension();}

  /**
   * it calculates the normalizing constant of the bases of the element 
   */
  REAL NormalizeConst();

  /**
   * it returns the connect index from the element 
   */
  int ConnectIndex(int side);
  void  SetConnectIndex(int /*inode*/, int index) {fConnectIndex = index;}

  /**
   * it returns the shapes number of the element 
   */
  int  NShapeF();

  void CreateInterface();

  void RemoveInterfaces();

  int ExistsInterface(TPZCompElSide cpsd);
  
};

inline TPZCompEl *TPZCompElDisc::CreateC3Disc(TPZGeoElC3d *geo, TPZCompMesh &mesh, int &index) {
  //  TPZGeoEl *cop = geo;
  return new TPZCompElDisc(mesh,(TPZGeoEl *) geo,index);
}

inline TPZCompEl *TPZCompElDisc::CreateT3Disc(TPZGeoElT3d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZCompElDisc(mesh,(TPZGeoEl *) geo,index);
}

inline TPZCompEl *TPZCompElDisc::CreatePi3Disc(TPZGeoElPi3d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZCompElDisc(mesh,(TPZGeoEl *) geo,index);
}

inline TPZCompEl *TPZCompElDisc::CreatePr3Disc(TPZGeoElPr3d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZCompElDisc(mesh,(TPZGeoEl *) geo,index);
}

inline TPZCompEl *TPZCompElDisc::CreateQ2Disc(TPZGeoElQ2d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZCompElDisc(mesh,(TPZGeoEl *) geo,index);
}
//Acessar com -> TPZGeoElXXd::SetCreateFunction(createCompXXDisc);

inline TPZCompEl *TPZCompElDisc::CreateT2Disc(TPZGeoElT2d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZCompElDisc(mesh,(TPZGeoEl *)geo,index);
}

#endif
  /**
   * pointer for the geometric element associate 
   */
//  TPZGeoEl *fGeoEl;

//$Id: TPZInterfaceEl.h,v 1.19 2003-12-02 12:37:58 tiago Exp $

#ifndef ELEMINTERFACEHH
#define ELEMINTERFACEHH


#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
class TPZGeoElQ2d;

class TPZInterfaceElement : public TPZCompEl {

  /**
   * element the left of the normal a interface 
   */
  TPZCompElDisc *fLeftEl;//este podia ser um TPZCompElSide

  /**
   * element the right of the normal a interface 
   */
  TPZCompElDisc *fRightEl;//este podia ser um TPZCompElSide

  /**
   * Normal to the face element
   */
  TPZVec<REAL> fNormal;

  /**
   * Geometric element to which this element refers
   */
  TPZGeoEl *fReference;

  /**
   * Material object of this element
   */
  TPZMaterial *fMaterial;//this variable can be gotten of the element of associated volume

 public:

  static TPZCompEl *CreateInterfaceQEl(TPZGeoElQ2d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreateInterfaceTEl(TPZGeoElT2d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreateInterface1dEl(TPZGeoEl1d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreateInterfacePointEl(TPZGeoElPoint *geo, TPZCompMesh &mesh, int &index);

  //construtor do descontínuo
  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index);
  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElDisc *left,TPZCompElDisc *right,int leftside);
  TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy);
  TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy, int &index);
  TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy, TPZVec<int> &destindex,int &index);

  ~TPZInterfaceElement(){};

  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    return new TPZInterfaceElement(mesh, *this);
  }
  TPZCompEl * CloneInterface(TPZCompMesh &aggmesh,int &index) const;
  //TPZCompEl * CloneInterface(TPZCompMesh &aggmesh, TPZVec<int> &destindex,int &index) const;

  /**return the geometric element to which this element references*/
  virtual TPZGeoEl *Reference() const { return fReference;}

  TPZMaterial *Material() const { return fMaterial;}

  void SetMaterial(TPZMaterial *mat) { fMaterial = mat;}

  /**
   * it identifies the elements of left and right volume of the interface 
   */
  void VolumeEls(TPZCompEl &thirdel);

  void GetTransformsLeftAndRight(TPZTransform &tl,TPZTransform &tr);

  /**
   * it returns the right element from the element interface 
   */
  TPZCompElDisc *RightElement() {return fRightEl;}

      void SetRightElement( TPZCompElDisc* el )
      {
	 fRightEl = el;
      }

  /**
   * it returns the left element from the element interface 
   */
  TPZCompElDisc *LeftElement() {return fLeftEl;}

      void SetLeftElement( TPZCompElDisc* el )
      {
	 fLeftEl = el;
      }

  /**
   * it returns the normal one to the face from the element
   */
  void Normal(TPZVec<REAL> &normal);

/*   void SetNormal(TPZVec<REAL> &normal); */

  /**
   * it returns the number from connectivities of the element 
   */
  int NConnects();

  /**
   * Its return the connects of the left and right element associates
   */
  int ConnectIndex(int i);

  /**
   * This function should not be called
   */
  void SetConnectIndex(int node, int index);

  /**
   * it returns the dimension from the element interface
   */
  int Dimension() const {return TPZCompElDisc::gInterfaceDimension;}

  /**
   * Type of the element 
   */
  MElementType Type() { return EInterface; }

  /**declare the element as interpolated or not.
   * You may not redefine this method, because a lot of "unsafe" casts depend
   * on the result of this method\n
   * Wherever possible, use dynamic_cast instead of this method
   * @return 0 if the element is not interpolated
   */
  virtual int IsInterpolated() {return 0;}

  /**
   * it returns the shapes number of the element 
   * the associated space of interpolation is gotten 
   * of the elements left and right
   */
  int  NShapeF() {return 0;}

  /**
   * CalcStiff computes the element stiffness matrix and right hand side
   * @param ek element matrix
   * @param ef element right hand side
   */
  virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * Print attributes of the object
   */
  void Print(ostream &out = cout);

  /**
   * it verifies the existence of interfaces associates
   * with the side of an element
   * case to interface should exist and exists only a returns 1
   * case to interface should not exist and does not exist returns 1
   * otherwise returns 0
   */
  static int ExistInterfaces(TPZCompElSide &comp);

  //it returns the normal one to the face from element
  void NormalToFace(TPZVec<REAL> &normal,int leftside);

  static int FreeInterface(TPZCompMesh &cmesh);

  /**
   * reproduz na malha aglomerada aggmesh uma copia da interface da malha fina
   */
  void CloneInterface(TPZCompMesh *aggmesh);

  static int main(TPZCompMesh &cmesh);

  void EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
		     TPZVec<REAL> &errors, TPZBlock * /*flux */);

};

inline TPZCompEl *TPZInterfaceElement::CreateInterfaceQEl(TPZGeoElQ2d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZInterfaceElement(mesh,(TPZGeoEl *) geo,index);
}
inline TPZCompEl *TPZInterfaceElement::CreateInterfaceTEl(TPZGeoElT2d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZInterfaceElement(mesh,(TPZGeoEl *) geo,index);
}
inline TPZCompEl *TPZInterfaceElement::CreateInterface1dEl(TPZGeoEl1d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZInterfaceElement(mesh,(TPZGeoEl *) geo,index);
}
inline TPZCompEl *TPZInterfaceElement::CreateInterfacePointEl(TPZGeoElPoint *geo, TPZCompMesh &mesh, int &index) {
  return new TPZInterfaceElement(mesh,(TPZGeoEl *) geo,index);
}
//Acessar com -> TPZGeoElXXd::SetCreateFunction(createInterfaceXXEl);
#endif


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
  TPZCompElDisc *fLeftEl;

  /**
   * element the right of the normal a interface 
   */
  TPZCompElDisc *fRightEl;

  /**
   * Geometric element to which this element refers
   */
  TPZGeoEl *fReference;

  /**
   * Material object of this element
   */
  TPZMaterial *fMaterial;

 public:

  static TPZCompEl *CreateInterfaceQEl(TPZGeoElQ2d *geo, TPZCompMesh &mesh, int &index);
  static TPZCompEl *CreateInterfaceTEl(TPZGeoElT2d *geo, TPZCompMesh &mesh, int &index);

  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index);
  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompEl &thirdel);
  //  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoElT2d *geo,int &index);
  ~TPZInterfaceElement(){};

  /**return the geometric element to which this element references*/
  virtual TPZGeoEl *Reference() { return fReference;}

  TPZMaterial *Material() { return fMaterial;}

  void SetMaterial(TPZMaterial *mat) { fMaterial = mat;}

  /**
   * it identifies the elements of left and right volume of the interface 
   */
  void VolumeEls(TPZCompEl &thirdel);

  /**
   * it returns the right element from the element interface 
   */
  TPZCompElDisc *RightElement(){return fRightEl;}

  /**
   * it returns the left element from the element interface 
   */
  TPZCompElDisc *LeftElement(){return fLeftEl;}

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
  int Dimension(){return 2;}

  /**
   * Type of the element 
   */
  MElementType Type() { return EInterfaceDisc; }

   /**
   * it returns the associated conectivity to the element 
   */
  //int ConnectIndex(int i){return fConnectIndex;}

  /**
   * it returns the shapes number of the element 
   */
  int  NShapeF() {return 0;}

  /**
   * Computes the element stiffness matrix and right hand side
   * param ek element matrix
   * param ef element right hand side
   */
  void CalcStiff(TPZFMatrix &ek,TPZFMatrix &ef);

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

  static int FreeInterface(TPZCompMesh &cmesh);

  static int main(TPZCompMesh &cmesh);

};

inline TPZCompEl *TPZInterfaceElement::CreateInterfaceQEl(TPZGeoElQ2d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZInterfaceElement(mesh,(TPZGeoEl *) geo,index);
}
//Acessar com -> TPZGeoElXXd::SetCreateFunction(createInterfaceEl);
inline TPZCompEl *TPZInterfaceElement::CreateInterfaceTEl(TPZGeoElT2d *geo, TPZCompMesh &mesh, int &index) {
  return new TPZInterfaceElement(mesh,(TPZGeoEl *) geo,index);
}
#endif


// -*- c++ -*-

//$Id: TPZInterfaceEl.h,v 1.50 2007-09-04 12:33:16 tiago Exp $

#ifndef ELEMINTERFACEHH
#define ELEMINTERFACEHH

#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "TPZCompElDisc.h"

/// This class computes the contribution over an interface between two discontinuous elements
/**
@ingroup CompElement
*/
class TPZInterfaceElement : public TPZCompEl {

   private :

  /**
   * Element the left of the normal a interface
   */
  TPZCompElSide fLeftElSide;

  /**
   * element the right of the normal a interface
   */
  TPZCompElSide fRightElSide;

  /**
   * Normal to the face element
   */
  TPZManVector<REAL,3> fNormal;

  /** Informs the connect that this element is no longer connected to it.
   */
  void DecreaseElConnected();

  /** Informs the connect that this element is connected to it.
   */
  void IncrementElConnected();

 /**
  * Extract connects from element el.
  */
  void GetConnects(TPZCompElSide &elside, TPZVec<TPZConnect*> &connects, TPZVec<int> &connectindex);

 protected:

  /** Initialize a material data and its attributes based on element dimension, number
   * of state variables and material definitions
   */
  void InitMaterialData(TPZMaterialData &data, TPZInterpolationSpace *left, TPZInterpolationSpace *right);

  /** Compute and fill data with requested attributes */
  void ComputeRequiredData(TPZMaterialData &data,
                           TPZInterpolationSpace *left, TPZInterpolationSpace *right,
                           TPZVec<REAL> &qsi,
                           TPZVec<REAL> &LeftIntPoint, TPZVec<REAL> &RightIntPoint);

  public:

  /** Compute solution at neighbour element in a given master coordinate qsi. It returns the axes
   * at which respect derivatives are computed.
   * @param [in] Neighbor
   * @param [in] qsi
   * @param [out] sol
   * @param [out] dsol
   * @param [out] NeighborAxes
   */
  void NeighbourSolution(TPZCompElSide & Neighbor, TPZVec<REAL> & qsi, TPZVec<REAL> &sol, TPZFMatrix &dsol, TPZFMatrix &NeighborAxes);
  
  protected:

  /** Check consistency of mapped qsi performed by method TPZInterfaceElement::MapQsi by
   * comparing the X coordinate of qsi and the correspondent NeighIntPoint.
   * It return true if everything is ok or false otherwise.
   */
  bool CheckConsistencyOfMappedQsi(TPZCompElSide &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL>&NeighIntPoint);

  void ComputeSideTransform(TPZCompElSide &Neighbor, TPZTransform &transf);

  /** Computes normal.
   */
  void ComputeNormal();

 public:

  /**
   * Maps qsi coordinate at this master element to qsi coordinate at neighbor master element.
   * @param Neighbor [in] may be this->LeftElementSide() or this->RightElementSide()
   * @param qsi [in] is the point at this element master
   * @param NeighIntPoint [out] is the point at neighbor element master. X[qsi] is equal to X[NeighIntPoint]
   */
  void MapQsi(TPZCompElSide &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL> &NeighIntPoint);

  enum CalcStiffOptions{ENone = -1, EStandard /*Deprecated*/ = 0, EPenalty, EContDisc,EReferred};

  /** Constuctor to continuous and/or discontinuous neighbours.
   */
  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);

  /**
   * Copy constructor.
   */
  TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy);


  /**
   * Clone constructor to a patch mesh
   * @param mesh reference to a clone mesh
   * @param copy element to be copied
   * @param gl2lcIdx map between global(original) and local (patch) connect indexes
   */
  TPZInterfaceElement(TPZCompMesh &mesh,
                      const TPZInterfaceElement &copy,
                      std::map<int,int> &gl2lcConIdx,
                      std::map<int,int> &gl2lcElIdx);


  /**
   *
   */
  TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy, int &index);

  /**
   * Empty constructor.
   */
  TPZInterfaceElement();

  /** Default TPZCompEl constructor. SetLeftRightElements must be called
   * before any computation.
   */
  TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index);

  /** Class destructor */
  ~TPZInterfaceElement();

  /** Set neighbors.
  */
  void SetLeftRightElements(TPZCompElSide & left, TPZCompElSide & right);

  /** Makes a clone of this */
  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    return new TPZInterfaceElement(mesh, *this);
  }

  /**
   * @see class TPZCompEl
   */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> &gl2lcConMap, std::map<int,int> &gl2lcElMap) const
  {
    return new TPZInterfaceElement(mesh, *this, gl2lcConMap,gl2lcElMap);
  }

  /** Method used in TPZAgglomerateElement::CreateAgglomerateMesh
   */
  TPZCompEl * CloneInterface(TPZCompMesh &aggmesh,int &index, /*TPZCompElDisc **/TPZCompElSide & left, /*TPZCompElDisc **/ TPZCompElSide &right) const;

  /**
   * it identifies the elements of left and right volume of the interface
   */
  void VolumeEls(TPZCompEl &thirdel);

  /**
   * it returns the right element from the element interface
   */
  TPZCompEl *RightElement() {return fRightElSide.Element();}

  /**
   * it returns the left element from the element interface
   */
  TPZCompEl *LeftElement() {return fLeftElSide.Element();}

  /**
   * Returns left neighbor
   */
  TPZCompElSide &LeftElementSide(){ return this->fLeftElSide; }

  /**
   * Returns right neighbor
   */
  TPZCompElSide &RightElementSide(){ return this->fRightElSide; }

  /**
   * it returns the normal of this interface which goes from left to right neighbors
   */
  void Normal(TPZVec<REAL> &normal);

  /**
   * it returns the number from connectivities of the element
   */
  virtual int NConnects() const;

  /**
   * it returns the number from connectivities of the element related to right neighbour
   */
  int NRightConnects() const;

  /**
   * it returns the number from connectivities of the element related to left neighbour
   */
  int NLeftConnects() const;

  /**
   * Its return the connects of the left and right element associates
   */
  int ConnectIndex(int i) const;

  /**
   * This function should not be called
   */
  void SetConnectIndex(int node, int index);

  /**
   * it returns the dimension from the element interface
   */
  int Dimension() const {
     return this->Reference()->Dimension();
  }

  /**
   * Type of the element
   */
  MElementType Type() { return EInterface; }

  /**
   * Loads the solution within the internal data structure of the element
   * Is used to initialize the solution of connect objects with dependency
   * Is also used to load the solution within SuperElements*/
  virtual void LoadSolution(){
    //NOTHING TO BE DONE HERE
  }

  /**
   * CalcStiff computes the element stiffness matrix and right hand side
   * @param ek element matrix
   * @param ef element right hand side
   */
  virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * CalcResidual only computes the element residual
   * @param ef element residual
   */
  virtual void CalcResidual(TPZElementMatrix &ef);

 /**
  * Computes solution and its derivatives in the local coordinate qsi.
  * @param [in] qsi master element coordinate
   * @param [out] leftsol left finite element solution
   * @param [out] rightsol right finite element solution
   * @param [out] dleftsol left solution derivatives
   * @param [out] drightsol right solution derivatives
   * @param [out] leftaxes axes associated with the derivative of the left element
   * @param [out] rightaxes axes associated with the derivative of the right element
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi,
                               TPZVec<REAL> &normal,
                               TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                               TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes);

  /**
  * Computes solution and its derivatives in local coordinate qsi
  * @param qsi master element coordinate
  * @param phi matrix containing shape functions compute in qsi point
  * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
  * @param [in] axes indicating the direction of the derivatives
  * @param sol finite element solution
  * @param dsol solution derivatives
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                               const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol);

 /**
  * Computes solution and its derivatives in the local coordinate qsi.
  * @param qsi master element coordinate
  * @param sol finite element solution
  * @param dsol solution derivatives
  * @param axes axes associated with the derivative of the solution
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi,
                               TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes);

  void VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet);

  /**
   * Print attributes of the object
   */
  void Print(std::ostream &out = std::cout);

  /**
   * it verifies the existence of interfaces associates
   * with the side of an element
   * case to interface should exist and exists only a returns 1
   * case to interface should not exist and does not exist returns 1
   * otherwise returns 0
   */
  static int ExistInterfaces(TPZCompElSide &comp);

  static int FreeInterface(TPZCompMesh &cmesh);

  /**
   * reproduz na malha aglomerada aggmesh uma copia da interface da malha fina
   */
  void CloneInterface(TPZCompMesh *aggmesh);

  static int main(TPZCompMesh &cmesh);

  void EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
		     TPZVec<REAL> &errors, TPZBlock * /*flux */);

  /**
   * ComputeError computes the element error estimator
  */
  virtual void ComputeError(int errorid,
                              TPZVec<REAL> &errorL,
                              TPZVec<REAL> &errorR);

  /**
   * Integrate a variable over the element.
   */
   virtual void Integrate(int variable, TPZVec<REAL> & value);
   
   void IntegrateInterface(int variable, TPZVec<REAL> & value);

  void EvaluateInterfaceJumps(TPZVec<REAL> &errors);

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

#endif


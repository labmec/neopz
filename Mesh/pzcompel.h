// -*- c++ -*-
// $Id: pzcompel.h,v 1.38 2007-06-20 12:37:03 tiago Exp $

#ifndef COMPELEMHPP
#define COMPELEMHPP

#include "pzreal.h"
//#include "pzshapelinear.h"
#include <iostream>
#include <fstream>
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzsave.h"
#include "pzmaterial.h"


class TPZBlockDiagonal;

struct TPZElementMatrix;
class TPZCompMesh;
class TPZBndCond;
class TPZInterpolatedElement;
class TPZInterfaceElement;
class TPZConnect;
class TPZMaterial;
class TPZGeoEl;
class TPZGeoNode;
class TPZMatrix;
class TPZFMatrix;
class TPZBlock;

template<class T>
class TPZVec;
template<class T, int N>
class TPZManVector;
template<class T, int N>
class TPZStack;

class TPZGraphMesh;
class TPZIntPoints;

class TPZTransform;
class TPZTransfer;
#include "pzeltype.h"



/**
 * @brief Class TPZCompEl defines the interface of a computational element

 * @ingroup CompElement
 */
class TPZCompEl : public virtual TPZSaveable {

protected:

  /**
   * Computational mesh to which the element belongs
   */
  TPZCompMesh 	*fMesh;

  /**
   * Element index into mesh element vector
   */
  int fIndex;

private:
  /**
  * Reference to geometric element
  */
  //TPZGeoEl *fReference;
  /**
  * Index of reference element
  */
  int fReferenceIndex;

public:

  /**
   * Simple Constructor
   */
  TPZCompEl();

  /**
   * Simple destructor
   */
  virtual ~TPZCompEl();

  /**
   * Method for creating a copy of the element
   */
  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const = 0;

  /**
   * Method for creating a copy of the element in a patch mesh
   * Otherwise of the previous clone function, this method don't
   * copy entire mesh. Therefore it needs to map the connect index
   * from the both meshes - original and patch
   * @param mesh Patch clone mesh
   * @param gl2lcMap map the connects indexes from global element (original) to the local copy.
     */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                  std::map<int,int> & gl2lcConMap,
                                  std::map<int,int> & gl2lcElMap) const = 0;


  /**
   * put a copy of the element in the referred mesh
   */
  TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy);

  /**
   * put a copy of the element in the patch mesh
   */
  TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<int,int> &gl2lcElMap);

  /**
   * copy of the element in the new mesh whit alocated index
   */
  TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, int &index);

  /**
   * Create a computational element within mesh
   * Inserts the element within the data structure of the mesh
   * @param mesh mesh wher will be created the element
   * @param index new elemen index
   */
  TPZCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

  /**
   * Set & Get the value of gOrder
   */
  static void SetgOrder( int order );

  static int GetgOrder();

  /**
   * Returns the volume of the geometric element associated.
   */
  virtual  REAL VolumeOfEl()
 {
   if(fReferenceIndex == 0) return 0.;
   return Reference()->Volume();
 }

  /**
   * Loads the geometric element referece
   */
  virtual void LoadElementReference();

  /**
   * This method verifies if the element has the given data characteristics
   * @param var state variable number
   * @param matname pointer to material name
   **/
  virtual REAL CompareElement(int var, char *matname);

  //  /**pointer to function which returns num orthogonal functions at the point x*/
  //  static void (*fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi);


  /**
   * @name Access
     Access to private data
   */
  //@{
  /**
   * Return the type of the element
   */
  //   * the types are listed in parameter : ENoType, EOned, ETriangle, EQuadrilateral, ESubstructure*/
  virtual MElementType Type();

  virtual int IsInterface() { return 0; }

  /**
   * Return a pointer to the corresponding geometric element if such exists
   * return 0 otherwise
   */
  TPZGeoEl *Reference() const
  {
    if ( fMesh->Reference() == NULL ) return NULL;
    return (fReferenceIndex == -1) ? 0 : fMesh->Reference()->ElementVec()[fReferenceIndex];
  }



  void SetReference(int referenceindex)
  {
    fReferenceIndex = referenceindex;
//    fReference = (referenceindex == -1) ? 0 : fMesh->Reference()->ElementVec()[fReferenceIndex];
  }

//   void SetReference(TPZGeoEl *ref)
//   {
//     fReference = ref;
//     if(ref)
//     {
//      fReferenceIndex = ref->Index();
//     } else {
//       fReferenceIndex = -1;
//     }
//   }

  /**
   * Return the number of nodes of the element
   */
  virtual int NConnects() const =0;

  /**
   * Return the number of equations of the element
   */
  virtual int NEquations();

  /**
   * Return element index of the mesh fELementVec list
   */
  int Index();

  /**
   * Set element index of the mesh fELementVec list
   */
  void SetIndex(int index);

  /**
   * Return the index of the ith connectivity of the element
   * @param i connectivity index who want knows
   */
  virtual int ConnectIndex(int i) const = 0;

  /**
   * Return a pointer to the ith node
   * @param i node index
   */
  virtual TPZConnect &Connect(int i) const;

  /**
   * Dimension of the element
   */
  virtual int Dimension() const = 0;

  /**
   * Identify the material object associated with the element
   */
  virtual TPZAutoPointer<TPZMaterial> Material() const
  {
    TPZAutoPointer<TPZMaterial> result;
    if(fMesh && Reference()) result = fMesh->FindMaterial(Reference()->MaterialId());
    return result;
  }

  /**
   * Set the material associated with the object
   * @param mat new element material
   */
//  virtual void SetMaterial(TPZAutoPointer<TPZMaterial> mat) = 0;

  /**
   * Returns the reference geometric element patch
   */
  TPZGeoEl * GetRefElPatch();

	//void SetIntegrationRule(int order);
  //@}

  /**
   * @name MODIFICATION_OF_PRIVATE_DATA
   * Methods that modify private data
   */
  //@{
  /**
   * Creates corresponding graphical element(s) if the dimension matches
   * graphical elements are used to generate output files
   * @param graphmesh graphical mesh where the element will be created
   * @param dimension target dimension of the graphical element
   */
  virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension);

  /**
   * Loads the solution within the internal data structure of the element
   * Is used to initialize the solution of connect objects with dependency
   * Is also used to load the solution within SuperElements*/
  virtual void LoadSolution();

  /**
   * Set the grid of the element
   * @param mesh new reference mesh
   */
  void SetMesh(TPZCompMesh *mesh);

  /**
   * Return a pointer to the grid of the element
   */
  TPZCompMesh *Mesh() const;
  //@}

  /**
   * @name Print
   * Methods for print data structure
   */
  //@{
  /**
   * Prints element data
   * @param out indicates the device where the data will be printed
   */
  virtual void Print(std::ostream & out = std::cout);

  /**
   * Output device operator
   * @param out indicates the device where the data will be printed
   * @param el element to print
   */
  friend std::ostream& operator<<(std::ostream &out,TPZCompEl &el);

  /**
   * Prints the solution - sol - for the variable "VarName"
   * at point specified in terms of the master element coordinates
   * @param point master element coordinate to print
   * @param VarName name of variable to print
   * @param out indicates the device where the data will be printed
   */
  virtual void PrintSolution(TPZVec<REAL> &point,char *VarName,std::ostream &out);

  /**
   * Prints one coordinate index corresponding to the point to the output stream
   * @param point master element coordinate to print
   * @param CoordinateIndex index of the coordinate corresponding to the point
   * @param out indicates the device where the data will be printed
   */
  virtual void PrintCoordinate(TPZVec<REAL> &point,int CoordinateIndex,std::ostream &out);

  /**
   * Prints the variables names associated with the element material
   * @param VarName pointer to variable parameter wha want to print
   * @param out indicates the device where the data will be printed
   */
  virtual void PrintTitle(char *VarName,std::ostream &out);
  //@}

  /**
   * Sets the orthogonal function which will be used throughout the program
   * by default this function is the Chebyshev function
   * @param orthogonal pointer to a function which will be used to generate the shape functions
   */
  static void SetOrthogonalFunction(void (*orthogonal)(REAL x,int num,TPZFMatrix & phi,
						       TPZFMatrix & dphi));

  //  /**
  //   * Coarsen the group of elements in elements
  //   */
  //   virtual void Coarsen(TPZVec<int> &elementindexes);

  /**
   * Compute the element stifness matrix
   * @param ek element stiffness matrix
   * @param ef element loads matrix
   */
  virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);

  /**
   * Compute the element right hand side
   * @param ef element load vector(s)
   */
  virtual void CalcResidual(TPZElementMatrix &ef);

  /**
   * Implementation of the orthogonal Chebyshev functions
   * @param x point where the Chebyshev function will be evaluated
   * @param num number of functions
   * @param phi values of the function
   * @param dphi values of derivative of the function
   */
  static void Chebyshev(REAL x,int num,TPZFMatrix &phi,TPZFMatrix &dphi);

  /**
   * Divide the computational element
   * if interpolate = 1, the solution is interpolated to the sub elements
   * This method needs to be implemented in the derived classes
   * @param index  index of the element which is being divided
   * @param subindex element vector where will be created the divided elements
   * @param interpolate boolean variable to indicates if the solution will be interpolated to the sub elements
   */
  virtual void Divide(int index, TPZVec<int> &subindex, int interpolate = 0);

  /**
   * Projects the flux function on the finite element space
   * @param ekmat element stiffness matrix
   * @param efmat element loads matrix
   */
  virtual void ProjectFlux(TPZElementMatrix &ek,TPZElementMatrix &ef);

  /**
   * Perform an error estimate on the elemen
   * @param fp function pointer which computes the exact solution
   * @param true_error (output)  the true error of the solution
   * @param L2_error (output) the L2 norm of the error of the solution
   * @param flux (input) value of the interpolated flux values
   * @param estimate (output) estimated error based on the implemented criterium
   */
  virtual void EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
			     TPZVec<REAL> &errors,TPZBlock *flux);

  /**
   * ComputeError computes the element error estimator
  */
  virtual void ComputeError(int errorid, TPZVec<REAL> &error){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented.\n";
  }

  /**
   * Integrate a variable over the element.
   */
   virtual void Integrate(int variable, TPZVec<REAL> & value){
    value.Fill(0.);
    PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented.\n";
   }

  /**
   * Calculates the solution - sol - for the variable var
   * at point qsi, where qsi is expressed in terms of the
   * master element coordinates
   * @param qsi master element coordinate
   * @param var variable name
   * @param sol vetor for the solution
   */
  virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol);

 /**
  * Computes solution and its derivatives in the local coordinate qsi.
  * @param qsi master element coordinate
  * @param sol finite element solution
  * @param dsol solution derivatives
  * @param axes axes associated with the derivative of the solution
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi,
                               TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes) = 0;

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
                               TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes) = 0;

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
                               const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol) = 0;

  /**
   * Builds the list of all connectivities related to the element including the
   * connects pointed to by dependent connects
   * Note : this method does not reset the stack to zero. The calling
   * method should do this
   * @param connectlist stack to receive the list
   */
  virtual void BuildConnectList(TPZStack<int> &connectlist);

  /**
   * Returns 1 if the element has at least one dependent node
   * returns 0 otherwise
   */
  virtual int HasDependency();

  /**
   * Domain Decomposition: Misael
   * This method will eliminate the nodes which are internal to the element from
   * the datastructure of the grid
   * After calling this method, the superelement will statically condense the
   * internal equations
   */
  virtual void ReduceInternalNodes() { }

  /**
   * Set the index i to node inode
   * @param inode node to set index
   * @param index index to be seted
   */
  virtual void SetConnectIndex(int inode, int index) = 0;

  /**
   * Calculates the diagonal block
   * @param connectlist stack list to calculates the diagonal block
   * @param block object to receive the diagonal block
   */
  virtual void CalcBlockDiagonal(TPZStack<int> &connectlist, TPZBlockDiagonal & block);

  REAL MaximumRadiusOfEl();

  REAL LesserEdgeOfEl();

  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);

private:
  /**
   * Default interpolation order
   */
    static int gOrder;


};


/**
 * Class TPZCompElSide implemments computational element side
 * @ingroup CompElement
 */

class TPZGeoElSide;

/// this class represents a computational element and a side
/**
This class was created to implement all algorithms associated with element/sides
Objects of this class are mostly temporary
@ingroup interpolation
*/
class TPZCompElSide {

  /**
   * Pointer to the computational element
   */
  TPZCompEl *fEl;

  /**
   * Index of the object side
   */
  int fSide;


public:

  //    /*PARA TESTES*/ TPZGeoEl *georeftest;

  /**
   * Simple Constructor
   */
  TPZCompElSide();

  /**
   * Creates a computational element side from an object TPZCompElSide
   * @param celside reference computational element side
   */
  TPZCompElSide(const TPZCompElSide &celside);

  /**
   * Creates an computational element side given a computational element
   * and the side index
   * @param cel pointer to the computational element
   * @param side index of the side
   */
  TPZCompElSide(TPZCompEl *cel,int side);

  /**
   * Gives a pointer to the reference computational element
   */
  TPZCompEl *Element() const {return fEl;}

  /**
   *  Set computational element pointer.
   */
  void SetElement(TPZCompEl* el){ fEl = el;}

  /**
   * Return the side index
   */
  int Side() const {return fSide;}

  /**
   * Set the side index
   * @param side new side index
   */
  void SetSide(int side);

  /**
   * Verifies if the object is non null (initialized)
   */
  int Exists() const {return fEl != 0;}

  /**
   * Reference to the geometric element
   */
  TPZGeoElSide Reference() const;

  /**
   * Returns all connected elements which have level higher to the current element
   * @param elsidevec side elements vector
   * @param onlyinterpolated if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
   * @param removeduplicates if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
   */
  void HigherLevelElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);

  /**
   * pushes all element/sides which have higher dimension than the current element/side
   * @param elsidevec side elements vector where elements/side will be put
   * @param onlyinterpolated if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
   * @param removeduplicates if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
   */
  void HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
  /**
   * Returns all connected elements to the current element
   * @param elsidevec side elements vector
   * @param onlyinterpolated if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
   * @param removeduplicates if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
   */
  void ConnectedElementList(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated, int removeduplicates);

  /**
   * Returns all connected elements which have equal level to the current element
   * This method will not put this on the stack
   * @param elsidevec side elements vector
   * @param onlyinterpolated  if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
   * @param removeduplicates  if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
   */
  void EqualLevelElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);

  /**
   * Returns all connected elements which have level lower to the current element
   * @param onlyinterpolated if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
    // if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
  */
  TPZCompElSide LowerLevelElementList(int onlyinterpolated);

  /**
   * Will remove elements which are direct neighbours from elvec (and elsides)
   * the method checks between any two elements in the list whether they are of equal level
   * and whether they are neighbours. If they are neighbours, one of the elements will be removed
   * from the list
   * the method NeighbourExists between any two elements of equal level will return 0
   * @param elvec computational element side vector
   */
  static void RemoveDuplicates(TPZStack<TPZCompElSide> &elvec);


  /**
   * Remove entries of the vector which share a connect along the side This should be
   * equivalent to RemoveDuplicates.
   * @param expandvec vector of TPZCompElSide objects
   */
  static void RemoveConnectDuplicates(TPZStack<TPZCompElSide> &expandvec);

  /**
   * Find the list element/side of the current element restrict nodes and elements
   * @param onlyinterpolated if ==1 only elements derived from TPZInterpolated will be put on the stack
   */
   // acha uma lista elemento/lado de elementos e n� restritos ao atual elemento*/
  void ExpandConnected(TPZStack<TPZCompElSide> &expandvec,int onlyinterpolated);

  /**
   * returns the element with lowest id of all direct neighbours of expandvec
   * @param expandvec elements whose neighbours will be checked
   * @param onlyinterpolated if == 1, only elements derived from TPZInterpolatedElement will be checked
   */
  TPZCompElSide LowerIdElementList(TPZCompElSide &expandvec,int onlyinterpolated);

  //inline////////////////////////////////////////////////////////////////////////

  /**
   * Return the index of the middle side connect alon fSide
   */
  int ConnectIndex() const{
    if(fEl) return fEl->ConnectIndex(fSide);
    else return -1;
  }

  bool operator != (const TPZCompElSide &other);
  bool operator == (const TPZCompElSide &other);


};
//  std::ostream & operator << (std::ostream &out,const TPZCompElSide &celside);

inline void TPZCompEl::CreateGraphicalElement(TPZGraphMesh &, int) {
  std::cout << "TPZCompEl::CreateGrafEl called\n";
}


inline void TPZCompEl::CalcStiff(TPZElementMatrix &,TPZElementMatrix &){
  std::cout << "TPZCompEl::CalcStiff(*,*) is called." << std::endl;
}

inline void TPZCompEl::ProjectFlux(TPZElementMatrix &ek,TPZElementMatrix &ef) {
  std::cout << "TPZCompEl::ProjectFlux is called." << std::endl;
}

inline bool TPZCompElSide::operator != (const TPZCompElSide &other)
{
  return (other.Element() != Element() || other.Side() != Side());
}

inline bool TPZCompElSide::operator == (const TPZCompElSide &other)
{
  return (other.Element() == Element() && other.Side() == Side());
}

inline std::ostream &operator << (std::ostream &out,const TPZCompElSide &celside)
{
  out << "Side = " << celside.Side()
      << " element: " << celside.Element()->Index()
      << std::endl;
  return out;
}

inline int TPZCompEl::Index() {
  return fIndex;
}

inline void TPZCompEl::SetgOrder( int order )
{
  gOrder = order;
}

inline int TPZCompEl::GetgOrder( )
{
  return gOrder;
}


#endif

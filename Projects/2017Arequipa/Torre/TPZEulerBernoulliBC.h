#ifndef TPZEulerBernoulliBCH
#define TPZEulerBernoulliBCH

#include "pzeltype.h"
//#include "pzinterpolationspace.h"
#include "pzcompel.h"
#include "TPZEulerBernoulliBeamData.h"
#include "pzquad.h"

/**
  * Implements boundary condition (external forces, masses and supports) fot the
  * Euler-Bernoulli beam element (TPZEulerBernoulliBeam).
  * @author Tiago
  * @since Feb 2017
 */
class TPZEulerBernoulliBC : public TPZCompEl {

public:

  enum EBCType{ ENone = 0,
                EForce = 1,
                ESupport = 2 };

  struct BCVal{
    BCVal(){
      fType = ENone;
      fVal = 0.;
    }
    BCVal(int type, REAL val){
      fType = type;
      fVal = val;
    }
    int fType; //EBCType
    REAL fVal;
  };

private:


  ////////////// material and section data //////////////////

  TPZAutoPointer< TPZEulerBernoulliBeamData > fPropertyData;

  ///Connect
  int fConnectIndex;

  ///Array of boundary condition values
  TPZManVector<BCVal,6> fBCVal;

  ///Array of nodal masses in x,y,z directions
  TPZManVector<REAL,3> fMasses;

  //////////////////////////////////////////////////////////

  //Search for a TPZCompEl neighbour
  TPZCompElSide FindNeighbourCompEl() const;

  //Gives element stiffness matrix and load vector
  void StiffnessMatrix(TPZFMatrix<STATE> &K, TPZFMatrix<STATE> &F);

  //Gives element mass matrix
  void MassMatrix(TPZFMatrix<STATE> &M) const;

public:

  /**
   * Simple Constructor
   */
  TPZEulerBernoulliBC();

  /**
   * Simple destructor
   */
  virtual ~TPZEulerBernoulliBC();

  /**
   * Method for creating a copy of the element
   */
  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const{
    DebugStop();
    return NULL;
  }

    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const {
		DebugStop();
	}

	/**
   * Method for creating a copy of the element in a patch mesh
   * Otherwise of the previous clone function, this method don't
   * copy entire mesh. Therefore it needs to map the connect index
   * from the both meshes - original and patch
   * @param mesh Patch clone mesh
   * @param gl2lcMap map the connects indexes from global element (original) to the local copy.
     */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                  std::map<int64_t,int64_t> & gl2lcConMap,
                                  std::map<int64_t,int64_t> & gl2lcElMap) const
  {
    DebugStop();
    return NULL;
  }

  /**
   * Create a computational element within mesh
   * Inserts the element within the data structure of the mesh
   * @param mesh mesh wher will be created the element
   * @param index new elemen index
   */
  TPZEulerBernoulliBC(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);

  //Copy constructor
  TPZEulerBernoulliBC(const TPZEulerBernoulliBC &cp);

  /**
   * Return the type of the element
   */
  //   * the types are listed in parameter : ENoType, EOned, ETriangle, EQuadrilateral, ESubstructure*/
  virtual MElementType Type(){
    return EPoint;
  }

  /**
   * Return the number of nodes of the element
   */
  virtual int NConnects() const{
    return 1;
  }

  /**
   * Return the number of equations of the element
   */
  virtual int NEquations(){
    return 6;
  }

  /**
   * Return the index of the ith connectivity of the element
   * @param i connectivity index who want knows
   */
  virtual int64_t ConnectIndex(int i) const{
    if(i != 0) DebugStop();
    return this->fConnectIndex;
  }

  /**
   * Set the index i to node inode
   * @param inode node to set index
   * @param index index to be seted
   */
  virtual void SetConnectIndex(int inode, int64_t index){
    if(inode != 0) DebugStop();
    this->fConnectIndex = index;
  }

  /**
   * Dimension of the element
   */
  virtual int Dimension() const{
    return 0;
  }

  ///Returns equations associated to this element
  void GetEquationIndices(TPZVec<int> &indices) const;

  void GetSolutionVector(TPZFMatrix<STATE> &u);

  /**
   * No TPZMaterial is associated to this computational element
   */
//	virtual TPZMaterial* Material() const{
//    DebugStop();
//    return NULL;
//  }

  /**
   * Prints element data
   * @param out indicates the device where the data will be printed
   */
  virtual void Print(std::ostream & out = std::cout) const;

  /**
   * Compute the element stifness matrix
   * @param ek element stiffness matrix
   * @param ef element loads matrix
   */
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);

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
                               TPZVec<REAL> &sol, TPZFMatrix<STATE> &dsol,TPZFMatrix<STATE> &axes, const int whichSol)
  {
    DebugStop();
  }

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
                               TPZVec<REAL> &leftsol, TPZFMatrix<STATE> &dleftsol,TPZFMatrix<STATE> &leftaxes,
                               TPZVec<REAL> &rightsol, TPZFMatrix<STATE> &drightsol,TPZFMatrix<STATE> &rightaxes, const int whichSol)
  {
    DebugStop();
  }

  /**
  * Computes solution and its derivatives in local coordinate qsi
  * @param qsi master element coordinate
  * @param phi matrix containing shape functions compute in qsi point
  * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
  * @param [in] axes indicating the direction of the derivatives
  * @param sol finite element solution
  * @param dsol solution derivatives
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphix,
                               const TPZFMatrix<STATE> &axes, TPZVec<REAL> &sol, TPZFMatrix<STATE> &dsol, const int whichSol)
  {
    DebugStop();
  }

  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid) const;

  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);

  ///Defines nodal masses
  void SetMasses( const TPZVec< REAL > & m ){
    if(m.NElements() != 0 && m.NElements() != 3) DebugStop();
    this->fMasses = m;
  }

  ///Defines boundary condition data
  void SetBCVal(const TPZVec< BCVal > &val){
    if(val.NElements() != 0 && val.NElements() != 6) DebugStop();
    this->fBCVal = val;
  }

  ///Sets fPropertyData attribute and defines boundary condition data
  void SetData(TPZAutoPointer< TPZEulerBernoulliBeamData > PropertyData){
    if(!PropertyData) DebugStop();
    fPropertyData = PropertyData;
  }

  //Const access do bc values
  const BCVal & GetBCVal( int i ) const { return fBCVal[i]; }

  //Returns true if at least one equation is a support boundary condition
  bool HasSupport() const{
    for(int i = 0; i < 6; i++){
      if(fBCVal[i].fType == ESupport) return true;
    }
    return false;
  }

    
    /** 2017 To make derived from TPZInterpolationSpace */
    virtual void SetPreferredOrder ( int order ) {
        order = 1;
    }
    virtual void PRefine ( int order ) {
        order += 1;
    }
    virtual const TPZIntPoints &GetIntegrationRule() const {
        static TPZInt1Point gin;
        return gin;
    }
    
    virtual TPZIntPoints &GetIntegrationRule() {
        static TPZInt1Point gin;
        return gin;
    }
    virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi) {
        int r=0;
        r+=1;
    }
    virtual int NShapeF() const {
        return 1;
    }
    virtual int NConnectShapeF(int icon, int order) const {
        return 1;
    }
    virtual int NSideConnects(int iside) const {
        return 1;
    }
    virtual int SideConnectLocId(int icon,int is) const {
        return 1;
    }

};


#endif


#ifndef TPZEulerBernoulliBeamH
#define TPZEulerBernoulliBeamH

#include "pzcompel.h"
#include "pzquad.h"
//#include "pzinterpolationspace.h"
#include "TPZEulerBernoulliBeamData.h"

/**
  * Implements a Euler-Bernoulli beam element.
  * To compute selfweight loads, Z is adopted as the vertical upward direction, i.e. gravity = -Z
  * @author Tiago
  * @since Feb 2017
 */
class TPZEulerBernoulliBeam : public TPZCompEl {
  //  class TPZEulerBernoulliBeam : public TPZInterpolationSpace {

private:

  bool SimplifiedCorotational() const;

  void GetInitialPosVec(TPZFMatrix<STATE> &InitialPos);

  void Angles(const TPZVec<REAL> &v, TPZVec<REAL> &angles) const;
  void ComputeRigidBodyMotionRotations(const TPZFMatrix<STATE> & uGlobal, TPZVec<REAL> & angles) const;
  REAL ComputeTwistRigidBodyAngle(const TPZFMatrix<STATE> & RotOrig,
	  const TPZFMatrix<STATE> & RotDeformed) const;
  void ComputeRotationBetweenTwoVectors(const TPZVec<REAL> &a, const TPZVec<REAL> &b,
	  TPZFMatrix<STATE> & RotBetweenVectors) const;

  void TransformDisplFromGlobalToLocalRefSystem(const TPZFMatrix<STATE> & uGlobal,
	  const TPZFMatrix<STATE> & Rot,
	  TPZFMatrix<STATE> &uLocal,
                                                REAL &twistAngle,
												TPZFMatrix<STATE> &twistRot);

  ////////////// material and section data //////////////////

  //Autopointer to material and section library
  TPZAutoPointer< TPZEulerBernoulliBeamData > fPropertyData;

  //Materail and section ids
  int fMaterialId, fSectionId;

  //The section is twisted over its axis of fAlfa (in rad)
  REAL fAlfa;

  ///Connects
  TPZManVector<int,2> fConnectIndexes;

  /** The beam may be longer (positive value) or shorter then designed. That fabrication error
  * results in normal forces in load vector.
  */
  REAL fFabricationErrorStrain;

  //////////////////////////////////////////////////////////

  //Search for a TPZCompEl neighbour
  TPZCompElSide FindNeighbourCompEl(int myside) const;

  //Returns element length
  void L(const TPZFMatrix<STATE> &u, REAL & Lorig, REAL & Ldef) const;

public:

  REAL LOriginal() const;

private:

  //Gives element stiffness matrix and load vector in global reference system
	void StiffnessMatrix(TPZFMatrix<STATE> &K, TPZFMatrix<STATE> &F);

  //Gives elements stiffness matrix and load vector in local reference system. It also returns roation matrix
	void LocalStiffnessMatrix(const TPZFMatrix<STATE> & uGlobal, TPZFMatrix<STATE> &K, TPZFMatrix<STATE> &F, TPZFMatrix<STATE> &Rot);

  //Gives element mass matrix in global reference system
	void MassMatrix(TPZFMatrix<STATE> &M);

  //Gives element mass matrix in local reference system
	void LocalMassMatrix(TPZFMatrix<STATE> &M) const;

  //Gets rotation matrix from local to global coordinate systems
	void RotateMatrix(TPZFMatrix<STATE> &Rotate, const TPZFMatrix<STATE> & u, bool uComposeRotation);

public:

  /**
   * Simple Constructor
   */
  TPZEulerBernoulliBeam();

  /**
   * Simple destructor
   */
  virtual ~TPZEulerBernoulliBeam();

  /**
   * Method for creating a copy of the element
   */
  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const{
    DebugStop();
    return NULL;
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
  TPZEulerBernoulliBeam(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);

  //Copy constructor
  TPZEulerBernoulliBeam(const TPZEulerBernoulliBeam &cp);

  //Returns section id
  int SectionId() const { return fSectionId; }

  /**
   * Return the type of the element
   */
  //   * the types are listed in parameter : ENoType, EOned, ETriangle, EQuadrilateral, ESubstructure*/
  virtual MElementType Type(){
    return EOned;
  }

  /**
   * Return the number of nodes of the element
   */
  virtual int NConnects() const{
    return 2;
  }

  /**
   * Return the number of equations of the element
   */
  virtual int NEquations(){
    return 12;
  }

  /**
   * Return the index of the ith connectivity of the element
   * @param i connectivity index who want knows
   */
  virtual int64_t ConnectIndex(int i) const{
    return this->fConnectIndexes[i];
  }

  /**
   * Set the index i to node inode
   * @param inode node to set index
   * @param index index to be seted
   */
  virtual void SetConnectIndex(int inode, int index){
    this->fConnectIndexes[inode] = index;
  }

  /**
   * Dimension of the element
   */
  virtual int Dimension() const{
    return 1;
  }

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
   * Post process forces and moments at both ends (nodes) of this beam element
   * in local reference system i.e. normal and shear forces and moments
   */
	void GetStaticNodalForces(TPZFMatrix<STATE> & NodalF);

  REAL GetStaticNormalForce();

  /**
   * Gather element solution vector.
   */
  void GetSolutionVector(TPZFMatrix<STATE> &u);

  ///Returns equations associated to this element
  void GetEquationIndices(TPZVec<int64_t> &indices) ;

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
	  TPZVec<REAL> &sol, TPZFMatrix<STATE> &dsol, TPZFMatrix<STATE> &axes, const int whichSol)
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
							   TPZVec<REAL> &leftsol, TPZFMatrix<STATE> &dleftsol, TPZFMatrix<STATE> &leftaxes,
							   TPZVec<REAL> &rightsol, TPZFMatrix<STATE> &drightsol, TPZFMatrix<STATE> &rightaxes, const int whichSol)
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

  //Sets element data
  void SetPropertyData( TPZAutoPointer< TPZEulerBernoulliBeamData > PropertyData,
                        int MaterialId, int SectionId, REAL alfa,
                        REAL fabErrorStrain  );

  /**
   * Divide the computational element
   * if interpolate = 1, the solution is interpolated to the sub elements
   * This method needs to be implemented in the derived classes
   * @param index  index of the element which is being divided
   * @param subindex element vector where will be created the divided elements
   * @param interpolate boolean variable to indicates if the solution will be interpolated to the sub elements
   */
  virtual void Divide(int64_t index, TPZVec<int64_t> &subindex, int interpolate = 0);

  /** @brief adds the connect indexes associated with base shape functions to the set */
  virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const {
  }

  /**
  * @brief Set the index i to node inode
  * @param inode node to set index
  * @param index index to be seted
  */
  virtual void SetConnectIndex(int inode, int64_t index) { }

/** 2017 To make derived from TPZInterpolationSpace */
        virtual void SetPreferredOrder ( int order ) {
            order = 1;
        }
        virtual void PRefine ( int order ) {
            order += 1;
        }
        virtual const TPZIntPoints &GetIntegrationRule() const {
            TPZInt1d gin;
            return gin;
        }
        
        virtual TPZIntPoints &GetIntegrationRule() {
            TPZInt1d gin;
            if (this->fIntegrationRule) {
                return *fIntegrationRule;
            }
            else
            {
                return gin;
            }
        }
        virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi) {
            int r=0;
            r+=1;
        }
        virtual int NShapeF() const {
            return 2;
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


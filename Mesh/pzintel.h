/**
 * @file
 * @brief Contains declaration of TPZInterpolatedElement class which implements computational element of the interpolation space.
 */

#ifndef PZINTEL_H
#define PZINTEL_H

#include "pzcompel.h"
#include "pzinterpolationspace.h"
//#include "pzmaterial.h"
struct TPZElementMatrix;

class TPZIntPoints;
#include "TPZCompElDisc.h"
#include "TPZOneShapeRestraint.h"
#include "pzblockdiag.h"

/**
 * @brief Implements computational element based on an interpolation space. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/** All h and p adaptive methods are implemented in this class\n
 * derived classes need to implement the behaviour (interface) defined in this class\n
 * this makes the adaptive process extremely general\n
 */
class TPZInterpolatedElement : public TPZInterpolationSpace {
	
protected:
	
	/**
	 * @brief Updates the interpolation order of all neighbouring elements along side
	 * to have side order equal to the side order of the current element
	 * @param side of the the element
	 * @param elvec vector of elements whose side order will be set
	 */
	void UpdateNeighbourSideOrder(int side, TPZVec<TPZCompElSide> &elvec);
	
	/**
	 * @brief Computes the minimum interpolation order of the elements contained in elementset
	 * this method is used to identify the side order of a set of equal level elements
	 * @param elementset vector of element/sides which will be used to compute the maximum order
	 * @return resulting side order
	 */
	static int ComputeSideOrder(TPZVec<TPZCompElSide> &elementset);
	
	
public:
	/**
	 * @brief Constructor with a mesh and geometric element as arguments
	 * @param mesh mesh object into which the element will insert itself
	 * @param reference reference object to which this element will refer
	 * @param index index in the vector of elements of mesh where this element was inserted
	 */
	TPZInterpolatedElement(TPZCompMesh &mesh, TPZGeoEl *reference, int64_t &index);
	
	/**
	 * @brief Constructor aimed at creating a copy of an interpolated element within a new mesh
	 */
	TPZInterpolatedElement(TPZCompMesh &mesh, const TPZInterpolatedElement &copy);
	
	/**
	 * @brief Copy the given element into a new patch mesh
	 * @param mesh patch mesh
	 * @param copy element to be copied
	 * @param gl2lcElMap map the indexes of the orginal mesh to the pacht mesh
	 */
	TPZInterpolatedElement ( TPZCompMesh &mesh,
							const TPZInterpolatedElement &copy,
							std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZInterpolatedElement();
	/** @brief Destructor, does nothing */
	virtual ~TPZInterpolatedElement();
	
	/** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh);
        
        int ClassId() const;
	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
	/**
	 * @name data access methods
	 * @brief Methods which allow to access the internal data structure of the element
	 * @{
	 */
	
	/** @brief Prints the relevant data of the element to the output stream */
	virtual void Print(std::ostream &out = std::cout) const;
	
	/** @brief Returns the total number of shapefunctions*/
	int NShapeF() const;
	
	/** @brief Returns the number of shape functions on a side*/
	int NSideShapeF(int side) const;
	
	/** @brief Returns the number of dof nodes along side iside*/
	virtual int NSideConnects(int iside) const = 0;
	
	/**
	 * @brief Returns the local node number of icon along is
	 * @param icon connect number along side is
	 * @param is side which is being queried
	 */
	virtual int SideConnectLocId(int icon,int is) const = 0;
	
	/**
	 * @brief Returns the local id of the connect in the middle of the side
	 * @param is side which is being queried
	 */
	virtual int MidSideConnectLocId(int is) const;
	
    /**
     * @brief Returns a reference to the connect in the middle of the side
     * @param is side which is being queried
     */
    virtual TPZConnect &MidSideConnect(int is) const;
    
	
	/** @brief Returns the index of the c th connect object along side is*/
	int64_t SideConnectIndex(int icon,int is) const;
	
	/** @brief Returns a pointer to the icon th connect object along side is */
	TPZConnect &SideConnect(int icon,int is);
	
	/** @brief Returns the dimension of the element */
	virtual int Dimension() const = 0;
	
	/** @brief Returns the number of corner connects of the element*/
	virtual int NCornerConnects() const = 0;
	
	/** @brief Returns the number of connect objects of the element*/
	virtual int NConnects() const = 0;
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const;

	
	/** @brief Identifies the interpolation order of all connects of the element different from the corner connects */
	/** Note there is a diference between the actual side order returned by this method
	 * and the side order which the element naturally would have, which is returned by
	 * PreferredSideOrder
	 * @param ord vector which will be filled with the side orders of the element */
	virtual void GetInterpolationOrder(TPZVec<int> &ord) = 0;
	
	/** @brief Returns the preferred order of the polynomial along side iside */
	virtual int PreferredSideOrder(int iside) = 0;
	
	/** 
	 * @brief Adjusts the preferredSideOrder for faces
	 * @param side : side for which the order needs adjustment
	 * @param order : original order which has to be compared with the sides
	 */
	int AdjustPreferredSideOrder(int side, int order);
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int EffectiveSideOrder(int side) const = 0;
	
    /// return true if the connects associated with the side have dependency with large and if the dependency dimensions match
    bool VerifyConstraintConsistency(int side, TPZCompElSide large) const;
	/** @} */
	
	/**
	 * @name Data modification methods
	 * @brief These methods which will modify the local datastructure of the element
	 * @{
	 */
	
	/** @brief Sets the node pointer of node i to nod */
	virtual void SetConnectIndex(int i, int64_t connectindex)=0;
	
	virtual void SetIntegrationRule(int order) {
		std::cout << "TPZInterpolatedElement::SetIntegrationRule called\n";
	}
	
	/** @brief Sets the preferred interpolation order along a side */
	/** 
	 * This method only updates the datastructure of the element
	 * In order to change the interpolation order of an element, use the method PRefine
	 */
	virtual void SetPreferredOrder(int order) = 0;
	
public:
	/** @brief Sets the interpolation order of side to order */
	/** 
	 * This method only updates the datastructure of the element and
	 * updates the blocksize of the associated connect object
	 * @note DO NOT CALL THIS METHOD
	 */
	virtual void SetSideOrder(int side, int order) = 0;
	
	/** @brief Impose an interpolation order on a given side (without using computesideorder) */
	virtual void ForceSideOrder(int side, int order);
	
	/** @} */
	
public:
	
	
	/**
	 * @name Computational methods
	 * @brief Methods used to perform computations on the interpolated element
	 * @{
	 */
	
	/** @brief Compute the values of the shape function along the side*/
	virtual void SideShapeFunction(int side, TPZVec<REAL> &point, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) = 0;
	
	/** @} */
	
	/**
	 * @name Post processing methods
	 * @brief The TPZInterpolatedElement class provides the user with a variety of methods for post-processing*/
	/**
	 * Methods for error evaluation\n
	 * Methods for computing derived post processed values (depending on the variational statement)\n
	 * @{
	 */
	
public:
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes indicating the direction of the derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
    
    public:
    
    /** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] qsi point in master element coordinates 
	 * @param[in] data stores all input data
	 */
    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data);
	
	/**
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions with respect of global coordinates: D[phi,x], D[phi,y], D[phi,z]
	 * @param axes axes indicating the direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
								 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.\n
	 * This method will function for both volumetric and interface elements
	 * @param qsi master element coordinate of the interface element
	 * @param normal unitary normal vector
	 * @param leftsol finite element solution
	 * @param dleftsol solution derivatives
	 * @param leftaxes axes associated with the left solution
	 * @param rightsol finite element solution
	 * @param drightsol solution derivatives
	 * @param rightaxes axes associated with the right solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZVec<REAL> &normal,
								 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
								 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes);
	
	/**
	 * @brief Compare the L2 norm of the difference between the švarš solution of the current element with
	 * the švarš solution of the element which is pointed to by the geometric element.
	 * @param var variable index indicating which difference is being integrated
	 * @param matname reference material name
	 */
	/**
	 * In order to use this method, call ResetReference on the geometric mesh to which the geometric reference element belongs
	 * and call LoadReference on the mesh of the element with which to compare the solution\n
	 * This method only computes the error is the name of the material is matname
	 */
	virtual REAL CompareElement(int var, char *matname);
	
	/** @} */
	
	/**
	 * @name Adaptivity methods
	 * @brief Methods which help to implement the adaptive process
	 * @{
	 */
	
	/**
	 * @brief Implement the refinement of an interpolated element
	 * @param index [in] index of the current element in the mesh
	 * @param sub [out] indices of the subelements (they are already inserted in the mesh)
	 * @param interpolatesolution [in] if == 1 the solution of the original element is projected on the
	 * solution of the subelements
	 * @see PRefine
	 * @note This is the user interface to adaptive refinement of this class
	 */
	/**
	 * Divides the current element into subelements. Inserts the subelements in the mesh of the element
	 * and returns their indices
	 */
	void Divide(int64_t index,TPZVec<int64_t> &sub,int interpolatesolution = 0);
	
	/**
	 * @brief Changes the interpolation order of a side. Updates all constraints and block sizes\n
	 * @param order interpolation order which the user requests
	 * @see ComputeSideOrder
	 * @note This is the user interface to adaptive refinement of this class
	 */
	/**
	 * This call will not šnecessarilyš modify the interpolation order of the side. The interpolation
	 * order of neighbouring elements need to remain compatible. The actual order is obtained by calling ComputeSideOrder
	 */
	void PRefine(int order);
	
	/** @brief Compute the shapefunction restraints which need to be applied to the shape functions
     on the side of the element*/
	virtual void RestrainSide(int side, TPZInterpolatedElement *neighbour, int neighbourside);
	
	/**
     * @enum MInsertMode
	 * @brief Defines a flag indicating the state of creation/deletion of the element
	 * This has an impact on how constraints are being computed
	 * @param EInsert The element is being inserted
	 * @param EDelete The element is being deleted
	 */
	enum MInsertMode {EInsert,EDelete};
	
	/** @brief Delete the restraints on the nodes of the connected elements if necessary
	 * @param mode indicates insertion or deletion of element
	 * @note THIS IS A VERY TRICKY METHOD
	 */
	/**
	 * Depending on the insert mode, the neighbouring elements will need to recompute their constraints \n
	 * insert indicates whether the element will be refined, coarsened, inserted or deleted \n
	 * This method is only called during deletion
	 */
	virtual void RemoveSideRestraintsII(MInsertMode mode);
	
	/**
	 * @brief Removes the side restraints of the current element along side with respect to neighbour/side
	 * @param side side of the current element which contains the constrained connect
	 * @param neighbour element/side with respect to which the connect is restrained
	 */
	/** 
	 * This method checks (extensively) if the relative positions between both elements makes sense
	 */
	virtual void RemoveSideRestraintWithRespectTo(int side, const TPZCompElSide &neighbour);
	
	/**
	 * @brief Will recompute the restraints of all connects which are restrained by this side
	 * @param side side of the large element
	 */
	/** This method will be called for a side when a connected lower dimension side is changing order */
	void RecomputeRestraints(int side);
	
    /// Add a shape restraint (meant to fit the pyramid to restraint
    virtual void AddShapeRestraint(TPZOneShapeRestraint restraint)
    {
        DebugStop();
    }
    
	/**
	 * @brief Accumulates the transfer coefficients between the current element and the
	 * coarse element into the transfer matrix, using the transformation t
	 * @param coarsel larger element with respect to which the transfer matrix is computed
	 * @param t transformation which maps the master element space of the current element into the master element space of the coarse element
	 * @param transfer transfer matrix mapping the solution of the coarse mesh into the fine mesh
	 */
	/** This method forms the basis for the multigrid method */
	virtual void BuildTransferMatrix(TPZInterpolatedElement &coarsel, TPZTransform<> &t, TPZTransfer<STATE> &transfer);
	
	/**
	 * @brief Verify the neighbours of the element and create a node along this side
	 * @note This is the central method for h-p adaptivity : the constructor of the element
	 * simply calls CreateMidSideConnect which does all the necessary adjustments
	 */
	/** If necessary. \n This method returns the index of the newly created node */
	virtual int64_t CreateMidSideConnect(int side);
	
	/**
	 * @brief Checks if the side order is consistent with the preferred side order and
	 * with the constraints and recomputes the constraints if necessary
	 */
	/** Calls IdentifySideOrder on higher level (i.e. smaller) connected elements recursively */
	virtual void IdentifySideOrder(int side);
	
	/** @} */
	
	/**
	 * @name Data consistency methods
	 * @brief These methods are used during the debugging phase and check the consistency of the data structures
	 * @{
	 */
	
	/** @brief Check the consistency of the constrained connects along a side*/
	void CheckConstraintConsistency(int side);
	
	/** @brief Check the consistency of the constrained connects for all sides*/
	void CheckConstraintConsistency();
	
	/** @brief Checks element data structure consistancy */
	virtual  int CheckElementConsistency();
	
	/**
	 * @brief Compare the shape functions of sides of an element
	 * @param sides small side
	 * @param sidel large side
	 * @param phis small side function values
	 * @param dphis small side gradient function values
	 * @param phil large side function values
	 * @param dphil large side gradient function values
	 * @param transform transformation matrix from large side to small side
	 */
	int CompareShapeF(int sides, int sidel, TPZFMatrix<REAL> &phis, TPZFMatrix<REAL> &dphis, TPZFMatrix<REAL> &phil, TPZFMatrix<REAL> &dphil, TPZTransform<> &transform);
	
	/**
	 * @brief Returns the transformation which transform a point from the side to the interior of the element
	 * @param side side from which the point will be tranformed (0<=side<=2)
	 * @return TPZTransform<> object
	 * @see the class TPZTransform
	 */
	virtual TPZTransform<> TransformSideToElement(int side) = 0;
	
	/** @} */
	
public:
	
	/** @brief To enable to work with discontinuous element that can have interface elements*/
	virtual void SetInterface(int /*side*/, int64_t /*index*/) { }
	virtual int Interface(int /*side*/) { return -1; }
	virtual int CanHaveInterface() { return 0; }
	virtual void DeleteInterfaces() { }
	/** @brief Returns total mass contained into the element */
	REAL MeanSolution(int var);
	/** @brief Computes the integral over the finite element */
	void CalcIntegral(TPZElementMatrix &ef);
	
};

#endif

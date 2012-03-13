/**
 * @file
 * @brief Contains declaration of TPZCompEl class which defines the interface of a computational element.
 */
// $Id: pzcompel.h,v 1.47 2011-05-11 02:27:20 phil Exp $

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
#include "pzcreateapproxspace.h"
#include "pzmaterialdata.h"
//#include "pzmaterial.h"


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

class TPZMaterialData;

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

#include <set>

/**
 * @brief Defines the interface of a computational element. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
class TPZCompEl : public virtual TPZSaveable {
	
protected:
	
	/** @brief Computational mesh to which the element belongs */
	TPZCompMesh 	*fMesh;
	
	/** @brief Element index into mesh element vector */
	int fIndex;
	
private:	
	/** @brief Index of reference element */
	int fReferenceIndex;
	
public:
	
	/** @brief Simple Constructor */
	TPZCompEl();
	
	/** @brief Simple destructor */
	virtual ~TPZCompEl();
	
	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const = 0;
	
	/**
	 * @brief Method for creating a copy of the element in a patch mesh
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the computational elements
     */
	/**
	 * Otherwise of the previous clone function, this method don't
	 * copy entire mesh. Therefore it needs to map the connect index
	 * from the both meshes - original and patch
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
									std::map<int,int> & gl2lcConMap,
									std::map<int,int> & gl2lcElMap) const = 0;
	
	
	/** @brief Put a copy of the element in the referred mesh */
	TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy);
	
	/** @brief Put a copy of the element in the patch mesh */
	TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<int,int> &gl2lcElMap);
	
	/** @brief Copy of the element in the new mesh whit alocated index */
	TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, int &index);
	
	/**
	 * @brief Creates a computational element within mesh. Inserts the element within the data structure of the mesh
	 * @param mesh mesh wher will be created the element
	 * @param gel geometric element for which the computational element will be created
	 * @param index new elemen index
	 */
	TPZCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);
	
	/** @brief Sets & Gets the value of gOrder */
	static void SetgOrder( int order );
	
	static int GetgOrder();
	
	/** @brief Sets create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh){
		mesh->SetAllCreateFunctionsContinuous();
	}
	
	/** @brief Returns the volume of the geometric element associated. */
	virtual  REAL VolumeOfEl()
	{
		if(fReferenceIndex == 0) return 0.;
		return Reference()->Volume();
	}
	
	/** @brief Loads the geometric element referece */
	virtual void LoadElementReference();
	
	/**
	 * @brief This method verifies if the element has the given data characteristics
	 * @param var state variable number
	 * @param matname pointer to material name
	 **/
	virtual REAL CompareElement(int var, char *matname);
	
	/**
	 * @name Access
	 * @brief Access to private data
	 * @{
	 */
	
	/** @brief Return the type of the element */
	/* the types are listed in parameter : ENoType, EOned, ETriangle, EQuadrilateral, ESubstructure */
	virtual MElementType Type();
	
	virtual int IsInterface() { return 0; }
	
	/** @brief Return a pointer to the corresponding geometric element if such exists, return 0 otherwise */
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
	
	/** @brief Returns the number of nodes of the element */
	virtual int NConnects() const =0;
	
	/** @brief Returns the number of equations of the element */
	virtual int NEquations();
	
	/** @brief Returns element index of the mesh fELementVec list */
	int Index() const;
	
	/** @brief Sets element index of the mesh fELementVec list */
	void SetIndex(int index);
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int ConnectIndex(int i) const = 0;
	
	/**
	 * @brief Returns a pointer to the ith node
	 * @param i node index
	 */
	virtual TPZConnect &Connect(int i) const;
	
	/** @brief Dimension of the element */
	virtual int Dimension() const = 0;
	
	/** @brief Identify the material object associated with the element */
	virtual TPZAutoPointer<TPZMaterial> Material() const;
	
	/**
	 * Sets the material associated with the object
	 * param mat new element material
	 */
	//  virtual void SetMaterial(TPZAutoPointer<TPZMaterial> mat) = 0;
	
	/** @brief Returns the reference geometric element patch 
     * Look for a geometric element which refers to a computational element and
     * Is neighbour of the current element AND is larger than the current element
     */
	TPZGeoEl * GetRefElPatch();
	
	//void SetIntegrationRule(int order);
	/**
	* @}
	 */
	
	/**
	 * @name MODIFICATION_OF_PRIVATE_DATA
	 * @brief Methods that modify private data
	 * @{
	 */
	
	/**
	 * @brief Creates corresponding graphical element(s) if the dimension matches
	 * graphical elements are used to generate output files
	 * @param graphmesh graphical mesh where the element will be created
	 * @param dimension target dimension of the graphical element
	 */
	virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension);
	
	/**
	 * @brief Loads the solution within the internal data structure of the element
	 */ 
	/** Is used to initialize the solution of connect objects with dependency
	 * Is also used to load the solution within SuperElements
	 */
	virtual void LoadSolution();
	
	/**
	 * @brief Sets the grid of the element
	 * @param mesh new reference mesh
	 */
	void SetMesh(TPZCompMesh *mesh);
	
	/** @brief Return a pointer to the grid of the element */
	TPZCompMesh *Mesh() const;

	/**
	 * @}
	 */
	
	/**
	 * @name Print
	 * @brief Methods for print data structure
	 * @{
	 */
	
	/**
	 * @brief Prints element data
	 * @param out Indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const;
	
	/**
	 * @brief Output device operator
	 * @param s Indicates the device where the data will be printed
	 * @param el Element to print
	 */
	friend std::ostream& operator<<(std::ostream &s,TPZCompEl &el);
	
	/**
	 * @brief Prints the solution - sol - for the variable "VarName"
	 * at point specified in terms of the master element coordinates
	 * @param point master element coordinate to print
	 * @param VarName name of variable to print
	 * @param out indicates the device where the data will be printed
	 */
	virtual void PrintSolution(TPZVec<REAL> &point,char *VarName,std::ostream &out);
	
	/**
	 * @brief Prints one coordinate index corresponding to the point to the output stream
	 * @param point master element coordinate to print
	 * @param CoordinateIndex index of the coordinate corresponding to the point
	 * @param out indicates the device where the data will be printed
	 */
	virtual void PrintCoordinate(TPZVec<REAL> &point,int CoordinateIndex,std::ostream &out);
	
	/**
	 * @brief Prints the variables names associated with the element material
	 * @param VarName pointer to variable parameter wha want to print
	 * @param out indicates the device where the data will be printed
	 */
	virtual void PrintTitle(char *VarName,std::ostream &out);

	/**
	 * @}
	 */
	
	/**
	 * @brief Sets the orthogonal function which will be used throughout the program
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
	 * @brief Computes the element stifness matrix and right hand side
	 * @param ek element stiffness matrix
	 * @param ef element load vector
	 */
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
	
	/** @brief Verifies if the material associated with the element is contained in the set */
	virtual bool HasMaterial(const std::set<int> &materialids);
	
	/**
	 * @brief Computes the element right hand side
	 * @param ef element load vector(s)
	 */
	virtual void CalcResidual(TPZElementMatrix &ef);
	
	/**
	 * @brief Implements of the orthogonal Chebyshev functions
	 * @param x point where the Chebyshev function will be evaluated
	 * @param num number of functions
	 * @param phi values of the function
	 * @param dphi values of derivative of the function
	 */
	static void Chebyshev(REAL x,int num,TPZFMatrix &phi,TPZFMatrix &dphi);
	
	/**
	 * @brief Divide the computational element. If interpolate = 1, the solution is interpolated to the sub elements
	 * @param index  index of the element which is being divided
	 * @param subindex element vector where will be created the divided elements
	 * @param interpolate boolean variable to indicates if the solution will be interpolated to the sub elements
	 */
	/** This method needs to be implemented in the derived classes */
	virtual void Divide(int index, TPZVec<int> &subindex, int interpolate = 0);
	
	/**
	 * @brief Projects the flux function on the finite element space
	 * @param ek element stiffness matrix
	 * @param ef element loads matrix
	 */
	virtual void ProjectFlux(TPZElementMatrix &ek,TPZElementMatrix &ef);
	
	/**
	 * @brief Performs an error estimate on the elemen
	 * @param fp function pointer which computes the exact solution
	 * @param errors [out] the L2 norm of the error of the solution
	 * @param flux [in] value of the interpolated flux values
	 */
	virtual void EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
							   TPZVec<REAL> &errors,TPZBlock *flux);
	
	/** @brief ComputeError computes the element error estimator */
	virtual void ComputeError(int errorid, TPZVec<REAL> &error){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented.\n";
	}
	
	/** @brief Integrates a variable over the element. */
	virtual void Integrate(int variable, TPZVec<REAL> & value){
		value.Fill(0.);
		PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented.\n";
	}
	
	/**
	 * @brief Calculates the solution - sol - for the variable var
	 * at point qsi, where qsi is expressed in terms of the
	 * master element coordinates
	 * @param qsi master element coordinate
	 * @param var variable name
	 * @param sol vetor for the solution
	 */
	virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol);
	
	virtual void ComputeSolution(TPZManVector<REAL,10> &qsi, TPZMaterialData &data)	{
		std::cout <<"Imposed for Hdiv solution ";
		DebugStop();
	};
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix &axes) = 0;
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
	 * This method will function for both volumetric and interface elements
	 * @param qsi master element coordinate of the interface element
	 * @param normal vector
	 * @param leftsol finite element solution
	 * @param dleftsol solution derivatives
	 * @param leftaxes axes associated with the left solution
	 * @param rightsol finite element solution
	 * @param drightsol solution derivatives
	 * @param rightaxes axes associated with the right solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZVec<REAL> &normal,
								 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix &leftaxes,
								 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix &rightaxes) = 0;
	
	/**
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
	 * @param axes [in] axes indicating the direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
								 const TPZFMatrix &axes, TPZSolVec &sol, TPZGradSolVec &dsol) = 0;
	/**
	 * @brief Builds the list of all connectivities related to the element including the
	 * connects pointed to by dependent connects
	 * @param indepconnectlist set of independent connect indices
	 * @param depconnectlist set of dependent connect indices
	 * @note Note : this method does not reset the set to zero. The calling
	 * method should do this
	 */
	virtual void BuildConnectList(std::set<int> &indepconnectlist, std::set<int> &depconnectlist);
	/**
	 * @brief Builds the list of all connectivities related to the element including the
	 * connects pointed to by dependent connects
	 * @param connectlist stack to receive the list
	 * @note Note : this method does not reset the stack to zero. The calling
	 * method should do this
	 */
	virtual void BuildConnectList(TPZStack<int> &connectlist);
	/**
	 * @brief Builds the list of all connectivities related to the element including the
	 * connects pointed to by dependent connects
	 * @param connectlist stack to receive the list
	 * @note Note : this method ADDS connects to the set.
	 */
	virtual void BuildConnectList(std::set<int> &connectlist);
	
	/** @brief Returns 1 if the element has at least one dependent node. Returns 0 otherwise */
	virtual int HasDependency();
	
    /** @brief returns the index of the pressure connect
     * returns -1 if their is no pressure connect
     */
    virtual int PressureConnectIndex() const
    {
        return -1;
    }

	/**
	 * @brief Domain Decomposition.\n
	 * This method will eliminate the nodes which are internal to the element from
	 * the datastructure of the grid \n
	 * After calling this method, the superelement will statically condense the
	 * internal equations
	 */
	virtual void ReduceInternalNodes() { }
	
	/**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, int index) = 0;
	
	/**
	 * @brief Calculates the diagonal block
	 * @param connectlist stack list to calculates the diagonal block
	 * @param block object to receive the diagonal block
	 */
	virtual void CalcBlockDiagonal(TPZStack<int> &connectlist, TPZBlockDiagonal & block);
	
	REAL MaximumRadiusOfEl();
	
	REAL LesserEdgeOfEl();
	
	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
private:
	/** @brief Default interpolation order */
    static int gOrder;

};


class TPZGeoElSide;

/**
 * @brief Implements computational element and a side. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */

/**
 This class was created to implement all algorithms associated with element/sides. \n
 Objects of this class are mostly temporary
 */
class TPZCompElSide {
	
	/** @brief Pointer to the computational element */
	TPZCompEl *fEl;
	
	/** @brief Index of the object side */
	int fSide;
	
	
public:
	
	//    /*PARA TESTES*/ TPZGeoEl *georeftest;
	
	/** @brief Simple Constructor */
	TPZCompElSide();
	
	/**
	 * @brief Creates a computational element side from an object TPZCompElSide
	 * @param celside reference computational element side
	 */
	TPZCompElSide(const TPZCompElSide &celside);
	
	/**
	 * @brief Creates an computational element side given a computational element
	 * and the side index
	 * @param cel pointer to the computational element
	 * @param side index of the side
	 */
	TPZCompElSide(TPZCompEl *cel,int side);
    
    /// constructor which allows us to create a vector of objects
    TPZCompElSide(int zero) : fEl(0), fSide(-1)
    {
    }
	
	/** @brief Gives a pointer to the reference computational element */
	TPZCompEl *Element() const {return fEl;}
    
    /**
     * @brief The conversion to bool indicates whether the object has an associated element
     */
    operator bool() const
    {
        return fEl != 0;
    }

	
	/** @brief Sets computational element pointer. */
	void SetElement(TPZCompEl* el){ fEl = el;}
	
	/** @brief Returns the side index */
	int Side() const {return fSide;}
	
	/**
	 * @brief Sets the side index
	 * @param side new side index
	 */
	void SetSide(int side);
	
	/** @brief Verifies if the object is non null (initialized) */
	int Exists() const {return fEl != 0;}
	
	/** @brief Reference to the geometric element */
	TPZGeoElSide Reference() const;
	
	/**
	 * @brief Returns all connected elements which have level higher to the current element
	 * @param elsidevec side elements vector
	 * @param onlyinterpolated if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
	 * @param removeduplicates if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
	 */
	void HigherLevelElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
	
	/**
	 * @brief Pushes all element/sides which have higher dimension than the current element/side
	 * @param elsidevec side elements vector where elements/side will be put
	 * @param onlyinterpolated if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
	 * @param removeduplicates if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
	 */
	void HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
	/**
	 * @brief Returns all connected elements to the current element
	 * @param elsidevec side elements vector
	 * @param onlyinterpolated if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
	 * @param removeduplicates if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
	 */
	void ConnectedElementList(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated, int removeduplicates);
	
	/**
	 * @brief Returns all connected elements which have equal level to the current element\n
	 * This method will not put this on the stack
	 * @param elsidevec side elements vector
	 * @param onlyinterpolated  if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
	 * @param removeduplicates  if removeduplicates == 1 no elements which are direct neighbours will be put on the stack
	 */
	void EqualLevelElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
	
	/**
	 * @brief Returns all connected elements which have level lower to the current element
	 * @param onlyinterpolated if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack
	 */
	/** If removeduplicates == 1 no elements which are direct neighbours will be put on the stack */
	TPZCompElSide LowerLevelElementList(int onlyinterpolated);
	
	/**
	 * @brief Will remove elements which are direct neighbours from elvec (and elsides)
	 * @param elvec computational element side vector
	 */
	/**
	 * The method checks between any two elements in the list whether they are of equal level
	 * and whether they are neighbours. If they are neighbours, one of the elements will be removed
	 * from the list \n
	 * The method NeighbourExists between any two elements of equal level will return 0
	 */
	static void RemoveDuplicates(TPZStack<TPZCompElSide> &elvec);
	
	
	/**
	 * @brief Remove entries of the vector which share a connect along the side This should be
	 * equivalent to RemoveDuplicates.
	 * @param expandvec vector of TPZCompElSide objects
	 */
	static void RemoveConnectDuplicates(TPZStack<TPZCompElSide> &expandvec);
	
	/**
	 * @brief Find the list element/side of the current element restrict nodes and elements
	 * @param expandvec Vector of the PZCompElSide
	 * @param onlyinterpolated if ==1 only elements derived from TPZInterpolated will be put on the stack
	 */
	void ExpandConnected(TPZStack<TPZCompElSide> &expandvec,int onlyinterpolated);
	
	/**
	 * @brief Returns the element with lowest id of all direct neighbours of expandvec
	 * @param expandvec elements whose neighbours will be checked
	 * @param onlyinterpolated if == 1, only elements derived from TPZInterpolatedElement will be checked
	 */
	TPZCompElSide LowerIdElementList(TPZCompElSide &expandvec,int onlyinterpolated);
	
	//inline//
	
	/** @brief Returns the index of the middle side connect alon fSide */
    int ConnectIndex() const;
	
	bool operator != (const TPZCompElSide &other);
	bool operator == (const TPZCompElSide &other);
	
	
};
//  std::ostream & operator << (std::ostream &out,const TPZCompElSide &celside);

inline void TPZCompEl::CreateGraphicalElement(TPZGraphMesh &, int) {
	std::cout << "TPZCompEl::CreateGrafEl called\n";
	this->Print(std::cout);
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

/** @brief Overload operator << to write computational element side data */
inline std::ostream &operator << (std::ostream &out,const TPZCompElSide &celside)
{
	out << "Side = " << celside.Side()
	<< " element: " << celside.Element()->Index()
	<< std::endl;
	return out;
}

inline int TPZCompEl::Index() const {
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

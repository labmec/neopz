/**
 * @file
 * @brief Contains declaration of TPZCompEl class which defines the interface of a computational element.
 */

#ifndef COMPELEMHPP
#define COMPELEMHPP

#include "pzreal.h"
#include <iostream>
#include <fstream>
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "TPZSavable.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzblock.h"
#include "pzblockdiag.h"
#include "pzcreateapproxspace.h"
#include "TPZOneShapeRestraint.h"
//#include "pztransfer.h"
#include <functional>


struct TPZElementMatrix;
template<class TVar>
struct TPZElementMatrixT;
class TPZCompMesh;
class TPZBndCond;
class TPZInterpolatedElement;
class TPZInterfaceElement;
class TPZConnect;
class TPZMaterial;
class TPZGeoEl;
class TPZGeoNode;

template<class T>
class TPZVec;
template<class T, int N>
class TPZManVector;
template<class T, int N>
class TPZStack;

class TPZGraphMesh;
class TPZIntPoints;

template<class T>
class TPZTransform;

#include "pzeltype.h"

#include <set>

/**
 * @brief Defines the interface of a computational element. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
class TPZCompEl : public virtual TPZSavable {
	
protected:
	
	/** @brief Computational mesh to which the element belongs */
	TPZCompMesh 	*fMesh;
	
	/** @brief Element index into mesh element vector */
	int64_t fIndex;
	
private:	
	/** @brief Index of reference element */
	int64_t fReferenceIndex;
	
public:
	
    static int StaticClassId();
    
    virtual int ClassId() const override;

    
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
									std::map<int64_t,int64_t> & gl2lcConMap,
									std::map<int64_t,int64_t> & gl2lcElMap) const = 0;
	
	
	/** @brief Put a copy of the element in the referred mesh */
	TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy);
	
	/** @brief Put a copy of the element in the patch mesh */
	TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<int64_t,int64_t> &gl2lcElMap);
		
	/**
	 * @brief Creates a computational element within mesh. Inserts the element within the data structure of the mesh
	 * @param mesh mesh wher will be created the element
	 * @param gel geometric element for which the computational element will be created
	 * @param index new elemen index
	 */
	TPZCompEl(TPZCompMesh &mesh, TPZGeoEl *gel);
	
	/** @brief Sets the value of the default interpolation order */
	static void SetgOrder( int order );
	
    /** @brief Set the default value of the interpolation order */
	static int GetgOrder();
	
	/** @brief Sets create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh);
	
	/** @brief Returns the volume of the geometric element associated. */
	virtual  REAL VolumeOfEl()
	{
		if(fReferenceIndex == -1) return 0.;
		return Reference()->Volume();
	}
	
	/** @brief Loads the geometric element reference */
	virtual void LoadElementReference();
	
	/**
	 * @brief This method computes the norm of the difference of a post processed variable with
     @ the same post processed variable of the element pointed to by the geometric element
	 * @param var state variable number
	 * @param matname only contribute to the norm if the material name matches the name of the material of the element
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
	TPZGeoEl *Reference() const;

    int64_t ReferenceIndex()
    {
        return fReferenceIndex;
    }
    
	void SetReference(int64_t referenceindex)
	{
		fReferenceIndex = referenceindex;
	}
	
    /// return true if the element has a variational statement associated with the material ids
    virtual bool NeedsComputing(const std::set<int> &materialids)
    {
        TPZGeoEl *gel = Reference();
        if (!gel) {
            DebugStop();
        }
        return materialids.find(gel->MaterialId()) != materialids.end();
    }
	/** @brief Returns the number of nodes of the element */
	virtual int NConnects() const = 0;
	
	/** @brief Returns the number of equations of the element */
	virtual int NEquations();
	
	/** @brief Returns element index of the mesh fELementVec list */
	int64_t Index() const;
	
	/** @brief Sets element index of the mesh fELementVec list */
	void SetIndex(int64_t index);
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int64_t ConnectIndex(int i) const = 0;
	
	/**
	 * @brief Returns a pointer to the ith node
	 * @param i node index
	 */
	virtual TPZConnect &Connect(int i) const;
	
	/** @brief Dimension of the element */
	virtual int Dimension() const = 0;
	
	/** @brief Identify the material object associated with the element */
	virtual TPZMaterial * Material() const;

	/** 
	 * @brief Returns the reference geometric element patch. \n
     * Look for a geometric element which refers to a computational element and
     * is neighbour of the current element AND is larger than the current element
     */
	TPZGeoEl * GetRefElPatch();
	
	/** @} */
	
	/**
	 * @name MODIFICATION_OF_PRIVATE_DATA
	 * @brief Methods that modify private data
	 * @{
	 */
	
	/** @brief Loads the solution within the internal data structure of the element */ 
	/** 
	 * Is used to initialize the solution of connect objects with dependency. \n
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
	 * @name Print methods
	 * @brief Methods for print data structure
	 * @{
	 */

	/**
	 * @brief Prints element data
	 * @param out Indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const;

    /**
     * @brief ShortPrint element data
     * @param out Indicates the device where the data will be printed
     * @note currently (04/03/20) this method is being implemented only for pzmultiphysicscompel::ShortPrint
     * @note If a short print is useful for your application, please implement as you see fit and modify these notes.
     */
    virtual void ShortPrint(std::ostream &out = std::cout) const;
	
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
	virtual void PrintSolution(TPZVec<REAL> &point, const char *VarName,std::ostream &out);
	
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
	virtual void PrintTitle(const char *VarName,std::ostream &out);

	/** @} */
	
	/**
	 * @brief Sets the orthogonal function which will be used throughout the program
	 * by default this function is the Chebyshev function
	 * @param orthogonal pointer to a function which will be used to generate the shape functions
	 */
	static void SetOrthogonalFunction(void (*orthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,
														 TPZFMatrix<REAL> & dphi));
    /**
     * @brief Computes the element stifness matrix and right hand side
     * in an internal data structure. Used for initializing condensed element data structures
     */
    virtual void Assemble();
    
    
    /** @brief Initialize element matrix in which is computed CalcStiff */
    virtual void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
    
    /** @brief Initialize element matrix in which is computed in CalcResidual */
    virtual void InitializeElementMatrix(TPZElementMatrix &ef);
    

	/**
	 * @brief Computes the element stifness matrix and right hand side
	 * @param ek element stiffness matrix
	 * @param ef element load vector
	 */
	virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef);
    virtual void CalcStiff(TPZElementMatrixT<CSTATE> &ek,TPZElementMatrixT<CSTATE> &ef);
    
	
	/** @brief Verifies if the material associated with the element is contained in the set */
	virtual bool HasMaterial(const std::set<int> &materialids) const;
	
	/**
	 * @brief Computes the element right hand side
	 * @param ef element load vector(s)
	 */
	virtual void CalcResidual(TPZElementMatrixT<STATE> &ef);
    virtual void CalcResidual(TPZElementMatrixT<CSTATE> &ef);
	
	/**
	 * @brief Implements of the orthogonal Chebyshev functions
	 * @param x point where the Chebyshev function will be evaluated
	 * @param num number of functions
	 * @param phi values of the function
	 * @param dphi values of derivative of the function
	 */
	static void Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	/**
	 * @brief Divide the computational element. If interpolate = 1, the solution is interpolated to the sub elements
	 * @param index  index of the element which is being divided
	 * @param subindex element vector where will be created the divided elements
	 * @param interpolate boolean variable to indicates if the solution will be interpolated to the sub elements
	 */
	/** This method needs to be implemented in the derived classes */
	virtual void Divide(int64_t index, TPZVec<int64_t> &subindex, int interpolate = 0);
	

  /**
	 * @brief Performs an error estimate on the element based on materials solution
	 * @param errors [out] the L2 norm of the error of the solution
	 * @param flux [in] value of the interpolated flux values
	 */
    virtual void EvaluateError(TPZVec<REAL> &errors, bool store_error);

	/** @brief ComputeError computes the element error estimator */
	virtual void ComputeError(int errorid, TPZVec<REAL> &error){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented.\n";
	}
	
	/**
	 * @brief Creates corresponding graphical element(s) if the dimension matches
	 * graphical elements are used to generate output files
	 * @param graphmesh graphical mesh where the element will be created
	 * @param dimension target dimension of the graphical element
	 */
	virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension);
	
	/** @brief Integrates a variable over the element. */
	virtual void Integrate(int variable, TPZVec<STATE> & value){
		value.Fill(0.);
		PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented.\n";
	}
    
    /** @brief Get the indices of the vector of element memory associated with the integration points */
    /**
     * Will return an empty vector if no memory is associated with the integration point
     * Is implemented in TPZCompElWithMem
     */
    virtual void GetMemoryIndices(TPZVec<int64_t> &indices) const
    {
        indices.resize(0);
    }
    
    /** @brief Set the indices of the vector of element memory associated with the integration points */
    /**
     * Will return an empty vector if no memory is associated with the integration point
     * Is implemented in TPZCompElWithMem
     */
    virtual void SetMemoryIndices(TPZVec<int64_t> &indices)
    {
        if(indices.size() != 0)
        {
            DebugStop();
        }
    }
	
    /** @brief Prepare the vector of the material withmem with the correct integration point indexes */
    virtual void PrepareIntPtIndices(){
        
    }
  
  /** @brief PrepareIntPtIndices initializes the material damage varibles memory in the proper material class. */
	virtual void ForcePrepareIntPtIndices(){
    
  }
  
  /** @brief Frees the material damage varibles memory in the proper material class. */
	virtual void SetFreeIntPtIndices(){
    
  }


	/** @brief Return the size of the elementvec in multiphysics, if it is not multiphysics, just return 1 */
	virtual int NumberOfCompElementsInsideThisCompEl(){
		return 1;
	}
    
    virtual void TransferMultiphysicsElementSolution()
    {
        // Nothing to be done here
    }
    
    virtual void SetMultiphysicsElementSolution()
    {
        DebugStop();
    }
	
    /// Add a shape restraint (meant to fit the pyramid to restraint
    virtual void AddShapeRestraint(TPZOneShapeRestraint restraint)
    {
        DebugStop();
    }

    /// Return a list with the shape restraints
    virtual std::list<TPZOneShapeRestraint> GetShapeRestraints() const
    {
        std::list<TPZOneShapeRestraint> loc;
        return loc;
    }

    /// Return a list with the shape restraints
    virtual void ResetShapeRestraints()
    {
    }
	/**
	 * @brief Calculates the solution - sol - for the variable var
	 * at point qsi, where qsi is expressed in terms of the
	 * master element coordinates
	 * @param qsi master element coordinate
	 * @param var variable name
	 * @param sol vetor for the solution
	 */
	virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol){
        SolutionInternal(qsi,var,sol);
    }
    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<CSTATE> &sol){
        SolutionInternal(qsi,var,sol);
    }
    
    /**
     * @brief Compute the integral of a variable
     */
    virtual TPZVec<STATE> IntegrateSolution(int var) const;
    
    /**
     * @brief Compute the integral of a variable defined by the string if the material id is included in matids
     */
    virtual TPZVec<STATE> IntegrateSolution(const std::string &varname, const std::set<int> &matids);
    //@}
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const = 0;

	/**
	 * @brief Builds the list of all connectivities related to the element including the
	 * connects pointed to by dependent connects
	 * @param indepconnectlist set of independent connect indices
	 * @param depconnectlist set of dependent connect indices
	 * @note Note : this method does not reset the set to zero. The calling
	 * method should do this
	 */
//protected:
	virtual void BuildConnectList(std::set<int64_t> &indepconnectlist, std::set<int64_t> &depconnectlist);
public:
	/**
	 * @brief Builds the list of all connectivities related to the element including the
	 * connects pointed to by dependent connects
	 * @param connectlist stack to receive the list
	 * @note Note : this method does not reset the stack to zero. The calling
	 * method should do this
	 */
	virtual void BuildConnectList(TPZStack<int64_t> &connectlist) const;
	/**
	 * @brief Builds the list of all connectivities related to the element including the
	 * connects pointed to by dependent connects
	 * @param connectlist stack to receive the list
	 * @note Note : this method ADDS connects to the set.
	 */
//protected:
	virtual void BuildConnectList(std::set<int64_t> &connectlist);
public:
	
	/** @brief Returns 1 if the element has at least one dependent node. Returns 0 otherwise */
	virtual int HasDependency();
	
    /** 
	 * @brief Returns the index of the pressure connect
     * @note Returns -1 if their is no pressure connect
     */
    virtual int PressureConnectIndex() const;

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
	virtual void SetConnectIndex(int inode, int64_t index) = 0;
	
	/**
	 * @brief Calculates the diagonal block
	 * @param connectlist stack list to calculates the diagonal block
	 * @param block object to receive the diagonal block
	 */
	virtual void CalcBlockDiagonal(TPZStack<int64_t> &connectlist, TPZBlockDiagonal<STATE> & block){
        CalcBlockDiagonalInternal(connectlist,block);
    }
    
    virtual void CalcBlockDiagonal(TPZStack<int64_t> &connectlist, TPZBlockDiagonal<CSTATE> & block){
        CalcBlockDiagonalInternal(connectlist,block);
    }
	
    /// Will return the maximum distance between the nodes of the reference element
	REAL MaximumRadiusOfEl();
	
    /// Will return the smallest distance between two nodes of the reference element
	REAL LesserEdgeOfEl();
	
	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;
	 
private:
    template<class TVar>
    void LoadSolutionInternal();
    template<class TVar>
    void SolutionInternal(TPZVec<REAL> &qsi,int var,TPZVec<TVar> &sol);
    template<class TVar>
    void CalcBlockDiagonalInternal(TPZStack<int64_t> &connectlist, TPZBlockDiagonal<TVar> & block);
	/** @brief Default interpolation order */
    static int gOrder;
    
protected:
    
    /// Integration rule established by the user
    TPZIntPoints *fIntegrationRule;
    
public:
    
    virtual void InitializeIntegrationRule(){
        std::cout << "TPZCompEl::InitializeIntegrationRule should not be called\n";
        DebugStop();
    }
    
    
    virtual int ComputeIntegrationOrder() const;
    
    /// Method to set a dynamically allocated integration rule
    virtual void SetIntegrationRule(TPZIntPoints *intrule);
    
    virtual void SetIntegrationRule(int order) {
        std::cout << "TPZCompEl::SetIntegrationRule should not be called\n";
    }
    
    virtual const TPZIntPoints &GetIntegrationRule() const
    {
        if(fIntegrationRule)
			return *fIntegrationRule;
        DebugStop();
		return *fIntegrationRule;
    }


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

	/** @brief Simple Constructor */
    TPZCompElSide() : fEl(0), fSide(-1)
    {
        
    }
	
	/**
	 * @brief Creates a computational element side from an object TPZCompElSide
	 * @param celside reference computational element side
	 */
    TPZCompElSide(const TPZCompElSide &celside) : fEl(celside.fEl), fSide(celside.fSide)
    {
        
    }
	
	/**
	 * @brief Creates an computational element side given a computational element
	 * and the side index
	 * @param cel pointer to the computational element
	 * @param side index of the side
	 */
    TPZCompElSide(TPZCompEl *cel,int side) : fEl(cel), fSide(side)
    {
        
    }
    
    /** @brief Constructor which allows us to create a vector of objects */
    TPZCompElSide(int zero) : fEl(0), fSide(-1)
    {
    }
	
	/** @brief Gives a pointer to the reference computational element */
	TPZCompEl *Element() const {return fEl;}
    
    /** @brief The conversion to bool indicates whether the object has an associated element */
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
    
    /**
     * @brief Split the connect shared between this and right compelsides
     * @param right compelside that shares the same connects as in this
     */
    void SplitConnect(const TPZCompElSide& right) const;

	/** @brief Returns the index of the middle side connect alon fSide */
    int64_t ConnectIndex() const;
	
	/** @brief Overlapping operator not equal */
	bool operator != (const TPZCompElSide &other);
	/** @brief Overlapping operator equal */
	bool operator == (const TPZCompElSide &other);

};

inline void TPZCompEl::CreateGraphicalElement(TPZGraphMesh &, int) {
	std::cout << "TPZCompEl::CreateGraphicalElement called\n";
	this->Print(std::cout);
}

inline void TPZCompEl::Assemble(){
    std::cout << "TPZCompEl::Assemble is called." << std::endl;
}

inline void TPZCompEl::CalcStiff(TPZElementMatrixT<STATE> &,TPZElementMatrixT<STATE> &){
	PZError << "TPZCompEl::CalcStiff(*,*) is called." << std::endl;
    DebugStop();
}

inline void TPZCompEl::CalcStiff(TPZElementMatrixT<CSTATE> &,TPZElementMatrixT<CSTATE> &){
	PZError << "TPZCompEl::CalcStiff(*,*) is called." << std::endl;
    DebugStop();
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

inline int64_t TPZCompEl::Index() const {
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


#define INSTANTIATE(TVar) \
extern template \
void TPZCompEl::SolutionInternal<TVar>(TPZVec<REAL> &qsi,int var,TPZVec<TVar> &sol); \
extern template \
void TPZCompEl::CalcBlockDiagonalInternal<TVar>(TPZStack<int64_t> &connectlist, TPZBlockDiagonal<TVar> & block);

INSTANTIATE(STATE)
INSTANTIATE(CSTATE)
#undef INSTANTIATE

#endif

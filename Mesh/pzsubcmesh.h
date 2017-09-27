/**
 * @file
 * @brief Contains declaration of TPZSubCompMesh class which implements a group of computational elements as a mesh and an element.
 */

#ifndef SUBCMESH_H
#define SUBCMESH_H

#if _MSC_VER > 1000
#pragma once
#endif

#include <stdlib.h>

#include "pzcompel.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzreal.h"
#include "pzanalysis.h"

#include <list>

class TPZSubMeshFrontalAnalysis;
class TPZSubMeshAnalysis;
class TPZAnalysis;
class TPZGuiInterface;

/**
 * @brief Implements a group of computational elements as a mesh and an element. \ref CompMesh "Computational Mesh"
 * @ingroup CompMesh
 * @ingroup CompElement
 */
/** Class TPZSubCompMesh is derived from Computacional mesh and
 * computacional element classes.
 */
class TPZSubCompMesh :
public TPZCompMesh,
public TPZCompEl
{
protected:
	/** @brief Pointer to submesh analysis object. Defines the resolution type. */
	TPZAutoPointer<TPZAnalysis> fAnalysis;
	
	/** @brief Pointer to external location index of the connection. */ 
	/** If the connection hasn't external location return the local id. */
	TPZManVector<long> fConnectIndex;
	
	/** @brief Indexes of the external connections. */ 
	/** If the connection isn't external id is -1! */
	TPZManVector<long> fExternalLocIndex;
	/** @brief Maps indicating the correspondence between the connect index of the father mesh and de local connect id */
	std::map<long,long> fFatherToLocal;
    
    /// Number of rigid body modes expected by the internal matrix inversion
    long fSingularConnect;
    
    /// List of connect/degree of freedom which are used as Lagrange multipliers
    std::list<std::pair<long,int> > fLagrangeEquations;
	
private:
	/** @brief Transfers one element from a submesh to another mesh. */
	long TransferElementTo(TPZCompMesh * mesh, long elindex);
	/** @brief Transfers one element from a specified mesh to the current submesh. */
	long TransferElementFrom(TPZCompMesh *mesh, long elindex);
	
	/** @brief Marks the connect to be local */
	void MakeInternalFast(long local);
    
	/**
	 * @brief Transfer the dependency list of a connect. This will
	 * make the dependency disappear for the corresponding father mesh
	 */
	/** It is necessary that the number of elements connected to the connect be equal one */
	void TransferDependencies(long local);
	
public:
	/**
	 * @brief Constructor.
	 * @param mesh reference mesh
	 * @param index reference mesh element index to transfer to submesh
	 */
	TPZSubCompMesh(TPZCompMesh &mesh, long &index);
	/** @brief Default constructor */
	TPZSubCompMesh();
	/** @brief Destructor. */
	virtual ~TPZSubCompMesh();
	
	virtual int NConnectShapeF(int inod, int order){
		PZError << "\nPLEASE IMPLEMENT ME: " << __PRETTY_FUNCTION__ << "\n";
		return 0;
	}
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		std::cout << "TPZSubCompMesh::Clone should be implemented\n";
		return 0;
	}
	
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
									std::map<long,long> & gl2lcConMap,
									std::map<long,long> & gl2lcElMap) const
	{
		std::cout << "TPZSubCompMesh::Clone should be implemented\n";
		return 0;
	}
	
	/** @brief Sets the analysis type. */
	void SetAnalysisFrontal(int numThreads, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /**
     * @brief Condense the internal equations using a skyline symetric matrix 
     * the preconditioned argument indicates whether the equations are condensed with a direct method (0) or 
     * with a GMRes solver preconditioned by the decomposed matrix
     */
	void SetAnalysisSkyline(int numThreads, int preconditioned, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    /**
     * @brief Condense the internal equations using a skyline symetric matrix
     * the preconditioned argument indicates whether the equations are condensed with a direct method (0) or
     * with a GMRes solver preconditioned by the decomposed matrix
     */
    void SetAnalysisSkyline(int numThreads, int preconditioned, TPZAutoPointer<TPZRenumbering> renumber);
    
	TPZAutoPointer<TPZAnalysis> Analysis()
	{
		return fAnalysis;
	}
	
	/** @brief This method will load the elements of the mesh in their corresponding geometric elements */
	virtual void LoadElementReference();
	
	/**
	 * @brief This method will initiate the comparison between the current computational
	 * mesh and the mesh which is referenced by the geometric mesh
	 * @param var state variable number
	 * @param matname pointer to material name
	 */
	virtual REAL CompareElement(int var, char *matname);
	
	/**
	 * @brief Verifies the transfer possibility of the connection elindex from
	 * the mesh to the current submesh.
	 * @param mesh  pointer to given mesh
	 * @param elindex given mesh element index
	 */
	int IsAllowedElement(TPZCompMesh *mesh, long elindex);
	
	/** @brief Computes the number of internal equations */
	long NumInternalEquations();
	
	/**
	 * @brief This method computes the skyline of the system of equations
	 * @param skyline vector where the skyline will be computed
	 */
	virtual void SkylineInternal(TPZVec<long> &skyline);
	
	/// Puts the nodes which can be transferred in an ordered list
	void PotentialInternal(std::list<long> &connectindices) const;

	/**
	 * @brief Changes an local internal connection to a external connection in the father mesh.
	 * @param local makes the connect with index local an external node
	 */
	void MakeExternal(long local);
	
	/**
	 * @name Mesh
	 * @brief Virtual methods derived from TPZCompMesh
	 * @{
	 */

	/** @brief Transfer one element form a submesh to another mesh. */
	virtual long TransferElement(TPZCompMesh *mesh, long elindex);
	
	/**
	 * @brief Makes all mesh connections internal mesh connections.
	 * @note This method is not working right!!! See comments in pzsubcmesh.cpp
	 */
	virtual void MakeAllInternal();
	
	/**
	 * @brief Returns the rootmesh who have the specified connection.
	 * @param local connection local index
	 */
	virtual TPZCompMesh * RootMesh(long local);
	
	/**
	 * @brief Makes a specified connection a internal mesh connection.
	 * @param local connection local number to be processed
	 */
	virtual void MakeInternal(long local);
	
	/**
	 * @brief Puts an local connection in the supermesh - Supermesh is one
	 * mesh who contains the analised submesh.
	 * @param local local index of the element to be trasnfered
	 * @param super pointer to the destination mesh
	 */
	virtual long PutinSuperMesh(long local, TPZCompMesh *super);
	
	/**
	 * @brief Gets an external connection from the supermesh - Supermesh is one
	 * mesh who contains the analised submesh.
	 * @param superind index of the element to be trasnfered
	 * @param super pointer to the destination mesh
	 */
	virtual long GetFromSuperMesh(long superind, TPZCompMesh *super);
	
	/**
	 * @brief Gives the commom father mesh of the specified mesh and the current submesh.
	 * @param mesh pointer to other mesh whosw want to know the commom father mesh
	 */
	virtual TPZCompMesh * CommonMesh (TPZCompMesh *mesh);
	
	/** @brief Returns the current submesh father mesh. */
	virtual TPZCompMesh * FatherMesh() const;

	/** @} */
	
	/** @brief Optimizes the connections positions on block. */
	void PermuteExternalConnects();
	
	/**
	 * @brief Computes the permutation vector which puts the internal connects to the first on the list \n
	 * Respect the previous order of the connects
	 */
	void ComputePermutationInternalFirst(TPZVec<long> &permute) const;
	
	/**
	 * @brief Permutes the potentially internal connects to the first on the list \n
	 * Respect the previous order of the connects
	 */
	void PermuteInternalFirst(TPZVec<long> &permute);
	
	/**
	 * @brief Prints the submesh information on the specified device/file out. 
	 * @param out indicates the device where the data will be printed
	 */ 
	/** This method use the virtual method from Computacional Mesh class. */
	virtual void Print(std::ostream &out = std::cout) const;
	
	/** @brief Verifies if any element needs to be computed corresponding to the material ids */
	bool NeedsComputing(const std::set<int> &matids);
	/** @brief Virtual Method to allocate new connect */
	virtual long AllocateNewConnect(int nshape, int nstate, int order);
	
	/** @brief Virtual Method to allocate new connect */
	virtual long AllocateNewConnect(const TPZConnect &connect);
	
	/** @brief Gives the id node  of one local node in containing mesh. */
	long NodeIndex(long nolocal, TPZCompMesh *super);
    
    /// return the index in the subcompmesh of a connect with index within the father
    long InternalIndex(long IndexinFather);
    	
	/** @brief Virtual Method! See declaration in TPZCompEl class. */ 
	/** The use of this method in submesh class return -1 == Error! */
	virtual int Dimension() const;
	
	/** @brief Computes the number of elements connected to each connect object */
	virtual void ComputeNodElCon();
	
	/** @brief Computes the number of elements connected to each connect object */
	virtual void ComputeNodElCon(TPZVec<int> &nelconnected) const;

	/**
	 * @name Element
	 * @brief Virtual methods derived from TPZCompEl
	 * @{
	 */

	/**
	 * @brief Changes the current node index -inode- to the specified node- index.
	 * @param inode node index to change
	 * @param index new node index for connect
	 */
	virtual void SetConnectIndex(int inode, long index);
	
  	/** @brief Calculates the submesh stiffness matrix */
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
		
	/**
	 * @brief Creates corresponding graphical element(s) if the dimension matches
	 * graphical elements are used to generate output files
	 * @param graphmesh graphical mesh where the element will be created
	 * @param dimension target dimension of the graphical element
	 */
	virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension);

	/** @brief Returns the connection index i. */
	virtual long ConnectIndex(int i) const;
	
	/** @brief Returns the number of connections. */
	virtual int NConnects() const;
	
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<long> &connectindexes) const;

    /**
	 * @brief Load the father mesh solution to all submesh connects -
  	 * (internal and external).
	 */
	virtual void LoadSolution();
    
    /**
     * @brief Transfer multiphysics solution
     */
    virtual void TransferMultiphysicsElementSolution();
    
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
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
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
	 * @param axes [in] indicating the direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
								 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
  
  virtual void EvaluateError(  void (*fp)(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv),
                                     TPZVec<REAL> &errors,TPZBlock<REAL> * /*flux */);
	
	/** @} */
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	static int ClassId();
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
    /// Method to verify that the datastructures are consistent
	bool VerifyDatastructureConsistency();
    
    /// Set the connect/degree of freedom as Lagrange variable
    void AddLagrangeDOF(long connectindex, int idf)
    {
        fLagrangeEquations.push_back(std::pair<long,int>(connectindex,idf));
        SetNumberRigidBodyModes(fLagrangeEquations.size());
    }
    
    /** @brief Set the number of rigid body modes associated with the internal degrees of freedom */
    void SetNumberRigidBodyModes(int nrigid, unsigned char lagrange = 0);
    
    /** @brief Return the number of rigid body modes associated with the internal degrees of freedom */
    int NumberRigidBodyModes();

	/** @brief Static function for validation tests. */
	static int main();
};

#endif

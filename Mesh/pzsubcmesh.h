/**
 * @file
 * @brief Contains declaration of TPZSubCompMesh class which implements a group of computational elements as a mesh and an element.
 */

#ifndef SUBCMESH_H
#define SUBCMESH_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <stdlib.h>
#include "pzcompel.h"
#include "pzcmesh.h"
//#include "pzsmanal.h"
#include "pzvec.h"
#include "pzreal.h"
#include "pzanalysis.h"

class TPZSubMeshFrontalAnalysis;
class TPZSubMeshAnalysis;
class TPZAnalysis;
class TPZGuiInterface;

/**
 * @brief Implements a group of computational elements as a mesh and an element. \ref CompMesh "Computational Mesh"
 * @ingroup CompMesh
 * @ingroup CompElement
 */
/** Class TPZSubCompMesh derived from Computacional mesh and
 * computacional element classes.
 */
class TPZSubCompMesh :
public TPZCompMesh,
public TPZCompEl
{
protected:
	/**
	 * @brief Pointer to submesh analysis object. Defines the resolution type.
	 */
	TPZAutoPointer<TPZAnalysis> fAnalysis;
	
	/**
	 * @brief Pointer to external location index of the connection.
	 * 
	 * If the connection hasn't external location return the local id.
	 */
	TPZManVector<int> fConnectIndex;
	
	/**
	 * @brief Indexes of the external connections.
	 * 
	 * If the connection isn't external id is -1!
	 */
	TPZManVector<int> fExternalLocIndex;
	/**
	 * @brief Maps indicating the correspondence between the connect index of the father mesh and de local connect id
	 */
	std::map<int,int> fFatherToLocal;
    
    /// Number of rigid body modes expected by the internal matrix inversion
    int fSingularConnect;
	
	
private:
	/** @brief Transfers one element from a submesh to another mesh. */
	int TransferElementTo(TPZCompMesh * mesh, int elindex);
	/** @brief Transfers one element from a specified mesh to the current submesh. */
	int TransferElementFrom(TPZCompMesh *mesh, int elindex);
	
	/** @brief Marks the connect to be local */
	void MakeInternalFast(int local);

	/**
	 * @brief Transfer the dependency list of a connect. This will
	 * make the dependency disappear for the corresponding father mesh
	 * 
	 * It is necessary that the number of elements connected to the connect be equal one
	 */
	void TransferDependencies(int local);
	
public:
	/**
	 * @brief Constructor.
	 * @param mesh reference mesh
	 * @param index reference mesh element index to transfer to submesh
	 */
	TPZSubCompMesh(TPZCompMesh &mesh, int &index);
	TPZSubCompMesh();
	
	/**
	 * @brief Destructor.
	 */
	virtual ~TPZSubCompMesh();
	virtual int NConnectShapeF(int inod){
		PZError << "\nPLEASE IMPLEMENT ME: " << __PRETTY_FUNCTION__ << "\n";
		return 0;
	}
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		std::cout << "TPZSubCompMesh::Clone should be implemented\n";
		return 0;
	}
	
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
									std::map<int,int> & gl2lcConMap,
									std::map<int,int> & gl2lcElMap) const
	{
		std::cout << "TPZSubCompMesh::Clone should be implemented\n";
		return 0;
	}
	
	/**
	 * @brief Static function for validation tests.
	 */
	static int main();
	
	/**
	 * @brief Sets the analysis type.
	 */
	void SetAnalysisFrontal(int numThreads, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /**
     * @brief Condense the internal equations using a skyline symetric matrix 
     * the preconditioned argument indicates whether the equations are condensed with a direct method (0) or 
     * with a GMRes solver preconditioned by the decomposed matrix
     */
	void SetAnalysisSkyline(int numThreads, int preconditioned, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	TPZAutoPointer<TPZAnalysis> Analysis()
	{
		return fAnalysis;
	}
	
	/**
	 * @brief This method will load the elements of the grid in their corresponding geometric
	 * elements
	 **/
	virtual void LoadElementReference();
	
	/**
	 * @brief This method will initiate the comparison between the current computational
	 * mesh and the mesh which is referenced by the geometric mesh
	 * @param var state variable number
	 * @param matname pointer to material name
	 **/
	virtual REAL CompareElement(int var, char *matname);
	
	/**
	 * @brief Verifies the transfer possibility of the connection elindex from
	 * the mesh to the current submesh.
	 * @param mesh  pointer to given mesh
	 * @param elindex given mesh element index
	 */
	int IsAllowedElement(TPZCompMesh *mesh, int elindex);
	
	/**
	 * @brief Computes the number of internal equations
	 */
	int NumInternalEquations();
	
	/**
	 * @brief This method computes the skyline of the system of equations
	 * @param skyline vector where the skyline will be computed
	 */
	virtual void SkylineInternal(TPZVec<int> &skyline);
	
	/// Puts the nodes which can be transferred in an ordered list
	void PotentialInternal(std::list<int> &connectindices) const;

	/**
	 * @brief Changes an local internal connection to a external connection in the father mesh.
	 * @param local makes the connect with index local an external node
	 */
	void MakeExternal(int local);
	
	/**
	 * @name Mesh
	 * @brief Virtual methods derived from TPZCompMesh
	 * @{
	 */

	/**
	 * @brief Transfer one element form a submesh to another mesh.
	 */
	virtual int TransferElement(TPZCompMesh *mesh, int elindex);
	
	/**
	 * @brief Makes all mesh connections internal mesh connections.
	 * @note This method is not working right!!! See comments in pzsubcmesh.cpp
	 */
	virtual void MakeAllInternal();
	
	/**
	 * @brief Returns the rootmesh who have the specified connection.
	 * @param local connection local index
	 */
	virtual TPZCompMesh * RootMesh(int local);
	
	/**
	 * @brief Makes a specified connection a internal mesh connection.
	 * @param local connection local number to be processed
	 */
	virtual void MakeInternal(int local);
	
	/**
	 * @brief Puts an local connection in the supermesh - Supermesh is one
	 * mesh who contains the analised submesh.
	 * @param local local index of the element to be trasnfered
	 * @param super pointer to the destination mesh
	 */
	virtual int PutinSuperMesh(int local, TPZCompMesh *super);
	
	/**
	 * @brief Gets an external connection from the supermesh - Supermesh is one
	 * mesh who contains the analised submesh.
	 * @param superind index of the element to be trasnfered
	 * @param super pointer to the destination mesh
	 */
	virtual int GetFromSuperMesh(int superind, TPZCompMesh *super);
	
	/**
	 * @brief Gives the commom father mesh of the specified mesh and the current submesh.
	 * @param mesh pointer to other mesh whosw want to know the commom father mesh
	 */
	virtual TPZCompMesh * CommonMesh (TPZCompMesh *mesh);
	
	/**
	 * @brief Returns the current submesh father mesh .
	 */
	virtual TPZCompMesh * FatherMesh() const;

	/** @} */
	
	/**
	 * @brief Optimizes the connections positions on block.
	 * void TPZSubCompMesh::PermuteExternalConnects(){
	 */
	void PermuteExternalConnects();
	
	/**
	 * @brief Computes the permutation vector which puts the internal connects to the first on the list
	 * Respect the previous order of the connects
	 */
	void ComputePermutationInternalFirst(TPZVec<int> &permute) const;
	
	/**
	 * @brief Permutes the potentially internal connects to the first on the list
	 * Respect the previous order of the connects
	 */
	void PermuteInternalFirst(TPZVec<int> &permute);
	
	/**
	 * @brief Prints the submesh information on the specified device/file out.
	 * 
	 * This method use the virtual method from Computacional Mesh class.
	 * @param out indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const;
	
	/**
	 * @brief Verifies if any element needs to be computed corresponding to the material ids
	 */
	bool NeedsComputing(const std::set<int> &matids);
	/**
	 * @brief Virtual Method to allocate new connect
	 */
	virtual int AllocateNewConnect(int nshape, int nstate, int order);
	
	/**
	 * @brief Gives the id node  of one local node in containing mesh.
	 */
	int NodeIndex(int nolocal, TPZCompMesh *super);
	
	
	/**
	 * @brief Virtual Method! See declaration in TPZCompEl class.
	 * 
	 * The use of this method in submesh class return -1 == Error!
	 */
	virtual int Dimension() const;
	
	/**
	 * @brief Computes the number of elements connected to each connect object
	 */
	virtual void ComputeNodElCon();
	
	/**
	 * @brief Computes the number of elements connected to each connect object
	 */
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
	virtual void SetConnectIndex(int inode, int index);
	
  	/**
	 * @brief Calculates the submesh stiffness matrix
	 */
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
		
	/**
	 * @brief Creates corresponding graphical element(s) if the dimension matches
	 * graphical elements are used to generate output files
	 * @param graphmesh graphical mesh where the element will be created
	 * @param dimension target dimension of the graphical element
	 */
	virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension);

	/**
	 * @brief Returns the connection index i.
	 */
	virtual int ConnectIndex(int i) const;
	
	/**
	 * @brief Returns the number of connections.
	 */
	virtual int NConnects() const;
	
    /**
	 * @brief Load the father mesh solution to all submesh connects -
  	 * (internal and external).
	 */
	virtual void LoadSolution();
		
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
	 * @param [in] axes indicating the direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
								 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
	
	/** @} */
	
	/**
	 * @brief Returns the unique identifier for reading/writing objects to streams
	 */
	virtual int ClassId() const;
	/**
	 @brief Saves the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/**
	 @brief Reads the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context);
	
	bool VerifyDatastructureConsistency();
    
    /// Set the number of rigid body modes associated with the internal degrees of freedom
    void SetNumberRigidBodyModes(int nrigid);
    
    /// Return the number of rigid body modes associated with the internal degrees of freedom
    int NumberRigidBodyModes();
	
};

#endif

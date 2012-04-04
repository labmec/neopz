/**
 * @file
 * @brief Contains the TPZDohrStructMatrix class which implements structural matrix divided in sub structures.
 */

#ifndef PZDOHRSTRUCTMATRIX 
#define PZDOHRSTRUCTMATRIX


#include "pzstrmatrix.h"
#include "tpzdohrassembly.h"
#include "pzsubcmesh.h"

/**
 * @ingroup substructure structural
 * @brief Implements structural matrix divided in sub structures. \ref structural "Structural Matrix" \ref substructure "Sub structure"
 * @author Philippe Devloo
 * @since 28/06/2010
 */
class TPZDohrStructMatrix : public TPZStructMatrix
{
	
	int fNumThreadsDecompose;
	
public:
	/** @brief We assume that the mesh consists of subcompmeshes */
	TPZDohrStructMatrix(TPZAutoPointer<TPZCompMesh> compmesh, int numthreads_compute, int numthreads_decompose);
	
	/** @brief Copy constructors */
	TPZDohrStructMatrix(const TPZDohrStructMatrix &copy);
	
	virtual ~TPZDohrStructMatrix();
	
	/** @brief Partition the mesh in submeshes */
	void SubStructure(int nsub);
	
	/** @brief This will create a DohrMatrix */
	virtual TPZMatrix<REAL> * Create();
	
	/**
	 * @brief This will return the pointer to the preconditioner AND abandon the pointer
	 * @warning This method can only be called once
	 */
	TPZAutoPointer<TPZMatrix<REAL> > Preconditioner()
	{
		TPZAutoPointer<TPZMatrix<REAL> > result = fDohrPrecond;
		fDohrPrecond = 0;
		return result;
	}
	
	/** @brief This will create a DohrMatrix and compute its matrices */
	virtual TPZMatrix<REAL> * CreateAssemble(TPZFMatrix<REAL> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/**
	 * @brief Assemble the global right hand side
	 */
	virtual void Assemble(TPZFMatrix<REAL> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/** @brief Creates a copy of itself */
	virtual TPZStructMatrix * Clone()
	{
		return new TPZDohrStructMatrix(*this);
	}
	
	/** 
	 * @brief Verifies if the subdomains are connected by sides of connectdimension and separate them if not
	 * @param domain_index
	 * @param nsub number of subdomains
	 * @param connectdimension
	 * @return returns the new number of subdomains
	 */
	int SeparateUnconnected(TPZVec<int> &domain_index, int nsub, int connectdimension);
	
	/**
	 * @brief Eliminates subdomains who are embedded in other subdomains
	 * @return returns the number of subdomains
	 */
	int ClusterIslands(TPZVec<int> &domain_index,int nsub,int connectdimension);
	
	
protected:
	
	TPZAutoPointer<TPZDohrAssembly> fDohrAssembly;
	
	TPZAutoPointer<TPZMatrix<REAL> > fDohrPrecond;
	
	/* @brief Get the global equation numbers of a substructure (and their inverse) */
	void IdentifyEqNumbers(TPZSubCompMesh *sub, std::map<int,int> &global, std::map<int,int> &globinv);
	
	/** @brief Computes the permutation vectors from the subcompmesh ordening to the "internal first" ordering */ 
	/** The mesh is modified during this method but is returned to its original state at the end of execution */
	void ComputeInternalEquationPermutation(TPZSubCompMesh *sub,
											TPZVec<int> &scatterpermute, TPZVec<int> &gatherpermute);
	
	/** @brief Identify the corner equations associated with a substructure */
	void IdentifySubCornerEqs(std::map<int,int> &globaltolocal, TPZVec<int> &cornereqs,
							  TPZVec<int> &coarseindex);
	
public:
	/** @brief Identify the external connects */
	void IdentifyExternalConnectIndexes();
	
private:
	/** @brief Identify cornernodes */
	void IdentifyCornerNodes();
	
	/** @brief The connect indexes which are external */
	TPZManVector<int> fExternalConnectIndexes;
	
	/** @brief A self administred pointer to the computational mesh */
	TPZAutoPointer<TPZCompMesh> fMesh;
	
	/** @brief The global equations defining the coarse matrix */
	std::set<int> fCornerEqs;
	
	/** @brief Mutexes (to choose which submesh is next) */
	pthread_mutex_t fAccessElement;
	
	friend struct ThreadDohrmanAssembly;
	
};

#include "tpzdohrsubstructCondense.h"

/**
 * @ingroup substructure
 * @brief Implements assembling by Dohrman algorithm.
 */
struct ThreadDohrmanAssembly {
	
	enum MTask {ENone, EComputeMatrix, EDecomposeInternal, EDecomposeBig};
	
	MTask fTask;
	TPZAutoPointer<TPZCompMesh> fMesh;
	int fSubMeshIndex;
	TPZAutoPointer<TPZDohrSubstructCondense> fSubstruct;
	TPZAutoPointer<TPZDohrAssembly> fAssembly;
	
	ThreadDohrmanAssembly(TPZAutoPointer<TPZCompMesh> mesh, int submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct,
						  TPZAutoPointer<TPZDohrAssembly> assembly) : 
	fTask(ENone), fMesh(mesh), fSubMeshIndex(submesh), fSubstruct(substruct), fAssembly(assembly)
	{
		
	}
	
	ThreadDohrmanAssembly(const ThreadDohrmanAssembly &copy) : fTask(copy.fTask), fMesh(copy.fMesh), fSubMeshIndex(copy.fSubMeshIndex),
	fSubstruct(copy.fSubstruct),fAssembly(copy.fAssembly)
	{
	}
	
	ThreadDohrmanAssembly(TPZAutoPointer<ThreadDohrmanAssembly> copy) : fTask(copy->fTask), fMesh(copy->fMesh), fSubMeshIndex(copy->fSubMeshIndex),
	fSubstruct(copy->fSubstruct),fAssembly(copy->fAssembly)
	{
	}
	
	ThreadDohrmanAssembly &operator=(const ThreadDohrmanAssembly &copy)
	{
		fTask = copy.fTask;
		fMesh = copy.fMesh;
		fSubMeshIndex = copy.fSubMeshIndex;
		fSubstruct = copy.fSubstruct;
		fAssembly = copy.fAssembly;
		return *this;
	}
	
	void AssembleMatrices(pthread_mutex_t &testthread);
};

/**
 * @ingroup substructure
 * @brief Implements a list of Dohrman assembling and control thread and semaphores.
 */
struct ThreadDohrmanAssemblyList {
	
	ThreadDohrmanAssemblyList();
	
	~ThreadDohrmanAssemblyList();
	
	std::list<TPZAutoPointer<ThreadDohrmanAssembly> > fList;
	
	void Append(TPZAutoPointer<ThreadDohrmanAssembly> object);
	
	/** @brief Returns an object and removes it from the list in a thread safe way */
	TPZAutoPointer<ThreadDohrmanAssembly> NextObject();
	
	static void *ThreadWork(void *voidptr);
	
	/** @brief Mutexes (to choose which submesh is next) */
	pthread_mutex_t fAccessElement;
	
	/** @brief mutex to debug the assembly process */
	pthread_mutex_t fTestThreads;
};

#endif

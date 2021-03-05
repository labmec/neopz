/**
 * @file
 * @brief Contains the TPZDohrStructMatrix class which implements structural matrix divided in sub structures.
 */

#ifndef PZDOHRSTRUCTMATRIX 
#define PZDOHRSTRUCTMATRIX

#include "pzstrmatrix.h"
#include "tpzdohrassembly.h"
#include "pzsubcmesh.h"

template<class TVar>
struct ThreadDohrmanAssembly;

/**
 * @ingroup substructure structural
 * @brief Implements structural matrix divided in sub structures. \ref structural "Structural Matrix" \ref substructure "Sub structure"
 * @author Philippe Devloo
 * @since 28/06/2010
 */
class TPZDohrStructMatrix : public TPZStructMatrix
{
		
public:
	
	/** @brief We assume that the mesh consists of subcompmeshes */
        TPZDohrStructMatrix(TPZAutoPointer<TPZCompMesh> compmesh);
	
	/** @brief Copy constructors */
	TPZDohrStructMatrix(const TPZDohrStructMatrix &copy);
	
	virtual ~TPZDohrStructMatrix();
	
	/** @brief Partition the mesh in submeshes */
	void SubStructure(int nsub);
	
	/** @brief This will create a DohrMatrix */
	virtual TPZMatrix<STATE> * Create() override;
	
	/**
	 * @brief This will return the pointer to the preconditioner AND abandon the pointer
	 * @warning This method can only be called once
	 */
	TPZAutoPointer<TPZMatrix<STATE> > Preconditioner()
	{
		TPZAutoPointer<TPZMatrix<STATE> > result = fDohrPrecond;
		//fDohrPrecond = 0; Essa linha me ferra pois nao posso pegar as subestruturas depois
		return result;
	}
	
	/** @brief This will create a DohrMatrix and compute its matrices */
	virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose) override;
	
    /**
	 * @brief Assemble the global system of equations into the matrix which has already been created
	 */
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                          unsigned numthreads_assemble, unsigned numthreads_decompose) override;

    void AssembleTBB(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

	/**
	 * @brief Assemble the global right hand side
	 */
	virtual void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
	
	/** @brief Creates a copy of itself */
	virtual TPZStructMatrix * Clone() override
	{
		return new TPZDohrStructMatrix(*this);
	}
    
    /** @brief Return the number of cornereqs */
    int NumberCornerEqs() const
    {
        return fCornerEqs.size();
    }
    
    // FOR DEBUG PURPOSES
    TPZAutoPointer<TPZDohrAssembly<STATE> > Assembly() {
        return fDohrAssembly;
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
    
    void Write(TPZStream &str, int withclassid) const override;
    
    void Read(TPZStream &str, void *context) override;
	
	
protected:
	
	TPZAutoPointer<TPZDohrAssembly<STATE> > fDohrAssembly;
	
	TPZAutoPointer<TPZMatrix<STATE> > fDohrPrecond;
	
	/* @brief Get the global equation numbers of a substructure (and their inverse) */
	void IdentifyEqNumbers(TPZSubCompMesh *sub, std::map<int,int> &global, std::map<int,int> &globinv);
	
	/** @brief Computes the permutation vectors from the subcompmesh ordening to the "internal first" ordering */ 
	/** The mesh is modified during this method but is returned to its original state at the end of execution */
	void ComputeInternalEquationPermutation(TPZSubCompMesh *sub,
											TPZVec<int> &scatterpermute, TPZVec<int> &gatherpermute);
	
	/** @brief Identify the corner equations associated with a substructure */
	void IdentifySubCornerEqs(std::map<int,int> &globaltolocal, TPZVec<int> &cornereqs,
							  TPZVec<int> &coarseindex);
    
    /** @brief Set the domain index of the lower dimension elements equal to the domain index of their neighbour */
    void CorrectNeighbourDomainIndex(TPZCompMesh *cmesh, TPZVec<int> &domainindex);
	
public:
	/** @brief Identify the external connects */
	void IdentifyExternalConnectIndexes();
	
private:
	TPZDohrStructMatrix();
        
	/** @brief Identify cornernodes */
	void IdentifyCornerNodes();
	
	/** @brief The connect indexes which are external */
	TPZManVector<int> fExternalConnectIndexes;
	
	/** @brief The global equations defining the coarse matrix */
	std::set<int> fCornerEqs;
	
	/** @brief Mutexes (to choose which submesh is next) */
	std::mutex fAccessElement;
	
	friend struct ThreadDohrmanAssembly<STATE>;
        
        friend TPZPersistenceManager;
	
};

#include "tpzdohrsubstructCondense.h"

/**
 * @ingroup substructure
 * @brief Implements assembling by Dohrman algorithm.
 */
template<class TVar>
struct ThreadDohrmanAssembly {
	
	enum MTask {ENone, EComputeMatrix, EDecomposeInternal, EDecomposeBig};
	
	MTask fTask;
	TPZAutoPointer<TPZCompMesh> fMesh;
	int fSubMeshIndex;
	TPZAutoPointer<TPZDohrSubstructCondense<TVar> > fSubstruct;
	TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembly;
	
	ThreadDohrmanAssembly(TPZAutoPointer<TPZCompMesh> mesh, int submesh, TPZAutoPointer<TPZDohrSubstructCondense<TVar> > substruct,
						  TPZAutoPointer<TPZDohrAssembly<TVar> > assembly) : 
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
	
	void AssembleMatrices(std::mutex &testthread, int numa_node);
};

/**
 * @ingroup substructure
 * @brief Implements a list of Dohrman assembling and control thread and semaphores.
 */
template<class TVar>
struct ThreadDohrmanAssemblyList {
	
	ThreadDohrmanAssemblyList();
	ThreadDohrmanAssemblyList(ThreadDohrmanAssemblyList<TVar> &cpy);
	
	~ThreadDohrmanAssemblyList();
	
	std::list<TPZAutoPointer<ThreadDohrmanAssembly<TVar> > > fList;
	
	void Append(TPZAutoPointer<ThreadDohrmanAssembly<TVar> > object);
	
	/** @brief Returns an object and removes it from the list in a thread safe way */
	TPZAutoPointer<ThreadDohrmanAssembly<TVar> > NextObject();
	
	static void *ThreadWork(void *voidptr);
	
	/** @brief Mutexes (to choose which submesh is next) */
	std::mutex fAccessElement;
	
	/** @brief mutex to debug the assembly process */
	std::mutex fTestThreads;
};

#endif

/*
 *  pzdohrstructmatrix.h
 *  SubStruct
 *
 *  Created by Philippe Devloo on 28/06/10.
 *  Copyright 2010 UNICAMP. All rights reserved.
 *
 */
#ifndef PZDOHRSTRUCTMATRIX 
#define PZDOHRSTRUCTMATRIX


#include "pzstrmatrix.h"
#include "tpzdohrassembly.h"
#include "pzsubcmesh.h"

class TPZDohrStructMatrix : public TPZStructMatrix
{
	
	int fNumThreadsDecompose;
	
public:
	// we assume that the mesh consists of subcompmeshes
	TPZDohrStructMatrix(TPZAutoPointer<TPZCompMesh> compmesh, int numthreads_assemble, int numthreads_decompose);
	
	// copy constructors
	TPZDohrStructMatrix(const TPZDohrStructMatrix &copy);
	
	virtual ~TPZDohrStructMatrix();
	
	// partition the mesh in submeshes
	void SubStructure(int nsub);
	
	// this will create a DohrMatrix
	virtual TPZMatrix * Create();
	
	/// this will return the pointer to the preconditioner AND abandon the pointer
	// @warning This method can only be called once
	TPZAutoPointer<TPZMatrix> Preconditioner()
	{
		TPZAutoPointer<TPZMatrix> result = fDohrPrecond;
		fDohrPrecond = 0;
		return result;
	}
	
	// this will create a DohrMatrix and compute its matrices
	virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	// create a copy of itself
	virtual TPZStructMatrix * Clone()
	{
		return new TPZDohrStructMatrix(*this);
	}
	
	/// Verifies if the subdomains are connected by sides of connectdimension and separate them if not
	// nsub : number of subdomains
	// returns the new number of subdomains
	int SeparateUnconnected(TPZVec<int> &domain_index, int nsub, int connectdimension);
	
	/// Eliminate subdomains who are embedded in other subdomains
	// returns the number of subdomains
	int ClusterIslands(TPZVec<int> &domain_index,int nsub,int connectdimension);

	
protected:
	
	TPZAutoPointer<TPZDohrAssembly> fDohrAssembly;
	
	TPZAutoPointer<TPZMatrix > fDohrPrecond;
	
	/// get the global equation numbers of a substructure (and their inverse)
	void IdentifyEqNumbers(TPZSubCompMesh *sub, std::map<int,int> &global, std::map<int,int> &globinv);
	
	// computes the permutation vectors from the subcompmesh ordening to the "internal first" ordering
	// the mesh is modified during this method but is returned to its original state at the end of execution
	void ComputeInternalEquationPermutation(TPZSubCompMesh *sub,
											TPZVec<int> &scatterpermute, TPZVec<int> &gatherpermute);
	
	/// Identify the corner equations associated with a substructure
	void IdentifySubCornerEqs(std::map<int,int> &globaltolocal, TPZVec<int> &cornereqs,
							  TPZVec<int> &coarseindex);
	
public:
	/// Identify the external connects
	void IdentifyExternalConnectIndexes();

private:
	/// identify cornernodes
	void IdentifyCornerNodes();
	
	/// The connect indexes which are external
	TPZManVector<int> fExternalConnectIndexes;
	
	/// A self administred pointer to the computational mesh
	TPZAutoPointer<TPZCompMesh> fMesh;

	// the global equations defining the coarse matrix
	std::set<int> fCornerEqs;
	
	// mutexes (to choose which submesh is next)
	pthread_mutex_t fAccessElement;
	
	friend class ThreadDohrmanAssembly;

};

#include "tpzdohrsubstructCondense.h"

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

struct ThreadDohrmanAssemblyList {

	ThreadDohrmanAssemblyList();
	
	~ThreadDohrmanAssemblyList();
	
	std::list<TPZAutoPointer<ThreadDohrmanAssembly> > fList;
	
	void Append(TPZAutoPointer<ThreadDohrmanAssembly> object);
	
	// returns an object and removes it from the list in a thread safe way
	TPZAutoPointer<ThreadDohrmanAssembly> NextObject();
	
	static void *ThreadWork(void *voidptr);
	
	// mutexes (to choose which submesh is next)
	pthread_mutex_t fAccessElement;
	
	/// mutex to debug the assembly process
	pthread_mutex_t fTestThreads;
};

#endif
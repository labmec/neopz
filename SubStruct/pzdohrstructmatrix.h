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
	
public:
	// we assume that the mesh consists of subcompmeshes
	TPZDohrStructMatrix(TPZAutoPointer<TPZCompMesh> compmesh);
	
	// copy constructors
	TPZDohrStructMatrix(const TPZStructMatrix &copy);
	
	virtual ~TPZDohrStructMatrix();
	
	// partition the mesh in submeshes
	static void SubStructure(TPZAutoPointer<TPZCompMesh> cmesh, int nsub);
	
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

	/// identify cornernodes
	void IdentifyCornerNodes();
	
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
	
	TPZAutoPointer<TPZCompMesh> fMesh;
	int fSubMeshIndex;
	TPZAutoPointer<TPZDohrSubstructCondense> fSubstruct;
	TPZAutoPointer<TPZDohrAssembly> fAssembly;
	
	ThreadDohrmanAssembly(TPZAutoPointer<TPZCompMesh> mesh, int submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct,
						  TPZAutoPointer<TPZDohrAssembly> assembly) : 
		fMesh(mesh), fSubMeshIndex(submesh), fSubstruct(substruct), fAssembly(assembly)
	{
		
	}
	
	void AssembleMatrices(pthread_mutex_t &testthread);
};

struct ThreadDohrmanAssemblyList {

	ThreadDohrmanAssemblyList();
	
	~ThreadDohrmanAssemblyList();
	
	std::list<TPZAutoPointer<ThreadDohrmanAssembly> > fList;
	
	void Append(TPZAutoPointer<ThreadDohrmanAssembly> object);
	
	TPZAutoPointer<ThreadDohrmanAssembly> NextObject();
	
	static void *ThreadWork(void *voidptr);
	
	// mutexes (to choose which submesh is next)
	pthread_mutex_t fAccessElement;
	
	/// mutex to debug the assembly process
	pthread_mutex_t fTestThreads;
};

#endif
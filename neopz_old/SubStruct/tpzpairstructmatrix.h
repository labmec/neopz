/*
 *  tpzpairstructmatrix.h
 *  SubStruct
 *
 *  Created by Philippe Devloo on 20/04/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */
#ifndef PAIRSTRUCTMATRIX
#define PAIRSTRUCTMATRIX

#include "pzcmesh.h"
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzmatrix.h"
#include "TPZGuiInterface.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"

class TPZPairStructMatrix
{
	TPZCompMesh *fMesh;
	TPZVec<int> fPermuteScatter;
	std::set<int> fMaterialIds;
	int fNumThreads;
	TPZVec<int> fElementOrder;
	
	void OrderElement();
	
	void PermuteScatter(TPZVec<int> &index);
	
public:
	
	static int gNumThreads;
	
	TPZPairStructMatrix(TPZCompMesh *mesh, TPZVec<int> &permutescatter);

	void SetNumThreads(int numthreads)
	{
		fNumThreads = numthreads;
	}
	
	/// Set the set of material ids which will be considered when assembling the system
	void SetMaterialIds(const std::set<int> &materialids);

	void Assemble(int mineq, int maxeq, TPZMatrix *first, TPZMatrix *second, TPZFMatrix &rhs);

	void SerialAssemble(int mineq, int maxeq, TPZMatrix *first, TPZMatrix *second, TPZFMatrix &rhs);
	
	void MultiThread_Assemble(int mineq, int maxeq, TPZMatrix *first, TPZMatrix *second, TPZFMatrix &rhs);
	
	struct ThreadData
	{
		// Initialize the mutex semaphores and others
		ThreadData(TPZCompMesh &mesh,TPZMatrix &mat1, TPZMatrix &mat2, TPZFMatrix &rhs, int mineq, int maxeq);
		// Destroy the mutex semaphores and others
		~ThreadData();
		// current structmatrix object
		TPZCompMesh *fMesh;
		// mutexes (to choose which element is next)
		pthread_mutex_t fAccessElement;
		// semaphore (to wake up the first assembly thread)
		TPZSemaphore fAssembly1;
		// semaphore (to wake up the second assembly thread)
		TPZSemaphore fAssembly2;
		// global matrix1
		TPZMatrix *fGlobMatrix1;
		// global matrix2
		TPZMatrix *fGlobMatrix2;
		// global rhs
		TPZFMatrix *fGlobRhs;
		// minimum equation to be assembled
		int fMinEq;
		// maximum equation to be assembled
		int fMaxEq;
		/// material identifiers which need to be computed
		std::set<int> fMaterialIds;
		/// vector which defines the permutation of all equations to internal equations
		TPZVec<int> fPermuteScatter;
		// list of computed element matrices (autopointers?)
		std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted1;
		// list of computed element matrices (autopointers?)
		std::map<int, TPZAutoPointer<TPZElementMatrix> > fSubmitted2;
		// elements which are being processed maintained by the first global matrix
		std::set<int> fProcessed1;
		// elements which are being processed maintained by the second global matrix
		std::set<int> fProcessed2;
		// current element
		int fNextElement;
		// look for an element index which needs to be computed and put it on the stack
		int NextElement();
		// put the computed element matrices in the map
		void ComputedElementMatrix(int iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef);
		// The function which will compute the matrices
		static void *ThreadWork(void *threaddata);
		// The function which will compute the assembly
		static void *ThreadAssembly1(void *threaddata);
		// The function which will compute the assembly
		static void *ThreadAssembly2(void *threaddata);
		
		/// Establish whether the element should be computed
		bool ShouldCompute(int matid)
		{
			return fMaterialIds.size()==0 || fMaterialIds.find(matid) != fMaterialIds.end();
		}
		void PermuteScatter(TPZVec<int> &index);
		
	};
	
};

#endif

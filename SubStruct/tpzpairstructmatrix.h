/**
 * @file
 * @brief Contains the TPZPairStructMatrix class.
 */
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

/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub Structure"
 */
class TPZPairStructMatrix
{
	TPZCompMesh *fMesh;
	TPZVec<int> fPermuteScatter;
	std::set<int> fMaterialIds;
	int fNumThreads;
	
	void PermuteScatter(TPZVec<int> &index);
	
public:
	
	static int gNumThreads;
	
	TPZPairStructMatrix(TPZCompMesh *mesh, TPZVec<int> &permutescatter);
	
	void SetNumThreads(int numthreads)
	{
		fNumThreads = numthreads;
	}
	
	/** @brief Set the set of material ids which will be considered when assembling the system */
	void SetMaterialIds(const std::set<int> &materialids);
	
	void Assemble(int mineq, int maxeq, TPZMatrix<REAL> *first, TPZMatrix<REAL> *second, TPZFMatrix<REAL> &rhs);
	
	void SerialAssemble(int mineq, int maxeq, TPZMatrix<REAL> *first, TPZMatrix<REAL> *second, TPZFMatrix<REAL> &rhs);
	
	void MultiThread_Assemble(int mineq, int maxeq, TPZMatrix<REAL> *first, TPZMatrix<REAL> *second, TPZFMatrix<REAL> &rhs);
	
	/** @brief Contains the thread data for matrices divided in sub structures. */
	struct ThreadData
	{
		/** @brief Initialize the mutex semaphores and others */
		ThreadData(TPZCompMesh &mesh,TPZMatrix<REAL> &mat1, TPZMatrix<REAL> &mat2, TPZFMatrix<REAL> &rhs, int mineq, int maxeq);
		/** @brief Destroy the mutex semaphores and others */
		~ThreadData();
		/** @brief Current structmatrix object */
		TPZCompMesh *fMesh;
		/** @brief Mutexes (to choose which element is next) */
		pthread_mutex_t fAccessElement;
		/** @brief Semaphore (to wake up the first assembly thread) */
		TPZSemaphore fAssembly1;
		/** @brief Semaphore (to wake up the second assembly thread) */
		TPZSemaphore fAssembly2;
		/** @brief Global matrix1 */
		TPZMatrix<REAL> *fGlobMatrix1;
		/** @brief Global matrix2 */
		TPZMatrix<REAL> *fGlobMatrix2;
		/** @brief Global rhs */
		TPZFMatrix<REAL> *fGlobRhs;
		/** @brief Minimum equation to be assembled */
		int fMinEq;
		/** @brief Maximum equation to be assembled */
		int fMaxEq;
		/** @brief Material identifiers which need to be computed */
		std::set<int> fMaterialIds;
		/** @brief Vector which defines the permutation of all equations to internal equations */
		TPZVec<int> fPermuteScatter;
		/** @brief List of computed element matrices (autopointers?) */
		std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted1;
		/** @brief List of computed element matrices (autopointers?) */
		std::map<int, TPZAutoPointer<TPZElementMatrix> > fSubmitted2;
		/** @brief Elements which are being processed maintained by the first global matrix */
		std::set<int> fProcessed1;
		/** @brief Elements which are being processed maintained by the second global matrix */
		std::set<int> fProcessed2;
		/** @brief Current element */
		int fNextElement;
		/** @brief Look for an element index which needs to be computed and put it on the stack */
		int NextElement();
		/** @brief Put the computed element matrices in the map */
		void ComputedElementMatrix(int iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef);
		/** @brief The function which will compute the matrices */
		static void *ThreadWork(void *threaddata);
		/** @brief The function which will compute the assembly */
		static void *ThreadAssembly1(void *threaddata);
		/** @brief The function which will compute the assembly */
		static void *ThreadAssembly2(void *threaddata);
		
		/** @brief Establish whether the element should be computed */
		bool ShouldCompute(int matid)
		{
			return fMaterialIds.size()==0 || fMaterialIds.find(matid) != fMaterialIds.end();
		}
		void PermuteScatter(TPZVec<int> &index);
		
	};
	
};

#endif

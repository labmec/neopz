/**
 * @file
 * @brief Contains the TPZPairStructMatrix class.
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
#include "TPZStructMatrix.h"
//class TPZStructMatrix;


//TODO: NEEDS REVIEW
/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub Structure"
 */
class TPZPairStructMatrix
{
	TPZVec<int> fPermuteScatter;
	
    TPZAutoPointer<TPZStructMatrix>fStrMatrix;
	
	void PermuteScatter(TPZVec<int> &index);
	void PermuteScatter(TPZVec<int64_t> &index);
	
public:
	
	static int gNumThreads;
	
	TPZPairStructMatrix(TPZCompMesh *mesh, TPZVec<int> &permutescatter);
    
    ~TPZPairStructMatrix()
    {
        
    }
	
	void SetNumThreads(int numthreads)
	{
		fStrMatrix->SetNumThreads(numthreads);
	}
	
	/** @brief Set the set of material ids which will be considered when assembling the system */
	void SetMaterialIds(const std::set<int> &materialids);
	
	void Assemble(TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs);
	void TBBAssemble(TPZMatrix<STATE> *first, 
                     TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs);

	void SerialAssemble(TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs);
	
	void MultiThread_Assemble(TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs);
	
	/** @brief Contains the thread data for matrices divided in sub structures. */
	struct ThreadData
	{
		/** @brief Initialize the mutex semaphores and others */
		ThreadData(TPZStructMatrix *strmatrix,TPZMatrix<STATE> &mat1, TPZMatrix<STATE> &mat2, TPZFMatrix<STATE> &rhs);
		/** @brief Destroy the mutex semaphores and others */
		~ThreadData();
		/** @brief Current structmatrix object */
		TPZStructMatrix *fStrMatrix;
		/** @brief Mutexes (to choose which element is next) */
		std::mutex fAccessElement;
		/** @brief Semaphore (to wake up the first assembly thread) */
		TPZSemaphore fAssembly1;
		/** @brief Semaphore (to wake up the second assembly thread) */
		TPZSemaphore fAssembly2;
		/** @brief Global matrix1 */
		TPZMatrix<STATE> *fGlobMatrix1;
		/** @brief Global matrix2 */
		TPZMatrix<STATE> *fGlobMatrix2;
		/** @brief Global rhs */
		TPZFMatrix<STATE> *fGlobRhs;
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
		void ComputedElementMatrix(int iel, TPZAutoPointer<TPZElementMatrixT<STATE>> &ek, TPZAutoPointer<TPZElementMatrixT<STATE>> &ef);
		/** @brief The function which will compute the matrices */
		static void *ThreadWork(void *threaddata);
		/** @brief The function which will compute the assembly */
		static void *ThreadAssembly1(void *threaddata);
		/** @brief The function which will compute the assembly */
		static void *ThreadAssembly2(void *threaddata);
		
		/** @brief Establish whether the element should be computed */
		bool ShouldCompute(int matid)
		{
            return fStrMatrix->ShouldCompute(matid);
		}
		void PermuteScatter(TPZVec<int> &index);
		void PermuteScatter(TPZVec<int64_t> &index);
		
	};
	
};

#endif

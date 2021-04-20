/**
 * @file
 * @brief Contains the TPZStructMatrixOT class which is a TPZStructMatrixBase using graph coloring to define the order
to process the elements and each color is processed and 
synchronized.
 */

#ifndef TPZStructMatrixOT_H
#define TPZStructMatrixOT_H

#include <set>
#include <map>
#include <semaphore.h>
#include <mutex>
#include <condition_variable>
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"
#include "TPZEquationFilter.h"
#include "TPZGuiInterface.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

class TPZStructMatrixOT;
#include "TPZStructMatrixBase.h"

/**
 * @brief Refines geometrical mesh (all the elements) num times
 * @ingroup geometry
 */
//void UniformRefine(int num, TPZGeoMesh &m);

/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrixOT : public TPZStructMatrixBase{
    
public:
    
    TPZStructMatrixOT(): TPZStructMatrixBase(){}
    
    TPZStructMatrixOT(TPZCompMesh *);
    
    TPZStructMatrixOT(TPZAutoPointer<TPZCompMesh> cmesh);
    
    TPZStructMatrixOT(const TPZStructMatrixOT &copy);
    virtual ~TPZStructMatrixOT() = default;
        
    TPZBaseMatrix * Create() override;
    
    TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose) {
        SetNumThreads(numthreads_assemble);
        return CreateAssemble(rhs, guiInterface);
    }
    
    TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override final;
    
    TPZStructMatrixOT * Clone() override final;
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    void Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override final;
    
    /** @brief Assemble the global right hand side */
    void Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override final;

    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;
    
protected:
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Serial_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void Serial_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
public:
    
    /** @brief Find the order to assemble the elements */
    static void OrderElement(TPZCompMesh *cmesh, TPZVec<int64_t> &ElementOrder);
    
    /** @brief Create blocks of elements to parallel processing */
    static void ElementColoring(TPZCompMesh *cmesh, TPZVec<int64_t> &elSequence, TPZVec<int64_t> &elSequenceColor, TPZVec<int64_t> &elBlocked, TPZVec<int64_t> &NumelColors);    
    
protected:
    
    /** @brief Structure to manipulate thread to solve system equations */
    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrixOT *strmat,int seqnum, TPZBaseMatrix &mat, TPZBaseMatrix &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrixOT *strmat, int seqnum, TPZBaseMatrix &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief The function which will compute the assembly */
        bool ShouldCompute(int matid)
        {
            return fStruct->ShouldCompute(matid);
        }
        
        /** @brief Current structmatrix object */
        TPZStructMatrixOT *fStruct;
        /** @brief Gui interface object */
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        /** @brief Global matrix */
        TPZBaseMatrix *fGlobMatrix;
        /** @brief Global rhs vector */
        TPZBaseMatrix *fGlobRhs;
        
        std::atomic<int64_t> *fCurrentIndex;
        /** @brief sequence number of the thread */
        int fThreadSeqNum;

        /** @brief vector indicating whether an element has been computed */
        TPZVec<int64_t> *fComputedElements;
        /** @brief Mutexes (to choose which element is next) */
        std::mutex fMutexAccessElement;
        std::condition_variable fConditionVar;
        
        int *fSomeoneIsSleeping;
        
        /// Vector for mesh coloring
        TPZVec<int64_t> *fElBlocked, *fElSequenceColor;
        
        /// All elements below or equal this index have been computed
        int64_t *fElementCompleted;
        
        static void *ThreadWorkResidual(void *datavoid);
    };
    friend struct ThreadData;
protected:
    
    /** @brief Vectors for mesh coloring */
    TPZVec<int64_t> fElBlocked, fElSequenceColor;
    
    /// vector of the size of the elements containing 0 or 1 if the element has been computed (in the order of computation sequence)
    TPZVec<int64_t> fElementsComputed;
    
    /// All elements below or equal this index have been computed
    int64_t fElementCompleted;
    
    /// variable indicating if a thread is sleeping
    int fSomeoneIsSleeping;
    
    std::atomic<int64_t> fCurrentIndex;
};

#endif

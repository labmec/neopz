/**
 * @file
 * @brief Contains the TPZStructMatrixOT class which is a TPZStructMatrixBase using graph coloring to define the order
to process the elements and each color is processed and 
synchronized.
 */

#ifndef TPZStructMatrixOT_H
#define TPZStructMatrixOT_H

#include "TPZStrMatParInterface.h"
#include "TPZSemaphore.h"
#include <mutex>

//forward declarations
class TPZElementMatrix;
class TPZBaseMatrix;
class TPZStructMatrix;


/**
 * @brief Parallel layer for struct matrices
 * using graph coloring to define the order 
 * in which to process the elements (Devloo-Gilvan).
 * Each color is processed and synchronized.
 * @ingroup structural
 */
template<class TVar>
class TPZStructMatrixOT : public virtual TPZStrMatParInterface{
public:
    //! Default constructor
    TPZStructMatrixOT();
    //! Copy constructor
    TPZStructMatrixOT(const TPZStructMatrixOT &copy);
    //! Move constructor
    TPZStructMatrixOT(TPZStructMatrixOT &&copy);
    //! Virtual destructor
    virtual ~TPZStructMatrixOT() = default;
    //! Copy assignment operator
    TPZStructMatrixOT& operator=(const TPZStructMatrixOT &);
    //! Move assignment operator
    TPZStructMatrixOT& operator=(TPZStructMatrixOT &&);

    /**
     * Functions overriden from TPZStrMatParInterface
     */
    //@{
    
    void Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    
    void Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    void InitCreateAssemble() override;
    //@}

    
    //@{
    //!Read and Write methods
    int ClassId() const override;
    
    void Read(TPZStream &buf, void *context) override;
    
    void Write(TPZStream &buf, int withclassid) const override;
    //@}
protected:
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Serial_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void Serial_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface); 

    /** @brief Find the order to assemble the elements */
    static void OrderElement(TPZCompMesh *cmesh, TPZVec<int64_t> &ElementOrder);
    
    /** @brief Create blocks of elements to parallel processing */
    static void ElementColoring(TPZCompMesh *cmesh, TPZVec<int64_t> &elSequence, TPZVec<int64_t> &elSequenceColor, TPZVec<int64_t> &elBlocked, TPZVec<int64_t> &NumelColors);   
    
    /** @brief Structure to manipulate thread to solve system equations */
    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrix *strmat,int seqnum, TPZBaseMatrix &mat, TPZBaseMatrix &rhs, const std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface, std::atomic<int64_t> *fCurrentIndex, bool computeRhs);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrix *strmat, int seqnum, TPZBaseMatrix &rhs, const std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface, std::atomic<int64_t> *fCurrentIndex);
        ~ThreadData();
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief The function which will compute the assembly */
        bool ShouldCompute(int matid) const;
        
        /** @brief Current structmatrix object */
        TPZStructMatrix *fStruct;
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
        /// Whether the rhs is being computed
        bool fComputeRhs;
    };
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


extern template class TPZStructMatrixOT<STATE>;
extern template class TPZStructMatrixOT<CSTATE>;
#endif

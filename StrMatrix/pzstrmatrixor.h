/**
 * @file
 * @brief Contains the TPZStructMatrixOR class which implements a
 * parallel strategy for a TPZStructMatrix using threads and mutexes
 * without mesh coloring.
 */

#ifndef TPZStructMatrixOR_H
#define TPZStructMatrixOR_H

#include "TPZStrMatParInterface.h"
#include "TPZSemaphore.h"
#include <mutex>

//forward declarations
class TPZElementMatrix;
template<class T>
class TPZElementMatrixT;
class TPZBaseMatrix;
class TPZStructMatrix;


/**
 * @brief Parallel layer for struct matrices using 
 * the producer-consumer paradigm.
 * This is the original and stable version of the multi threaded assemble.
 * @ingroup structural
 */
template<class TVar>
class TPZStructMatrixOR : public virtual TPZStrMatParInterface {    
public:
    //! Default constructor.
    TPZStructMatrixOR();
    //! Copy constructor
    TPZStructMatrixOR(const TPZStructMatrixOR &copy) = default;
    //! Move constructor
    TPZStructMatrixOR(TPZStructMatrixOR &&copy) = default;
    //! Virtual destructor
    virtual ~TPZStructMatrixOR() = default;
    //! Copy assignment operator
    TPZStructMatrixOR& operator=(const TPZStructMatrixOR &) = default;
    //! Move assignment operator
    TPZStructMatrixOR& operator=(TPZStructMatrixOR &&) = default;

    /**
     * Functions overriden from TPZStrMatParInterface
     */
    //@{
    
    void Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    
    void Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
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

    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrix *strmat,TPZBaseMatrix &mat, TPZBaseMatrix &rhs, const std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface, bool computeRhs);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrix *strmat, TPZBaseMatrix &rhs, const std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface,bool computeRhs);
        /** @brief Destructor: Destroy the mutex semaphores and others */
        ~ThreadData();
        /** @brief Look for an element index which needs to be computed and put it on the stack */
        int64_t NextElement();
        /** @brief Put the computed element matrices in the map */
        void ComputedElementMatrix(int64_t iel, TPZAutoPointer<TPZElementMatrixT<TVar>> &ek, TPZAutoPointer<TPZElementMatrixT<TVar>> &ef);
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief The function which will compute the assembly */
        static void *ThreadAssembly(void *threaddata);
        /** @brief Establish whether the element should be computed */
        bool ShouldCompute(int matid) const;
        /** @brief Current structmatrix object */
        TPZStructMatrix *fStruct;
        /** @brief Gui interface object */
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        /** @brief Global matrix */
        TPZBaseMatrix *fGlobMatrix;
        /** @brief Global rhs vector */
        TPZBaseMatrix *fGlobRhs;
        /** @brief List of computed element matrices (autopointers?) */
        std::map<int, std::pair< TPZAutoPointer<TPZElementMatrixT<TVar>>, TPZAutoPointer<TPZElementMatrixT<TVar>> > > fSubmitted;
        /** @brief Elements which are being processed */
        std::set<int> fProcessed;
        /** @brief  Current element */
        int64_t fNextElement;
        /** @brief Mutexes (to choose which element is next) */
        std::mutex fMutexAccessElement;
        /** @brief Semaphore (to wake up assembly thread) */
//        std::condition_variable fAssembly;
        TPZSemaphore fAssembly;
        bool fComputeRhs{true};
    };
};


extern template class TPZStructMatrixOR<STATE>;
extern template class TPZStructMatrixOR<CSTATE>;
#endif

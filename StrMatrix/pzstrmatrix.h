/**
 * @file
 * @brief Contains the TPZStructMatrixOR class which is a TPZStructMatrix using threads and mutexes without mesh coloring for parallelization.
 */

#ifndef TPZStructMatrixOR_H
#define TPZStructMatrixOR_H


#include <mutex>

class TPZStructMatrixOR;
#include "TPZStructMatrixBase.h"
#include "TPZSemaphore.h"

//forward declarations
class TPZElementMatrix;
class TPZBaseMatrix;
/**
 * @brief Refines geometrical mesh (all the elements) num times
 * @ingroup geometry
 */
//void UniformRefine(int num, TPZGeoMesh &m);

/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrixOR : public TPZStructMatrixBase {    
public:
    //! Default constructor.
    TPZStructMatrixOR() = default;
    //! Constructor setting a computational mesh via raw pointer.
    TPZStructMatrixOR(TPZCompMesh *);
    //! Constructor setting a computational mesh via TPZAutoPointer.
    TPZStructMatrixOR(TPZAutoPointer<TPZCompMesh> cmesh);
    //! Copy constructor
    TPZStructMatrixOR(const TPZStructMatrixOR &copy) = default;
    //! Move constructor (deleted)
    TPZStructMatrixOR(const TPZStructMatrixOR &&copy) = delete;
    //! Virtual destructor
    virtual ~TPZStructMatrixOR() = default;

    TPZStructMatrixOR& operator=(const TPZStructMatrixOR &) = default;

    TPZStructMatrixOR& operator=(const TPZStructMatrixOR &&) = delete;

    /**
     * Functions overriden from TPZStructMatrixBase
     */
    //@{
    
    virtual TPZMatrix<STATE> * Create() override;
    
    virtual TPZStructMatrixOR * Clone() override;

    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    
    virtual void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    //@}
    
    using TPZStructMatrixBase::CreateAssemble;
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose);
    
    //virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                          unsigned numthreads_assemble, unsigned numthreads_decompose) {
        std::cout << "Nothing to do." << std::endl;
    }

    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;


protected:
    //forward declaration
    struct ThreadData;
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Serial_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    
    friend struct ThreadData;

    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrixOR *strmat,TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrixOR *strmat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        /** @brief Destructor: Destroy the mutex semaphores and others */
        ~ThreadData();
        /** @brief Look for an element index which needs to be computed and put it on the stack */
        int64_t NextElement();
        /** @brief Put the computed element matrices in the map */
        void ComputedElementMatrix(int64_t iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef);
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief The function which will compute the assembly */
        static void *ThreadAssembly(void *threaddata);
        /** @brief Establish whether the element should be computed */
        bool ShouldCompute(int matid)
        {
            return fStruct->ShouldCompute(matid);
        }
        
        /** @brief Current structmatrix object */
        TPZStructMatrixOR *fStruct;
        /** @brief Gui interface object */
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        /** @brief Global matrix */
        TPZBaseMatrix *fGlobMatrix;
        /** @brief Global rhs vector */
        TPZBaseMatrix *fGlobRhs;
        /** @brief List of computed element matrices (autopointers?) */
        std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted;
        /** @brief Elements which are being processed */
        std::set<int> fProcessed;
        /** @brief  Current element */
        int64_t fNextElement;
        /** @brief Mutexes (to choose which element is next) */
        std::mutex fMutexAccessElement;
        /** @brief Semaphore (to wake up assembly thread) */
//        std::condition_variable fAssembly;
        TPZSemaphore fAssembly;
    };
};

#endif

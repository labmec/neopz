/**
 * @file
 * @brief Contains the TPZStructMatrixGC class which responsible for a interface among Matrix and Finite Element classes.
 */

#ifndef TPZStructMatrixGC_H
#define TPZStructMatrixGC_H

#include <set>
#include <map>
#include <semaphore.h>
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"
#include "TPZEquationFilter.h"
#include "TPZGuiInterface.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
class TPZStructMatrixGC;
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
class TPZStructMatrixGC : public TPZStructMatrixBase {
    
public:
    
    TPZStructMatrixGC(): TPZStructMatrixBase() {}
    
    TPZStructMatrixGC(TPZCompMesh *);
    
    TPZStructMatrixGC(TPZAutoPointer<TPZCompMesh> cmesh);
    
    TPZStructMatrixGC(const TPZStructMatrixGC &copy);
    
    virtual ~TPZStructMatrixGC(){};
    
    virtual TPZMatrix<STATE> * Create();
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose) {
        SetNumThreads(numthreads_assemble);
        return CreateAssemble(rhs, guiInterface);
    }
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    virtual TPZStructMatrixGC * Clone();
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                          unsigned numthreads_assemble, unsigned numthreads_decompose) {
        std::cout << "Nothing to do." << std::endl;
    }
    
    /** @brief Assemble the global right hand side */
    virtual void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    public:
virtual int ClassId() const;
    void Read(TPZStream& buf, void* context);
    void Write(TPZStream& buf, int withclassid) const;

    
protected:
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Serial_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
public:    
    /** @brief Find the order to assemble the elements */
    static void OrderElement(TPZCompMesh *cmesh, TPZVec<int64_t> &ElementOrder);
    
    /** @brief Create blocks of elements to parallel processing */
    static void ElementColoring(TPZCompMesh *cmesh, TPZVec<int64_t> &elSequence, TPZVec<int64_t> &elSequenceColor, TPZVec<int64_t> &elBlocked);
    
protected:
    
    /** @brief Structure to manipulate thread to solve system equations */
    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrixGC *strmat,TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrixGC *strmat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        /** @brief Destructor: Destroy the mutex semaphores and others */
        ~ThreadData();
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief The function which will compute the assembly */
        bool ShouldCompute(int matid)
        {
            return fStruct->ShouldCompute(matid);
        }
        
        /** @brief Current structmatrix object */
        TPZStructMatrixGC *fStruct;
        /** @brief Gui interface object */
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        /** @brief Global matrix */
        TPZMatrix<STATE> *fGlobMatrix;
        /** @brief Global rhs vector */
        TPZFMatrix<STATE> *fGlobRhs;
        /** @brief List of computed element matrices (autopointers?) */
        std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted;
        /** @brief Elements which are being processed */
        std::set<int> fProcessed;
        /** @brief  Current element */
        int64_t fNextElement;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t fAccessElement;
        /** @brief Semaphore (to wake up assembly thread) */
        TPZSemaphore fAssembly;
        
        pthread_cond_t fCondition;
        bool fSleeping;
        
        // Vectors for mesh coloring
        std::map<int,int> felBlocked;
        /// Vector for mesh coloring
        TPZVec<int64_t> *fnextBlocked, *felSequenceColor;
        
        static void *ThreadWorkResidual(void *datavoid);
    };
    
    friend struct ThreadData;
protected:
    /** @brief Vectors for mesh coloring */
    TPZVec<int64_t> fnextBlocked, felSequenceColor;
};

#endif

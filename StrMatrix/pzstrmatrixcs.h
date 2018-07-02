/**
 * @file
 * @brief Contains the TPZStructMatrixCS class which responsible for a interface among Matrix and Finite Element classes.
 */

#ifndef TPZStructMatrixCS_H
#define TPZStructMatrixCS_H

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

class TPZStructMatrixCS;
#include "TPZStructMatrixBase.h"

#ifdef USING_TBB
#include <tbb/tbb.h>
#endif

/**
 * @brief Refines geometrical mesh (all the elements) num times
 * @ingroup geometry
 */
//void UniformRefine(int num, TPZGeoMesh &m);

/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrixCS : public TPZStructMatrixBase{
    
public:
    
    TPZStructMatrixCS() : TPZStructMatrixBase(){}
    
    TPZStructMatrixCS(TPZCompMesh *);
    
    TPZStructMatrixCS(TPZAutoPointer<TPZCompMesh> cmesh);
    
    TPZStructMatrixCS(const TPZStructMatrixCS &copy);
    
    virtual ~TPZStructMatrixCS(){};
    
    /** @brief Sets number of threads in Assemble process */
    
    virtual TPZMatrix<STATE> * Create();
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose) {
        SetNumThreads(numthreads_assemble);
        return CreateAssemble(rhs, guiInterface);
    }
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    virtual TPZStructMatrixCS * Clone();
    
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
    
protected:
    
    /** @brief Structure to manipulate thread to solve system equations */
    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrixCS *strmat,TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZStructMatrixCS *strmat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        /** @brief Destructor: Destroy the mutex semaphores and others */
        ~ThreadData();
        /** @brief Look for an element index which needs to be computed and put it on the stack */
        int64_t NextElement();
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief Establish whether the element should be computed */
        bool ShouldCompute(int matid)
        {
            return fStruct->ShouldCompute(matid);
        }
        
        /** @brief Current structmatrix object */
        TPZStructMatrixCS *fStruct;
        /** @brief Gui interface object */
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        /** @brief Global matrix */
        TPZMatrix<STATE> *fGlobMatrix;
        /** @brief Global rhs vector */
        TPZFMatrix<STATE> *fGlobRhs;
        /** @brief List of computed element matrices (autopointers?) */
        std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted;
        /** @brief  Current element */
        int64_t fNextElement;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t fAccessElement;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t fAccessElementK;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t fAccessElementF;
        /** @brief Semaphore (to wake up assembly thread) */
        TPZSemaphore fAssembly;
    };
    
#ifdef USING_TBB
    struct AssembleTask {
        ThreadData *data;
        AssembleTask(ThreadData *dt) : data(dt) {};
        void operator()(const tbb::blocked_range<size_t>& range) const;
    };
#endif
    
    friend struct ThreadData;
};

#endif

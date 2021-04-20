/**
 * @file
 * @brief Contains the TPZStructMatrixTBBFlow class which responsible for a interface among Matrix and Finite Element classes using TBB library. Usage of this class without USING_TBB=ON when configuring the NeoPZ library will result in runtime errors.
 */

#ifndef TPZStructMatrixTBBFlow_H
#define TPZStructMatrixTBBFlow_H

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

class TPZStructMatrixTBBFlow;
#include "TPZStructMatrixBase.h"

class TPZFlowGraph;

/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrixTBBFlow : public TPZStructMatrixBase
{
    
public:
    
    TPZStructMatrixTBBFlow();
    
    TPZStructMatrixTBBFlow(TPZCompMesh *);
    
    TPZStructMatrixTBBFlow(TPZAutoPointer<TPZCompMesh> cmesh);
    
    TPZStructMatrixTBBFlow(const TPZStructMatrixTBBFlow &copy);
    
    virtual ~TPZStructMatrixTBBFlow();
    
    virtual TPZBaseMatrix * Create() override;
    
    virtual TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose) {
        SetNumThreads(numthreads_assemble);
        return CreateAssemble(rhs, guiInterface);
    }
    
    virtual TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
    virtual TPZStructMatrixTBBFlow * Clone() override;
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    virtual void Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                          unsigned numthreads_assemble, unsigned numthreads_decompose) {
        std::cout << "Nothing to do." << std::endl;
    }
    
    /** @brief Assemble the global right hand side */
    virtual void Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
    public:
int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;


protected:
    
    //    /** @brief Assemble the global system of equations into the matrix which has already been created */
    //    virtual void Serial_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    //
    //    /** @brief Assemble the global right hand side */
    //    virtual void Serial_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
protected:
    
    /** @brief Pointer to the computational mesh from which the matrix will be generated */
    TPZCompMesh * fMesh;
    /** @brief Autopointer control of the computational mesh */
    TPZAutoPointer<TPZCompMesh> fCompMesh;
    /** @brief Object which will determine which equations will be assembled */
    TPZEquationFilter fEquationFilter;

    TPZFlowGraph *fFlowGraph;
    
protected:
    
    /** @brief Set of material ids to be considered. It is a private attribute. */
    /** Use ShouldCompute method to know if element must be assembled or not    */
    std::set<int> fMaterialIds;
    
    /** @brief Number of threads in Assemble process */
    int fNumThreads;
};

#endif

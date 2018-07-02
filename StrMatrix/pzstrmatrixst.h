/**
 * @file
 * @brief Contains the TPZStructMatrixST class which responsible for a interface among Matrix and Finite Element classes.
 */

#ifndef TPZStructMatrixST_H
#define TPZStructMatrixST_H

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

class TPZStructMatrixST;
#include "TPZStructMatrixBase.h"


/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrixST : public TPZStructMatrixBase {
    
public:
    
    TPZStructMatrixST(): TPZStructMatrixBase(){}
    
    TPZStructMatrixST(TPZCompMesh *);
    
    TPZStructMatrixST(TPZAutoPointer<TPZCompMesh> cmesh);
    
    TPZStructMatrixST(const TPZStructMatrixST &copy);
    
    virtual ~TPZStructMatrixST(){};
        
    /** @brief */
    virtual TPZMatrix<STATE> * Create();
    
    /** @brief */
    virtual TPZStructMatrixST * Clone();

    /** @brief */
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief */
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose) {
        SetNumThreads(numthreads_assemble);
        return CreateAssemble(rhs, guiInterface);
    }
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Assemble(TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief */
    virtual void Assemble(TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                          unsigned numthreads_assemble, unsigned numthreads_decompose) {
        std::cout << "Nothing to do." << std::endl;
    }
    
    /** @brief Assemble the global right hand side */
    virtual void Assemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    public:
virtual int ClassId() const;
    void Read(TPZStream& buf, void* context);
    void Write(TPZStream& buf, int withclassid) const;

    
protected:
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void OnlyAssemble(TPZMatrix<STATE> *mat, TPZFMatrix<STATE> *rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void OnlyAssemble(TPZFMatrix<STATE> *rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief */
    virtual void ExecuteAssemble(TPZMatrix<STATE> *fGlobMatrix, TPZFMatrix<STATE> *fGlobRhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
public:
    
    /** @brief Determine that the assembly refers to a range of equations */
    void SetEquationRange(int64_t mineq, int64_t maxeq)
    {
        fEquationFilter.Reset();
        fEquationFilter.SetMinMaxEq(mineq, maxeq);
    }
    
    /** @brief Verify if a range has been specified */
    virtual bool HasRange() const
    {
        return fEquationFilter.IsActive();
    }
    
    /** @brief access method for the equation filter */
    TPZEquationFilter &EquationFilter()
    {
        return fEquationFilter;
    }
    
    /** @brief number of equations after applying the filter */
    int64_t NReducedEquations() const
    {
        return fEquationFilter.NActiveEquations();
    }
    
    /** @brief Access method for the mesh pointer */
    TPZCompMesh *Mesh() const
    {
        return fMesh;
    }
    
    /** @brief Filter out the equations which are out of the range */
    virtual void FilterEquations(TPZVec<int64_t> &origindex, TPZVec<int64_t> &destindex) const;
    
    /** @brief Set the set of material ids which will be considered when assembling the system */
    void SetMaterialIds(const std::set<int> &materialids);
    
    /** @brief Establish whether the element should be computed */
    bool ShouldCompute(int matid) const
    {
        const size_t size = fMaterialIds.size();
        return size == 0 || fMaterialIds.find(matid) != fMaterialIds.end();
    }
    /** @brief Returns the material ids */
    const std::set<int> &MaterialIds()
    {
        return fMaterialIds;
    }
    
protected:
    
    /** @brief Pointer to the computational mesh from which the matrix will be generated */
    TPZCompMesh * fMesh;
    /** @brief Autopointer control of the computational mesh */
    TPZAutoPointer<TPZCompMesh> fCompMesh;
    /** @brief Object which will determine which equations will be assembled */
    TPZEquationFilter fEquationFilter;
    
protected:
    
    /** @brief Set of material ids to be considered. It is a private attribute. */
    /** Use ShouldCompute method to know if element must be assembled or not    */
    std::set<int> fMaterialIds;
    
    /** @brief Number of threads in Assemble process */
    int fNumThreads;
};

#endif

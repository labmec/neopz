#ifndef TPZSTRUCTMATRIXBASE_H
#define TPZSTRUCTMATRIXBASE_H

#include "TPZEquationFilter.h"
#include "tpzautopointer.h"

class TPZCompMesh;

template <class T> class TPZMatrix;

template <class T> class TPZFMatrix;

class TPZGuiInterface;

class TPZStructMatrixBase : public TPZSavable {
public:
    virtual void SetMesh(TPZCompMesh *);

    virtual void SetMesh(TPZAutoPointer<TPZCompMesh>);

    virtual TPZStructMatrixBase *Clone() = 0;

    virtual TPZMatrix<STATE> * Create() = 0;

    virtual void Assemble(TPZMatrix<STATE> &stiffness, TPZFMatrix<STATE> &rhs,
                          TPZAutoPointer<TPZGuiInterface> guiInterface) = 0;

    virtual void Assemble(TPZFMatrix<STATE> &rhs,
                          TPZAutoPointer<TPZGuiInterface> guiInterface) = 0;

    virtual TPZMatrix<STATE> *
    CreateAssemble(TPZFMatrix<STATE> &rhs,
                   TPZAutoPointer<TPZGuiInterface> guiInterface);

    /** @brief Filter out the equations which are out of the range */
    virtual void FilterEquations(TPZVec<int64_t> &origindex, TPZVec<int64_t> &destindex) const;
    
    /** @brief Set the set of material ids which will be considered when assembling the system */
    virtual void SetMaterialIds(const std::set<int> &materialids);

    inline virtual void SetNumThreads(int n) {
        this->fNumThreads = n;
    }

    inline virtual int GetNumThreads() const {
        return this->fNumThreads;
    }

    inline virtual void SetEquationRange(int64_t mineq, int64_t maxeq) {
        fEquationFilter.Reset();
        fEquationFilter.SetMinMaxEq(mineq, maxeq);
    }

    /** @brief Verify if a range has been specified */
    inline virtual bool HasRange() const {
        return fEquationFilter.IsActive();
    }

    /** @brief access method for the equation filter */
    inline virtual TPZEquationFilter &EquationFilter() {
        return fEquationFilter;
    }

    /** @brief number of equations after applying the filter */
    inline virtual int64_t NReducedEquations() const {
        return fEquationFilter.NActiveEquations();
    }

    /** @brief Access method for the mesh pointer */
    inline virtual TPZCompMesh *Mesh() const {
        return fMesh;
    }
    
    /** @brief Establish whether the element should be computed */
    virtual bool ShouldCompute(int matid) const
    {
        const unsigned int size = fMaterialIds.size();
        return size == 0 || fMaterialIds.find(matid) != fMaterialIds.end();
    }
    /** @brief Returns the material ids */
    virtual const std::set<int> &MaterialIds()
    {
        return fMaterialIds;
    }
    
    public:
virtual int ClassId() const;
    void Read(TPZStream& buf, void* context);
    void Write(TPZStream& buf, int withclassid) const;
    virtual ~TPZStructMatrixBase() {}
  protected:
    TPZStructMatrixBase();
    TPZStructMatrixBase(const TPZStructMatrixBase &);
    TPZStructMatrixBase(TPZCompMesh *);
    TPZStructMatrixBase(TPZAutoPointer<TPZCompMesh>);

  protected:
    /** @brief Pointer to the computational mesh from which the matrix will be
     * generated */
    TPZCompMesh *fMesh;
    /** @brief Autopointer control of the computational mesh */
    TPZAutoPointer<TPZCompMesh> fCompMesh;
    /** @brief Object which will determine which equations will be assembled */
    TPZEquationFilter fEquationFilter;
    /** @brief Set of material ids to be considered. It is a private attribute.
     */
    /** Use ShouldCompute method to know if element must be assembled or not */
    std::set<int> fMaterialIds;
    /** @brief Number of threads in Assemble process */
    int fNumThreads;
};


/** This is the original and stable version of multi_thread_assemble (producer-consumer) */
#include "pzstrmatrix.h"
typedef TPZStructMatrixOR TPZStructMatrix;

/** This version has a clean code with openmp parallelism */
//#include "pzstrmatrixst.h"
//typedef TPZStructMatrixST TPZStructMatrix;

/** This version uses locks in the assemble contribution with tbb (Nathan-Borin) */
//#include "pzstrmatrixcs.h"
//typedef TPZStructMatrixCS TPZStructMatrix;

/** This version uses graph coloring to define the order to process the elements (Devloo-Gilvan) */
//#include "pzstrmatrixgc.h"
//typedef TPZStructMatrixGC TPZStructMatrix;

/** This version uses graph coloring to define the order to process the elements (Devloo-Gilvan) and
 * each color is processed and synchronized */
//#include "pzstrmatrixot.h"
//typedef TPZStructMatrixOT TPZStructMatrix;

/** This version uses graph coloring to define the order to process the elements (Devloo-Quinelato) and
 * each color is processed separately using TPZThreadPool */
//#include "TPZStrMatrixGCTP.h"
//typedef TPZStructMatrixGCTP TPZStructMatrix;

/** This version uses the graph coloring and creates a tbb::flow::graph to process in parallel */
//https://trac.macports.org/wiki/MigrationTBB
//#include "pzstrmatrixtbb.h"
//typedef TPZStructMatrixTBB TPZStructMatrix;

/** This version uses the graph coloring and creates a tbb::flow::graph to process in parallel 
 *  every node of the tbb flow graph computes calc and the assemble
 */
//#include "pzstrmatrixflowtbb.h"
//typedef TPZStructMatrixTBBFlow TPZStructMatrix;

#endif // TPZSTRUCTMATRIXBASE_H

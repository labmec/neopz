/**
    \file TPZStructMatrixBase.h
    Class defining common behaviour for the TPZStructMatrix hierarchy
*/

#ifndef TPZSTRUCTMATRIXBASE_H
#define TPZSTRUCTMATRIXBASE_H

#include "TPZEquationFilter.h"
#include "tpzautopointer.h"

#include <set>

//forward declarations
class TPZCompMesh;
class TPZGuiInterface;
template <class T> class TPZMatrix;
template <class T> class TPZFMatrix;

class TPZStructMatrixBase : public TPZSavable {
public:
    /*Setting the destructor as default would work. however,
      it would request the deletion of the TPZAutoPointer<TPZCompMesh>
      and, since at this point TPZCompMesh is an incomplete type,
      it would lead to a -Wdelete-incomplete warning.
     */
    virtual ~TPZStructMatrixBase();

    /**
     *  Functions to be overriden in child classes
     */
    //@{
    //! Clone method
    virtual TPZStructMatrixBase *Clone() = 0;
    //!Creates a matrix for assembling
    virtual TPZMatrix<STATE> * Create() = 0;
    //! Assemble the global system of equations into a matrix that has already been created.
    virtual void Assemble(TPZMatrix<STATE> &stiffness, TPZFMatrix<STATE> &rhs,
                          TPZAutoPointer<TPZGuiInterface> guiInterface) = 0;
    //! Assemble the global right hand side vector.
    virtual void Assemble(TPZFMatrix<STATE> &rhs,
                          TPZAutoPointer<TPZGuiInterface> guiInterface) = 0;

    virtual TPZMatrix<STATE> *
    CreateAssemble(TPZFMatrix<STATE> &rhs,
                   TPZAutoPointer<TPZGuiInterface> guiInterface);
    //@}

    //! Sets the computational mesh via raw pointer.
    void SetMesh(TPZCompMesh *);
    //! Sets the computational mesh via TPZAutoPointer.
    void SetMesh(TPZAutoPointer<TPZCompMesh>);
    
    //! Filter out the equations which are out of the range
    void FilterEquations(TPZVec<int64_t> &origindex, TPZVec<int64_t> &destindex) const;
    
    //! Set the set of material ids which will be considered when assembling the system
    void SetMaterialIds(const std::set<int> &materialids);
    //! Set number of threads to be used in the assembly.
    inline void SetNumThreads(int n) {
        this->fNumThreads = n;
    }
    //! Get number of threads to be used in the assembly.
    inline int GetNumThreads() const {
        return this->fNumThreads;
    }
    //! Set range of equations to be solved.
    inline void SetEquationRange(int64_t mineq, int64_t maxeq) {
        fEquationFilter.Reset();
        fEquationFilter.SetMinMaxEq(mineq, maxeq);
    }

    //! Verify if a range of equations has been specified.
    inline bool HasRange() const {
        return fEquationFilter.IsActive();
    }

    //! Access method for the equation filter.
    inline TPZEquationFilter &EquationFilter() {
        return fEquationFilter;
    }

    //! Number of equations after applying the filter.
    inline int64_t NReducedEquations() const {
        return fEquationFilter.NActiveEquations();
    }

    //! Access method for the mesh pointer.
    inline TPZCompMesh *Mesh() const {
        return fMesh;
    }
    
    //! Establish whether a given material id should be computed.
    inline bool ShouldCompute(int matid) const
    {
        const unsigned int size = fMaterialIds.size();
        return size == 0 || fMaterialIds.find(matid) != fMaterialIds.end();
    }
    //! Returns the material ids to be computed.
    inline const std::set<int> &MaterialIds()
    {
        return fMaterialIds;
    }
    
    /*! Computes a color for each element.
      \return Number of colors for parallel assembly and the color -1 when the element should not be computed */
    int ComputeElementColors(TPZVec<int> &elementcolors);

    //@{
    //!Read and Write methods
    int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
  protected:
    TPZStructMatrixBase();
    TPZStructMatrixBase(const TPZStructMatrixBase &);
    TPZStructMatrixBase(TPZStructMatrixBase &&) = delete;
    TPZStructMatrixBase& operator=(const TPZStructMatrixBase &) = default;
    TPZStructMatrixBase& operator=(TPZStructMatrixBase &&) = delete;
    TPZStructMatrixBase(TPZCompMesh *);
    TPZStructMatrixBase(TPZAutoPointer<TPZCompMesh>);

    //! Non-managed pointer to the computational mesh from which the matrix will be generated.
    TPZCompMesh *fMesh;
    //! Autopointer of the computational mesh
    TPZAutoPointer<TPZCompMesh> fCompMesh;
    //! Object which will determine which equations will be assembled.
    TPZEquationFilter fEquationFilter;
    //! Set of material ids to be considered in the Assemble process.
    std::set<int> fMaterialIds;
    //! Number of threads in Assemble process.
    int fNumThreads;
};


/** This is the original and stable version of multi_thread_assemble (producer-consumer) */
#include "pzstrmatrix.h"
typedef TPZStructMatrixOR TPZStructMatrix;

/** This version uses graph coloring to define the order to process the elements (Devloo-Gilvan) and
 * each color is processed and synchronized */
//#include "pzstrmatrixot.h"
//typedef TPZStructMatrixOT TPZStructMatrix;

/** This version uses the graph coloring and creates a tbb::flow::graph to process in parallel 
 *  every node of the tbb flow graph computes calc and the assemble
 */
//#include "pzstrmatrixflowtbb.h"
//typedef TPZStructMatrixTBBFlow TPZStructMatrix;

#endif // TPZSTRUCTMATRIXBASE_H

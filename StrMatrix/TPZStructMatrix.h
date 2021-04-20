/**
    \file TPZStructMatrix.h
    Class defining common behaviour for the TPZStructMatrix hierarchy
*/

#ifndef TPZSTRUCTMATRIX_H
#define TPZSTRUCTMATRIX_H

#include "TPZStrMatParInterface.h"
#include "tpzautopointer.h"

#include <set>

//forward declarations
class TPZCompMesh;
class TPZGuiInterface;
class TPZBaseMatrix;

/*!
  Describes the interface that should be implemented for a given structural matrix.

  It is expected that the child classes will be created as:
  template <typename TVar=STATE, typename TPar=TPZStructMatrixOR>
  class TPZStrDerived{};

  Meaning that it will have one template parameter corresponding to its
  type and one template parameter corresponding to its Parallel Strategy,
  a class derived from \ref TPZStrMatParInterface.
*/
class TPZStructMatrix : public virtual TPZStrMatParInterface{
public:
    /**
     *  Methods to be overriden in child classes
     */
    //@{
    //! Clone method
    virtual TPZStructMatrix *Clone() = 0;
    //!Creates a matrix for assembling
    virtual TPZBaseMatrix * Create() = 0;
    //@}
    
    /*
      Other methods
    */

    //@}
    /*Setting the destructor as default would work. however,
      it would request the deletion of the TPZAutoPointer<TPZCompMesh>
      and, since at this point TPZCompMesh is an incomplete type,
      it would lead to a -Wdelete-incomplete warning.
     */
    //! Destructor
    ~TPZStructMatrix();
    //! Sets the computational mesh via raw pointer.
    void SetMesh(TPZCompMesh *);
    //! Sets the computational mesh via TPZAutoPointer.
    void SetMesh(TPZAutoPointer<TPZCompMesh>);
    //! Access method for the mesh pointer.
    inline TPZCompMesh *Mesh() const {
        return fMesh;
    }
    
    //! Set the set of material ids which will be considered when assembling the system
    void SetMaterialIds(const std::set<int> &materialids);
    
    //! Returns the material ids to be computed.
    inline const std::set<int> &MaterialIds()
    {
        return fMaterialIds;
    }

    //! Establish whether a given material id should be computed.
    inline bool ShouldCompute(int matid) const
    {
        const unsigned int size = fMaterialIds.size();
        return size == 0 || fMaterialIds.find(matid) != fMaterialIds.end();
    }

    //! Filter out the equations which are out of the range
    inline void FilterEquations(TPZVec<int64_t> &origindex,
                                TPZVec<int64_t> &destindex) const
    {
        //destindex = origindex;
        fEquationFilter.Filter(origindex, destindex);
    
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
    /* See comment on destructor on why this is not default
     */
    TPZStructMatrix();
    TPZStructMatrix(const TPZStructMatrix &);
    TPZStructMatrix(TPZStructMatrix &&) = delete;
    TPZStructMatrix& operator=(const TPZStructMatrix &) = default;
    TPZStructMatrix& operator=(TPZStructMatrix &&) = delete;
    TPZStructMatrix(TPZCompMesh *);
    TPZStructMatrix(TPZAutoPointer<TPZCompMesh>);

    //! Non-managed pointer to the computational mesh from which the matrix will be generated.
    TPZCompMesh *fMesh{nullptr};
    //! Autopointer of the computational mesh
    TPZAutoPointer<TPZCompMesh> fCompMesh{nullptr};
    //! Set of material ids to be considered in the Assemble process.
    std::set<int> fMaterialIds;
    //! Object which will determine which equations will be assembled.
    TPZEquationFilter fEquationFilter;
};


/**
 * This is the original and stable version of multi_thread_assemble
 * as producer-consumer
 * TPZStructMatrixOR
 * 
 * This version uses graph coloring to define the order to process the
 * elements (Devloo-Gilvan) and each color is processed and synchronized
 * TPZStructMatrixOT
 *
 * This version uses the graph coloring and creates a tbb::flow::graph
 * to process in parallel 
 * every node of the tbb flow graph computes calc and the assemble
 * TPZStructMatrixTBBFlow
 **/
#endif

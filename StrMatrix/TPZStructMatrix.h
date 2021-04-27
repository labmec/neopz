/**
    \file TPZStructMatrix.h
    Class defining common behaviour for the TPZStructMatrix hierarchy
*/

#ifndef TPZSTRUCTMATRIX_H
#define TPZSTRUCTMATRIX_H

#include "TPZStrMatParInterface.h"
#include "tpzautopointer.h"
#include "TPZGuiInterface.h"//avoiding warnings in external projects
#include <set>

//forward declarations
class TPZCompMesh;
class TPZBaseMatrix;

/**
  @brief Describes the type-agnostic interface that should be implemented 
  for a given structural matrix. 
  
  @note  The class \ref TPZStructMatrixT is the one that 
  structural matrices should inherit from.
  @ingroup structural
*/
class TPZStructMatrix : public virtual TPZStrMatParInterface{
public:
    //! Default constructor
    TPZStructMatrix() = default;

    /*
      @orlandini:
      Setting the destructor as default would work.
      However, it would request the deletion of the
      
      TPZAutoPointer<TPZCompMesh>
      
      and, since at this point TPZCompMesh is an incomplete type,
      it would lead to a -Wdelete-incomplete warning.
      The same applies for the constructor taking the
      
      TPZAutoPointer<TPZCompMesh>

     */
    
    //! Destructor
    ~TPZStructMatrix();
    //! Copy constructor
    TPZStructMatrix(const TPZStructMatrix &);
    //! Move constructor (deleted)
    TPZStructMatrix(TPZStructMatrix &&) = delete;
    //! Copy assignment operator
    TPZStructMatrix& operator=(const TPZStructMatrix &) = default;
    //! Move assignment operator
    TPZStructMatrix& operator=(TPZStructMatrix &&) = delete;
    //! Constructor taking the non-managed mesh as a raw pointer
    TPZStructMatrix(TPZCompMesh *);
    //! Constructor taking the mesh as a TPZAutoPointer
    TPZStructMatrix(TPZAutoPointer<TPZCompMesh>);
    
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
    int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
  protected:

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

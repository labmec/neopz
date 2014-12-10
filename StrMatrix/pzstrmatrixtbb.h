/**
 * @file
 * @brief Contains the TPZStructMatrixTBB class which responsible for a interface among Matrix and Finite Element classes.
 */

#ifndef TPZStructMatrixTBB_H
#define TPZStructMatrixTBB_H

#include <set>
#include <map>
#include <semaphore.h>
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"
#include "pzequationfilter.h"
#include "TPZGuiInterface.h"

class TPZCompMesh;
template<class TVar>
class TPZMatrix;
template<class TVar>
class TPZFMatrix;

#ifdef USING_TBB
#include "tbb/tbb.h"
#include "tbb/flow_graph.h"
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
class TPZStructMatrixTBB {
    
public:
    
    TPZStructMatrixTBB(TPZCompMesh *);
    
    TPZStructMatrixTBB(TPZAutoPointer<TPZCompMesh> cmesh);
    
    TPZStructMatrixTBB(const TPZStructMatrixTBB &copy);
    
    virtual ~TPZStructMatrixTBB();
    
    /** @brief Sets number of threads in Assemble process */
    void SetNumThreads(int n){
        this->fNumThreads = n;
    }
    
    int GetNumThreads() const{
        return this->fNumThreads;
    }
    
    virtual TPZMatrix<STATE> * Create();
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose) {
        SetNumThreads(numthreads_assemble);
        return CreateAssemble(rhs, guiInterface);
    }
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    virtual TPZStructMatrixTBB * Clone();
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                          unsigned numthreads_assemble, unsigned numthreads_decompose) {
        std::cout << "Nothing to do." << std::endl;
    }
    
    /** @brief Assemble the global right hand side */
    virtual void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
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
    
    /** @brief Determine that the assembly refers to a range of equations */
    void SetEquationRange(long mineq, long maxeq)
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
    long NReducedEquations() const
    {
        return fEquationFilter.NActiveEquations();
    }
    
    /** @brief Access method for the mesh pointer */
    TPZCompMesh *Mesh() const
    {
        return fMesh;
    }
    
    /** @brief Find the order to assemble the elements */
    static void OrderElement(TPZCompMesh *cmesh, TPZVec<int> &ElementOrder);
    
    /** @brief Create blocks of elements to parallel processing */
    static void ElementColoring(TPZCompMesh *cmesh, TPZVec<int> &elSequence, TPZVec<int> &elSequenceColor, TPZVec<int> &elBlocked);
    
    
    /** @brief Filter out the equations which are out of the range */
    virtual void FilterEquations(TPZVec<long> &origindex, TPZVec<long> &destindex) const;
    
    /** @brief Set the set of material ids which will be considered when assembling the system */
    void SetMaterialIds(const std::set<int> &materialids);
    
    /** @brief Establish whether the element should be computed */
    bool ShouldCompute(int matid) const
    {
        const unsigned int size = fMaterialIds.size();
        return size == 0 || fMaterialIds.find(matid) != fMaterialIds.end();
    }
    /** @brief Returns the material ids */
    const std::set<int> &MaterialIds()
    {
        return fMaterialIds;
    }
    
protected:
    
#ifdef USING_TBB
    struct TPZGraphThreadData {
        // copy constructor
        TPZGraphThreadData(TPZVec<int> fnextBlocked, TPZVec<int> felSequenceColor);
        // destructor
        ~TPZGraphThreadData();
        // tbb tasks graph
        tbb::flow::graph fAssembleGraph;
        // initial node
        tbb::flow::broadcast_node<tbb::flow::continue_msg> fStart;
        // store all the nodes
        std::vector<tbb::flow::continue_node<tbb::flow::continue_msg>* > fGraphNodes;
        // vector for coloring mesh
        TPZVec<int> felSequenceColor;
        
        
        
        /// current structmatrix object
        TPZStructMatrixTBB *fStruct;
        /// gui interface object
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        /// global matrix
        TPZMatrix<STATE> *fGlobMatrix;
        /// global rhs vector
        TPZFMatrix<STATE> *fGlobRhs;
    };
    
    struct TPZGraphThreadNode {
        TPZGraphThreadData *data;
        int iel;
        TPZGraphThreadNode(TPZGraphThreadData *data, int el)
        : data(data), iel(el) {}
        void operator()(tbb::flow::continue_msg) const;
    };
    
    /** @brief Pointer to the tbb::parallel graph structure */
    TPZGraphThreadData *fAssembleThreadGraph;
    
#endif
    

protected:
    
    /** @brief Pointer to the computational mesh from which the matrix will be generated */
    TPZCompMesh * fMesh;
    /** @brief Autopointer control of the computational mesh */
    TPZAutoPointer<TPZCompMesh> fCompMesh;
    /** @brief Object which will determine which equations will be assembled */
    TPZEquationFilter fEquationFilter;
    /** @brief Vectors for mesh coloring */
    TPZVec<int> fnextBlocked, felSequenceColor;
    
protected:
    
    /** @brief Set of material ids to be considered. It is a private attribute. */
    /** Use ShouldCompute method to know if element must be assembled or not    */
    std::set<int> fMaterialIds;
    
    /** @brief Number of threads in Assemble process */
    int fNumThreads;
};

#endif

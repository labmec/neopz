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
#include "TPZEquationFilter.h"
#include "TPZGuiInterface.h"
#include <list>
#include "pzmatrix.h"
#include "pzfmatrix.h"

class TPZStructMatrixTBB;
#include "TPZStructMatrixBase.h"
#ifdef USING_TBB
#include "tbb/tbb.h"
#include "tbb/flow_graph.h"


/**
 * @brief Refines geometrical mesh (all the elements) num times
 * @ingroup geometry
 */
//void UniformRefine(int num, TPZGeoMesh &m);

/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrixTBB : public TPZStructMatrixBase{
    
public:
    
    TPZStructMatrixTBB();
    
    TPZStructMatrixTBB(TPZCompMesh *, bool onlyrhs = false);
    
    TPZStructMatrixTBB(TPZAutoPointer<TPZCompMesh> cmesh, bool onlyrhs = false);
    
    TPZStructMatrixTBB(const TPZStructMatrixTBB &copy);
    
    virtual ~TPZStructMatrixTBB();
    
       
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
    
    public:
    virtual int ClassId() const;
    void Read(TPZStream& buf, void* context);
    void Write(TPZStream& buf, int withclassid) const;
    
protected:
    
//    /** @brief Assemble the global system of equations into the matrix which has already been created */
//    virtual void Serial_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
//    
//    /** @brief Assemble the global right hand side */
//    virtual void Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
protected:
    
    
    class TPZFlowGraph {
        
        struct TPZPointers
        {
            TPZElementMatrix *fEk;
            TPZElementMatrix *fEf;
            
            TPZPointers() : fEk(0), fEf(0){}
        };
        
    public:
        TPZFlowGraph(TPZStructMatrixTBB *strmat, bool onlyrhs);
        ~TPZFlowGraph();
        TPZFlowGraph(TPZFlowGraph const &copy);
        
    protected:
        
        TPZStack<int64_t> fFirstElColor;

        // vectors for mesh coloring
        TPZVec<int64_t> fnextBlocked, felSequenceColor, felSequenceColorInv;
        TPZManVector<int64_t> fElementOrder;
        
        TPZManVector<int,20> fNodeDest;

        
        /// Computational mesh
        TPZCompMesh *fCMesh;
        
        /// tbb graph
        tbb::flow::graph fGraph;
        
        /// pointer to the nodes of the graph
        std::vector<tbb::flow::continue_node<tbb::flow::continue_msg> *> fNodes;
        
        /// variable to store pointers to element matrices
        TPZVec<TPZPointers> fElMatPointers;
        
        /// variable which is true if only the rhs will be assembled
        bool fOnlyRhs;
        
    public:
        
//        std::vector<std::list<int64_t> > nextTasks;
//        std::vector<int64_t> predecessorTasks;
        
        void OrderElements();
        void ElementColoring();
        
        void CreateGraph();
        
        void CreateGraphRhs();
        
        void ExecuteGraph(TPZFMatrix<STATE> *rhs, TPZMatrix<STATE> *matrix);
        
        void ExecuteGraph(TPZFMatrix<STATE> *rhs);

		friend class TPZStructMatrixTBB;
        
    protected:
        /// current structmatrix object
        TPZStructMatrixTBB *fStruct;
        /// gui interface object
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        
        /// matrix for accumulating the rhs during the rhs assembly process
        TPZFMatrix<STATE> fRhsFat;
        /// global matrix
        TPZMatrix<STATE> *fGlobMatrix;
        /// global rhs vector
        TPZFMatrix<STATE> *fGlobRhs;
        
#ifdef USING_TBB
        
        class TPZAssembleTask  {
        public:
            TPZAssembleTask() : fOrigin(0), fElMat(0), fIel(-1){
            }
            
            TPZAssembleTask(int64_t iel, TPZFlowGraph *origin) : fOrigin(origin), fIel(iel)
            {
                fElMat = &origin->fElMatPointers[iel];
            }
            
            ~TPZAssembleTask() {}
            
            TPZAssembleTask(const TPZAssembleTask &copy) : fOrigin(copy.fOrigin), fElMat(copy.fElMat), fIel(copy.fIel)
            {
                
            }
            
            void operator=(const TPZAssembleTask &copy)
            {
                fOrigin = copy.fOrigin;
                fElMat = copy.fElMat;
                fIel = copy.fIel;
            }
            
            tbb::flow::continue_msg operator()(const tbb::flow::continue_msg &msg);
            
            /// FlowGraph Object which spawned the node
            TPZFlowGraph *fOrigin;
            
            TPZPointers *fElMat;
            
            int64_t fIel;
            
            
        };
        
        class TPZCalcTask {
        public:
            
            TPZCalcTask(TPZFlowGraph *flowgraph) : fFlowGraph(flowgraph) {}
            
            void operator()(const tbb::blocked_range<int64_t>& range) const;
            
            
        private:
            TPZFlowGraph *fFlowGraph;
        };
        
        class TAssembleOneColor {
        public:
            TAssembleOneColor(TPZFlowGraph *graph) : fFlowGraph(graph)
            {
            }
            
            void operator()(const tbb::blocked_range<int64_t> &range) const;
            
        private:
            TPZFlowGraph *fFlowGraph;
        };
        
        class TComputeElementRange {
            
        public:
            TComputeElementRange(TPZFlowGraph *graph, int64_t color) : fFlowGraph(graph), fColor(color)
            {
                
            }
            void operator()(const tbb::blocked_range<int64_t> &range) const;
        private:
            TPZFlowGraph *fFlowGraph;
            int64_t fColor;
            
        };
        
        class TSumTwoColors
        {
        public:
            TSumTwoColors(int64_t firstcolumn, int64_t secondcolumn, TPZFMatrix<STATE> *rhs) : fFirstColumn(firstcolumn), fSecondColumn(secondcolumn), fRhs(rhs)
            {
            }
            
            tbb::flow::continue_msg operator()(const tbb::flow::continue_msg &msg) const
            {
                int64_t nrow = fRhs->Rows();
                for (int64_t ir=0; ir<nrow; ir++) {
                    (*fRhs)(ir,fFirstColumn) += (*fRhs)(ir,fSecondColumn);
                }
                return tbb::flow::continue_msg();

            }
        private:
            int64_t fFirstColumn;
            int64_t fSecondColumn;
            TPZFMatrix<STATE> *fRhs;
            
        };

    };

#endif


protected:

#ifdef USING_TBB
    TPZFlowGraph *fFlowGraph;
#endif
};


#endif
#endif

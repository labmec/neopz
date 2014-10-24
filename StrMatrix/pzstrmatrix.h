/**
 * @file
 * @brief Contains the TPZStructMatrix class which responsible for a interface among Matrix and Finite Element classes.
 */

#ifndef TPZSTRUCTMATRIX_H
#define TPZSTRUCTMATRIX_H

#include <set>
#include <map>
#include <semaphore.h>

#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"

#include "pzequationfilter.h"

class TPZCompMesh;
template<class TVar>
class TPZMatrix;
template<class TVar>
class TPZFMatrix;

#include "TPZGuiInterface.h"

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
class TPZStructMatrix {
	
public:    
	
	TPZStructMatrix(TPZCompMesh *);
    
    TPZStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);
	
	TPZStructMatrix(const TPZStructMatrix &copy);
	
	virtual ~TPZStructMatrix();
	
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
        return NULL;
    }
	
	virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	virtual TPZStructMatrix * Clone();
	
	/** @brief Assemble the global system of equations into the matrix which has already been created */
	virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                          unsigned numthreads_assemble, unsigned numthreads_decompose) {
        std::cout << "Nothing to do." << std::endl;
    }
	
	/** @brief Assemble the global right hand side */
	virtual void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    /** @brief Find the order to assemble the elements */
    static void OrderElement(TPZCompMesh *cmesh, TPZVec<int> &ElementOrder);
    
    /** @brief Create blocks of elements to parallel processing */
    
    static void ElementColoring(TPZCompMesh *cmesh, TPZVec<int> &elSequence, TPZVec<int> &elSequenceColor, TPZVec<int> &elBlocked);
    
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
    
    /// access method for the equation filter
    TPZEquationFilter &EquationFilter()
    {
        return fEquationFilter;
    }
    
    /// number of equations after applying the filter
    long NReducedEquations() const
    {
        return fEquationFilter.NActiveEquations();
    }
    
    /// Access method for the mesh pointer
    TPZCompMesh *Mesh() const
    {
        return fMesh;
    }
	
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
	
	/** @brief Structure to manipulate thread to solve system equations */
	struct ThreadData
	{
		/// Initialize the mutex semaphores and others
		ThreadData(TPZStructMatrix *strmat,TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
		/// Initialize the mutex semaphores and others
		ThreadData(TPZStructMatrix *strmat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
		/// Destructor: Destroy the mutex semaphores and others
		~ThreadData();
		/// Current structmatrix object
		TPZStructMatrix *fStruct;
		/// Gui interface object
		TPZAutoPointer<TPZGuiInterface> fGuiInterface;
		/// Mutexes (to choose which element is next)
		pthread_mutex_t fAccessElement;
		/// Semaphore (to wake up assembly thread)
		TPZSemaphore fAssembly;
		/// Global matrix
		TPZMatrix<STATE> *fGlobMatrix;
		/// Global rhs vector
		TPZFMatrix<STATE> *fGlobRhs;
		/// List of computed element matrices (autopointers?)
		std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted;
		/// Elements which are being processed
		std::set<int> fProcessed;
		/// Current element
		long fNextElement;
		/// Look for an element index which needs to be computed and put it on the stack
		long NextElement();
		/// Put the computed element matrices in the map
		void ComputedElementMatrix(long iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef);
		/// The function which will compute the matrices
		static void *ThreadWork(void *threaddata);
		/// The function which will compute the assembly
		static void *ThreadAssembly(void *threaddata);
		
		/// Establish whether the element should be computed
		bool ShouldCompute(int matid)
		{
            return fStruct->ShouldCompute(matid);
		}
        // Vectors for mesh coloring
        std::map<int,int> felBlocked;
        /// Vector for mesh coloring
        TPZVec<int> *fnextBlocked, *felSequenceColor;
        
        // Condition Variables
        pthread_cond_t fCondition;
        bool fSleeping;
        
        static void *ThreadWorkResidual(void *datavoid);
		
	};
#ifdef USING_TBB
    struct GraphThreadData {
        // create tbb::flow::graph
        GraphThreadData(TPZStructMatrix *strmat, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
        // destructor
        ~GraphThreadData();
        // tbb tasks graph
        tbb::flow::graph fAssembleGraph;
        // initial node
        tbb::flow::broadcast_node<tbb::flow::continue_msg> fStart;
        // store all the nodes
        std::vector<tbb::flow::continue_node<tbb::flow::continue_msg>* > fGraphNodes;
        // vector for coloring mesh
        TPZVec<int> fnextBlocked, felSequenceColor;
        
        
        /// current structmatrix object
        TPZStructMatrix *fStruct;
        /// gui interface object
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        /// global matrix
        TPZMatrix<STATE> *fGlobMatrix;
        /// global rhs vector
        TPZFMatrix<STATE> *fGlobRhs;
    };
    
    struct GraphThreadNode {
        GraphThreadData *data;
        int iel;
        GraphThreadNode(GraphThreadData *data, int el)
        : data(data), iel(el) {}
        void operator()(tbb::flow::continue_msg) const;
    };
    
#endif
    
	friend struct ThreadData;
protected:
	
	/** @brief Pointer to the computational mesh from which the matrix will be generated */
	TPZCompMesh * fMesh;
#ifdef USING_TBB
    /** @brief Pointer to the tbb::parallel graph structure */
    GraphThreadData *fAssembleThreadGraph;
#endif
    
    /** @brief Autopointer control of the computational mesh */
    TPZAutoPointer<TPZCompMesh> fCompMesh;
	
    /// Object which will determine which equations will be assembled
    TPZEquationFilter fEquationFilter;
    
    /// Vector for mesh coloring
    TPZVec<int> fnextBlocked, felSequenceColor;
    
protected:
	
	/** @brief Set of material ids to be considered. It is a private attribute. */
	/** Use ShouldCompute method to know if element must be assembled or not */
	std::set<int> fMaterialIds;
	
	/** @brief Number of threads in Assemble process */
	int fNumThreads;
};

#endif

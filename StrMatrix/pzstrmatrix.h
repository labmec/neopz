/**
 * @file
 * @brief Contains the TPZStructMatrix class which responsible for a interface among Matrix and Finite Element classes.
 */
//$Id: pzstrmatrix.h,v 1.16 2011-05-23 19:00:50 fortiago Exp $

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

class TPZCompMesh;
class TPZMatrix;
class TPZFMatrix;

#include "TPZGuiInterface.h"

/**
 * @brief Refines geometrical mesh (all the elements) num times 
 * @ingroup geometry
 */
void UniformRefine(int num, TPZGeoMesh &m);

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
	
	/** @brief Sets number of threads in Assemble process
	 */
	void SetNumThreads(int n){
		this->fNumThreads = n;
	}
	
	int GetNumThreads() const{
		return this->fNumThreads;
	}
	
	/** @brief Calling this method indicates that only internal equations whould be assembled
	 *
	 * This is the default behaviour for objects of type TPZSubCompMesh
	 */
	void AssembleOnlyInternalEquations()
	{
		fOnlyInternal = true;
	}
	
	/** @brief Calling this method indicates all equations should be assembled
	 *
	 * this behaviour is default for objects of type TPZCompMesh
	 */
	void AssembleAllEquations()
	{
		fOnlyInternal = false;
	}
	
	virtual TPZMatrix * Create();
	
	virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	virtual TPZStructMatrix * Clone();
	
	/**
	 * @brief Assemble the global system of equations into the matrix which has already been created
	 */
	virtual void Assemble(TPZMatrix & mat, TPZFMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/**
	 * @brief Assemble the global right hand side
	 */
	virtual void Assemble(TPZFMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
protected:
	
	/**
	 * @brief Assemble the global system of equations into the matrix which has already been created
	 */
	virtual void Serial_Assemble(TPZMatrix & mat, TPZFMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/**
	 * @brief Assemble the global right hand side
	 */
	virtual void Serial_Assemble(TPZFMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	
	/**
	 * @brief Assemble the global right hand side
	 */
	virtual void MultiThread_Assemble(TPZFMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/**
	 * @brief Assemble the global system of equations into the matrix which has already been created
	 */
	virtual void MultiThread_Assemble(TPZMatrix & mat, TPZFMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
public:
	
	/** @brief Determine that the assembly refers to a range of equations */
	void SetEquationRange(int mineq, int maxeq)
	{
		fMinEq = mineq;
		fMaxEq = maxeq;
	}
	
	/** @brief Verify if a range has been specified */
	bool HasRange()
	{
		return (fMinEq != -1 && fMaxEq != -1);
	}
	
	/** @brief Filter out the equations which are out of the range */
	static void FilterEquations(TPZVec<int> &origindex, TPZVec<int> &destindex, int mineq, int upeq);
	
	/** @brief Set the set of material ids which will be considered when assembling the system */
	void SetMaterialIds(const std::set<int> &materialids);
	
	/** @brief Establish whether the element should be computed */
	bool ShouldCompute(int matid)
	{
		const unsigned int size = fMaterialIds.size();
		return size == 0 || fMaterialIds.find(matid) != fMaterialIds.end();
	}
	/** @brief Returns the material ids */
	std::set<int> MaterialIds()
	{
		return fMaterialIds;
	}
	
protected:
	
	/** @brief Structure to manipulate thread to solve system equations */
	struct ThreadData
	{
		/// Initialize the mutex semaphores and others
		ThreadData(TPZCompMesh &mesh,TPZMatrix &mat, TPZFMatrix &rhs, int mineq, int maxeq, std::set<int> &MaterialIds, TPZAutoPointer<TPZGuiInterface> guiInterface);
		/// Destructor: Destroy the mutex semaphores and others
		~ThreadData();
		/// Current structmatrix object
		TPZCompMesh *fMesh;
		/// Gui interface object
		TPZAutoPointer<TPZGuiInterface> fGuiInterface;
		/// Mutexes (to choose which element is next)
		pthread_mutex_t fAccessElement;
		/// Semaphore (to wake up assembly thread)
		TPZSemaphore fAssembly;
		/// Global matrix
		TPZMatrix *fGlobMatrix;
		/// Global rhs vector
		TPZFMatrix *fGlobRhs;
		/// Minimum equation to be assembled
		int fMinEq;
		/// Maximum equation to be assembled
		int fMaxEq;
		/// Material identifiers which need to be computed
		std::set<int> fMaterialIds;
		/// List of computed element matrices (autopointers?)
		std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted;
		/// Elements which are being processed
		std::set<int> fProcessed;
		/// Current element
		int fNextElement;
		/// Look for an element index which needs to be computed and put it on the stack
		int NextElement();
		/// Put the computed element matrices in the map
		void ComputedElementMatrix(int iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef);
		/// The function which will compute the matrices
		static void *ThreadWork(void *threaddata);
		/// The function which will compute the assembly
		static void *ThreadAssembly(void *threaddata);
		
		/// Establish whether the element should be computed
		bool ShouldCompute(int matid)
		{
			return fMaterialIds.size()==0 || fMaterialIds.find(matid) != fMaterialIds.end();
		}
		
	};
	
	friend struct ThreadData;
protected:
	
	/**
	 * @brief Pointer to the computational mesh from which the matrix will be generated
	 */
	TPZCompMesh * fMesh;
    
    /**
     * @brief Autopointer control of the computational mesh
     */
    TPZAutoPointer<TPZCompMesh> fCompMesh;
	
	/** @brief Equation range for assembly of the global matrix
	 */
	int fMinEq, fMaxEq;
	
	bool fOnlyInternal;
	
protected:
	
	/** @brief Set of material ids to be considered. It is a private attribute.
	 *
	 * Use ShouldCompute method to know if element must be assembled or not
	 */
	std::set<int> fMaterialIds;
	
	/** @brief Number of threads in Assemble process
	 */
	int fNumThreads;
};
#endif //TPZSTRUCTMATRIX_H

/**
 * @file
 * @brief Contains the TPZParFrontStructMatrix class which is a structural matrix with parallel techniques included.
 */

#ifndef TPZPARFRONTSTRUCTMATRIX_H
#define TPZPARFRONTSTRUCTMATRIX_H

#include "TPZFrontStructMatrix.h"
#include "pzstack.h"
//#define STACKSTORAGE

/**
 * @brief Is a structural matrix with parallel techniques included. \ref structural "Structural Matrix"
 * @ingroup frontal structural
 */
/** 
 * TPZParFrontStructMatrix is derived from TPZFrontStructMatrix.
 * It uses TPZParFrontMatrix as its FrontalMatrix
 */

template<class TFront,
         class TVar = STATE,
         class TPar = TPZStructMatrixOR<TVar>>
class TPZParFrontStructMatrix :
	public TPZFrontStructMatrix<TFront,TVar,TPar> {
	
private:
	
	TPZAutoPointer<TPZGuiInterface> fGuiInterface;
	
public:     
		
	/** @brief It clones a TPZStructMatrix */
	/** Virtual function must return same type */
	TPZStructMatrix *Clone() override;
	
	/** @brief Constructor passing as parameter a TPZCompMesh */
	TPZParFrontStructMatrix(
							TPZCompMesh *mesh //! Mesh to refer to
							);
	
	int ClassId() const override;
	
	/**
	* @brief Destructor
	*/
	virtual ~TPZParFrontStructMatrix();
	
	/**
	 * @brief Returns a pointer to TPZMatrix
	 * @param rhs Load matrix
	 * @param guiInterface pointer to user interface
	 */
	virtual TPZMatrix<TVar> * CreateAssemble( TPZFMatrix<TVar> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	
	virtual void Assemble(TPZMatrix<TVar> & mat, TPZFMatrix<TVar> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/** @brief It computes element matrices in an independent thread. */
	/** 
	 * It is passed as a parameter to the thread creation function. \n
	 * It is a 'static void *'
	 */
	static void *ElementAssemble(void *t);
	/** @brief It assembles element matrices in the global stiffness matrix, it is also executed in an independent thread. */
	/**
	 * It is passed as a parameter to the thread creation function. \n
	 * It is a 'static void *'
	 */
	static void *GlobalAssemble(void *t);
	
private:
	/** @brief Current computed element*/
	int64_t fCurrentElement;
	/** @brief Current assembled element in the global stiffness matrix*/
	int64_t fCurrentAssembled;
	/** @brief Total number of elements*/
	int64_t fNElements;
	/** @brief Maximum stack size allowed. */
	/** Whenever this value is reached a execution of element computing is suspended */
	int fMaxStackSize;
	/** @brief Local pointer to stiffness matrix*/
	TPZMatrix<TVar> * fStiffness;
	/** @brief Local pointer to load matrix*/
	TPZFMatrix<TVar> * fRhs;
	
	/** @brief Stack containing elements to be assembled on Stiffness matrix. */
	/** 
	 * ElementAssemble pushes elements on the stack. \n
	 * GlobalAssemble pops elements from the stack.
	 */
	TPZStack <int64_t> felnum;
	TPZStack <TPZElementMatrix *> fekstack;
	TPZStack <TPZElementMatrix *> fefstack;
    
    
	
};

#endif //TPZPARFRONTSTRUCTMATRIX_H

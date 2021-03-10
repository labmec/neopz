/**
 * @file
 * @brief Contains the TPZParFrontMatrix class which implements FrontMatrix with parallel techniques.
 */

#ifndef TPZPARFRONTMATRIX_H
#define TPZPARFRONTMATRIX_H
#include "TPZFrontMatrix.h"
#include "TPZFront.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"


#include "TPZStackEqnStorage.h"
#include "pzvec.h"
#include "pzstack.h"
#include "TPZFileEqnStorage.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

struct TPZElementMatrix;
class TPZCompMesh;

template<class TVar>
class TPZEqnArray;


/**
 * @brief FrontMatrix with parallel techniques included. \ref frontal "Frontal"
 * @ingroup frontal
 */
/**
 * As its base class it is also a template class. The parameters store and front can \n
 * assume the values TPZFileEqnStorage or TPZStackEqnStorage for store and TPZFrontSym or TPZFrontNonSym \n
 * for front.
 */ 
template <class TVar, class store, class front>
class TPZParFrontMatrix : public TPZFrontMatrix<TVar, store, front> 
{
    public:
int ClassId() const override;
	/** @brief Used in an independent thread to write decomposed equations to a binary file */
	static void * WriteFile(void *t);
	
	/** @brief Simple Destructor*/
	~TPZParFrontMatrix();
	/** @brief Simple Constructor */
	TPZParFrontMatrix();
	/** 
	 * @brief Constructor with a globalsize parameter 
	 * @param globalsize Indicates initial global size
	 */
	TPZParFrontMatrix(int64_t globalsize);
	
	TPZParFrontMatrix(const TPZParFrontMatrix &cp) : TPZRegisterClassId(&TPZParFrontMatrix::ClassId),
    TPZFrontMatrix<TVar, store,front>(cp), fFinish(0), fwritelock(), fwritecond()
	{
        fEqnStack.Resize(0);
	}
	
	//CLONEDEF(TPZParFrontMatrix)
	virtual TPZMatrix<TVar>*Clone() const override { return new TPZParFrontMatrix(*this); }
	
    /** 
	 * @brief Add a contribution of a stiffness matrix 
	 * @param elmat Member stiffness matrix beeing added
	 * @param sourceindex Source position of values on member stiffness matrix
	 * @param destinationindex Positioning of such members on global stiffness matrix
	 */
    virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec < int64_t > & sourceindex, TPZVec < int64_t > & destinationindex) override;
	
    /**
	 * @brief Add a contribution of a stiffness matrix putting it on destination indexes position
	 * @param elmat Member stiffness matrix beeing added
	 * @param destinationindex Positioning of such members on global stiffness matrix
     */
    virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec < int64_t > & destinationindex) override;
	
	/** @brief Sets the flag fFinish to its true value*/    		
	void FinishWriting();    		
	
private:    
	/** @brief Buffer of pointers to decomposed equations. Stored in a Stack.*/
	TPZStack<TPZEqnArray<TVar> *> fEqnStack;
	/** @brief Boolean responsibility. Assumes values 0 and 1*/
	int fFinish;
	/** @brief Mutual exclusion locks used in management of writeing to disk and decomposition.*/
	std::mutex fwritelock;
	/** @brief Condition variable used in management of writeing to disk and decomposition.*/
	std::condition_variable fwritecond;
	
	/** @link dependency */
	/*#  TPZFrontNonSym lnkTPZFrontNonSym; */ 
	
	/** @link dependency */
	/*#  TPZFrontSym lnkTPZFrontSym; */
};

template <class TVar, class store, class front>
int TPZParFrontMatrix<TVar,store,front>::ClassId() const{
    return Hash("TPZParFrontMatrix") ^ TPZFrontMatrix<TVar, store, front>::ClassId() << 1;
}
#endif //TPZPARFRONTMATRIX_H

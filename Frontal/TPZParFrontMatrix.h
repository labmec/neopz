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

#include <pthread.h>

#include "TPZStackEqnStorage.h"
#include "pzvec.h"
#include "pzstack.h"
#include "TPZFileEqnStorage.h"
#include "pzmatrix.h"

struct TPZElementMatrix;
class TPZCompMesh;

template<class TVar>
class TPZFMatrix;

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
template <class store, class front>
class TPZParFrontMatrix : public TPZFrontMatrix<store, front> 
{
public: 
	
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
	TPZParFrontMatrix(int globalsize);
	
	TPZParFrontMatrix(const TPZParFrontMatrix &cp) : TPZFrontMatrix<store,front>(cp), fFinish(0)
	{
        fEqnStack.Resize(0);
        pthread_mutex_t mlocal = PTHREAD_MUTEX_INITIALIZER;
        fwritelock = mlocal;
        pthread_cond_t clocal = PTHREAD_COND_INITIALIZER;
        fwritecond = clocal;
	}
	
	//CLONEDEF(TPZParFrontMatrix)
	virtual TPZMatrix<REAL>*Clone() const { return new TPZParFrontMatrix(*this); }
	
    /** 
	 * @brief Add a contribution of a stiffness matrix 
	 * @param elmat Member stiffness matrix beeing added
	 * @param sourceindex Source position of values on member stiffness matrix
	 * @param destinationindex Positioning of such members on global stiffness matrix
	 */
    void AddKel(TPZFMatrix<REAL> & elmat, TPZVec < int > & sourceindex, TPZVec < int > & destinationindex);
	
    /**
	 * @brief Add a contribution of a stiffness matrix putting it on destination indexes position
	 * @param elmat Member stiffness matrix beeing added
	 * @param destinationindex Positioning of such members on global stiffness matrix
     */
    virtual void AddKel(TPZFMatrix<REAL> & elmat, TPZVec < int > & destinationindex);
	
	/** @brief Sets the flag fFinish to its true value*/    		
	void FinishWriting();    		
	
private:    
	/** @brief Buffer of pointers to decomposed equations. Stored in a Stack.*/
	TPZStack<TPZEqnArray *> fEqnStack;
	/** @brief Boolean responsibility. Assumes values 0 and 1*/
	int fFinish;
	/** @brief Mutual exclusion locks used in management of writeing to disk and decomposition.*/
	pthread_mutex_t fwritelock;
	/** @brief Condition variable used in management of writeing to disk and decomposition.*/
	pthread_cond_t fwritecond;
	
	/** @link dependency */
	/*#  TPZFrontNonSym lnkTPZFrontNonSym; */ 
	
	/** @link dependency */
	/*#  TPZFrontSym lnkTPZFrontSym; */
};
#endif //TPZPARFRONTMATRIX_H

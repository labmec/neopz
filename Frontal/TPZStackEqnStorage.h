/**
 * @file
 * @brief Contains the TPZStackEqnStorage class responsible for storing arrays of equations.
 */

#ifndef TPZSTACKEQNSTORAGE_H
#define TPZSTACKEQNSTORAGE_H
//#include "TPZFrontMatrix.h"

#include "tpzeqnarray.h"
#include "pzstack.h"

/** 
 * @brief Responsible for storing arrays of equations (mostly in a decomposed form). \ref frontal "Frontal"
 * @ingroup frontal
 */
/** 
 * It has methods for operating over a set of equations
 * The arrays of equations are in the form of a Stack of EqnArrays
 */
class TPZStackEqnStorage {
public:
	
	/** @brief Method to make this class "template compatible" with the file equation storage */
	void ReOpen();
	
	/**  @brief Reinitialize the object */
	void Zero();
    /** @brief It closes the opened binary file. */
	void FinishWriting();
	
    /** @brief Only to make both possible templates similar in terms of methods and constructors */
	TPZStackEqnStorage(char option, const char * name);
    /** @brief Static main for testing */
	static void main();
    /** @brief Simple Destructor */
    ~TPZStackEqnStorage();
    /** @brief Simple Constructor */
    TPZStackEqnStorage();
    
    TPZStackEqnStorage(const TPZStackEqnStorage &cp) : fEqnStack(cp.fEqnStack)
    {
    }
	
    /** 
	 * @brief Adds an EqnArray to EqnStack object
	 * @param *EqnArray Pointer to EqnArray to be added to the Stack
	 */
    void AddEqnArray(TPZEqnArray *EqnArray);
	
    /** 
	 * @brief Prints TPZEqnStorage data. 
	 * @param name file title to print to
	 * @param out object type file
	 */
    void Print(const char *name, std::ostream& out) const;
	
    /** @brief Resets data structure */
    void Reset();
    /** 
	 * @brief Executes a Backward substitution Stack object
	 * @param f Matrix to apply Backward substitution on
	 * @param dec Decomposition type of f, depends on what decomposition method was used to decompose f
	 */
    void Backward(TPZFMatrix<REAL> &f, DecomposeType dec) const ;
	
    /** 
	 * @brief Executes a Forward substitution Stack object 
	 * @param f Matrix to apply Forward substitution on
	 * @param dec Decomposition type of f. Depends on what decomposition method was used to decompose f
	 */
    void Forward(TPZFMatrix<REAL> &f, DecomposeType dec) const;
	
	//Only to make it compatible with FIleEqnStorage

    /**
     * Only to make both possible templates similar in terms of methods and constructors 
     */
	//void SetFileName(const char *name);
	
    /** @brief Only to make both possible templates similar in terms of methods and constructors */
	void OpenGeneric(char option, const char * name);
	
    /** @brief Only to make both possible templates similar in terms of methods and constructors */
	void ReadBlockPositions();
	
	/** @brief Name of Storage */
	std::string GetStorage();
private:
    /** @brief Sets the block size to be used */
    void SetBlockSize();
    /** @brief Defines a stack of EqnArrays */
    TPZStack<TPZEqnArray> fEqnStack;
	
    /** @label Several objects are stored within a stack object
     * @ directed
     * @ link association*/
    /*#  TPZEqnArray lnkTPZEqnArray; */
};

#endif //TPZSTACKEQNSTORAGE_H

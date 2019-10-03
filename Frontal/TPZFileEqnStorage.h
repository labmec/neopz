/**
 * @file
 * @brief Contains the TPZFileEqnStorage class which implements an equation array and stores the EqnArrays.
 */

#ifndef TPZFILEEQNSTORAGE_H
#define TPZFILEEQNSTORAGE_H

#include "tpzeqnarray.h"
#include "pzstack.h"
#include "pzstack.h"
#include "tpzeqnarray.h"
#include "pzfmatrix.h"
/** 
 * @brief Has the same porpouse of EqnStack but stores the EqnArrays in a different form (binary files). \ref frontal "Frontal"
 * @ingroup frontal
 */
/*
 * It has methods for operating over a set of equations. \n
 * The arrays of equations are in the form of a binary files of EqnArrays.
 */
template<class TVar>
class TPZFileEqnStorage : public TPZSavable {
public:
	/** @brief Reopens an binary file with its current fFileName*/
	void ReOpen();
	
	/** @brief Reinitialize the object */
	void Zero();
    /** @brief Method used for binary input/output */
	void ReadBlockPositions();
    /** @brief Method used for binary input/output */
	void FinishWriting();
    /** Sets file name and if it is for input or output, the second term can be either 'r' for input and 'w' for output.*/
	//void SetFileName(const char *name);
	
    /** 
	 * @brief Sets file name and if it is for input or output, the second term can be either 'r' for input and 'w' for output. 
	 * @param option 'w' means writing and 'r' reading
	 * @param name The file name to print to
	 */
	void OpenGeneric(char option, const char * name);
	
    /** @brief Simple constructor */
	TPZFileEqnStorage();
	/** @brief Copy constructor */
	TPZFileEqnStorage(const TPZFileEqnStorage &cp);
	
    /** Static main for testing */
	static void main();
	
	/** @brief Simple destructor */
    ~TPZFileEqnStorage();
    /** 
	 * @brief Constructor with option (can assume "w" or "r") for writing and reading respectively 
	 * @param option 'w' means writing and 'r' reading
	 * @param name the file name to print to
	 */
    TPZFileEqnStorage(char option, const std::string &name);
	
    /** 
	 * @brief Adds an EqnArray 
	 * @param EqnArray added to the binary file
	 */
    void AddEqnArray(TPZEqnArray<TVar> *EqnArray);
	
    /**
     * @brief It prints TPZEqnStorage data.
	 * @param name File name to print to
	 * @param out ofstream object name
     */
    void Print(const char *name, std::ostream& out) const; 
	
    /** @brief Resets data */
    void Reset();
    /** 
	 * @brief Executes a Backward substitution
	 * @param f Full matrix already decomposed
	 * @param dec Decomposition type of f
	 */
    void Backward(TPZFMatrix<TVar> &f, DecomposeType dec) const;
	
    /**
	 * @brief Executes a Forward substitution 
	 * @param f Full matrix already decomposed
	 * @param dec Decomposition type of f
	 */
    void Forward(TPZFMatrix<TVar> &f, DecomposeType dec) const;

    /**
	 * @brief Stores from ieq to jeq equations on a binary file 
	 * @param ieq Initial equation to be added to EqnArray
	 * @param jeq Final equation to be added to EqnArray
	 * @param name Binary file name
	 */
    void Store(int ieq, int jeq, const char *name);    
	
    /** @brief Writes the header of the binary file */
    void WriteHeaders();
	
	/** @brief Type of Storage */
	std::string GetStorage();
	
        public:
        int ClassId() const override;
	
private:
	/** In blocks position */
	//TPZStack<int64_t> fSubBlockIndex;
	
    /** @brief Indicates the number of headers for the object */
	int fNumHeaders;
	
    /** @brief Stack containing block positions */
    TPZStack<int64_t> fBlockPos;

    /** @brief file name containing binary data */
    std::string fFileName;
	
    /** @brief binary file itself */
    FILE *fIOStream;
    /** @brief Used with binary input/output aritimethics */
    int fCurrentBlock;
    /** @brief Used with binary input/output aritimethics */
    int fCurBlockPosition;
    /** @brief Used with binary input/output aritimethics */
    int fNumBlocks;
};

template<class TVar>
int TPZFileEqnStorage<TVar>::ClassId() const{
    return Hash("TPZFileEqnStorage") ^ ClassIdOrHash<TVar>() << 1;
}

#endif //TPZFILEEQNSTORAGE_H

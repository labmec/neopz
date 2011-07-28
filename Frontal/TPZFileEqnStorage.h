/* Generated by Together */

#ifndef TPZFILEEQNSTORAGE_H
#define TPZFILEEQNSTORAGE_H

#include "tpzeqnarray.h"
#include "pzstack.h"
#include "pzstack.h"
#include "tpzeqnarray.h"
#include "pzfmatrix.h"
/** 
 * @brief Has the same purpouse of EqnStack but stores the EqnArrays \n
 * in a different form (binary files). \n
 * 
 * It has methods for operating over a set of equations. \n
 * @ingroup frontal
 */
/**
 * The arrays of equations are in the form of a binary files of EqnArrays.
 */
class TPZFileEqnStorage {
public:
	/** @brief Reopens an binary file with its current fFileName*/
	void ReOpen();
	
	/**
	 * @brief Reinitialize the object
	 */
	void Zero();
    /** @brief Method used for binary input/output */
	void ReadBlockPositions();
    /** @brief Method used for binary input/output */
	void FinishWriting();
    /** Sets file name and if it is for input or output, the second term can be either 'r' for input and 'w' for output.*/
	//void SetFileName(const char *name);
	
    /**
     * @brief Sets file name and if it is for input or output, the second term can be either 'r' for input and 'w' for output. 
     */
	void OpenGeneric(
					 char option //! 'w' means writing and 'r' reading
					 , const char * name //! The file name to print to
					 );
    /** @brief Simple constructor */
	TPZFileEqnStorage();
	
	TPZFileEqnStorage(const TPZFileEqnStorage &cp);
	
    /** Static main for testing */
	static void main();
	
	/** @brief Simple destructor */
    ~TPZFileEqnStorage();
    /** @brief Constructor with
	 option (can assume "w" or "r") for writing and reading respectively
	 */
    TPZFileEqnStorage(
					  char option //! 'w' means writing and 'r' reading
					  , const std::string &name //! the file name to print to
					  );
    /** @brief Adds an EqnArray */
    void AddEqnArray(
					 TPZEqnArray *EqnArray //! EqnArray added to the binary file
					 );
	
    /**
     * @brief It prints TPZEqnStorage data. 
     */
    void Print(
			   const char *name //!File name to print to
			   , std::ostream& out //!ofstream object name
			   ) const; 
    /** @brief Resets data */
    void Reset();
    /** @brief Executes a Backward substitution */
    void Backward(
				  TPZFMatrix &f //!Full matrix already decomposed.
				  , DecomposeType dec //!Decomposition type of f
				  ) const;
	
    /** @brief Executes a Forward substitution */
    void Forward(
				 TPZFMatrix &f //!Full matrix already decomposed.
				 , DecomposeType dec //!Decomposition type of f
				 ) const;
	
    /**
     * @stereotype void 
     */
    /** @brief Stores from ieq to jeq equations on a binary file */
    void Store(
			   int ieq //!Initial equation to be added to EqnArray
			   , int jeq //!Final equation to be added to EqnArray
			   , const char *name //!Binary file name
			   );    
	
    /**
     * @brief Writes the header of the binary file 
     */
    void WriteHeaders();
    /** Sets the block size for writing */
	/* void SetBlockSize(int bs);*/
	
	/**
	 * @brief Type of Storage
	 */
	std::string GetStorage();
	
	
private:
	/** In blocks position */
	//TPZStack<long int> fSubBlockIndex;
	
    /** @brief Indicates the number of headers for the object */
	int fNumHeaders;
    /** Indicates blocksize */
	//int fBlockSize;
	
    /**
     * @brief Stack containing block positions 
     */
    TPZStack<long int> fBlockPos;
	
    /** @label Several objects are stored within a stack object
     * @directed
     * @link association*/
    /*#  TPZEqnArray lnkTPZEqnArray; */
    /** @brief file name containing binary data */
    std::string fFileName;
	
    /**
     * @brief binary file itself 
     */
    FILE *fIOStream;
    /** @brief Used with binary input/output aritimethics */
    int fCurrentBlock;
    /** @brief Used with binary input/output aritimethics */
    int fCurBlockPosition;
    /** @brief Used with binary input/output aritimethics */
    int fNumBlocks;
};

#endif //TPZFILEEQNSTORAGE_H

/* Generated by Together */

//class TPZFMatrix;

#ifndef TPZFRONTMATRIX_H
#define TPZFRONTMATRIX_H


#include "TPZFront.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"


#include "TPZStackEqnStorage.h"
#include "pzvec.h"
#include "TPZFileEqnStorage.h"
#include "pzmatrix.h"

class TPZFMatrix;
/**
 @brief Implements a matrix stored in a frontal decomposition scheme
 */
class TPZAbstractFrontMatrix : public TPZMatrix
{
public:

	TPZAbstractFrontMatrix() : TPZMatrix()
	{
	}
	
	TPZAbstractFrontMatrix(int ieq, int jeq) : TPZMatrix(ieq,jeq)
	{
	}
	
	virtual TPZFront & GetFront() = 0;
		
};
/**
 * @brief Class responsible for the frontal method as a whole. 
 * Manages the remaining classes connecting them
 * @ingroup frontal
 */
template <class store, class front>
class TPZFrontMatrix : public TPZAbstractFrontMatrix {
protected:

    /**
     * Indicates storage schema. Assumes values TPZFileEqnStorage for binary file and TPZStackEqnStorage for a TPZStack storage
     */
     store fStorage;

    /**
     * Indicates Front matrix type. Assumes values TPZFrontSym for symmetric front and TPZFrontNonSym for non symmetric front matrix
     */
     front fFront;
public:
 	int Work();
	/**Finishes writing of a binary file and closes it*/
	void FinishWriting();
	/**Reopen the binary file*/
	void ReOpen();

  /**
   * reinitialize the structure of the frontal matrix
   
   */
   virtual int Zero();
    /** Allocates data for the FrontMatrix */
	void AllocData();

	/** returns a pointer to the front matrix */
	front &GetFront() { return fFront;}
    /** Chacks if FrontMatrix needs a compression,
    if so calls Compress method */
	void CheckCompress();
    /** Static main for testing */
	static void main();
    /** Prints a FrontMatrix object */
	void Print(const char * name, std::ostream & out);
    /** Simple Destructor */
    ~TPZFrontMatrix();
    /** Simple Constructor */
    TPZFrontMatrix();
    /** Constructor with a globalsize parameter */
	TPZFrontMatrix(
		       int globalsize //! Indicates initial global size
		       );
 
 TPZFrontMatrix(const TPZFrontMatrix &cp) : TPZAbstractFrontMatrix(cp), fStorage(cp.fStorage),
 fFront(cp.fFront),fNumEq(cp.fNumEq),fLastDecomposed(cp.fLastDecomposed), fNumElConnected(cp.fNumElConnected),fNumElConnectedBackup(fNumElConnectedBackup)
    {
    }
    
    CLONEDEF(TPZFrontMatrix)

    /**
     * Sends a message to decompose equations from lower_eq to upper_eq,
     * according to destination index
     */
    void EquationsToDecompose(
			      TPZVec<int> &destinationindex //! Contains destination indexes of equations
			      , int &lower_eq //! Starting index
			      , int &upper_eq //! Finishing index
			      );


    /** Add a matrix to the frontal matrix */
//    void AssembleMatrix(TPZVec < int > & eqnumbers, TPZFMatrix & ek, TPZFMatrix & ef);

    /** Initializes the number of elements connected to each equation */
    void SetNumElConnected(
			   TPZVec < int > &numelconnected //! Indicates number of elements connected to that equation
			   );

    /**
     * Add a contribution of a stiffness matrix
	* putting it on destination indexes position
     */
virtual void AddKel(
		TPZFMatrix & elmat //! Member stiffness matrix beeing added
		, TPZVec < int > & destinationindex //! Positioning of such members on global stiffness matrix
		);

    /** Add a contribution of a stiffness matrix using the indexes to compute the frontwidth
      * It does it symbolicaly
      */
    void SymbolicAddKel(
        TPZVec < int > & destinationindex //! Array containing destination indexes
        );

    /** Add a contribution of a stiffness matrix */
    void AddKel(
		TPZFMatrix & elmat //! Member stiffness matrix beeing added
		, TPZVec < int > & sourceindex //! Sorce position of values on member stiffness matrix
		, TPZVec < int > & destinationindex //! Positioning of such members on global stiffness matrix
		);
    /**Forward substitution and result is on b*/
	int Subst_Forward(
        TPZFMatrix *b //! Result of the substitution
        ) const;

    /**Backward substitution and result is on b*/
	int Subst_Backward(TPZFMatrix *b) const;
    /**Executes a substitution on a TPZFmatrix object
       applies both forward and backward substitution automaticaly*/
	int Substitution(TPZFMatrix *) const;
    /*
    void SetFileName(
        const char *name = SetTempFileName() //! Name of the file
        );
    */
    /**
     * Sets a file name for generic input or output
     */
    void SetFileName(
        char option //! It can be either 'w' or 'r' for writeing or reading respectively
        , const char *name //! Name of the file.
        );

protected:
    /** Indicates number of equations */
	int fNumEq;
    /** Indicates last decomposed equation */
	int fLastDecomposed;

    /** @link aggregationByValue */
//    TPZFront fFront;

    /** @link aggregationByValue */
    //TPZStackEqnStorage fStorage;

    /** Contains the number of elements which still need to contribute to a given equation */
    TPZVec <int> fNumElConnected;

    /** Contains the number of elements which still need to contribute to a given equation
     * This backup copy is needed to be able to reinitialize the matrix through the
     * Zero() method
    */
    TPZVec <int> fNumElConnectedBackup;

    /**
     * @link aggregationByValue
     */
    //TPZFileEqnStorage fFileEqnStorage;
};
#endif //TPZFRONTMATRIX_H

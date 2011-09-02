/**
 * @file
 * @brief Contains the TPZFront class which implements decomposition process of the frontal matrix.
 */

#ifndef TPZFRONT_H
#define TPZFRONT_H

#include "pzmatrix.h"
#include "pzstack.h"
#include "pzvec.h"

class TPZEqnArray;

/** 
 * The Front matrix itself. \n
 * It is controled by TPZFrontMatrix.
 */
/**
 * @ingroup frontal
 * @brief Abstract class implements storage and decomposition process of the frontal matrix. \ref frontal "Frontal"
 */
class TPZFront {
	
public:
	
    /** @brief Static main used for testing */	
	static void main();
	
	int NElements();
    /** @brief Simple destructor */
    virtual ~TPZFront();
    /** @brief Simple constructor */
    TPZFront();
    /** @brief Constructor with a initial size parameter */
	TPZFront(
			 int GlobalSize //! Initial size of the Frontal Matrix
			 );
	
	TPZFront(const TPZFront &cp);
	
    /** 
	 * @brief Decompose these equations in a symbolic way and store freed indexes in fFree 
	 * @param mineq Initial equation
	 * @param maxeq Final equation
	 */
    void SymbolicDecomposeEquations(int mineq, int maxeq);
	
	/** 
	 * @brief Add a contribution of a stiffness matrix using the indexes to compute the frontwidth 
	 * @param destinationindex Destination index of each element added
	 */
	void SymbolicAddKel(TPZVec < int > & destinationindex);

	int Work(){
		return fWork;
	}
	
	/** @brief Indicate the first equation dedicated to rigid body modes */
	void SetNumRigidBodyModes(int nrigid)
	{
		fNextRigidBodyMode = fLocal.NElements()-nrigid;
		
		std::cout << " fNextRigidBody Mode neste ponto " << fNextRigidBodyMode<<std::endl;
	}
	
protected:
	int fWork;
private:    
	
    /**
     * @brief Sets the global equation as freed, allowing the space used by this equation to be used
	 * @param global Equation number to be freed
	 */
	/** By future assembly processes. 
     */
    void FreeGlobal(int global);
    
	/** 
	 * @brief return a local index corresponding to a global equation number 
	 * @param global Global equation index which has a local indexation
     */
    int Local(int global);
	
public:
	/** @brief Extracts the so far condensed matrix */
	virtual	void ExtractFrontMatrix(TPZFMatrix &front) {
		std::cout << "TPZFront ExtractFrontMatrix should never be called\n";
		DebugStop();
	}
    /** @brief Returns the number of free equations */
	virtual int NFree();
    /** Resets data structure */
	void Reset(int GlobalSize=0);
	
    /** @brief It prints TPZFront data */
	void Print(const char *name, std::ostream& out) const;
	void PrintGlobal(const char *name, std::ostream& out = std::cout);
	
	/** @brief returns the actual front size */
	int FrontSize()
	{
		return fFront;
	}

protected:
    
    /** @brief Maximum size of the front */
    int fMaxFront;

    /**
     * @brief Global equation associated to each front equation.
	 */
	/**
	 * If we need a position in globalmatrix of a equation "i" in the frontmatrix \n
	 * then we can use fGlobal[i]. If the global equation "i" is not used \f$ then fGlobal[i]==-1 \f$
     */
    TPZManVector <int> fGlobal;
	
    /** @brief Front equation to each global equation */
    /**
	 * If we need a position in frontmatrix of a global equation "i" \n
	 * then we can use fLocal[i]. If the global equation is not represented in the front then \f$ fLocal[i]==-1 \f$.
     */
    TPZVec<int> fLocal;
	
    /** @brief Actual front size */
    int fFront;
	
	/** @brief Equation where rigid body modes can be stored */
	int fNextRigidBodyMode;
	
    /** @brief Colection of already decomposed equations still on the front */
    TPZStack <int> fFree;
	
    /** @brief Frontal matrix data */
    TPZVec<REAL> fData;
    
    /** @brief Expansion Ratio of frontal matrix */
    int fExpandRatio;
};

#endif //TPZFRONT_H

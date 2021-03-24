/**
 * @file
 * @brief Contains the TPZFrontSym class which implements decomposition process of the frontal matrix (case symmetric).
 */
template<class TVar>
class TPZEqnArray;

#ifndef TPZFRONTSYM_H
#define TPZFRONTSYM_H

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "pzmatrix.h"
#include "pzstack.h"
#include "pzvec.h"
#include "TPZFront.h"
#include "TPZFileEqnStorage.h"
#include "TPZStackEqnStorage.h"

/** 
 * The Front matrix itself.
 * It is controled by TPZFrontMatrix.\n
 * TPZFrontSym is a symmetrical matrix. It uses a Cholesky decomposition scheme.
 */
/**
 * @brief Abstract class implements storage and decomposition process of the frontal matrix involving simmetry characteristics. \ref frontal "Frontal"
 * @ingroup frontal
 */
template <class TVar>
class TPZFrontSym : public TPZFront<TVar> {
public:
	/** @brief Returns its type*/
	std::string GetMatrixType();
	
    /** Static main used for testing */
	static void main();
	/** @brief Simple destructor */
    ~TPZFrontSym();
    /** @brief Simple constructor */
    TPZFrontSym();
    
    TPZFrontSym(const TPZFrontSym<TVar> &cp) : TPZRegisterClassId(&TPZFrontSym<TVar>::ClassId),
    TPZFront<TVar>(cp)
    {
    }
    /** @brief Constructor with a initial size parameter */
	TPZFrontSym(int64_t GlobalSize);
	
    /// Set the decomposition type
    virtual void SetDecomposeType(DecomposeType dectype) override
    {
        if (dectype == ECholesky || dectype == ELDLt) {
            this->fDecomposeType = dectype;
        }
        else
        {
            DebugStop();
        }
    }
    

    /**
     * @brief Decompose these equations and put the result in eqnarray \n
     * Default decompose method is Cholesky
	 * @param mineq index of equations to be decomposed
	 * @param maxeq index of equations to be decomposed
	 * @param result result of decomposition
     */
    void DecomposeEquations(int64_t mineq, int64_t maxeq, TPZEqnArray<TVar> & result);
	
    /**
     * @brief Decompose these equations in a symbolic way and store freed indexes in fFree 
	 * @param mineq Initial equation index
	 * @param maxeq Final equation index
     */
    void SymbolicDecomposeEquations(int64_t mineq, int64_t maxeq);
	
	/** @brief Add a contribution of a stiffness matrix using the indexes to compute the frontwidth */
	void SymbolicAddKel(TPZVec < int64_t > & destinationindex);
	
    /** @brief Compress data structure */
    void Compress();
	
	/** @brief Expand the front matrix */
	void Expand(int largefrontsize);
	
    /** @brief Returns ith, jth element of matrix. \f$ (sourceindex[i],sourceindex[j]) \f$ */
	TVar & Element(int64_t i, int64_t j){
		if(i>j){
			int64_t i_temp=i;
			i=j;
			j=i_temp;
		}
		return this->fData[(j*(j+1))/2+i];
	}
    /** @brief Returns ith, jth element of matrix. \f$ (sourceindex[i],sourceindex[j]) \f$ */
    const TVar & Element(int64_t i, int64_t j) const {
        if(i>j){
            int64_t i_temp=i;
            i=j;
            j=i_temp;
        }
        return this->fData[(j*(j+1))/2+i];
    }
    /** @brief Add a contribution of a stiffness matrix*/
    void AddKel(TPZFMatrix<TVar> &elmat, TPZVec<int64_t> &destinationindex);
	
    /**@brief Add a contribution of a stiffness matrix*/
    virtual void AddKel(TPZFMatrix<TVar> &elmat, TPZVec<int64_t> &sourceindex,  TPZVec<int64_t> &destinationindex);
	
	/** @brief Reorders the elements of the frontmatrix into the full matrix */
	virtual void ExtractFrontMatrix(TPZFMatrix<TVar> &front) override;
        
        public:
int ClassId() const override;
	
private:    

	TVar & Element4JGreatEqualI(int64_t i, int64_t j){
#ifdef PZDEBUG
    if(i>j){
      DebugStop();
    }
#endif
	  return this->fData[(j*(j+1))/2+i];
  }
	
	
	/** 
	 * @brief Decomposes ieq equation and add the result to EqnArray 
	 * @param ieq Index of equation to be decomposed 
	 * @param eqnarray EqnArray to store resulting members
	 */
    void DecomposeOneEquation(int64_t ieq, TPZEqnArray<TVar> &eqnarray);
	
    /**
     * @brief Sets the global equation as freed, allowing the space 
	 */
	/** 
	 * Used by this equation to be used by future assembly processes 
     */
    void FreeGlobal(int64_t global);
    /** @brief return a local index corresponding to a global equation number */
    int Local(int64_t global);
public:
    /** @brief Returns the number of free equations */
	virtual int64_t NFree() override;
    /** @brief Resets data structure */
	void Reset(int64_t GlobalSize=0);
    /** @brief Allocates data for Front */
	void AllocData();
	
	/** @brief Prints TPZFront data */
	void Print(const char *name, std::ostream& out=std::cout) const;
	void PrintGlobal(const char *name, std::ostream& out = std::cout);
	
	/** @brief Returns decomposition type*/
	DecomposeType GetDecomposeType() const;
	
	/** @brief Does the tensor product betweem two vectors in the positions dependent of ithread*/
	virtual void TensorProductIJ(int ithread, typename TPZFront<TVar>::STensorProductMTData *data) override;
	
    /** @link dependency */
    /*#  TPZFileEqnStorage lnkTPZFileEqnStorage; */
	
    /** @link dependency */
    /*#  TPZStackEqnStorage lnkTPZStackEqnStorage; */
};

#endif //TPZFRONTSYM_H

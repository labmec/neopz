/* Generated by Together */
class TPZEqnArray;

#ifndef TPZFRONTNONSYM_H
#define TPZFRONTNONSYM_H


#include "pzmatrix.h"
#include "pzstack.h"
#include "pzvec.h"
#include "TPZFront.h"

#include <math.h>
#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
#include <fstream>
#include "TPZStackEqnStorage.h"
#include "TPZFileEqnStorage.h"

#ifdef USING_BLAS
extern "C"{
  #include "cblas.h"
//#include "g2c.h"
//#include "fblaswr.h"
};
#endif

#ifdef USING_ATLAS
extern "C"{
  #include <cblas.h>
//#include "g2c.h"
//#include "fblaswr.h"
};
#endif


/**
 * The Front matrix itself. \n
 * It is controled by TPZFrontMatrix. \n
 * TPZFrontNonSym is a non symmetrical matrix. It uses LU decomposition scheme.
 * @ingroup frontal
 */
/** Jorge
 * @brief Abstract class implements storage and decomposition process of the frontal matrix involving non-simetry characteristics
 */
class TPZFrontNonSym : public TPZFront {
public:
     /**
      *Type of matrix
      */
     std::string GetMatrixType();

    /** Static main used for testing */
				 static void main();
    /** Simple destructor */
    ~TPZFrontNonSym();
    /** Simple constructor */
    TPZFrontNonSym();
    /** Constructor with a initial size parameter */
        TPZFrontNonSym(int GlobalSize);

        TPZFrontNonSym(const TPZFrontNonSym &cp);
    /**
     * Decompose these equations and put the result in eqnarray. \n
     * Default decompose method is LU
     */
    void DecomposeEquations(
                            int mineq //! Starting index of equations to be decomposed
                            , int maxeq //! Finishing index of equations to be decomposed
                            , TPZEqnArray & result //! Result of decomposition
                            );

    /**
     * Decompose these equations in a symbolic way and store freed indexes in fFree 
     */
    void SymbolicDecomposeEquations(int mineq, int maxeq);
        
        /** Add a contribution of a stiffness matrix using the indexes to compute the frontwidth */
        void SymbolicAddKel(TPZVec < int > & destinationindex);

    /**
     * Compress data structure 
     */
    void Compress();

        /**
        * Expand the front matrix
        */
        void Expand(
             int largefrontsize//! New size of front
             );

    /**
     * Returns the ith,jth element of the matrix.
     * @associates <{mat(sourceindex[i],sourceindex[j])}>
     * @semantics += 
     */
     REAL & Element(int i, int j){
             return fData[fMaxFront*j + i];
     }
    
    /**
     * Returns the ith,jth element of the matrix.
     * @associates <{mat(sourceindex[i],sourceindex[j])}>
     * @semantics += 
     */
    const REAL & Element(int i, int j) const{
        return fData[fMaxFront*j + i];
    }
    
    /**Add a contribution of a stiffness matrix*/
    void AddKel(
          TPZFMatrix &elmat //! Already formed element matrix
          , TPZVec<int> &destinationindex //! Destine index on the global matrix
          );

    /**Add a contribution of a stiffness matrix*/
    void AddKel(
          TPZFMatrix &elmat //! Already formed element matrix
          , TPZVec<int> &sourceindex //! Source index 
          ,  TPZVec<int> &destinationindex //! Destine index on the global matrix
          );    

	/// Extract the front matrix
	virtual void ExtractFrontMatrix(TPZFMatrix &front);
	
private:    
	
    /**
     * Decomposes ieq equation and add the result to EqnArray 
     */
    void DecomposeOneEquation(
                              int ieq //! index of equation to be decomposed
                              , TPZEqnArray &eqnarray //! EqnArray to store resulting members
                              );
    /**
     * Sets the global equation as freed, allowing the space \n
     * used by this equation to be used by future assembly processes. 
     */
    void FreeGlobal(
          int global //! global index to be freed.
          );
    /**
     * return a local index corresponding to a global equation number 
     */
    int Local(
          int global //! global index inquired.
          );
public:
    /** Returns the number of free equations */
        virtual int NFree();
    /** Resets data structure */
        void Reset(
             int GlobalSize=0 //! Initial global size to be used in reseting.
             );
    /** Allocates data for Front */
        void AllocData();

    /**
     * It prints TPZFront data 
     */
     void Print(const char *name, std::ostream& out) const;
     void PrintGlobal(const char *name, std::ostream& out);
     /**Returns decomposition type. \n Default LU*/
     DecomposeType GetDecomposeType() const;



private:
    /**
     * Used Decomposition method 
     */
    DecomposeType fDecomposeType;

    /** @link dependency */
    /*#  TPZStackEqnStorage lnkTPZStackEqnStorage; */

    /** @link dependency */
    /*#  TPZFileEqnStorage lnkTPZFileEqnStorage; */
};
#endif //TPZFRONTNONSYM_H

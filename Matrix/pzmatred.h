/**
 * @file
 * @brief Contains TPZMatRed class which implements a simple substructuring of a linear system of equations, composed of 4 submatrices.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tmatred.h
//
// Class:  TPZMatRed
//
// Obs.:   Subestruturacao simples de um sistema de equacoes.
//
//				[K00][U0] + [K01][U1] = [F0]
//				[K10][U0] + [K11][U1] = [F1]
//
// Versao: 04 / 1996.
//

#ifndef _TMATREDHH_
#define _TMATREDHH_

//#include "tintvec.h"

/** @brief Id of TPZMATRED matrix very sparse */
#define TPZMATRED_VERYSPARSE_ID 28291103
/** @brief If of TPZMATRED full matrix */
#define TPZMATRED_FMATRIX_ID 28291102

#include "pzmatrix.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzsolve.h"
#include "tpzverysparsematrix.h"

#ifdef OOPARLIB
#include "pzsaveable.h"
#include "pzmatdefs.h"
#endif

template<class TVar>
class TPZFMatrix;

/**
 * @brief Implements a simple substructuring of a linear system of equations, composed of 4 submatrices. \ref matrix "Matrix"
 * @ingroup matrix
 */
/**
 * Implements a matrix composed of 4 submatrices:
 *			\f[ [K00] [U0] + [K01] [U1] = [F0] \f]
 *			\f[ [K10] [U0] + [K11] [U1] = [F1] \f]
 */

template<class TVar, class TSideMatrix >
class TPZMatRed: public TPZMatrix<TVar>
{
public:
	
	friend class TPZMatRed<TVar, TPZFMatrix<TVar> >;
	friend class TPZMatRed<TVar ,TPZVerySparseMatrix<TVar> >;
	/** @brief Simple constructor */
	TPZMatRed();
	
	/**
	 * @brief Constructor with 2 parameters
	 * @param dim assumes the value of n1+n2
	 * @param dim00 equals n1
	 */
	TPZMatRed(const int dim, const int dim00);
	
	template<class TSideCopy>
	TPZMatRed<TVar ,TSideMatrix>(const TPZMatRed<TVar, TSideCopy> &cp): TPZMatrix<TVar>(cp), fK11(cp.fK11), fK01(cp.fK01), fK10(cp.fK10), fF0(cp.fF0), fF1(cp.fF1),fMaxRigidBodyModes(cp.fMaxRigidBodyModes),fNumberRigidBodyModes(cp.fNumberRigidBodyModes)
	{
		fDim0=cp.fDim0;
		fDim1=cp.fDim1;
		fF0IsComputed=cp.fF0IsComputed;
        fK00IsDecomposed = cp.fK00IsDecomposed;
		fK11IsReduced=cp.fK11IsReduced;
		fK01IsComputed = cp.fK01IsComputed;
		fF1IsReduced=cp.fF1IsReduced;
		fIsReduced = cp.fIsReduced;
		fSolver = cp.fSolver;
		
		if(cp.fK00) fK00 = cp.fK00;
	}
	
	CLONEDEF(TPZMatRed)
	/** @brief Simple destructor */
	~TPZMatRed();
	
	/** @brief returns 1 or 0 depending on whether the fK00 matrix is zero or not */
	virtual int IsSimetric() const;
	
	/** @brief changes the declared dimension of the matrix to fDim1 */
	void SetReduced()
	{
		TPZMatrix<TVar>::Resize(fDim1, fDim1);
		fIsReduced = 1;
	}
	
	/**
	 * @brief Put and Get values without bounds checking
	 * these methods are faster than "Put" e "Get" if DEBUG is defined
	 */
	virtual int PutVal(const int row, const int col, const TVar& value);
	virtual const TVar &GetVal(const int row, const int col) const;
	virtual TVar &s(int row, int col);
	
	/** @brief This method will zero all submatrices associated with this reducable matrix class */
	virtual int Zero();
	
	/**
	 * @brief Sets the matrix pointer of the upper left matrix to K00
	 * @param K00 pointer to an upper left matrix
	 */
	void SetK00(TPZAutoPointer<TPZMatrix<TVar> > K00);
	
	TPZAutoPointer<TPZMatrix<TVar> > K00()
	{
		return fK00;
	}
	TSideMatrix &K01()
	{
		return fK01;
	}
	TSideMatrix &K10()
	{
		return fK10;
	}
	void SetSolver(TPZAutoPointer<TPZMatrixSolver<TVar> > solver);
    TPZAutoPointer<TPZMatrixSolver<TVar> > Solver()
    {
        return fSolver;
    }
	/**
	 * @brief Copies the F vector in the internal data structure
	 * @param F vector containing data to stored in current object
	 */
	void SetF(const TPZFMatrix<TVar> & F);
    
    /** @brief indicate how many degrees of freedom are reserved for rigid body modes */
    void SetMaxNumberRigidBodyModes(int maxrigid)
    {
        fMaxRigidBodyModes = maxrigid;
    }
    
    /** @brief return the number of rigid body modes detected during decomposition */
    int NumberRigidBodyModes()
    {
        return fNumberRigidBodyModes;
    }
	
	/** @brief Indicate whether F0 needs to be reduced or not */
	void SetF0IsComputed(bool value)
	{
		fF0IsComputed = value;
	}
	/** @brief Indicate that the value of F1 has been reduced */
	void SetF1IsReduced(bool value)
	{
		fF1IsReduced = value;
	}
	
	/** @brief Computes the reduced version of the right hand side \f$ [F1]=[F1]-[K10][A00^-1][F0] \f$ */
	const TPZFMatrix<TVar> & F1Red();
	
	/** @brief Computes the K11 reduced \f$ [K11]=[K11]-[K10][A00^-1][A01] \f$ */
	const TPZFMatrix<TVar> & K11Red();
	
	/**
	 * @brief Returns the second vector, inverting K11
	 * @param F contains second vector
	 */
	void U1(TPZFMatrix<TVar> & F);
	
	/**
	 * @brief Computes the complete vector based on the solution u1.
	 * @param U1 right hand side
	 * @param result contains the result of the operation
	 */
	void UGlobal(const TPZFMatrix<TVar> & U1, TPZFMatrix<TVar> & result);
	void UGlobal2(TPZFMatrix<TVar> & U1, TPZFMatrix<TVar> & result);
	//  TPZFMatrix<>U0(TPZMatrix<>*u1 = NULL);
	
	
	/** @brief Prints the object data structure */
	void Print(const char *name = NULL, std::ostream &out = std::cout,
			   const MatrixOutputFormat = EFormatted) const;
	
	
	/** @brief Redim: Set the dimension of the complete matrix and reduced matrix */
	int Redim(int dim, int dim00); //Cesar 19/12/00
	
	
	/**
	 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
	 * @param x Is x on the above operation
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 * @param stride Indicates n/N where n is dimension of the right hand side
	 * vector and N is matrix dimension
	 */
	void MultAdd(const TPZFMatrix<TVar> &x, const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
				 const TVar alpha, const TVar beta, const int opt, const int stride) const;
	
	/** @brief If fK00 is simetric, only part of the matrix is accessible to external objects. */
	/** Simetrizes copies the data of the matrix to make its data simetric */
	void Simetrize();
	
	/**
	 * template class TPZMatRed<TPZVerySparseMatrix>;
	 template class TPZMatRed<TPZFMatrix>;
	 *
	 */
	
	/** @brief Saveable methods */
	int ClassId() const;
	
	//TPZMATRED_FMATRIX_ID
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
	
private:
	
	//static int Error(const char *msg ,const char *msg2 = "");
	
	/**
	 * @brief Swaps the row and column index
	 * @param row Row number
	 * @param col Column number
	 */
	static void Swap(int *row, int *col);
    
    /** @brief Decompose K00 and adjust K01 and K10 to reflect rigid body modes */
    void DecomposeK00();
	
	/** @brief Stiffnes matrix */
	TPZAutoPointer<TPZMatrix<TVar> > fK00;
	
	/** @brief Solution method for inverting \f$ fK00 \f$ */
	TPZAutoPointer<TPZMatrixSolver<TVar> > fSolver;
	/** @brief Full Stiffnes matrix */
	TSideMatrix fK11;
	TSideMatrix fK01, fK10;
	
	/** @brief Right hand side or force matrix */
	TPZFMatrix<TVar> fF0, fF1;

	/** @brief Stores matricess \f$ fKij \f$ dimensions */
	int fDim0, fDim1;
	
	/** @brief Is true if the declared dimension of the matrix is fDim0 */
	char fIsReduced;
	
	/** @brief Is true if \f$ [(K00)^-1][F0] \f$ has been calculated and overwritten on \f$ [F0] \f$ */
	char fF0IsComputed;
    
    /** @brief Flag indicating that K00 was decomposed */
    char fK00IsDecomposed;
	
	/** @brief Is true if \f$ [(K00)^-1][KO1] \f$ has been computed and overwritten \f$ [K01] \f$ */
	char fK01IsComputed;
	
	/** @brief fK11IsReduced is true if \f$ [K11]=[K11]-[K10][(A00)^-1][A01] \f$ exists */
	char fK11IsReduced;
	
	/** @brief fF1IsReduced is true if  \f$ [F1]=[F1]-[K10][(A00)^-1][F0] \f$ exists */
	char fF1IsReduced;
    
    /** @brief Number of rigid body modes foreseen in the computational mesh */
    int fMaxRigidBodyModes;
    
    /** @brief Number of rigid body modes identified during the decomposition of fK00 */
    int fNumberRigidBodyModes;
};

/************/
/*** Swap ***/
/* @note Modificacao por Philippe Devloo (insercao de inline )*/
template<class TVar, class TSideMatrix>
inline void TPZMatRed<TVar, TSideMatrix>::Swap(int *a, int *b)
{
	int aux = *a;
	*a = *b;
	*b = aux;
}


#endif

/**
 * @file
 * @brief Contains TPZMatRed class which implements a simple substructuring of a linear system of equations, composed of 4 submatrices.
 */

#ifndef _TMATREDHH_
#define _TMATREDHH_

//#include "tintvec.h"
#include "pzmatrixid.h"


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
	TPZMatRed(const long dim, const long dim00);
	
	template<class TSideCopy>
	TPZMatRed<TVar ,TSideMatrix>(const TPZMatRed<TVar, TSideCopy> &cp): TPZMatrix<TVar>(cp), fK11(cp.fK11), fK01(cp.fK01), fK10(cp.fK10), fF0(cp.fF0), fF1(cp.fF1),fMaxRigidBodyModes(cp.fMaxRigidBodyModes),fNumberRigidBodyModes(cp.fNumberRigidBodyModes)
	{
		fDim0=cp.fDim0;
		fDim1=cp.fDim1;
		fK01IsComputed = cp.fK01IsComputed;
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
	
    void ReallocSolver() {
        fSolver->ReallocMatrix();
    }
	/**
	 * @brief Put and Get values without bounds checking
	 * these methods are faster than "Put" e "Get" if DEBUG is defined
	 */
	virtual int PutVal(const long row, const long col, const TVar& value);
	virtual const TVar &GetVal(const long row, const long col) const;
	virtual TVar &s(const long row, const long col);
	
	/** @brief This method will zero all submatrices associated with this reducable matrix class */
	virtual int Zero();
	
	/**
	 * @brief Sets the matrix pointer of the upper left matrix to K00
	 * @param K00 pointer to an upper left matrix
	 */
	void SetK00(TPZAutoPointer<TPZMatrix<TVar> > K00);
    
    /**
     * @brief Sets K01 as computed
     */
    
    void SetK01IsComputed (int n)
    {
        fK01IsComputed=n;
    }
	
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
    
    long Dim0()
    {
        return fDim0;
    }
    
    long Dim1()
    {
        return fDim1;
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
	
	
	/** @brief Computes the reduced version of the right hand side \f$ [F1]=[F1]-[K10][A00^-1][F0] \f$ */
	void F1Red(TPZFMatrix<TVar> &F1);
	
	/** @brief Computes the K11 reduced \f$ [K11]=[K11]-[K10][A00^-1][A01] \f$ */
	void K11Reduced(TPZFMatrix<TVar> &K11, TPZFMatrix<TVar> &F1);
	
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
	
	/** @brief Prints the object data structure */
	void Print(const char *name = NULL, std::ostream &out = std::cout,
			   const MatrixOutputFormat = EFormatted) const;
	
	/** @brief Redim: Set the dimension of the complete matrix and reduced matrix */
	int Redim(const long dim,const long dim00); //Cesar 19/12/00
	
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
	void SimetrizeMatRed();
	
	/** @brief Saveable methods */
	int ClassId() const;
	
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
	
private:
	
	/**
	 * @brief Swaps the row and column index
	 * @param row Row number
	 * @param col Column number
	 */
	static void Swap(long *row, long *col);
    
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
	long fDim0, fDim1;
	
	/** @brief Is true if the declared dimension of the matrix is fDim0 */
	bool fIsReduced;
	
	/** @brief Is true if \f$ [(K00)^-1][KO1] \f$ has been computed and overwritten \f$ [K01] \f$ */
	bool fK01IsComputed;
	
    /** @brief Number of rigid body modes foreseen in the computational mesh */
    int fMaxRigidBodyModes;
    
    /** @brief Number of rigid body modes identified during the decomposition of fK00 */
    int fNumberRigidBodyModes;
};

/************/
/*** Swap ***/
/* @note Modificacao por Philippe Devloo (insercao de inline )*/
template<class TVar, class TSideMatrix>
inline void TPZMatRed<TVar, TSideMatrix>::Swap(long *a, long *b)
{
	long aux = *a;
	*a = *b;
	*b = aux;
}

#endif

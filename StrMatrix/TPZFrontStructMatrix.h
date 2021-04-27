/**
 * @file
 * @brief Contains the TPZFrontStructMatrix class which responsible for a interface among Finite Element Package and Matrices package to frontal method.
 */

#ifndef TPZFRONTSTRUCTMATRIX_H
#define TPZFRONTSTRUCTMATRIX_H

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"
/**
 * @brief Implements a structural matrix using the frontal method.
 * @ingroup structural frontal
 * @note Type parameter for TPZFrontStructMatrix frontal matrix. It can assume values TPZFrontSym and TPZFrontNonSym for symmetric and non symmetric matrices.*/
template<class TFront,
         class TVar = STATE,
         class TPar = TPZStructMatrixOR<TVar>> 
class TPZFrontStructMatrix : public TPZStructMatrixT<TVar>,
                             public TPar
{
	
protected:
	/** @brief This vector contains an ordered list */
	/** The elements must be asssembled in that order so the frontal works on its best performance */
	TPZVec<int> fElementOrder;
    int f_quiet{0};
	
	/**
	 * @brief Returns a vector containing all elements connected to a degree of freedom.
	 * @param numelconnected Vector containing the number of connections for every ith dof
	 */
	void GetNumElConnected(TPZVec <int> &numelconnected);

	/** @brief It is applied over fElementOrder putting it in the correct order. */
	void OrderElement();//TPZVec <int> &elorder);
	
	/** @brief Resequence the connects according to the element order*/
	void AdjustSequenceNumbering();
	
    /** @brief Used Decomposition method */
    DecomposeType fDecomposeType{ENoDecompose};

public:
    using TPZStructMatrixT<TVar>::TPZStructMatrixT;

    /// Set the decomposition type
    virtual void SetDecomposeType(DecomposeType dectype)
    {
        fDecomposeType = dectype;
    }
    
    
	/** @brief Returns a pointer to TPZBaseMatrix */
	TPZMatrix<TVar> * Create() override;
	
	/** @brief Clones a TPZFrontStructMatrix */
	TPZStructMatrix * Clone() override;
	
	/**
	 * @brief Assemble a stiffness matrix according to rhs
	 * @param stiffness Stiffness matrix to assembled
	 * @param rhs Vector containing loads
	 * @param guiInterface pointer to user interface
	 */ 	
	void AssembleNew(TPZMatrix<TVar> & stiffness
					 , TPZFMatrix<TVar> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/**
	 * @brief Assemble a stiffness matrix.
	 * @param stiffness Stiffness matrix to assembled
	 * @param rhs Vector containing loads
	 * @param guiInterface pointer to user interface
	 */ 	
	void Assemble(TPZBaseMatrix & stiffness
				  , TPZBaseMatrix & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface) override;
	
	/**
	 * @brief Computes element matrices.
	 * @param el Actual element being computed
	 * @param ek Formed element matrix
	 * @param ef Global element load matrix
	 * @param stiffness Global stiffness matrix
	 * @param rhs Global load matrix
	 */
	/** 
	 * Each computed element matrices would then be added to Stiffness matrix
	 */
	void AssembleElement(TPZCompEl *el, TPZElementMatrix & ek
						 , TPZElementMatrix & ef, TPZMatrix<TVar> & stiffness, TPZFMatrix<TVar> & rhs); 
	/**
	 * @brief Returns a pointer to TPZMatrix.
	 * @param rhs Load matrix
	 * @param guiInterface pointer to user interface
	 */
	/** 
	 * This is a mandatory function, it is neded by all StructMatrix. \n
	 * Except in frontal matrices, the returned matrix is not in its decomposed form.
	 */
	TPZMatrix<TVar> *CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
	
    void SetQuiet(int quiet);
    //@{
    int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
private:
	
    friend TPZPersistenceManager;
	
};

#endif //TPZFRONTSTRUCTMATRIX_H

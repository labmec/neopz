/**
 * @file
 * @brief Contains declaration of TPZElementMatrixT struct which associates an element matrix with the coeficients of its contribution in the global stiffness matrix.
 */

#ifndef TPZELEMENTMATRIXT_H
#define TPZELEMENTMATRIXT_H

#include "pzelmat.h"

template<class TVar>
struct TPZElementMatrixT : public TPZElementMatrix {
	void Reset(TPZCompMesh *mesh = NULL, MType type=Unknown)
	{
      TPZElementMatrix::Reset(mesh,type);
      fBlock.SetNBlocks(0);
      fConstrBlock.SetNBlocks(0);
      fMat.Resize(0,0);
      fConstrMat.Resize(0,0);
    }
	
	TPZElementMatrixT(TPZCompMesh *mesh, MType type) :
        TPZElementMatrix(mesh,type),
        fMat(0,0), fBlock(),
        fConstrMat(0,0), fConstrBlock()
    {
        fBlock.SetMatrix(&fMat);
        fConstrBlock.SetMatrix(&fConstrMat);
    }

    TPZElementMatrixT() :
        TPZElementMatrix(), fMat(0,0), fConstrMat(0,0),
        fBlock(&fMat), fConstrBlock(&fConstrMat)
    {
    }
    
    TPZElementMatrixT &operator=(const TPZElementMatrixT &copy);
	
    TPZElementMatrixT(const TPZElementMatrixT &copy);
    
	void Print(std::ostream &out) override;
	
	void SetMatrixSize(short NumBli, short NumBlj,
                       short BlSizei, short BlSizej) override;
	
	void SetMatrixMinSize(short NumBli, short NumBlj,
                          short BlMinSizei, short BlMinSizej) override;

    /** @brief permute the order of the connects */
    void PermuteGather(TPZVec<int64_t> &permute) override;
	/** @brief Apply the constraints applied to the nodes by transforming the tangent matrix and right hand side */
	void ApplyConstraints() override;
    
    /// Apply the constraint of the one shape restraints
    void ApplyOneShapeConstraints(int constraintindex) override;

    TVar &at(int64_t ibl, int64_t jbl, int idf, int jdf)
    {
        return fMat.at(fBlock.at(ibl,jbl,idf,jdf));
    }
    TVar &at(int64_t ibl, int idf)
    {
        return fMat(fBlock.Index(ibl,idf));
    }

    TPZFMatrix<TVar> & Matrix() override{
        return fMat;
    }

    TPZFMatrix<TVar> & ConstrMatrix() override{
        return fConstrMat;
    }

    TPZBlock & Block() override{
        return fBlock;
    }

    TPZBlock & ConstrBlock() override{
        return fConstrBlock;
    }
    /** @brief Pointer to a blocked matrix object*/
	TPZFNMatrix<1000, TVar> fMat;
    /** @brief Block structure associated with fMat*/
	TPZBlock fBlock;
	/** @brief Pointer to the constrained matrix object*/
	TPZFNMatrix<1000, TVar> fConstrMat;
    /** @brief Block structure associated with fConstrMat*/
	TPZBlock fConstrBlock;
};

extern template class TPZElementMatrixT<STATE>;
extern template class TPZElementMatrixT<CSTATE>;
#endif

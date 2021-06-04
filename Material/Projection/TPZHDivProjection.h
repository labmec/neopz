/**
 * @file TPZHDivProjection.h
 * @brief Header file for TPZHDivProjection.
 * It implements a projection of a given solution using the HDiv norm
 */

#ifndef TPZHDIVPROJECTION_H
#define TPZHDIVPROJECTION_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"

/**
 * @brief Implements a projection of a given solution on a HDiv approximation space.
 * This material uses a HDiv approximation space in a geometric domain of
 * dimension defined by the user (2 or 3D) and projects a given solution using
 * the HDiv norm.
 * The solution to be projected is set by the Forcing Function of the 
 * TPZMaterialT.
 * It is expected to be a TPZVec<TVar> of size 4, as the solution is always in
 * a 3D space and the last position is the divergence.
 * For the error analysis, the solution is set through the ExactSol of
 * TPZMatErrorSingleSpace.
 */
template<class TVar=STATE>
class  TPZHDivProjection :
    public TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>,TPZMatErrorSingleSpace<TVar>>
{
	using TBase = TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>,TPZMatErrorSingleSpace<TVar>>;
public:
    //! Default constructor
    TPZHDivProjection();
    /** @brief Creates the TPZHDivProjection material.
     * @param id material identifier
     * @param dim dimension
     */
    TPZHDivProjection(int id, int dim);
    
    std::string Name() const override { return "TPZHDivProjection"; }

    /** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, ESolution = 1 , EDivergence = 2};
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {fScale = scale;}
    
    [[nodiscard]] REAL ScaleFactor() const
    {return fScale;}

	int Dimension() const  override { return this->fDim; }
    
    /** @brief Sets problem dimension */
    virtual void SetDimension(int dim) { this->fDim = dim; }
	
    int NStateVariables() const override { return 1; }
    void Contribute(const TPZMaterialDataT<TVar> &data,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    /** @brief To create another material of the same type */
	TPZMaterial * NewMaterial() const override;
	
	/** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name) const override;
	
	int NSolutionVariables(int var) const override;

    void Solution(const TPZMaterialDataT<TVar> &data,
                  int var, TPZVec<TVar> &solOut) override;
    
    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override;

    void Errors(const TPZMaterialDataT<TVar> &data,
                TPZVec<REAL> &errors) override;

    [[nodiscard]] int IntegrationRuleOrder(const int elPMaxOrder) const override;
    virtual int ClassId() const override;
protected:
    /** @brief Problem dimension */
	int fDim;
    
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale{1};
};

extern template class TPZHDivProjection<STATE>;
extern template class TPZHDivProjection<CSTATE>;
#endif


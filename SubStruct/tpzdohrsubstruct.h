/**
 * @file
 * @brief Contains the TPZDohrSubstruct class which implements sub structure matrices using Dohrman algorithm.
 */

#ifndef TPZDOHRSUBSTRUCT_H
#define TPZDOHRSUBSTRUCT_H

#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "pzstepsolver.h"

/**
 * @ingroup substructure
 * @author Philippe Devloo
 * @warning Hardcoded definition of MAKEINTERNAL  ! 
 */
#ifndef MAKEINTERNAL
#define MAKEINTERNAL
#endif

/**
 * @ingroup substructure
 * @author Philippe Devloo
 * @brief Implements sub structure matrices using Dohrman algorithm. \ref substructure "Sub Structure"
 */
template<class TVar>
class TPZDohrSubstruct : public TPZSavable {
	// @TODO implement the interface to make the substruct class actually saveable
public:
    TPZDohrSubstruct();
	
    ~TPZDohrSubstruct();
    
    enum EWeightType {CorrectWeight, Uniform};
    
    int ClassId() const override  {
        return Hash("TPZDohrSubstruct") ^ ClassIdOrHash<TVar>() << 1;
    }
    
    static EWeightType fWeightType;
	
	/** @TODO Remove this - For PerfTests Only */
	void ReallocMatRed() {};

    /**
     * @brief ContributeResidual
     * @param u is the global solution
     * @param r is the residual of the solution "u"
     */
	/**
	 * It makes \f$ u(i)=R(i)*u \f$, \f$ r(i)=F(i)-K(i)*u(i) \f$ and \f$ R(i)*r(novo)=R(i)*r+r(i) \f$ where r(novo) is the new global residual
     * It puts r(novo) in r
	 */
    void ContributeResidual(TPZFMatrix<TVar> &u, TPZFMatrix<TVar> &r);
	
    void ContributeTestV1(TPZFMatrix<TVar> &testV1, int NumCoarse);
    /**
     * @brief It computes W(i)*R(i)*r and stores it in fLocalWeightedResidual;
     * @param r_global is the residual to the global solution
     */
    void LoadWeightedResidual(const TPZFMatrix<TVar> &r_global);
	
    /** @brief Adjust the residual to reflect a static condensation */
	/** The residual corresponding to the internal nodes will be zeroed */
    void AdjustResidual(TPZFMatrix<TVar> &r_global);
    /** @brief It computes the local contribution to r(c). */
	/** The method LoadWeightedResidual must be called before this one. */
    void Contribute_rc(TPZFMatrix<TVar> &rc);
    /** @brief It computes the local contribution to r(c). */
	/** The method LoadWeightedResidual must be called before this one. */
    void Contribute_rc_local(TPZFMatrix<TVar> &residual_local, TPZFMatrix<TVar> &rc_local) const;
    /** @brief It computes the local contribution to K(c) */
    void Contribute_Kc(TPZMatrix<TVar> &Kc, TPZVec<int> &coarseindex);
    /**
	 * @brief It computes the local contribution to v1.
	 * @param v1 
	 * @param invKc_rc is the product K(c)_inverted*r(c) 
	 */
	/** Of course r(c) must be computed, using Contribute_rc(), before calling this method */
    void Contribute_v1(TPZFMatrix<TVar> &v1, TPZFMatrix<TVar> &invKc_rc);
    /**
	 * @brief It computes the local contribution to v1.
	 * @param v1_local 
	 * @param invKc_rc is the product K(c)_inverted*r(c)
	 */
	/** Of course r(c) must be computed, using Contribute_rc(), before calling this method */
    void Contribute_v1_local(TPZFMatrix<TVar> &v1_local, TPZFMatrix<TVar> &invKc_rc) const;
    /** @brief It computes the local contribution to v2. */
    void Contribute_v2(TPZFMatrix<TVar> &v2);
    /** @brief It computes the local contribution to v2. */
    void Contribute_v2_local(TPZFMatrix<TVar> &residual_local, TPZFMatrix<TVar> &v2_local);
    /**
	 * @brief It computes the local contribution to v(3)
	 * @param v3 
	 * @param r is the global residual
	 * @param v1Plusv2 is the sum "v1 + v2"
	 */
    void Contribute_v3(TPZFMatrix<TVar> &v3, const TPZFMatrix<TVar> &r, TPZFMatrix<TVar> &v1Plusv2) const;
    /**
	 * @brief It computes the local contribution to v(3)
	 * @param v3 
	 * @param v1Plusv2 is the sum "v1 + v2"
	 */
    void Contribute_v3_local(TPZFMatrix<TVar> &v3, TPZFMatrix<TVar> &v1Plusv2) const;
	
    void Print(std::ostream &out) const;
    /**
     * @brief Contribute to the global matrix vector multiplication. It's needed to the MultAdd method of TPZDohrMatrix.
     */
    void ContributeKU(const TVar alpha, const TPZFMatrix<TVar> &uglobal, TPZFMatrix<TVar> &z) const;
	
	/** @brief Compute the multiplication of the local stiffness matrix with the vector u */
	void ContributeKULocal(const TVar alpha, const TPZFMatrix<TVar> &u, TPZFMatrix<TVar> &z) const;
    /**
	 * @brief Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
	 * @param StiffnessDiag is the diagonal of the stiffness matrix
	 */
    void ContributeGlobalDiagonal(TPZFMatrix<TVar> &StiffnessDiag);
    /**
	 * @brief Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
	 * @param StiffnessDiag is the diagonal of the stiffness matrix
	 */
    void ContributeDiagonalLocal(TPZFMatrix<TVar> &StiffnessDiag);
    /**
	 * @brief Computes the weight matrix.
	 * @param StiffnessDiag is the diagonal of the global matrix (global numbering)
	 */
    void ComputeWeights(TPZFMatrix<TVar> &StiffnessDiag);
    /**
	 * @brief Computes the weight matrix.
	 * @param StiffnessDiagLocal is the diagonal of the global matrix (local numbering)
	 */
    void ComputeWeightsLocal(TPZFMatrix<TVar> &StiffnessDiagLocal);
    /**
	 * @brief Initializes the substructure.
	 */
    void Initialize();
    
    //Daqui pra baixo são funcões cujos nomes foram dados pelo Phil
    /**
     * @brief Assembles the contribution to the coarse residual
     */
    void GetCoarseResidual(TPZFMatrix<TVar> &rc, TPZVec<int> &indices);
    
    /**
     * @brief Assembles the coarse dof stiffness matrix
     */
    void GetCoarseStiffness(TPZMatrix<TVar> &stiff, TPZVec<int> &indices);
public:
    /**
	 * @brief It prepares the datas for solving systems for phi and zi
	 */
    void PrepareSystems();
    /** @brief Solves the system for Phi and for v2 */
	/** It stores the results in fPhiC and fzi */
    void SolveSystemPhi();
    /** @brief Solves the system for zi */
	/** It stores the results in fzi */
    void SolveSystemZi();
    /** @brief Computes K(ci) and stores it in fKCi */
    void ComputeCoarseStiffness();
	/** @brief Add the internal solution to the final result */
	void AddInternalSolution(TPZFMatrix<TVar> &sol);
	
    /** @brief Variables needed to solve the systems for phi and zi */
	TPZFMatrix<TVar> fC_star;
    TPZFMatrix<TVar> fKeC_star; //K_star_inv*C_star_trans
    TPZStepSolver<TVar> finv;
    TPZFMatrix<TVar> fNullPivots;

	/** @brief Number of equations of the substructure */
	int fNEquations; 
	
	int NumEquations() const
	{
		return fNEquations;
	}
    /**
	 * @brief 
	 */
    TPZVec<int> fCoarseNodes;
	
    /** @brief Constraint definition associated with this substructure */
	/** Input data */
    TPZFMatrix<TVar> fC; //C(i)
    /** @brief Vectors associated with each constraint */
    TPZFMatrix<TVar> fPhiC; //Phí(i)
	/** @brief Phi * W matrix and condensed to the equations which are part of the global system */
	TPZFMatrix<TVar> fPhiC_Weighted_Condensed;

    /** @brief Global index associated with each coarse degree of freedom */
	/** Input data */
    TPZVec<int> fCoarseIndex; //R(ci)
    /** @brief Weights associated with each variable/equation */ 
	/** Computed */
    TPZManVector<TVar> fWeights; //W(i) = diagonal(fWeights)
    /** @brief Stiffness matrix associated with the substructure. */
	/** Input data */
    TPZAutoPointer<TPZMatrix<TVar> > fStiffness; //K(i)
    /** @brief Stiffness matrix associated with the constraints */
	/** Computed */
    TPZFMatrix<TVar> fKCi; //K(ci)
    /** @brief Global vector indices of the equations/degrees of freedom */
	/** Input data */
    TPZVec<int> fGlobalIndex; //R(i)
	
	
	/** @brief Indices of the equations in the global system corresponding to the local equations */
	TPZVec<std::pair<int,int> > fGlobalEqs;
    /** @brief Inverted (LU or Cholesky) stiffness matrix */
	/** Computed */
    TPZStepSolver<TVar> fInvertedStiffness; //Ke, or K_star
    /** @brief Inverted (LU or Cholesky or LDLt) stiffness matrix for the internal degrees of freedom */
	/** Computed */
    mutable TPZStepSolver<TVar> fInvertedInternalStiffness; //R(Ii)*K(i)*R(Ii)invertida
    /** @brief Inverted restraint matrix */
	/** Computed */
    TPZStepSolver<TVar> fKCInvert; //K(c) invertida
    /** @brief Internal nodes */
	/** Data */
    TPZVec<int> fInternalEqs; //R(Ii)
	/** @brief Equations corresponding to boundary of the substructure */
	TPZVec<int> fBoundaryEqs;
    /** @brief Local load vector */
	/** Data */
    TPZFMatrix<TVar> fLocalLoad; //F(i) - Actually it's a vector
    /** @brief Local weighted residual - \f$ W(i)*R(i)*r \f$ */
    TPZFMatrix<TVar> fLocalWeightedResidual; //W(i)*R(i)*r - Actually it's a vector
    /** @brief Needed to compute v2. Is the solution of a system */
    TPZFMatrix<TVar> fzi; //z(i)
	
	/** @brief Solution vector which needs to be added to the converged system solution */
	/** This variable is initialized and set in the AdjustResidual method */
	TPZFMatrix<TVar> fAdjustSolution;
	
};

#endif

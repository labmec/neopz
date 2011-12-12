/**
 * @file
 * @brief Contains the TPZDohrSubstruct class which implements sub structure matrices using Dohrman algorithm.
 */
/***************************************************************************
 *   Copyright (C) 2006 by Philippe Devloo   *
 *   phil@fec.unicamp.br   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TPZDOHRSUBSTRUCT_H
#define TPZDOHRSUBSTRUCT_H

#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "pzstepsolver.h"

/**
 @ingroup substructure
 @author Philippe Devloo
 */

// #warning Hardcoded definition of MAKEINTERNAL !!!!
/** @warning Hardcoded definition of MAKEINTERNAL  ! */
#define MAKEINTERNAL

/**
 * @ingroup substructure
 * @author Philippe Devloo
 * @brief Implements sub structure matrices using Dohrman algorithm. \ref substructure "Sub Structure"
 */
class TPZDohrSubstruct : public TPZSaveable {
	// @TODO implement the interface to make the substruct class actually saveable
public:
    TPZDohrSubstruct();
	
    ~TPZDohrSubstruct();
    
    enum EWeightType {CorrectWeight, Uniform};
    
    static EWeightType fWeightType;
	
    /**
     * @brief ContributeResidual
     * @param u is the global solution
     * @param r is the residual of the solution "u"
     */
	/**
	 * It makes \f$ u(i)=R(i)*u \f$, \f$ r(i)=F(i)-K(i)*u(i) \f$ and \f$ R(i)*r(novo)=R(i)*r+r(i) \f$ where r(novo) is the new global residual
     * It puts r(novo) in r
	 */
    void ContributeResidual(TPZFMatrix &u, TPZFMatrix &r);
	
    void ContributeTestV1(TPZFMatrix &testV1, int NumCoarse);
    /**
     * @brief It computes W(i)*R(i)*r and stores it in fLocalWeightedResidual;
     * @param r_global is the residual to the global solution
     */
    void LoadWeightedResidual(const TPZFMatrix &r_global);
	
    /** @brief Adjust the residual to reflect a static condensation */
	/** The residual corresponding to the internal nodes will be zeroed */
    void AdjustResidual(TPZFMatrix &r_global);
    /** @brief It computes the local contribution to r(c). */
	/** The method LoadWeightedResidual must be called before this one. */
    void Contribute_rc(TPZFMatrix &rc);
    /** @brief It computes the local contribution to r(c). */
	/** The method LoadWeightedResidual must be called before this one. */
    void Contribute_rc_local(TPZFMatrix &residual_local, TPZFMatrix &rc_local);
    /** @brief It computes the local contribution to K(c) */
    void Contribute_Kc(TPZMatrix &Kc, TPZVec<int> &coarseindex);
    /**
	 * @brief It computes the local contribution to v1.
	 * @param v1 
	 * @param invKc_rc is the product K(c)_inverted*r(c) 
	 */
	/** Of course r(c) must be computed, using Contribute_rc(), before calling this method */
    void Contribute_v1(TPZFMatrix &v1, TPZFMatrix &invKc_rc);
    /**
	 * @brief It computes the local contribution to v1.
	 * @param v1_local 
	 * @param invKc_rc is the product K(c)_inverted*r(c)
	 */
	/** Of course r(c) must be computed, using Contribute_rc(), before calling this method */
    void Contribute_v1_local(TPZFMatrix &v1_local, TPZFMatrix &invKc_rc);
    /** @brief It computes the local contribution to v2. */
    void Contribute_v2(TPZFMatrix &v2);
    /** @brief It computes the local contribution to v2. */
    void Contribute_v2_local(TPZFMatrix &residual_local, TPZFMatrix &v2_local
							 );
    /**
	 * @brief It computes the local contribution to v(3)
	 * @param v3 
	 * @param r is the global residual
	 * @param v1Plusv2 is the sum "v1 + v2"
	 */
    void Contribute_v3(TPZFMatrix &v3, const TPZFMatrix &r, TPZFMatrix &v1Plusv2) const;
    /**
	 * @brief It computes the local contribution to v(3)
	 * @param v3 
	 * @param v1Plusv2 is the sum "v1 + v2"
	 */
    void Contribute_v3_local(TPZFMatrix &v3, TPZFMatrix &v1Plusv2) const;
	
    void Print(std::ostream &out) const;
    /**
     * @brief Contribute to the global matrix vector multiplication. It's needed to the MultAdd method of TPZDohrMatrix.
     */
    void ContributeKU(const REAL alpha, const TPZFMatrix &uglobal, TPZFMatrix &z) const;
	
	/** @brief Compute the multiplication of the local stiffness matrix with the vector u */
	void ContributeKULocal(const REAL alpha, const TPZFMatrix &u, TPZFMatrix &z) const;
    /**
	 * @brief Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
	 * @param StiffnessDiag is the diagonal of the stiffness matrix
	 */
    void ContributeGlobalDiagonal(TPZFMatrix &StiffnessDiag);
    /**
	 * @brief Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
	 * @param StiffnessDiag is the diagonal of the stiffness matrix
	 */
    void ContributeDiagonalLocal(TPZFMatrix &StiffnessDiag);
    /**
	 * @brief Computes the weight matrix.
	 * @param StiffnessDiag is the diagonal of the global matrix (global numbering)
	 */
    void ComputeWeights(TPZFMatrix &StiffnessDiag);
    /**
	 * @brief Computes the weight matrix.
	 * @param StiffnessDiagLocal is the diagonal of the global matrix (local numbering)
	 */
    void ComputeWeightsLocal(TPZFMatrix &StiffnessDiagLocal);
    /**
	 * @brief Initializes the substructure.
	 */
    void Initialize();
    
    //Daqui pra baixo são funcões cujos nomes foram dados pelo Phil
    /**
     * @brief Assembles the contribution to the coarse residual
     */
    void GetCoarseResidual(TPZFMatrix &rc, TPZVec<int> &indices);
    
    /**
     * @brief Assembles the coarse dof stiffness matrix
     */
    void GetCoarseStiffness(TPZMatrix &stiff, TPZVec<int> &indices);
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
	void AddInternalSolution(TPZFMatrix &sol);
	
    /** @brief Variables needed to solve the systems for phi and zi */
	TPZFMatrix fC_star;
    TPZFMatrix fKeC_star; //K_star_inv*C_star_trans
    TPZStepSolver finv;
    TPZFMatrix fNullPivots;

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
    /**
     * @brief Constraint definition associated with this substructure
     */
	/** Input data */
    TPZFMatrix fC; //C(i)
    /** @brief Vectors associated with each constraint */
    TPZFMatrix fPhiC; //Phí(i)
	/** @brief Phi * W matrix and condensed to the equations which are part of the global system */
	TPZFMatrix fPhiC_Weighted_Condensed;
    /**
     * @brief Global index associated with each coarse degree of freedom
	 */
	/** Input data */
    TPZVec<int> fCoarseIndex; //R(ci)
    /**
     * @brief Weights associated with each variable/equation
     */ 
	/** Computed */
    TPZManVector<REAL> fWeights; //W(i) = diagonal(fWeights)
    /**
     * @brief Stiffness matrix associated with the substructure.
	 */
	/** Input data */
    TPZAutoPointer<TPZMatrix> fStiffness; //K(i)
    /**
     * @brief Stiffness matrix associated with the constraints
	 */
	/** Computed */
    TPZFMatrix fKCi; //K(ci)
    /**
     * @brief Global vector indices of the equations/degrees of freedom
	 */
	/** Input data */
    TPZVec<int> fGlobalIndex; //R(i)
	
	
	/** @brief Indices of the equations in the global system corresponding to the local equations */
	TPZVec<std::pair<int,int> > fGlobalEqs;
    /**
     * @brief Inverted (LU or Cholesky) stiffness matrix
	 */
	/** Computed */
    TPZStepSolver fInvertedStiffness; //Ke, or K_star
    /**
     * @brief Inverted (LU or Cholesky or LDLt) stiffness matrix for the internal degrees of freedom
     */
	/** Computed */
    mutable TPZStepSolver fInvertedInternalStiffness; //R(Ii)*K(i)*R(Ii)invertida
    /**
     * @brief Inverted restraint matrix
     */
	/** Computed */
    TPZStepSolver fKCInvert; //K(c) invertida
    /**
     * @brief Internal nodes
	 */
	/** Data */
    TPZVec<int> fInternalEqs; //R(Ii)
	/**
	 * @brief Equations corresponding to boundary of the substructure
	 */
	TPZVec<int> fBoundaryEqs;
    /**
     * @brief Local load vector
     */
	/** Data */
    TPZFMatrix fLocalLoad; //F(i) - Actually it's a vector
    /**
	 * @brief Local weighted residual - W(i)*R(i)*r
	 */
    TPZFMatrix fLocalWeightedResidual; //W(i)*R(i)*r - Actually it's a vector
    /**
	 * The columns of fEigenVectors are the eigenvectors of fStiffness, or K(i), associated with the null eigenvalue
	 */
	//    TPZFMatrix fEigenVectors;
    /**
	 * @brief Needed to compute v2. Is the solution of a system
	 */
    TPZFMatrix fzi; //z(i)
	
	/**
	 * @brief Solution vector which needs to be added to the converged system solution
	 */
	/** This variable is initialized and set in the AdjustResidual method */
	TPZFMatrix fAdjustSolution;
	
};

#endif

/**
 * @file
 * @brief Contains the TPZDohrSubstructCondense class which condenses matrix divided in sub structures.
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
#ifndef TPZDOHRSUBSTRUCTCONDENSE_H
#define TPZDOHRSUBSTRUCTCONDENSE_H

#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "pzstepsolver.h"
#include "pzysmp.h"
#include "pzmatred.h"
#include "tpzverysparsematrix.h"

/**
 * @brief To condense matrix divided in sub structures. \ref substructure "Sub Structure"
 * @author Philippe Devloo
 * @ingroup substructure
 */
class TPZDohrSubstructCondense : public TPZSaveable
{
public:
	TPZDohrSubstructCondense();
	
	~TPZDohrSubstructCondense();
	
	enum EWeightType {CorrectWeight, Uniform};
	
	enum ENumbering {Submesh, InternalFirst, ExternalFirst};
	
	static EWeightType fWeightType;
	
	/** @brief Return the number of equations which contribute to the global system */
	int NumEquations()
	{
		return fNumExternalEquations;
	}
	
	/**
	 * @brief Adjust the residual to reflect a static condensation
	 *
	 * The residual corresponding to the internal nodes will be zeroed
	 */
	void AdjustResidual(TPZFMatrix<REAL> &r_global);
	/**
	 * @brief It computes the local contribution to r(c).
	 *
	 * The method LoadWeightedResidual must be called before this one.
	 */
	void Contribute_rc_local(TPZFMatrix<REAL> &residual_local, TPZFMatrix<REAL> &rc_local);
	/**
	 * @brief It computes the local contribution to K(c)
	 */
	void Contribute_Kc(TPZMatrix<REAL> &Kc, TPZVec<int> &coarseindex);
	/**
	 * @brief It computes the local contribution to v1.
	 * @param v1_local 
	 * @param invKc_rc is the product K(c)_inverted*r(c)
	 * Of course r(c) must be computed, using Contribute_rc(), before calling this method
	 */
	void Contribute_v1_local(TPZFMatrix<REAL> &v1_local, TPZFMatrix<REAL> &invKc_rc);
	/**
	 * @brief It computes the local contribution to v2.
	 */
	void Contribute_v2_local(TPZFMatrix<REAL> &residual_local, TPZFMatrix<REAL> &v2_local
							 );
	/**
	 * @brief It computes the local contribution to v(3)
	 * @param v3
	 * @param v1Plusv2 is the sum "v1 + v2"
	 */
	void Contribute_v3_local(TPZFMatrix<REAL> &v3, TPZFMatrix<REAL> &v1Plusv2) const;
	/**
	 * @brief 
	 */
	void Print(std::ostream &out) const;
	
	/** @brief compute the multiplication of the local stiffness matrix with the vector u */
	void ContributeKULocal(const REAL alpha, const TPZFMatrix<REAL> &u, TPZFMatrix<REAL> &z) const;
	/**
	 * @brief Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
	 * @param StiffnessDiag is the diagonal of the stiffness matrix
	 */
	void ContributeDiagonalLocal(TPZFMatrix<REAL> &StiffnessDiag);
	/**
	 * @brief Computes the weight matrix.
	 * @param StiffnessDiagLocal is the diagonal of the global matrix (local numbering)
	 */
	void ComputeWeightsLocal(TPZFMatrix<REAL> &StiffnessDiagLocal);
	/**
	 * @brief Computes the condensed right hand side for the substructure
	 * @param rhs the right hand side ordered external first
	 */
	void ContributeRhs(TPZFMatrix<REAL> &rhs);
	/**
	 * @brief Computes the global solution based on the interface solution
	 */
	void UGlobal(TPZFMatrix<REAL> &Uext, TPZFMatrix<REAL> &UGlobal);
	/**
	 * @brief Initializes the substructure.
	 */
	void Initialize();
	
	//Daqui pra baixo são funcões cujos nomes foram dados pelo Phil
	/**
	 * @brief Assemble the contribution to the coarse residual
	 */
	void GetCoarseResidual(TPZFMatrix<REAL> &rc, TPZVec<int> &indices);
	
	/**
	 * @brief Assemble the coarse dof stiffness matrix
	 */
	void GetCoarseStiffness(TPZMatrix<REAL> &stiff, TPZVec<int> &indices);
	
	/**
	 * @brief Apply a scatter permutation to the input vector using a scatter permutation
	 * output[permute[i]] = input[i-first], first <= i < last
	 */
	/** This method does not resize the elements */
	static void PermuteScatter(const TPZVec<int> &permute, const TPZFMatrix<REAL> &input, TPZFMatrix<REAL> &output, int first, int last);
	/**
	 * @brief Apply a gather permutation to the input vector using a scatter permutation
	 * output[i-first] = input[permute[i]], first <= i < last
	 */
	/** This method does not resize the elements */
	static void PermuteGather(const TPZVec<int> &permute, const TPZFMatrix<REAL> &input, TPZFMatrix<REAL> &output, int first, int last);
	
public:
	/**
	 * @brief It prepares the datas for solving systems for phi and zi
	 */
	void PrepareSystems();
	/**
	 * @brief Solves the system for Phi and for v2 \n
	 * It stores the results in fPhiC and fzi
	 */
	void SolveSystemPhi();
	/**
	 * @brief Add the internal solution to the final result
	 */
	void AddInternalSolution(TPZFMatrix<REAL> &sol);
	/**
	 * Variables needed to solve the systems for phi and zi
	 */
	//	TPZFMatrix<REAL> fC_star;
	//    TPZFMatrix<REAL> fKeC_star; //K_star_inv*C_star_trans
	//    TPZStepSolver finv;
	//    TPZFMatrix<REAL> fNullPivots;
	
	/** @brief Matrix problem which solves the phi and zi problems */
	TPZAutoPointer<TPZMatRed<REAL, TPZFMatrix<REAL> > > fMatRedComplete;
	/** @brief Number of equations of the substructure */
	int fNEquations; 
	
	/** @brief Number of internal equations of the substructure */
	int fNumInternalEquations;
	
	/** @brief Number of equations which connect to the global structure */
	int fNumExternalEquations;
	/**
	 * @brief 
	 */
	TPZVec<int> fCoarseNodes;
	/**
	 * Constraint definition associated with this substructure
	 *
	 * Input data
	 */
	//    TPZFMatrix<REAL> fC; //C(i)
	/**
	 * @brief Vectors associated with each constraint
	 *
	 * Computed
	 */
	TPZFMatrix<REAL> fPhiC; //Phí(i)
	/** @brief Phi * W matrix and condensed to the equations which are part of the global system */
	TPZFMatrix<REAL> fPhiC_Weighted_Condensed;
	/**
	 * Global index associated with each coarse degree of freedom
	 * Input data
	 */
	//    TPZVec<int> fCoarseIndex; //R(ci)
	/**
	 * @brief Weights associated with each variable/equation
	 * 
	 * Computed
	 */
	TPZManVector<REAL> fWeights; //W(i) = diagonal(fWeights)
	/**
	 * Stiffness matrix associated with the substructure
	 * Input data
	 */
	//    TPZAutoPointer<TPZMatrix> fStiffness; //K(i)
	/**
	 * @brief Stiffness matrix associated with the constraints
	 *
	 * Computed
	 */
	TPZFMatrix<REAL> fKCi; //K(ci)
	/**
	 * Global vector indices of the equations/degrees of freedom
	 * Input data
	 */
	//TPZVec<int> fGlobalIndex; //R(i)
	
	/** @brief Permutation vectors
	 *
	 * This map holds all permutation vectors for the equations
	 * The format of a permutation is that permute[ind] corresponds to the destination index of the element index ind of the original vector
	 * solout[permute[ind]] = solin[ind], ind = 0, n, 1
	 */
	std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> > fPermutationsScatter;
	
	TPZVec<int> &GatherVec(ENumbering origin, ENumbering destination)
	{
		return fPermutationsScatter[std::pair<ENumbering,ENumbering>(destination,origin)];
	}
	
	TPZVec<int> &ScatterVec(ENumbering origin, ENumbering destination)
	{
		return fPermutationsScatter[std::pair<ENumbering,ENumbering>(origin,destination)];
	}
	
	const TPZVec<int> &GatherVec(ENumbering origin, ENumbering destination) const;
	
	const TPZVec<int> &ScatterVec(ENumbering origin, ENumbering destination) const;
	
	/**
	 * @brief Inverted (LU or Cholesky or LDLt) stiffness matrix for the internal degrees of freedom
	 * Calculado
	 */
	mutable TPZAutoPointer<TPZMatRed<REAL, TPZFMatrix<REAL> > > fMatRed; //R(Ii)*K(i)*R(Ii)invertida
	/**
	 * @brief Local load vector
	 * 
	 * Data
	 */
	TPZFMatrix<REAL> fLocalLoad; //F(i) - Actually it's a vector
	/**
	 * @brief Local weighted residual - W(i)*R(i)*r
	 */
	TPZFMatrix<REAL> fLocalWeightedResidual; //W(i)*R(i)*r - Actually it's a vector
	
	/**
	 * @brief Solution vector which needs to be added to the converged system solution
	 * 
	 * This variable is initialized and set in the AdjustResidual method
	 */
	TPZFMatrix<REAL> fAdjustSolution;
	
};

#endif

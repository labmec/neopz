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
 @author Philippe Devloo
 */
class TPZDohrSubstructCondense
	{
	public:
		TPZDohrSubstructCondense();
		
		~TPZDohrSubstructCondense();
		
		enum EWeightType {CorrectWeight, Uniform};
		
		enum ENumbering {Submesh, InternalFirst, ExternalFirst};
		
		static EWeightType fWeightType;
		
		/// return the number of equations which contribute to the global system
		int NumEquations()
		{
			return fNumExternalEquations;
		}
		
		/**
		 * Adjust the residual to reflect a static condensation
		 * The residual corresponding to the internal nodes will be zeroed
		 */
		void AdjustResidual(TPZFMatrix &r_global);
		/**
		 * It computes the local contribution to r(c).
		 * The method LoadWeightedResidual must be called before this one.
		 */
		void Contribute_rc_local(TPZFMatrix &residual_local, TPZFMatrix &rc_local);
		/**
		 * It computes the local contribution to K(c)
		 */
		void Contribute_Kc(TPZMatrix &Kc, TPZVec<int> &coarseindex);
		/**
		 * It computes the local contribution to v1.
		 * @param invKc_rc is the product K(c)_inverted*r(c)
		 * Of course r(c) must be computed, using Contribute_rc(), before calling this method
		 */
		void Contribute_v1_local(TPZFMatrix &v1_local, TPZFMatrix &invKc_rc);
		/**
		 * It computes the local contribution to v2.
		 */
		void Contribute_v2_local(TPZFMatrix &residual_local, TPZFMatrix &v2_local
								 );
		/**
		 * It computes the local contribution to v(3)
		 * @param r is the global residual
		 * @param v1Plusv2 is the sum "v1 + v2"
		 */
		void Contribute_v3_local(TPZFMatrix &v3, TPZFMatrix &v1Plusv2) const;
		/**
		 *
		 */
		void Print(std::ostream &out) const;
		
		/// compute the multiplication of the local stiffness matrix with the vector u
		void ContributeKULocal(const REAL alpha, const TPZFMatrix &u, TPZFMatrix &z) const;
		/**
		 * Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
		 * @param StiffnessDiag is the diagonal of the stiffness matrix
		 */
		void ContributeDiagonalLocal(TPZFMatrix &StiffnessDiag);
		/**
		 * Computes the weight matrix.
		 * @param StiffnessDiag is the diagonal of the global matrix (local numbering)
		 */
		void ComputeWeightsLocal(TPZFMatrix &StiffnessDiagLocal);
		/**
		 * Initializes the substructure.
		 */
		void Initialize();
		
		//Daqui pra baixo são funcões cujos nomes foram dados pelo Phil
		/**
		 * Assemble the contribution to the coarse residual
		 */
		void GetCoarseResidual(TPZFMatrix &rc, TPZVec<int> &indices);
		
		/**
		 * Assemble the coarse dof stiffness matrix
		 */
		void GetCoarseStiffness(TPZMatrix &stiff, TPZVec<int> &indices);
		
		/**
		 * Apply a scatter permutation to the input vector using a scatter permutation
		 * output[permute[i]] = input[i-first], first <= i < last
		 * this method does not resize the elements
		 */
		static void PermuteScatter(const TPZVec<int> &permute, const TPZFMatrix &input, TPZFMatrix &output, int first, int last);
		/**
		 * Apply a gather permutation to the input vector using a scatter permutation
		 * output[i-first] = input[permute[i]], first <= i < last
		 * this method does not resize the elements
		 */
		static void PermuteGather(const TPZVec<int> &permute, const TPZFMatrix &input, TPZFMatrix &output, int first, int last);
		
	public:
		/**
		 * It prepares the datas for solving systems for phi and zi
		 */
		void PrepareSystems();
		/**
		 * Solves the system for Phi and for v2
		 * It stores the results in fPhiC and fzi
		 */
		void SolveSystemPhi();
		/**
		 * Add the internal solution to the final result
		 */
		void AddInternalSolution(TPZFMatrix &sol);
		/**
		 * Variables needed to solve the systems for phi and zi
		 */
		//	TPZFMatrix fC_star;
		//    TPZFMatrix fKeC_star; //K_star_inv*C_star_trans
		//    TPZStepSolver finv;
		//    TPZFMatrix fNullPivots;
		
		/// Matrix problem which solves the phi and zi problems
		TPZAutoPointer<TPZMatRed<TPZFMatrix> > fMatRedComplete;
		/// Number of equations of the substructure
		int fNEquations; 
		
		/// Number of internal equations of the substructure
		int fNumInternalEquations;
		
		/// Number of equations which connect to the global structure
		int fNumExternalEquations;
		/**
		 *
		 */
		TPZVec<int> fCoarseNodes;
		/**
		 * Constraint definition associated with this substructure
		 * Dado de entrada
		 */
		//    TPZFMatrix fC; //C(i)
		/**
		 * Vectors associated with each constraint
		 * Calculado
		 */
		TPZFMatrix fPhiC; //Phí(i)
		/// Phi * W matrix and condensed to the equations which are part of the global system
		TPZFMatrix fPhiC_Weighted_Condensed;
		/**
		 * Global index associated with each coarse degree of freedom
		 * Dado de entrada
		 */
		//    TPZVec<int> fCoarseIndex; //R(ci)
		/**
		 * Weights associated with each variable/equation
		 * Calculado
		 */
		TPZManVector<REAL> fWeights; //W(i) = diagonal(fWeights)
		/**
		 * Stiffness matrix associated with the substructure
		 * Dado de entrada
		 */
		//    TPZAutoPointer<TPZMatrix> fStiffness; //K(i)
		/**
		 * Stiffness matrix associated with the constraints
		 * Calculado
		 */
		TPZFMatrix fKCi; //K(ci)
		/**
		 * Global vector indices of the equations/degrees of freedom
		 * Dado de entrada
		 */
		//TPZVec<int> fGlobalIndex; //R(i)
		
		/// Permutation vectors
		/**
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
		
		/// Indices of the equations in the global system corresponding to the local equations
		//	TPZVec<std::pair<int,int> > fGlobalEqs;
		/**
		 * Inverted (LU or Cholesky) stiffness matrix
		 * Calculado
		 */
		//    TPZStepSolver fInvertedStiffness; //Ke, or K_star
		/**
		 * Inverted (LU or Cholesky or LDLt) stiffness matrix for the internal degrees of freedom
		 * Calculado
		 */
		mutable TPZAutoPointer<TPZMatRed<TPZVerySparseMatrix> > fMatRed; //R(Ii)*K(i)*R(Ii)invertida
		/**
		 * Sparse matrices corresponding to the matrices after the internal equations were reordered as first
		 */
		//	TPZFYsmpMatrix fK01,fK10,fK11;
		/**
		 * Inverted restraint matrix
		 * Calculado
		 */
		//    TPZStepSolver fKCInvert; //K(c) invertida
		/**
		 * internal nodes
		 * Dado
		 */
		//    TPZVec<int> fInternalEqs; //R(Ii)
		/**
		 * equations corresponding to boundary of the substructure
		 */
		//	TPZVec<int> fBoundaryEqs;
		/**
		 * local load vector
		 * Dado
		 */
		TPZFMatrix fLocalLoad; //F(i) - Actually it's a vector
		/**
		 * Local weighted residual - W(i)*R(i)*r
		 */
		TPZFMatrix fLocalWeightedResidual; //W(i)*R(i)*r - Actually it's a vector
		/**
		 * The columns of fEigenVectors are the eigenvectors of fStiffness, or K(i), associated with the null eigenvalue
		 */
		//    TPZFMatrix fEigenVectors;
		/**
		 * Needed to compute v2. Is the solution of a system
		 */
		//    TPZFMatrix fzi; //z(i)
		
		/**
		 * Solution vector which needs to be added to the converged system solution
		 * this variable is initialized and set in the AdjustResidual method
		 */
		TPZFMatrix fAdjustSolution;
		
	};

#endif

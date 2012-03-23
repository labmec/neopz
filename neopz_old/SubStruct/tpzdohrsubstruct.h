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
@author Philippe Devloo
*/


#warning Hardcoded definition of MAKEINTERNAL !!!!
#define MAKEINTERNAL




class TPZDohrSubstruct{
public:
    TPZDohrSubstruct();

    ~TPZDohrSubstruct();
    
    enum EWeightType {CorrectWeight, Uniform};
    
    static EWeightType fWeightType;

    /**
     * ContributeResidual
     * It makes u(i)=R(i)*u, r(i)=F(i)-K(i)*u(i) and R(i)*r(novo)=R(i)*r+r(i) where r(novo) is the new global residual
     * It puts r(novo) in r
     * @param u is the global solution
     * @param r is the residual of the solution "u"
     */
    void ContributeResidual(TPZFMatrix &u, TPZFMatrix &r);
    /**
    *
    */
    void ContributeTestV1(TPZFMatrix &testV1, int NumCoarse);
    /**
     * It computes W(i)*R(i)*r and stores it in fLocalWeightedResidual;
     * @param r_global is the residual to the global solution
     */
    void LoadWeightedResidual(const TPZFMatrix &r_global);

    /**
     * Adjust the residual to reflect a static condensation
     * The residual corresponding to the internal nodes will be zeroed
     */
    void AdjustResidual(TPZFMatrix &r_global);
    /**
    * It computes the local contribution to r(c).
    * The method LoadWeightedResidual must be called before this one.
    */
    void Contribute_rc(TPZFMatrix &rc);
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
    void Contribute_v1(TPZFMatrix &v1, TPZFMatrix &invKc_rc);
    /**
	 * It computes the local contribution to v1.
	 * @param invKc_rc is the product K(c)_inverted*r(c)
	 * Of course r(c) must be computed, using Contribute_rc(), before calling this method
	 */
    void Contribute_v1_local(TPZFMatrix &v1_local, TPZFMatrix &invKc_rc);
    /**
    * It computes the local contribution to v2.
    */
    void Contribute_v2(TPZFMatrix &v2);
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
    void Contribute_v3(TPZFMatrix &v3, const TPZFMatrix &r, TPZFMatrix &v1Plusv2) const;
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
    /**
     * Contribute to the global matrix vector multiplication. It's needed to the MultAdd method of TPZDohrMatrix.
     */
    void ContributeKU(const REAL alpha, const TPZFMatrix &uglobal, TPZFMatrix &z) const;
	
	/// compute the multiplication of the local stiffness matrix with the vector u
	void ContributeKULocal(const REAL alpha, const TPZFMatrix &u, TPZFMatrix &z) const;
    /**
    * Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
    * @param StiffnessDiag is the diagonal of the stiffness matrix
    */
    void ContributeGlobalDiagonal(TPZFMatrix &StiffnessDiag);
    /**
	 * Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
	 * @param StiffnessDiag is the diagonal of the stiffness matrix
	 */
    void ContributeDiagonalLocal(TPZFMatrix &StiffnessDiag);
    /**
    * Computes the weight matrix.
    * @param StiffnessDiag is the diagonal of the global matrix (global numbering)
    */
    void ComputeWeights(TPZFMatrix &StiffnessDiag);
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
    * Solves the system for zi
    * It stores the results in fzi
    */
    void SolveSystemZi();
    /**
    * Computes K(ci) and stores it in fKCi
    */
    void ComputeCoarseStiffness();
	/**
	 * Add the internal solution to the final result
	 */
	void AddInternalSolution(TPZFMatrix &sol);
    /**
     * Variables needed to solve the systems for phi and zi
     */
	TPZFMatrix fC_star;
    TPZFMatrix fKeC_star; //K_star_inv*C_star_trans
    TPZStepSolver finv;
    TPZFMatrix fNullPivots;
	/// Number of equations of the substructure
	int fNEquations; 

	int NumEquations() const
	{
		return fNEquations;
	}
    /**
    *
    */
    TPZVec<int> fCoarseNodes;
    /**
     * Constraint definition associated with this substructure
     * Dado de entrada
     */
    TPZFMatrix fC; //C(i)
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
    TPZVec<int> fCoarseIndex; //R(ci)
    /**
     * Weights associated with each variable/equation
     * Calculado
     */
    TPZManVector<REAL> fWeights; //W(i) = diagonal(fWeights)
    /**
     * Stiffness matrix associated with the substructure
     * Dado de entrada
     */
    TPZAutoPointer<TPZMatrix> fStiffness; //K(i)
    /**
     * Stiffness matrix associated with the constraints
     * Calculado
     */
    TPZFMatrix fKCi; //K(ci)
    /**
     * Global vector indices of the equations/degrees of freedom
     * Dado de entrada
     */
    TPZVec<int> fGlobalIndex; //R(i)
	
	
	/// Indices of the equations in the global system corresponding to the local equations
	TPZVec<std::pair<int,int> > fGlobalEqs;
    /**
     * Inverted (LU or Cholesky) stiffness matrix
     * Calculado
     */
    TPZStepSolver fInvertedStiffness; //Ke, or K_star
    /**
     * Inverted (LU or Cholesky or LDLt) stiffness matrix for the internal degrees of freedom
     * Calculado
     */
    mutable TPZStepSolver fInvertedInternalStiffness; //R(Ii)*K(i)*R(Ii)invertida
    /**
     * Inverted restraint matrix
     * Calculado
     */
    TPZStepSolver fKCInvert; //K(c) invertida
    /**
     * internal nodes
     * Dado
     */
    TPZVec<int> fInternalEqs; //R(Ii)
	/**
	 * equations corresponding to boundary of the substructure
	 */
	TPZVec<int> fBoundaryEqs;
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
    TPZFMatrix fzi; //z(i)
	
	/**
	 * Solution vector which needs to be added to the converged system solution
	 * this variable is initialized and set in the AdjustResidual method
	 */
	TPZFMatrix fAdjustSolution;
	
};

#endif

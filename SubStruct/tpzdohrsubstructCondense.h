/**
 * @file
 * @brief Contains the TPZDohrSubstructCondense class which condenses matrix divided in sub structures.
 */

#ifndef TPZDOHRSUBSTRUCTCONDENSE_H
#define TPZDOHRSUBSTRUCTCONDENSE_H

#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "pzstepsolver.h"
#include "pzysmp.h"
#include "pzmatred.h"
#include "TPZSavable.h"

/**
 * @brief To condense matrix divided in sub structures. \ref substructure "Sub Structure"
 * @author Philippe Devloo
 * @ingroup substructure
 */
template<class TVar>
class TPZDohrSubstructCondense : public TPZSavable
{
	public:
            
            public:
int ClassId() const override;

            
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

		/** @brief Adjust the residual to reflect a static condensation */
		/** The residual corresponding to the internal nodes will be zeroed */
		void AdjustResidual(TPZFMatrix<TVar> &r_global);
		/** @brief It computes the local contribution to r(c). */
		/** The method LoadWeightedResidual must be called before this one. */
		void Contribute_rc_local(TPZFMatrix<TVar> &residual_local, TPZFMatrix<TVar> &rc_local) const;
		/** @brief It computes the local contribution to K(c) */
		void Contribute_Kc(TPZMatrix<TVar> &Kc, TPZVec<int> &coarseindex);
		/**
		 * @brief It computes the local contribution to v1.
		 * @param v1_local 
		 * @param invKc_rc is the product K(c)_inverted*r(c)
		 * Of course r(c) must be computed, using Contribute_rc(), before calling this method
		 */
		void Contribute_v1_local(TPZFMatrix<TVar> &v1_local, TPZFMatrix<TVar> &invKc_rc) const;
		/** @brief It computes the local contribution to v2. */
		void Contribute_v2_local(TPZFMatrix<TVar> &residual_local, TPZFMatrix<TVar> &v2_local);
		/**
		 * @brief It computes the local contribution to v(3)
		 * @param v3
		 * @param v1Plusv2 is the sum "v1 + v2"
		 */
		void Contribute_v3_local(TPZFMatrix<TVar> &v3, TPZFMatrix<TVar> &v1Plusv2) const;
		/**
		 * @brief 
		 */
		void Print(std::ostream &out) const;

		/** @brief compute the multiplication of the local stiffness matrix with the vector u */
		void ContributeKULocal(const TVar alpha, const TPZFMatrix<TVar> &u, TPZFMatrix<TVar> &z) const;
		/**
		 * @brief Computes the contribution of each substructure node to global Stiffness diagonal (or something like that).
		 * @param StiffnessDiag is the diagonal of the stiffness matrix
		 */
		void ContributeDiagonalLocal(TPZFMatrix<TVar> &StiffnessDiag);
		/**
		 * @brief Computes the weight matrix.
		 * @param StiffnessDiagLocal is the diagonal of the global matrix (local numbering)
		 */
		void ComputeWeightsLocal(TPZFMatrix<TVar> &StiffnessDiagLocal);
		/**
		 * @brief Computes the condensed right hand side for the substructure
		 * @param rhs the right hand side ordered external first
		 */
		void ContributeRhs(TPZFMatrix<TVar> &rhs);
		/** @brief Computes the global solution based on the interface solution */
		void UGlobal(TPZFMatrix<TVar> &Uext, TPZFMatrix<TVar> &UGlobal);
		/** @brief Initializes the substructure. */
		void Initialize();

		/** @brief Assemble the contribution to the coarse residual */
		void GetCoarseResidual(TPZFMatrix<TVar> &rc, TPZVec<int> &indices);

		/** @brief Assemble the coarse dof stiffness matrix */
		void GetCoarseStiffness(TPZMatrix<TVar> &stiff, TPZVec<int> &indices);

		/**
		 * @brief Apply a scatter permutation to the input vector using a scatter permutation
		 * output[permute[i]] = input[i-first], first <= i < last
		 */
		/** This method does not resize the elements */
		static void PermuteScatter(const TPZVec<int> &permute, const TPZFMatrix<TVar> &input, TPZFMatrix<TVar> &output, int first, int last);
		/**
		 * @brief Apply a gather permutation to the input vector using a scatter permutation
		 * output[i-first] = input[permute[i]], first <= i < last
		 */
		/** This method does not resize the elements */
		static void PermuteGather(const TPZVec<int> &permute, const TPZFMatrix<TVar> &input, TPZFMatrix<TVar> &output, int first, int last);

		/** @brief method for streaming the object to a stream */
		void Write(TPZStream &buf, int withclassid) const override;

		/** @brief method for reading the object for a stream */
		void Read(TPZStream &input, void *context) override;

	public:
		/** @brief It prepares the datas for solving systems for phi and zi */
		void PrepareSystems();
		/**
		 * @brief Solves the system for Phi and for v2 \n
		 * It stores the results in fPhiC and fzi
		 */
		void SolveSystemPhi();
		/** @brief Add the internal solution to the final result */
		void AddInternalSolution(TPZFMatrix<TVar> &sol);

		/** @brief Matrix problem which solves the phi and zi problems */
		TPZAutoPointer<TPZMatRed<TVar, TPZFMatrix<TVar> > > fMatRedComplete;
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

		/** @brief Vectors associated with each constraint */
		/** Computed */
		TPZFMatrix<TVar> fPhiC; //Phí(i)
		/** @brief Phi * W matrix and condensed to the equations which are part of the global system */
		TPZFMatrix<TVar> fPhiC_Weighted_Condensed;

		/**
		 * @brief Weights associated with each variable/equation
		 * 
		 * Computed
		 */
		TPZManVector<TVar> fWeights; //W(i) = diagonal(fWeights)
		/**
		 * @brief Stiffness matrix associated with the constraints
		 *
		 * Computed
		 */
		TPZFMatrix<TVar> fKCi; //K(ci)

		/** @brief Permutation vectors
		 *
		 * This map holds all permutation vectors for the equations
		 * The format of a permutation is that permute[ind] corresponds to the destination index of the element index ind of the original vector
		 * solout[permute[ind]] = solin[ind], ind = 0, n, 1
		 */
		std::map<std::pair<ENumbering, ENumbering> , TPZVec<int> > fPermutationsScatter;

		bool wasRealloc;

		/** @TODO Remove this - PerfTests Only */
		void ReallocMatRed() 
		{
			if(!wasRealloc) {
				std::cout << "TPZDohrSubstructCondense::ReallocMatRed\n";
				fMatRed.ReallocForNuma(0);
                fMatRed->ReallocSolver();
				fMatRedComplete.ReallocForNuma(0);
                
				wasRealloc = true;
			}
		}

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
		mutable TPZAutoPointer<TPZMatRed<TVar, TPZFMatrix<TVar> > > fMatRed; //R(Ii)*K(i)*R(Ii)invertida
		/** @brief Local load vector */
		TPZFMatrix<TVar> fLocalLoad; //F(i) - Actually it's a vector
		/** @brief Local weighted residual - \f$ W(i)*R(i)*r \f$ */
		TPZFMatrix<TVar> fLocalWeightedResidual; //W(i)*R(i)*r - Actually it's a vector

		/** @brief Solution vector which needs to be added to the converged system solution */ 
		/** This variable is initialized and set in the AdjustResidual method */
		TPZFMatrix<TVar> fAdjustSolution;
};

template<class TVar>
int TPZDohrSubstructCondense<TVar>::ClassId() const{
    return Hash("TPZDohrSubstructCondense") ^ ClassIdOrHash<TVar>() << 1;
}

#endif

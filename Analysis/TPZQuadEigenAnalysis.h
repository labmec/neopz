//
// Created by Francisco Teixeira Orlandini on 14/06/23.
//

#ifndef PZ_TPZQUADEIGENANALYSIS_H
#define PZ_TPZQUADEIGENANALYSIS_H
#include "TPZEigenAnalysisBase.h"     //For TPZEigenBanalysisBase
#include "tpzautopointer.h" //For TPZAutoPointer
#include "pzmatrix.h"       //For TPZFMatrix

template<typename TVar>
class TPZQuadEigenSolver;

/**
 * @ingroup Analysis
 * @brief Performs the Finite Element Analysis of a quadratic eigenvalue problem.

 The problem is of the form:
 -b^2 Mu + i b Lu + Ku = 0

 See TPZQuadEigenSolver for details.
 For now, this is the only available solver.
*/
class TPZQuadEigenAnalysis : public virtual TPZEigenAnalysisBase{
public:
    enum class Mat{K,L,M};//< Which matrix should be assembled
    /** @name Constructors*/
    /** @{ */
    /** @brief Create an empty TPZQuadEigenAnalysis object*/
    TPZQuadEigenAnalysis();
    /** @brief Create an TPZQuadEigenAnalysis object from one mesh pointer */
    TPZQuadEigenAnalysis(TPZCompMesh *mesh, const RenumType &rt = RenumType::EDefault, std::ostream &out = std::cout);
    /** @brief Create an TPZQuadEigenAnalysis object from one mesh auto pointer object */
    TPZQuadEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh, const RenumType &rt = RenumType::EDefault, std::ostream &out = std::cout);

    /** @} */
    /** @name FEM*/
    /** @{ */
    /** @brief Gets the eigensolver */
    template<class TVar>
    TPZQuadEigenSolver<TVar> &EigenSolver();
    /** @brief Set the solver
      @note In this function it will be checked if the solver is a TPZEigenSolver*/
    void SetSolver(const TPZSolver &solver) override;
    /** @brief Assemble the matrices associated with the EVP*/
    void Assemble() override;
    /** @brief Assemble one of the matrices associated with the EVP*/
    void AssembleMat(const TPZQuadEigenAnalysis::Mat mat);
    /** @brief Solve the EVP problem*/
    void Solve() override;
    /** @} */
    /** @name ReadWrite*/
    /**  @{ */
    //! Class unique identifier
    int ClassId() const override;
    /** @} */
protected:
    template<class TVar, Mat MAT>
    void AssembleT();
    template<class TVar>
    void SolveT();
    template<class TVar, Mat MAT>
    TPZAutoPointer<TPZMatrix<TVar>> GetSolverMat();
};

#define INSTANTIATE_TEMPLATES(TVar)                                            \
  extern template TPZQuadEigenSolver<TVar> &TPZQuadEigenAnalysis::EigenSolver<TVar>();

INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)
#undef INSTANTIATE_TEMPLATES

#endif //PZ_TPZQUADEIGENANALYSIS_H

//
// Created by Francisco Teixeira Orlandini on 11/17/17.
//

#ifndef PZ_TPZEIGENANALYSIS_H
#define PZ_TPZEIGENANALYSIS_H
#include "TPZEigenAnalysisBase.h"     //For TPZEigenBanalysisBase
#include "tpzautopointer.h" //For TPZAutoPointer
#include "pzmatrix.h"       //For TPZFMatrix

template<typename TVar>
class TPZLinearEigenSolver;

/**
 * @ingroup Analysis
 * @brief Performs the Finite Element Analysis of a (standard or generalised) eigenvalue problem.
 */
class TPZEigenAnalysis : public virtual TPZEigenAnalysisBase{
public:
    enum class Mat{A,B};//< Which matrix should be assembled
    /** @name Constructors*/
    /** @{ */
    /** @brief Create an empty TPZEigenAnalysis object*/
    TPZEigenAnalysis();

  /** @brief Create an TPZEigenAnalysis object from one mesh pointer */
  TPZEigenAnalysis(TPZCompMesh *mesh, const RenumType& renumtype = RenumType::EDefault, std::ostream &out = std::cout);
  /** @brief Create an TPZEigenAnalysis object from one mesh auto pointer object */
  TPZEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh, const RenumType& renumtype = RenumType::EDefault, std::ostream &out = std::cout);


    /** @} */
    /** @name FEM*/
    /** @{ */
    /** @brief Gets the eigensolver */
    template<class TVar>
    TPZLinearEigenSolver<TVar> &EigenSolver();
    /** @brief Set the solver
      @note In this function it will be checked if the solver is a TPZEigenSolver*/
    void SetSolver(const TPZSolver &solver) override;
    /** @brief Assemble the matrices associated with the EVP*/
    void Assemble() override;
    /** @brief Assemble one of the matrices associated with the EVP*/
    void AssembleMat(const TPZEigenAnalysis::Mat mat);

    int ClassId() const override;
private:
    template<class TVar, Mat MAT>
    void AssembleT();
    template<class TVar>
    void ConfigAssemble();
    bool fIsSetUp{false};
};

#define INSTANTIATE_TEMPLATES(TVar)                                            \
  extern template TPZLinearEigenSolver<TVar> &TPZEigenAnalysis::EigenSolver<TVar>();

INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)
#undef INSTANTIATE_TEMPLATES

#endif //PZ_TPZEIGENANALYSIS_H

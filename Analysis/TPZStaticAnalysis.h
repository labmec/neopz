#ifndef TPZEQSYSANALYSIS_H
#define TPZEQSYSANALYSIS_H
#include "TPZAnalysis.h"

template<class TVar>
class TPZMatrixSolver;
/**
 * @ingroup analysis
 * @brief Performs the Finite Element Analysis of a equation system.
 */
class TPZStaticAnalysis : public TPZAnalysis{
public:

  /** @brief Create an empty TPZStaticAnalysis object */
	TPZStaticAnalysis();

	/** @brief Create an TPZStaticAnalysis object from one mesh pointer */
	TPZStaticAnalysis(TPZCompMesh *mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);
    	
	/** @brief Create an TPZStaticAnalysis object from one mesh auto pointer object */
	TPZStaticAnalysis(TPZAutoPointer<TPZCompMesh> mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);
  /** @brief Returns the load vector */
	TPZSolutionMatrix &Rhs() { return fRhs;}

  /** @name Graphical */
  /** @{ */
  /** @brief Run and print the solution step by step */
	void AnimateRun(int64_t num_iter, int steps,
                  TPZVec<std::string> &scalnames,
                  TPZVec<std::string> &vecnames,
                  const std::string &plotfile);
  
  using TPZAnalysis::PostProcess;
  void PostProcess(int resolution, int dimension) override;
  /** @} */
  /** @brief Get the matrix solver */
    template<class TVar>
	TPZMatrixSolver<TVar> & MatrixSolver();

  /** @brief Assemble the stiffness matrix and load vector */
	 void Assemble() override;
	
	/** @brief Assemble the load vector */
	virtual void AssembleResidual();
	
	/** @brief Invert the stiffness matrix */
	void Solve() override;
  /** @brief Set the solver
      @note In this function it will be checked if the solver is a TPZMatrixSolver*/
  void SetSolver(const TPZSolver &solver) override;

  /** @name ReadWrite
      @{ */
  int ClassId() const override;
  
  void Write(TPZStream &buf, int withclassid) const override;

  void Read(TPZStream &buf, void *context) override;
  /** @} */
protected:
  /** @brief Load vector */
	TPZSolutionMatrix fRhs;
private:
  template <class TVar> void AssembleT();
  template <class TVar> void AssembleResidualT();
  template <class TVar> void SolveT();
  template <class TVar>
  void AnimateRunT(int64_t num_iter, int steps,
                   TPZVec<std::string> &scalnames,
                   TPZVec<std::string> &vecnames,
                   const std::string &plotfile);
};


#define INSTANTIATE_TEMPLATES(TVar)                                     \
  extern template                                                       \
  TPZMatrixSolver<TVar> &TPZStaticAnalysis::MatrixSolver<TVar>();

INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)
#undef INSTANTIATE_TEMPLATES

#endif
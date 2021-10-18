#ifndef TPZEQSYSANALYSIS_H
#define TPZEQSYSANALYSIS_H
#include "TPZAnalysis.h"

template<class TVar>
class TPZMatrixSolver;
/**
 * @ingroup analysis
 * @brief Performs the Finite Element Analysis of a equation system.
 */
class TPZLinearAnalysis : public TPZAnalysis{
public:

  /** @name Constructors */
  /** @{ */
  /** @brief Create an empty TPZLinearAnalysis object */
	TPZLinearAnalysis();

	/** @brief Create an TPZLinearAnalysis object from one mesh pointer */
	TPZLinearAnalysis(TPZCompMesh *mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);
    	
	/** @brief Create an TPZLinearAnalysis object from one mesh auto pointer object */
	TPZLinearAnalysis(TPZAutoPointer<TPZCompMesh> mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);
  /** @} */
  
  /** @name FEM */
  /** @{ */
  /** @brief Get the matrix solver */
    template<class TVar>
	TPZMatrixSolver<TVar> & MatrixSolver();

  /** @brief Returns the load vector */
	TPZSolutionMatrix &Rhs() { return fRhs;}

  /** @brief Assemble the stiffness matrix and load vector */
	 void Assemble() override;
	
	/** @brief Assemble the load vector */
	virtual void AssembleResidual();
	
	/** @brief Invert the stiffness matrix */
	void Solve() override;
  /** @brief Set the solver
      @note In this function it will be checked if the solver is a TPZMatrixSolver*/
  void SetSolver(const TPZSolver &solver) override;
  /** @} */

  /** @name Graphical */
  /** @{ */
  /** @brief Run and print the solution step by step */
	void AnimateRun(int64_t num_iter, int steps,
                  TPZVec<std::string> &scalnames,
                  TPZVec<std::string> &vecnames,
                  const std::string &plotfile);
  

  void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout) override{
    TPZAnalysis::PostProcess(loc,out);
  }//just to avoid hiding the parent overloaded method
  void PostProcess(int resolution) override{
    TPZAnalysis::PostProcess(resolution);
  }//just to avoid hiding the parent overloaded method
  
  void PostProcess(int resolution, int dimension) override;

  /** @brief Sets time used in OpenDX files */
	inline void SetTime(REAL time) {
    this->fTime = time;
  }
	/** @brief Gets time used in OpenDX files */
	inline REAL GetTime() const{
    return this->fTime;
  }
  /** @} */
  
  /** @name ReadWrite
      @{ */
  int ClassId() const override;
  
  void Write(TPZStream &buf, int withclassid) const override;

  void Read(TPZStream &buf, void *context) override;
  /** @} */
protected:
  /** @brief Load vector */
	TPZSolutionMatrix fRhs;
  /** @brief Time variable used for post-processing*/
	REAL fTime{0.};
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
  TPZMatrixSolver<TVar> &TPZLinearAnalysis::MatrixSolver<TVar>();

INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)
#undef INSTANTIATE_TEMPLATES

#endif
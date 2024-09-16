//
// Created by Francisco Teixeira Orlandini on 14/06/23.
//

#ifndef PZ_TPZEIGENANALYSISBASE_H
#define PZ_TPZEIGENANALYSISBASE_H
#include "TPZAnalysis.h"     //For TPZAnalysis
#include "tpzautopointer.h" //For TPZAutoPointer
#include "pzmatrix.h"       //For TPZFMatrix

template<typename TVar>
class TPZLinearEigenSolver;

/**
 * @ingroup Analysis
 * @brief Performs the Finite Element Analysis of a (standard or generalised) eigenvalue problem.
 */
class TPZEigenAnalysisBase : public TPZAnalysis{
public:
    /** @brief Create an empty TPZEigenAnalysisBase object*/
    TPZEigenAnalysisBase();
    /** @brief Create an TPZEigenAnalysisBase object from one mesh pointer */
    TPZEigenAnalysisBase(TPZCompMesh *mesh, const RenumType &rt = RenumType::EDefault, std::ostream &out = std::cout);
    /** @brief Create an TPZEigenAnalysisBase object from one mesh auto pointer object */
    TPZEigenAnalysisBase(TPZAutoPointer<TPZCompMesh> mesh, const RenumType &rt = RenumType::EDefault, std::ostream &out = std::cout);
    /** @brief Set the solver: derived classes should check the type*/
    void SetSolver(const TPZSolver &solver) override = 0;
    /** @brief Gets the eigenvectors calculated by the Solve method*/
    TPZFMatrix<CSTATE> &GetEigenvectors()
    {return fEigenvectors;}
    /** @brief Gets the eigenvalues  by the Solve method*/
    TPZVec<CSTATE> &GetEigenvalues()
    {return fEigenvalues;}
    /** @brief Gets the eigenvectors calculated by the Solve method*/
    const TPZFMatrix<CSTATE> &GetEigenvectors() const
    {return fEigenvectors;}
    /** @brief Gets the eigenvalues  by the Solve method*/
    const TPZVec<CSTATE> &GetEigenvalues() const
    {return fEigenvalues;}
    /** @brief Sets the eigenvectors*/
    void SetEigenvectors(const TPZFMatrix<CSTATE> &ev){fEigenvectors = ev;}
    /** @brief Sets the eigenvalues*/
    void SetEigenvalues(const TPZVec<CSTATE> &ev){fEigenvalues = ev;}
    /** @brief Set to compute eigenvectors or just eigenvalues*/
    inline void SetComputeEigenvectors(const bool opt){fCalcVectors=opt;}
    /** @brief Whether to compute eigenvectors or just eigenvalues*/
    inline bool ComputeEigenvectors() const{return fCalcVectors;}

    /** @brief Solve the EVP problem*/
    void Solve() override;
    
    /** @} */
    /** @name ReadWrite*/
    /**  @{ */
    //! Class unique identifier
    int ClassId() const override;
    //! Write attributes to TPZStream
    void Write(TPZStream &buf, int withclassid) const override;
    //! Read attributes from TPZStream to replicate instance from file
    void Read(TPZStream &buf, void *context) override;
    /** @} */
protected:
    /**
    * @brief Stores the computed eigenvalues
    */
    TPZManVector<CSTATE,10> fEigenvalues;
    /**
     * @brief Stores the computed eigenvectors
     */
    TPZFMatrix<CSTATE> fEigenvectors;
    //! Whether to compute eigenvectors
    bool fCalcVectors{true};
private:
    template<class TVar>
    void SolveT();
    /*the following methods are set to private until they are properly implemented*/
    using TPZAnalysis::PrePostProcessTable;
    using TPZAnalysis::PostProcessTable;
    using TPZAnalysis::PostProcessError;
    using TPZAnalysis::PostProcessErrorSerial;
    using TPZAnalysis::PostProcessErrorParallel;
};
#endif //PZ_TPZEIGENANALYSISBASE_H

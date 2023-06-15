#include "TPZQuadEigenAnalysis.h"
#include "TPZQuadEigenSolver.h"
#include "TPZMatQuadraticEigenVal.h"
#include "TPZSpStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZMaterial.h"
#include "pzcmesh.h"

TPZQuadEigenAnalysis::TPZQuadEigenAnalysis()
    : TPZRegisterClassId(&TPZQuadEigenAnalysis::ClassId),TPZEigenAnalysisBase()
{

}

TPZQuadEigenAnalysis::TPZQuadEigenAnalysis(TPZCompMesh *mesh,
                                           const RenumType &rt,
                                           std::ostream &out)
    : TPZEigenAnalysisBase(mesh, rt,out)
{

}

TPZQuadEigenAnalysis::TPZQuadEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh,
                                           const RenumType &rt,
                                           std::ostream &out)
    : TPZEigenAnalysisBase(mesh, rt, out)
{

}

template<class TVar>
TPZQuadEigenSolver<TVar> &TPZQuadEigenAnalysis::EigenSolver()
{
    const auto tmp = dynamic_cast<TPZQuadEigenSolver<TVar>*>(fSolver);
    if(fSolver && !tmp){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" incompatible Solver type! Aborting\n";
        DebugStop();
    }
    return *tmp;
}

void TPZQuadEigenAnalysis::SetSolver(const TPZSolver &solver)
{
    if(fSolver) delete fSolver;
    auto *tmpState =
        dynamic_cast<const TPZEigenSolver<STATE>*>(&solver);
    auto *tmpCState =
        dynamic_cast<const TPZEigenSolver<CSTATE>*>(&solver);
    if(tmpState && fSolType == EReal){
        fSolver = (TPZEigenSolver<STATE> *) solver.Clone();
        return;
    }
    else if(tmpCState && fSolType == EComplex){
        fSolver = (TPZEigenSolver<CSTATE> *) solver.Clone();
        return;
    }
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" Incompatible types!\n";
    PZError<<" Aborting...\n";
    DebugStop();
}

void TPZQuadEigenAnalysis::Assemble()
{
    if(fSolType == EReal){
        AssembleT<STATE,Mat::K>();
        //we want to avoid computing sparsity pattern of other matrices,
        //so we assign them (knowing that they will be zeroed)
        auto matK = this->GetSolverMat<STATE, Mat::K>();
        auto &eigSolver = this->EigenSolver<STATE>();
        auto matL = matK->Clone();
        eigSolver.SetMatrixL(matL);
        auto matM = matK->Clone();
        eigSolver.SetMatrixM(matM);
        
        AssembleT<STATE,Mat::L>();
        AssembleT<STATE,Mat::M>();
    }
    else{
        AssembleT<CSTATE,Mat::K>();
        //we want to avoid computing sparsity pattern of other matrices,
        //so we assign them (knowing that they will be zeroed)
        auto matK = this->GetSolverMat<CSTATE, Mat::K>();
        auto &eigSolver = this->EigenSolver<CSTATE>();
        auto matL = matK->Clone();
        eigSolver.SetMatrixL(matL);
        auto matM = matK->Clone();
        eigSolver.SetMatrixM(matM);
        
        AssembleT<CSTATE,Mat::L>();
        AssembleT<CSTATE,Mat::M>();
    }
}

void TPZQuadEigenAnalysis::AssembleMat(const TPZQuadEigenAnalysis::Mat mat){
    if(fSolType == EReal){
        switch(mat){
        case Mat::K:
            this->AssembleT<STATE,Mat::K>();
            break;
        case Mat::L:
            this->AssembleT<STATE,Mat::L>();
            break;
        case Mat::M:
            this->AssembleT<STATE,Mat::M>();
            break;
        }
    }
    else{
        switch(mat){
        case Mat::K:
            this->AssembleT<CSTATE,Mat::K>();
            break;
        case Mat::L:
            this->AssembleT<CSTATE,Mat::L>();
            break;
        case Mat::M:
            this->AssembleT<CSTATE,Mat::M>();
            break;
        }
    }
}

template<class TVar, TPZQuadEigenAnalysis::Mat MAT>
void TPZQuadEigenAnalysis::AssembleT()
{

    if (!this->fStructMatrix) {
#ifdef USING_MKL
        std::cout << "Setting default struct matrix: sparse(non-symmetric)"
                  << std::endl;
        TPZSpStructMatrix<TVar> defaultMatrix(fCompMesh);
#else
        std::cout << "Setting default struct matrix: skyline(non-symmetric)"
                  << std::endl;
        TPZSkylineNSymStructMatrix<TVar> defaultMatrix(fCompMesh);
#endif
        this->SetStructuralMatrix(defaultMatrix);
    }

    fStructMatrix->SetComputeRhs(false);
    //it wont be resized or anything.
    TPZFMatrix<TVar> dummyRhs;
    const auto sz = fStructMatrix->EquationFilter().NActiveEquations();
    auto &eigSolver = this->EigenSolver<TVar>();


    auto &materialVec = fCompMesh->MaterialVec();
    for (auto &&item : materialVec) {
        auto mat = (item.second);
        auto eigmat = dynamic_cast<TPZMatQuadraticEigenVal *>(mat);
        if (eigmat){
            if constexpr (MAT==TPZQuadEigenAnalysis::Mat::K){
                eigmat->SetMatrixK();
            }else if constexpr (MAT==TPZQuadEigenAnalysis::Mat::L){
                eigmat->SetMatrixL();
            }else if constexpr (MAT==TPZQuadEigenAnalysis::Mat::M){
                eigmat->SetMatrixM();
            }
        } else {
            PZError << __PRETTY_FUNCTION__;
            PZError << "\nERROR: cannot solve quadratic EVP if\n"
                    << "the materials do not have the interface:\n"
                    << "TPZMatQuadraticEigenVal\n"
                    << "Aborting...\n";
            DebugStop();
        }
    }
    auto mat = GetSolverMat<TVar,MAT>();
    
    if (mat && mat->Rows() == sz) {
        mat->Zero();
        fStructMatrix->Assemble(*mat, dummyRhs);
    } else {
        TPZAutoPointer<TPZMatrix<TVar>> mat =
            dynamic_cast<TPZMatrix<TVar>*>(
              fStructMatrix->CreateAssemble(dummyRhs));
    
        if constexpr (MAT==TPZQuadEigenAnalysis::Mat::K){
            eigSolver.SetMatrixK(mat);
        }else if constexpr (MAT==TPZQuadEigenAnalysis::Mat::L){
            eigSolver.SetMatrixL(mat);
        }else{
            eigSolver.SetMatrixM(mat);
        }
    }
}


void TPZQuadEigenAnalysis::Solve()
{
    if(fSolType == EReal)
        SolveT<STATE>();
    else
        SolveT<CSTATE>();
}


template<class TVar>
void TPZQuadEigenAnalysis::SolveT()
{
    const auto nEq = fCompMesh->NEquations();
    const auto nReducedEq = fStructMatrix->NReducedEquations();
    const bool isReduced = nReducedEq == nEq ? false : true;
    auto &eigSolver = this->EigenSolver<TVar>();
    const auto nev = eigSolver.NEigenpairs();
    fEigenvalues.Resize(nev);
    if (ComputeEigenvectors()) {
      TPZFMatrix<CTVar> *eigvectors = &fEigenvectors;
      if(isReduced){
        eigvectors = new TPZFMatrix<CTVar>;
      }
      //if there is no equation filter nReducedEq == numeq
      eigvectors->Redim(nReducedEq, nev);
      eigSolver.Solve(fEigenvalues, *eigvectors);
      if(isReduced){
        fEigenvectors.Redim(nEq,nev);
        fStructMatrix->EquationFilter().Scatter(*eigvectors,fEigenvectors);
        delete eigvectors;
      }
    } else {
      eigSolver.Solve(fEigenvalues);
    }
}




template<class TVar, TPZQuadEigenAnalysis::Mat MAT>
TPZAutoPointer<TPZMatrix<TVar>> TPZQuadEigenAnalysis::GetSolverMat()
{
    auto &eigSolver = this->EigenSolver<TVar>();
    if constexpr (MAT==TPZQuadEigenAnalysis::Mat::K){
        return eigSolver.MatrixK();
    }else if constexpr (MAT==TPZQuadEigenAnalysis::Mat::L){
        return eigSolver.MatrixL();
    }else{
        return eigSolver.MatrixM();
    }
}


int TPZQuadEigenAnalysis::ClassId() const
{
  return Hash("TPZQuadEigenAnalysis") ^
    TPZEigenAnalysisBase::ClassId() << 1;
}

#define INSTANTIATE_TEMPLATES(TVar)                                     \
  template TPZQuadEigenSolver<TVar> &TPZQuadEigenAnalysis::EigenSolver<TVar>(); \

INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)
#undef INSTANTIATE_TEMPLATES

template class TPZRestoreClass<TPZQuadEigenAnalysis>;
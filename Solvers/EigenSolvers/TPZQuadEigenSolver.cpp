#include "TPZQuadEigenSolver.h"
#include "TPZFMatrixRef.h"
#include "TPZLapackEigenSolver.h"
#include "TPZPardisoSolver.h"
#include "TPZYSMPPardiso.h"
#include "TPZSYSMPPardiso.h"
#include "TPZSimpleTimer.h"


template<class TVar>
void TPZQuadEigenSolver<TVar>::SetMatrixK(TPZAutoPointer<TPZMatrix<TVar>> K)
{
    bool e1{false}, e2{false};
    if(fL){
        if(fL->Rows() != K->Rows()){e1=true;}
    }
    if(fM){
        if(fM->Rows() != K->Rows()){e2=true;}
    }
    const bool error = e1 | e2;
    if(error){
        PZError<<__PRETTY_FUNCTION__
               <<"\nIncompatible dimensions! errors: "
               <<e1<<" "
               <<e2<<"\n"
               <<"Aborting..."<<std::endl;
        DebugStop();
    }
    fK=K;
}

template<class TVar>
void TPZQuadEigenSolver<TVar>::SetMatrixL(TPZAutoPointer<TPZMatrix<TVar>> L)
{
    bool e1{false}, e2{false};
    if(fK){
        if(fK->Rows() != L->Rows()){e1=true;}
    }
    if(fM){
        if(fM->Rows() != L->Rows()){e2=true;}
    }
    const bool error = e1 | e2;
    if(error){
        PZError<<__PRETTY_FUNCTION__
               <<"\nIncompatible dimensions! errors: "
               <<e1<<" "
               <<e2<<"\n"
               <<"Aborting..."<<std::endl;
        DebugStop();
    }
    fL=L;
}

template<class TVar>
void TPZQuadEigenSolver<TVar>::SetMatrixM(TPZAutoPointer<TPZMatrix<TVar>> M)
{
    bool e1{false}, e2{false};
    if(fK){
        if(fK->Rows() != M->Rows()){e1=true;}
    }
    if(fL){
        if(fL->Rows() != M->Rows()){e2=true;}
    }
    const bool error = e1 | e2;
    if(error){
        PZError<<__PRETTY_FUNCTION__
               <<"\nIncompatible dimensions! errors: "
               <<e1<<" "
               <<e2<<"\n"
               <<"Aborting..."<<std::endl;
        DebugStop();
    }
    fM=M;
}

template<class TVar>
void TPZQuadEigenSolver<TVar>::SetMatrices(TPZAutoPointer<TPZMatrix<TVar>> K,
                                           TPZAutoPointer<TPZMatrix<TVar>> L,
                                           TPZAutoPointer<TPZMatrix<TVar>> M)
{
    SetMatrixK(K);
    SetMatrixL(L);
    SetMatrixM(M);
}

template<class TVar>
void TPZQuadEigenSolver<TVar>::ComputeAndDecomposeDelta()
{

    //delta = -(s^2M+sL+K)

    //delta = K
    fDelta = fK->Clone();
    const auto s = this->Shift();

    auto tmp = fL->Clone();
    *tmp *= s;
    //delta = K + sL
    fDelta->Storage() += tmp->Storage();
    //@orlandini: at least for now, we cannot use
    // *tmp = *fM;
    tmp->Zero();
    tmp->Storage() += fM->Storage();
    *tmp *= s*s;
    //delta = K + sL+s^2M
    fDelta->Storage() += tmp->Storage();
    //delta = -(s^2M+sL+K)
    *fDelta *= -1;
    delete tmp;
    fDelta->Decompose(ELU);
}

template<class TVar>
void TPZQuadEigenSolver<TVar>::ApplyOperator(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &res) const
{
    /**dinv = delta^-1

       (Q-sP)^-1 P (u)  = ( dinv (Lu + sM u + Mw)
                   (w)    ( dinv (-Ku + sMw)
    **/

    const int n = this->SystemSize()/2;
    //let us split x into u and w
    TPZFMatrix<TVar> u(n, 1, &x.g(0, 0), n);
    TPZFMatrix<TVar> w(n, 1, &x.g(n, 0), n);

    res.Zero();
    //let us split res into u and w
    TPZFMatrix<TVar> res_u(n, 1, &res.g(0, 0), n);
    TPZFMatrix<TVar> res_w(n, 1, &res.g(n, 0), n);

    const auto s = this->Shift();
    //first we compute res u
    
    fM->Multiply(u,res_u);
    res_u *= s;
    fM->Multiply(w,fScratch);
    res_u += fScratch;

    fL->Multiply(u,fScratch);
    res_u+=fScratch;
    //now res w
    fM->Multiply(w,res_w);
    res_w *= s;
    fK->Multiply(u,fScratch);
    fScratch *= -1;
    res_w += fScratch;

    const auto dectype = fDelta->IsDecomposed();
    fDelta->SolveDirect(res_u,dectype);
    fDelta->SolveDirect(res_w,dectype);
}

template<class TVar>
void
TPZQuadEigenSolver<TVar>::PreSolve()
{
    ComputeAndDecomposeDelta();
    if(this->KrylovDim() == -1){
        this->SetKrylovDim(this->NEigenpairs()*10);
    }
}

template<class TVar>
void
TPZQuadEigenSolver<TVar>::TransformEigenvalues(TPZVec<CTVar> &w)
{
    const auto &s = this->fShift;
    for(auto &mappedw : w) mappedw = ((CTVar)1.0)/mappedw + s;
}


template<class TVar>
TPZQuadEigenSolver<TVar> *
TPZQuadEigenSolver<TVar>::Clone() const
{
    return new TPZQuadEigenSolver<TVar>(*this);
}
template<class TVar>
void TPZQuadEigenSolver<TVar>::ResetMatrix()
{
    TPZAutoPointer<TPZMatrix<TVar>> newK, newL, newM;
    fK = newK;
    fL = newL;
    fM = newM;
}
template<class TVar>
void TPZQuadEigenSolver<TVar>::SetNEigenpairs(int n)
{
    if(n < 1) n = 1;
    this->fNEigenpairs = n;
    if(n > this->fKrylovDim){
        this->fKrylovDim = 10 * n;
        std::cout<< "Adjusted Krylov dim to "<<this->fKrylovDim<<std::endl;
    }
}

template class TPZQuadEigenSolver<float>;
template class TPZQuadEigenSolver<double>;

template class TPZQuadEigenSolver<std::complex<float>>;
template class TPZQuadEigenSolver<std::complex<double>>;
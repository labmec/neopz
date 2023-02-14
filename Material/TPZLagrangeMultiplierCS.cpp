#include "TPZLagrangeMultiplierCS.h"
#include "pzaxestools.h"

#ifdef USING_LAPACK
#include "TPZLapack.h"
#endif

template<class TVar>
void TPZLagrangeMultiplierCS<TVar>::ContributeInterface(
    const TPZMaterialDataT<TVar> &data,
    const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
    const std::map<int, TPZMaterialDataT<TVar>> &dataright,
    REAL weight, TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef)
{
#ifdef PZDEBUG
    if(dataleft.size() != 1 || dataright.size() != 1) DebugStop();
#endif
    

    const auto *phiLPtr = &dataleft.begin()->second.phi;
    const auto *phiRPtr = &dataright.begin()->second.phi;
        
    TPZManVector<TVar,3> solLvec = dataleft.begin()->second.sol[0];
    TPZManVector<TVar,3> solRvec = dataright.begin()->second.sol[0];
    TVar solL = solLvec[0];
    TVar solR = solRvec[0];

    const TPZFMatrix<REAL> &phiL = *phiLPtr;
    const TPZFMatrix<REAL> &phiR = *phiRPtr;
    
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
    static int count  = 0;

    if((nrowl+nrowr)*fNStateVariables != ek.Rows() && count < 20)
    {
        std::cout<<"ek.Rows() "<< ek.Rows()<<
        " nrowl " << nrowl <<
        " nrowr " << nrowr << " may give wrong result " << std::endl;
        count++;
    }

    int secondblock = ek.Rows()-phiR.Rows()*fNStateVariables;
    int il,jl,ir,jr;
    
    // 3) phi_I_left, phi_J_right
    for(il=0; il<nrowl; il++) {
        for (int ist=0; ist<fNStateVariables; ist++) {
            ef(fNStateVariables*il+ist,0) += weight * fMultiplier * phiL.GetVal(il,0) * solR;
        }
        for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL.GetVal(il,0) * phiR.GetVal(jr,0));
            }
        }
    }
    
    //    // 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        for (int ist=0; ist<fNStateVariables; ist++) {
            ef(ir*fNStateVariables+ist+secondblock,0) += weight * fMultiplier * phiR.GetVal(ir,0) * solL;
        }
        for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR.GetVal(ir,0) * phiL.GetVal(jl,0));
            }
        }
    }
}

template<>
void TPZLagrangeMultiplierCS<STATE>::ContributeInterface(
                                                         const TPZMaterialDataT<STATE> &data,
                                                         const std::map<int, TPZMaterialDataT<STATE>> &dataleft,
                                                         const std::map<int, TPZMaterialDataT<STATE>> &dataright,
                                                         REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
#ifdef PZDEBUG
    if(dataleft.size() != 1 || dataright.size() != 1)
        DebugStop();
#endif


        //const TPZFMatrix<REAL> *phiL = &(dataleft.begin()->second.phi);
        TPZFMatrix<REAL> phiLdummy = dataleft.begin()->second.phi;
        TPZFMatrix<REAL> *phiL = &phiLdummy;

        TPZFMatrix<REAL> phiRdummy = dataright.begin()->second.phi;
        TPZFMatrix<REAL> *phiR = &phiRdummy;
    
        TPZManVector<REAL,3> solLvec = dataleft.begin()->second.sol[0];
        TPZManVector<REAL,3> solRvec = dataright.begin()->second.sol[0];
        REAL solL = solLvec[0];
        REAL solR = solRvec[0];

        int nrowl = phiL->Rows();
        int nrowr = phiR->Rows();
        static int count  = 0;

        if((nrowl+nrowr)*fNStateVariables != ek.Rows() && count < 20)
        {
            std::cout<<"ek.Rows() "<< ek.Rows()<<
                     " nrowl " << nrowl <<
                     " nrowr " << nrowr << " may give wrong result "<< "fNStateVariables\t"<< fNStateVariables << "count " << count  << std::endl;
            count++;
        }
  
#ifdef USING_LAPACK2
    TPZFMatrix<REAL> phiLBlas(fNStateVariables,fNStateVariables*nrowl,0.);
    TPZFMatrix<REAL> phiRBlas(fNStateVariables,fNStateVariables*nrowr,0.);
    for (int i = 0; i < nrowl; i++) {
        for (int j = 0; j < fNStateVariables; j++) {
            phiLBlas(j,j+i*fNStateVariables) = phiLdummy(i,0);
        }
    }
    for (int i = 0; i < nrowr; i++) {
        for (int j = 0; j < fNStateVariables; j++) {
            phiRBlas(j,j+i*fNStateVariables) = phiRdummy(i,0);
        }
    }
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k;
        m = phiLBlas.Cols(); // number of rows of mat op(A)
        n = phiRBlas.Cols(); // number of cols of mat op(B)
        k = phiLBlas.Rows(); // number of cols of mat op(A)
        alpha = weight * fMultiplier;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDC = nrowl*fNStateVariables+nrowr*fNStateVariables; // first dimension of C
        LDA = k; // if not transpose max( 1, m ), otherwise max( 1, k ).
        LDB = k; // if not transpose max( 1, k ), otherwise max( 1, n ).
        C = &ek(0,nrowl*fNStateVariables);
        A = &phiLBlas(0,0);
        B = &phiRBlas(0,0);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    {

        double *A, *B, *C;
        double alpha, beta;
        int m,n,k;
        m = phiRBlas.Cols(); // number of rows of mat op(A)
        n = phiLBlas.Cols(); // number of cols of mat op(B)
        k = phiRBlas.Rows(); // number of cols of mat op(A)
        alpha = weight * fMultiplier ;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDC = nrowl*fNStateVariables+nrowr*fNStateVariables;
        LDA = k; // if not transpose max( 1, m ), otherwise max( 1, k ).
        LDB = k; // if not transpose max( 1, k ), otherwise max( 1, n ).
        C = &ek(nrowl*fNStateVariables,0);
        B = &phiLBlas(0,0);
        A = &phiRBlas(0,0);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }

#else
    
  
        
        int secondblock = ek.Rows()-phiR->Rows()*fNStateVariables;
        int il,jl,ir,jr;

        // 3) phi_I_left, phi_J_right
        for(il=0; il<nrowl; il++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ef(fNStateVariables*il+ist,0) += weight *(-1.0)* fMultiplier * phiL->GetVal(il,0) * solR;
            }
            for(jr=0; jr<nrowr; jr++) {
                for (int ist=0; ist<fNStateVariables; ist++) {
                    double L = *(&phiLdummy(0,0)+il);
                    double R = *(&phiRdummy(0,0)+jr);
                    ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * L * R;
                }
            }
        }

        //    // 4) phi_I_right, phi_J_left
        for(ir=0; ir<nrowr; ir++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ef(ir*fNStateVariables+ist+secondblock,0) += weight * (-1.0)*fMultiplier * phiR->GetVal(ir,0) * solL;
            }
            for(jl=0; jl<nrowl; jl++) {
                for (int ist=0; ist<fNStateVariables; ist++) {
                    double L = *(&phiLdummy(0,0)+jl);
                    double R = *(&phiRdummy(0,0)+ir);
                    ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (R * L);

                }

            }
        }
    
#endif
}


template<class TVar>
void TPZLagrangeMultiplierCS<TVar>::FillDataRequirementsInterface(
    TPZMaterialDataT<TVar> &data,
    std::map<int, TPZMaterialDataT<TVar>> &datavec_left,
    std::map<int, TPZMaterialDataT<TVar>> &datavec_right)
{
    data.SetAllRequirements(false);
    for(auto it = datavec_left.begin() ; it != datavec_left.end() ; it++){
        it->second.fNeedsSol = true;
    }
    for(auto it = datavec_right.begin() ; it != datavec_right.end() ; it++){
        it->second.fNeedsSol = true;
    }
}
// print the data in human readable form
template<class TVar>
void TPZLagrangeMultiplierCS<TVar>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TBase::Print(out);
    out << "NStateVariables " << this->fNStateVariables << std::endl;
    out << "fDimension " << this->fDimension << std::endl;
    out << "fMultiplier " << this->fMultiplier << std::endl;
}

template<class TVar>
int TPZLagrangeMultiplierCS<TVar>::ClassId() const{
    return Hash("TPZLagrangeMultiplierCS") ^
        TBase::ClassId() << 1;
}


template<class TVar>
void TPZLagrangeMultiplierCS<TVar>::Write(TPZStream &buf, int withclassid) const
{
    TBase::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
}


template<class TVar>
void TPZLagrangeMultiplierCS<TVar>::Read(TPZStream &buf, void *context)
{
    TBase::Read(buf, context);
    buf.Read(&fNStateVariables);
}

template class TPZLagrangeMultiplierCS<STATE>;
template class TPZLagrangeMultiplierCS<CSTATE>;

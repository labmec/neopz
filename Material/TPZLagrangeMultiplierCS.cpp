#include "TPZLagrangeMultiplierCS.h"
#include "pzaxestools.h"


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
        for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL.GetVal(il,0) * phiR.GetVal(jr,0));
            }
        }
    }
    
    //	// 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR.GetVal(ir,0) * phiL.GetVal(jl,0));
            }
        }
    }

}

template<class TVar>
void TPZLagrangeMultiplierCS<TVar>::FillDataRequirementsInterface(
    TPZMaterialDataT<TVar> &data,
    std::map<int, TPZMaterialDataT<TVar>> &datavec_left,
    std::map<int, TPZMaterialDataT<TVar>> &datavec_right)
{
    data.SetAllRequirements(false);
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
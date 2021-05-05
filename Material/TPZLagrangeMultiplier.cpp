#include "TPZLagrangeMultiplier.h"
#include "pzaxestools.h"


template<class TVar>
void TPZLagrangeMultiplier<TVar>::ContributeInterface(
    const TPZMaterialDataT<TVar> &data,
    const TPZMaterialDataT<TVar> &dataleft,
    const TPZMaterialDataT<TVar> &dataright,
    REAL weight, TPZFMatrix<TVar> &ek,
    TPZFMatrix<TVar> &ef)
{
	const TPZFMatrix<REAL> &phiL = dataleft.phi;
	const TPZFMatrix<REAL> &phiR = dataright.phi;

	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
#ifdef PZDEBUG
    if(phiL.Rows()*fNStateVariables+phiR.Rows()*fNStateVariables != ek.Rows())
    {
        DebugStop();
    }
#endif
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

// print the data in human readable form
template<class TVar>
void TPZLagrangeMultiplier<TVar>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMaterial::Print(out);
    out << "NStateVariables " << this->fNStateVariables << std::endl;
    out << "fDimension " << this->fDimension << std::endl;
    out << "fMultiplier " << this->fMultiplier << std::endl;
}

template<class TVar>
int TPZLagrangeMultiplier<TVar>::ClassId() const{
    return Hash("TPZLagrangeMultiplier") ^
        TBase::ClassId() << 1;
}


template<class TVar>
void TPZLagrangeMultiplier<TVar>::Write(TPZStream &buf, int withclassid) const
{
    TBase::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
}


template<class TVar>
void TPZLagrangeMultiplier<TVar>::Read(TPZStream &buf, void *context)
{
    TBase::Read(buf, context);
    buf.Read(&fNStateVariables);
    
}

template class TPZLagrangeMultiplier<STATE>;
template class TPZLagrangeMultiplier<CSTATE>;
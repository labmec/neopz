#include "TPZHDivProjection.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

template<class TVar>
TPZHDivProjection<TVar>::TPZHDivProjection(int id, int dim) :
    TPZRegisterClassId(&TPZHDivProjection::ClassId),
    TBase(id), fDim(dim)
{
}

template<class TVar>
TPZMaterial * TPZHDivProjection<TVar>::NewMaterial() const{
	return new TPZHDivProjection(*this);
}

template<class TVar>
void TPZHDivProjection<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                       REAL weight,
                                       TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	
	const int nshape = data.phi.Rows();
    //last position of solLoc is the divergence
    TPZManVector<TVar,4> solLoc(4);
    if(!this->HasForcingFunction()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" a forcing function (solution to be projected)\n";
        PZError<<"has not been set! Aborting...\n";
        DebugStop();
    }
    this->fForcingFunction(data.x,solLoc);
    const TVar divSol = solLoc[3];

    /*
     * HDiv approximation spaces in NeoPZ are created 
     * from the combination of scalar functions with 
     * constant vector fields
     */
    
    // Setting the scalar phis
    const TPZFMatrix<REAL> &phiQ = data.phi;
    const TPZFMatrix<REAL> &dphiQ = data.dphix;
    
    const int phrq = data.fVecShapeIndex.NElements();
    
    //Calculate the matrix contribution for flux. Matrix A
    TPZFNMatrix<3,REAL> ivec(3,1,0.);
    TPZFNMatrix<3,REAL> jvec(3,1,0.);
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        const int ivecind = data.fVecShapeIndex[iq].first;
        const int ishapeind = data.fVecShapeIndex[iq].second;
        
        for(int id=0; id<3; id++){
            ivec(id,0) = data.fDeformedDirections.GetVal(id,ivecind);
        }
        TVar ff = 0.;
        
        for (int i=0; i<3; i++) {
            ff += ivec(i,0)*solLoc[i];
        }
        
        ef(iq,0) += weight*ff*phiQ.GetVal(ishapeind,0);
        REAL divqi = 0.;
        TPZFNMatrix<3,REAL> axesvec(3,1,0.);
        data.axes.Multiply(ivec,axesvec);
        //computin div(qj)
        for(int iloc=0; iloc<fDim; iloc++){
            divqi += axesvec(iloc,0)*dphiQ.GetVal(iloc,ishapeind);
        }
        for (int jq=0; jq<phrq; jq++)
        {
            int jvecind = data.fVecShapeIndex[jq].first;
            int jshapeind = data.fVecShapeIndex[jq].second;
            
            for(int id=0; id<3; id++){
                jvec(id,0) = data.fDeformedDirections.GetVal(id,jvecind);
            }
            
            //jvecZ.Print("mat1 = ");
            const REAL prod1 = ivec(0,0)*jvec(0,0) + ivec(1,0)*jvec(1,0) + ivec(2,0)*jvec(2,0);
            ek(iq,jq) +=
                    weight*phiQ.GetVal(ishapeind,0)*phiQ.GetVal(jshapeind,0)*prod1;

            REAL divqj = 0.;
            TPZFNMatrix<3,REAL> axesvec(3,1,0.);
            data.axes.Multiply(jvec,axesvec);
            //computin div(qj)
            for(int jloc=0; jloc<fDim; jloc++){
                divqj += axesvec(jloc,0)*dphiQ.GetVal(jloc,jshapeind);
            }
            ek(iq,jq) += weight*divqi*divqj;
        }
    }
}

template<class TVar>
void TPZHDivProjection<TVar>::ContributeBC(const TPZMaterialDataT<TVar> &data,
                                         REAL weight,
                                         TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                                         TPZBndCondT<TVar> &bc)
{
	const auto &phi = data.phi;
	const int phr = phi.Rows();

    const TVar v2 = [&bc = std::as_const(bc),
                     &data = std::as_const(data)](){
        if(bc.HasForcingFunctionBC()){
            TPZManVector<TVar> res(3);
            TPZFNMatrix<9,TVar> dummy(3,3,0.);
            bc.ForcingFunctionBC()(data.x,res,dummy);
            return res[0];
        }else {
            return bc.Val2()[0];
        }
    }();
    
	switch (bc.Type()){
			// Dirichlet condition
		case 0 : 
            for(auto in = 0; in < phr; in++) {
                ef(in,0) +=
                    (TVar)TPZMaterial::fBigNumber *
                    v2 *phi.GetVal(in,0) * (TVar)weight;
					for (auto jn = 0 ; jn < phr; jn++) {
						ek(in,jn) +=
                            TPZMaterial::fBigNumber *
                            phi.GetVal(in,0) * phi.GetVal(jn,0) * weight;
					}//jn
				}//in
			break;
			// Neumann condition
		case 1 :
            for(auto in = 0 ; in < phr; in++) {
                ef(in,0) +=
                    v2 * (TVar)fScale * (TVar)phi.GetVal(in,0) * (TVar)weight;
            }//in
			break;
        case 2 :			// mixed condition
            for(int iq = 0; iq < phr; iq++) {
                ef(iq,0) += v2*phi.GetVal(iq,0)*weight;
                for (int jq = 0; jq < phr; jq++) {
                    ek(iq,jq) += weight*bc.Val1().GetVal(0,0)*phi.GetVal(iq,0)*phi.GetVal(jq,0);
                }
            }
		default:{
			std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
		}
	}//switch
	
}
template<class TVar>
void TPZHDivProjection<TVar>::GetSolDimensions(uint64_t &u_len,
                                             uint64_t &du_row,
                                             uint64_t &du_col) const
{
    u_len=3;
    du_row=1;
    du_col=1;
}


template<class TVar>
int TPZHDivProjection<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return ESolution;
    if(!strcmp("Divergence",name.c_str())) return EDivergence;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int TPZHDivProjection<TVar>::NSolutionVariables(int var) const{
	if(var == ESolution) return fDim;
    if (var == EDivergence) {
        return 1;
    }
	
    return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void TPZHDivProjection<TVar>::Solution(const TPZMaterialDataT<TVar> &data,
                                     int var, TPZVec<TVar> &solOut)
{
    const auto &sol = data.sol[0];
    const auto &dsol = data.dsol[0];
	if (var == ESolution){
        solOut.Resize(sol.size());
        for (int i=0; i<sol.size(); i++) {
            solOut[i] = sol[i];
        }
		return;
	}
    if (var == EDivergence) {
        solOut.Resize(fDim);
        for (int i=0; i<fDim; i++) {
            solOut[i] = dsol.GetVal(i,0);
        }
        return;
    }
}

template<class TVar>
void TPZHDivProjection<TVar>::Errors(const TPZMaterialDataT<TVar> &data,
                                     TPZVec<REAL> &values)
{
#ifdef PZDEBUG
    if(!this->HasExactSol()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe material has no associated exact solution. Aborting...\n";
        DebugStop();
    }
#endif
    const auto &x = data.x;
    const auto &u = data.sol[0];
    const auto &du = data.divsol[0];
    
    TPZManVector<TVar,3> u_exact={0.,0.,0.};
    TPZFNMatrix<3,TVar> du_exact(1,1,0.);
    
    this->ExactSol()(x,u_exact,du_exact);
    
    values.Resize(this->NEvalErrors());
    values.Fill(0.0);
    
    //values[0] : error in HDiv norm
    //values[1] : error in L2 norm
    //values[2] : error in HDiv semi-norm
    TVar diff = (u[0] - u_exact[0]);
    if constexpr (is_complex<TVar>::value){
        values[1]  = std::real((diff*std::conj(diff)));
    }else{
        values[1]  = diff*diff;
    }

    if constexpr (is_complex<TVar>::value){
        values[2] = std::real(std::conj(du_exact(0,0)) * du[0]);
    }else{
        values[2] = du_exact(0,0) *du[0];
    }
    
    values[0]  = values[1]+values[2];
}

template<class TVar>
int TPZHDivProjection<TVar>::IntegrationRuleOrder(const int elPMaxOrder) const
{
    //for order k functions of order k+1
    return  TBase::IntegrationRuleOrder(elPMaxOrder+1);
}

template<class TVar>
int TPZHDivProjection<TVar>::ClassId() const{
    return Hash("TPZHDivProjection") ^ TBase::ClassId() << 1;
}


template class TPZHDivProjection<STATE>;
template class TPZHDivProjection<CSTATE>;
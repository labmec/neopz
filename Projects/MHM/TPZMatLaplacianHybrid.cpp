/**
 * @file
 * @brief Contains implementations of the TPZMatLaplacian methods.
 */

#include "TPZMatLaplacianHybrid.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"
#include "pzaxestools.h"

#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.TPZMatLaplacian"));
#endif


using namespace std;

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid(int nummat, int dim) : TPZMatLaplacian(nummat,dim)
{
}

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid() : TPZMatLaplacian()
{
}


TPZMatLaplacianHybrid & TPZMatLaplacianHybrid::operator=(const TPZMatLaplacianHybrid &copy){
	TPZMatLaplacian::operator = (copy);
	return *this;
}


TPZMatLaplacianHybrid::~TPZMatLaplacianHybrid() {
}


void TPZMatLaplacianHybrid::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
    TPZMatLaplacian::Print(out);
}

void TPZMatLaplacianHybrid::Contribute(TPZVec<TPZMaterialData> &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
	
    
    TPZFMatrix<REAL>  &phi = data[0].phi;
    TPZFMatrix<REAL> &dphi = data[0].dphix;
    TPZVec<REAL>  &x = data[0].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE XfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        TPZFMatrix<STATE> dres(Dimension(),1);
        fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
        XfLoc = res[0];
    }
    
    STATE KPerm = fK;
    if (fPermeabilityFunction) {
        TPZFNMatrix<9,STATE> perm, invperm;
        TPZManVector<STATE,3> func;
        TPZFNMatrix<18,STATE> dfunc(6,3,0.);
        fPermeabilityFunction->Execute(x, func, dfunc);
        KPerm = dfunc(0,0);
    }
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * XfLoc * (STATE)phi(in,0);
        for( int jn = 0; jn < phr; jn++ ) {
            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight * (+KPerm * (STATE)( dphi(kd,in) * dphi(kd,jn)));
            }
        }
    }
    //constantes 1
    for (int in=0; in<phr; in++) {
        ek(phr,in) -= (STATE)phi(in,0)*weight;
        ek(in,phr) -= (STATE)phi(in,0)*weight;
    }
    
    //constante 2
    ek(phr+1,phr) += weight;
    ek(phr,phr+1) += weight;
    //ek(phr,phr) -= 1.*(STATE)weight;
    
#ifdef PZDEBUG
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
#endif
    
}

int TPZMatLaplacianHybrid::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str()))        return  31;
	if(!strcmp("ExactSolution",name.c_str()))   return  32;
    if(!strcmp("PressureConstante",name.c_str())) return 33;
	return TPZMatLaplacian::VariableIndex(name);
}

int TPZMatLaplacianHybrid::NSolutionVariables(int var){
	if(var==31 || var==32 || var==33) return 1;
	return TPZMatLaplacian::NSolutionVariables(var);
}

void TPZMatLaplacianHybrid::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
#ifdef PZDEBUG
    if (bc.Type() > 2 ) {
        std::cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
        DebugStop();
    }
#endif
    
    TPZFMatrix<REAL>  &phiQ = data.phi;
    int phrq = phiQ.Rows();
    
    REAL v2;
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(this->fDim,1);
        bc.ForcingFunction()->Execute(data.x,res,gradu);
        v2 = res[0];
    }else
    {
        v2 = bc.Val2()(0,0);
    }
    
    switch (bc.Type()) {
        case 0 :		// Dirichlet condition
                        //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2*phiQ(iq,0)*weight;
            }
            break;
            
        case 1 :			// Neumann condition
                            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                ef(iq,0)+= gBigNumber*v2*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                }
            }
            break;
            
        case 2 :			// mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
                ef(iq,0) += v2*phiQ(iq,0)*weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq,jq) += weight*bc.Val1()(0,0)*phiQ(iq,0)*phiQ(jq,0);
                }
            }
            
            break;
    }
    
}



void TPZMatLaplacianHybrid::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    if (var == 31) {
        Solout[0] = datavec[0].sol[0][0];
        return;
    }
    
    
    if (var == 33) {
        Solout[0] = datavec[2].sol[0][0];
        return;
    }
    
    TPZVec<STATE> solExata(1);
	TPZFMatrix<STATE> flux(2,1);
    //Exact soluion
	if(var == 32){
		fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}
    TPZMatLaplacian::Solution(datavec[0], var, Solout);
}


void TPZMatLaplacianHybrid::Write(TPZStream &buf, int withclassid) const{
	TPZMatLaplacian::Write(buf, withclassid);
}

void TPZMatLaplacianHybrid::Read(TPZStream &buf, void *context){
	TPZMatLaplacian::Read(buf, context);
}

int TPZMatLaplacianHybrid::ClassId() const{
    return Hash("TPZMatLaplacianHybrid") ^ TPZMatLaplacian::ClassId() << 1;
}

template class TPZRestoreClass<TPZMatLaplacianHybrid>;

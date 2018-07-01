/**
 * @file
 * @brief Contains implementations of the TPZMatLaplacian methods.
 */

#include "TPZMatLaplacianLagrange.h"
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

TPZMatLaplacianLagrange::TPZMatLaplacianLagrange(int nummat, int dim) : TPZMatLaplacian(nummat,dim)
{
}

TPZMatLaplacianLagrange::TPZMatLaplacianLagrange() : TPZMatLaplacian()
{
}


TPZMatLaplacianLagrange & TPZMatLaplacianLagrange::operator=(const TPZMatLaplacianLagrange &copy){
	TPZMatLaplacian::operator = (copy);
	return *this;
}


TPZMatLaplacianLagrange::~TPZMatLaplacianLagrange() {
}


void TPZMatLaplacianLagrange::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
    TPZMatLaplacian::Print(out);
}

void TPZMatLaplacianLagrange::Contribute(TPZVec<TPZMaterialData> &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
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

int TPZMatLaplacianLagrange::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str()))        return  31;
	if(!strcmp("ExactSolution",name.c_str()))   return  32;
    if(!strcmp("PressureConstante",name.c_str())) return 33;
	return TPZMatLaplacian::VariableIndex(name);
}

int TPZMatLaplacianLagrange::NSolutionVariables(int var){
	if(var==31 || var==32 || var==33) return 1;
	return TPZMatLaplacian::NSolutionVariables(var);
}



void TPZMatLaplacianLagrange::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
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

void TPZMatLaplacianLagrange::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    TPZMaterial::Errors(data[0],u_exact,du_exact,errors);
}



void TPZMatLaplacianLagrange::Write(TPZStream &buf, int withclassid){
	TPZMatLaplacian::Write(buf, withclassid);
}

void TPZMatLaplacianLagrange::Read(TPZStream &buf, void *context){
	TPZMatLaplacian::Read(buf, context);
}

template class TPZRestoreClass<TPZMatLaplacianLagrange,TPZMatLaplacianLagrangeID>;

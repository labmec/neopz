/*
 * @file Helmhotz1D.cpp
 * @brief Contains the implementation of the TPZHelmholtzComplex1D methods.
 */

#include "pzhelmholtzcomplex1D.h"
#include "pzlog.h"

#ifdef LOG4CXX
        static LoggerPtr logger(Logger::getLogger("pz.helmholtz1D"));
#endif

TPZHelmholtzComplex1D::TPZHelmholtzComplex1D(int id,int dim) : TPZMat1dLin(id) {
}

TPZHelmholtzComplex1D::TPZHelmholtzComplex1D(const TPZHelmholtzComplex1D &helm) : TPZMat1dLin(helm) {
}

TPZHelmholtzComplex1D::~TPZHelmholtzComplex1D() {
}

void TPZHelmholtzComplex1D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

    TPZManVector<STATE,2> alphaval(1), betaval(1), phiaval(1);
    fAlpha->Execute(data.x, alphaval);
    fBeta->Execute(data.x, betaval);
    fPhi->Execute(data.x, phiaval);
    
    #ifdef LOG4CXX 
    {
        std::stringstream sout;        
        sout << "Coordenate x: " << data.x << " alpha = " << alphaval << " beta = " << betaval << " phi = " << phiaval;
        LOGPZ_DEBUG(logger, sout.str());
    }
    #endif
    
    TPZFNMatrix<4, STATE> xk(1, 1), xb(1, 1), xc(1, 1, 0.), xf(1, 1);
    xk(0,0) = alphaval[0];    
    xb(0,0) = betaval[0];    
    xf(0,0) = -phiaval[0];   
       
    SetMaterial(xk, xc, xb, xf);
    TPZMat1dLin::Contribute(data, weight, ek, ef);
}

/** @brief Returns the variable index associated with the name */
int TPZHelmholtzComplex1D::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(), "state")) return 0; 
	return TPZMat1dLin::VariableIndex(name);
}
    
/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZHelmholtzComplex1D::NSolutionVariables(int var) {
	if(var == 0) return 2 * NStateVariables();
	return TPZMat1dLin::NSolutionVariables(var);
}
    
/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZHelmholtzComplex1D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) {
	if(var == 0) {
		Solout[0] = data.sol[0][0].real();
		Solout[1] = data.sol[0][0].imag();
		return;
	}
	TPZMaterial::Solution(data, var, Solout);
}


/*
void TPZHelmholtz1D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc) {

    
}
*/
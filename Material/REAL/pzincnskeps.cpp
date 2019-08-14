/**
 * @file
 * @brief Contains implementations of the TPZIncNavierStokesKEps methods.
 */

#include "pzincnskeps.h"    

TPZIncNavierStokesKEps::TPZIncNavierStokesKEps(int id, int dimension):
TPZRegisterClassId(&TPZIncNavierStokesKEps::ClassId),
TPZMaterial(id){
	this->fDimension = dimension;
}

TPZIncNavierStokesKEps::~TPZIncNavierStokesKEps(){}

void TPZIncNavierStokesKEps::SetParameters(STATE MU, STATE RHO, STATE Cmu, STATE SigmaK, STATE SigmaEps, STATE Cepsilon1, STATE Cepsilon2, TPZVec<STATE> &BodyForce ){
	fMU = MU;
	fRHO = RHO;
	fCmu = Cmu;
	fSigmaK = SigmaK;
	fSigmaEps = SigmaEps;
	fCepsilon1 = Cepsilon1;
	fCepsilon2 = Cepsilon2;
	fBodyForce = BodyForce;
}

void TPZIncNavierStokesKEps::GetParameters(STATE &MU, STATE &RHO, STATE &Cmu, STATE &SigmaK, STATE &SigmaEps, STATE &Cepsilon1, STATE &Cepsilon2, TPZVec<STATE> &BodyForce ){
	MU = fMU;
	RHO = fRHO;
	Cmu = fCmu;
	SigmaK = fSigmaK;
	SigmaEps = fSigmaEps;
	Cepsilon1 = fCepsilon1;
	Cepsilon2 = fCepsilon2;
	BodyForce = fBodyForce;
}

int TPZIncNavierStokesKEps::Dimension() const {
	return this->fDimension;
}

int TPZIncNavierStokesKEps::NStateVariables() const {
	return this->Dimension() + 3; //Vi + p + K + Eps
}

void TPZIncNavierStokesKEps::Print(std::ostream &out){
	TPZMaterial::Print(out);    
}

int TPZIncNavierStokesKEps::NSolutionVariables(int var){
	if ( (var == EK) || (var == EEpsilon) || (var == EPressure) ) return 1;  
	if (var == EVvector) return this->Dimension();
	PZError << "ERROR: " << __PRETTY_FUNCTION__ << " - Variable var " << var << " not found";
	return 0;
}

void TPZIncNavierStokesKEps::Solution(TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol,
									  TPZFMatrix<REAL> &axes, int var, TPZVec<STATE> &Solout){
	
	if (var == EK){
		Solout[0] = Sol[EK];  
		return;
	}//if  
	
	if (var == EEpsilon){
		Solout[0] = Sol[EEpsilon];  
		return;
	}//if    
	
	if (var == EPressure){
		Solout[0] = Sol[EPressure];  
		return;
	}//if
	
	if (var == EVvector){
		int i, n = this->Dimension();   
		for(i = 0; i < n; i++){
			Solout[i] = Sol[EVx + i];
		}//for
		return;   
	}//if
	
	PZError << "Variable " << var << " does not have post processing or does not exist" << std::endl;
	
}

void TPZIncNavierStokesKEps::Contribute(TPZMaterialData &data,
                                        REAL weight,
                                        TPZFMatrix<STATE> &ek,
                                        TPZFMatrix<STATE> &ef){
	
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZVec<STATE> &sol=data.sol[0];
	TPZFMatrix<STATE> &dsol=data.dsol[0];
	
	STATE valor;
	
	const int dim = this->Dimension();
	const int nstate = this->NStateVariables();
	//Getting state variables
	STATE K = sol[EK];
	STATE Eps = sol[EEpsilon];
	STATE Pressure = sol[EPressure];  
	TPZManVector<STATE,3> V(dim);
	int i,j;
	for(i = 0; i < dim; i++) V[i] = sol[EVx+i];
	
	//Getting Grad[state variables]
	TPZManVector<STATE,3> GradK(dim,1);
	for(i = 0; i < dim; i++) GradK[i] = dsol(i, EK);
	TPZManVector<STATE,3> GradEps(dim,1);
	for(i = 0; i < dim; i++) GradEps[i] = dsol(i, EEpsilon);
	TPZManVector<STATE,3> GradPressure(dim,1);
	for(i = 0; i < dim; i++) GradPressure[i] = dsol(i, EPressure);
	TPZFNMatrix<9,STATE> GradV(dim,dim);
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			GradV(i,j) = dsol(j,EVx+i);
		}
	}
	
	//Constants:
	STATE Rt = K*K /(Eps*fMU/fRHO);
	STATE muT = fRHO * fCmu * (K*K/Eps) * exp(-2.5/(1.+Rt/50.));
	
	//TURBULENCE RESIDUALS
	//CONSERVATION OF K                              
	const int nShape = phi.Rows(); 
	TPZFNMatrix<9,STATE> S(dim,dim);
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			S(i,j) = 0.5 * (GradV(i,j) + GradV(j,i) );
		}
	}
	STATE Diss = 2. * fMU * (1. / (4. * K) ) * this->Dot(GradK, GradK);
	
	TPZManVector<STATE,3> GradPhi(dim);
	for(i = 0; i < nShape; i++){
		int k;
		for(k = 0; k < dim; k++) GradPhi[k] = dphi(k,i);
		ef(i*nstate+EK) += -1. * fRHO * this->Dot(V, GradPhi) * K 
		+ (fMU + muT / fSigmaK) * Dot(GradK,GradPhi) 
		-2.0 * muT * this->Dot(S,GradV)*phi[i]
		+ fRHO * Eps * phi[i] 
		-1. * Diss;
		valor = -1. * fRHO * this->Dot(V, GradPhi) * K 
		+ (fMU + muT / fSigmaK) * Dot(GradK,GradPhi) 
		-2.0 * muT * this->Dot(S,GradV)*phi[i]
		+ fRHO * Eps * phi[i] 
		-1. * Diss;                
	}
	
	//CONSERVATION OF EPSILON
	for(i = 0; i < nShape; i++){
		int k;
		for(k = 0; k < dim; k++) GradPhi[k] = dphi(k,i);
		ef(i*nstate+EEpsilon) += -1. * fRHO * this->Dot(V, GradPhi) * Eps
		+ (fMU + muT / fSigmaEps) * this->Dot(GradEps, GradPhi)
		-1. * fCepsilon1 * (Eps/K) * 2. * muT * this->Dot(S,GradV) * phi[i]
		+ fCepsilon2 * (Eps*Eps/K) * phi[i];
		valor = -1. * fRHO * this->Dot(V, GradPhi) * Eps
		+ (fMU + muT / fSigmaEps) * this->Dot(GradEps, GradPhi)
		-1. * fCepsilon1 * (Eps/K) * 2. * muT * this->Dot(S,GradV) * phi[i]
		+ fCepsilon2 * (Eps*Eps/K) * phi[i];                  
	}
	
	//INCOMPRESSIBLE NAVIERS-STOKES RESIDUALS  
	//CONTINUITY EQUATION
	for(i = 0; i < nShape; i++){
		STATE trGradV = 0.;
		int k;
		for(k = 0; k < dim; k++) trGradV += GradV(k,k);
		ef(i*nstate+EPressure) += trGradV * phi[i];
		valor = trGradV * phi[i];
	}  
	
	//CONSERVATION OF LINEAR MOMENTUM
	TPZFNMatrix<9,STATE> T(dim,dim);
	//T = -p I + (mu + muT) * 2 * S
	T = S;
	T *= (fMU + muT ) *2.;
	for(i = 0; i < dim; i++) T(i,i) += -1. * Pressure;
	
	for(j = 0; j < dim; j++){
		for(i = 0; i < nShape; i++){
			int k;
			for(k = 0; k < dim; k++) GradPhi[k] = dphi(k,i);
			valor = fRHO * this->Dot(V, GradV, j) * phi[i]
			-1.  * this->Dot(GradPhi, T, j)
			-1.  * fBodyForce[j] * phi[i];                      
			ef(i*nstate+ EVx+j) += valor;
		}
	}
	
}//method

void TPZIncNavierStokesKEps::Contribute(TPZMaterialData &data,
                                        REAL weight,
                                        TPZFMatrix<STATE> &ef){
	
}


void TPZIncNavierStokesKEps::ContributeBC(TPZMaterialData &data,
                                          REAL weight,
                                          TPZFMatrix<STATE> &ek,
                                          TPZFMatrix<STATE> &ef,
                                          TPZBndCond &bc){
	
}

STATE TPZIncNavierStokesKEps::Dot(TPZFMatrix<STATE> &A, TPZFMatrix<STATE> &B){
	STATE sum = 0.;
	int i, j, rows, cols;
	rows = A.Rows();
	cols = A.Cols();
	for(i = 0; i < rows; i++){
		for(j = 0; j < cols; j++){
			sum += A(i,j) * B(i,j);
		}
	}
	return sum;
}

STATE TPZIncNavierStokesKEps::Dot(TPZVec<STATE> &A, TPZVec<STATE> &B){
	STATE sum = 0.;
	int i, dim;
	dim = A.NElements();
	for(i = 0; i < dim; i++){
		sum += A[i] * B[i];
	}
	return sum;
}

STATE TPZIncNavierStokesKEps::Dot(TPZVec<STATE> &A, TPZFMatrix<STATE> &B, int BRow){
	STATE sum = 0.;
	int i, dim;
	dim = A.NElements();
	for(i = 0; i < dim; i++){
		sum += A[i] * B(BRow,i);
	}
	return sum;
}

int TPZIncNavierStokesKEps::ClassId() const{
    return Hash("TPZIncNavierStokesKEps") ^ TPZMaterial::ClassId() << 1;
}

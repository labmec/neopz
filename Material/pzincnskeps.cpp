/**
 * \file
 * @brief Contains implementations of the TPZIncNavierStokesKEps methods.
 */
#include "pzincnskeps.h"    

TPZIncNavierStokesKEps::TPZIncNavierStokesKEps(int id, int dimension):TPZMaterial(id){
	this->fDimension = dimension;
}

TPZIncNavierStokesKEps::~TPZIncNavierStokesKEps(){}

void TPZIncNavierStokesKEps::SetParameters(REAL MU, REAL RHO, REAL Cmu, REAL SigmaK, REAL SigmaEps, REAL Cepsilon1, REAL Cepsilon2, TPZVec<REAL> &BodyForce ){
	fMU = MU;
	fRHO = RHO;
	fCmu = Cmu;
	fSigmaK = SigmaK;
	fSigmaEps = SigmaEps;
	fCepsilon1 = Cepsilon1;
	fCepsilon2 = Cepsilon2;
	fBodyForce = BodyForce;
}

void TPZIncNavierStokesKEps::GetParameters(REAL &MU, REAL &RHO, REAL &Cmu, REAL &SigmaK, REAL &SigmaEps, REAL &Cepsilon1, REAL &Cepsilon2, TPZVec<REAL> &BodyForce ){
	MU = fMU;
	RHO = fRHO;
	Cmu = fCmu;
	SigmaK = fSigmaK;
	SigmaEps = fSigmaEps;
	Cepsilon1 = fCepsilon1;
	Cepsilon2 = fCepsilon2;
	BodyForce = fBodyForce;
}

int TPZIncNavierStokesKEps::Dimension(){
	return this->fDimension;
}

int TPZIncNavierStokesKEps::NStateVariables(){
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

void TPZIncNavierStokesKEps::Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
									  TPZFMatrix &axes, int var, TPZVec<REAL> &Solout){
	
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
                                        TPZFMatrix &ek,
                                        TPZFMatrix &ef){
	
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	TPZVec<REAL> &sol=data.sol[0];
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	TPZFMatrix &dsol=data.dsol[0];
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	REAL valor;
	
	const int dim = this->Dimension();
	const int nstate = this->NStateVariables();
	//Getting state variables
	REAL K = sol[EK];
	REAL Eps = sol[EEpsilon];
	REAL Pressure = sol[EPressure];  
	TPZManVector<REAL,3> V(dim);
	int i,j;
	for(i = 0; i < dim; i++) V[i] = sol[EVx+i];
	
	//Getting Grad[state variables]
	TPZManVector<REAL,3> GradK(dim,1);
	for(i = 0; i < dim; i++) GradK[i] = dsol(i, EK);
	TPZManVector<REAL,3> GradEps(dim,1);
	for(i = 0; i < dim; i++) GradEps[i] = dsol(i, EEpsilon);
	TPZManVector<REAL,3> GradPressure(dim,1);
	for(i = 0; i < dim; i++) GradPressure[i] = dsol(i, EPressure);
	TPZFNMatrix<9> GradV(dim,dim);
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			GradV(i,j) = dsol(j,EVx+i);
		}
	}
	
	//Constants:
	REAL Rt = K*K /(Eps*fMU/fRHO);
	REAL muT = fRHO * fCmu * (K*K/Eps) * exp(-2.5/(1.+Rt/50.));
	
	//TURBULENCE RESIDUALS
	//CONSERVATION OF K                              
	const int nShape = phi.Rows(); 
	TPZFNMatrix<9> S(dim,dim);
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			S(i,j) = 0.5 * (GradV(i,j) + GradV(j,i) );
		}
	}
	REAL Diss = 2. * fMU * (1. / (4. * K) ) * this->Dot(GradK, GradK);
	
	TPZManVector<REAL,3> GradPhi(dim);
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
		REAL trGradV = 0.;
		int k;
		for(k = 0; k < dim; k++) trGradV += GradV(k,k);
		ef(i*nstate+EPressure) += trGradV * phi[i];
		valor = trGradV * phi[i];
	}  
	
	//CONSERVATION OF LINEAR MOMENTUM
	TPZFNMatrix<9> T(dim,dim);
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
                                        TPZFMatrix &ef){
	
}


void TPZIncNavierStokesKEps::ContributeBC(TPZMaterialData &data,
                                          REAL weight,
                                          TPZFMatrix &ek,
                                          TPZFMatrix &ef,
                                          TPZBndCond &bc){
	
}


REAL TPZIncNavierStokesKEps::Dot(TPZFMatrix &A, TPZFMatrix &B){
	REAL sum = 0.;
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

REAL TPZIncNavierStokesKEps::Dot(TPZVec<REAL> &A, TPZVec<REAL> &B){
	REAL sum = 0.;
	int i, dim;
	dim = A.NElements();
	for(i = 0; i < dim; i++){
		sum += A[i] * B[i];
	}
	return sum;
}

REAL TPZIncNavierStokesKEps::Dot(TPZVec<REAL> &A, TPZFMatrix &B, int BRow){
	REAL sum = 0.;
	int i, dim;
	dim = A.NElements();
	for(i = 0; i < dim; i++){
		sum += A[i] * B(BRow,i);
	}
	return sum;
}



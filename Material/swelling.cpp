/**
 * \file
 * @brief Contains implementations of the TPZSwelling methods.
 */
#include "swelling.h"
#include "pzelmat.h"
#include "pzbndcond.h" 
#include "pzmatrix.h"
#include "pzfmatrix.h"
//#include "pztempmat.h"
#include "pzerror.h"
#include "checkconv.h"
#include <math.h>

using namespace std;

#include <cmath>

REAL TPZSwelling::gFaraday = 96.4853;
REAL TPZSwelling::gVPlus = 2.3;
REAL TPZSwelling::gVMinus = 15.17;
REAL TPZSwelling::gRGas = 8.3145;
REAL TPZSwelling::gTemp = 293.;
REAL TPZSwelling::gMuRef[3] = {0.,0.,0.};

TPZFMatrix<REAL> TPZSwelling::gState;
TPZFMatrix<REAL> TPZSwelling::gphi(1,1,1.);
TPZFMatrix<REAL> TPZSwelling::gdphi(3,1,0.34);

#ifdef _AUTODIFF

void ToMatrix(TPZVec<FADREAL> &vec, TPZFMatrix<REAL> &ek);

#endif

TPZSwelling::TPZSwelling(int matindex, REAL lambda, REAL shear, REAL alfa, REAL M, REAL Gamma, REAL Kperm, REAL DPlus, REAL DMinus,
						 REAL rHinder, REAL Cfc, REAL Nf0, REAL NPlus0, REAL NMinus0) : 
TPZMaterial(matindex) {
	fComputationMode = 0;
	fLambda = lambda;
	fShear = shear;
	fAlfa = alfa;
	fM = M; 
	fGamma = Gamma;
	fKperm.Redim(3,3);
	fDPlus = DPlus;
	fDMinus = DMinus;
	frHinder = rHinder;
	fInitDeform = 1.;
	fCfc = Cfc;
	fNf0 = Nf0;
	fNPlus0 = NPlus0;
	fNMinus0 = NMinus0;
	int ic;
	for(ic=0; ic<3; ic++) {
		fKperm(ic,ic) = Kperm;
	}
	fTheta = 1.;
	fDelt = 0.5;
}

TPZSwelling::~TPZSwelling() {
}


void TPZSwelling::Print(std::ostream &out) {
	TPZMaterial::Print(out);
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "Computation mode : " << fComputationMode << endl;
	out << "Compression modulus " << fLambda << endl;
	out << "Shear modulus " << fShear << endl;
	out << "Biot coupling coeficient " <<  fAlfa << endl;
	out << "Storage modulus " << fM << endl;
	out << "Osmotic coeficient " << fGamma << endl;
	out << "Hydraulic permeability " << fKperm << endl;
	out << "Diffusion coeficient for cations " << fDPlus << endl;
	out << "Diffusion coeficient for anions " << fDMinus << endl;
	out << "Hindrance factor " << frHinder << endl;
	out << "Initial deformation " << fInitDeform << endl;
	out << "Fixed charge density " << fCfc << endl;
	out << "Initial fluid volume fraction " << fNf0 << endl;
	out << "Initial cation volume fraction " << fNPlus0 << endl;
	out << "Initial anion volume fraction " << fNMinus0 << endl;
	out << "Weight factor for time integration of ionic conservation law " << fTheta << endl;
	out << "Timestep " << fDelt << endl;
	out << "Faraday constant " << gFaraday << endl;
	out << "Molar volume cation " <<gVPlus << endl;
	out << "Molar volume anions " << gVMinus << endl;
	out << "Gas constant " << gRGas << endl;
	out << "Absolute temperature " << gTemp << endl;
	int ieq;
	for(ieq=0; ieq<3; ieq++) {
		out << "Reference chemical potentials [" << ieq << "] " << gMuRef[ieq] << endl;
	}
}



/** returns the variable index associated with the name*/
int TPZSwelling::VariableIndex(const std::string &name){
	return TPZMaterial::VariableIndex(name);
}

int TPZSwelling::NSolutionVariables(int var){
	
	return TPZMaterial::NSolutionVariables(var);
}

void TPZSwelling::Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout){
	
	TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}




void TPZSwelling::ContributeBC(TPZMaterialData &data,
                               REAL weight,
                               TPZFMatrix<REAL> &ek,
                               TPZFMatrix<REAL> &ef,
                               TPZBndCond &bc) {
	
	
	// TPZFMatrix<REAL> &dphi = data.dphix;
	// TPZFMatrix<REAL> &dphiL = data.dphixl;
	// TPZFMatrix<REAL> &dphiR = data.dphixr;
	TPZFMatrix<REAL> &phi = data.phi;
	// TPZFMatrix<REAL> &phiL = data.phil;
	// TPZFMatrix<REAL> &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZVec<REAL> &sol=data.sol[0];
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix<REAL> &dsol=data.dsol;
	// TPZFMatrix<REAL> &dsolL=data.dsoll;
	// TPZFMatrix<REAL> &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix<REAL> &daxesdksi=data.daxesdksi;
	// TPZFMatrix<REAL> &axes=data.axes;
	
	if(bc.Material().operator ->() != this){
		PZError << "TPZMatHyperElastic.ContributeBC : this material don't exists \n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "ContributeBC.aplybc, unknown boundary condition type : "<<bc.Type() << endl;
	}
	
	int ndof = NStateVariables();
	int nnod = ek.Rows()/ndof;
	int r = ndof;
	
	int idf,jdf,in,jn;
	switch(bc.Type()){
		case 0:
			for(in=0 ; in<nnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += gBigNumber*phi(in,0)*(bc.Val2()(idf,0)-sol[idf])*weight;
				}
				for(jn=0 ; jn<nnod ; ++jn) {
					for(idf = 0;idf<r;idf++) {
						ek(in*r+idf,jn*r+idf) += gBigNumber*phi(in,0)*phi(jn,0)*weight;
					}
				}
			}
			break;
			
		case 1:
			for(in=0 ; in<nnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					//(ef)(in*r+idf,0) += weight*phi(in,0)*(bc.Val2()(idf,0)-sol[idf]);
					(ef)(in*r+idf,0) += weight*phi(in,0)*(bc.Val2()(idf,0));
				}
			}
			break;
			
		case 2:
			for(in=0 ; in<nnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					for (jdf=0; jdf<r; jdf++){
						(ef)(in*r+idf,0) += phi(in,0)*bc.Val1()(idf,jdf)*(bc.Val2()(jdf,0)-sol[jdf])*weight;
					}
					for(jn=0 ; jn<nnod ; ++jn) {
						for(idf = 0;idf<r;idf++) {
							for(jdf = 0;jdf<r;jdf++) {
								ek(in*r+idf,jn*r+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
							}
						}
					}
				}
				
			}
			
	}//fim switch
}


#ifdef _AUTODIFF


/* The function below makes the correspondence between the dsol vector
 and a matrix ordered F operator
 */
inline int ith(const int i, const int j)
{
	return i*3+j;
}


void TPZSwelling::ContributeElastEnergy(TPZVec<FADFADREAL> &dsol,
										FADFADREAL &U, REAL weight)
{
	FADFADREAL J, TrC; // J = det(F); TrC = Trace(C)
	int numeq = dsol[0].size();
	FADREAL defaultFAD(numeq, REAL(0.), REAL(0.));
	FADFADREAL defaultFADFAD(numeq, defaultFAD, defaultFAD);
	TPZManVector<FADFADREAL,9> GradMap(9);
	int iel;
	for(iel=0; iel<9; iel++) {
		GradMap[iel] = dsol[iel]*FADREAL(fInitDeform);
	}
	
	for(iel=0; iel<3; iel++) GradMap[ith(iel,iel)].val().val() += fInitDeform;
	
	TrC = GradMap[0]*GradMap[0]+GradMap[1]*GradMap[1]+GradMap[2]*GradMap[2]+
    GradMap[3]*GradMap[3]+GradMap[4]*GradMap[4]+GradMap[5]*GradMap[5]+
    GradMap[6]*GradMap[6]+GradMap[7]*GradMap[7]+GradMap[8]*GradMap[8];
	
	J =  GradMap[ith(0,0)] * GradMap[ith(1,1)] * GradMap[ith(2,2)] + 
    GradMap[ith(0,1)] * GradMap[ith(1,2)] * GradMap[ith(2,0)] +
    GradMap[ith(0,2)] * GradMap[ith(1,0)] * GradMap[ith(2,1)] -
    GradMap[ith(0,2)] * GradMap[ith(1,1)] * GradMap[ith(2,0)] -
    GradMap[ith(0,1)] * GradMap[ith(1,0)] * GradMap[ith(2,2)] -
    GradMap[ith(0,0)] * GradMap[ith(1,2)] * GradMap[ith(2,1)]; //  J = det(F)
	
	// expression of the elastic energy
	U += (J*J - FADREAL(1.)) * FADREAL(weight*fLambda/4.) -
    log( J ) * FADREAL(weight*(fLambda/2.+fShear)) +
    (TrC - FADREAL(3.)) * FADREAL(weight*fShear/2.);
	
}

/**
 * Computes the residual vector at an integration point and its tangent matrix by automatic differentiation
 */
void TPZSwelling::ContributeResidual(TPZVec<REAL> & x,
									 TPZVec<FADREAL> & sol, 
									 TPZVec<FADREAL> &dsol,
									 TPZFMatrix<REAL> &phi,
									 TPZFMatrix<REAL> &dphi,
									 TPZVec<FADREAL> &RES,
									 REAL weight){
	const int nstate = 8;
	
	if(fComputationMode == 0) {
		ContributePrevResidual(x,sol,dsol,phi,dphi,RES,weight);
		return;
	}
	int numeq = sol[0].size();
	FADREAL defaultFAD(numeq, REAL(0.), REAL(0.));
	FADFADREAL defaultFADFAD(numeq, defaultFAD, defaultFAD);
	FADFADREAL deform(defaultFADFAD);
	TPZManVector<FADFADREAL> dsolFADFAD(9,defaultFADFAD);
	
	// extend the gradient of the solution to a variable which contains second derivatives
	int ieq,der;
	for(der=0; der<9; der++) {
		dsolFADFAD[der].val().val() = dsol[der].val();
		for(ieq=0; ieq<numeq; ieq++) {
			dsolFADFAD[der].val().fastAccessDx(ieq) = dsol[der].fastAccessDx(ieq);
			dsolFADFAD[der].fastAccessDx(ieq).val() = dsol[der].fastAccessDx(ieq);
		}
	}
	// The first derivative of the deform variable contains the residual vector
	// The second derivative contains the tangent matrix
	ContributeElastEnergy(dsolFADFAD, deform, weight);
	
	int jeq;
	for(ieq=0; ieq<numeq; ieq++) {
		RES[ieq].val() += deform.fastAccessDx(ieq).val();
		for(jeq=0; jeq<numeq; jeq++) {
			RES[ieq].fastAccessDx(jeq) += deform.fastAccessDx(ieq).fastAccessDx(jeq);
		}
	}
	
	//return;
	//  cout << "dsol " << dsol;
	
	// compute the gradient of the map induced by the displacement of the element
	// include the derivative of the map with respect to the solution
	FADREAL GradMap[3][3];
	GradMap[0][0] = dsol[0]+REAL(1.);
	GradMap[0][1] = dsol[1];
	GradMap[0][2] = dsol[2];
	GradMap[1][0] = dsol[3];
	GradMap[1][1] = dsol[4]+REAL(1.);
	GradMap[1][2] = dsol[5];
	GradMap[2][0] = dsol[6];
	GradMap[2][1] = dsol[7];
	GradMap[2][2] = dsol[8]+REAL(1.);
	
	// compute the determinant of the map
	// this computation will carry the derivative w.r.t the solution
	FADREAL J = GradMap[0][0] * GradMap[1][1] * GradMap[2][2] + 
    GradMap[0][1] * GradMap[1][2] * GradMap[2][0] +
    GradMap[0][2] * GradMap[1][0] * GradMap[2][1] -
    GradMap[0][2] * GradMap[1][1] * GradMap[2][0] -
    GradMap[0][1] * GradMap[1][0] * GradMap[2][2] -
    GradMap[0][0] * GradMap[1][2] * GradMap[2][1]; //  J = det(F)
	
	// compute the inverse of the map, including its derivatives
	FADREAL GradMapInv[3][3];
	for(ieq=0; ieq<3; ieq++) {
		int ieqp = (ieq+1)%3;
		int ieqpp = (ieq+2)%3;
		for(jeq=0; jeq<3; jeq++) {
			int jeqp = (jeq+1)%3;
			int jeqpp = (jeq+2)%3;
			GradMapInv[ieq][jeq] = (GradMap[ieqp][jeqp]*GradMap[ieqpp][jeqpp]-GradMap[ieqpp][jeqp]*GradMap[ieqp][jeqpp])/J;
			//      cout << ieq << ' ' << jeq << endl << ieqp << ' ' << jeqp << '+' << ieqpp << ' ' <<  jeqpp << endl 
			//	   << ieqpp << ' ' << jeqp << '-' << ieqp << ' ' << jeqpp <<  endl << GradMapInv[ieq][jeq];
		}
	}
	
	// compute the Lagrangian volume fractions in function of mu, pres and ksi
	TPZVec<FADREAL> N(3);
	ComputeN(sol,N);
	
	// uncomment to observe that the numerical procedure gives the same result as the analytic procedure
	/*
	 TPZVec<REAL> N2(3);
	 NResidual(sol,N2);
	 
	 for(ieq=0; ieq<3; ieq++) {
	 N2[ieq] -= N[ieq];
	 }
	 
	 cout << "the difference between both approaches\n" << N2;
	 */
	
	int ishape,nshape = phi.Rows();
	for(ishape=0; ishape<nshape; ishape++) {
		// compute the gradient of the weighting function in the deformed configuration
		FADREAL GradPhi[3];
		for(ieq=0; ieq<3; ieq++) {
			GradPhi[ieq] = GradMapInv[0][ieq]*dphi(0,ishape)+GradMapInv[1][ieq]*dphi(1,ishape)+GradMapInv[2][ieq]*dphi(2,ishape);
			// Add the contribution of the pressure to the linear momentum equation
			RES[ishape*nstate+ieq] += sol[3]*(GradPhi[ieq]*weight)*J;
		}
		// Add the contribution of the rate of change of the volume to the mass balance of the mixture equation
		RES[ishape*nstate+3]+= J*(phi(ishape,0)*weight);
		// Add the contribution of the Lagrange volume fractions to the mass balance of the mixture equation
		RES[ishape*nstate+3] -= (N[0]+N[1]+N[2])*(phi(ishape,0)*weight);
		// Add the contribution to electro neutrality equation
		RES[ishape*nstate+7] -= (N[1]/gVPlus-N[2]/gVMinus)*(phi(ishape,0)*gFaraday*weight);
		for(ieq=0; ieq<3; ieq++) {
			// Add the contribution to the mass balance of the constituents
			RES[ishape*nstate+4+ieq] += (N[ieq]*phi(ishape,0)*weight);
			REAL KGradPhi = 0.;
			for(jeq=0; jeq<3; jeq++) {
				KGradPhi += dphi(jeq,ishape)*fKperm(ieq,jeq)*fTheta*fDelt*weight;
			}
			// Add the contribution of the difusion of the constituents
			RES[ishape*nstate+4+ieq] += KGradPhi*dsol[ieq+3*(4+ieq)];
		}
	}
}

void TPZSwelling::ContributePrevResidual(TPZVec<REAL> & x,
										 TPZVec<FADREAL> & sol,
										 TPZVec<FADREAL> &dsol,
										 TPZFMatrix<REAL> &phi,
										 TPZFMatrix<REAL> &dphi,
										 TPZVec<FADREAL> &RES,
										 REAL weight){
	const int nstate = 8;
	
	if(fComputationMode != 0) {
		cout << "TPZSwelling::ContributePrevResidual should not be called\n";
		return;
	}
	int numeq = sol[0].size();
	FADREAL defaultFAD(numeq, REAL(0.), REAL(0.));
	FADFADREAL defaultFADFAD(numeq, defaultFAD, defaultFAD);
	FADFADREAL deform(defaultFADFAD);
	TPZManVector<FADFADREAL> dsolFADFAD(9,defaultFADFAD);
	
	// extend the gradient of the solution to a variable which contains second derivatives
	int ieq,der;
	for(der=0; der<9; der++) {
		dsolFADFAD[der].val().val() = dsol[der].val();
		for(ieq=0; ieq<numeq; ieq++) {
			dsolFADFAD[der].val().fastAccessDx(ieq) = dsol[der].fastAccessDx(ieq);
			dsolFADFAD[der].fastAccessDx(ieq).val() = dsol[der].fastAccessDx(ieq);
		}
	}
	// The first derivative of the deform variable contains the residual vector
	// The second derivative contains the tangent matrix
	
	int jeq;
	
	//return;
	//  cout << "dsol " << dsol;
	
	// compute the gradient of the map induced by the displacement of the element
	// include the derivative of the map with respect to the solution
	REAL GradMap[3][3];
	GradMap[0][0] = dsol[0].val()+1.;
	GradMap[0][1] = dsol[1].val();
	GradMap[0][2] = dsol[2].val();
	GradMap[1][0] = dsol[3].val();
	GradMap[1][1] = dsol[4].val()+1.;
	GradMap[1][2] = dsol[5].val();
	GradMap[2][0] = dsol[6].val();
	GradMap[2][1] = dsol[7].val();
	GradMap[2][2] = dsol[8].val()+1.;
	
	// compute the determinant of the map
	// this computation will carry the derivative w.r.t the solution
	REAL J = GradMap[0][0] * GradMap[1][1] * GradMap[2][2] +
    GradMap[0][1] * GradMap[1][2] * GradMap[2][0] +
    GradMap[0][2] * GradMap[1][0] * GradMap[2][1] -
    GradMap[0][2] * GradMap[1][1] * GradMap[2][0] -
    GradMap[0][1] * GradMap[1][0] * GradMap[2][2] -
    GradMap[0][0] * GradMap[1][2] * GradMap[2][1]; //  J = det(F)
	
	// compute the inverse of the map, including its derivatives
	/*
	 REAL GradMapInv[3][3];
	 for(ieq=0; ieq<3; ieq++) {
	 int ieqp = (ieq+1)%3;
	 int ieqpp = (ieq+2)%3;
	 for(jeq=0; jeq<3; jeq++) {
	 int jeqp = (jeq+1)%3;
	 int jeqpp = (jeq+2)%3;
	 GradMapInv[ieq][jeq] = (GradMap[ieqp][jeqp]*GradMap[ieqpp][jeqpp]-GradMap[ieqpp][jeqp]*GradMap[ieqp][jeqpp])/J;
	 //      cout << ieq << ' ' << jeq << endl << ieqp << ' ' << jeqp << '+' << ieqpp << ' ' <<  jeqpp << endl
	 //	   << ieqpp << ' ' << jeqp << '-' << ieqp << ' ' << jeqpp <<  endl << GradMapInv[ieq][jeq];
	 }
	 }
	 */
	// compute the Lagrangian volume fractions in function of mu, pres and ksi
	TPZVec<REAL> N(3);
	ComputeN(sol,N);
	
	// uncomment to observe that the numerical procedure gives the same result as the analytic procedure
	/*
	 TPZVec<REAL> N2(3);
	 NResidual(sol,N2);
	 
	 for(ieq=0; ieq<3; ieq++) {
	 N2[ieq] -= N[ieq];
	 }
	 
	 cout << "the difference between both approaches\n" << N2;
	 */
	
	int ishape,nshape = phi.Rows();
	for(ishape=0; ishape<nshape; ishape++) {
		// Add the contribution of the rate of change of the volume to the mass balance of the mixture equation
		RES[ishape*nstate+3] -= J*(phi(ishape,0)*weight);
		// Add the contribution of the Lagrange volume fractions to the mass balance of the mixture equation
		RES[ishape*nstate+3] += (N[0]+N[1]+N[2])*(phi(ishape,0)*weight);
		// Add the contribution to electro neutrality equation
		RES[ishape*nstate+7] += (N[1]/gVPlus-N[2]/gVMinus)*(phi(ishape,0)*gFaraday*weight);
		for(ieq=0; ieq<3; ieq++) {
			// Add the contribution to the mass balance of the constituents
			RES[ishape*nstate+4+ieq] += (N[ieq]*phi(ishape,0)*weight);
			REAL KGradPhi = 0.;
			for(jeq=0; jeq<3; jeq++) {
				KGradPhi += dphi(jeq,ishape)*fKperm(ieq,jeq)*(1.-fTheta)*fDelt*weight;
			}
			// Add the contribution of the difusion of the constituents
			RES[ishape*nstate+4+ieq] += KGradPhi*dsol[ieq+3*(4+ieq)];
		}
	}
}



void TPZSwelling::ComputeW(FADFADREAL &W, TPZVec<REAL> &N) {
	FADREAL defaultFAD(3,REAL(0.),REAL(0.));
	FADFADREAL defaultFADFAD(3,defaultFAD,defaultFAD);
	W = defaultFADFAD;
	FADFADREAL NFAD[3] = {defaultFADFAD,defaultFADFAD,defaultFADFAD};
	
	// transform the N (volume fractions) in second order derivatives
	int ieq;
	for(ieq=0; ieq<3; ieq++) {
		NFAD[ieq].val().val() = N[ieq];
		NFAD[ieq].val().fastAccessDx(ieq) = 1.;
		NFAD[ieq].fastAccessDx(ieq).val() = 1.;
	}
	W += FADREAL(gMuRef[0])*NFAD[0]+FADREAL(gMuRef[1])*NFAD[1]+
    FADREAL(gMuRef[2])*NFAD[2]-
    FADREAL(gRGas*gTemp*fGamma)*(NFAD[1]*FADREAL(1./gVPlus)+NFAD[2]*FADREAL(1./gVMinus))*log(NFAD[0])+
	FADREAL(gRGas*gTemp/gVPlus)*NFAD[1]*(log(NFAD[1]/FADREAL(gVPlus))-FADREAL(1.))+
	FADREAL(gRGas*gTemp/gVMinus)*NFAD[2]*(log(NFAD[2]/FADREAL(gVMinus))-FADREAL(1.));
}

void TPZSwelling::ComputeN(TPZVec<REAL> &mu, REAL ksi, REAL pressure, TPZVec<REAL> &N) {
	
	TPZFMatrix<REAL> res,tangent;
	ExactSolution(mu,ksi,pressure,N);
	NResidual(mu,ksi,pressure,N,res,tangent);
	// compute the residual of the current N values
	REAL resnorm = Norm(res);
	int ieq;
	while(resnorm > 1.e-6) {
		// Compute the correction as suggested by the Newton method
		// This only works if the residual is suficiently small
		tangent.SolveDirect(res,ELU);
		for(ieq=0; ieq<3; ieq++) {
			N[ieq] -= res(ieq,0);
		}
		NResidual(mu,ksi,pressure,N,res,tangent);
		resnorm = Norm(res);
	}
}

void TPZSwelling::NResidual(TPZVec<REAL> &mu, REAL ksi, REAL pressure, TPZVec<REAL> &N, TPZFMatrix<REAL> &res, TPZFMatrix<REAL> &tangent) {
	
	FADREAL defaultFAD(3,REAL(0.),REAL(0.));
	FADFADREAL defaultFADFAD(3,defaultFAD,defaultFAD),W;
	W = defaultFADFAD;
	
	// Compute the mixing energy as a second order derivative variable
	// the gradient contributes to the residual
	// the second derivatives contribute to the tangent matrix
	ComputeW(W,N);
	
	//alternative direct formula
	//  REAL gradw[3];
	//  gradw[0] = gRGas*gTemp*fGamma*(N[1]/gVPlus + N[2]/gVMinus)/N[0];
	//  gradw[1] = (gRGas*gTemp/gVPlus)*log(N[1]/(pow(N[0],fGamma)*gVPlus));
	//  gradw[2] = (gRGas*gTemp/gVMinus)*log(N[2]/(pow(N[0],fGamma)*gVMinus));
	
	res.Redim(3,1);
	tangent.Redim(3,3);
	REAL z[3] = {0.,1.,-1.};
	REAL v[3] = {1.,gVPlus,gVMinus};
	int ieq,jeq;
	for(ieq=0; ieq<3; ieq++) {
		res(ieq,0) = -mu[ieq]+W.val().fastAccessDx(ieq)+gFaraday*z[ieq]*ksi/v[ieq]+pressure;
		for(jeq=0; jeq<3; jeq++) {
			tangent(ieq,jeq) = W.fastAccessDx(ieq).fastAccessDx(jeq);
		}
	}
}

void TPZSwelling::NResidual(TPZVec<FADREAL> &sol, TPZVec<FADREAL> &N) {
	
	// compute the values of N numerically and its derivatives by applying a single newton iteration using FAD variables
	TPZManVector<REAL,3> muloc(3),Nloc(3);
	// create the residual vector and tangent matrix for solving for N
	TPZVec<FADREAL> res(3);
	TPZVec<TPZVec<FADREAL> > tangent(3,res);
	
	FADFADREAL W;
	
	int ieq;
	REAL z[3] = {0.,1.,-1.};
	REAL v[3] = {1.,gVPlus,gVMinus};
	FADREAL &ksi = sol[7];
	FADREAL &pressure = sol[3];
	for(ieq=0; ieq<3; ieq++) {
		muloc[ieq]= sol[ieq+4].val();
		Nloc[ieq] = N[ieq].val();
	}
	// compute the value of N by any method (without carrying the derivatives)
	ExactSolution(muloc,ksi.val(),pressure.val(),Nloc);
	ComputeW(W,Nloc);
	
	// initialize the N variable with the correct solution but with zero derivatives
	for(ieq=0; ieq<3; ieq++) {
		N[ieq] = FADREAL(sol[0].size(),Nloc[ieq],0.);
	}
	// perform a single Newton iteration
	// build the tangent matrix and residual (including the derivatives)
	int jeq;
	for(ieq=0; ieq<3; ieq++) {
		res[ieq] = -sol[ieq+4]+W.val().fastAccessDx(ieq)+gFaraday*z[ieq]*ksi/v[ieq]+pressure;
		for(jeq=0; jeq<3; jeq++) {
			tangent[ieq][jeq] = FADREAL(sol[0].size(),W.fastAccessDx(ieq).fastAccessDx(jeq),0);
		}
	}
	// update by a single Newton iteration
	Solve(tangent,res);
	for(ieq=0; ieq<3; ieq++) {
		N[ieq] -= res[ieq];
	}
}

void TPZSwelling::Solve(TPZVec<TPZVec<FADREAL> > &tangent, TPZVec<FADREAL> &res) {
	// implements a simple LU decomposition using FADREAL variables
	int neq = res.NElements();
	int k,i,j;
	for(k=0; k<neq; k++) {
		FADREAL pivot = tangent[k][k];
		for ( i = k+1; i < neq; i++ ) {
			tangent[i][k] /= pivot;
			for (j = k+1; j < neq; j++ ) tangent[i][j] -= tangent[i][k]*tangent[k][j];
			res[i] -= tangent[i][k]*res[k];
		}
	}
	for ( i = neq-1; i >= 0; i-- ) {
		for (j = i+1; j < neq ; j++ ) {
			res[i] -= tangent[i][j]*res[j];
		}
		res[i] /= tangent[i][i];
	}
}

int TPZSwelling::main() {
	
	// program created for testing purposes
	int matindex = 1;
	REAL lambda = 1;
	REAL shear = 0.5;
	REAL alfa = 1.;
	REAL M = 0.2;
	REAL Gamma = 0.9;
	REAL Kperm = 1.e-4;
	REAL DPlus = 5.e-4;
	REAL DMinus = 5.e-4;
	REAL rHinder = 0.4;
	REAL Cfc = -2.e-4;
	REAL Nf0 = 0.8;
	// initiele externe zoutconcentratie
	REAL C0 = 0.15e-3;
	REAL cplus0 = (-Cfc+sqrt(Cfc*Cfc+4.*C0*C0))/2.;
	REAL cminus0 = (Cfc+sqrt(Cfc*Cfc+4.*C0*C0))/2.;
	REAL NPlus0 = gVPlus*cplus0*Nf0;
	REAL NMinus0 = gVMinus*cminus0*Nf0;
	
	TPZSwelling test(matindex, lambda, shear, alfa, M, Gamma, Kperm, DPlus, DMinus,
					 rHinder, Cfc, Nf0, NPlus0, NMinus0);
	
	test.Print(cout);
	
	// initialize a set of variables as they would be setup by the element
	int ieq,d;
	TPZFMatrix<REAL> phi(1,1),dphi(3,1);
	phi(0,0) = 1.;
	for(ieq=0; ieq<3; ieq++) {
		dphi(ieq,0) = 0.34;
	}
	FADREAL defaultFAD(8,REAL(0.),REAL(0.));
	TPZVec<FADREAL> sol(8,defaultFAD),dsol(24,defaultFAD);
	TPZVec<REAL> x(3,1.);
	REAL weight = 0.33;
	TPZVec<FADREAL> RES(8,defaultFAD);
	REAL values[8] = {1.,1.,1.,0.001,-0.7243,-9258.,-1401.,-15};
	
	TPZVec<REAL> mu(3),N(3,0.);
	for(ieq=0; ieq<3; ieq++) mu[ieq] = values[ieq+4];
	REAL ksi,J = 1.0062,pres;
	test.ComputeInitialGuess(mu,J,pres,ksi,N);
	values[3] = pres;
	values[7] = ksi;
	cout << "Based on J " << J << " and mu \n" << mu << "\nThe initial guesses pres " << pres << " ksi " << ksi << "\n N\n" << N << endl;
	
	
	for(ieq=0; ieq<8; ieq++) {
		sol[ieq].val() = values[ieq];
		sol[ieq].fastAccessDx(ieq) = phi(0);
		for(d=0; d<3; d++) {
			dsol[ieq*3 + d].val() = dphi(d,0)*values[ieq];
			dsol[ieq*3 + d].fastAccessDx(ieq) = dphi(d,0);
		}
	}
	
	// procedure to check whether the stiffness matrix is tangent to the residual vector
	
	REAL rangeval[11] = {0.1,0.1,0.1,0.0001,0.01,1.,1.,0.1,0.,0.,0.};
	TPZFMatrix<REAL> state(11,1),range(11,1,0.0);
	for(ieq=0; ieq<8; ieq++) {
		state(ieq,0) = values[ieq];
		range(ieq,0) = rangeval[ieq];
	}
	for(ieq=8; ieq<11; ieq++) {
		state(ieq,0) = N[ieq-8];
		range(ieq,0) = rangeval[ieq];
	}
	TPZVec<REAL> coefs(1,1.);
	CheckConvergence(test,state,range,coefs);
	
	
	// sample computation of a contribution to a stiffness matrix
	// this procedure was built to check whether the stiffness matrix is symetric
	
	test.ContributeResidual(x,sol,dsol,phi,dphi,RES,weight);
	TPZFMatrix<REAL> ek;
	ToMatrix(RES,ek);
	ek.Print("stiffness matrix");
	
	cout << RES;
	return 0;
	
}

void ToMatrix(TPZVec<FADREAL> &vec, TPZFMatrix<REAL> &ek) {
	int nel = vec.NElements();
	ek.Redim(nel,nel);
	int ieq,jeq;
	for(ieq=0; ieq<nel; ieq++) {
		for(jeq=0; jeq<nel; jeq++) {
			ek(ieq,jeq) = vec[ieq].fastAccessDx(jeq);
		}
	}
}

/*
 Methods to implement automatic conformity checks
 
 */
int TPZSwelling::NumCases() {
	return 1;
}

void TPZSwelling::LoadState(TPZFMatrix<REAL> &state) {
	if(state.Rows() != 11) {
		cout << "TPZSwelling LoadState wrong number of variables expect bad things\n";
	}
	gState = state;
}
void TPZSwelling::ComputeTangent(TPZFMatrix<REAL> &tangent,TPZVec<REAL> &coefs, int cases) {
	tangent.Redim(11,11);
	int ieq,jeq,d;
	FADREAL defaultFAD(8,REAL(0.),REAL(0.));
	TPZVec<FADREAL> sol(8,defaultFAD),dsol(24,defaultFAD);
	TPZVec<REAL> x(3,1.);
	REAL weight = 0.33;
	TPZVec<FADREAL> RES(8,defaultFAD);
	for(ieq=0; ieq<8; ieq++) {
		sol[ieq].val() = gphi(0,0)*gState(ieq,0);
		sol[ieq].fastAccessDx(ieq) = gphi(0);
		for(d=0; d<3; d++) {
			dsol[ieq*3 + d].val() = gdphi(d,0)*gState(ieq,0);
			dsol[ieq*3 + d].fastAccessDx(ieq) = gdphi(d,0);
		}
	}
	ContributeResidual(x,sol,dsol,gphi,gdphi,RES,weight);
	for(ieq=0; ieq<8; ieq++) {
		for(jeq=0; jeq<8; jeq++) {
			tangent(ieq,jeq) = RES[ieq].fastAccessDx(jeq);
		}
	}
	TPZVec<REAL> N(3,0.);
	for(ieq=0; ieq<3; ieq++) N[ieq] = gState(8+ieq,0);
	FADREAL defwFAD(3,REAL(0.),REAL(0.));
	FADFADREAL defwFADFAD(3,defwFAD,defwFAD);
	FADFADREAL W(defwFADFAD);
	ComputeW(W,N);
	for(ieq=0; ieq<3; ieq++) {
		for(jeq=0; jeq<3; jeq++) {
			tangent(8+ieq,8+jeq) = W.fastAccessDx(ieq).fastAccessDx(jeq);
		}
	}
}
void TPZSwelling::Residual(TPZFMatrix<REAL> &res, int cases) {
	res.Redim(11,1);
	int ieq,d;
	FADREAL defaultFAD(8,REAL(0.),REAL(0.));
	TPZVec<FADREAL> sol(8,defaultFAD),dsol(24,defaultFAD);
	TPZVec<REAL> x(3,1.);
	REAL weight = 0.33;
	TPZVec<FADREAL> RES(8,defaultFAD);
	for(ieq=0; ieq<8; ieq++) {
		sol[ieq].val() = gphi(0,0)*gState(ieq,0);
		sol[ieq].fastAccessDx(ieq) = gphi(0);
		for(d=0; d<3; d++) {
			dsol[ieq*3 + d].val() = gdphi(d,0)*gState(ieq,0);
			dsol[ieq*3 + d].fastAccessDx(ieq) = gdphi(d,0);
		}
	}
	ContributeResidual(x,sol,dsol,gphi,gdphi,RES,weight);
	for(ieq=0; ieq<8; ieq++) {
		res(ieq,0) = RES[ieq].val();
	}
	TPZVec<REAL> N(3,0.);
	for(ieq=0; ieq<3; ieq++) N[ieq] = gState(8+ieq,0);
	FADREAL defwFAD(3,REAL(0.),REAL(0.));
	FADFADREAL defwFADFAD(3,defwFAD,defwFAD);
	FADFADREAL W(defwFADFAD);
	ComputeW(W,N);
	for(ieq=0; ieq<3; ieq++) {
		res(8+ieq,0) = W.fastAccessDx(ieq).val();
	}
}

#endif

void TPZSwelling::ComputeInitialGuess(TPZVec<REAL> &mu, REAL J, REAL &pres, REAL &ksi, TPZVec<REAL> &N) {
	
	pres = 0.;
	REAL expcontr = exp((gVPlus*(mu[1]-gMuRef[1]-pres)+gVMinus*(mu[2]-gMuRef[2]-pres))/(gRGas*gTemp));
	REAL expcontr2 = expcontr*REAL(4.);
	REAL sqrtval = sqrt(fCfc*fCfc+expcontr2);
	REAL cplus = (-fCfc + sqrtval)/2.;
	REAL cmin = (fCfc+sqrtval)/2.;
	
	int iter;
	for(iter=0; iter<5; iter++) {
		pres = (mu[0]-gMuRef[0])+fGamma*gRGas*gTemp*(cplus+cmin);
		cout.precision(12);
		cout << "pres " << pres << " cplus " << cplus << " cmin " << cmin << " expcontr " << expcontr << endl;
		REAL expcontr = exp((gVPlus*(mu[1]-gMuRef[1]-pres)+gVMinus*(mu[2]-gMuRef[2]-pres))/(gRGas*gTemp));
		REAL expcontr2 = expcontr*REAL(4.);
		sqrtval = sqrt(fCfc*fCfc+expcontr2);
		cplus = (-fCfc + sqrtval)/2.;
		cmin = (fCfc+sqrtval)/2.;
	}
    cout.precision(12);
	cout << "cplus " << cplus << " cmin " << cmin << endl;
	ksi= (gRGas*gTemp*log(cmin/cplus)+gVPlus*(mu[1]-gMuRef[1]-pres)-gVMinus*(mu[2]-gMuRef[2]-pres))/(2.*gFaraday);
	N[0] = J-1+fNf0;
	N[1] = N[0]*gVPlus*cplus;
	N[2] = N[0]*gVMinus*cmin;
	
}

void TPZSwelling::ExactSolution(TPZVec<REAL> &mu, REAL ksi, REAL pres, TPZVec<REAL> &N) {
	
	// computes the analytic solution of N in function of mu, ksi and pressure
	/*
	 REAL test[3];
	 test[0] = pres-mu[0]-gRGas*gTemp*fGamma*(N[1]/(N[0]*gVPlus)+N[2]/(N[0]*gVMinus));
	 test[1] = pres-mu[1]+gFaraday*ksi/gVPlus+gRGas*gTemp*log(N[1]/(pow(N[0],fGamma)*gVPlus))/gVPlus;
	 test[2] = pres-mu[2]-gFaraday*ksi/gVMinus+gRGas*gTemp*log(N[2]/(pow(N[0],fGamma)*gVMinus))/gVMinus;
	 */
	
	REAL c1,c2;
	c1 = (mu[1]-gFaraday*ksi/gVPlus - pres)*gVPlus/(gRGas*gTemp);
	c2 = (mu[2]+gFaraday*ksi/gVMinus - pres)*gVMinus/(gRGas*gTemp);
	
	//  REAL test1 = pow(N[0],fGamma)*exp(c1)*gVPlus - N[1];
	//  REAL test2 = pow(N[0],fGamma)*exp(c2)*gVMinus - N[2];
	
	REAL fac = exp(c1)+exp(c2);
	REAL fac2 = (-mu[0]+pres)/(gRGas*gTemp*fGamma);
	N[0] = pow(fac2/fac,1./(fGamma-1));
	
	N[1] = pow(N[0],fGamma)*exp(c1)*gVPlus;
	N[2] = pow(N[0],fGamma)*exp(c2)*gVMinus;
	
}

#ifdef _AUTODIFF

void TPZSwelling::ComputeN(TPZVec<FADREAL> &sol, TPZVec<FADREAL> &N) {
	
	// computes the analytic solution of N carrying the derivatives
	
	FADREAL &pres = sol[3];
	FADREAL &mu0 = sol[4];
	FADREAL &mu1 = sol[5];
	FADREAL &mu2 = sol[6];
	FADREAL &ksi = sol[7];
	/*
	 REAL test[3];
	 test[0] = pres-mu[0]-gRGas*gTemp*fGamma*(N[1]/(N[0]*gVPlus)+N[2]/(N[0]*gVMinus));
	 test[1] = pres-mu[1]+gFaraday*ksi/gVPlus+gRGas*gTemp*log(N[1]/(pow(N[0],fGamma)*gVPlus))/gVPlus;
	 test[2] = pres-mu[2]-gFaraday*ksi/gVMinus+gRGas*gTemp*log(N[2]/(pow(N[0],fGamma)*gVMinus))/gVMinus;
	 */
	
	FADREAL expc1,expc2;
	expc1 = exp((mu1-gFaraday*ksi/gVPlus - pres)*gVPlus/(gRGas*gTemp));
	
	expc2 = exp((mu2+gFaraday*ksi/gVMinus - pres)*gVMinus/(gRGas*gTemp));
	
	//  REAL test1 = pow(N[0],fGamma)*exp(c1)*gVPlus - N[1];
	//  REAL test2 = pow(N[0],fGamma)*exp(c2)*gVMinus - N[2];
	
	
	FADREAL fac2 = (-mu0+pres)/(gRGas*gTemp*fGamma);
	N[0] = pow(fac2/(expc1+expc2),1./(fGamma-1));
	FADREAL N0Gamma = pow(N[0],fGamma);
	
	N[1] = N0Gamma*expc1*gVPlus;
	N[2] = N0Gamma*expc2*gVMinus;
}

void TPZSwelling::ComputeN(TPZVec<FADREAL> &sol, TPZVec<REAL> &N) {
	
	// computes the analytic solution of N carrying the derivatives
	
	REAL &pres = sol[3].val();
	REAL &mu0 = sol[4].val();
	REAL &mu1 = sol[5].val();
	REAL &mu2 = sol[6].val();
	REAL &ksi = sol[7].val();
	/*
	 REAL test[3];
	 test[0] = pres-mu[0]-gRGas*gTemp*fGamma*(N[1]/(N[0]*gVPlus)+N[2]/(N[0]*gVMinus));
	 test[1] = pres-mu[1]+gFaraday*ksi/gVPlus+gRGas*gTemp*log(N[1]/(pow(N[0],fGamma)*gVPlus))/gVPlus;
	 test[2] = pres-mu[2]-gFaraday*ksi/gVMinus+gRGas*gTemp*log(N[2]/(pow(N[0],fGamma)*gVMinus))/gVMinus;
	 */
	
	REAL expc1,expc2;
	expc1 = exp((mu1-gFaraday*ksi/gVPlus - pres)*gVPlus/(gRGas*gTemp));
	
	expc2 = exp((mu2+gFaraday*ksi/gVMinus - pres)*gVMinus/(gRGas*gTemp));
	
	//  REAL test1 = pow(N[0],fGamma)*exp(c1)*gVPlus - N[1];
	//  REAL test2 = pow(N[0],fGamma)*exp(c2)*gVMinus - N[2];
	
	
	REAL fac2 = (-mu0+pres)/(gRGas*gTemp*fGamma);
	N[0] = pow(fac2/(expc1+expc2),(float)1./(fGamma-1));
	REAL N0Gamma = pow(N[0],fGamma);
	
	N[1] = N0Gamma*expc1*gVPlus;
	N[2] = N0Gamma*expc2*gVMinus;
}

#endif

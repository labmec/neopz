/**
 * @file
 * @brief Contains implementations of the TPZSpaceTimeRichardsEq methods.
 */

#include "pzspacetimerichardseq.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzaxestools.h"

#include <fstream>
#include <cmath>
using namespace std;

/** @brief Inicializing local variable TCoeff */
STATE TPZSpaceTimeRichardsEq::TCoeff = 1./60.;
/** @brief Inicializing local variable LCoeff */
STATE TPZSpaceTimeRichardsEq::LCoeff = 1000.;
/** @brief Inicializing loval variable deltaDerivada */
STATE TPZSpaceTimeRichardsEq::deltaDerivada = 1.e-3;

TPZSpaceTimeRichardsEq::TPZSpaceTimeRichardsEq():
TPZRegisterClassId(&TPZSpaceTimeRichardsEq::ClassId),
TPZMaterial()
{
}

TPZSpaceTimeRichardsEq::TPZSpaceTimeRichardsEq(int id): 
TPZRegisterClassId(&TPZSpaceTimeRichardsEq::ClassId),
TPZMaterial(id)
{
}

TPZSpaceTimeRichardsEq::TPZSpaceTimeRichardsEq(int matid, STATE Alpha, STATE N, STATE ThetaS, STATE ThetaR, STATE Ks) 
: TPZRegisterClassId(&TPZSpaceTimeRichardsEq::ClassId),
TPZMaterial(matid){
	this->Set(Alpha, N, ThetaS, ThetaR, Ks);
}

TPZSpaceTimeRichardsEq::~TPZSpaceTimeRichardsEq()
{
}

void TPZSpaceTimeRichardsEq::Set(STATE Alpha, STATE N, STATE ThetaS, STATE ThetaR, STATE Ks){
	this->fAlpha = Alpha;
	this->fN = N;
	this->fThetaS = ThetaS;
	this->fThetaR = ThetaR;
	this->fKs = Ks;
}

int TPZSpaceTimeRichardsEq::Dimension() const {
	return 2;
}

int TPZSpaceTimeRichardsEq::NStateVariables() const{
	return 1;
}

void TPZSpaceTimeRichardsEq::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &dphi = data.dphix;
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	const STATE sol = data.sol[0][0];
	
	const STATE BetaBarT = 0*LCoeff*data.detjac/2.; //beta=(0,1)
	
	TPZFNMatrix<2,STATE> dsol(2,1,0.);
    /*
    TPZFNMatrix<9,STATE> axes(data.axes.Rows(),data.axes.Cols());
    for (int r=0; r<axes.Rows(); r++) {
        for (int c=0; c<axes.Cols(); c++) {
            axes(r,c) = data.axes(r,c);
        }
    }
     */
	TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsol, data.axes);
	
	const int phr = phi.Rows();
	int i, j;
	
	for(i = 0; i < phr; i++){
		const STATE BetaBarGradV = BetaBarT*dphi(1,i);
		ef(i,0) += -1.*weight *( -1.*sol*dphi(1,i) +1.*dsol(1,0)*BetaBarGradV + /*(K/C)*/(dsol(0,0))*dphi(0,i));
		for(j = 0; j < phr; j++){
			ek(i,j) += weight * ( -1.*phi(j,0)*dphi(1,i)+dphi(1,j)*BetaBarGradV + dphi(0,i)*dphi(0,j) );
		}
	}//for i
}//Contribute

void TPZSpaceTimeRichardsEq::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	
	const STATE v2 = bc.Val2()(0,0);
	TPZFMatrix<REAL> &phi = data.phi;
	const int phr = phi.Rows();
	int in, jn;
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	
	switch (bc.Type()){
			
			// Dirichlet condition
		case 0 : {
			for(in = 0 ; in < phr; in++) {
				ef(in,0) += weight * ( gBigNumber * phi(in,0) * (v2 - data.sol[0][0]) );
				for (jn = 0 ; jn < phr; jn++) {
					ek(in,jn) +=  gBigNumber * phi(in,0) * phi(jn,0) * weight;
				}
			}
			break;
		}
			
			// Neumann condition
		case 1:{
			// please implement me
		}
			break;
			
			// outflow condition
		case 3 : { 
			
			const STATE sol = data.sol[0][0];
			STATE ConvDir[2] = {0., 1.}; 
			STATE normal[2];
			normal[0] = data.axes(0,1);
			normal[1] = -1.*data.axes(0,0);
			
			STATE ConvNormal = ConvDir[0]*normal[0] + ConvDir[1]*normal[1];
			if(ConvNormal > 0.) {
				for(int il = 0; il < phr; il++) {
					for(int jl = 0; jl < phr; jl++) {
						ek(il,jl) += weight * ConvNormal * phi(il)*phi(jl);
					}
					ef(il,0) += -1. * weight * ConvNormal * phi(il) * sol;
				}
			}
			else{
				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
			}  
		}
			break;
			
		default:{
			std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
		}
	}//switch
	
}//ContributeBC

STATE TPZSpaceTimeRichardsEq::C_Coef(STATE sol){
	
	sol = sol/LCoeff;
	STATE n = this->fN;
	STATE m = 1.-(1./n);
	STATE TR = this->fThetaR;
	STATE TS = this->fThetaS;
	STATE a = this->fAlpha;
	STATE result = (m*n*pow(pow(a,(STATE)2.)*pow(sol,(STATE)2.),n/2.)*pow(1./(1. + pow(pow(a,(STATE)2.)*pow(sol,(STATE)2.),n/2.)),1. + m)*(TR - TS))/sol;
	return result/LCoeff;
}

STATE TPZSpaceTimeRichardsEq::Se(STATE sol){
	STATE n = this->fN;
	STATE m = 1.-(1./n);
	STATE result = pow((STATE)(1./(1.+pow(fabs(this->fAlpha*sol),n))),m);
	return result;
}

STATE TPZSpaceTimeRichardsEq::Theta(STATE Se){
	STATE result = Se*(fThetaS-fThetaR)+fThetaR;
	return result;
}

/** @brief Return sign of the STATE A. If A is closed to zero up to tolerance tol returns zero (in absolute value) */
int Sign(double A, double tol){
	if(fabs(A) < tol) return 0;
	if(A > 0.) return +1;
	return -1;
}

STATE TPZSpaceTimeRichardsEq::DCDsol(STATE sol){
	STATE sol1 = sol*(1.-deltaDerivada);
	STATE antes =   this->C_Coef(sol1);
	STATE sol2 = sol * (1.+deltaDerivada);
	STATE depois = this->C_Coef(sol2);
	STATE result = (depois-antes)/(sol2-sol1);
	return result;
}

void TPZSpaceTimeRichardsEq::AnalysisOfParameters(STATE sol0, STATE solL, char* filename){
	int np = 100;
	TPZFMatrix<STATE> C(np,2), K(np,2), dCdSol(np,2), dKdsol(np,2);
	double delta = (solL - sol0)/(np-1);
	double sol;
	for(int i = 0; i < np; i++){
		sol = sol0+i*delta;
		C(i,0) = sol;
		C(i,1) = this->C_Coef(sol);
		K(i,0) = sol;
		K(i,1) = this->K_Coef(sol);
		dCdSol(i,0) = sol;
		dCdSol(i,1) = this->DCDsol(sol);
		dKdsol(i,0) = sol;
		dKdsol(i,1) = this->DKDsol(sol);
	}
	
	std::ofstream file(filename);
	C.Print("C=",file,EMathematicaInput);
	K.Print("K=",file,EMathematicaInput);
	dCdSol.Print("dCdSol=",file,EMathematicaInput);
	dKdsol.Print("dKdsol=",file,EMathematicaInput);
	
	K.Print("K=",file);	
}

STATE TPZSpaceTimeRichardsEq::K_Coef(STATE sol) {
	
	sol = sol / LCoeff;
	
	STATE n = this->fN;
	STATE m = 1.-(1./n);
	STATE Se = this->Se(sol);
	STATE Ks = this->fKs;
	STATE result = Ks*sqrt(Se)*pow(((STATE)1.)-pow(((STATE)1.)-pow(Se,((STATE)1.)/m),m),(STATE)2.);
	
	return result*LCoeff/TCoeff;
}


STATE TPZSpaceTimeRichardsEq::DKDsol(STATE sol){
	
    STATE sol1 = sol*(1.-deltaDerivada);
    STATE antes = this->K_Coef(sol1);
    STATE sol2 = sol * (1.+deltaDerivada);
    STATE depois = this->K_Coef(sol2);
    STATE resposta = (depois-antes)/(sol2-sol1);
    return resposta;
	
	sol = sol / LCoeff;
	STATE n = this->fN;
	STATE m = 1.-(1./n);  
	STATE Se = this->Se(sol);
	STATE DkDSe = 0.5 * this->fKs * (1.-pow(((STATE)1.)-pow(Se,((STATE)1.)/m),m))*(4.*pow(Se,((STATE)-0.5)+(((STATE)1.)/m))*pow(((STATE)1.)-pow(Se,((STATE)1.)/m),m-((STATE)1.))-(-1.+pow(((STATE)1.)-pow(Se,((STATE)1.)/m),m))/sqrt(Se));
	STATE DSeDsol = -(1./sol)*m*n*pow(fabs(this->fAlpha*sol),n)*pow(1./(1.+pow(fabs(this->fAlpha*sol),n)),m+1.);
	STATE dkdsol = DkDSe * DSeDsol;
	return dkdsol*LCoeff/(TCoeff*LCoeff);
}

int TPZSpaceTimeRichardsEq::ClassId() const{
    return Hash("TPZSpaceTimeRichardsEq") ^ TPZMaterial::ClassId() << 1;
}

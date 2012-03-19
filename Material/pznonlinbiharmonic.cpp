/**
 * \file
 * @brief Contains implementations of the TPZNonLinBiharmonic methods.
 */
//$Id: pznonlinbiharmonic.cpp,v 1.6 2008-10-08 02:09:28 phil Exp $

#include "pznonlinbiharmonic.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>
#include <cmath>

//	NIPG	SIPG	SSIPG1	SSIPG2
REAL TPZNonLinBiharmonic::gLambda1 = 1.0; //	-1	1	-1	1
REAL TPZNonLinBiharmonic::gLambda2 = 1.0; //	-1	1	1	-1
REAL TPZNonLinBiharmonic::gSigmaA  = 10.0;//	 10	10	10	10
REAL TPZNonLinBiharmonic::gSigmaB  = 10.0;//	 10	10	10	10
REAL TPZNonLinBiharmonic::gL_alpha = 6.0; //		  [0, 6
REAL TPZNonLinBiharmonic::gM_alpha = 3.0; //            IGUAL
REAL TPZNonLinBiharmonic::gL_betta = 4.0; //            [-2, 4
REAL TPZNonLinBiharmonic::gM_betta = 1.0; //            IGUAL
REAL TPZNonLinBiharmonic::g_teta = 0.5; // Parametro da parte advectiva.
REAL TPZNonLinBiharmonic::Re = 50.0; // 
int TPZNonLinBiharmonic::NorP = 1; // Constante. Se for 1, entao Metodo de Newton
//            Se for 0, entao Metodo de Picard

using namespace std;

TPZNonLinBiharmonic::TPZNonLinBiharmonic(int nummat, REAL f) : TPZDiscontinuousGalerkin(nummat),
fXf(f){}

TPZNonLinBiharmonic::~TPZNonLinBiharmonic() {
}

void TPZNonLinBiharmonic::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	TPZMaterial::Print(out);
	
}


// Fazer outro metodo igual a esse sem os ... ek ... e com os parametros adequados

void TPZNonLinBiharmonic::Contribute(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix<REAL> &ek,
                                     TPZFMatrix<REAL> &ef) {
	TPZFMatrix<REAL> &dphi = data.dphix;
	// TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	// TPZFMatrix<REAL> &dphiR = dataright.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	// TPZFMatrix<REAL> &phiL = dataleft.phi;
	// TPZFMatrix<REAL> &phiR = dataright.phi;
	// TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=dataleft.sol;
	// TPZVec<REAL> &solR=dataright.sol;
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZFMatrix<REAL> &dsol=data.dsol[0];
	// TPZFMatrix<REAL> &dsolL=dataleft.dsol;
	// TPZFMatrix<REAL> &dsolR=dataright.dsol;
	// REAL &faceSize=data.HSize;
	
	int phr = phi.Rows();
	REAL Re_1 = 1./Re;
	
	if(fForcingFunction) {            // phi(in, 0) = phi_in
		TPZManVector<REAL> res(1);
		fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
		fXf = res[0];
	}
	//Equa�o de non linear biharmonic
	for( int in = 0; in < phr; in++ ) {
		ef(in, 0) +=  weight * fXf*phi(in,0)
		- weight*(Re_1*dsol(2,0)*dphi(2,in)
				  + dsol(2,0)*(dsol(1,0)*dphi(0,in) - dsol(0,0)*dphi(1,in)) ) ;
		
		for( int jn = 0; jn < phr; jn++ ) {
			ek(in,jn) +=  weight * ( Re_1*dphi(2,in) * dphi(2,jn)   
									
									+   dphi(2,jn)*(dsol(1,0)*dphi(0,in) - dsol(0,0)*dphi(1,in)) ); 
			if(NorP==1)
				ek(in,jn) +=  weight * ( dsol(2,0)*(dphi(1,jn)*dphi(0,in) - dphi(0,jn)*dphi(1,in) ) );
		}
	}
}


void TPZNonLinBiharmonic::ContributeBC(TPZMaterialData &data,
                                       REAL weight,
                                       TPZFMatrix<REAL> &ek,
                                       TPZFMatrix<REAL> &ef,
                                       TPZBndCond &bc) {
	
	//NOT TO BE DONE HERE
	PZError << "TPZBiHarminic::ContributeBC - It should never be called.";
}

/** returns the variable index associated with the name*/
int TPZNonLinBiharmonic::VariableIndex(const std::string &name){
	if(!strcmp("Displacement6",name.c_str()))   return  0;
	if(!strcmp("Solution",name.c_str()))        return  1;
	if(!strcmp("Derivate",name.c_str()))        return  2;
	if(!strcmp("POrder",name.c_str()))          return 10;
	cout << "TPZNonLinBiharmonic::VariableIndex Error\n";
	return -1;
}

int TPZNonLinBiharmonic::NSolutionVariables(int var){
	if(var == 0) return 6;
	if(var == 1) return 1;
	if(var == 2) return 2;
	if(var == 10) return 1;
	return TPZMaterial::NSolutionVariables(var);
	//  cout << "TPZNonLinBiharmonic::NSolutionVariables Error\n";
	//  return 0;
}

void TPZNonLinBiharmonic::Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &/*axes*/,
								   int var,TPZVec<REAL> &Solout){
	if(var == 0 || var == 1) Solout[0] = Sol[0];//function
	if(var == 2) {
		Solout.Resize(DSol.Rows());
		int id;
		for(id=0 ; id  < DSol.Rows(); id++) {
			Solout[id] = DSol(id,0);//derivate
		}
	}
}

void TPZNonLinBiharmonic::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/,
							   TPZFMatrix<REAL> &/*DSol*/, TPZFMatrix<REAL> &/*axes*/,
							   TPZVec<REAL> &/*flux*/) {
	//Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux)
}

void TPZNonLinBiharmonic::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u, TPZFMatrix<REAL> &dudx,
								 TPZFMatrix<REAL> &axes, TPZVec<REAL> &/*flux*/,
								 TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,
								 TPZVec<REAL> &values) {
	
	TPZVec<REAL> sol(1), dsol(8,0.);
	Solution(u,dudx,axes,1,sol);
	Solution(u,dudx,axes,2,dsol);
    //values[1] : error em norma L2
	values[1]  = (sol[0] - u_exact[0])*(sol[0] - u_exact[0]);
	
	//values[2] : erro em semi norma H1
	values[2] = 0.;
	for(int id=0; id<2; id++) {
		values[2]  += (dsol[id] - du_exact(id,0))*(dsol[id] - du_exact(id,0));
	}
	//values[3] : erro em semi norma Laplace
	values[3]  = (dsol[2] - du_exact(2,0))*(dsol[2] - du_exact(2,0));
	
	//values[0] : erro em norma H1
	values[0]  = values[1]+values[2];
	// dxx
	values[5] = (dsol[5] - du_exact(5,0))*(dsol[5] - du_exact(5,0));
	// dyy
	values[6] = (dsol[6] - du_exact(6,0))*(dsol[6] - du_exact(6,0));
	// dxy
	values[7] = (dsol[7] - du_exact(7,0))*(dsol[7] - du_exact(7,0));
	//values[4] : erro em norma H2
	values[4]  = values[5]+values[6]+values[7]+values[0];
	
}





void TPZNonLinBiharmonic::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                              REAL weight,
                                              TPZFMatrix<REAL> &ek,
                                              TPZFMatrix<REAL> &ef){
	// TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &dphiR = dataright.dphix;
	// TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	int LeftPOrder=dataleft.p;
	int RightPOrder=dataright.p;
	// TPZVec<REAL> &sol=data.sol;
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZVec<REAL> &solL=dataleft.sol[0];
	TPZVec<REAL> &solR=dataright.sol[0];
	// TPZFMatrix<REAL> &dsol=data.dsol;
	TPZFMatrix<REAL> &dsolL=dataleft.dsol[0];
	TPZFMatrix<REAL> &dsolR=dataright.dsol[0];
	REAL faceSize=data.HSize;
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr,id;
	REAL Re_1 = 1./Re;
	
	REAL alpha=gSigmaA*(pow(((REAL)LeftPOrder),gL_alpha)+pow(((REAL)RightPOrder), gL_alpha) ) /
    (2. * pow(faceSize, gM_alpha) );
	// cout <<  "faceSize em ContributeInterface = " <<faceSize << endl;
	//  cout <<  "LeftPOrder em ContributeInterface = " <<LeftPOrder << endl;
	//  cout <<  "alpha em ContributeInterface = " <<alpha << endl;
	
	
	REAL betta=gSigmaB*(pow(((REAL)LeftPOrder),gL_betta)+pow(((REAL)RightPOrder), gL_betta) ) /
    (2. * pow(faceSize, gM_betta) );
	// cout <<  "betta em ContributeInterface = " <<betta << endl;
	
	// 
	// advectivo
	//
	// 'Unica Integral 
	REAL norFnorF = normal[0]*normal[0];
	REAL norSnorS = normal[1]*normal[1];
    
	REAL uLnorF=  dsolL(1,0)*normal[0];
	REAL vLnorS= -dsolL(0,0)*normal[1];
	
	REAL uRnorF=  dsolR(1,0)*normal[0];  //
	REAL vRnorS= -dsolR(0,0)*normal[1];
    
	REAL absUnorL = fabs(uLnorF + vLnorS);
	REAL absUnorR = fabs(uRnorF + vRnorS);
	
	
	for(il=0; il<nrowl; il++) {
		ef(il,0) +=  weight*0.5*( dsolL(2,0)*phiL(il,0)*(uLnorF + vLnorS))
		+weight*g_teta*absUnorL*(dsolL(2,0)*phiL(il,0)*norFnorF + 
								 dsolL(2,0)*phiL(il,0)*norSnorS);
		
		for(jl=0; jl<nrowl; jl++) {
			REAL deriv_uLnorF =  dphiL(1,jl)*normal[0];
			REAL deriv_vLnorS = -dphiL(0,jl)*normal[1];
			
			REAL deriv_absUnorL = 0.;
			if(uLnorF+vLnorS > 0.) deriv_absUnorL += (deriv_uLnorF+deriv_vLnorS);
			else deriv_absUnorL -= (deriv_uLnorF+deriv_vLnorS);
			
			ek(il,jl) -= weight*0.5*( dphiL(2,jl)*phiL(il,0)*(uLnorF + vLnorS))
			+weight*g_teta*absUnorL*(dphiL(2,jl)*phiL(il,0)*norFnorF + 
									 dphiL(2,jl)*phiL(il,0)*norSnorS) ; 
			
			if(NorP==1) 					  
				ek(il,jl) -= + weight*0.5*(dsolL(2,0)*phiL(il,0)*(deriv_uLnorF+deriv_vLnorS)) +
				weight*g_teta*deriv_absUnorL*(dsolL(2,0)*phiL(il,0)*norFnorF +
											  dsolL(2,0)*phiL(il,0)*norSnorS );
		}
	}
	
	
	for(ir=0; ir<nrowr; ir++) {
		//
		ef(ir+nrowl,0) += weight*0.5*(
									  - dsolR(2,0)*phiR(ir,0)*(uRnorF + vRnorS))
		+weight*g_teta*absUnorR*(dsolR(2,0)*phiR(ir,0)*norFnorF + 
								 dsolR(2,0)*phiR(ir,0)*norSnorS);
		
		for(jr=0; jr<nrowr; jr++) {
			REAL deriv_uRnorF =  dphiR(1,jr)*normal[0];
			REAL deriv_vRnorS = -dphiR(0,jr)*normal[1];
			
			REAL deriv_absUnorR = 0.;
			if(uRnorF + vRnorS > 0.) deriv_absUnorR += (deriv_uRnorF+deriv_vRnorS);
			else deriv_absUnorR -= (deriv_uRnorF+deriv_vRnorS);
			
			ek(ir+nrowl,jr+nrowl) -= weight*0.5*(
												 - dphiR(2,jr)*phiR(ir,0)*(uRnorF + vRnorS))
			+weight*g_teta*absUnorR*(dphiR(2,jr)*phiR(ir,0)*norFnorF + 
									 dphiR(2,jr)*phiR(ir,0)*norSnorS );
			if(NorP==1)					  
				ek(ir+nrowl,jr+nrowl) -=   + weight*0.5*(
														 - dsolR(2,0)*phiR(ir,0)*(deriv_uRnorF + deriv_vRnorS))
				+weight*g_teta*deriv_absUnorR*(dsolR(2,0)*phiR(ir,0)*norFnorF + 
											   dsolR(2,0)*phiR(ir,0)*norSnorS);
		}
	}
	
	
	for(il=0; il<nrowl; il++) {
		//
		ef(il,0) += weight*0.5*(
								dsolR(2,0)*phiL(il,0)*(uRnorF + vRnorS))      
		-weight*g_teta*absUnorR*(dsolR(2,0)*phiL(il,0)*norFnorF + 
								 dsolR(2,0)*phiL(il,0)*norSnorS);
		
		
		for(jr=0; jr<nrowr; jr++) {
			REAL deriv_uRnorF =  dphiR(1,jr)*normal[0];
			REAL deriv_vRnorS = -dphiR(0,jr)*normal[1];
			
			REAL deriv_absUnorR = 0.;
			if(uRnorF + vRnorS > 0.) deriv_absUnorR += (deriv_uRnorF+deriv_vRnorS);
			else deriv_absUnorR -= (deriv_uRnorF+deriv_vRnorS);
			
			ek(il,jr+nrowl) -= weight*0.5*(
										   dphiR(2,jr)*phiL(il,0)*(uRnorF + vRnorS))      
			-weight*g_teta*absUnorR*(dphiR(2,jr)*phiL(il,0)*norFnorF + 
									 dphiR(2,jr)*phiL(il,0)*norSnorS) ;
			
			if(NorP==1)
				ek(il,jr+nrowl) -=   +weight*0.5*(
												  dsolR(2,0)*phiL(il,0)*(deriv_uRnorF + deriv_vRnorS))      
				-weight*g_teta*deriv_absUnorR*(dsolR(2,0)*phiL(il,0)*norFnorF + 
											   dsolR(2,0)*phiL(il,0)*norSnorS) ;
		}
	}
	
	
	for(ir=0; ir<nrowr; ir++) {
		ef(ir+nrowl,0) += weight*0.5*(
									  - dsolL(2,0)*phiR(ir,0)*(uLnorF + vLnorS)) 
		-weight*g_teta*absUnorL*(dsolL(2,0)*phiR(ir,0)*norFnorF + 
								 dsolL(2,0)*phiR(ir,0)*norSnorS);
		
		for(jl=0; jl<nrowl; jl++) {
			REAL deriv_uLnorF =  dphiL(1,jl)*normal[0];
			REAL deriv_vLnorS = -dphiL(0,jl)*normal[1];
			
			REAL deriv_absUnorL = 0.;
			if(uLnorF+vLnorS > 0.) deriv_absUnorL += (deriv_uLnorF+deriv_vLnorS);
			else deriv_absUnorL -= (deriv_uLnorF+deriv_vLnorS);
			
			ek(ir+nrowl,jl) -= weight*0.5*(
										   - dphiL(2,jl)*phiR(ir,0)*(uLnorF + vLnorS))            
			-weight*g_teta*absUnorL*(dphiL(2,jl)*phiR(ir,0)*norFnorF + 
									 dphiL(2,jl)*phiR(ir,0)*norSnorS );
			if(NorP==1)
				ek(ir+nrowl,jl) -=  +weight*0.5*(
												 - dsolL(2,0)*phiR(ir,0)*(deriv_uLnorF + deriv_vLnorS))
				
				-weight*g_teta*deriv_absUnorL*(dsolL(2,0)*phiR(ir,0)*norFnorF + 
											   dsolL(2,0)*phiR(ir,0)*norSnorS);
			
		}
	}
	
	
	/*
	 * biharmonic 
	 */
	
	/* Primeira Integral */
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(3+id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(3+id,jl)*normal[id];
			}
			
			ek(il,jl) += Re_1*weight*0.5*(+gLambda1*dphiLinormal*phiL(jl,0) + dphiLjnormal*phiL(il,0));
		}
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(3+id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*0.5*(+gLambda1*dphiLinormal*solL[0] + dsolLnormal*phiL(il,0));
	}
	
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(3+id,ir)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(3+id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += Re_1*weight*0.5*(
													  -gLambda1*dphiRinormal*phiR(jr,0) - dphiRjnormal*phiR(ir,0));
		}
		REAL dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(3+id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*0.5*(
										   -gLambda1*dphiRinormal*solR[0] - dsolRnormal*phiR(ir,0));
		
	}
	
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(3+id,il)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(3+id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += Re_1*weight*0.5*(
												- gLambda1*dphiLinormal*phiR(jr,0) + dphiRjnormal*phiL(il,0));
		}
		REAL dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(3+id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*0.5*(
									 - gLambda1*dphiLinormal*solR[0] + dsolRnormal*phiL(il,0));
		
	}
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(3+id,ir)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(3+id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += Re_1*weight*0.5*(
												+ gLambda1*dphiRinormal*phiL(jl,0) - dphiLjnormal*phiR(ir,0));
		}
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(3+id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*0.5*(
										   + gLambda1*dphiRinormal*solL[0] - dsolLnormal*phiR(ir,0));
	}
	
	/* Segunda Integral */
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			
			ek(il,jl) += Re_1*weight*0.5*(
										  - dphiLinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiL(2,il));
		}
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*0.5*(
									 - dphiLinormal*dsolL(2,0) - gLambda2*dsolLnormal*dphiL(2,il) );
		
	}
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += Re_1*weight*0.5*(
													  + dphiRinormal*dphiR(2,jr) + gLambda2*dphiRjnormal*dphiR(2,ir));
		}
		REAL dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*0.5*(
										   + dphiRinormal*dsolR(2,0) + gLambda2*dsolRnormal*dphiR(2,ir));
	}
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += Re_1*weight*0.5*(
												- dphiLinormal*dphiR(2,jr) + gLambda2*dphiRjnormal*dphiL(2,il));
		}
		REAL dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*0.5*(
									 - dphiLinormal*dsolR(2,0) + gLambda2*dsolRnormal*dphiL(2,il));
	}
	
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += Re_1*weight*0.5*(
												+ dphiRinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiR(2,ir));
		}
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*0.5*(
										   + dphiRinormal*dsolL(2,0) - gLambda2*dsolLnormal*dphiR(2,ir));    
		
	}
	
	
	/* Terceira Integral */
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			
			ek(il,jl) += Re_1*weight*(
									  alpha * phiL(jl,0)*phiL(il,0) +
									  betta * dphiLinormal*dphiLjnormal);
		}
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*(
								 alpha * solL[0]*phiL(il,0) +
								 betta * dphiLinormal*dsolLnormal);
	}
	
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += Re_1*weight*(
												  alpha * phiR(jr,0)*phiR(ir,0) +
												  betta * dphiRinormal*dphiRjnormal );
		}
		REAL dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*(
									   alpha * solR[0]*phiR(ir,0) +
									   betta * dphiRinormal*dsolRnormal );
	}
	
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += Re_1*weight*(
											-  alpha * phiR(jr,0)*phiL(il,0)
											-  betta * dphiLinormal*dphiRjnormal);
		}
		REAL dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(id,0)*normal[id];
		}    
		ef(il,0) -= Re_1*weight*(
								 -  alpha * solR[0]*phiL(il,0)
								 -  betta * dphiLinormal*dsolRnormal);    
	}
	
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += Re_1*weight*(
											- alpha * phiL(jl,0)*phiR(ir,0)
											- betta * dphiRinormal*dphiLjnormal);
		}
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*(
									   - alpha * solL[0]*phiR(ir,0)
									   - betta * dphiRinormal*dsolLnormal);
	}
}

void TPZNonLinBiharmonic::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                                REAL weight, 
                                                TPZFMatrix<REAL> &ek,
                                                TPZFMatrix<REAL> &ef,
                                                TPZBndCond &bc) {
	// TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	// TPZFMatrix<REAL> &dphiR = dataright.dphix;
	// TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	// TPZFMatrix<REAL> &phiR = dataright.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	int POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZVec<REAL> &solL=dataleft.sol[0];
	// TPZVec<REAL> &solR=dataright.sol;
	// TPZFMatrix<REAL> &dsol=data.dsol;
	TPZFMatrix<REAL> &dsolL=dataleft.dsol[0];
	// TPZFMatrix<REAL> &dsolR=dataright.dsol;
	REAL faceSize=data.HSize;
	
	REAL alpha = gSigmaA*pow(((REAL)POrder), gL_alpha) /  pow(faceSize, gM_alpha);
	REAL betta = gSigmaB*pow(((REAL)POrder), gL_betta) /  pow(faceSize, gM_betta);
	//   cout <<  "faceSize em ContributeBCInterface = " <<faceSize << endl;
	//   cout <<  "POrder em ContributeBCInterface = " <<POrder << endl;
	//   cout <<  "alpha em ContributeBCInterface = " <<alpha << endl;
	//   cout <<  "betta em ContributeBCInterface = " <<betta << endl;
	
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	REAL Re_1 = 1./Re;
	
	/* Primeira Integral */
	for(il=0; il<nrowl; il++) {
		
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(3+id,il)*normal[id];
		}
		
		
		// Termos de Dirichlet     -     em Val2()(0,0)
		ef(il,0) += + Re_1*gLambda1*weight*dphiLinormal*bc.Val2()(0,0)
		+ Re_1*alpha * weight*bc.Val2()(0,0)*phiL(il) ;
		
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(3+id,jl)*normal[id];
			}
			ek(il,jl) += Re_1*weight*(+ gLambda1*dphiLinormal*phiL(jl,0)+ dphiLjnormal*phiL(il,0)); //2 1
		}
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(3+id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*(+ gLambda1*dphiLinormal*solL[0]+ dsolLnormal*phiL(il,0));
		
	}
	
	/* Segunda Integral */
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		
		// Termos de Neuwmann     -      em Val2()(1,0)
		ef(il,0) += - Re_1*gLambda2*weight*bc.Val2()(1,0)*dphiL(2,il)
		+ Re_1*betta * weight*bc.Val2()(1,0)*dphiLinormal ;
		
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(il,jl) += Re_1*weight*(- dphiLinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiL(2,il) );
		}
		
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*(- dphiLinormal*dsolL(2,0) - gLambda2*dsolLnormal*dphiL(2,il) );
		
	}
	
	
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(il,jl) += Re_1*weight*(
									  alpha * phiL(jl,0)*phiL(il,0) +
									  betta * dphiLinormal*dphiLjnormal);
		}
		REAL dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*(
								 alpha * solL[0]*phiL(il,0) +
								 betta * dphiLinormal*dsolLnormal);
		
	}
	
	/*
	 *  Termo advectivo
	 */
	
	/* REAL ULnormal = dsolL(1,0)*normal[0] - dsolL(0,0)*normal[1];
	 for(il=0; il<nrowl; il++) 
	 for(jl=0; jl<nrowl; jl++) 
	 ek(il,jl) -= weight*dphiL(2,jl)*phiL(il,0)*ULnormal;
	 
	 */
	
	REAL ULnormal = dsolL(1,0)*normal[0] - dsolL(0,0)*normal[1];
	for(il=0; il<nrowl; il++){ 
		for(jl=0; jl<nrowl; jl++){ 
			ek(il,jl) -= weight*dphiL(2,jl)*phiL(il,0)*ULnormal;
			
			if(NorP==1)
				ek(il,jl) -=   +weight*dsolL(2,0)*phiL(il,0)*( dphiL(1,jl)*normal[0] - dphiL(0,jl)*normal[1]);
		}
		
		ef(il,0) += weight*dsolL(2,0)*phiL(il,0)*ULnormal;
		
	}

}

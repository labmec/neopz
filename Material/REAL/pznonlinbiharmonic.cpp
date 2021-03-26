/**
 * @file
 * @brief Contains implementations of the TPZNonLinBiharmonic methods.
 */

#include "pznonlinbiharmonic.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>
#include <cmath>

//	NIPG	SIPG	SSIPG1	SSIPG2
STATE TPZNonLinBiharmonic::gLambda1 = 1.0; //	-1	1	-1	1
STATE TPZNonLinBiharmonic::gLambda2 = 1.0; //	-1	1	1	-1
STATE TPZNonLinBiharmonic::gSigmaA  = 10.0;//	 10	10	10	10
STATE TPZNonLinBiharmonic::gSigmaB  = 10.0;//	 10	10	10	10
STATE TPZNonLinBiharmonic::gL_alpha = 6.0; //		  [0, 6
STATE TPZNonLinBiharmonic::gM_alpha = 3.0; //            IGUAL
STATE TPZNonLinBiharmonic::gL_betta = 4.0; //            [-2, 4
STATE TPZNonLinBiharmonic::gM_betta = 1.0; //            IGUAL
STATE TPZNonLinBiharmonic::g_teta = 0.5; // Parametro da parte advectiva.
STATE TPZNonLinBiharmonic::Re = 50.0; // 
int TPZNonLinBiharmonic::NorP = 1; // Constante. Se for 1, entao Metodo de Newton
//            Se for 0, entao Metodo de Picard

using namespace std;

TPZNonLinBiharmonic::TPZNonLinBiharmonic(int nummat, STATE f) : TPZRegisterClassId(&TPZNonLinBiharmonic::ClassId),
TPZMaterial(nummat),
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
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef) {
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZManVector<REAL,3> &x = data.x;
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	TPZFMatrix<STATE> &dsol=data.dsol[0];
	
	int phr = phi.Rows();
	STATE Re_1 = 1./Re;
	
	if(fForcingFunction) {            // phi(in, 0) = phi_in
		TPZManVector<STATE> res(1);
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
                                       TPZFMatrix<STATE> &ek,
                                       TPZFMatrix<STATE> &ef,
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
}

void TPZNonLinBiharmonic::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &/*axes*/,
								   int var,TPZVec<STATE> &Solout){
	if(var == 0 || var == 1) Solout[0] = Sol[0];//function
	if(var == 2) {
		Solout.Resize(DSol.Rows());
		int id;
		for(id=0 ; id  < DSol.Rows(); id++) {
			Solout[id] = DSol(id,0);//derivate
		}
	}
}


void TPZNonLinBiharmonic::Errors(TPZVec<REAL> &/*x*/,TPZVec<STATE> &u, TPZFMatrix<STATE> &dudx,
								 TPZFMatrix<REAL> &axes,
								 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,
								 TPZVec<REAL> &values) {
	
	TPZVec<STATE> sol(1), dsol(8,0.);
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
                                              TPZFMatrix<STATE> &ek,
                                              TPZFMatrix<STATE> &ef){
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &dphiR = dataright.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	
	int LeftPOrder=dataleft.p;
	int RightPOrder=dataright.p;
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	TPZVec<STATE> &solL=dataleft.sol[0];
	TPZVec<STATE> &solR=dataright.sol[0];
	TPZFMatrix<STATE> &dsolL=dataleft.dsol[0];
	TPZFMatrix<STATE> &dsolR=dataright.dsol[0];
	STATE faceSize=data.HSize;
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr,id;
	STATE Re_1 = 1./Re;
	
	STATE alpha=gSigmaA*(pow(((STATE)LeftPOrder),gL_alpha)+pow(((STATE)RightPOrder), gL_alpha) ) /
    (2. * pow(faceSize, gM_alpha) );
	
	
	STATE betta=gSigmaB*(pow(((STATE)LeftPOrder),gL_betta)+pow(((STATE)RightPOrder), gL_betta) ) /
    (2. * pow(faceSize, gM_betta) );
	
	// 
	// advectivo
	// 'Unica Integral 
	STATE norFnorF = normal[0]*normal[0];
	STATE norSnorS = normal[1]*normal[1];
    
	STATE uLnorF=  dsolL(1,0)*normal[0];
	STATE vLnorS= -dsolL(0,0)*normal[1];
	
	STATE uRnorF=  dsolR(1,0)*normal[0];  //
	STATE vRnorS= -dsolR(0,0)*normal[1];
    
	STATE absUnorL = fabs(uLnorF + vLnorS);
	STATE absUnorR = fabs(uRnorF + vRnorS);
	
	
	for(il=0; il<nrowl; il++) {
		ef(il,0) +=  weight*0.5*( dsolL(2,0)*phiL(il,0)*(uLnorF + vLnorS))
		+weight*g_teta*absUnorL*(dsolL(2,0)*phiL(il,0)*norFnorF + 
								 dsolL(2,0)*phiL(il,0)*norSnorS);
		
		for(jl=0; jl<nrowl; jl++) {
			STATE deriv_uLnorF =  dphiL(1,jl)*normal[0];
			STATE deriv_vLnorS = -dphiL(0,jl)*normal[1];
			
			STATE deriv_absUnorL = 0.;
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
			STATE deriv_uRnorF =  dphiR(1,jr)*normal[0];
			STATE deriv_vRnorS = -dphiR(0,jr)*normal[1];
			
			STATE deriv_absUnorR = 0.;
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
			STATE deriv_uRnorF =  dphiR(1,jr)*normal[0];
			STATE deriv_vRnorS = -dphiR(0,jr)*normal[1];
			
			STATE deriv_absUnorR = 0.;
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
			STATE deriv_uLnorF =  dphiL(1,jl)*normal[0];
			STATE deriv_vLnorS = -dphiL(0,jl)*normal[1];
			
			STATE deriv_absUnorL = 0.;
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
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(3+id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(3+id,jl)*normal[id];
			}
			
			ek(il,jl) += Re_1*weight*0.5*(+gLambda1*dphiLinormal*phiL(jl,0) + dphiLjnormal*phiL(il,0));
		}
		STATE dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(3+id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*0.5*(+gLambda1*dphiLinormal*solL[0] + dsolLnormal*phiL(il,0));
	}
	
	for(ir=0; ir<nrowr; ir++) {
		STATE dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(3+id,ir)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			STATE dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(3+id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += Re_1*weight*0.5*(
													  -gLambda1*dphiRinormal*phiR(jr,0) - dphiRjnormal*phiR(ir,0));
		}
		STATE dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(3+id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*0.5*(
										   -gLambda1*dphiRinormal*solR[0] - dsolRnormal*phiR(ir,0));
		
	}
	
	for(il=0; il<nrowl; il++) {
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(3+id,il)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			STATE dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(3+id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += Re_1*weight*0.5*(
												- gLambda1*dphiLinormal*phiR(jr,0) + dphiRjnormal*phiL(il,0));
		}
		STATE dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(3+id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*0.5*(
									 - gLambda1*dphiLinormal*solR[0] + dsolRnormal*phiL(il,0));
		
	}
	for(ir=0; ir<nrowr; ir++) {
		STATE dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(3+id,ir)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(3+id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += Re_1*weight*0.5*(
												+ gLambda1*dphiRinormal*phiL(jl,0) - dphiLjnormal*phiR(ir,0));
		}
		STATE dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(3+id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*0.5*(
										   + gLambda1*dphiRinormal*solL[0] - dsolLnormal*phiR(ir,0));
	}
	
	/* Segunda Integral */
	for(il=0; il<nrowl; il++) {
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			
			ek(il,jl) += Re_1*weight*0.5*(
										  - dphiLinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiL(2,il));
		}
		STATE dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*0.5*(
									 - dphiLinormal*dsolL(2,0) - gLambda2*dsolLnormal*dphiL(2,il) );
		
	}
	for(ir=0; ir<nrowr; ir++) {
		STATE dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			STATE dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += Re_1*weight*0.5*(
													  + dphiRinormal*dphiR(2,jr) + gLambda2*dphiRjnormal*dphiR(2,ir));
		}
		STATE dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*0.5*(
										   + dphiRinormal*dsolR(2,0) + gLambda2*dsolRnormal*dphiR(2,ir));
	}
	for(il=0; il<nrowl; il++) {
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			STATE dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += Re_1*weight*0.5*(
												- dphiLinormal*dphiR(2,jr) + gLambda2*dphiRjnormal*dphiL(2,il));
		}
		STATE dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*0.5*(
									 - dphiLinormal*dsolR(2,0) + gLambda2*dsolRnormal*dphiL(2,il));
	}
	
	for(ir=0; ir<nrowr; ir++) {
		STATE dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += Re_1*weight*0.5*(
												+ dphiRinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiR(2,ir));
		}
		STATE dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*0.5*(
										   + dphiRinormal*dsolL(2,0) - gLambda2*dsolLnormal*dphiR(2,ir));    
		
	}
	
	/* Terceira Integral */
	for(il=0; il<nrowl; il++) {
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			
			ek(il,jl) += Re_1*weight*(
									  alpha * phiL(jl,0)*phiL(il,0) +
									  betta * dphiLinormal*dphiLjnormal);
		}
		STATE dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*(
								 alpha * solL[0]*phiL(il,0) +
								 betta * dphiLinormal*dsolLnormal);
	}
	
	for(ir=0; ir<nrowr; ir++) {
		STATE dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			STATE dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += Re_1*weight*(
												  alpha * phiR(jr,0)*phiR(ir,0) +
												  betta * dphiRinormal*dphiRjnormal );
		}
		STATE dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(id,0)*normal[id];
		}
		ef(ir+nrowl,0) -= Re_1*weight*(
									   alpha * solR[0]*phiR(ir,0) +
									   betta * dphiRinormal*dsolRnormal );
	}
	
	for(il=0; il<nrowl; il++) {
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			STATE dphiRjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += Re_1*weight*(
											-  alpha * phiR(jr,0)*phiL(il,0)
											-  betta * dphiLinormal*dphiRjnormal);
		}
		STATE dsolRnormal = 0.;
		for(id=0; id<2; id++) {
			dsolRnormal += dsolR(id,0)*normal[id];
		}    
		ef(il,0) -= Re_1*weight*(
								 -  alpha * solR[0]*phiL(il,0)
								 -  betta * dphiLinormal*dsolRnormal);    
	}
	
	for(ir=0; ir<nrowr; ir++) {
		STATE dphiRinormal = 0.;
		for(id=0; id<2; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += Re_1*weight*(
											- alpha * phiL(jl,0)*phiR(ir,0)
											- betta * dphiRinormal*dphiLjnormal);
		}
		STATE dsolLnormal = 0.;
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
                                                TPZFMatrix<STATE> &ek,
                                                TPZFMatrix<STATE> &ef,
                                                TPZBndCond &bc) {
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	int POrder=data.p;
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	TPZVec<STATE> &solL=dataleft.sol[0];
	TPZFMatrix<STATE> &dsolL=dataleft.dsol[0];
	STATE faceSize=data.HSize;
	
	STATE alpha = gSigmaA*pow(((STATE)POrder), gL_alpha) /  pow(faceSize, gM_alpha);
	STATE betta = gSigmaB*pow(((STATE)POrder), gL_betta) /  pow(faceSize, gM_betta);
	
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	STATE Re_1 = 1./Re;
	
	/* Primeira Integral */
	for(il=0; il<nrowl; il++) {
		
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(3+id,il)*normal[id];
		}
		
		// Termos de Dirichlet     -     em Val2()(0,0)
		ef(il,0) += + Re_1*gLambda1*weight*dphiLinormal*bc.Val2()(0,0)
		+ Re_1*alpha * weight*bc.Val2()(0,0)*phiL(il) ;
		
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(3+id,jl)*normal[id];
			}
			ek(il,jl) += Re_1*weight*(+ gLambda1*dphiLinormal*phiL(jl,0)+ dphiLjnormal*phiL(il,0)); //2 1
		}
		STATE dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(3+id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*(+ gLambda1*dphiLinormal*solL[0]+ dsolLnormal*phiL(il,0));
		
	}
	
	/* Segunda Integral */
	for(il=0; il<nrowl; il++) {
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		
		// Termos de Neuwmann     -      em Val2()(1,0)
		ef(il,0) += - Re_1*gLambda2*weight*bc.Val2()(1,0)*dphiL(2,il)
		+ Re_1*betta * weight*bc.Val2()(1,0)*dphiLinormal ;
		
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(il,jl) += Re_1*weight*(- dphiLinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiL(2,il) );
		}
		
		STATE dsolLnormal = 0.;
		for(id=0; id<2; id++) {
			dsolLnormal += dsolL(id,0)*normal[id];
		}
		ef(il,0) -= Re_1*weight*(- dphiLinormal*dsolL(2,0) - gLambda2*dsolLnormal*dphiL(2,il) );
		
	}
	
	for(il=0; il<nrowl; il++) {
		STATE dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(il,jl) += Re_1*weight*(
									  alpha * phiL(jl,0)*phiL(il,0) +
									  betta * dphiLinormal*dphiLjnormal);
		}
		STATE dsolLnormal = 0.;
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
	
	STATE ULnormal = dsolL(1,0)*normal[0] - dsolL(0,0)*normal[1];
	for(il=0; il<nrowl; il++){ 
		for(jl=0; jl<nrowl; jl++){ 
			ek(il,jl) -= weight*dphiL(2,jl)*phiL(il,0)*ULnormal;
			
			if(NorP==1)
				ek(il,jl) -=   +weight*dsolL(2,0)*phiL(il,0)*( dphiL(1,jl)*normal[0] - dphiL(0,jl)*normal[1]);
		}
		
		ef(il,0) += weight*dsolL(2,0)*phiL(il,0)*ULnormal;
		
	}
	
}

int TPZNonLinBiharmonic::ClassId() const{
    return Hash("TPZNonLinBiharmonic") ^ TPZMaterial::ClassId() << 1;
}

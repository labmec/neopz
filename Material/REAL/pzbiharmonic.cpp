/**
 * \file
 * @brief Contains implementations of the TPZBiharmonic methods.
 */

#include "pzbiharmonic.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

#include <cmath>

//	NIPG	SIPG	SSIPG1	SSIPG2
REAL TPZBiharmonic::gLambda1 = 1.0; //	-1	1	-1	1
REAL TPZBiharmonic::gLambda2 = 1.0; //	-1	1	1	-1
REAL TPZBiharmonic::gSigmaA  = 10.0;//	 10	10	10	10
REAL TPZBiharmonic::gSigmaB  = 10.0;//	 10	10	10	10
REAL TPZBiharmonic::gL_alpha = 6.0; //		  [0, 6]
REAL TPZBiharmonic::gM_alpha = 3.0; //            IGUAL
REAL TPZBiharmonic::gL_betta = 4.0; //            [-2, 4]
REAL TPZBiharmonic::gM_betta = 1.0; //            IGUAL
using namespace std;


TPZBiharmonic::TPZBiharmonic(int nummat, REAL f) 
: TPZRegisterClassId(&TPZBiharmonic::ClassId), TPZMaterial(nummat),
fXf(f){}

TPZBiharmonic::~TPZBiharmonic() {
}

void TPZBiharmonic::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	TPZMaterial::Print(out);
	
}

void TPZBiharmonic::Contribute(TPZMaterialData &data,
                               REAL weight,
                               TPZFMatrix<STATE> &ek,
							   TPZFMatrix<STATE> &ef) {
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZManVector<REAL,3> &x = data.x;
	
	int phr = phi.Rows();
	
	if(fForcingFunction) {            // phi(in, 0) = phi_in
		TPZManVector<STATE> res(1);
        TPZFMatrix<STATE> grad;
		fForcingFunction->Execute(x,res,grad);       // dphi(i,j) = dphi_j/dxi
		fXf = res[0];
	}
	//Equa�o de Poisson
	for( int in = 0; in < phr; in++ ) {
		ef(in, 0) +=  weight * fXf*phi(in,0);
		for( int jn = 0; jn < phr; jn++ ) {
			ek(in,jn) +=  weight * ( dphi(2,in) * dphi(2,jn) );
		}
	}
}


void TPZBiharmonic::ContributeBC(TPZMaterialData &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,
                                 TPZFMatrix<STATE> &ef,
                                 TPZBndCond &bc) {
	
	//NOT TO BE DONE HERE
	PZError << "TPZBiHarminic::ContributeBC - It should never be called.";
}

/** returns the variable index associated with the name*/
int TPZBiharmonic::VariableIndex(const std::string &name){
	if(!strcmp("Displacement6",name.c_str()))   return  0;
	if(!strcmp("Solution",name.c_str()))        return  1;
	if(!strcmp("Derivate",name.c_str()))        return  2;
	if(!strcmp("POrder",name.c_str()))          return 10;
	cout << "TPZBiharmonic::VariableIndex Error\n";
	return -1;
}

int TPZBiharmonic::NSolutionVariables(int var){
	if(var == 0) return 6;
	if(var == 1) return 1;
	if(var == 2) return 2;
	if(var == 10) return 1;
	return TPZMaterial::NSolutionVariables(var);
	//  cout << "TPZBiharmonic::NSolutionVariables Error\n";
	//  return 0;
}

void TPZBiharmonic::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &/*axes*/,
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


void TPZBiharmonic::Errors(TPZVec<REAL> &/*x*/,TPZVec<STATE> &u, TPZFMatrix<STATE> &dudx,
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

void TPZBiharmonic::ContributeInterface(TPZMaterialData &data , TPZMaterialData &dataleft, TPZMaterialData &dataright,
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
	REAL faceSize=data.HSize;	
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr,id;
	
	REAL alpha=gSigmaA*(pow(((REAL)LeftPOrder),gL_alpha)+pow(((REAL)RightPOrder), gL_alpha) ) /
    (2. * pow(faceSize, gM_alpha) );
	
	REAL betta=gSigmaB*(pow(((REAL)LeftPOrder),gL_betta)+pow(((REAL)RightPOrder), gL_betta) ) /
    (2. * pow(faceSize, gM_betta) );
	
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
			
			ek(il,jl) += weight*0.5*(+gLambda1*dphiLinormal*phiL(jl,0) + dphiLjnormal*phiL(il,0));
		}
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
			ek(ir+nrowl,jr+nrowl) += weight*0.5*(
												 -gLambda1*dphiRinormal*phiR(jr,0) - dphiRjnormal*phiR(ir,0)
												 );
		}
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
			ek(il,jr+nrowl) += weight*0.5*(
										   - gLambda1*dphiLinormal*phiR(jr,0) + dphiRjnormal*phiL(il,0)
										   );
		}
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
			ek(ir+nrowl,jl) += weight*0.5*(
										   + gLambda1*dphiRinormal*phiL(jl,0) - dphiLjnormal*phiR(ir,0)
										   );
		}
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
			
			ek(il,jl) += weight*0.5*(
									 - dphiLinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiL(2,il)
									 );
		}
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
			ek(ir+nrowl,jr+nrowl) += weight*0.5*(
												 + dphiRinormal*dphiR(2,jr) + gLambda2*dphiRjnormal*dphiR(2,ir)
												 );
		}
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
			ek(il,jr+nrowl) += weight*0.5*(
										   - dphiLinormal*dphiR(2,jr) + gLambda2*dphiRjnormal*dphiL(2,il)
										   );
		}
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
			ek(ir+nrowl,jl) += weight*0.5*(
										   + dphiRinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiR(2,ir)
										   );
		}
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
			
			ek(il,jl) += weight*(
								 alpha * phiL(jl,0)*phiL(il,0) +
								 betta * dphiLinormal*dphiLjnormal
								 );
		}
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
			ek(ir+nrowl,jr+nrowl) += weight*(
											 alpha * phiR(jr,0)*phiR(ir,0) +
											 betta * dphiRinormal*dphiRjnormal
											 );
		}
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
			ek(il,jr+nrowl) += weight*(
									   -  alpha * phiR(jr,0)*phiL(il,0)
									   -  betta * dphiLinormal*dphiRjnormal
									   );
		}
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
			ek(ir+nrowl,jl) += weight*(
									   - alpha * phiL(jl,0)*phiR(ir,0)
									   - betta * dphiRinormal*dphiLjnormal
									   );
		}
	}
}

void TPZBiharmonic::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                          REAL weight,
                                          TPZFMatrix<STATE> &ek,
                                          TPZFMatrix<STATE> &ef,
                                          TPZBndCond &bc) {
	
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	int POrder=data.p;                            // I need some explains, why you use reference & - Jorge
	REAL faceSize=data.HSize;
	
	REAL alpha = gSigmaA*pow(((REAL)POrder), gL_alpha) /  pow(faceSize, gM_alpha);
	REAL betta = gSigmaB*pow(((REAL)POrder), gL_betta) /  pow(faceSize, gM_betta);
	
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	
	/* Primeira Integral */
	for(il=0; il<nrowl; il++) {
		
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(3+id,il)*normal[id];
		}
		
		// Termos de Dirichlet     -     em Val2()(0,0)
		ef(il,0) += + gLambda1*weight*dphiLinormal*bc.Val2()(0,0)
		+ alpha * weight*bc.Val2()(0,0)*phiL(il) ;
		
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(3+id,jl)*normal[id];
			}
			ek(il,jl) += weight*(+ gLambda1*dphiLinormal*phiL(jl,0)+ dphiLjnormal*phiL(il,0)); //2 1
		}
	}
	
	/* Segunda Integral */
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<2; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		
		// Termos de Neuwmann     -      em Val2()(1,0)
		ef(il,0) += - gLambda2*weight*bc.Val2()(1,0)*dphiL(2,il)
		+ betta * weight*bc.Val2()(1,0)*dphiLinormal ;
		
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<2; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(il,jl) += weight*(- dphiLinormal*dphiL(2,jl) - gLambda2*dphiLjnormal*dphiL(2,il) );
		}
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
			ek(il,jl) += weight*(
								 alpha * phiL(jl,0)*phiL(il,0) +
								 betta * dphiLinormal*dphiLjnormal
								 );
		}
	}
}

int TPZBiharmonic::ClassId() const{
    return Hash("TPZBiharmonic") ^ TPZMaterial::ClassId() << 1;
}

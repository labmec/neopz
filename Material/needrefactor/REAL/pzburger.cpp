/**
 * \file
 * @brief Contains implementations of the TPZBurger methods.
 */

#include "pzburger.h"
#include "pzbndcond.h"

using namespace std;

int TPZBurger::gStabilizationScheme = ESUPG;

TPZBurger::TPZBurger(int nummat, int dim):TPZRegisterClassId(&TPZBurger::ClassId),
TPZMatPoisson3dReferred(nummat, dim){
	this->fIsReferred = true;
	this->fSolRef = 1.;
	TPZBurger::gStabilizationScheme = ESUPG;
}

TPZBurger::TPZBurger(const TPZBurger &cp):TPZRegisterClassId(&TPZBurger::ClassId),
TPZMatPoisson3dReferred(cp){
	this->fIsReferred = cp.fIsReferred;
	this->fSolRef = cp.fSolRef;
}

TPZBurger::~TPZBurger(){
	
}

void TPZBurger::ContributeGradStab(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,REAL weight,
								   TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
								   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
	STATE fK = AVGK();

    if (this->IsReferred()){
		this->SetConvectionTerm(dsol, axes);
	}
	
	int phr = phi.Rows();
	
	
	
	if(fForcingFunction) {            // phi(in, 0) = phi_in
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
		fXf = res[0];
	}
	REAL delx = 0.;
	REAL ConvDirAx[3] = {0.};
	if(fC != 0.0) {
		int di,dj;
		delx = 0.;
		for(di=0; di<fDim; di++) {
			for(dj=0; dj<fDim; dj++) {
				delx = (delx<fabs(jacinv(di,dj))) ? fabs(jacinv(di,dj)) : delx;
			}
		}
		delx = 2./delx;
		
		
		switch(fDim) {
			case 1:
				ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
				break;
			case 2:
				ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
				ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
				break;
			case 3:
				ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
				ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
				ConvDirAx[2] = axes(2,0)*fConvDir[0]+axes(2,1)*fConvDir[1]+axes(2,2)*fConvDir[2];
				break;
			default:
				PZError << "TPZBurger::Contribute dimension error " << fDim << endl;
		}
	}
    
	REAL signsol = (sol[0] < 0) ? -1. : +1.;
	for( int in = 0; in < phr; in++ ) {
		int kd;
		REAL dphiic = 0;
		
		ef(in, 0) += - weight * ( fXf*phi(in,0) + 0.5*fSD*delx*fC*dphiic*fXf * fabs(sol[0])/fSolRef );
		for(kd = 0; kd < fDim; kd++){
			ef(in, 0) += -1. * weight * ( +fK * ( dphi(kd,in) * dsol(kd,0) )
										 -fC * ( ConvDirAx[kd]* dphi(kd,in) * sol[0]*sol[0]/fSolRef )
										 +0.5 * fSD * delx * fC * dphiic * dsol(kd,0) * ConvDirAx[kd] * fabs(sol[0])/fSolRef );
		}//kd
		
		for( int jn = 0; jn < phr; jn++ ) {
			ek(in, jn) += weight * (0.5*fSD*delx*fC*dphiic*fXf*phi(jn,0)*signsol/fSolRef + 0.5*fSD*delx*fC*dphiic*fXf * fabs(sol[0])/fSolRef);
			for(kd=0; kd<fDim; kd++) {
				ek(in,jn) += weight * (
									   +fK * ( dphi(kd,in) * dphi(kd,jn) ) 
									   -fC * ( ConvDirAx[kd]* dphi(kd,in) * phi(jn) * 2.*sol[0]/fSolRef )
									   +0.5 * fSD * delx * fC * dphiic *  ConvDirAx[kd] * signsol * (dphi(kd,jn)*sol[0]+dsol(kd,0)*phi(jn,0))/fSolRef
									   );
			}
		}
	}//in
    
	if (this->fC == 0.){    
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}                          
}

void TPZBurger::ContributeSUPG(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,REAL weight,
							   TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
	if (this->IsReferred()){
		this->SetConvectionTerm(dsol, axes);
	}

    STATE fK = AVGK();
	int phr = phi.Rows();
	
	if(fForcingFunction) {            // phi(in, 0) = phi_in
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
		fXf = res[0];
	}
	REAL delx = 0.;
	REAL ConvDirAx[3] = {0.};
	if(fC != 0.0) {
		int di,dj;
		delx = 0.;
		for(di=0; di<fDim; di++) {
			for(dj=0; dj<fDim; dj++) {
				delx = (delx<fabs(jacinv(di,dj))) ? fabs(jacinv(di,dj)) : delx;
			}
		}
		delx = 2./delx;
		
		
		switch(fDim) {
			case 1:
				ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
				break;
			case 2:
				ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
				ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
				break;
			case 3:
				ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
				ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
				ConvDirAx[2] = axes(2,0)*fConvDir[0]+axes(2,1)*fConvDir[1]+axes(2,2)*fConvDir[2];
				break;
			default:
				PZError << "TPZBurger::Contribute dimension error " << fDim << endl;
		}
	}
    
	REAL signsol = (sol[0] < 0) ? -1. : +1.;
	for( int in = 0; in < phr; in++ ) {
		int kd;
		REAL dphiic = 0;
		for(kd = 0; kd<fDim; kd++) dphiic += ConvDirAx[kd]*dphi(kd,in);
		REAL norm = 0.;
		for(kd = 0; kd<fDim; kd++) norm += (ConvDirAx[kd]*fC) * (ConvDirAx[kd]*fC);
		norm = sqrt(norm);
		
		ef(in, 0) += - weight * ( fXf*phi(in,0) + 0.5*fSD*delx*fC*dphiic*fXf * fabs(sol[0])/fSolRef );
		for(kd = 0; kd < fDim; kd++){
			ef(in, 0) += -1. * weight * ( +fK * ( dphi(kd,in) * dsol(kd,0) )
										 -fC * ( ConvDirAx[kd]* dphi(kd,in) * sol[0]*sol[0]/fSolRef )
										 +0.5 * fSD * delx * fC * dphiic * dsol(kd,0) * ConvDirAx[kd] * fabs(sol[0])/fSolRef );
		}//kd
		
		for( int jn = 0; jn < phr; jn++ ) {
			ek(in, jn) += weight * (0.5*fSD*delx*fC*dphiic*fXf*phi(jn,0)*signsol/fSolRef);
			for(kd=0; kd<fDim; kd++) {
				ek(in,jn) += weight * (
									   +fK * ( dphi(kd,in) * dphi(kd,jn) ) 
									   -fC * ( ConvDirAx[kd]* dphi(kd,in) * phi(jn) * 2.*sol[0]/fSolRef )
									   +0.5 * fSD * delx * fC * dphiic *  ConvDirAx[kd] * signsol * (dphi(kd,jn)*sol[0]+dsol(kd,0)*phi(jn,0))/fSolRef
									   );
			}
		}
	}//in
    
	if (fC == 0.){    
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}                          
}

void TPZBurger::ContributeBC(TPZMaterialData &data,
                             REAL weight,
                             TPZFMatrix<STATE> &ek,
                             TPZFMatrix<STATE> &ef,
                             TPZBndCond &bc){
    int numbersol = data.sol.size();
    if(numbersol != 1)
    {
        DebugStop();
    }
	TPZFMatrix<REAL> &phi = data.phi;
	TPZVec<STATE> &sol=data.sol[0];
	TPZFMatrix<REAL> &axes=data.axes;
	
	int phr = phi.Rows();
	short in,jn;
	REAL v2[1];
	v2[0] = bc.Val2()(0,0);
	
	switch (bc.Type()) {
		case 0 : {      // Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				ef(in,0) += weight * ( gBigNumber * v2[0] * phi(in,0) - gBigNumber * phi(in,0) * sol[0] );
				for (jn = 0 ; jn < phr; jn++) {
					ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
				}
			}
		}
			break;
			
		case 1 : {      // Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(in,0) += v2[0] * phi(in,0) * weight;
			}
		}
			break;
			
		case 2 :{    // condicao mista
			cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
		}
			break;
			
		case 3 : { // outflow condition
			
			int id, il, jl;
			REAL normal[3];
			if (fDim == 1) PZError << __PRETTY_FUNCTION__ << " - ERROR! The normal vector is not available for 1D TPZInterpolatedElement\n";
			if (fDim == 2){
				normal[0] = axes(0,1);
				normal[1] = axes(1,1);
			}
			if (fDim == 3){
				normal[0] = axes(0,2);
				normal[1] = axes(1,2);
				normal[2] = axes(2,2);
			}
			REAL ConvNormal = 0.;    
			for(id=0; id<fDim; id++) ConvNormal += fC*fConvDir[id]*normal[id];  
			if(ConvNormal > 0.) {
				for(il=0; il<phr; il++) {
					for(jl=0; jl<phr; jl++) {
						ek(il,jl) += weight * ConvNormal * phi(il)*phi(jl)*2.*sol[0]/fSolRef;
					}
					ef(il,0) += -1. * weight * ConvNormal * phi(il) * sol[0]*sol[0]/fSolRef;
				}
			}
			else{
				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
			}
		}
			break;
			
		default :{
			PZError << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - Error! Wrong boundary condition type\n";
		}
			break;
	}
	
	if (this->IsSymetric()) {//only 1.e-3 because of bignumbers.
		if ( !ek.VerifySymmetry( 1.e-3 ) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
}

void TPZBurger::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                    REAL weight,
                                    TPZFMatrix<STATE> &ek,
                                    TPZFMatrix<STATE> &ef){
    STATE fK = AVGK();

    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &dphiR = dataright.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	TPZVec<STATE> &solL=dataleft.sol[0];
	TPZVec<STATE> &solR=dataright.sol[0];
	TPZFMatrix<STATE> &dsolL=dataleft.dsol[0];
	TPZFMatrix<STATE> &dsolR=dataright.dsol[0];
	
	if (this->IsReferred()){
		this->SetConvectionTermInterface(dsolL, dsolR);
	}
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr,id;
	
	//Convection term
	REAL ConvNormal = 0.;
	for(id=0; id<fDim; id++) ConvNormal += fC * fConvDir[id]*normal[id];
	if(ConvNormal > 0.) {
		for(il=0; il<nrowl; il++) {
			ef(il, 0) += -1. * weight * ConvNormal * phiL(il) * solL[0]*solL[0]/fSolRef;
			for(jl=0; jl<nrowl; jl++) {
				ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl) *2.*solL[0]/fSolRef;
			}
		}
		for(ir=0; ir<nrowr; ir++) {
			ef(ir+nrowl,0) += -1. * (-1. * weight * ConvNormal * phiR(ir) * solL[0]*solL[0]/fSolRef);
			for(jl=0; jl<nrowl; jl++) {
				ek(ir+nrowl,jl) -= weight * ConvNormal * phiR(ir) * phiL(jl) * 2.*solL[0]/fSolRef;
			}
		}
	} else {
		for(ir=0; ir<nrowr; ir++) {
			ef(ir+nrowl,0) += -1. * (-1. * weight * ConvNormal * phiR(ir) * solR[0]*solR[0]/fSolRef );
			for(jr=0; jr<nrowr; jr++) {
				ek(ir+nrowl,jr+nrowl) -= weight * ConvNormal * phiR(ir) * phiR(jr) *2.*solR[0]/fSolRef;
			}
		}
		for(il=0; il<nrowl; il++) {
			ef(il,0) += -1. * weight * ConvNormal * phiL(il) * solR[0]*solR[0]/fSolRef;
			for(jr=0; jr<nrowr; jr++) {
				ek(il,jr+nrowl) += weight * ConvNormal * phiL(il) * phiR(jr) * 2.*solR[0]/fSolRef;
			}
		}
	}
	
	
	//diffusion term
	REAL leftK, rightK;
	leftK  = fK;
	rightK = fK;
	
	//Compute GradSol . normal
	REAL DSolLNormal = 0.;
	REAL DSolRNormal = 0.;
	for(id=0; id<fDim; id++) {
		DSolLNormal += dsolL(id,0)*normal[id];
		DSolRNormal += dsolR(id,0)*normal[id];
	}//for
	
	// 1) phi_I_left, phi_J_left
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		
		ef(il,0) += -1. * (weight * leftK * (this->fSymmetry * 0.5 * dphiLinormal*solL[0]-0.5*DSolLNormal*phiL(il,0)));
		
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(il,jl) += weight * leftK * (
										   this->fSymmetry * 0.5*dphiLinormal*phiL(jl,0)-0.5*dphiLjnormal*phiL(il,0)
										   );
		}
	}
	
	// 2) phi_I_right, phi_J_right
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		
		//ef = F - K u
		ef(ir+nrowl,0) += -1. * weight * rightK * ( this->fSymmetry * (-0.5 * dphiRinormal * solR[0] ) + 0.5 * DSolRNormal * phiR(ir) );
		
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += weight * rightK * (
														this->fSymmetry * (-0.5 * dphiRinormal * phiR(jr) ) + 0.5 * dphiRjnormal * phiR(ir)
														);
		}
	}
	
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL(id,il)*normal[id];
		}
		
		//ef = F - K u
		ef(il,0) += -1. * weight * ( this->fSymmetry * (-0.5 * dphiLinormal * leftK * solR[0] ) - 0.5 * DSolRNormal * rightK * phiL(il) );
		
		for(jr=0; jr<nrowr; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR(id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += weight * (
										 this->fSymmetry * (-0.5 * dphiLinormal * leftK * phiR(jr) ) - 0.5 * dphiRjnormal * rightK * phiL(il)
										 );
		}
	}
	
	// 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR(id,ir)*normal[id];
		}
		
		//ef = F - K u
		ef(ir+nrowl,0) += -1. * weight * (this->fSymmetry * 0.5 * dphiRinormal * rightK * solL[0] + 0.5 * DSolLNormal * leftK * phiR(ir));
		
		for(jl=0; jl<nrowl; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL(id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += weight * (
										 this->fSymmetry * 0.5 * dphiRinormal * rightK * phiL(jl) + 0.5 * dphiLjnormal * leftK * phiR(ir)
										 );
		}
	}
	
	if (this->IsSymetric()){
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
}

void TPZBurger::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                      REAL weight,
                                      TPZFMatrix<STATE> &ek,
                                      TPZFMatrix<STATE> &ef,
                                      TPZBndCond &bc) {
    STATE fK = AVGK();

    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	TPZVec<STATE> &solL=dataleft.sol[0];
	TPZFMatrix<STATE> &dsolL=dataleft.dsol[0];
	
	if (this->IsReferred()){
		this->SetConvectionTermInterface(dsolL, dsolL);
	}
	
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	REAL ConvNormal = 0.;
	for(id=0; id<fDim; id++) ConvNormal += fC * fConvDir[id]*normal[id];
	
	//Compute GradSol . normal
	REAL DSolLNormal = 0.;
	for(id=0; id<fDim; id++) {
		DSolLNormal += dsolL(id,0)*normal[id];
	}//for
	
	switch(bc.Type()) {
		case 0: // DIRICHLET
			
			//Diffusion
			for(il=0; il<nrowl; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL(id,il)*normal[id];
				}
				ef(il,0) += weight*fK*dphiLinormal*bc.Val2()(0,0) * this->fSymmetry;
				
				//ef = F - K u
				ef(il,0) += -1. * weight*fK*(this->fSymmetry * dphiLinormal * solL[0] - DSolLNormal * phiL(il,0));
				
				for(jl=0; jl<nrowl; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL(id,jl)*normal[id];
					}
					ek(il,jl) += weight*fK*(this->fSymmetry * dphiLinormal * phiL(jl,0) - dphiLjnormal * phiL(il,0));
				}
			}
			
			//Convection
			if(ConvNormal > 0.) {
				for(il=0; il<nrowl; il++) {
					ef(il,0) += -1. * weight * ConvNormal * phiL(il) * solL[0]*solL[0]/fSolRef;
					for(jl=0; jl<nrowl; jl++) {
						ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl) * 2.*solL[0]/fSolRef;
					}
				}
			} else {
				for(il=0; il<nrowl; il++) {
					ef(il,0) -= weight * ConvNormal * bc.Val2()(0,0) * phiL(il);
				}
			}
			
			break;
		case 1: // Neumann
			for(il=0; il<nrowl; il++) {
				ef(il,0) += weight*phiL(il,0)*bc.Val2()(0,0);
			}
			break;
			
		case 3: // outflow condition
			if(ConvNormal > 0.) {
				for(il=0; il<nrowl; il++) {
					for(jl=0; jl<nrowl; jl++) {
						ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl) *2.*solL[0]/fSolRef;
					}
					ef(il,0) += -1. * weight * ConvNormal * phiL(il) * solL[0]*solL[0]/fSolRef;
				}
			}
			else {
				if (ConvNormal < 0.){
					std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
				}
			}
			break;    
			
		default:
			PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
			break;
	}
    if (this->IsSymetric()){
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
}

int TPZBurger::ClassId() const{
    return Hash("TPZBurger") ^ TPZMatPoisson3dReferred::ClassId() << 1;
}

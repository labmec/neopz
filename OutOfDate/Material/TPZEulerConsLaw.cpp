/**
 * \file
 * @brief DEPRECATED FILE. This file contains implementations of the TPZEulerConsLawDEP methods
 */
#include "TPZEulerConsLaw.h" 
#include "TPZDiffusionConsLaw.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzreal.h"
#include <math.h>

using namespace std;
TPZEulerConsLawDEP::~TPZEulerConsLawDEP(){
	
}

TPZEulerConsLawDEP::TPZEulerConsLawDEP(int nummat,REAL delta_t,REAL gamma,int dim,const std::string &artdiff) :
TPZConservationLawDEP(nummat,delta_t,dim) {
	
	if( strcmp("SUPG",artdiff.c_str()) && strcmp("LS",artdiff.c_str()) && strcmp("Bornhaus",artdiff.c_str()) ){
		PZError << "TPZEulerConsLawDEP::TPZEulerConsLawDEP artificial diffusion parameter, default LS\n";
		fArtificialDiffusion = "LS";
	}
	else
	{
		fArtificialDiffusion = artdiff;//SUPG, LS, Bornhauss
	}
	fGamma = gamma;
}

TPZEulerConsLawDEP::TPZEulerConsLawDEP(TPZEulerConsLawDEP & copy) : TPZConservationLawDEP(copy){
	fArtificialDiffusion = copy.fArtificialDiffusion;
	fGamma = copy.fGamma;
}

TPZAutoPointer<TPZMaterial> TPZEulerConsLawDEP::NewMaterial(){
	TPZEulerConsLawDEP *result = new TPZEulerConsLawDEP(*this);
	return result;
}

void TPZEulerConsLawDEP::Print(std::ostream &out) {
	
	TPZDiffusionConsLaw *diff;
	out << "\nName of material : " << Name() << "\n";
	out << "Properties : \n";
	out << "Gamma : " << fGamma << endl;
	if(TPZDiffusionConsLaw::fCFL != 0.) out << "CFL  : " << TPZDiffusionConsLaw::fCFL << endl;
	else  out << "CFL  : 1/(2p+1) "  << endl;
	if(TPZDiffusionConsLaw::fDelta > 0) out << "Delta : " << TPZDiffusionConsLaw::fDelta << endl;
	else out << "Delta otimo : " << diff->DeltaOtimo() << endl;
	out << "Time step : " << TimeStep() << endl;
	out << "Difusao artificial : " << fArtificialDiffusion << endl;
	out << "Dimensao do problema : " << Dimension() << endl;
	out << "Numero de variaveis de estado : " << NStateVariables() << endl;
	out << "Numero de fluxos : " << NFluxes() << endl;
	
	TPZMaterial::Print(out);
}

void TPZEulerConsLawDEP::Contribute(TPZMaterialData &data, REAL weight,
                                 TPZFMatrix &ek,TPZFMatrix &ef) {
	
	TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	int phr = phi.Rows();// phi(in, 0) = phi_in  ,  dphi(i,jn) = dphi_jn/dxi
	int i,k,l,p,nstate = NStateVariables();//3, 4 ou 5
	if(fForcingFunction) {
		//na 2a itera�o deve-se ter fForcingFunction = 0
		TPZManVector<REAL> res(nstate);
		fForcingFunction(x,res);
		for(i=0;i<nstate;i++) sol[i] = res[i];
	}
	int dim = dphi.Rows();//dx, dy ou dz
	if(Dimension() != dim)
		PZError << "TPZEulerConsLawDEP::Contribute dimension error, dimension = " << dim << endl;
	
	//neste passo �calculada (�x/�,�y/�,�z/�) (no construtor)
	TPZDiffusionConsLaw diffusion(sol,fGamma,dim,fArtificialDiffusion);
	TPZVec<REAL> Fx(nstate),Fy(nstate),Fz(nstate);
	Flux(sol,Fx,Fy,Fz);
	TPZVec<REAL> gradphi(3,0.);
	TPZFMatrix Tx(nstate,nstate),Ty(nstate,nstate),Tz(nstate,nstate);
	TPZFMatrix DF1(nstate,nstate),DF2(nstate,nstate),DF3(nstate,nstate);
	diffusion.Tau(Tx,Ty,Tz);
	TPZFMatrix Trx(nstate,nstate),Try(nstate,nstate),Trz(nstate,nstate);
	int permite = 0;//SEM TRANSPOStA �A TEORIA
	if(permite){
		Tx.Transpose(&Trx);
		Ty.Transpose(&Try);
		Tz.Transpose(&Trz);
	}
	diffusion.GradientOfTheFlow(DF1,DF2,DF3);
	REAL timestep = TimeStep();
	REAL delta = diffusion.Delta(),Lpl,Hkp,Pkl;
	if(fDelta != 0.) delta = diffusion.DeltaOtimo();
	TPZFMatrix divF(nstate,1),prodpoint(nstate,nstate);
	TPZVec<REAL> sum1(nstate,0.),sum2(nstate,0.);
	
	for( int in = 0; in < phr; in++ ) {
		
		// w * Un
		for(i=0;i<nstate;i++) sum1[i] = phi(in, 0) * sol[i];
		
		// grad(w) . F
		for(i=0;i<nstate;i++){
			if(dim>0) sum2[i]  = Fx[i] * dphi(0,in);
			if(dim>1) sum2[i] += Fy[i] * dphi(1,in);
			if(dim>2) sum2[i] += Fz[i] * dphi(2,in);
		}
		
		//EF : w * Un + deltaT * (grad(w) . F)
		for(i=0;i<nstate;i++)
			ef(in * nstate + i, 0) += weight * (sum1[i] + timestep * sum2[i]);
		
		//EK
		// w * Un+1 + (grad(w) @ T) * div F(Un+1)    
		for( int jn = 0; jn < phr; jn++ ) {
			
			// DIFUS� + w * Un+1
			for(k=0;k<nstate;k++){
				for(l=0;l<nstate;l++){
					
					Pkl = 0.0;
					for(p=0;p<nstate;p++){
						
						if(permite){
							if(dim>0) {Hkp  = dphi(0,jn)*Trx(k,p); Lpl  = dphi(0,in)*DF1(p,l);}
							if(dim>1) {Hkp += dphi(1,jn)*Try(k,p); Lpl += dphi(1,in)*DF2(p,l);}
							if(dim>2) {Hkp += dphi(2,jn)*Trz(k,p); Lpl += dphi(2,in)*DF3(p,l);}
						} else {
							if(dim>0) {Hkp  = dphi(0,jn)*Tx(k,p); Lpl  = dphi(0,in)*DF1(p,l);}
							if(dim>1) {Hkp += dphi(1,jn)*Ty(k,p); Lpl += dphi(1,in)*DF2(p,l);}
							if(dim>2) {Hkp += dphi(2,jn)*Tz(k,p); Lpl += dphi(2,in)*DF3(p,l);}
						}
						Pkl += Hkp * Lpl;
						
					}
					
					// Dt * delta * (grad(w) o T) * div F (nstatex1)
					ek(nstate  * in + k, nstate  * jn + l) += weight * timestep * delta * Pkl;
					
					// w * Un+1 (nstatex1)
					if(l == k) ek(nstate * in + k, nstate * jn + l) += weight * phi(in,0) * phi(jn,0);
				}
			}
		}//jn
	}//in
}

void TPZEulerConsLawDEP::ContributeInterface(TPZMaterialData &data,
                                          REAL weight,
                                          TPZFMatrix &ek,
                                          TPZFMatrix &ef){
	// TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	// TPZFMatrix &phi = data.phi;
	TPZFMatrix &phiL = data.phil;
	TPZFMatrix &phiR = data.phir;
	TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
	TPZVec<REAL> &solL=data.soll;
	TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	int phrl = phiL.Rows();
	int phrr = phiR.Rows();
	int i,nstate = NStateVariables();//3, 4 ou 5
	if(fForcingFunction) {
		TPZManVector<REAL> res(nstate);// phi(in, 0) = phi_in , dphi(i,j) = dphi_j/dxi
		fForcingFunction(x,res);//fXf(0,0) = res[0];
		//si o elemento �BC ent� n� tem espao de interpola� associado => phi = {} (vacio) e a CC �aplicada
		if(phrl) for(i=0;i<nstate;i++) solL[i] = res[i];//solu�o inicial (interior, t=0) fForcingFunction != 0
		if(phrr) for(i=0;i<nstate;i++) solR[i] = res[i];//nas itera�es > 0 fForcingFunction = 0
	}//caso phr=0 a CC j�est�em sol
	
	//No que se segue a defini� das vari�eis pode ser otimizada
	REAL flux_rho,flux_rhou,flux_rhov,flux_rhow,flux_rhoE;
	TPZVec<REAL> flux_Roe(nstate,0.);
	
	if(nstate == 5){//dimens� 3
		TPZDiffusionConsLaw::Roe_Flux(solL[0],solL[1],solL[2],solL[3],solL[4],
									  solR[0],solR[1],solR[2],solR[3],solR[4],
									  normal[0],normal[1],normal[2],fGamma,
									  flux_rho,flux_rhou,flux_rhov,flux_rhow,flux_rhoE);
		flux_Roe[0] = flux_rho;
		flux_Roe[1] = flux_rhou;
		flux_Roe[2] = flux_rhov;
		flux_Roe[3] = flux_rhow;
		flux_Roe[4] = flux_rhoE;
		
	} else if(nstate == 4){//dimens� 2
		
		if(1){//TESTE : utiliza Roe 3D para o c�culo de Roe 2D
			TPZVec<REAL> left(2),right(2);
			left[0] = 0.0;
			left[1] = solL[3];
			right[0] = 0.0;
			right[1] = solR[3];          // ro     ro_u    ro_v    0.0     E
			TPZDiffusionConsLaw::Roe_Flux(solL[0],solL[1],solL[2],left[0],left[1],          //
										  solR[0],solR[1],solR[2],right[0],right[1],        // dimens� 1 : Roe 3D -> Roe 2D
										  normal[0],normal[1],normal[2],fGamma,             //
										  flux_rho,flux_rhou,flux_rhov,flux_rhow,flux_rhoE);//
		}
		if(0){
			TPZDiffusionConsLaw::Roe_Flux(solL[0],solL[1],solL[2],solL[3],
										  solR[0],solR[1],solR[2],solR[3],
										  normal[0],normal[1],fGamma,
										  flux_rho,flux_rhou,flux_rhov,flux_rhoE);
		}
		
		flux_Roe[0] = flux_rho;
		flux_Roe[1] = flux_rhou;
		flux_Roe[2] = flux_rhov;
		flux_Roe[3] = flux_rhoE;
		
	} else if(nstate == 3){//dimens� 1
		TPZVec<REAL> left(2),right(2);
		left[0] = 0.0;
		left[1] = solL[2];
		right[0] = 0.0;
		right[1] = solR[2];          // ro     ro_u    0.0     E
		TPZDiffusionConsLaw::Roe_Flux(solL[0],solL[1],left[0],left[1],
									  solR[0],solR[1],right[0],right[1],
									  normal[0],normal[1],fGamma,
									  flux_rho,flux_rhou,flux_rhov,flux_rhoE);
		flux_Roe[0] = flux_rho;
		flux_Roe[1] = flux_rhou;
		flux_Roe[2] = flux_rhoE;
	}
	//contribui� do elemento interface para a carga
	REAL timestep = TimeStep();
	int pos = 0;
	//a normal �interface sempre aponta do seu elemento esquerdo para o seu direito
	
	for( int in = 0; in < phrl; in++ ) {
		for(i=0;i<nstate;i++)
			ef(pos * nstate + i , 0) += -timestep * weight * phiL(in, 0) * flux_Roe[i];//left element
		pos++;
	}
	for( int in = 0; in < phrr; in++ ) {
		for(i=0;i<nstate;i++)
			ef(pos * nstate + i , 0) +=   timestep * weight * phiR(in, 0) * flux_Roe[i];//right element
		pos++;
	}
}


void TPZEulerConsLawDEP::ContributeBC(TPZMaterialData &data,REAL weight,
                                   TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {
	TPZFMatrix dphi = data.dphix;
	TPZFMatrix dphiL = data.dphixl;
	TPZFMatrix dphiR = data.dphixr;
	TPZFMatrix phi = data.phi;
	TPZFMatrix phiL = data.phil;
	TPZFMatrix phiR = data.phir;
	TPZManVector<REAL,3> normal = data.normal;
	TPZManVector<REAL,3> x = data.x;
	// int POrder=data.p;
	// int LeftPOrder=data.leftp;
	// int RightPOrder=data.rightp;
	TPZVec<REAL> sol=data.sol;
	TPZVec<REAL> solL=data.soll;
	TPZVec<REAL> solR=data.solr;
	TPZFMatrix dsol=data.dsol;
	TPZFMatrix dsolL=data.dsoll;
	TPZFMatrix dsolR=data.dsolr;
	// REAL faceSize=data.HSize;
	
	int phr = phi.Rows();
	short in,jn,i,j;
	int nstate = NStateVariables();
	REAL v2[5];//m�imo nstate
	for(i=0;i<nstate;i++) v2[i] = bc.Val2()(i,0);
	
	switch (bc.Type()) {
		case 0 :// Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				for(i = 0 ; i < nstate; i++)
					ef(in*nstate+i,0) += gBigNumber * weight * v2[i] * phi(in,0);
				for (jn = 0 ; jn < phr; jn++) {
					for(i = 0 ; i < nstate; i++)
						ek(in*nstate+i,jn*nstate+i) += gBigNumber * weight * phi(in,0) * phi(jn,0);
				}
			}
			break;
		case 1 :// Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {
				for(i = 0 ; i < nstate; i++)
					ef(in*nstate+i,0) += v2[i] * phi(in,0) * weight;
			}
			break;
		case 2 :// condi�o mista
			for(in = 0 ; in < phi.Rows(); in++) {
				for(i = 0 ; i < nstate; i++)
					ef(in*nstate+i, 0) += weight * v2[i] * phi(in, 0);
				for (jn = 0 ; jn < phi.Rows(); jn++) {
					for(i = 0 ; i < nstate; i++) for(j = 0 ; j < nstate; j++)
						ek(in*nstate+i,jn*nstate+j) += weight * bc.Val1()(i,j) * phi(in,0) * phi(jn,0);
				}
			}
	}
}

void TPZEulerConsLawDEP::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
							 TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
							 TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {
	
	cout << "\nTPZEulerConsLawDEP::Errors not implemented yet, program exit\n";
	exit(-1);
	
	//   TPZVec<REAL> sol(1),dsol(3);
	//   Solution(u,dudx,axes,1,sol);
	//   Solution(u,dudx,axes,2,dsol);
	//   REAL dx = dsol[0]*axes(0,0)+dsol[1]*axes(1,0)+dsol[2]*axes(2,0);
	//   REAL dy = dsol[0]*axes(0,1)+dsol[1]*axes(1,1)+dsol[2]*axes(2,1);
	//   REAL dz = dsol[0]*axes(0,2)+dsol[1]*axes(1,2)+dsol[2]*axes(2,2);
	//   //values[1] : eror em norma L2
	//   values[1]  = pow(sol[0] - u_exact[0],2.0);
	//   //values[2] : erro em semi norma H1
	//   values[2]  = pow(dx - du_exact(0,0),2.0);
	//   values[2] += pow(dy - du_exact(1,0),2.0);
	//   values[2] += pow(dz - du_exact(2,0),2.0);
	//   //values[0] : erro em norma H1 <=> norma Energia
	//   values[0]  = values[1]+values[2];
}

int TPZEulerConsLawDEP::NStateVariables() {
	return (2 + Dimension());//U = (ro, rou, rov, row, roe)
}

int TPZEulerConsLawDEP::NSolutionVariables(int var){
	
	if(var == 1 || var == 3 || var == 4 || var == 6) return 1;
	if(var == 2) return Dimension();
	if(var == 5) return NStateVariables();
	
	cout << "TPZEulerConsLawDEP::NSolutionVariables not defined\n";
	return 0;
}

/** returns the variable index associated with the name*/
int TPZEulerConsLawDEP::VariableIndex(const std::string &name) {
	if( !strcmp(name.c_str(),"density")  )     return 1;//rho
	if( !strcmp(name.c_str(),"velocity") )     return 2;//(u,v,w)
	if( !strcmp(name.c_str(),"energy")   )     return 3;//E
	if( !strcmp(name.c_str(),"pressure") )     return 4;//p
	if( !strcmp(name.c_str(),"solution") )     return 5;//(ro,u,v,w,E)
	if( !strcmp(name.c_str(),"normvelocity") ) return 6;//sqrt(u+v+w)
	cout << "TPZEulerConsLawDEP::VariableIndex not defined\n";
	return TPZMaterial::VariableIndex(name);
}

void TPZEulerConsLawDEP::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){
	
	if(fabs(Sol[0]) < 1.e-10) {
		cout << "\nTPZEulerConsLawDEP::Solution: Densidade quase nula\n";
		cout << "Densidade = " << Sol[0] << endl;
	}
	
	if(var == 1) {
		Solout.Resize(1);
		Solout[0] = Sol[0];//density
		return;
	} else if(var == 2) {
		int dim = Dimension();
		Solout.Resize(dim);
		for(int i=0;i<dim;i++) Solout[i] = Sol[i+1]/Sol[0];//velocity vector
		return;
	} else if(var == 3) {
		Solout.Resize(1);
		int pos = Dimension() + 1;
		Solout[0] = Sol[pos];//energy
		return;
	} else if(var == 4) {
		Solout.Resize(1);
		Solout[0] = Pressure(Sol);//pressure
		return;
	} else if(var == 5) {
		int nstate = NStateVariables();
		Solout.Resize(nstate);
		for(int i=0;i<nstate;i++) Solout[i] = Sol[i];//(ro,ro*u,ro*v,ro*w,E)
		return;
	} else if(var == 6) {
		int nstate = NStateVariables();
		Solout.Resize(1);
		REAL ro2 = Sol[0]*Sol[0];
		REAL veloc = 0.0;
		for(int i=1;i<nstate-1;i++) veloc += Sol[i]*Sol[i];//velocity vector
		Solout[0] = sqrt(veloc/ro2);
		return;
	} else {
		//cout << "TPZEulerConsLawDEP::Solution variable in the base class\n";
		TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
	}
}


REAL TPZEulerConsLawDEP::Pressure(TPZVec<REAL> &U) {
	
	if(fabs(U[0]) < 1.e-6) {
		cout << "\nTPZEulerConsLawDEP::Pressure: Densidade quase nula\n";
		cout << "Densidade = " << U[0] << endl;
		//exit(-1);
	}
	// Press� = (gam-1)*(E - ro*||(u,v,w)||/2)
	// onde aqui ro_e = E (nota�)
	
	int nstate = NStateVariables();
	REAL press = 0.0;
	if(U.NElements() != nstate) U.Resize(nstate);
	if(nstate == 5){    
		//U = (U0,U1,U2,U3,U4) = (ro , ro u , ro v , ro w , ro e)
		REAL rho_velocity = ( U[1]*U[1] + U[2]*U[2] + U[3]*U[3] )/U[0];
		press = ((fGamma-1.)*( U[4] - 0.5 * rho_velocity ));
	} else
		if(nstate == 4){    
			//U = (U0,U1,U2,U3,U4) = (ro , ro u , ro v , ro e)
			REAL rho_velocity = ( U[1]*U[1] + U[2]*U[2] )/U[0];
			press = ((fGamma-1.)*( U[3] - 0.5 * rho_velocity ));
		} else
			if(nstate == 3){
				//U = (U0,U1,U2,U3,U4) = (ro , ro u , ro e)
				REAL rho_velocity = ( U[1]*U[1] )/U[0];
				press = ((fGamma-1.)*( U[2] - 0.5 * rho_velocity ));
			} else { 
				cout << "\nTPZEulerConsLawDEP::Pressure caso nao tratado retorna nulo\n";
				return 0.0;
			}
	if(press < 0){
		PZError << "TPZEulerConsLawDEP::Pressure pressao negativa: pressao = " << press << endl;
		press = (fGamma-1.)*U[nstate-1];
		PZError << "TPZEulerConsLawDEP::Pressure pressao substituta (gama-1)*E = " << press << endl;
	}
	return press;
}

void TPZEulerConsLawDEP::Flux(TPZVec<REAL> &U,TPZVec<REAL> &Fx,TPZVec<REAL> &Fy,TPZVec<REAL> &Fz) {
	
	REAL press = Pressure(U);
	int nstate = NStateVariables();
	if(nstate < 3 && nstate > 5){
		cout << "TPZEulerConsLawDEP::Flux case not implemented\n";
		Fx.Resize(0);
		Fy.Resize(0);
		Fz.Resize(0);
		exit(-1);
	}
	
	Fx.Resize(5,0.0);//v�ido
	Fy.Resize(5,0.0);//para
	Fz.Resize(5,0.0);//R , R , R
	
	if(nstate == 5){
		Fx[0] = U[1];//ro u
		Fx[1] = (U[1]/U[0])*U[1] + press;//ro u2 + p
		Fx[2] = U[1]*(U[2]/U[0]);//ro u v
		Fx[3] = U[1]*(U[3]/U[0]);//ro u w
		Fx[4] = (U[4]+press)*(U[1]/U[0]);//(ro e + p) u
		
		Fy[0] = U[2];//ro v
		Fy[1] = U[2]*(U[1]/U[0]);//ro v u
		Fy[2] = (U[2]/U[0])*U[2] + press;//ro v2 + p
		Fy[3] = U[2]*(U[3]/U[0]);//ro v w
		Fy[4] = (U[4] + press)*(U[2]/U[0]);//(ro e + p) v
		
		Fz[0] = U[3];//ro w
		Fz[1] = U[3]*(U[1]/U[0]);//ro w u
		Fz[2] = U[3]*(U[2]/U[0]);//ro w v
		Fz[3] = (U[2]/U[0])*U[2] + press;//ro w2 + p  
		Fz[4] = (U[4] + press)*(U[2]/U[0]);//(ro e + p) w
		return;
	}
	
	if(nstate == 4){
		Fx[0] = U[1];//ro u
		Fx[1] = U[1]*U[1] / U[0] + press;//ro u2 + p
		Fx[2] = U[1]*U[2] / U[0];//ro u v
		Fx[3] = (U[3]+press)*U[1] / U[0];//(E + p) u
		
		Fy[0] = U[2];//ro v
		Fy[1] = U[1]*U[2] / U[0];//ro u v
		Fy[2] = U[2]*U[2] / U[0] + press;//ro v2 + p
		Fy[3] = (U[3] + press)*U[2] / U[0];//(E + p) v
		return;
	}
	
	if(nstate == 3){
		Fx[0] = U[1];//ro u
		Fx[1] = (U[1]/U[0])*U[1] + press;//ro u2 + p
		Fx[2] = (U[2]+press)*(U[1]/U[0]);//(ro e + p) u
	}
}

void TPZEulerConsLawDEP::Flux(TPZVec<REAL> &x,TPZVec<REAL> &Sol,TPZFMatrix &DSol,
						   TPZFMatrix &axes,TPZVec<REAL> &flux) {
	TPZVec<REAL> Fx,Fy,Fz;
	Flux(Sol,Fx,Fy,Fz);
	int cap = Sol.NElements();
	int nstate = NStateVariables(),i;
	if(cap != nstate){
		PZError << "\nTPZEulerConsLawDEP::Flux data size error\n";
		flux.Resize(0);
		return;
	}
	if(nstate == 3){
		flux.Resize(3);
		for(i=0;i<3;i++) flux[i] = Fx[i];
		return;
	} else
		if(nstate == 4){
			flux.Resize(8);
			for(i=0;i<4;i++) flux[i] = Fx[i];
			for(i=4;i<8;i++) flux[i] = Fy[i];
			return;
		} else
			if(nstate == 5){
				flux.Resize(15);
				for(i=00;i<05;i++) flux[i] = Fx[i];
				for(i=05;i<10;i++) flux[i] = Fy[i];
				for(i=10;i<15;i++) flux[i] = Fz[i];
			}
}

void TPZEulerConsLawDEP::Contribute(TPZMaterialData &data,
                                 REAL weight,TPZFMatrix &ef) {
	
	TPZFMatrix dphi = data.dphix;
	TPZFMatrix dphiL = data.dphixl;
	TPZFMatrix dphiR = data.dphixr;
	TPZFMatrix phi = data.phi;
	TPZFMatrix phiL = data.phil;
	TPZFMatrix phiR = data.phir;
	TPZManVector<REAL,3> normal = data.normal;
	TPZManVector<REAL,3> x = data.x;
	// int POrder=data.p;
	// int LeftPOrder=data.leftp;
	// int RightPOrder=data.rightp;
	TPZVec<REAL> sol=data.sol;
	TPZVec<REAL> solL=data.soll;
	TPZVec<REAL> solR=data.solr;
	TPZFMatrix dsol=data.dsol;
	TPZFMatrix dsolL=data.dsoll;
	TPZFMatrix dsolR=data.dsolr;
	// REAL faceSize=data.HSize;
	
	int phr = phi.Rows();// phi(in, 0) = phi_in  ,  dphi(i,jn) = dphi_jn/dxi
	int i,nstate = NStateVariables();//3, 4 ou 5
	if(fForcingFunction) {
		//na 2a itera�o deve-se ter fForcingFunction = 0
		TPZManVector<REAL> res(nstate);
		fForcingFunction(x,res);
		for(i=0;i<nstate;i++) sol[i] = res[i];
	}
	int dim = dphi.Rows();//dx, dy ou dz
	if(Dimension() != dim)
		PZError << "TPZEulerConsLawDEP::Contribute dimension error, dimension = " << dim << endl;
	
	//neste passo �calculada (�x/�,�y/�,�z/�) (no construtor)
	TPZDiffusionConsLaw diffusion(sol,fGamma,dim,fArtificialDiffusion);
	TPZVec<REAL> Fx(nstate),Fy(nstate),Fz(nstate);
	Flux(sol,Fx,Fy,Fz);
	TPZVec<REAL> gradphi(3,0.);
	TPZFMatrix Tx(nstate,nstate),Ty(nstate,nstate),Tz(nstate,nstate);
	TPZFMatrix DF1(nstate,nstate),DF2(nstate,nstate),DF3(nstate,nstate);
	diffusion.Tau(Tx,Ty,Tz);
	TPZFMatrix Trx(nstate,nstate),Try(nstate,nstate),Trz(nstate,nstate);
	Tx.Transpose(&Trx);
	Ty.Transpose(&Try);
	Tz.Transpose(&Trz);
	diffusion.GradientOfTheFlow(DF1,DF2,DF3);
	REAL timestep = TimeStep();
	//REAL delta = diffusion.Delta();
	//if(fDelta!= 0.) delta = diffusion.DeltaOtimo();
	TPZFMatrix divF(nstate,1),prodpoint(nstate,nstate);
	TPZVec<REAL> sum1(nstate,0.),sum2(nstate,0.);
	
	for( int in = 0; in < phr; in++ ) {
		
		// w * Un
		for(i=0;i<nstate;i++) sum1[i] = phi(in, 0) * sol[i];
		
		// grad(w) . F
		for(i=0;i<nstate;i++){
			if(dim>0) sum2[i]  = Fx[i] * dphi(0,in);
			if(dim>1) sum2[i] += Fy[i] * dphi(1,in);
			if(dim>2) sum2[i] += Fz[i] * dphi(2,in);
		}
		
		//EF : w * Un + deltaT * (grad(w) . F)
		for(i=0;i<nstate;i++)
			ef(in * nstate + i, 0) += weight * (sum1[i] + timestep * sum2[i]);
	}//in
}

//m�odo para testes de implementa�, verifica� das integrais elementares
void TPZEulerConsLawDEP::ContributeTESTE(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,
									  REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,
									  TPZFMatrix &ek,TPZFMatrix &ef) {
	
	//   est�chamada deve ser feita do Contribute(..)
	//   if(0){
	//     ContributeTESTE(x,jacinv,sol,dsol,weight,axes,phi,dphi,ek,ef);//PARA TESTES
	//     return;
	//   }
	
	static int firsttime = 1;
	static REAL EK[2],EF[3];
	if(firsttime && 0){
		EF[0]=1.; EF[1]=1.; EF[2]=1.;
		cout << "TPZEulerConsLawDEP::Contribute EK1+EK2: [0:ambas][1:EK 1][2:EK 2]";
		int par;
		cin >> par;
		if(par==0){EK[0] = 1.0; EK[1] = 1.0;} 
		if(par==1){EK[0] = 1.0; EK[1] = 0.0;}
		if(par==2){EK[0] = 0.0; EK[1] = 1.0;}
		firsttime = 0;
	}
	if(firsttime && 1){
		EK[0] = 1.0; EK[1] = 1.0;
		cout << "TPZEulerConsLawDEP::Contribute EF1+EF2+EF3:"
		<< "[0:todas][1:EF 1][2:EF 2][3:EF 3]";
		int par;
		cin >> par;
		if(par==0){EF[0] = 1.0; EF[1] = 1.0; EF[2] = 1.0;} 
		if(par==1){EF[0] = 1.0; EF[1] = 0.0; EF[2] = 0.0;}
		if(par==2){EF[0] = 0.0; EF[1] = 1.0; EF[2] = 0.0;}
		if(par==3){EF[0] = 0.0; EF[1] = 0.0; EF[2] = 1.0;}
		firsttime = 0;
	}
	
	//int phr = phi.Rows();// phi(in, 0) = phi_in  ,  dphi(i,jn) = dphi_jn/dxi
	int nstate = NStateVariables();//3, 4 ou 5
	if(fForcingFunction) {
		//na 2a itera�o deve-se ter fForcingFunction = 0
		TPZManVector<REAL> res(nstate);
		fForcingFunction(x,res);
		for(int i=0;i<nstate;i++) sol[i] = res[i];
	}
	int dim = dphi.Rows();//dx, dy ou dz
	if(Dimension() != dim)
		PZError << "TPZEulerConsLawDEP::Contribute dimension error, dimension = " << dim << endl;
	
	//   TPZDiffusionConsLaw diffusion(sol,fGamma,dim,fArtificialDiffusion);
	//   TPZVec<REAL> Fx(nstate),Fy(nstate),Fz(nstate);
	//   Flux(sol,Fx,Fy,Fz);
	//   TPZVec<REAL> gradphi(3,0.),prodpoint(nstate),divF(nstate);
	//   diffusion.GradientOfTheFlow(DF1,DF2,DF3);
	//   REAL timestep = TimeStep();
	//   //REAL delta = diffusion.Delta();
	//   REAL delta = diffusion.DeltaOtimo(),Lpl,Hkp,Pkl;
	//   TPZFMatrix divF(nstate,1),prodpoint(nstate,nstate);
	//   TPZVec<REAL> sum1(nstate,0.),sum2(nstate,0.);
	
	//   for( int in = 0; in < phr; in++ ) {
    
	//     // w * Un
	//     for(i=0;i<nstate;i++) sum1[i] = phi(in, 0) * sol[i];
	
	//     // grad(w) . F
	//     for(i=0;i<nstate;i++){
	//       if(dim>0) sum2[i]  = Fx[i] * dphi(0,in);
	//       if(dim>1) sum2[i] += Fy[i] * dphi(1,in);
	//       if(dim>2) sum2[i] += Fz[i] * dphi(2,in);
	//     }
	
	//     //EF : w * Un + deltaT * (grad(w) . F)
	//     for(i=0;i<nstate;i++)
	//       ef(in * nstate + i, 0) += weight * (EF[0] * sum1[i] + EF[1] * timestep * sum2[i]);
	
	//     //EK
	//     // w * Un+1 + (grad(w) @ T) * div F(Un+1)    
	//     for( int jn = 0; jn < phr; jn++ ) {
	
	//       // DIFUS� + w * Un+1
	//       for(k=0;k<nstate;k++){
	// 	for(l=0;l<nstate;l++){
	
	// 	  Pkl = 0.0;
	// 	  for(p=0;p<nstate;p++){
	
	// 	    if(dim>0) {Hkp  = dphi(0,jn)*Tx(k,p); Lpl  = dphi(0,in)*DF1(p,l);}
	// 	    if(dim>1) {Hkp += dphi(1,jn)*Ty(k,p); Lpl += dphi(1,in)*DF2(p,l);}
	// 	    if(dim>2) {Hkp += dphi(2,jn)*Tz(k,p); Lpl += dphi(2,in)*DF3(p,l);}
	// 	    Pkl += Hkp * Lpl;
	
	// 	  }
	
	// 	  // Dt * delta * (grad(w) o T) * div F (nstatex1)
	// 	  ek(nstate  * in + k, nstate  * jn + l) += EK[1] * weight * timestep * delta * Pkl;
	
	// 	  // w * Un+1 (nstatex1)
	// 	  if(l == k) ek(nstate * in + k, nstate * jn + l) += EK[0] * weight * phi(in,0) * phi(jn,0);
	// 	}
	//       }
	//     }//jn
	//   }//in
}

void TPZEulerConsLawDEP::ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft){
	
	if(!bcleft){
		PZError << "TPZEulerConsLawDEP::ComputeSolLeft null bundary condition return\n";
		return;
	}
	int i,nstate = NStateVariables();
	REAL vpn=0.;
	switch (bcleft->Type()){
		case 0://Dirichlet
		case 1://Neumann
		case 2://Mista
			PZError << "TPZEulerConsLawDEP::ComputeSolLeft boundary condition error\n";
			break;
		case 3://Dirichlet: nada a fazer a CC �a correta
			for(i=0;i<nstate;i++) soll[i] = bcleft->Val2()(i,0);
			break;
		case 4://recuperar valor da solu� MEF direita: saida livre
			for(i=0;i<nstate; i++) soll[i] = solr[i];
			break;
		case 5://parede
			for(i=1;i<nstate-1;i++) vpn += solr[i]*normal[i-1];//v.n
			for(i=1;i<nstate-1;i++) soll[i] = solr[i] - 2.0*vpn*normal[i-1];
			soll[0] = solr[0];
			soll[nstate-1] = solr[nstate-1];
			break;
		case 6://n� refletivas (campo distante)
			for(i=0;i<nstate;i++) soll[i] = solr[i];
			break;
		default:
			PZError << "TPZEulerConsLawDEP::ContributeInterface Boundary Condition Type Not Exists\n";
	}
}


void TPZEulerConsLawDEP::ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright){
	
	if(!bcright){
		PZError << "TPZEulerConsLawDEP::ComputeSolLeft null bundary condition return\n";
		return;
	}
	int i,nstate = NStateVariables();
	REAL vpn=0.;
	switch (bcright->Type()){
		case 0://Dirichlet
		case 1://Neumann
		case 2://Mista
			PZError << "TPZEulerConsLawDEP::ComputeSolLeft boundary condition error\n";
			break;
		case 3://Dirichlet: nada a fazer a CC �a correta
			break;
		case 4://recuperar valor da solu� MEF esquerda: saida livre
			for(i=0;i<nstate; i++) solr[i] = soll[i];
			break;
		case 5://condi� de parede
			for(i=1;i<nstate-1;i++) vpn += soll[i]*normal[i-1];//v.n
			for(i=1;i<nstate-1;i++) solr[i] = soll[i] - 2.0*vpn*normal[i-1];
			solr[0] = soll[0];
			solr[nstate-1] = soll[nstate-1];
			break;
		case 6://n� refletivas (campo distante)
			for(i=0;i<nstate;i++) solr[i] = soll[i];
			break;
		default:
			PZError << "TPZEulerConsLawDEP::ContributeInterface Boundary Condition Type Not Exists\n";
	}
}

void TPZEulerConsLawDEP::SetDeltaTime(REAL maxveloc,REAL deltax,int degree){
	
	REAL CFL = 1./((2.0*(REAL)degree) + 1.0);
	//  TPZDiffusionConsLaw *diff;
	//  if(diff->fCFL) CFL = diff->fCFL;//o primeiro valor n� nulo �mantido
	REAL deltaT = CFL*deltax/maxveloc;
	cout << "TPZCompMesh::Delta Time : " << deltaT << endl;
	SetTimeStep(deltaT);
}

void TPZEulerConsLawDEP::TestOfRoeFlux(REAL &tetainit,REAL &tetamax,REAL &tol,REAL &increment){
	//Problema teste choque refletido estacion�io de tr� estados constantes 
	//OS VALORES ENCONTRADOS S�:
	//61 GRAUS 
	//-66.7169 GRAUS (-1.16443 RADIANOS)
	//regi� R1
	//  REAL r1M = 2.9;//Mach
	REAL r1ro = 1.0;
	REAL r1u = 2.9;
	REAL r1v = 0.0;
	REAL r1w = 0.0;
	REAL r1p = 0.714286;
	REAL r1vel2 = r1u*r1u+r1v*r1v+r1w*r1w;
	//regi� R2
	//REAL r2M = 2.3781;//Mach
	REAL r2ro = 1.7;
	REAL r2u = 2.61934;
	REAL r2v = -0.50632;
	REAL r2w = 0.0;
	REAL r2p = 1.52819;
	REAL r2vel2 = r2u*r2u+r2v*r2v+r2w*r2w;
	//regi� R3
	//REAL r3M = 1.94235;//Mach
	REAL r3ro = 2.68728;
	REAL r3u = 2.40140;
	REAL r3v = 0.0;
	REAL r3w = 0.0;
	REAL r3p = 2.93407;
	REAL r3vel2 = r3u*r3u+r3v*r3v+r3w*r3w;
	//procurando a normal �frente estacionaria: aproximadamente 23.28 graus ~ 0.406313 radianos
	//entre as regi�s R1 e R2
	//reta (cos  , sen ) de �gulo 
	//(sen  , -cos ) normal �reta de �gulo  apontando da regi� R1 para a regi� R2
	REAL teta=0.0,flux_rho,flux_rhou,flux_rhov,flux_rhoE,gama = 1.4;//flux_rhow,
	//REAL increment=0.00001;
	TPZVec<REAL> U1(5,0.),U2(5,0.),U3(5,0.),r1Fx(0),r1Fy(0),r1Fz(0),r2Fx(0),r2Fy(0),r2Fz(0),r3Fx(0),r3Fy(0),r3Fz(0),n(3,0.);
	int i,enter=0;
	U1[0] = r1ro;
	U1[1] = r1ro*r1u;
	U1[2] = r1ro*r1v;
	U1[3] = r1p/(gama-1.0)+r1ro*r1vel2*0.5;
	U1[4] = 0.0;//
	U2[0] = r2ro;
	U2[1] = r2ro*r2u;
	U2[2] = r2ro*r2v;
	U2[3] = r2p/(gama-1.0)+r2ro*r2vel2*0.5;
	U2[4] = 0.0;//
	U3[0] = r3ro;
	U3[1] = r3ro*r3u;
	U3[2] = r3ro*r3v;
	U3[3] = r3p/(gama-1.0)+r3ro*r3vel2*0.5;
	U3[4] = 0.0;
	//REAL tetamax = 2.0*asin(1.0)+0.1;
	REAL soma = 1.0;
	REAL suma = 1.0;
	teta = tetainit;
	while(teta < tetamax){
		teta += increment;
		n[0] = cos(teta);
		n[1] = sin(teta);
		Flux(U1,r1Fx,r1Fy,r1Fz);
		Flux(U2,r2Fx,r2Fy,r2Fz);
		Flux(U3,r3Fx,r3Fy,r3Fz);
		soma = 0.0;
		suma = 0.0;
		for(i=0;i<4;i++){
			soma += (r1Fx[i] - r2Fx[i])*n[0] + (r1Fy[i] - r2Fy[i])*n[1];
		}
		for(i=0;i<4;i++){
			suma += (r2Fx[i] - r3Fx[i])*n[0] + (r2Fy[i] - r3Fy[i])*n[1];
		}
		if(fabs(soma) < tol){
			cout << "TPZEulerConsLawDEP::TestOfRoeFlux found angle, angle = " << teta << "\tradians\n";
			//a normal aponta da regi� R1 para a regi� R2: U1 �esquerdo e U2 direito
			TPZDiffusionConsLaw::Roe_Flux(U1[0],U1[1],U1[2],U1[3],
										  U2[0],U2[1],U2[2],U2[3],
										  n[0],n[1],gama,
										  flux_rho,flux_rhou,flux_rhov,flux_rhoE);
			cout << "TPZEulerConsLawDEP::TestOfRoeFlux flow in the R1 region\n";
			cout << "density : " << r1Fx[0]*n[0]+r1Fy[0]*n[1] << endl
			<< "u*ro    : " << r1Fx[1]*n[0]+r1Fy[1]*n[1] << endl
			<< "v*ro    : " << r1Fx[2]*n[0]+r1Fy[2]*n[1] << endl
			<< "energy  : " << r1Fx[3]*n[0]+r1Fy[3]*n[1] << endl << endl;
			cout << "TPZEulerConsLawDEP::TestOfRoeFlux flow in the R2 region\n";
			cout << "density : " << r2Fx[0]*n[0]+r2Fy[0]*n[1] << endl
			<< "u*ro    : " << r2Fx[1]*n[0]+r2Fy[1]*n[1] << endl
			<< "v*ro    : " << r2Fx[2]*n[0]+r2Fy[2]*n[1] << endl
			<< "energy  : " << r2Fx[3]*n[0]+r2Fy[3]*n[1] << endl << endl;
			cout << "TPZEulerConsLawDEP::TestOfRoeFlux the calculation of the flow of Roe is\n";
			cout << "density : " << flux_rho << endl
			<< "u*ro    : " << flux_rhou << endl
			<< "v*ro    : " << flux_rhov << endl
			<< "energy  : " << flux_rhoE << endl << endl;
			cout << "\nTPZEulerConsLawDEP::TestOfRoeFlux norma da diferenca |F1*n - F2*n| = " << fabs(soma) << "\n\n";
			enter = 1;
		}
		if(fabs(suma) < tol){
			cout << "TPZEulerConsLawDEP::TestOfRoeFlux found angle, angle = " << teta << "\tradians\n";
			//a normal aponta da regi� R1 para a regi� R2: 
			TPZDiffusionConsLaw::Roe_Flux(U2[0],U2[1],U2[2],U2[3],
										  U3[0],U3[1],U3[2],U3[3],
										  n[0],n[1],gama,
										  flux_rho,flux_rhou,flux_rhov,flux_rhoE);
			cout << "TPZEulerConsLawDEP::TestOfRoeFlux flow in the R1 region\n";
			cout << "density : " << r2Fx[0]*n[0]+r2Fy[0]*n[1] << endl
			<< "u*ro    : " << r2Fx[1]*n[0]+r2Fy[1]*n[1] << endl
			<< "v*ro    : " << r2Fx[2]*n[0]+r2Fy[2]*n[1] << endl
			<< "energy  : " << r2Fx[3]*n[0]+r2Fy[3]*n[1] << endl << endl;
			cout << "TPZEulerConsLawDEP::TestOfRoeFlux flow in the R2 region\n";
			cout << "density : " << r3Fx[0]*n[0]+r3Fy[0]*n[1] << endl
			<< "u*ro    : " << r3Fx[1]*n[0]+r3Fy[1]*n[1] << endl
			<< "v*ro    : " << r3Fx[2]*n[0]+r3Fy[2]*n[1] << endl
			<< "energy  : " << r3Fx[3]*n[0]+r3Fy[3]*n[1] << endl << endl;
			cout << "TPZEulerConsLawDEP::TestOfRoeFlux the calculation of the flow of Roe is\n";
			cout << "density : " << flux_rho << endl
			<< "u*ro    : " << flux_rhou << endl
			<< "v*ro    : " << flux_rhov << endl
			<< "energy  : " << flux_rhoE << endl << endl;
			cout << "\nTPZEulerConsLawDEP::TestOfRoeFlux norma da diferenca |F2*n - F3*n| = " << fabs(suma) << "\n\n";
			enter = 1;
			//break;
			//       if(fabs(soma) > 1.e-9){
			// 	tetamax = teta + 10.0*increment;
			// 	tetainit = teta - 10.0*increment;
			// 	tol /= 10.0;
			// 	increment /=1000.0;
			// 	TestOfRoeFlux(tetainit,tetamax,tol,increment);
			//       }
		}
	}
	if(enter) cout << "\nTPZEulerConsLawDEP::TestOfRoeFlux angle not found, The End\n";
	else cout << "\nTPZEulerConsLawDEP::TestOfRoeFlux norma da diferenca |F1*n - F2*n| = " << fabs(soma) << "\n\n";
	cout << "\nTPZEulerConsLawDEP::TestOfRoeFlux concluded test\n";
	//OS VALORES ENCONTRADOS S�:
	//61 GRAUS 
	//-66.7169 GRAUS (-1.16443 RADIANOS)
}

//       if(diffright){// DIFUS� NA CARGA : 2
// 	// w * Un+1 
// 	for(i=0;i<nstate;i++)
// 	  ek(nstate * in + i, nstate * jn + i) += weight * phi(in,0) * phi(jn,0);
//       }

//     if(diffright){// DIFUS� NA CARGA : 1
//       // grad(w)
//       for(i=0;i<dim;i++) gradphi[i] = dphi(i,in);
//       // grad(w) @ T (point product)
//       diffusion.PointOperator(gradphi,prodpoint);
//       for(i=0;i<dim;i++) gradphi[i] = dsol(i,0);
//       diffusion.Divergence(gradphi,divF);
//       TPZVec<REAL> diff_term(nstate,0.);
//       for(i=0;i<nstate;i++) for(j=0;j<nstate;j++) diff_term[i] += prodpoint(i,j) * divF(j,0);
//       for(i=0;i<nstate;i++)
// 	 ef(in * nstate + i, 0) += weight * timestep * delta * diff_term[i];
//     }


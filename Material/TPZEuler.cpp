/**
 * \file
 * @brief Contains implementations of the TPZEuler methods.
 */

#include "TPZEuler.h"
#include "pzerror.h"
#include "pzbndcond.h"

using namespace std;
TEulerDiffusivity TPZEuler::gEul;

void TPZEuler::SetData(istream &data) {
	TPZMaterial::SetData(data);
	data >> fDeltaT;
}
TPZAutoPointer<TPZMaterial> TPZEuler::NewMaterial() {
	TPZEuler *result = new TPZEuler(*this);
	
	return result;
}
void TPZEuler::Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,
						TPZVec<REAL> &Solout){
	if(var == 1) {
		Solout.Resize(1);
		Solout[0] = gEul.Pressure(Sol);
	} else if(var ==2) {
		Solout.Resize(1);
		Solout[0] = Sol[0];
	} else if(var == 3) {
		Solout.Resize(2);
		Solout[0] = Sol[1];
		Solout[1] = Sol[2];
	} else TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}
int TPZEuler::NSolutionVariables(int index) {
	if(index == 1) return 1;
	if(index == 2) return 1;
	if(index == 3) return 2;
	return TPZMaterial::NSolutionVariables(index);
}
int TPZEuler::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"pressure")) return 1;
	if(!strcmp(name.c_str(),"density")) return 2;
	if(!strcmp(name.c_str(),"velocity")) return 3;
	return TPZMaterial::VariableIndex(name);
}
void TPZEuler::Print(std::ostream & out) {
    TPZMaterial::Print(out);
}
void TPZEuler::ContributeBC(TPZMaterialData &data,REAL weight,
							TPZFMatrix<REAL> &ek,
							TPZFMatrix<REAL> &ef,TPZBndCond &bc) {
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
	TPZFMatrix<REAL> &axes=data.axes;
	
	if(fState == 0) return;
	if(bc.Material().operator ->() != this){
		PZError << "TPZMat1dLin.apply_bc warning : this material didn't create the boundary condition!\n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 3){
		PZError << "TPZEuler.aplybc, unknown boundary condition type :"  <<
		bc.Type() << " boundary condition ignored\n";
		return;
	}
	
	
	int numdof = NStateVariables();
	int numnod = ek.Rows()/numdof;
	int r = numdof;
	
	TPZVec<REAL> flux(8);
	gEul.Flux(sol,flux);
	REAL normal[2] = {axes(0,1),-axes(0,0)};
	/*
	 int i;
	 cout << " flux = " << x[0] << ' ' << x[1] << ' ' << "normal" << normal[0] << ' ' << normal[1] << ' ';
	 for(i=0; i<4; i++) cout << flux[i+4] << ' ';
	 cout << endl;
	 */
	int idf,jdf,in,jn;
	switch(bc.Type()){
			
		case 0:
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += gBigNumber*phi(in,0)*bc.Val2()(idf,0)*weight;
				}
				for(jn=0 ; jn<numnod ; ++jn) {
					for(idf = 0;idf<r;idf++) {
						ek(in*r+idf,jn*r+idf) += gBigNumber*phi(in,0)*phi(jn,0)*weight;
					}
				}
			}
			break;
			
		case 1:
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += phi(in,0)*bc.Val2()(idf,0)*weight;
				}
			}
			break;
			
		case 2:
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += phi(in,0)*bc.Val2()(idf,0)*weight;
				}
				for(jn=0 ; jn<numnod ; ++jn) {
					for(idf = 0;idf<r;idf++) {
						for(jdf = 0;jdf<r;jdf++) {
							ek(in*r+idf,jn*r+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
						}
					}
				}
			case 3: {
				TPZFMatrix<REAL> A(4,4),B(4,4);
				gEul.JacobFlux(sol,A,B);
				for(in=0; in<numnod; in++) {
					for(idf=0; idf<4; idf++) {
						/*
						 ef(4*in+idf) += weight*fDeltaT*(
						 -phi(in,0)*flux[idf]*normal[0]
						 -phi(in,0)*flux[idf+4]*normal[1]
						 );
						 */
						for(jn=0; jn<numnod; jn++) {
							for(jdf=0; jdf<4; jdf++) {      
								ek(4*in+idf,4*jn+jdf) += weight*fDeltaT*
								(phi(in,0)*A(idf,jdf)*phi(jn,0)*normal[0]
								 + phi(in,0)*B(idf,jdf)*phi(jn,0)*normal[1]);
								
							}
						}
					}
				}
			}
			}//fim switch
	}
}
void TPZEuler::Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<REAL> &ek,
						  TPZFMatrix<REAL> &ef) {
	TPZFMatrix<REAL> &dphi = data.dphix;
	// TPZFMatrix<REAL> &dphiL = data.dphixl;
	// TPZFMatrix<REAL> &dphiR = data.dphixr;
	TPZFMatrix<REAL> &phi = data.phi;
	// TPZFMatrix<REAL> &phiL = data.phil;
	// TPZFMatrix<REAL> &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
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
	TPZFMatrix<REAL> &daxesdksi=data.jacinv;
	TPZFMatrix<REAL> &axes=data.axes;
	
    int nshape = phi.Rows();
    REAL dphix[2];
    int in,jn,idf,jdf;
    for(in=0; in<nshape; in++) {
		dphix[0] = dphi(0,in)*axes(0,0)+dphi(1,in)*axes(1,0);
		dphix[1] = dphi(1,in)*axes(0,1)+dphi(1,in)*axes(1,1);
		dphi(0,in) = dphix[0];
		dphi(1,in) = dphix[1];
    }
	
    if(fState == 0 && ! fForcingFunction) {
		cout << "TPZEuler failed to suply initial conditions\n";
		return;
    }
    if(fState == 0) {
		TPZVec<REAL> force(4);
		fForcingFunction->Execute(x,force);
		for(in=0; in<nshape; in++) {
			for(idf=0; idf<4; idf++) {
				ef(4*in+idf) += weight*phi(in,0)*gBigNumber*force[idf];
				for(jn=0; jn<nshape; jn++) {
					ek(4*in+idf,4*jn+idf) += weight*phi(in,0)*phi(jn,0)*gBigNumber;
				}
			}
		}
		return;
    }
    TPZVec<REAL> flux(8);
    gEul.Flux(sol,flux);
    /*
	 int i;
	 cout << " flux = " << x[0] << ' ' << x[1] << ' ';
	 for(i=0; i<4; i++) cout << flux[i+4] << ' ';
	 cout << endl << "dsol";
	 for(i=0; i<4; i++) cout << dsol(1,i) << ' ';
	 cout << endl;
	 */
    TPZFMatrix<REAL> KXX(4,4),KXY(4,4),KYX(4,4),KYY(4,4),A(4,4),B(4,4);
    TPZFMatrix<REAL> jacinv(2,2);
    REAL jacdet = daxesdksi(0,0)*daxesdksi(1,1)-daxesdksi(0,1)*daxesdksi(1,0);
    jacinv(0,0) = daxesdksi(1,1)/jacdet;
    jacinv(1,1) = daxesdksi(0,0)/jacdet;
    jacinv(0,1) = -daxesdksi(0,1)/jacdet;
    jacinv(1,0) = -daxesdksi(1,0)/jacdet;
    gEul.MatrixDiff(sol,axes,jacinv,KXX,KXY,KYX,KYY);
    gEul.JacobFlux(sol,A,B);
    for(in=0; in<nshape; in++) {
		for(idf=0; idf<4; idf++) {
			ef(4*in+idf) += weight*phi(in,0)*sol[idf];
			
			//+fDeltaT*dphi(0,in)*flux[idf]
			//+fDeltaT*dphi(1,in)*flux[idf+4]);
			/*
			 for(jdf=0; jdf<4; jdf++) {
			 ef(in*4+idf,0) += weight*fDeltaT*(
			 -dphi(0,in)*KXX(idf,jdf)*dsol(0,jdf)
			 -dphi(0,in)*KXY(idf,jdf)*dsol(1,jdf)
			 -dphi(1,in)*KYX(idf,jdf)*dsol(0,jdf)
			 -dphi(1,in)*KYY(idf,jdf)*dsol(1,jdf)
			 );
			 }
			 */
			for(jn=0; jn<nshape; jn++) {
				for(jdf=0; jdf<4;jdf++) {
					
					ek(4*in+idf,4*jn+jdf) += fDeltaT*weight*(
															 dphi(0,in)*KXX(idf,jdf)*dphi(0,jn)
															 +dphi(0,in)*KXY(idf,jdf)*dphi(1,jn)
															 +dphi(1,in)*KYX(idf,jdf)*dphi(0,jn)
															 +dphi(1,in)*KYY(idf,jdf)*dphi(1,jn)
															 
															 -dphi(0,in)*A(idf,jdf)*phi(jn)
															 -dphi(1,in)*B(idf,jdf)*phi(jn)
															 ); 
				}
				ek(4*in+idf,4*jn+idf) += weight*phi(in,0)*phi(jn,0);
			}
		}
    }
}
int TPZEuler::NStateVariables()  {
	return 4;
}
int TPZEuler::Dimension() {
	return 2;
}
TPZEuler::TPZEuler(TPZEuler & copy) : TPZMaterial(copy){
	fDeltaT = copy.fDeltaT;
	fState = copy.fState;
}
TPZEuler::TPZEuler(int id, REAL deltat) : TPZMaterial(id){
	fDeltaT = deltat;
	fState = 0;
}

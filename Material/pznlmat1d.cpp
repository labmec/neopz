/**
 * \file
 * @brief Contains implementations of the TPZMat1dLin methods.
 */
#include "pzmat1dlin.h"
#include "pzmaterial.h"
//#include "pztempmat.h"
#include "pzconnect.h"
#include "pzbndcond.h"
#include "pzerror.h"
#include "pzvec.h"
//#include "pzmatrix.h"

#include <math.h>
using namespace std;

void TPZMat1dLin::Contribute(TPZMaterialData &data,
                             REAL weight,
                             TPZFMatrix<REAL> &ek, 
                             TPZFMatrix<REAL> &ef){
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
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix<REAL> &dsol=data.dsol;
	// TPZFMatrix<REAL> &dsolL=data.dsoll;
	// TPZFMatrix<REAL> &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix<REAL> &daxesdksi=data.daxesdksi;
	// TPZFMatrix<REAL> &axes=data.axes;
	
	// this method adds the contribution of the material to the stiffness
	// matrix and right hand side
	
	// check on the validity of the arguments
	
	if(phi.Cols() != 1 || dphi.Rows() != 1 || phi.Rows() != dphi.Cols()){
		PZError << "TPZMat1dLin.contr, inconsistent input data : phi.Cols() = "
	    << phi.Cols() << " dphi.Cols + " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\n";
	}
	
	if(fForcingFunction) {
		TPZManVector<REAL> xfloat(fXf.Rows());
		fForcingFunction->Execute(x,xfloat);//fXf = xfloat
		int i;
		for(i=0; i<fXf.Rows(); i++) fXf(i,0) = xfloat[i];
	}
	int r = fXk.Rows();
	int c = fXk.Cols();
	TPZFMatrix<REAL> submat(r,c);
	for(int in=0 ; in < phi.Rows() ; ++in){
		ef.AddSub(in*r, 0, (fXf*(phi(in,0)*weight)));
		for(int jn=0 ; jn<phi.Rows() ; ++jn){
			submat =  fXb*(phi(in,0)*phi(jn,0)*weight);
			submat += fXk*(dphi(0,in)*dphi(0,jn)*weight);
			submat += fXc*(phi(in,0)*dphi(0,jn)*weight);
			ek.AddSub(in*r,jn*c,submat);
		}
	}
}

void TPZMat1dLin::ContributeBC(TPZMaterialData &data,
                               REAL weight,
                               TPZFMatrix<REAL> &ek,
                               TPZFMatrix<REAL> &ef,
                               TPZBndCond &bc){
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
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix<REAL> &dsol=data.dsol;
	// TPZFMatrix<REAL> &dsolL=data.dsoll;
	// TPZFMatrix<REAL> &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix<REAL> &daxesdksi=data.daxesdksi;
	// TPZFMatrix<REAL> &axes=data.axes;
	
	//void TPZMat1dLin::ContributeBc(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*sol*/, TElementMatrix &ek, TElementMatrix &ef, TPZBndCond &bc, int nod) {
	
	// this method applies the boundary condition itype to ek and ef
	
	if(bc.Material().operator ->() != this){
		PZError << "TPZMat1dLin.apply_bc warning : this material didn't create the boundary condition!\n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "TPZMat1dLin.aplybc, unknown boundary condition type :"  <<
		bc.Type() << " boundary condition ignored\n";
	}
	int bcv1r,bcv1c,bcv2r,bcv2c;
	int r = fXk.Rows();
	int numnod = ek.Rows()/r;
	//	ekrsub = ek.mat->rowsub(0,0);
	bcv1r = bc.Val1().Rows();
	bcv1c = bc.Val1().Cols();
	bcv2r = bc.Val2().Rows();
	bcv2c = bc.Val1().Cols();
	if( bcv1r != r ||
	   bcv1c != r ||
	   bcv2r != r ||
	   bcv2c != 1 ) {
		PZError << "TPZMat1dLin.aplybc, incompatible number of degrees of " <<
		"freedom, \n"
		" val1.Rows =" << bc.Val1().Rows() << " xk.Rows = " << fXk.Rows() << "\n"
		" val2.Cols() = " << bc.Val2().Cols() << " val2.Rows() = " << bc.Val2().Rows() << "\n"
		" val1.Cols() = " << bc.Val1().Cols() << "\n";
		//		pzerror.show();
	}
	
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
			}
			break;
			
	}
	//	pzerror.show();
}

void TPZMat1dLin::Print(std::ostream & out){
	
	out << "Material type TPZMat1dLin -- number = " << Id() << "\n";
	out << "Matrix xk ->  "; fXk.Print("fXk",out);
	out << "Matrix xc ->  "; fXc.Print("fXc",out);
	out << "Matrix xb ->  "; fXb.Print("fXb",out);
	out << "Matrix xf ->  "; fXf.Print("fXf",out);
}

void TPZMat1dLin::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*u*/, TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &/*axes*/, TPZVec<REAL> &fl) {
	
	int row = NStateVariables();
	for(int i=0; i<row; i++){
		fl[i]  = 0.;
		for(int j=0; j<row; j++) {
			fl[i] += -fXk(i,j)*dudx(0,j);
		}
	}
}

void TPZMat1dLin::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,TPZFMatrix<REAL> &dudx,TPZFMatrix<REAL> &/*axes*/, TPZVec<REAL> &flux,
						 TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values) {
	
	TPZVec<REAL> udif(u);
	int nelem= udif.NElements(),i;
	for(i=0; i<nelem; i++) udif[i] -= u_exact[i];
	TPZFMatrix<REAL> dudif(dudx);
	
	int r = NStateVariables();
	TPZVec<REAL> flux_el( r );
	short idf;
	for(idf=0; idf<r; idf++) {
		dudif(0,idf) -= du_exact(0,idf);
	}
	
	values.Fill(0.);  //misael
	
	for (idf=0; idf<r; idf++) {
		values[1] += udif[idf]*udif[idf];
		for (short jdf=0; jdf<r; jdf++) {
			values[0] += dudif(0,idf)*fXk(idf,jdf)*dudif(0,jdf) + udif[idf]*fXb(idf,jdf)*udif[jdf];
			flux_el[idf] -= fXk(idf,jdf)*dudx(0,jdf);
		}
	}
	
	for (idf=0; idf<r; idf++) {
		REAL dif = flux[idf]-flux_el[idf];
		if(fabs(fXk(idf,idf)) >= 1.e-10)
		{ //Erico cout<<endl<<fXk(idf,idf)<<endl;
			values[2] += dif*dif/sqrt(fabs( fXk(idf,idf) ));
		}
	}
}


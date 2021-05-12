/**
 * @file
 * @brief Contains implementations of the TPZNLMat1d methods. (non linear one dimensional equation)
 */

#include "pznlmat1d.h"
#include "TPZMaterial.h"
#include "pzconnect.h"
#include "pzbndcond.h"
#include "pzerror.h"
#include "pzvec.h"
#include "Hash/TPZHash.h"
#include <math.h>
using namespace std;

TPZNLMat1d::TPZNLMat1d(int id) : TPZMaterial(id)
{}


TPZNLMat1d::~TPZNLMat1d()
{}

int TPZNLMat1d::ClassId() const{
return Hash("TPZNLMat1d") ^ TPZMaterial::ClassId() << 1;
}
void TPZNLMat1d::Contribute(TPZMaterialData &data, REAL weight,
                             TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
//	TPZManVector<REAL,3> &x = data.x;
	
	// this method adds the contribution of the material to the stiffness
	// matrix and right hand side
	
	// check on the validity of the arguments
	
	if(phi.Cols() != 1 || dphi.Rows() != 1 || phi.Rows() != dphi.Cols()){
		PZError << "TPZMat1dLin.contr, inconsistent input data : phi.Cols() = "
	    << phi.Cols() << " dphi.Cols + " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\n";
	}
	
	// IT IS NOT IMPLEMENTED 
	
/*	if(fForcingFunction) {
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
 */
}

void TPZNLMat1d::ContributeBC(TPZMaterialData &data, REAL weight,
                               TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                               TPZBndCond &bc) {
	
//	TPZFMatrix<REAL> &phi = data.phi;
	
	// this method applies the boundary condition itype to ek and ef
	
	if(bc.Material() != this){
		PZError << "TPZMat1dLin.apply_bc warning : this material didn't create the boundary condition!\n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "TPZMat1dLin.aplybc, unknown boundary condition type :"  <<
		bc.Type() << " boundary condition ignored\n";
	}
//	int bcv1r,bcv1c,bcv2r,bcv2c;

	// IT IS NOT IMPLEMENTED YET
	
/*	int r = fXk.Rows();
	int numnod = ek.Rows()/r;
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
 */
}

void TPZNLMat1d::Print(std::ostream & out){
	
	out << "Material type TPZNLMat1D -- number = " << Id() << "\n";
//	out << "Matrix xk ->  "; fXk.Print("fXk",out);
//	out << "Matrix xc ->  "; fXc.Print("fXc",out);
//	out << "Matrix xb ->  "; fXb.Print("fXb",out);
//	out << "Matrix xf ->  "; fXf.Print("fXf",out);
}


#include "MiViga1D.h"

#include "pzconnect.h"
#include "pzbndcond.h"

TPZMiViga1D::TPZMiViga1D(int nId,REAL Elasticity,REAL Poisson,REAL Inercia,REAL Ksi,REAL Area) : TPZMaterial(nId) {
	fE = Elasticity;

	if(Poisson < 0) fNi = 0.;
	else if(Poisson > 0.5) fNi = 0.5;
	else fNi = Poisson;

	fI = Inercia;
	fKsi = Ksi;
	fArea = Area;

	fG = fE/(2*(1.+fNi));
}

void TPZMiViga1D::Contribute(TPZMaterialData &data,
                             REAL weight,
                             TPZFMatrix<STATE> &ek, 
                             TPZFMatrix<STATE> &ef) {

	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
//	TPZManVector<REAL,3> &x = data.x;

	// check on the validity of the arguments
	if(phi.Cols() != 1 || dphi.Rows() != 1 || phi.Rows() != dphi.Cols()) {
		PZError << "TPZMiViga1D.contr, inconsistent input data : phi.Cols() = "
	    << phi.Cols() << " dphi.Cols + " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\n";
	}
	int r = phi.Rows();
	int c = phi.Cols();
	
	TPZFMatrix<STATE> submat(r,c);
	for(int in=0 ; in < phi.Rows() ; ++in) {
		// contribucion en el vector de carga
		ef(2*in,0) += weight*fCarga*phi(in,0);
		// contribucion en la matriz rigidez
		for(int jn=0 ; jn<phi.Rows() ; ++jn) {
			ek(2*in,2*jn) += weight*dphi(0,in)*dphi(0,jn)*fG*fArea*fKsi;
			ek(2*in,2*jn+1) += weight*dphi(0,in)*phi(jn,0)*(-fG*fArea*fKsi);
			ek(2*in+1,2*jn) += weight*phi(in,0)*dphi(0,jn)*(-fG*fArea*fKsi);
			ek(2*in+1,2*jn+1) += weight*(fE*fI*dphi(0,in)*dphi(0,jn)+fG*fArea*fKsi*phi(in,0)*phi(jn,0));
		}
	}
}

void TPZMiViga1D::ContributeBC(TPZMaterialData &data,
                               REAL weight,
                               TPZFMatrix<STATE> &ek,
                               TPZFMatrix<STATE> &ef,
                               TPZBndCond &bc) {
	TPZFMatrix<REAL> &phi = data.phi;
	
	// this method applies the boundary condition itype to ek and ef
	
	if(bc.Material() == this){
		PZError << "TPZMiViga1D.apply_bc warning : this material didn't create the boundary condition!\n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 1){
		PZError << "TPZMiViga1D.aplybc, unknown boundary condition type :"  <<
		bc.Type() << " boundary condition ignored\n";
	}

	int numnod = ek.Rows()/2;
	
	int idf, in, jn;
	switch(bc.Type()) {
		case 0:  // Dirichlet
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<2;idf++) {
					ef(2*in+idf,0) += gBigNumber*phi(in,0)*bc.Val2()(idf,0)*weight;
				}
				for(jn=0; jn<numnod; ++jn) {
					for(idf = 0;idf<2;idf++) {
						ek(2*in+idf,2*jn+idf) += gBigNumber*phi(in,0)*phi(jn,0)*weight;
					}
				}
			}
			break;
		case 1:   // Neumann
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<2;idf++) {
					ef(2*in+idf,0) += phi(in,0)*bc.Val2()(idf,0)*weight;
				}
			}
			break;			
	}
}

int TPZMiViga1D::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"desplazamiento")) return 1;
	if(!strcmp(name.c_str(),"giro")) return 2;
	return TPZMaterial::VariableIndex(name);
}

int TPZMiViga1D::NSolutionVariables(int var) { 
	switch(var) {
	case 1:
		return 2;
	case 2:
		return 2;
	}
	return -1;
}

void TPZMiViga1D::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) {
	
	switch(var) {
		case 1:   // desplazamiento
			Solout[0] = 0.;
			Solout[1] = Sol[0];
			break;
		case 2:   // giro
			Solout[0] = 0.;
			Solout[1] = Sol[1];
			break;
		default:
			std::cout << "TPZElasticityMaterial::Solution Error\n";
			TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
			break;
	}
}

/**
 * \file
 * @brief Contains implementations of the TPZMatHybrid methods.
 */
#include "TPZFVHybrid.h"

#include "pzbndcond.h"

#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

using namespace std;

TPZMatHybrid::TPZMatHybrid(int nummat) : TPZMaterial(nummat), fXf(1,1,0.) {
	
	fNumMat = nummat;//material id
}

TPZMatHybrid::~TPZMatHybrid() {
}

int TPZMatHybrid::NStateVariables() const {
	return 1;
}

void TPZMatHybrid::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "Material id   = " << fNumMat   << endl;
}

void TPZMatHybrid::Contribute(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef) {
	
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
	TPZFMatrix<REAL> &axes=data.axes;
	
	
	int phr = phi.Rows();
	int dim = phi.Cols();
	int k;
	if(fForcingFunction) {            // phi(in, 0) :  fun�o de forma associada ao n�in
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(x,res);       // dphi(i,j) :  derivada c/r a xi da fun�o de forma j
		fXf(0,0) = res[0];
	}
	
	TPZVec<REAL> normal(3);
	for(k=0;k<3;k++) normal[k] = axes(0,k);//primeira linha de axes
	
	for( int in = 0; in < phr; in++ ) {
		ef(in, 0) += weight * fXf(0,0) * phi(in, 0) ;
		for( int jn = 0; jn < phr; jn++ ) {
			REAL sum = 0.0;
			for(k=0;k<dim;k++) sum +=dphi(k,in) * dphi(k,jn);
			ek(in,jn) += weight * sum;
		}
	}
}

void TPZMatHybrid::ContributeBC(TPZMaterialData &data,
                                REAL weight,
                                TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef,
                                TPZBndCond &bc) {
	
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
	
	if(bc.Material() != this){
		PZError << "TPZMat1dLin.apply_bc warning : this material didn't create the boundary condition!\n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "TPZMat1dLin.aplybc, unknown boundary condition type :"  <<
		bc.Type() << " boundary condition ignored\n";
	}
	
	
	int numdof = NStateVariables();
	int numnod = ek.Rows()/numdof;
	int r = numdof;
	
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
			}//fim switch
	}
}

/** returns the variable index associated with the name*/
int TPZMatHybrid::VariableIndex(const std::string &name){
	if(!strcmp("Pressao",name.c_str()))       return 0;
	if(!strcmp("Solution",name.c_str()))      return 1;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMatHybrid::NSolutionVariables(int var){
	
	if(var == 0 || var == 1 || var == 2) return 1;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMatHybrid::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout){
	
	if(var == 1) {
		cout << "TPZMatHybrid::Solution implementar Pressao\n";
		return;
	}
	if(var == 1) {
		Solout[0] = Sol[0];
		return;
	}
	
	TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

//ofstream arq1("Simetria.dat");
void TPZMatHybrid::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
						  TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
	
	TPZVec<STATE> sol(1),dsol(2);
	Solution(u,dudx,axes,1,sol);
	Solution(u,dudx,axes,2,dsol);
	for (int i = 0; i < 3; i++)
		values[i] = 0.0;
	//values[0] : norma energia do erro
	values[0] = (dsol[0] - du_exact(0,0))*(dsol[0] - du_exact(0,0));
	values[0] += (dsol[1] - du_exact(1,0))*(dsol[1] - du_exact(1,0));
	//values[1] -> norma energia da solucao Uhp (aproximada)
	values[1]  = (dsol[0])*(dsol[0]);
	values[1] += (dsol[1])*(dsol[1]);
	//values[2] -> morma enaergia da solucao U (exata)
	values[2]  = (du_exact(0,0))*(du_exact(0,0));
	values[2] += (du_exact(1,0))*(du_exact(1,0));
}


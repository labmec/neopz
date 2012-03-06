/**
 * @file
 * @brief Contains the implementations of the TPZConsLawTest methods.
 */
#include "TPZConsLawTest.h" 
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>
#include <cmath>

using namespace std;


TPZConsLawTest::TPZConsLawTest(int nummat, TPZVec<REAL> B,int artdiff,REAL delta_t,int dim,REAL delta,int test) : TPZConservationLaw(nummat,delta_t,dim), fXf(1,1,0.), fB(dim)  {
	
	if(artdiff<0 || artdiff>3){
		PZError << "TPZConsLawTest::TPZConsLawTest artificial diffusion parameter, default 1\n";
		artdiff = 1;
	}
	if(B.NElements() != dim){
		PZError << "TPZConsLawTest::TPZConsLawTest error in dimension of B, default dimension " << dim << endl;
		B.Resize(dim);
		int i;
		for(i=0;i<dim;i++) B[i] = 0.;
		B[dim-1] = 1.0;
	}
	int i;
	for(i=0;i<dim;i++) fB[i] = B[i];
	fArtificialDiffusion = artdiff;//SUPG, SL or Bornhauss
	fDelta = delta;
	fTest = test;
}

TPZConsLawTest::~TPZConsLawTest() {
}

int TPZConsLawTest::NStateVariables() {
	return 1;
}

void TPZConsLawTest::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";           
	TPZMaterial::Print(out);
}

void TPZConsLawTest::Contribute(TPZMaterialData &data,
                                REAL weight,
                                TPZFMatrix &ek,
                                TPZFMatrix &ef) {
	
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
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZVec<REAL> &sol=data.sol[0];
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	int phr = phi.Rows();// phi(in, 0) = phi_in  ,  dphi(i,jn) = dphi_jn/dxi
	
	if(fForcingFunction) {
		TPZManVector<REAL> res(1);
		fForcingFunction->Execute(x,res);//fXf(0,0) = res[0];//if(!sol[0]) 
		sol[0] = res[0];//solu�o inicial, na itera�o 0 sol = 0 e res != 0 
	}                  //nas itera�es > 0 sol != 0 e res = 0
	int dim = dphi.Rows();
	if(Dimension() != dim)
		PZError << "TPZConsLawTest::Contribute dimension error, dimension = " << dim;
	
	REAL sum1 = 0.,sum2=0.;
	int i;
	REAL time = TimeStep();
	for( int in = 0; in < phr; in++ ) {
		//primeira parcela
		sum1 = phi(in, 0) * sol[0];
		//segunda parcela
		sum2 = 0.;
		for(i=0;i<dim;i++) sum2 += B(i,x) * dphi(i,in);//F = (B0*u.B1*u,B2*u)
		sum2 *= time * sol[0];
		//terceira parcela: contribui� do elemento interface: feita na classe TPZInterfaceElement
		//contribui� final
		ef(in, 0) += weight *(sum1+sum2);//* fXf(0,0);
		
		for( int jn = 0; jn < phr; jn++ ) {
			//primeira parcela
			ek(in,jn) += weight * phi(in,0) * phi(jn,0);
			for(int ki=0; ki<dim; ki++) {
				for(int kj=0; kj<dim; kj++) {
					//segunda parcela
					ek(in,jn) += weight * time * Delta() * T(ki,x) * B(kj,x) * ( dphi(kj,in) * dphi(ki,jn) );
				}
			}
		}
	}
}

REAL TPZConsLawTest::Delta(){
	
	return fDelta;
}

REAL TPZConsLawTest::DeltaOtimo(){
	
	//REAL cfl = CFL(1);
	//REAL delta = ( (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10. );
	return 0.0;//delta;
}

REAL TPZConsLawTest::CFL(int degree){
	
	return (1.0/(2.0*degree+1.0));
} 

REAL TPZConsLawTest::B(int i,TPZVec<REAL> &x){
	
	if(fTest==0) return fB[i];
	
	if(fTest==1){
		if(i==0) return -x[1];
		if(i==1) return  x[0];
		if(i==2) return  0.0;
	}
	if(fTest==2){
		if(i==0) return -x[1];
		if(i==1) return  x[0];
		if(i==2) return  0.;
	}
	return 0.;
}

REAL TPZConsLawTest::T(int jn,TPZVec<REAL> &x){
	
	//SUPG
	if(fArtificialDiffusion == 1){
		int i;
		REAL norm = 0;
		int dim = Dimension();
		for(i=0;i<dim;i++) norm += B(i,x)*B(i,x);
		norm = sqrt(norm);
		if(jn==0) return B(0,x)/norm;
		if(jn==1) return B(1,x)/norm;
		if(jn==2) return B(2,x)/norm;
	}
	//LS
	if(fArtificialDiffusion == 2){
		cout << "TPZConsLawTest::T artificial diffusion LS not implemented\n";
		return 0;
	}
	//Bornhaus
	if(fArtificialDiffusion == 3){
		cout << "TPZConsLawTest::T artificial diffusion Bornhaus not implemented\n";
		return  0;
	}
	cout << "TPZConsLawTest::T case not implemented, case = " << jn << endl;
	return 0.0;
}

void TPZConsLawTest::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                         REAL weight,
                                         TPZFMatrix &ek,
                                         TPZFMatrix &ef){
	
	// TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	// TPZFMatrix &phi = data.phi;
	TPZFMatrix &phiL = dataleft.phi;
	TPZFMatrix &phiR = dataright.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZVec<REAL> &solL=dataleft.sol[0];
	TPZVec<REAL> &solR=dataright.sol[0];
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	
	int phrl = phiL.Rows();
	int phrr = phiR.Rows();
	
	if(fForcingFunction) {      // phi(in, 0) = phi_in
		TPZManVector<REAL> res(1);// dphi(i,j) = dphi_j/dxi
		fForcingFunction->Execute(x,res);
		if(phrl) solL[0] = res[0];//solu�o inicial interior (t=0), na itera�o 0 sol = 0 e res != 0
		if(phrr) solR[0] = res[0];//nas itera�es > 0 sol != 0 e res = 0
	}
	
	REAL Bn=0.0;
	//a normal est�imersa no mesmo espao da malha: R, R, R
	int dim = Dimension(),i;
	for(i=0;i<dim;i++) Bn += B(i,x)*normal[i];
	//contribui� do elemento interface
	REAL time = TimeStep();
	int efc = 0;
	if(Bn > 0){
		//a normal �interface sempre aponta do seu elemento esquerdo para o seu direito
		for( int in = 0; in < phrl; in++ ) {
			ef(efc  , 0) += -time * weight * solL[0] * phiL(in, 0) * Bn;//0 2 4 ..
			efc++;
		}
		for( int in = 0; in < phrr; in++ ) {
			ef(efc, 0) += -time * weight * solL[0] * phiR(in, 0) * -Bn;//1 3 5 ..
			efc++;
		}
	} else {
		for( int in = 0; in < phrl; in++ ) {
			ef(efc  , 0) += -time * weight * solR[0] * phiL(in, 0) * Bn;//0 2 4 ..
			efc++;
		}
		for( int in = 0; in < phrr; in++ ) {
			ef(efc, 0) += -time * weight * solR[0] * phiR(in, 0) * -Bn;//1 3 5 ..
			efc++;
		}
	}
}

void TPZConsLawTest::ContributeBC(TPZMaterialData &data,
                                  REAL weight,
                                  TPZFMatrix &ek,
                                  TPZFMatrix &ef,
                                  TPZBndCond &bc) {
	
	// TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	int phr = phi.Rows();
	short in,jn;
	REAL v2[1];
	v2[0] = bc.Val2()(0,0);
	
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				ef(in,0) += gBigNumber * v2[0] * phi(in,0) * weight;
				for (jn = 0 ; jn < phr; jn++) {
					ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
				}
			}
			break;
		case 1 :			// Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(in,0) += v2[0] * phi(in,0) * weight;
			}
			break;
		case 2 :		// condi�o mista
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(in, 0) += v2[0] * phi(in, 0) * weight;
				for (jn = 0 ; jn < phi.Rows(); jn++) {
					ek(in,jn) += bc.Val1()(0,0) * phi(in,0) *
					phi(jn,0) * weight;     // peso de contorno => integral de contorno
				}
			}
	}
}

/** returns the variable index associated with the name*/
int TPZConsLawTest::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str()))        return  1;
	if(!strcmp("Derivate",name.c_str()))        return  2;
	cout << "TPZConsLawTest::VariableIndex Error\n";
	return -1;
}

int TPZConsLawTest::NSolutionVariables(int var){
	
	if(var == 1) return 1;
	if(var == 2) return Dimension();
	cout << "TPZConsLawTest::NSolutionVariables Error\n";
	return 0;
}

void TPZConsLawTest::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout){
	
	if(var == 0 || var == 1) Solout[0] = Sol[0];//function
	if(var == 2) {
		Solout[0] = DSol(0,0);//derivate
		Solout[1] = DSol(1,0);//derivate
		Solout[2] = DSol(2,0);//derivate
	}
}

void TPZConsLawTest::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
	//Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

void TPZConsLawTest::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
							TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
							TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {
	
	TPZVec<REAL> sol(1),dsol(3);
	Solution(u,dudx,axes,1,sol);
	Solution(u,dudx,axes,2,dsol);
	REAL dx = dsol[0]*axes(0,0)+dsol[1]*axes(1,0)+dsol[2]*axes(2,0);
	REAL dy = dsol[0]*axes(0,1)+dsol[1]*axes(1,1)+dsol[2]*axes(2,1);
	REAL dz = dsol[0]*axes(0,2)+dsol[1]*axes(1,2)+dsol[2]*axes(2,2);
	//values[1] : eror em norma L2
	values[1]  = pow(sol[0] - u_exact[0],(REAL)2.0);
	//values[2] : erro em semi norma H1
	values[2]  = pow(dx - du_exact(0,0),(REAL)2.0);
	values[2] += pow(dy - du_exact(1,0),(REAL)2.0);
	values[2] += pow(dz - du_exact(2,0),(REAL)2.0);
	//values[0] : erro em norma H1 <=> norma Energia
	values[0]  = values[1]+values[2];
}

void TPZConsLawTest::ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft){
	
	if(!bcleft){
		PZError << "TPZConsLawTest::ComputeSolLeft null bundary condition return\n";
		return;
	}
	//int nstate = NStateVariables();
	TPZFMatrix jacinv(0,0),axes(0,0);
	switch (bcleft->Type()){
		case 0://Dirichlet
		case 1://Neumann
		case 2://Mista
			PZError << "TPZConsLawTest::ComputeSolLeft boundary condition error\n";
			break;
		case 3://Dirichlet: nada a fazer a CC �a correta
			break;
		case 4://recuperar valor da solu� MEF direita: saida livre
			soll[0] = solr[0];
			break;
		default:
			PZError << "TPZConsLawTest::ContributeInterface Boundary Condition Type Not Exists\n";
	}
}


void TPZConsLawTest::ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright){
	
	if(!bcright){
		PZError << "TPZConsLawTest::ComputeSolLeft null bundary condition return\n";
		return;
	}
	//int nstate = NStateVariables();
	TPZFMatrix jacinv(0,0),axes(0,0);
	switch (bcright->Type()){
		case 0://Dirichlet
		case 1://Neumann
		case 2://Mista
			PZError << "TPZConsLawTest::ComputeSolLeft boundary condition error\n";
			break;
		case 3://Dirichlet: nada a fazer a CC �a correta
			break;
		case 4://recuperar valor da solu� MEF esquerda: saida livre
			solr[0] = soll[0];
			break;
		default:
			PZError << "TPZConsLawTest::ContributeInterface Boundary Condition Type Not Exists\n";
	}
}

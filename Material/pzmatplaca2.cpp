/**
 * \file
 * @brief Contains implementations of the TPZMatPlaca2 methods.
 */

#include "pzmatplaca2.h" 
#include "pzmaterial.h"

#include "pzbndcond.h"
#include <math.h>
#include <fstream>
using namespace std;
#include "pzvec.h"
#include "pzerror.h"


TPZMatPlaca2::TPZMatPlaca2(int num, REAL h, REAL f, REAL E1 , REAL E2 ,
						   REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
						   REAL G23 , TPZFMatrix &naxes, TPZVec<REAL> &xf) :
TPZMaterial(num),
fIdfMax(6),fE1(E1), fE2(E2), fG12(G12), fG13(G13), fG23(G23),
fh(h),ff(f),fmi(1./(-1.+ni1*ni2)),fni1(ni1),fni2(ni2),
fnaxes(naxes),
fRmat(6,6,0.),fRmatT(6,6,0.),
fKxxR(6,6,0.),fKyyR(6,6,0.),
fKxyR(6,6,0.),fKyxR(6,6,0.),
fBx0R(6,6,0.),fB0xR(6,6,0.),
fBy0R(6,6,0.),fB0yR(6,6,0.),
fB00R(6,6,0.),
fKxx(6,6,0.),fKyy(6,6,0.),
fKxy(6,6,0.),fKyx(6,6,0.),
fBx0(6,6,0.),fB0x(6,6,0.),
fBy0(6,6,0.),fB0y(6,6,0.),
fB00(6,6,0.), fXF(xf)
{
	REAL Small , k, mi;
	Small = E1*REAL(1.E-5);
	k = 5./6.; // coeficiente de cisalhamento
	mi = 1.0/(-1.0 + ni1 * ni2);
	
	fRmat(0,0) = fnaxes(0,0); fRmat(0,1) = fnaxes(0,1); fRmat(0,2) = fnaxes(0,2);
	fRmat(1,0) = fnaxes(1,0); fRmat(1,1) = fnaxes(1,1); fRmat(1,2) = fnaxes(1,2);
	fRmat(2,0) = fnaxes(2,0); fRmat(2,1) = fnaxes(2,1); fRmat(2,2) = fnaxes(2,2);
	
	fRmat(3,3) = fnaxes(0,0); fRmat(3,4) = fnaxes(0,1); fRmat(3,5) = fnaxes(0,2);
	fRmat(4,3) = fnaxes(1,0); fRmat(4,4) = fnaxes(1,1); fRmat(4,5) = fnaxes(1,2);
	fRmat(5,3) = fnaxes(2,0); fRmat(5,4) = fnaxes(2,1); fRmat(5,5) = fnaxes(2,2);
	
	fRmat.Transpose(&fRmatT);
	
	
	
	fKxx(0,0) = -E1*h*mi;
	fKxx(0,4) = -E1*f*h*mi;
	fKxx(1,1) =  G12*h/2. + h*Small/4.;
	fKxx(1,3) = -f*G12*h/2.;
	fKxx(2,2) =  G13*h*k/2.;
	fKxx(3,1) = -f*G12*h/2.;
	fKxx(3,3) =  f*f*G12*h/2. + G12*h*h*h/24.;
	fKxx(4,0) = -E1*f*h*mi;
	fKxx(4,4) = -E1*f*f*h*mi - E1*h*h*h*mi/12.;
	
	fKxy(0,1) = -E1*h*mi*ni2;
	fKxy(0,3) =  E1*f*h*mi*ni2;
	fKxy(1,0) =  G12*h/2. - h*Small/4.;
	fKxy(1,4) =  f*G12*h/2.;
	fKxy(3,0) = -f*G12*h/2.;
	fKxy(3,4) = -f*f*G12*h/2. - G12*h*h*h/24.;
	fKxy(4,1) = -E1*f*h*mi*ni2;
	fKxy(4,3) =  E1*f*f*h*mi*ni2 + E1*h*h*h*mi*ni2/12.;
	
	fKxy.Transpose(&fKyx);
	
	fKyy(0,0) =  G12*h/2. + h*Small/4.;
	fKyy(0,4) =  f*G12*h/2.;
	fKyy(1,1) = -E2*h*mi;
	fKyy(1,3) =  E2*f*h*mi;
	fKyy(2,2) =  G23*h*k/2.;
	fKyy(3,1) =  E2*f*h*mi;
	fKyy(3,3) = -E2*f*f*h*mi - E2*h*h*h*mi/12.;
	fKyy(4,0) =  f*G12*h/2.;
	fKyy(4,4) =  f*f*G12*h/2. + G12*h*h*h/24.;
	
	fB0x(4,2) =  G13*h*k/2.;
	fB0x(5,1) = -h*Small/2.;
	
	fB0x.Transpose(&fBx0);
	
	fB0y(3,2) = -G23*h*k/2.;
	fB0y(5,0) =  h*Small/2.;
	
	fB0y.Transpose(&fBy0);
	
	fB00(3,3) =  G23*h*k/2.;
	fB00(4,4) =  G13*h*k/2.;
	fB00(5,5) =  h*Small;
	
	fKxxR = fRmatT * (fKxx * fRmat);
	fKyxR = fRmatT * (fKyx * fRmat);
	fKxyR = fRmatT * (fKxy * fRmat);
	fKyyR = fRmatT * (fKyy * fRmat);
	fB0xR = fRmatT * (fB0x * fRmat);
	fB0yR = fRmatT * (fB0y * fRmat);
	fBx0R = fRmatT * (fBx0 * fRmat);
	fBy0R = fRmatT * (fBy0 * fRmat);
	fB00R = fRmatT * (fB00 * fRmat);
	
}

void TPZMatPlaca2::SetNAxes(TPZFMatrix &n) {
	
	fnaxes = n;
	int numl = fIdfMax/3;
	int i;
	for(i=0;i<numl;i++) {
		fRmat(3*i+0,3*i+0) = fnaxes(0,0); 
		fRmat(3*i+0,3*i+1) = fnaxes(0,1); 
		fRmat(3*i+0,3*i+2) = fnaxes(0,2);
		fRmat(3*i+1,3*i+0) = fnaxes(1,0); 
		fRmat(3*i+1,3*i+1) = fnaxes(1,1); 
		fRmat(3*i+1,3*i+2) = fnaxes(1,2);
		fRmat(3*i+2,3*i+0) = fnaxes(2,0); 
		fRmat(3*i+2,3*i+1) = fnaxes(2,1); 
		fRmat(3*i+2,3*i+2) = fnaxes(2,2);
	}
	
	fRmat.Transpose(&fRmatT);
	
	
	
	
	
	
	
	TPZFMatrix tmp;
	fKxx.Multiply(fRmat,tmp);
	fRmatT.Multiply(tmp,fKxxR);
	fKyxR = fRmatT * (fKyx * fRmat);
	fKxyR = fRmatT * (fKxy * fRmat);
	fKyyR = fRmatT * (fKyy * fRmat);
	fB0xR = fRmatT * (fB0x * fRmat);
	fB0yR = fRmatT * (fB0y * fRmat);
	fBx0R = fRmatT * (fBx0 * fRmat);
	fBy0R = fRmatT * (fBy0 * fRmat);
	fB00R = fRmatT * (fB00 * fRmat);
	
	
}

/** @brief Output file to write data of the test over shell (placa)*/
ofstream placatest("placatest.dat");
void TPZMatPlaca2::Contribute(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef) {
	// this method adds the contribution of the material to the stiffness
	// matrix and right hand side
	
	// check on the validity of the arguments
	//rows x cols
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
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	TPZFMatrix &axes=data.axes;
	
	if(phi.Cols() != 1 || dphi.Rows() != 2 || phi.Rows() != dphi.Cols()){
		PZError << "TPZMatPlaca2.contr, inconsistent input data : phi.Cols() = "
	    << phi.Cols() << " dphi.Cols + " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\n";
	}
	if(fForcingFunction) {
		fForcingFunction->Execute(x,fXF);//fXf = xfloat
	}
	
	
	REAL Dax1n1, Dax1n2, Dax2n1, Dax2n2,intern33;
	
	Dax1n1 = axes(0,0)* fnaxes(0,0) + axes(0,1)* fnaxes(0,1) + axes(0,2)* fnaxes(0,2);
	Dax1n2 = axes(0,0)* fnaxes(1,0) + axes(0,1)* fnaxes(1,1) + axes(0,2)* fnaxes(1,2);
	Dax2n1 = axes(1,0)* fnaxes(0,0) + axes(1,1)* fnaxes(0,1) + axes(1,2)* fnaxes(0,2);
	Dax2n2 = axes(1,0)* fnaxes(1,0) + axes(1,1)* fnaxes(1,1) + axes(1,2)* fnaxes(1,2);
	intern33 = axes(2,0)* fnaxes(2,0) + axes(2,1)* fnaxes(2,1) + axes(2,2)* fnaxes(2,2);
	if(fabs(fabs(intern33)-1.) > 1.e-6) {
		PZError << "third axes of the material not parallel to element axis " << intern33 << endl;
		PZError << axes(2,0) << ' ' << axes(2,1) << ' ' << axes(2,2) << endl;
		PZError << fnaxes(2,0) << ' ' << fnaxes(2,1) << ' ' << fnaxes(2,2) << endl;
	}
	
	TPZFMatrix Kn1n1(fIdfMax,fIdfMax,0.),Kn2n2(fIdfMax,fIdfMax,0.),
    Kn1n2(fIdfMax,fIdfMax,0.),Kn2n1(fIdfMax,fIdfMax,0.),
    Bn10(fIdfMax,fIdfMax,0.) ,B0n1(fIdfMax,fIdfMax,0.),
    Bn20(fIdfMax,fIdfMax,0.),B0n2(fIdfMax,fIdfMax,0.),
    B000(fIdfMax,fIdfMax,0.);
	
	
	Kn1n1 = fKxxR * Dax1n1 * Dax1n1 + fKyyR * Dax1n2 * Dax1n2 +
    fKxyR * Dax1n1 * Dax1n2 + fKyxR * Dax1n1 * Dax1n2 ;
	
	Kn2n2 = fKxxR * Dax2n1 * Dax2n1 + fKyyR * Dax2n2 * Dax2n2 +
    fKxyR * Dax2n2 * Dax2n1 + fKyxR * Dax2n2 * Dax2n1;
	
	Kn1n2 = fKxxR * Dax1n1 * Dax2n1 + fKyyR * Dax1n2 * Dax2n2 +
    fKxyR * Dax1n1 * Dax2n2 + fKyxR * Dax1n2 * Dax2n1;
	
	Kn2n1 = fKxxR * Dax2n1 * Dax1n1 + fKyyR * Dax2n2 * Dax1n2 +
    fKxyR * Dax2n1 * Dax1n2 + fKyxR * Dax2n2 * Dax1n1;
	
	
	B0n1  = fB0xR * Dax1n1 + fB0yR * Dax1n2 ;
	
	Bn10  = fBx0R * Dax1n1 + fBy0R * Dax1n2 ;
	
	B0n2  = fB0xR * Dax2n1 + fB0yR * Dax2n2;
	
	Bn20  = fBx0R * Dax2n1 + fBy0R * Dax2n2;
	
	B000  = fB00R ;
	
	/*
	 Kn1n1.Print("Kn1n1 = ",placatest);
	 Kn1n2.Print("Kn1n2 = ",placatest);
	 Kn2n2.Print("Kn2n2 = ",placatest);
	 Kn2n1.Print("Kn2n1 = ",placatest);
	 B0n1.Print("B0n1 = ",placatest);
	 Bn10.Print("Bn10 = ",placatest);
	 B0n2.Print("B0n2 = ",placatest);
	 Bn20.Print("B0n1 = ",placatest);
	 B000.Print("B00 = ",placatest);
	 placatest.flush();
	 */
	int nshape = phi.Rows();
	TPZFMatrix KIJ(fIdfMax,fIdfMax);
	
	
	
	int idf,jdf,i,j;
	REAL contrib[3];
	for(i=0;i<3;i++) {
		contrib[i]=0.;
		for(j=0;j<3;j++) {
			contrib[i] += fXF[j]*fnaxes(j,i);
		}
	}
	for(i=0; i<nshape; i++) {
		for(idf=0; idf<3; idf++) {
			ef(fIdfMax*i+idf) += weight*contrib[idf]*phi(i,0);
		}
		for(idf=3; idf<fIdfMax; idf++) ef(fIdfMax*i+idf) += weight*fXF[idf]*phi(i,0);
		for(j=0; j<nshape; j++) {
			KIJ = weight*(dphi(0,i)*dphi(0,j)*Kn1n1+
						  dphi(0,i)*dphi(1,j)*Kn1n2+
						  dphi(1,i)*dphi(0,j)*Kn2n1+
						  dphi(1,i)*dphi(1,j)*Kn2n2+
						  dphi(0,i)*phi(j)   *Bn10 +
						  dphi(1,i)*phi(j)   *Bn20 +
						  phi(i)   *dphi(0,j)*B0n1 +
						  phi(i)   *dphi(1,j)*B0n2 +
						  phi(i)   *phi(j)   *B000 );
			for(idf=0; idf<fIdfMax; idf++) for(jdf=0; jdf<fIdfMax; jdf++)
				ek(i*fIdfMax+idf,j*fIdfMax+jdf) += KIJ(idf,jdf);
		}
	}
}

void TPZMatPlaca2::ContributeBC(TPZMaterialData &data,
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
	
	if(bc.Material().operator ->() != this){
		PZError << "TPZMatPlaca2.ContributeBC warning : this material didn't create the boundary condition!\n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "TPZMatPlaca2.aplybc, unknown boundary condition type :"  <<
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

int TPZMatPlaca2::NFluxes() {return 1;}

void TPZMatPlaca2::Flux(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*u*/,TPZFMatrix &/*dudx*/,TPZFMatrix &/*axes*/,TPZVec<REAL> &/*fl*/) {
	PZError << "TPZMatPlaca2::Flux is called\n";
}

void TPZMatPlaca2::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*u*/,TPZFMatrix &/*dudx*/,TPZFMatrix &/*axes*/,TPZVec<REAL> &/*flux*/,
						  TPZVec<REAL> &/*u_exact*/,TPZFMatrix &/*du_exact*/,TPZVec<REAL> &/*values*/) {
	PZError << "TPZMatPlaca2::Errors is called\n";
}

void TPZMatPlaca2::Print(std::ostream & out) {
	//out << "Material type TPZMatPlaca2 -- number = " << Id() << "\n";
	//out << "Matrix xk ->  "; fXk.Print("fXk",out);
	//out << "Matrix xC ->  "; fCk.Print("fCf",out);
	//out << "Matrix xf ->  "; fXf.Print("fXf",out);
}

/*TPZBndCond *TPZMatPlaca2::CreateBC(int num,int typ,TPZFMatrix &val1,TPZFMatrix &val2) {
 PZError << "TPZMatPlaca2::CreateBC is called\n";
 return 0;
 }        */

/**returns the variable index associated with the name*/
int TPZMatPlaca2::VariableIndex(const std::string &name){
	if(!strcmp(name.c_str(),"Deslocamentos nodais")) return 0;
	if(!strcmp(name.c_str(),"Deslocx")) return 2;// Desloc. eixo x global
	if(!strcmp(name.c_str(),"Deslocy")) return 3;// Desloc. eixo y global
	if(!strcmp(name.c_str(),"Deslocz")) return 4;// Desloc. eixo z global
	if(!strcmp(name.c_str(),"Mn1 Mn2 e Mn1n2"))     return 5;// Momentos nas direcoes dos eixos n1 e n2
	if(!strcmp(name.c_str(),"Ma1 Ma2 e Ma1a2"))     return 50;// Momentos nas direcoes dos eixos a1 e a2
	//if(!strcmp(name,"Mn2"))     return 6;// Mom. fletor eixo n2 da fibra
	//if(!strcmp(name,"Mn1n2"))   return 7;// Mom. volvente eixos n1 e n2 da fibra
	if(!strcmp(name.c_str(),"Vn1"))     return 8;// forca cortante Vn1 (positiva se antihorario)
	if(!strcmp(name.c_str(),"Vn2"))     return 9;// forca cortante Vn2 (positiva se antihorario)
	if(!strcmp(name.c_str(),"Sign1"))   return 10;// tens� normal na dire�o n1
	if(!strcmp(name.c_str(),"Sign2"))   return 11;// tens� normal na dire�o n2
	if(!strcmp(name.c_str(),"Taun1n2")) return 12;// tens� cisalhamento eixos n1 e n2
	if(!strcmp(name.c_str(),"Na1, Na2 e Na1a2")) return 54;//Tensoes normais nas direcoes dos eixos a1,a2
	if(!strcmp(name.c_str(),"Taun1n3")) return 13;// tens� cisalhamento eixos n1 e n3
	if(!strcmp(name.c_str(),"Taun2n3")) return 14;// tens� cisalhamento eixos n2 e n3
	if(!strcmp(name.c_str(),"Translacoes na superficie de referencia (u,v,w)")) return 15;// translacoes u,v,w
	
	
	
	return TPZMaterial::VariableIndex(name);
}

/** returns the number of variables associated with the variable indexed by var.
 var is obtained by calling VariableIndex*/
int TPZMatPlaca2::NSolutionVariables(int var){
	if(var == 2) return 1;  // Desloc. na superficie de referencia na eixo x global
	if(var == 3) return 1;  // Desloc. na superficie de referencia na eixo y global
	if(var == 4) return 1;  // Desloc. na superficie de referencia na eixo z global
	if(var == 5) return 3;  // Momentos na superficie de referencia nas direcoes dos eixos n1 e n2
	if(var== 50) return 3;  // Momentos na superficie de referencia nas direcoes dos eixos a1 e a2
	//if(var == 6) return 1;
	//if(var == 7) return 1;
	if(var == 8) return 1;
	if(var == 9) return 1;
	if(var == 10) return 1;
	if(var == 11) return 1;
	if(var == 12) return 1;
	if(var == 54) return 3;
	if(var == 13) return 1;
	if(var == 14) return 1;
	if(var == 15) return 3;
	
	return TPZMaterial::NSolutionVariables(var);  // todos os Deslocamentos nodais e suas derivadas nas
	// direcoes dos eixos axes (a1 e a2)
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TPZMatPlaca2::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,
							TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){
	
	REAL k = 5./6.;
	
	if(var == 2) {
		Solout.Resize(1);
		Solout[0] = Sol[0];
		return;
	}
	if(var == 3) {
		Solout.Resize(1);
		Solout[0] = Sol[1];
		return;
	}
	if(var == 4) {
		Solout.Resize(1);
		Solout[0] = Sol[2];
		return;
	}
	if(var == 15) {
		Solout.Resize(3);
		Solout[0] = Sol[0];
		Solout[1] = Sol[1];
		Solout[2] = Sol[2];
		return;
	}
	
	
	TPZVec<REAL> Soln(6);
	
	TPZFMatrix DSolnax(2,6),DSolnn(2,6);
	
	
	int idf,jdf;
	for(idf=0; idf<6; idf++) {
		Soln[idf] = 0;
		DSolnax(0,idf) = 0;
		DSolnax(1,idf) = 0.;
		for(jdf=0; jdf<6; jdf++) {
			Soln[idf] += fRmat(idf,jdf)*Sol[jdf];
			DSolnax(0,idf) += fRmat(idf,jdf)*DSol(0,jdf);
			DSolnax(1,idf) += fRmat(idf,jdf)*DSol(1,jdf);
		}
	}
	
	
	//     REAL Dax1n1, Dax1n2, Dax2n1, Dax2n2;
	
    //         Dax1n1 = axes(0,0)* fnaxes(0,0) + axes(0,1)* fnaxes(0,1) + axes(0,2)* fnaxes(0,2);
    //         Dax1n2 = axes(0,0)* fnaxes(1,0) + axes(0,1)* fnaxes(1,1) + axes(0,2)* fnaxes(1,2);
    //         Dax2n1 = axes(1,0)* fnaxes(0,0) + axes(1,1)* fnaxes(0,1) + axes(1,2)* fnaxes(0,2);
    //         Dax2n2 = axes(1,0)* fnaxes(1,0) + axes(1,1)* fnaxes(1,1) + axes(1,2)* fnaxes(1,2);
	
	TPZFMatrix Rmatan(2,2);
	Rmatan(0,0)= axes(0,0)* fnaxes(0,0) + axes(0,1)* fnaxes(0,1) + axes(0,2)* fnaxes(0,2);
	Rmatan(1,0)= axes(0,0)* fnaxes(1,0) + axes(0,1)* fnaxes(1,1) + axes(0,2)* fnaxes(1,2);
	Rmatan(0,1)= axes(1,0)* fnaxes(0,0) + axes(1,1)* fnaxes(0,1) + axes(1,2)* fnaxes(0,2);
	Rmatan(1,1)= axes(1,0)* fnaxes(1,0) + axes(1,1)* fnaxes(1,1) + axes(1,2)* fnaxes(1,2);
	
	for(idf=0;idf<6;idf++) {
		DSolnn(0,idf) = Rmatan(0,0)*DSolnax(0,idf)+Rmatan(0,1)*DSolnax(1,idf);
		DSolnn(1,idf) = Rmatan(1,0)*DSolnax(0,idf)+Rmatan(1,1)*DSolnax(1,idf);
	}
	
	if (var==5){
		REAL Mn1;
		Mn1 = -fE1*fh*fh*fh*fmi*(-fni2*DSolnn(1,3)+DSolnn(0,4))/12.;
		Mn1 += fE1*ff*fh*fmi*(ff*fni2*DSolnn(1,3)-fni2*DSolnn(1,1)
							  -ff*DSolnn(0,4)-DSolnn(0,0));
		
		
		REAL Mn2;
		Mn2  =-(fE2*fh*fh*fh*fmi*(-DSolnn(1,3) + fni1*DSolnn(0,4)))/12.;
		Mn2 +=  fE2*ff*fh*fmi*(ff*DSolnn(1,3)  - DSolnn(1,1) -
							   ff*fni1*DSolnn(0,4) - fni1*DSolnn(0,0));
		
		REAL Mn1n2;
		Mn1n2  = fG12*fh*fh*fh*(DSolnn(1,4) - DSolnn(0,3))/24.;
		Mn1n2 += (ff*fG12*fh*(ff*DSolnn(1,4) + DSolnn(1,0) -
							  ff*DSolnn(0,3) + DSolnn(0,1)))/2.;
		
		Solout.Resize(3);
		Solout[0]=Mn1;
		Solout[1]=Mn2;
		Solout[2]=Mn1n2;
		return;
	}
	
	if (var==50){
		
		REAL Mn1;
		Mn1 = -fE1*fh*fh*fh*fmi*(-fni2*DSolnn(1,3)+DSolnn(0,4))/12.;
		Mn1 += fE1*ff*fh*fmi*(ff*fni2*DSolnn(1,3)-fni2*DSolnn(1,1)
							  -ff*DSolnn(0,4)-DSolnn(0,0));
		
		
		REAL Mn2;
		Mn2  =-(fE2*fh*fh*fh*fmi*(-DSolnn(1,3) + fni1*DSolnn(0,4)))/12.;
		Mn2 +=  fE2*ff*fh*fmi*(ff*DSolnn(1,3)  - DSolnn(1,1) -
							   ff*fni1*DSolnn(0,4) - fni1*DSolnn(0,0));
		
		REAL Mn1n2;
		Mn1n2  = fG12*fh*fh*fh*(DSolnn(1,4) - DSolnn(0,3))/24.;
		Mn1n2 += (ff*fG12*fh*(ff*DSolnn(1,4) + DSolnn(1,0) -
							  ff*DSolnn(0,3) + DSolnn(0,1)))/2.;
		REAL Ma1;
		Ma1 = Mn1 * Rmatan(0,0) * Rmatan(0,0) + Mn2*Rmatan(1,0) * Rmatan(1,0) + 2.* Rmatan(0,0) * Rmatan(1,0) * Mn1n2;
		REAL Ma2;
		Ma2 = Mn1 * Rmatan(0,1) * Rmatan(0,1) + Mn2*Rmatan(1,1) * Rmatan(1,1) + 2.* Rmatan(1,1) * Rmatan(0,1) * Mn1n2;
		
		REAL  Ma1a2;
		Ma1a2 = Mn1n2 * (Rmatan(0,0)*Rmatan(0,0) + Rmatan(1,0) *Rmatan(0,1))+ (Mn1-Mn2)*Rmatan(0,0)*Rmatan(0,1);
		
		Solout.Resize(3);
		Solout[0]=Ma1;
		Solout[1]=Ma2;
		Solout[2]=Ma1a2;
		return;
	}
	
	if (var == 8) {
		REAL Vn1;
		Vn1  =  (fG12*fh*k*(Soln[4]+DSolnn(0,2)));
		Solout[0] = Vn1;
		return;
	}
	if (var == 9) {
		REAL Vn2;
		Vn2  = (fG12*fh*k*(-Soln[3]+ DSolnn(1,2)))/2.;
		Solout[0] = Vn2;
		return;
	}
	
	REAL Sign1 = 0.;
	
	if (var == 10  || var == 54) {
		REAL z;
		z=0.;
		Sign1 = -fE1*fmi*(fni2*(-(ff + z)*DSolnn(1,3) + DSolnn(1,1))+
						  (ff + z)*DSolnn(0,4) + DSolnn(0,0));
		if(var == 10) {
			Solout[0] = Sign1;
			return;
		}
	}
	REAL Sign2 = 0.;
	if (var == 11 || var == 54) {
		REAL z;
		z=0.;
		Sign2=-fE2*fmi*(-(ff + z)*DSolnn(1,3) + DSolnn(1,1) +
						fni1*((ff + z)*DSolnn(0,4) + DSolnn(0,0)));
		if(var == 11) {
			Solout[0] = Sign2;
			return;
		}
	}
	REAL Taun1n2 = 0.;
	if (var == 12 || var == 54) {
		REAL z;
		z=0.;
		Taun1n2=(fG12*((ff + z)*DSolnn(1,4) + DSolnn(1,0) -
					   (ff + z)*DSolnn(0,3) + DSolnn(0,1)))/2.;
		if(var == 12) {
			Solout[0] = Taun1n2;
			return;
		}
	}
	
	if(var == 54) {
		REAL Siga1;
		Siga1 = Sign1 * Rmatan(0,0) * Rmatan(0,0) + Sign2*Rmatan(1,0) * Rmatan(1,0) + 2.* Rmatan(0,0) * Rmatan(1,0) * Taun1n2;
		REAL Siga2;
		Siga2 = Sign1 * Rmatan(0,1) * Rmatan(0,1) + Sign2*Rmatan(1,1) * Rmatan(1,1) + 2.* Rmatan(1,1) * Rmatan(0,1) * Taun1n2;
		
		REAL  Siga1a2;
		Siga1a2 = Taun1n2 * (Rmatan(0,0)*Rmatan(0,0) + Rmatan(1,0) *Rmatan(0,1))+ (Sign1-Sign2)*Rmatan(0,0)*Rmatan(0,1);
		
		Solout.Resize(3);
		Solout[0]=Siga1*fh;
		Solout[1]=Siga2*fh;
		Solout[2]=Siga1a2*fh;
		return;
	}
	if (var == 13 ) {
		REAL Taun1n3;
		Taun1n3=fG13*(Soln[4] + DSolnn(0,2))/2.;
		Solout[0] = Taun1n3;
		return;
	}
	if (var == 14 ) {
		REAL Taun2n3;
		Taun2n3=fG23*(-Soln[3] + DSolnn(1,2))/2.;
		Solout[0] = Taun2n3;
		return;
	}
	
	TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}



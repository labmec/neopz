/**
 * \file
 * @brief Contains implementations of the TPZMultPlaca methods.
 */
#include "pzmaterial.h"
//#include "pztempmat.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"
#include <math.h>
#include <fstream>
using namespace std;

#include "pzvec.h" 
#include "pzerror.h"

#include "pzmultplaca.h"

TPZMultPlaca::TPZMultPlaca(int num, REAL h, TPZVec<REAL> &esp, REAL f, REAL E1 ,
						   REAL E2 , REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
						   REAL G23 , TPZFMatrix<REAL> &naxes, TPZVec<REAL> &xf,
						   int camadaref, int camadaatual) :
TPZMatPlaca2(num, h, f, E1 , E2 , ni1 , ni2 , G12 , G13 ,
			 G23 , naxes, xf), fT(6,6,0.) {
	
	
	//variaveis de entrada
	/*  fIdfMax   = numero de graus de liberdade por ponto
	 icA       = camada atual
	 icR       = camada de referencia
	 icC       = camada corrente
	 esp       = vetor espessuras das camadas
	 f   = distancia do eixo de referencia ate o meio da camada atual
	 f>0 quando icA > icR     e    f<0 quando icA <icR
	 */
	fIdfMax = (3*(esp.NElements()+1));
	int icR(camadaref), icA(camadaatual), icC;
	
	
	
	// Preparacao da matriz T de ordem 6x(3(NCamadas+1))
	// que transforma as matrizes fKxx em KxxMC
	
	
	
	fT.Resize(6, fIdfMax);
	
	int imin,imax,idif,sign,i,j;
	
	if (icA != icR) {
		if (icA > icR) {
			idif=icA - icR -1;
			icC=icR;
			sign=1;
		}
		else {
			idif=icR-icA-1;
			icC=icA;
			sign=-1;
        };
		
		for (i=1; i<=idif; i++) {
			icC=icC+1;
			j=3*(icC+1);
			fT(0,j+1)  = sign*esp[icC];
			fT(1,j)=-sign*esp[icC];
		};
		
		j=3*(icA+1);
		fT(0,j+1)=-f+sign*esp[icA]/2.0;
		fT(1,j  )= f-sign*esp[icA]/2.0;
    };
	
    if (icA>icR){
		imin=3*(icR+1);
		imax=3*(icA+1);
	}
	else{
		imin=3*(icA+1);
		imax=3*(icR+1);
	};
    for (i=0;i<3;i++){
		fT(i,i)=1.0;
		fT(3+i,imin+i)=1.0;
		fT(3+i,imax+i)=1.0;
    };
	
	
	
	
	TPZFMatrix<REAL> TTransp;
	fT.Transpose(&TTransp);
	
	
	// Obtencao das matrizes do elemento no plano xy, com efeito multicamada
	
	TPZFMatrix<REAL> KxxMC,KyyMC,KxyMC,KyxMC,Bx0MC,B0xMC,By0MC,B0yMC,B00MC;
	
	KxxMC=TTransp*(fKxx*fT);
	KyyMC=TTransp*(fKyy*fT);
	KxyMC=TTransp*(fKxy*fT);
	Bx0MC=TTransp*(fBx0*fT);
	By0MC=TTransp*(fBy0*fT);
	B00MC=TTransp*(fB00*fT);
	
	KxyMC.Transpose(&KyxMC);
	Bx0MC.Transpose(&B0xMC);
	By0MC.Transpose(&B0yMC);
	
	
	
	// Geracao da Matriz de Rotacao para passar do plano para o espaco
	
	fRmat.Redim(fIdfMax,fIdfMax);
	fRmatT.Redim(fIdfMax,fIdfMax);
	
	imax= fIdfMax / 3;
	for (i=0; i<imax; i++){
		j=3*i;
		fRmat(j,j)   = fnaxes(0,0); fRmat(j,j+1)   = fnaxes(0,1); fRmat(j,j+2)   = fnaxes(0,2);
		fRmat(j+1,j) = fnaxes(1,0); fRmat(j+1,j+1) = fnaxes(1,1); fRmat(j+1,j+2) = fnaxes(1,2);
		fRmat(j+2,j) = fnaxes(2,0); fRmat(j+2,j+1) = fnaxes(2,1); fRmat(j+2,j+2) = fnaxes(2,2);
	}
	
	fRmat.Transpose(&fRmatT);
	
	
	// Obtencao da matriz auxiliar que ira entrar no calculo das matrizes Knn (do espaco)
	
	fKxxR.Redim(fIdfMax,fIdfMax); fKyyR.Redim(fIdfMax,fIdfMax);
	fKxyR.Redim(fIdfMax,fIdfMax); fKyxR.Redim(fIdfMax,fIdfMax);
	fBx0R.Redim(fIdfMax,fIdfMax); fB0xR.Redim(fIdfMax,fIdfMax);
	fBy0R.Redim(fIdfMax,fIdfMax); fB0yR.Redim(fIdfMax,fIdfMax);
	fB00R.Redim(fIdfMax,fIdfMax);
	
	fKxxR = fRmatT * (KxxMC * fRmat);
	fKyxR = fRmatT * (KyxMC * fRmat);
	fKxyR = fRmatT * (KxyMC * fRmat);
	fKyyR = fRmatT * (KyyMC * fRmat);
	fB0xR = fRmatT * (B0xMC * fRmat);
	fB0yR = fRmatT * (B0yMC * fRmat);
	fBx0R = fRmatT * (Bx0MC * fRmat);
	fBy0R = fRmatT * (By0MC * fRmat);
	fB00R = fRmatT * (B00MC * fRmat);
	
	// testes
	
	ofstream out("saida.dat");
	//  T.Print("Matriz T",out);
	//  KxxMC.Print("Matriz KxxMC",out);
	//  KyyMC.Print("Matriz KyyMC",out);
	//  KxyMC.Print("Matriz KxyMC",out);
	//  Bx0MC.Print("Matriz Bx0MC",out);
	//  By0MC.Print("Matriz By0MC",out);
	//  B00MC.Print("Matriz B00MC",out);
	
	//  KyxMC.Print("Matriz KyxMC",out);
	//  B0xMC.Print("Matriz B0xMC",out);
	//  B0yMC.Print("Matriz B0yMC",out);
	
	TPZFMatrix<REAL> Temp(fIdfMax,fIdfMax),Transp(fIdfMax,fIdfMax);
	
	fKxxR.Transpose(&Transp);
	Temp= fKxxR - Transp;
	Temp.Print("fKxxR - fKxxR^t",out);
	
	fKyyR.Transpose(&Transp);
	Temp= fKyyR - Transp;
	Temp.Print("fKyyR - fKyyR^t",out);
	
	fKyxR.Transpose(&Transp);
	Temp= fKxyR - Transp;
	Temp.Print("fKxyR - fKyxR^t",out);
	
	fKxyR.Transpose(&Transp);
	Temp= fKyxR - Transp;
	Temp.Print("fKyxR - fKxyR^t",out);
	
	fB0xR.Transpose(&Transp);
	Temp= fBx0R - Transp;
	Temp.Print("fBx0R - fB0xR^t",out);
	
	fBx0R.Transpose(&Transp);
	Temp= fB0xR - Transp;
	Temp.Print("fB0xR - fBx0R^t",out);
	
	fB0yR.Transpose(&Transp);
	Temp= fBy0R - Transp;
	Temp.Print("fBy0R - fB0yR^t",out);
	
	fBy0R.Transpose(&Transp);
	Temp= fB0yR - Transp;
	Temp.Print("fB0yR - fKy0R^t",out);
	
	fB00R.Transpose(&Transp);
	Temp= fB00R - Transp;
	Temp.Print("fB00R - fB00R^t",out);
};

/**returns the solution associated with the var index based on the finite element approximation*/
void TPZMultPlaca::Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,
							TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout){
	
	//        REAL k = 5./6.;
	
	if(var == 2 || var ==3 || var == 4) {
		TPZMatPlaca2::Solution(Sol,DSol,axes,var,Solout);
		return;
	}
	
	if(var > 4) {
		TPZVec<REAL> Soln(fIdfMax);
		
		TPZFMatrix<REAL> DSolnax(2,fIdfMax),DSolnn(2,fIdfMax);
		
		
		int idf,jdf;
		for(idf=0; idf<fIdfMax; idf++) {
			Soln[idf] = 0;
			DSolnax(0,idf) = 0;
			DSolnax(1,idf) = 0.;
			for(jdf=0; jdf<fIdfMax; jdf++) {
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
		
		TPZFMatrix<REAL> Rmatan(2,2);
		Rmatan(0,0)= axes(0,0)* fnaxes(0,0) + axes(0,1)* fnaxes(0,1) + axes(0,2)* fnaxes(0,2);
		Rmatan(0,1)= axes(0,0)* fnaxes(1,0) + axes(0,1)* fnaxes(1,1) + axes(0,2)* fnaxes(1,2);
		Rmatan(1,0)= axes(1,0)* fnaxes(0,0) + axes(1,1)* fnaxes(0,1) + axes(1,2)* fnaxes(0,2);
		Rmatan(1,1)= axes(1,0)* fnaxes(1,0) + axes(1,1)* fnaxes(1,1) + axes(1,2)* fnaxes(1,2);
		
		for(idf=0;idf<fIdfMax;idf++) {
			DSolnn(0,idf) = Rmatan(0,0)*DSolnax(0,idf)+Rmatan(0,1)*DSolnax(1,idf);
			DSolnn(1,idf) = Rmatan(1,0)*DSolnax(0,idf)+Rmatan(1,1)*DSolnax(1,idf);
		}
		
		TPZVec<REAL> SolStarnn(6);
		TPZFMatrix<REAL> DSolStarnn(2,6);
		
		for (idf=0; idf<6; idf++){
			SolStarnn[idf]=0.;
			DSolStarnn(0,idf)=0.;
			DSolStarnn(1,idf)=0.;
			for (jdf=0; jdf<fIdfMax; jdf++){
				SolStarnn[idf]    += fT(idf,jdf)*Soln[jdf];
				DSolStarnn(0,idf) += fT(idf,jdf)*DSolnn(0,jdf);
				DSolStarnn(1,idf) += fT(idf,jdf)*DSolnn(1,jdf);
			}
		}
		TPZVec<REAL> Sol6(6,0.);
		TPZFMatrix<REAL> DSol6(2,6,0.),DSoln6a(2,6,0.);
		for(idf=0;idf<6;idf++) {
			DSoln6a(0,idf) = Rmatan(0,0)*DSolStarnn(0,idf)+Rmatan(1,0)*DSolStarnn(1,idf);
			DSoln6a(1,idf) = Rmatan(0,1)*DSolStarnn(0,idf)+Rmatan(1,1)*DSolStarnn(1,idf);
		}
		for(idf=0; idf<6; idf++) {
			Sol6[idf] = 0;
			DSol6(0,idf) = 0;
			DSol6(1,idf) = 0.;
			for(jdf=0; jdf<6; jdf++) {
				Sol6[idf] += fRmat(jdf,idf)*SolStarnn[jdf];
				DSol6(0,idf) += fRmat(jdf,idf)*DSoln6a(0,jdf);
				DSol6(1,idf) += fRmat(jdf,idf)*DSoln6a(1,jdf);
			}
		}
		
		
		TPZMatPlaca2::Solution(Sol6,DSol6,axes,var,Solout);
    }
	return;
	
}

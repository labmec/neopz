#include "TPZShapeDisc.h"
#include "pzshapelinear.h"
#include "pzreal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void (*TPZShapeDisc::fOrthogonal)(REAL C, REAL x0, REAL x,int degree, TPZFMatrix & phi, TPZFMatrix & dphi, int n) = TPZShapeDisc::Polynomial;

TPZShapeDisc::TPZShapeDisc(){
}

void TPZShapeDisc::Polynomial(REAL C,REAL x0,REAL x,int degree,TPZFMatrix &phi,TPZFMatrix &dphi, int n){
   
   if(degree < 0){
      PZError << "TPZShapeDisc::Polynomial the degree of the polynomial cannot be minus, aborting\n";
      exit(-1);
   }
   phi.Redim(degree+1,1);
   dphi.Redim(n,degree+1);
   //grau zero ou constante
   phi(0,0) = 1.0;
   if(degree == 0) return;
   //grau 1
   REAL val = (x-x0)/C;
   phi(1,0) = val;
   dphi(0,1) = 1./C;
   //grau maior que 1
   int p;
   degree++;
   for(p = 2;p < degree; p++) {
      phi(p,0) = phi(p-1,0)*val;
      int ideriv;
      dphi(0,p) = dphi(0,p-1)*val+phi(p-1,0)/C;
      for(ideriv=2; ideriv<=n; ideriv++) {
	 dphi(ideriv-1,p) = val*dphi(ideriv-1,p-1)+(ideriv/C)*dphi(ideriv-2,p-1);
      }
   }
}

void TPZShapeDisc::Legendre(REAL C,REAL x0,REAL x,int degree,TPZFMatrix & phi,TPZFMatrix & dphi, int n){
   
   x = (x - x0) / C;
   phi.Redim(degree+1,1);
   dphi.Redim(n,degree+1);

   TPZShapeLinear::Legendre(x, degree+1, phi, dphi, n);

   REAL val;
   for(int ord = 0; ord < degree+1; ord++)
      for (int ideriv = 0; ideriv < n; ideriv++){
	 val = dphi(ideriv, ord) / pow(C, ideriv+1);
	 dphi.Put(ideriv, ord, val);
      }

} //end of method


void  TPZShapeDisc::Shape0D(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi){

  if(degree < 0){
    PZError << "TPZShapeDisc::Polynomial the degree of the polynomial cannot be minus, aborting\n";
    exit(-1);
  }
  int cap = degree+1;
  phi.Redim(cap,1);
  for(int i=0;i<cap;i++) phi.PutVal(i,0,1.0);//fun¢ão unitária
  dphi.Redim(0,0);//não existe a derivada em um ponto: dimensão nula
} 

void  TPZShapeDisc::Shape1D(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi){

  //suponha a componente não nula sendo a primeira
  REAL x0 = X0[0];
  REAL x = X[0];

  fOrthogonal(C,x0,x,degree,phi,dphi,1);

}

void TPZShapeDisc::Shape2D(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi, MShapeType type){

  REAL x0 = X0[0];
  REAL y0 = X0[1];
  REAL x = X[0];
  REAL y = X[1];

  TPZFMatrix phix,phiy,dphix,dphiy;

  fOrthogonal(C,x0,x,degree,phix,dphix);
  fOrthogonal(C,y0,y,degree,phiy,dphiy);

  int i,j,counter=0,num = degree+1;
  
  int nshape = NShapeF(degree,2,type);
//  int count=num,ind=0;
  phi.Redim(nshape,1);
  dphi.Redim(2,nshape);
  phi.Zero();
  dphi.Zero();

  //valor da função
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree && type == EOrdemTotal) break;
      phi(counter,0) = phix(i,0)*phiy(j,0);
      dphi(0,counter) = dphix(0,i)*phiy(j,0);
      dphi(1,counter) = phix(i,0)*dphiy(0,j);
      counter++;
    }
  }
  if(counter != nshape) {
    PZError << "TPZShapeDisc::Shape2D wrong shape count\n";
  }
  REAL phi0,dphi0[2];
  phi0 = phi(0);
  dphi0[0] = dphi(0,0);
  dphi0[1] = dphi(1,0);
  phi(0) = phi(nshape-1);
  dphi(0,0) = dphi(0,nshape-1);
  dphi(1,0) = dphi(1,nshape-1);

  phi(nshape-1) = phi0;
  dphi(0,nshape-1) = dphi0[0];
  dphi(1,nshape-1) = dphi0[1];
}

void  TPZShapeDisc::Shape3D(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi, MShapeType type){

  REAL x0 = X0[0];
  REAL y0 = X0[1];
  REAL z0 = X0[2];
  REAL x = X[0];
  REAL y = X[1];
  REAL z = X[2];

  TPZFMatrix phix,phiy,phiz,dphix,dphiy,dphiz;


  fOrthogonal(C,x0,x,degree,phix,dphix);
  fOrthogonal(C,y0,y,degree,phiy,dphiy);
  fOrthogonal(C,z0,z,degree,phiz,dphiz);
  
 
  int i,j,k,nshape = 0,num=degree+1, counter=0;
  for(i=0;i<num;i++) nshape += (i+1)*(i+2)/2;
//  int count=num,ind=0;
  phi.Redim(nshape,1);
  dphi.Redim(3,nshape);
  phi.Zero();
  dphi.Zero();

  //function value
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree && type == EOrdemTotal) break;
      for(k=0;k<num;k++){
	if(i+j+k > degree && type == EOrdemTotal) break;
	phi(counter,0) = phix(i,0)*phiy(j,0)*phiz(k,0);
	dphi(0, counter) = dphix(0,i)*phiy(j,0)*phiz(k,0);
	dphi(1, counter ) = phix(i,0)*dphiy(0,j)*phiz(k,0);
	dphi(2, counter ) = phix(i,0)*phiy(j,0)*dphiz(0,k);
	counter++;
      }
    }
  }

}


void  TPZShapeDisc::Shape2DFull(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi, MShapeType type){

  REAL x0 = X0[0];
  REAL y0 = X0[1];
  REAL x = X[0];
  REAL y = X[1];


  TPZFNMatrix<30> phix(degree+1,1),phiy(degree+1,1),dphix(3,degree+1),dphiy(3,degree+1);

  fOrthogonal(C,x0,x,degree,phix,dphix,3);
  fOrthogonal(C,y0,y,degree,phiy,dphiy,3);
 
  int i,j,counter=0,num = degree+1;
  
  int nshape = NShapeF(degree,2,type);
//  int count=num,ind=0;
  phi.Redim(nshape,1);
  dphi.Redim(5,nshape);
  phi.Zero();
  dphi.Zero();

  //valor da função
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree && type == EOrdemTotal) break;
      phi(counter,0) = phix(i,0)*phiy(j,0);
      dphi(0,counter) = dphix(0,i)*phiy(j,0);
      dphi(1,counter) = phix(i,0)*dphiy(0,j);
      dphi(2,counter) = dphix(1,i)*phiy(j,0)+phix(i,0)*dphiy(1,j); // Laplaciano
      dphi(3,counter) = dphix(2,i)*phiy(j,0)+dphix(0,i)*dphiy(1,j); // D/Dx Laplaciano
      dphi(4,counter) = dphix(1,i)*dphiy(0,j)+phix(i,0)*dphiy(2,j); // D/Dy Laplaciano
      counter++;
    }
  }
  if(counter != nshape) {
    PZError << "TPZShapeDisc::Shape2D wrong shape count\n";
  }
  REAL phi0,dphi0[5];
  phi0 = phi(0);
  for(i=0; i<5; i++) dphi0[i] = dphi(i,0);
  phi(0) = phi(nshape-1);
  for(i=0; i<5; i++) dphi(i,0) = dphi(i,nshape-1);

  phi(nshape-1) = phi0;
  for(i=0; i<5; i++) dphi(i,nshape-1) = dphi0[i];
}

int  TPZShapeDisc::NShapeF(int degree, int dimension, MShapeType type) {
  int sum =0,i;
  switch(dimension) {
  case 1:
    return degree+1;
    break;
  case 2:
    if(type == ETensorial) return (degree+1)*(degree+1);
    else return (degree+1)*(degree+2)/2;
    break;
  case 3:
    if(type == ETensorial) return (degree+1)*(degree+1)*(degree+1);
    for(i=0;i<(degree+1);i++) sum += (i+1)*(i+2)/2;
    return sum;
    break;
  default:
      PZError << "TPZShapeDisc::NShapeF case does not exists\n";
      return -1;
  }
}


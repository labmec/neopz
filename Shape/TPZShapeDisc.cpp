#include "TPZShapeDisc.h"
#include "pzreal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>




TPZShapeDisc::TPZShapeDisc(){

}

void TPZShapeDisc::Polynomial(REAL C,REAL x0,REAL x,int degree,TPZFMatrix &phi,TPZFMatrix &dphi){

  phi.Redim(degree+1,1);
  dphi.Redim(1,degree+1);

  if(degree < 0) return;
  //grau zero ou constante
  phi.Put(0,0,1.0);
  dphi.Put(0,0, 0.0);
  if(degree == 0) return;
  //grau 1
  REAL val = (x-x0)/C;
  REAL deriv = 1./C;
  phi.Put(1,0, val);
  dphi.Put(0,1, deriv);
  //grau maior que 1
  int p;
  degree++;
  for(p = 2;p < degree; p++) {
    REAL q = (REAL) p;
    val = pow((x-x0)/C,q);
    phi.Put(p,0, val);
    deriv = q*pow((x-x0),q-1.)/pow(C,q);
    dphi.Put(0,p, deriv);
  }
}


void  TPZShapeDisc::Shape1D(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi){

  //suponha a componente não nula sendo a primeira
  REAL x0 = X0[0];
  REAL x = X[0];
  Polynomial(C,x0,x,degree,phi,dphi);
}

void  TPZShapeDisc::Shape2D(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi){

  REAL x0 = X0[0];
  REAL y0 = X0[1];
  REAL x = X[0];
  REAL y = X[1];

  TPZFMatrix phix,phiy,dphix,dphiy;
  Polynomial(C,x0,x,degree,phix,dphix);
  Polynomial(C,y0,y,degree,phiy,dphiy);
  
  int i,j,num = degree+1;
  int nshape = num*(num+1)/2;
//  int count=num,ind=0;
  int k,l;
  phi.Redim(nshape,1);
  dphi.Redim(2,nshape);
  phi.Zero();
  dphi.Zero();

  //valor da função
  k = 0;
  l = 0;
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree) break;
      if(l == num){
	k++;
	l = 0;
      }
      phi(num*k + l++,0) = phix(i,0)*phiy(j,0);
    }
  }

  //derivada com respeito a x
  k = 0;
  l = 0;
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree) break;
      if(l == num){
	k++;
	l = 0;
      }
      dphi(0,num*k + l++) = dphix(0,i)*phiy(j,0);
    }
  }

  //derivada com respeito a y
  k = 0;
  l =  0;
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree) break;
      if(l == num){
	k++;
	l = 0;
      }
      dphi(1,num*k + l++) = phix(i,0)*dphiy(0,j);
    }
  }
}

void  TPZShapeDisc::Shape3D(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi){

  REAL x0 = X0[0];
  REAL y0 = X0[1];
  REAL z0 = X0[2];
  REAL x = X[0];
  REAL y = X[1];
  REAL z = X[2];

  TPZFMatrix phix,phiy,phiz,dphix,dphiy,dphiz;
  Polynomial(C,x0,x,degree,phix,dphix);
  Polynomial(C,y0,y,degree,phiy,dphiy);
  Polynomial(C,z0,z,degree,phiz,dphiz);
  
  int i,j,k,nshape = 0,num=degree+1;
  for(i=0;i<num;i++) nshape += (i+1)*(i+2)/2;
//  int count=num,ind=0;
  phi.Redim(nshape,1);
  dphi.Redim(3,nshape);
  phi.Zero();
  dphi.Zero();

  //function value
  int l=0,m=0,n=0;
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree) break;
      for(k=0;k<num;k++){
	if(i+j+k > degree) break;
	if(n == num){
	  m++;
	  n = 0;
	  if(m==num){
	    l++;
	    m = 0;
	  }	  
	}
	phi( num * (num*l + m) + n++ ,0) = phix(i,0)*phiy(j,0)*phiz(k,0);
      }
    }
  }

  //derivada com respeito a x
  l=0,m=0,n=0;
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree) break;
      for(k=0;k<num;k++){
	if(i+j+k > degree) break;
	if(n == num){
	  m++;
	  n = 0;
	  if(m==num){
	    l++;
	    m = 0;
	  }	  
	}
	dphi(0, num * (num*l + m) + n++ ) = dphix(0,i)*phiy(j,0)*phiz(k,0);
      }
    }
  }

  //derivada com respeito a y
  l=0,m=0,n=0;
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree) break;
      for(k=0;k<num;k++){
	if(i+j+k > degree) break;
	if(n == num){
	  m++;
	  n = 0;
	  if(m==num){
	    l++;
	    m = 0;
	  }	  
	}
	dphi(1, num * (num*l + m) + n++ ) = phix(i,0)*dphiy(0,j)*phiz(k,0);
      }
    }
  }

  //derivada com respeito a z
  l=0,m=0,n=0;
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree) break;
      for(k=0;k<num;k++){
	if(i+j+k > degree) break;
	if(n == num){
	  m++;
	  n = 0;
	  if(m==num){
	    l++;
	    m = 0;
	  }	  
	}
	dphi(2, num * (num*l + m) + n++ ) = phix(i,0)*phiy(j,0)*dphiz(0,k);
      }
    }
  } 
}

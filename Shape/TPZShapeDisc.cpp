/**
 * @file
 * @brief Contains the implementation of the TPZShapeDisc methods.
 */
#include "TPZShapeDisc.h"
#include "pzshapelinear.h"
#include "pzreal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>

using namespace std;

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {

	void (*TPZShapeDisc::fOrthogonal)(REAL C, REAL x0, REAL x,int degree, TPZFMatrix<REAL> & phi, TPZFMatrix<REAL> & dphi, int n) = TPZShapeDisc::Polynomial;

TPZShapeDisc::TPZShapeDisc(){

}

TPZShapeDisc::TPZShapeDisc(const TPZShapeDisc &copy){

}

TPZShapeDisc::~TPZShapeDisc(){

}


void TPZShapeDisc::IntegratedLegendre(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, int n){
#ifdef DEBUG
  if(IsZero(C)){
    DebugStop();
  }
#endif
  const REAL L = C;
  if(IsZero(L)){
    DebugStop();
  }
  const REAL b = 2.*(x - x0)/L;

  if(degree >= 17) DebugStop();

  REAL shapes[17] = {
    1,
    b,
    pow(b,2)/2.,
    -b/2. + pow(b,3)/2.,
    (-3*pow(b,2))/4. + (5*pow(b,4))/8.,
    (3*b)/8. - (5*pow(b,3))/4. + (7*pow(b,5))/8.,
    (15*pow(b,2))/16. - (35*pow(b,4))/16. + (21*pow(b,6))/16.,
   (-5*b)/16. + (35*pow(b,3))/16. - (63*pow(b,5))/16. + (33*pow(b,7))/16.,
   (-35*pow(b,2))/32. + (315*pow(b,4))/64. - (231*pow(b,6))/32. + (429*pow(b,8))/128.,
   (35*b)/128. - (105*pow(b,3))/32. + (693*pow(b,5))/64. - (429*pow(b,7))/32. + (715*pow(b,9))/128.,
   (315*pow(b,2))/256. - (1155*pow(b,4))/128. + (3003*pow(b,6))/128. - (6435*pow(b,8))/256. + (2431*pow(b,10))/256.,
   (-63*b)/256. + (1155*pow(b,3))/256. - (3003*pow(b,5))/128. + (6435*pow(b,7))/128. - (12155*pow(b,9))/256. + (4199*pow(b,11))/256.,
   (-693*pow(b,2))/512. + (15015*pow(b,4))/1024. - (15015*pow(b,6))/256. + (109395*pow(b,8))/1024. - (46189*pow(b,10))/512. + (29393*pow(b,12))/1024.,
   (231*b)/1024. - (3003*pow(b,3))/512. + (45045*pow(b,5))/1024. - (36465*pow(b,7))/256. + (230945*pow(b,9))/1024. - (88179*pow(b,11))/512. + (52003*pow(b,13))/1024.,
   (3003*pow(b,2))/2048. - (45045*pow(b,4))/2048. + (255255*pow(b,6))/2048. - (692835*pow(b,8))/2048. + (969969*pow(b,10))/2048. - (676039*pow(b,12))/2048. + (185725*pow(b,14))/2048.,
   (-429*b)/2048. + (15015*pow(b,3))/2048. - (153153*pow(b,5))/2048. + (692835*pow(b,7))/2048. - (1616615*pow(b,9))/2048. + (2028117*pow(b,11))/2048. - (1300075*pow(b,13))/2048. + (334305*pow(b,15))/2048.,
   (-6435*pow(b,2))/4096. + (255255*pow(b,4))/8192. - (969969*pow(b,6))/4096. + (14549535*pow(b,8))/16384. - (7436429*pow(b,10))/4096. + (16900975*pow(b,12))/8192. -(5014575*pow(b,14))/4096. + (9694845*pow(b,16))/32768.
   };

  REAL dshapes[17] = {
    0,1,b,-0.5 + (3*pow(b,2))/2.,(-3*b)/2. + (5*pow(b,3))/2.,0.375 - (15*pow(b,2))/4. + (35*pow(b,4))/8.,(15*b)/8. - (35*pow(b,3))/4. + (63*pow(b,5))/8.,
   -0.3125 + (105*pow(b,2))/16. - (315*pow(b,4))/16. + (231*pow(b,6))/16.,(-35*b)/16. + (315*pow(b,3))/16. - (693*pow(b,5))/16. + (429*pow(b,7))/16.,
   0.2734375 - (315*pow(b,2))/32. + (3465*pow(b,4))/64. - (3003*pow(b,6))/32. + (6435*pow(b,8))/128.,
   (315*b)/128. - (1155*pow(b,3))/32. + (9009*pow(b,5))/64. - (6435*pow(b,7))/32. + (12155*pow(b,9))/128.,
   -0.24609375 + (3465*pow(b,2))/256. - (15015*pow(b,4))/128. + (45045*pow(b,6))/128. - (109395*pow(b,8))/256. + (46189*pow(b,10))/256.,
   (-693*b)/256. + (15015*pow(b,3))/256. - (45045*pow(b,5))/128. + (109395*pow(b,7))/128. - (230945*pow(b,9))/256. + (88179*pow(b,11))/256.,
   0.2255859375 - (9009*pow(b,2))/512. + (225225*pow(b,4))/1024. - (255255*pow(b,6))/256. + (2078505*pow(b,8))/1024. - (969969*pow(b,10))/512. + (676039*pow(b,12))/1024.,
   (3003*b)/1024. - (45045*pow(b,3))/512. + (765765*pow(b,5))/1024. - (692835*pow(b,7))/256. + (4849845*pow(b,9))/1024. - (2028117*pow(b,11))/512. + (1300075*pow(b,13))/1024.,
   -0.20947265625 + (45045*pow(b,2))/2048. - (765765*pow(b,4))/2048. + (4849845*pow(b,6))/2048. - (14549535*pow(b,8))/2048. + (22309287*pow(b,10))/2048. - (16900975*pow(b,12))/2048. +
    (5014575*pow(b,14))/2048.,(-6435*b)/2048. + (255255*pow(b,3))/2048. - (2909907*pow(b,5))/2048. + (14549535*pow(b,7))/2048. - (37182145*pow(b,9))/2048. + (50702925*pow(b,11))/2048. -
    (35102025*pow(b,13))/2048. + (9694845*pow(b,15))/2048. };

   const REAL dbdx = 2./L;


  if(degree < 0 || n > 1){
    DebugStop();
  }
  phi.Redim(degree+1,1);
  dphi.Redim(n,degree+1);
  for(int i = 0; i < degree+1; i++){
    phi(i,0) = shapes[i];
    dphi(0,i) = dshapes[i]*dbdx;
  }

}//method

void TPZShapeDisc::Polynomial(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, int n){
#ifdef DEBUG
  if(IsZero(C)){
    DebugStop();
  }
#endif
   if(degree < 0){
      PZError << "TPZShapeDisc::Polynomial the degree of the polynomial cannot be minus, aborting\n";
      DebugStop();
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

void TPZShapeDisc::PolynomialWithoutScale(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n){
#ifdef DEBUG
  if(IsZero(C)){
    DebugStop();
  }
#endif
   if(degree < 0){
      PZError << "TPZShapeDisc::Polynomial the degree of the polynomial cannot be minus, aborting\n";
      DebugStop();
   }
   phi.Redim(degree+1,1);
   dphi.Redim(n,degree+1);
   //grau zero ou constante
   phi(0,0) = 1.0;
   if(degree == 0) return;
   //grau 1
   REAL val = x;
   phi(1,0) = val;
   dphi(0,1) = 1.; 
   //grau maior que 1
   int p;
   degree++;
   for(p = 2;p < degree; p++) {
      phi(p,0) = phi(p-1,0)*val;
      int ideriv;
      dphi(0,p) = dphi(0,p-1)*val + phi(p-1,0)* 1.0;
      for(ideriv=2; ideriv<=n; ideriv++) {
	 dphi(ideriv-1,p) = val*dphi(ideriv-1,p-1)+(ideriv)*dphi(ideriv-2,p-1);
      }
   }         
}

void TPZShapeDisc::Legendre(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n){
#ifdef DEBUG
  if(IsZero(C)){
    DebugStop();
  }
#endif

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

void TPZShapeDisc::ChebyshevWithoutScale(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n){
#ifdef DEBUG
  if(IsZero(C)){
    DebugStop();
  }
#endif

   phi.Redim(degree+1,1);
   dphi.Redim(n,degree+1);

   TPZShapeLinear::Chebyshev(x, degree+1, phi, dphi);

} //end of method


void TPZShapeDisc::LegendreWithoutScale(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n){
#ifdef DEBUG
  if(IsZero(C)){
    DebugStop();
  }
#endif

   phi.Redim(degree+1,1);
   dphi.Redim(n,degree+1);

   TPZShapeLinear::Legendre(x, degree+1, phi, dphi, n);

} //end of method

void TPZShapeDisc::Shape(int dimension, REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type){
  switch (dimension){
    case(0) :{
      TPZShapeDisc::Shape0D(C,X0,X,degree,phi,dphi);
      break;
    }
    case(1) : {
      TPZShapeDisc::Shape1D(C,X0,X,degree,phi,dphi);
      break;
    }
    case(2) : {
      if(type == ETensorial || type == EOrdemTotal){
        TPZShapeDisc::Shape2D(C,X0,X,degree,phi,dphi,type);
      }
      else{
        TPZShapeDisc::Shape2DFull(C,X0,X,degree,phi,dphi,type);
      }
      break;
    }
    case(3) : {
      TPZShapeDisc::Shape3D(C,X0,X,degree,phi,dphi,type);
      break;
    }
    default:{
      PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
      PZError.flush();
      DebugStop();
    }
  }
}
		
		void TPZShapeDisc::Shape(int &dimension,int &degree,TPZVec<REAL> &X, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,MShapeType type)
		{
				REAL C=1;//fator de escala utilizado neste metodo
				TPZManVector<REAL,3> X0(3,0.);//centro do elemento
				
				if(type == ETensorial ||type == EOrdemTotal){
						TPZShapeDisc::Shape2D(C,X0,X,degree,phi,dphi,type);
				}
				else {
						std::cout<<"Not implement"<<std::endl;
						DebugStop();
				}
				
				
				
		}
		
    


void  TPZShapeDisc::Shape0D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){

  if(degree < 0){
    PZError << "TPZShapeDisc::Polynomial the degree of the polynomial cannot be minus, aborting\n";
    DebugStop();
  }
  int cap = degree+1;
  phi.Redim(cap,1);
		for(int i=0;i<cap;i++) phi.PutVal(i,0,1.0);//funÔøΩ unitÔøΩia
		dphi.Redim(0,0);//nÔøΩ existe a derivada em um ponto: dimensÔøΩ nula
}

void  TPZShapeDisc::Shape1D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){

		//suponha a componente nÔøΩ nula sendo a primeira
  REAL x0 = X0[0];
  REAL x = X[0];

  fOrthogonal(C,x0,x,degree,phi,dphi,1);

}

void TPZShapeDisc::Shape2D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type){

  if(degree == 0) {
    phi(0,0) = 1.;
    dphi(0,0) = 0.;
    dphi(1,0) = 0.;
    return;
  }
  REAL x0 = X0[0];
  REAL y0 = X0[1];
  REAL x = X[0];
  REAL y = X[1];

  TPZFNMatrix<660> phix,phiy,dphix,dphiy;

  fOrthogonal(C,x0,x,degree,phix,dphix,1);
  fOrthogonal(C,y0,y,degree,phiy,dphiy,1);

  int i, j, ix, iy, counter=0, num = degree+1;

  int nshape = NShapeF(degree,2,type);
//  int count=num,ind=0;
  phi.Redim(nshape,1);
  dphi.Redim(2,nshape);
  phi.Zero();
  dphi.Zero();

  if(type==EOrdemTotal)
  {
  //mounts the shape function according to the max order:
  //  ____
  // |///
  // |//
  // |/

    for(i = 0; i <= degree; i++)
    {
       for(j = 0; j <= i; j++)
       {
          // notice that ix+iy is always constant for given i
          ix = j;
	  iy = i - j;
          phi(counter,0) = phix(ix,0)*phiy(iy,0);
          dphi(0,counter) = dphix(0,ix)*phiy(iy,0);
          dphi(1,counter) = phix(ix,0)*dphiy(0,iy);
          counter++;
       }
    }
  }else{
  // mounts the shape functions according to each
  // direction order
  // The functions are assembled in this sequence:
  // x functions, y functions, max xy
  // (horizontal, vertical lines and the 0 in the graph)
  // _____
  // |0|||
  // |-0||
  // |--0|
  // |---0
    // type == ETensorial
    for(i=0;i<num;i++)
    {
      for(j=0;j<i;j++)
      {
         ix = j;
	 iy = i;
         phi(counter,0) = phix(ix,0)*phiy(iy,0);
         dphi(0,counter) = dphix(0,ix)*phiy(iy,0);
         dphi(1,counter) = phix(ix,0)*dphiy(0,iy);
	 counter++;
      }
      for(j=0;j<i;j++)
      {
         ix = i;
	 iy = j;
         phi(counter,0) = phix(ix,0)*phiy(iy,0);
         dphi(0,counter) = dphix(0,ix)*phiy(iy,0);
         dphi(1,counter) = phix(ix,0)*dphiy(0,iy);
	 counter++;
      }

      ix = i;
      iy = i;
      phi(counter,0) = phix(ix,0)*phiy(iy,0);
      dphi(0,counter) = dphix(0,ix)*phiy(iy,0);
      dphi(1,counter) = phix(ix,0)*dphiy(0,iy);
      counter++;
    }
  }

  if(counter != nshape) {
    PZError << "TPZShapeDisc::Shape2D wrong shape count\n";
  }
  REAL phi0,dphi0[2];
  phi0 = phi(0,0);//here
  dphi0[0] = dphi(0,0);
  dphi0[1] = dphi(1,0);
  phi(0,0) = phi(nshape-1,0);//here
  dphi(0,0) = dphi(0,nshape-1);
  dphi(1,0) = dphi(1,nshape-1);

  phi(nshape-1,0) = phi0;//here
  dphi(0,nshape-1) = dphi0[0];
  dphi(1,nshape-1) = dphi0[1];
}

void  TPZShapeDisc::Shape3D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type){

  REAL x0 = X0[0];
  REAL y0 = X0[1];
  REAL z0 = X0[2];
  REAL x = X[0];
  REAL y = X[1];
  REAL z = X[2];

  TPZFNMatrix<1000> phix,phiy,phiz,dphix,dphiy,dphiz;

  fOrthogonal(C,x0,x,degree,phix,dphix,1);
  fOrthogonal(C,y0,y,degree,phiy,dphiy,1);
  fOrthogonal(C,z0,z,degree,phiz,dphiz,1);

  int i,j,k,nshape = 0,num=degree+1, counter=0;
  //for(i=0;i<num;i++) nshape += (i+1)*(i+2)/2;
  nshape = NShapeF(degree,3,type);
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

  REAL phi0,dphi0[3];
  phi0 = phi(0);
  dphi0[0] = dphi(0,0);
  dphi0[1] = dphi(1,0);
  dphi0[2] = dphi(2,0);
  phi(0) = phi(nshape-1);
  dphi(0,0) = dphi(0,nshape-1);
  dphi(1,0) = dphi(1,nshape-1);
  dphi(2,0) = dphi(2,nshape-1);

  phi(nshape-1) = phi0;
  dphi(0,nshape-1) = dphi0[0];
  dphi(1,nshape-1) = dphi0[1];
  dphi(2,nshape-1) = dphi0[2];

}


void  TPZShapeDisc::Shape2DFull(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type){

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
  dphi.Redim(8,nshape);
  phi.Zero();
  dphi.Zero();

  //valor da fun�o
  for(i=0;i<num;i++){
    for(j=0;j<num;j++){
      if(i+j > degree && type == EOrdemTotalFull) break;
      phi(counter,0) = phix(i,0)*phiy(j,0);
      dphi(0,counter) = dphix(0,i)* phiy(j,0); // dx 
      dphi(1,counter) =  phix(i,0)*dphiy(0,j); // dy
      dphi(2,counter) = dphix(1,i)* phiy(j,0)+ phix(i,0)*dphiy(1,j); // Laplaciano
      dphi(3,counter) = dphix(2,i)* phiy(j,0)+dphix(0,i)*dphiy(1,j); // D/Dx Laplaciano
      dphi(4,counter) = dphix(1,i)*dphiy(0,j)+ phix(i,0)*dphiy(2,j); // D/Dy Laplaciano
      dphi(5,counter) = dphix(1,i)* phiy(j,0);   // dxx 
      dphi(6,counter) =  phix(i,0)*dphiy(1,j);   // dyy 
      dphi(7,counter) = dphix(0,i)*dphiy(0,j);   // dxy 
      counter++;
    }
  }
  if(counter != nshape) {
    PZError << "TPZShapeDisc::Shape2D wrong shape count\n";
  }
  REAL phi0,dphi0[8];
  phi0 = phi(0);
  for(i=0; i<8; i++) dphi0[i] = dphi(i,0);
  phi(0) = phi(nshape-1);
  for(i=0; i<8; i++) dphi(i,0) = dphi(i,nshape-1);

  phi(nshape-1) = phi0;
  for(i=0; i<8; i++) dphi(i,nshape-1) = dphi0[i];
}

int  TPZShapeDisc::NShapeF(int degree, int dimension, MShapeType type) {
  int sum =0,i;
  switch(dimension) {
  case 1:
    return degree+1;
//    break;
  case 2:
    if(type == ETensorial || type == ETensorialFull) return (degree+1)*(degree+1);
    else return (degree+1)*(degree+2)/2;
//    break;
  case 3:
    if(type == ETensorial || type == ETensorialFull) return (degree+1)*(degree+1)*(degree+1);
    for(i=0;i<(degree+1);i++) sum += (i+1)*(i+2)/2;
    return sum;
 //   break;
  default:
      PZError << "TPZShapeDisc::NShapeF case does not exists\n";
      return -1;
  }
//  return -1;
}

};

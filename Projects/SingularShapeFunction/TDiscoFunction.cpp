//$Id: TDiscoFunction.cpp,v 1.2 2008-02-09 21:12:07 tiago Exp $

#include "TDiscoFunction.h"

 TDiscoFunction:: TDiscoFunction(){

}

 TDiscoFunction::~ TDiscoFunction(){

}
    
void  TDiscoFunction::Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix &df){
  f.Resize(1);
  df.Resize(2,1);
  double Xc = 0.;
  double Yc = 0.;
  double Xp = x[0];
  double Yp = x[1];
  double r = sqrt( (Xp-Xc)*(Xp-Xc) +(Yp-Yc)*(Yp-Yc) );
  f[0] = log(r);
  df(0,0) = (Xp - Xc)/(r*r);
  df(1,0) = (Yp - Yc)/(r*r);
}  

int  TDiscoFunction::NFunctions(){
  return 1;
}
  
int  TDiscoFunction::PolynomialOrder(){
  return 100;
}


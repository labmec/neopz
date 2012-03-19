//$Id: TExtFunction.cpp,v 1.1 2008-02-06 18:18:44 tiago Exp $

#include "TExtFunction.h"

TExtFunction::TExtFunction(){

}

TExtFunction::~TExtFunction(){

}
    
void TExtFunction::Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix<REAL> &df){
  f.Resize(1);
  df.Resize(2,1);
  double r = 2.-x[0];
  if(r < 1e-6) r = 1e-6;
  f[0] = log(r);
  df(0,0) = -1./(r);
  df(1,0) = 0.;
}  

int TExtFunction::NFunctions(){
  return 1;
}
  
int TExtFunction::PolynomialOrder(){
  return 1;
}


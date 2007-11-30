

#include "pzfmatrix.h"
//#include "pztempmat.h" 
#include "pzskylmat.h"
#include "pzsolve.h"
#include "pzmat1dlin.h"
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include <stdlib.h>
#include <iostream>
#include "pzvec.h"

#include "pzmat2dlin.h"
#include "pzanalysis.h"
#include "pzmetis.h"
#include "pzmatplaca2.h"
#include "pzmultplaca.h"
#include "tpzmultcamada.h"

#include <stdio.h>
#include <time.h>
//template<class T>
//class TPZVec;
//#define NOTDEBUG

TPZMultCamada *LerMaterial(char *filename);

int main(){ 
  
  TPZMultCamada *mat1=LerMaterial("teste1.ttt");

  TPZMultCamada *mat2=LerMaterial("teste2.ttt");

  
  TPZVec<REAL> sol(6), solout1(3),solout2(3);
  TPZFMatrix dsol(2,6), axes(3,3);
  // valores arbitrarios para teste
   sol[0]=0.1 ; dsol(0,0)=0.11; dsol(1,0)=0.12;
   sol[1]=0.2; dsol(0,1)=0.13; dsol(1,1)=0.15;
   sol[2]=0.3; dsol(0,2)=0.16; dsol(1,2)=0.17;
   sol[3]=0.4; dsol(0,3)=0.25; dsol(1,3)=-0.9;
   sol[4]=0.4; dsol(0,4)=0.01; dsol(1,4)=-0.07;
   sol[5]=0.5; dsol(0,5)=1,07; dsol(1,5)=-0.33;

  // direcoes de integracao
  axes(0,0)=1.0; axes(1,0)=0.0; axes(2,0)=0.0;
  axes(0,1)=0.0; axes(1,1)=1.0; axes(2,1)=0.0;
  axes(0,2)=0.0; axes(1,2)=0.0; axes(2,2)=1.0;


  
  mat1->Solution(sol,dsol,axes,2,solout1);
  mat2->Solution(sol,dsol,axes,2,solout2);
 
  ofstream saida("testesai2.dat");
  REAL dif;
  saida << "Diferencas > 1E-11 para vetor ej \n";
  saida << "Desloc.  MAT1    MAT2. \n";
  int numocorrencias=0;
  for (int i=0; i<3; i++){
    dif=solout1[i]-solout2[i];
    if (fabs(dif) > 1E-11) {
      saida << i  << "  " << solout1[i] << "  " << solout2[i] << endl;
      numocorrencias += 1;
    }
  }
  saida.flush();
  cout << "Num de ocorrencias" << numocorrencias;

  delete mat1;
  delete mat2;
  return 0;
}


TPZMultCamada *LerMaterial(char *filename) {
  ifstream input(filename);
  TPZFMatrix naxes(3,3);
  REAL ni1,ni2,h,E1,E2,G12,G13,G23,f;
  REAL n00,n01,n02,n10,n11,n12,n20,n21,n22;
  TPZVec<REAL> xf(6), esp(1);
  REAL g,st,ct,PI=3.1415926536;
  int matindex, numcamadas,i,numbc,camadaref, materialtype;
  input >> matindex >> numcamadas >> camadaref >> materialtype;
  TPZMultCamada *matcamadas = new TPZMultCamada(matindex);
  esp.Resize(numcamadas);
  for(i=0; i<numcamadas; i++) input >> esp[i];
  for(i=0; i< numcamadas; i++) {
   input >> f  >>  h  >>
            E1 >>  E2 >>
            G12 >> G13 >> G23 >>
            ni1 >> ni2;
    // input >> n00 >> n01 >> n02 >> n10 >> n11 >> n12 >> n20 >> n21 >> n22;
    cout << "Angulo (graus) do eixo x ao eixo n = ";
    cin >> g;
    g = g*PI/180.;
    ct = cos(g);
    st = sin(g);
    n00=ct; n01=st; n02=0; n10=-st; n11=ct; n12=0.; n20=0.; n21=0.; n22=1.;
 
    if(materialtype == 1) {
      xf.Resize(3*(numcamadas+1));
      int i;
      for(i=0; i< 3*(numcamadas+1); i++) input >> xf[i];
    } else {
      input >> xf[0] >> xf[1] >> xf[2] >> xf[3] >> xf[4] >> xf[5];
    }
    naxes(0,0) =  n00;    naxes(0,1) =  n01;    naxes(0,2) =  n02;
    naxes(1,0) =  n10;    naxes(1,1) =  n11;    naxes(1,2) =  n12;
    naxes(2,0) =  n20;    naxes(2,1) =  n21;    naxes(2,2) =  n22;
    if(materialtype == 1) {
      matcamadas->AddLayer(new TPZMultPlaca(-1,h,esp,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf,camadaref,i));
    } else {
      matcamadas->AddLayer(new TPZMatPlaca2(-1,h,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf));
    }
  }
  return matcamadas;
}

  

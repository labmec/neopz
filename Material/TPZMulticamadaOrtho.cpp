#include "TPZPlacaOrthotropic.h"

#include "TPZMulticamadaOrtho.h"
#include "pzmatorthotropic.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
// #include "pzelmat.h"
// #include "pzbndcond.h"
// #include "pzmatrix.h"
// #include "pzfmatrix.h"
// #include "pzerror"
// #include "pztempmat.h"
// #include "pzmanvector.h"
// #include <math.h>
// #include <fstream>
// using namespace std;




TPZMulticamadaOrthotropic::TPZMulticamadaOrthotropic(REAL z,REAL dx,REAL dy, TPZGeoMesh *gmesh, TPZCompMesh *cmesh){

  fGeoMesh = gmesh;
  fCompMesh = cmesh;
  fZ  =  z;
  fDx = dx;
  fDy = dy;


}


void TPZMulticamadaOrthotropic::GenerateMesh(){

}


void TPZMulticamadaOrthotropic::AddPlacaOrtho(TPZPlacaOrthotropic *placa){
  
  fPlacaOrth.Push(placa);
  
}

void  TPZMulticamadaOrthotropic::Print(ostream &out){

  int i, nplaca=fPlacaOrth.NElements();
  out << "TPZMulticamadaOrthotropic::Print\n";
  out << nplaca << endl;
  for (i=0; i<nplaca; i++){
    //TPZMulticamadaOrthotropic *multcam;
    TPZPlacaOrthotropic *placa = fPlacaOrth[i];
    out << "placa : " << i << endl;
    fPlacaOrth[i]->Print();
   //out << "quantidade de camadas :" << multcam->ZHight(placa); 
  }
}
  
REAL TPZMulticamadaOrthotropic::ZHight(TPZPlacaOrthotropic *placa){

  int i,j;
  fZ = 0.;
  int quantplacas;

  TPZMulticamadaOrthotropic * multcam; 
  j = multcam->RQPlacas();
  for (i=0; i<j; i++){
    fZ += multcam->RPlacaOrtho()[i]->FH();
    
  }

  return fZ;
}

int TPZMulticamadaOrthotropic::RQPlacas(){

  int quantplacas;
  cout << "digitar a quantidade de placas desejadas : ";
  cin >> quantplacas;
  fQuantPlacas = quantplacas;
  return fQuantPlacas;
}



  


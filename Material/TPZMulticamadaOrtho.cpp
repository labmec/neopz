
#include "TPZMulticamadaOrtho.h"
#include "pzmatorthotropic.h"


TPZMulticamadaOrthotropic::TPZMulticamadaOrthotropic(REAL z,REAL dx,REAL dy){

  fZ  =  z;
  fDx = dx;
  fDy = dy;


}


void TPZMulticamadaOrthotropic::GenerateMesh(){}


int TPZMulticamadaOrthotropic::main(){
/**
 * dica: devera ler as carasteristica de cada placa a partir
 * de um arquivo de entrada e criar o material ortotropico
 * com aqueles dados (classe TPZMatOrthotropic), logo chamar o
 * método AddPlacaOrth() para adicionar esta placa na 
 * lista de placas, variável fPlacaOrth
 */ 

}

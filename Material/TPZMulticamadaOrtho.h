#ifndef MULTICAMADAORTH
#define MULTICAMADAORTH

#include <iostream>
using namespace std;
#include "pzvec.h"
#include "pzstack.h"

class TPZGeoMesh;
class TPZCompMesh;
class TPZPlacaOrthotropic;


/**
 * Gerencia um conjunto de placas
 * dispostas em forma multicamada
 */

class TPZMulticamadaOrthotropic {

  /**malha geométrica onde são inseridas as placas geométricas*/
  TPZGeoMesh             *fGeoMesh;
  /**malha computacional: elementos computacionais correspondentes*/
   TPZCompMesh            *fCompMesh;
  /**Vetor de placas*/
  //TPZVec<TPZPlacaOrthotropic *>  fPlacaOrth;
   TPZStack<TPZPlacaOrthotropic *> fPlacaOrth;
  /**
   * fZ: altura máxima (camada mais longe do plano XY)
   * fDx,fDy: dimensões das placas (constantes para todas as placas)   
   */
  REAL fDx,fDy;
  REAL fZ;
  double fQuantPlacas;

 public:
  /**construtor*/
  TPZMulticamadaOrthotropic(REAL z,REAL dx,REAL dy, TPZGeoMesh *gmesh, TPZCompMesh *cmesh);
  /*destrutor*/
  ~TPZMulticamadaOrthotropic(){}

  /**Adiciona placas ao conjunto*/
  void AddPlacaOrtho(TPZPlacaOrthotropic *placa);
  /**gera a malha computacional do conjunto de placas*/
  void GenerateMesh();
  /**cria o conjunto de placas multicamada*/
  //  static int main();
  void Print(ostream &out = cout);
  /*criando método para retornar a altura da multicamada*/
  REAL ZHight(TPZPlacaOrthotropic *placa);
  /*criando método para retornar fPlacaOrth*/
  TPZVec<TPZPlacaOrthotropic *> &RPlacaOrtho(){return fPlacaOrth;}
  /*criando método para contar quant de placas*/
  int RQPlacas();

};
#endif

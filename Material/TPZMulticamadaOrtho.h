#ifndef MULTICAMADAORTH
#define MULTICAMADAORTH

#include <iostream>
using namespace std;
#include "pzvec.h"
#include "pzstack.h"
#include "TPZPlacaOrthotropic.h"

class TPZGeoMesh;
class TPZCompMesh;
class TPZPlacaOrthotropic;
class TPZMatOrthotropic;


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
   TPZStack<TPZPlacaOrthotropic> fPlacaOrth;
  /**
   * fZ: altura máxima (camada mais longe do plano XY)
   * fDx,fDy: dimensões das placas (constantes para todas as placas)   
   */
  REAL fDx,fDy;
  /**
   * fNelx, fNely : Numero de elementos na direcao x e y
   */
  int fNelx, fNely;
  REAL fZMin, fZMax;
  //  double fQuantPlacas;

 public:
  /**construtor*/
  TPZMulticamadaOrthotropic(REAL z,REAL dx,REAL dy, int nelx, int nely);
  /*destrutor*/
  ~TPZMulticamadaOrthotropic(){}

  /**Adiciona placas ao conjunto*/
  void AddPlacaOrtho(TPZMatOrthotropic *material, REAL height);
  /**gera a malha computacional do conjunto de placas*/
  void GenerateMesh();
  /**cria o conjunto de placas multicamada*/
  //  static int main();
  void Print(ostream &out = cout);
  /*criando método para retornar a altura da multicamada*/
  REAL Height();
  /*criando método para retornar fPlacaOrth*/
  TPZVec<TPZPlacaOrthotropic> &RPlacaOrtho(){return fPlacaOrth;}
  /*criando método para contar quant de placas*/
  int NPlacas();

};
#endif

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
  REAL fMX[3],fMY[3],fMXY[3],fQX[3],fQY[3],fNX[3],fNY[3],fNXY[3];



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
  /**
   * Compute a tension state corresponding to the difference between the target state
   * and tension state loaded in the solution
   */
  void AnalyticTensor(TPZVec<REAL> &co, TPZFMatrix &tensor);

  /**
   * Tensor which needs to be applied at the given coordinate
   */
  void Tensor(TPZVec<REAL> &x, int placa, TPZFMatrix &tensor);
  /**
   * Computes the global efforts of the finite element solution
   */
  void ComputeCenterForces();

  void ComputeSolution();

  void SetMX(REAL MX) { 
    fMX[0] = MX;
  }

  void SetNX(REAL NX) {
    fNX[0] = NX;
  }

};
#endif

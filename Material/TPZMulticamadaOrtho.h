// -*- c++ -*-

// $Id: TPZMulticamadaOrtho.h,v 1.8 2003-11-07 00:40:34 phil Exp $
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

  void ComputeSolution(ostream &out = cout,int print = 0);

  void SetMX(REAL MX) { 
    fMX[0] = MX;
    fMX[2] = MX;
  }

  void SetNX(REAL NX) {
    fNX[0] = NX;
    fNX[2] = NX;
  }

  void SetMY(REAL MX) { 
    fMY[0] = MX;
    fMY[2] = MX;
  }

  void SetNXY(REAL NX) {
    fNXY[0] = NX;
    fNXY[2] = NX;
  }

  void SetMXY(REAL MXY) {
    fMXY[0] = MXY;
    fMXY[2] = MXY;
  }

  void SetQX(REAL QX) {
    fQX[0] = QX;
    fQX[2] = QX;
  }

  void SetQY(REAL QY) {
    fQY[0] = QY;
    fQY[2] = QY;
  }

  TPZGeoMesh *GeoMesh(){return fGeoMesh;}

  TPZCompMesh *CompMesh(){return fCompMesh;}
};
#endif

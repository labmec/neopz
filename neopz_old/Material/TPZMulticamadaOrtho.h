// -*- c++ -*-

// $Id: TPZMulticamadaOrtho.h,v 1.15 2007-01-03 00:08:26 phil Exp $
#ifndef MULTICAMADAORTH
#define MULTICAMADAORTH

#include <iostream>
#include "pzvec.h"
#include "pzstack.h"
#include "TPZPlacaOrthotropic.h"
#include "tpzautopointer.h"

class TPZGeoMesh;
class TPZCompMesh;
class TPZPlacaOrthotropic;
class TPZMatOrthotropic;
class TPZAnalysis;
class TPZMaterial;
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
  REAL fdMXdX[3],fdMXdY[3],fdMYdX[3],fdMYdY[3],fdMXYdX[3],fdMXYdY[3],
     fdQXdX[3],fdQXdY[3],fdQYdX[3],fdQYdY[3],fdNXdX[3],fdNXdY[3],fdNYdX[3],fdNYdY[3],fdNXYdX[3],fdNXYdY[3];
  int fLinearX,fLinearY;
  TPZManVector<REAL,3> fDirx, fDiry;

  /**
   * Relaxation factor to correct resulting forces.
   * @since Feb 10, 2004
   */
  REAL fCorrect;



 public:
  /**construtor*/
  TPZMulticamadaOrthotropic(REAL z,REAL dx,REAL dy, int nelx, int nely, REAL Correct = 1.0);
  /*destrutor*/
  ~TPZMulticamadaOrthotropic(){}

  /**Adiciona placas ao conjunto*/
  void AddPlacaOrtho(TPZAutoPointer<TPZMaterial> material, REAL height);
  /**gera a malha computacional do conjunto de placas*/
  void GenerateMesh();
  /**cria o conjunto de placas multicamada*/
  //  static int main();
  void Print(std::ostream &out = std::cout);
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

  void ComputeSolution(std::ostream &out = std::cout,int print = 0);

  void ComputeSolution(TPZMaterial *mat,std::ofstream &out,int numiter);

  void SetMX(REAL MX) { 
    fMX[0] = MX;
    fMX[2] = MX;
  }

  void SetNX(REAL NX) {
    fNX[0] = NX;
    fNX[2] = NX;
  }

  void SetNY(REAL NY) {
    fNY[0] = NY;
    fNY[2] = NY;
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
    fdMXdX[0] = QX;
    fdMXdX[2] = QX;
    if(QX != 0.) fLinearX = 1;
    else fLinearX = 0;
  }

  void SetQY(REAL QY) {
    fQY[0] = QY;
    fQY[2] = QY;
    fdMYdY[0] = QY;
    fdMYdY[2] = QY;
    if(QY != 0.) fLinearY = 1;
    else fLinearY = 0;

  }

  TPZGeoMesh *GeoMesh(){return fGeoMesh;}

  TPZCompMesh *CompMesh(){return fCompMesh;}

  void SetCorrect(REAL Correct){ fCorrect = Correct;}

  REAL CorrectFactor(){return fCorrect;}

  void PrintTensors(std::ostream &out);

  void PrintTensors(std::ostream &out,TPZFMatrix &tensorin,TPZFMatrix &tensorout);

  void PrintCenterForces(std::ostream &out);
};
#endif

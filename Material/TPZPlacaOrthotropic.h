// -*- c++ -*-
// $Id: TPZPlacaOrthotropic.h,v 1.6 2003-11-11 20:52:48 cedric Exp $
#ifndef PLACAORTHOTROPIC
#define PLACAORTHOTROPIC
#include "pzvec.h"
template<class>
class TPZVec;
class TPZFMatrix;
class TPZInterpolatedElement;
class TPZCompEl;
class TPZGeoEl;

// o objeto desta classe representa uma placa do objeto multicamada; 
//este útlimo representa fisicamente a superposi¢ão de várias placas. 
//não é o material em si.


class TPZPlacaOrthotropic {

 private:

  TPZGeoEl *fGeoEl;
  /**elemento computacional da placa*/
  TPZInterpolatedElement *fIntel;
  /**Espessura da placa*/
  REAL fH;
  REAL fZMin, fZMax;
  int fTensorVar;
  //double dx, dy;
  
 public:

  TPZPlacaOrthotropic();
  /**construtor da placa*/
  TPZPlacaOrthotropic(TPZGeoEl *gel,REAL zmin, REAL zmax);
  /*destrutor*/
  ~TPZPlacaOrthotropic(){}
  /**devolve o tensor de tensões da placa
   * @param ksi ponto no espaco parametrico
   */
  void Tensor(TPZVec<REAL> &ksi, TPZFMatrix &T);
  /**Dados dois vetores n1 e n2 retorna o momento*/
  REAL Moment(REAL zref, TPZVec<REAL> &normal, TPZVec<REAL> &direction);
  /**Dados dois vetores n1 e n2 retorna a forca*/
  REAL Force(TPZVec<REAL> &normal, TPZVec<REAL> &direction);

  /**
   * gradiente do momento na direcao indicada por graddir (a direcao e dada em espaco parametrico)
   * a derivada e devolvida em espaco real
   */
  REAL GradMoment(REAL zref, TPZVec<REAL> &graddir, TPZVec<REAL> &normal, TPZVec<REAL> &direction);

  /**
   * gradiente da forca na direcao indicada por graddir (a direcao e dada em espaco parametrico)
   * a derivada e devolvida em espaco real
   */
  REAL GradForce(TPZVec<REAL> &graddir, TPZVec<REAL> &normal, TPZVec<REAL> &direction);

  /**
   * gradiente do tensor na direcao indicada por graddir (a direcao e dada em espaco parametrico)
   * a derivada e devolvida em espaco real
   */
  void GradTensor(TPZVec<REAL> &graddir, TPZVec<REAL> &ksi,  TPZFMatrix &gradtensor);

  void PrintTensors(ostream &out);

  void PrintTensors(ostream &out,TPZFMatrix &tensorin,TPZFMatrix &tensorout);
  
  REAL Height(){return fH;}

  REAL ZMin() { return fZMin;}

  REAL ZMax() { return fZMax;}

  void IdentifyCompEl();

  void Print();

  TPZInterpolatedElement *ComputEl() {return fIntel;}
  //usa * no nome ComptEl, porque 
};
#endif

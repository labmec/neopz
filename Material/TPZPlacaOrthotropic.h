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
  /**devolve o tensor de tensões da placa*/
  void Tensor(REAL ksi, TPZFMatrix &T);
  /**Dados dois vetores n1 e n2 retorna o momento*/
  REAL Moment(REAL zref, TPZVec<REAL> &normal, TPZVec<REAL> &direction);
  /**Dados dois vetores n1 e n2 retorna a forca*/
  REAL Force(TPZVec<REAL> &normal, TPZVec<REAL> &direction);
  
  REAL Height(){return fH;}

  REAL ZMin() { return fZMin;}

  REAL ZMax() { return fZMax;}

  void IdentifyCompEl();

  void Print();

  TPZInterpolatedElement *ComputEl() {return fIntel;}
  //usa * no nome ComptEl, porque 
};
#endif

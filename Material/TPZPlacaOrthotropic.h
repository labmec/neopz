#ifndef PLACAORTHOTROPIC
#define PLACAORTHOTROPIC
#include "pzvec.h"
template<class>
class TPZVec;
class TPZFMatrix;
class TPZInterpolatedElement;
class TPZCompEl;

// o objeto desta classe representa uma placa do objeto multicamada; 
//este útlimo representa fisicamente a superposi¢ão de várias placas. 
//não é o material em si.


class TPZPlacaOrthotropic {

 private:

  /**elemento computacional da placa*/
  TPZInterpolatedElement *fIntel;
  /**Espessura da placa*/
  REAL fH;
  //double dx, dy;
  
 public:

  /**construtor da placa*/
  TPZPlacaOrthotropic(TPZInterpolatedElement *cel,REAL hight = 0. );
  /*destrutor*/
  ~TPZPlacaOrthotropic(){}
  /**devolve o tensor de tensões da placa*/
  void Tensor(TPZFMatrix &T, REAL z);
  /**Dados dois vetores n1 e n2 retorna o momento*/
  REAL Moment(TPZVec<REAL> n1, TPZVec<REAL> n2);
  /**Dados dois vetores n1 e n2 retorna a forca*/
  REAL Force(TPZVec<REAL> n1, TPZVec<REAL> n2);
  
  REAL FH(){return fH;}

  void Print();

  TPZInterpolatedElement *ComputEl() {return fIntel;}
  //usa * no nome ComptEl, porque 
};
#endif

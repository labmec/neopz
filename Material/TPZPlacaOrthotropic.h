#ifndef PLACAORTHOTROPIC
#define PLACAORTHOTROPIC
#include "pzvec.h"
template<class>
class TPZVec;
class TPZFMatrix;
class TPZInterpolatedElement;

class TPZPlacaOrthotropic {

 private:

  /**elemento computacional da placa*/
  TPZInterpolatedElement *fIntel;
  /**Espessura da placa*/
  REAL fH;
  
 public:

  /**construtor da placa*/
  TPZPlacaOrthotropic(REAL hight);
  /**devolve o tensor de tensões da placa*/
  void Tensor(TPZFMatrix &T, REAL z);
  /**Dados dois vetores n1 e n2 retorna o momento*/
  REAL Moment(TPZVec<REAL> n1, TPZVec<REAL> n2);
  /**Dados dois vetores n1 e n2 retorna a forca*/
  REAL Force(TPZVec<REAL> n1, TPZVec<REAL> n2);


};
#endif


#include "TPZPlacaOrthotropic.h"

TPZPlacaOrthotropic::TPZPlacaOrthotropic(REAL hight){
	

  fH = hight;
}

REAL TPZPlacaOrthotropic::Moment(TPZVec<REAL> n1, TPZVec<REAL> n2){/**just for fix warning*/ return 0.;}

REAL TPZPlacaOrthotropic::Force(TPZVec<REAL> n1, TPZVec<REAL> n2){/**just for fix warning*/return 0.;}

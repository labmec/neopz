
#include "pzbctension.h"
#include "pzadmchunk.h"


TPZBCTension::TPZBCTension(TPZMaterial *material,int id,int type,
			   TPZFMatrix &val1,TPZFMatrix &val2,TPZFMatrix &esforcos) :
  TPZBndCond(material,id,type,val1,val2) {
  cout << "TPZBCTension::TPZBCTensio FALTA DEFINIR\n";  
}





// $Id: pzbctension.cpp,v 1.4 2003-11-04 16:45:07 phil Exp $

#include "pzbctension.h"
#include "pzadmchunk.h"
#include "pzintel.h"
#include "TPZMulticamadaOrtho.h"
//#include "TPZPlacaOrthotropic.h"

TPZBCTension::TPZBCTension(TPZMaterial *material,int id,int type,
			   TPZFMatrix &val1,TPZFMatrix &val2, TPZMulticamadaOrthotropic *mult, int camada) :
  TPZBndCond(material,id,type,val1,val2) {
  fCamada = camada;


}




// $Id: pzbctension.cc,v 1.5 2003-11-05 00:24:58 phil Exp $

#include "pzbctension.h"
#include "pzadmchunk.h"
#include "pzintel.h"
#include "TPZMulticamadaOrtho.h"
//#include "TPZPlacaOrthotropic.h"

TPZBCTension::TPZBCTension(TPZMaterial *material,int id,int type,
			   TPZFMatrix &val1,TPZFMatrix &val2, TPZMulticamadaOrthotropic *mult, int camada) :
  TPZBndCond(material,id,type,val1,val2) {
  fCamada = camada;
  fMultCam = mult;
  

}




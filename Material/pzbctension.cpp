// $Id: pzbctension.cpp,v 1.6 2003-11-07 00:40:34 phil Exp $

#include "pzbctension.h"
#include "pzadmchunk.h"
#include "pzintel.h"
#include "TPZMulticamadaOrtho.h"
//#include "TPZPlacaOrthotropic.h"

TPZBCTension::TPZBCTension(TPZMaterial *material,int id,int type,
			   TPZFMatrix &val1,TPZFMatrix &val2, REAL sign, TPZMulticamadaOrthotropic *mult, int camada) :
  TPZBndCond(material,id,type,val1,val2) {
  fCamada = camada;
  fMultCam = mult;
  fSign = sign;
  

}




/**
 * \file
 * @brief Contains implementations of the TPZBCTension methods.
 */
// $Id: pzbctension.cpp,v 1.8 2007-05-11 19:15:17 joao Exp $ 

#include "pzbctension.h"
#include "pzadmchunk.h"
#include "pzintel.h"
#include "TPZMulticamadaOrtho.h"
//#include "TPZPlacaOrthotropic.h"

TPZBCTension::TPZBCTension(TPZAutoPointer<TPZMaterial> &material,int id,int type,
						   TPZFMatrix &val1,TPZFMatrix &val2, REAL sign, TPZMulticamadaOrthotropic *mult, int camada) :
TPZBndCond(material,id,type,val1,val2) {
	fCamada = camada;
	fMultCam = mult;
	fSign = sign;
	
	
}

/**
 * @file
 * @brief Contains implementations of the TPZBCTension methods.
 */

#include "pzbctension.h"
#include "pzadmchunk.h"
#include "pzintel.h"
#include "TPZMulticamadaOrtho.h"

TPZBCTension::TPZBCTension(TPZMaterial * &material,int id,int type,
						   TPZFMatrix<REAL> &val1,TPZFMatrix<REAL> &val2, REAL sign, TPZMulticamadaOrthotropic *mult, int camada) :
TPZBndCond(material,id,type,val1,val2) {
	fCamada = camada;
	fMultCam = mult;
	fSign = sign;
	
	
}

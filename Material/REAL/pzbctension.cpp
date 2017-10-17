/**
 * @file
 * @brief Contains implementations of the TPZBCTension methods.
 */

#include "pzbctension.h"
#include "pzadmchunk.h"
#include "pzintel.h"
#include "TPZMulticamadaOrtho.h"

TPZBCTension::TPZBCTension(TPZMaterial * &material,int id,int type,
						   TPZFMatrix<STATE> &val1,TPZFMatrix<STATE> &val2, REAL sign, TPZMulticamadaOrthotropic *mult, int camada) :
TPZRegisterClassId(&TPZBCTension::ClassId),TPZBndCond(material,id,type,val1,val2) {
	fCamada = camada;
	fMultCam = mult;
	fSign = sign;
	
	
}
int TPZBCTension::ClassId() const{
    return Hash("TPZBCTension") ^ TPZBndCond::ClassId() << 1;
}

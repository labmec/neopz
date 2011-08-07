/**
 * @file
 * @brief Contains the implementation of the TPZGaussLobattoRule methods. 
 */
// $Id: tpzgausslobattorule.cpp,v 1.2 2011-03-28 18:19:41 fortiago Exp $
#include "tpzgausslobattorule.h"
#include "tpzintrulelist.h"
#include "pzerror.h"

//***************************************
//***************************************
TPZGaussLobattoRule::TPZGaussLobattoRule(int precision){
	
	if(precision < 1 && precision > NUMINT_RULES){
		PZError << "TPZGaussLobattoRule creation precision = " << precision << " not allowable\n";
		fNumInt = 0;
		fLocation = NULL;
		fWeight = NULL;
		return;
	}
	
	int intpoints[] = {2,2,2,3,4,5,6,7,8,9,10};
	
	int numpoints = (int)((precision+3.)/2.+0.5);
	
	if(numpoints > 10) DebugStop();
	fNumInt = (short) intpoints[numpoints];
	fLocation = new REAL[fNumInt];
	fWeight = new REAL[fNumInt];
	
	if(fLocation == NULL || fWeight == NULL){
		fNumInt = 0;
		return;
	}
	
	switch(fNumInt){
			
		case 1:
		case 2:
			fLocation[0] = -1.;					fWeight[0] = 1.;
			fLocation[0] =  1.;					fWeight[0] = 1.;
			break;
			
		case 3:
			fLocation[0] = -1.;					fWeight[0] = 1./3.;
			fLocation[1] =  0.;					fWeight[1] = 4./3.;
			fLocation[2] =  1.;					fWeight[2] = 1./3.;
			break;
			
		case 4:
			fLocation[0] = -1.;					fWeight[0] = 1./6.;
			fLocation[1] = -0.447213595499958;	fWeight[1] = 5./6.;
			fLocation[2] = -fLocation[1];			fWeight[2] = 5./6.;
			fLocation[3] =  1.;					fWeight[3] = 1./6.;
			break;
			
		case 5:
			fLocation[0] = -1.;					fWeight[0] =  1./10.;
			fLocation[1] = -0.654653670707977;	fWeight[1] = 49./90.;
			fLocation[2] =  0.;					fWeight[2] = 64./90.;
			fLocation[3] = -fLocation[1];			fWeight[3] = 49./90.;
			fLocation[4] =  1.;					fWeight[4] =  1./10.;
			break;
			
		case 6:
			fLocation[0] = -1.;					fWeight[0] =  1./15.;
			fLocation[1] = -0.765055323929465;	fWeight[1] =  0.378474956297847;
			fLocation[2] = -0.285231516480645;	fWeight[2] =  0.554858377035486;
			fLocation[3] = -fLocation[2];			fWeight[3] =  fWeight[2];
			fLocation[4] = -fLocation[1];			fWeight[4] =  fWeight[1];
			fLocation[5] =  1.;					fWeight[5] =  1./15.;
			break;
			
		case 7:
			fLocation[0] = -1.;					fWeight[0] =  1./21.;
			fLocation[1] = -0.830223896278567;	fWeight[1] =  0.276826047361566;
			fLocation[2] = -0.468848793470714;	fWeight[2] =  0.431745381209863;
			fLocation[3] =  0.;					fWeight[3] =  0.487619047619048;
			fLocation[4] = -fLocation[2];			fWeight[4] =  fWeight[2];
			fLocation[5] = -fLocation[1];			fWeight[5] =  fWeight[1];
			fLocation[6] =  1.;					fWeight[6] =  1./21.;
			break;
			
		case 8:
			fLocation[0] = -1.;					fWeight[0] =  1./28.;
			fLocation[1] = -0.871740148509607;	fWeight[1] =  0.210704227143506;
			fLocation[2] = -0.591700181433142;	fWeight[2] =  0.341122692483504;
			fLocation[3] = -0.209299217902479;	fWeight[3] =  0.412458794658704;
			fLocation[4] = -fLocation[3];			fWeight[4] =  fWeight[3];
			fLocation[5] = -fLocation[2];			fWeight[5] =  fWeight[2];
			fLocation[6] = -fLocation[1];			fWeight[6] =  fWeight[1];
			fLocation[7] =  1.;					fWeight[7] =  1./28.;
			break;
			
		case 9:
			fLocation[0] = -1.;					fWeight[0] =  1./36.;
			fLocation[1] = -0.899757995411460;	fWeight[1] =  0.165495361560805;
			fLocation[2] = -0.677186279510738;	fWeight[2] =  0.274538712500162;
			fLocation[3] = -0.363117463826178;	fWeight[3] =  0.346428510973046;
			fLocation[4] =  0.;					fWeight[4] =  0.371519274376417;
			fLocation[5] = -fLocation[3];			fWeight[5] =  fWeight[3];
			fLocation[6] = -fLocation[2];			fWeight[6] =  fWeight[2];
			fLocation[7] = -fLocation[1];			fWeight[7] =  fWeight[1];
			fLocation[8] =  1.;					fWeight[8] =  1./36.;
			break;
			
		case 10:
			fLocation[0] = -1.;					fWeight[0] =  1./45.;
			fLocation[1] = -0.919533908166459;	fWeight[1] =  0.133305990851070;
			fLocation[2] = -0.738773865105505;	fWeight[2] =  0.224889342063126;
			fLocation[3] = -0.477924949810444;	fWeight[3] =  0.292042683679684;
			fLocation[4] = -0.165278957666387;	fWeight[4] =  0.327539761183897;
			fLocation[5] =  -fLocation[4];		fWeight[5] =  fWeight[4];
			fLocation[6] =  -fLocation[3];		fWeight[6] =  fWeight[3];
			fLocation[7] =  -fLocation[2];		fWeight[7] =  fWeight[2];
			fLocation[8] =  -fLocation[1];		fWeight[8] =  fWeight[1];
			fLocation[9] =  1.;					fWeight[9] =  1./45.;
			break;
			
		default:
			PZError << "TPZGaussLobattoRule creation : invalid number of integration points "
			" specified\n";
			
			fNumInt = 0;
	}
	
	
}

TPZGaussLobattoRule::~TPZGaussLobattoRule(){
	
	if (fLocation) delete []fLocation;
	if (fWeight)   delete []fWeight;
	
}

REAL TPZGaussLobattoRule::Loc(int i) {
	
	if (fLocation && i>=0 && i<fNumInt)
		return fLocation[i];
	else {
		PZError << "ERROR(TPZGaussLobattoRule::loc) Out of bounds!!\n";
		return 0.0;
	}
}

REAL TPZGaussLobattoRule::W(int i) {
	
	if (fWeight && i>=0 && i<fNumInt)
		return fWeight[i];
	else {
		PZError << "ERROR(TPZGaussLobattoRule::w) Out of bounds!!\n";
		return 0.0;
	}
}


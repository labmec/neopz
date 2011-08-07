/**
 * @file
 * @brief Contains the implementation of the TPZIntRuleT3D methods. 
 */
//
// C++ Implementation: tpzintrulet3d
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzintrulet3d.h"
#include "pzerror.h"
#include "pzvec.h"

REAL TPZIntRuleT3D::W(int i) const {
	
	if (fWeight && i>=0 && i<fNumInt)
		return fWeight[i];
	else {
		PZError << "ERROR(TPZIntRuleT3D::w) Out of bounds!!\n";
		return 0.0;
	}
}

//------------------------------------------------------------------------------
TPZIntRuleT3D::~TPZIntRuleT3D(){
	
	if (fLocationKsi) delete []fLocationKsi;
	if (fLocationEta) delete []fLocationEta;
	if (fLocationZeta) delete []fLocationZeta;
	if (fWeight)   delete []fWeight;
	
}
//------------------------------------------------------------------------------
void TPZIntRuleT3D::Loc(int i, TPZVec<REAL> &Points) const {
	
	if (fLocationKsi && fLocationEta && fLocationZeta && i>=0 && i<fNumInt){
		Points[0] = fLocationKsi[i];
		Points[1] = fLocationEta[i];
		Points[2] = fLocationZeta[i];
		return;
	}
	else {
		PZError << "ERROR(TPZIntRuleT::loc) Out of bounds!!\n";
	}
}
//------------------------------------------------------------------------------
TPZIntRuleT3D::TPZIntRuleT3D(int precision){
	
	if(precision < 1 && precision > NUMINT_RULEST3D){
		PZError << "TPZIntRule creation precision = " << precision << " not available\n";
		precision = NUMINT_RULEST3D;
	}
	
	fNumInt = precision;
	
	if (precision <= 1)  fNumInt =   1;//integra constantes e lineares
	else if (precision == 2)  fNumInt =   4;//integra quadraticas
	else if (precision == 3)  fNumInt =  18;//integra cubicas
	else if (precision == 4)  fNumInt =  24;//integra de grau 4
	else if (precision <= 6)  fNumInt =  65;//integra de grau 6
	else if (precision == 7)  fNumInt = 125;//integra de grau 7
	else if (precision == 8)  fNumInt = 150;//integra de grau 8
	
	fLocationKsi  = new REAL[fNumInt];
	fLocationEta  = new REAL[fNumInt];
	fLocationZeta = new REAL[fNumInt];
	fWeight = new REAL[fNumInt];
	
	if(fLocationKsi == NULL || fLocationEta == NULL || fLocationZeta == NULL || fWeight == NULL){
		fNumInt = 0;
		return;
	}
	REAL p,q,w;
	
	switch(fNumInt){
			
		case 1://integra lineares
			fLocationKsi[0]  = 0.25;//1./4.
			fLocationEta[0]  = 0.25;//1./4.
			fLocationZeta[0] = 0.25;//1./4.
			fWeight[0]       = 1./6.;
			break;
			
		case 4://integra quadraticas
			p = 1./4. - sqrt(5.)/20.;
			q = 1./4. + 3 * sqrt(5.)/20.;
			w = 1./24.;
			
			fLocationKsi[0]  = p;
			fLocationEta[0]  = p;
			fLocationZeta[0] = p;
			fWeight[0]  = w;
			
			fLocationKsi[1]  = p;
			fLocationEta[1]  = q;
			fLocationZeta[1] = p;
			fWeight[1]       = w;
			
			fLocationKsi[2]  = q;
			fLocationEta[2]  = p;
			fLocationZeta[2] = p;
			fWeight[2]       = w;
			
			fLocationKsi[3]  = p;
			fLocationEta[3]  = p;
			fLocationZeta[3] = q;
			fWeight[3]       = w;
			break;
			
			/*  case 5://integra cubicas com peso negativo
			 p = 0.25;//1./4.
			 q = 1./6.;
			 r = 0.5;//1./2.
			 fLocationKsi[0]  = p;
			 fLocationEta[0]  = p;
			 fLocationZeta[0] = p;
			 fWeight     [0]  = -2./15.;
			 
			 fLocationKsi[1]  = q;
			 fLocationEta[1]  = q;
			 fLocationZeta[1] = q;
			 fWeight     [1]  = 0.075;//3./40.
			 
			 fLocationKsi[2]  = q;
			 fLocationEta[2]  = q;
			 fLocationZeta[2] = r;
			 fWeight     [2]  = 0.075;
			 
			 fLocationKsi[3]  = q;
			 fLocationEta[3]  = r;
			 fLocationZeta[3] = q;
			 fWeight     [3]  = 0.075;
			 
			 fLocationKsi[4]  = r;
			 fLocationEta[4]  = q;
			 fLocationZeta[4] = q;
			 fWeight     [4]  = 0.075;
			 break;               */
			
			
        case 18://integra cubicas com todos os pesos positivos
			
			fLocationKsi[0]  = 0.092060081835903;
			fLocationEta[0]  = 0.010320791771678;
			fLocationZeta[0] = 0.887298334620742;
			fWeight[0]       = 0.000193968091080;
			
			fLocationKsi[1]  = 0.408423786490229;
			fLocationEta[1]  = 0.045788106754886;
			fLocationZeta[1] = 0.500000000000000;
			fWeight[1]       = 0.006108430203073;
			
			fLocationKsi[2]  = 0.724787491144556;
			fLocationEta[2]  = 0.081255421738093;
			fLocationZeta[2] = 0.112701665379258;
			fWeight[2]       = 0.012022892315067;
			
			fLocationKsi[3]  = 0.010320791771678;
			fLocationEta[3]  = 0.092060081835903;
			fLocationZeta[3] = 0.887298334620742;
			fWeight[3]       = 0.000193968091080;
			
			fLocationKsi[4]  = 0.045788106754886;
			fLocationEta[4]  = 0.408423786490229;
			fLocationZeta[4] = 0.500000000000000;
			fWeight[4]       = 0.006108430203073;
			
			fLocationKsi[5]  = 0.081255421738093;
			fLocationEta[5]  = 0.724787491144556;
			fLocationZeta[5] = 0.112701665379258;
			fWeight[5]       = 0.012022892315067;
			
			fLocationKsi[6]  = 0.010320791771678;
			fLocationEta[6]  = 0.010320791771678;
			fLocationZeta[6] = 0.887298334620742;
			fWeight[6]       = 0.000193968091080;
			
			fLocationKsi[7]  = 0.045788106754886;
			fLocationEta[7]  = 0.045788106754886;
			fLocationZeta[7] = 0.500000000000000;
			fWeight[7]       = 0.006108430203073;
			
			fLocationKsi[8]  = 0.081255421738093;
			fLocationEta[8]  = 0.081255421738093;
			fLocationZeta[8] = 0.112701665379258;
			fWeight[8]       = 0.012022892315067;
			
			fLocationKsi[9]  = 0.012183390180066;
			fLocationEta[9]  = 0.050259137599596;
			fLocationZeta[9] = 0.887298334620742;
			fWeight[9]       = 0.000394071972775;
			
			fLocationKsi[10]  = 0.054051509084035;
			fLocationEta[10]  = 0.222974245457983;
			fLocationZeta[10] = 0.500000000000000;
			fWeight[10]       = 0.012410088315445;
			
			fLocationKsi[11]  = 0.095919627988004;
			fLocationEta[11]  = 0.395689353316369;
			fLocationZeta[11] = 0.112701665379258;
			fWeight[11]       = 0.024426104658116;
			
			fLocationKsi[12]  = 0.050259137599596;
			fLocationEta[12]  = 0.012183390180066;
			fLocationZeta[12] = 0.887298334620742;
			fWeight[12]       = 0.000394071972775;
			
			fLocationKsi[13]  = 0.222974245457983;
			fLocationEta[13]  = 0.054051509084035;
			fLocationZeta[13] = 0.500000000000000;
			fWeight[13]       = 0.012410088315445;
			
			fLocationKsi[14]  = 0.395689353316369;
			fLocationEta[14]  = 0.095919627988004;
			fLocationZeta[14] = 0.112701665379258;
			fWeight[14]       = 0.024426104658116;
			
			fLocationKsi[15]  = 0.050259137599596;
			fLocationEta[15]  = 0.050259137599596;
			fLocationZeta[15] = 0.887298334620742;
			fWeight[15]       = 0.000394071972775;
			
			fLocationKsi[16]  = 0.222974245457983;
			fLocationEta[16]  = 0.222974245457983;
			fLocationZeta[16] = 0.500000000000000;
			fWeight[16]       = 0.012410088315445;
			
			fLocationKsi[17]  = 0.395689353316369;
			fLocationEta[17]  = 0.395689353316369;
			fLocationZeta[17] = 0.112701665379258;
			fWeight[17]       = 0.024426104658116;
			
			
			break;
			
		case 24://integra at� grau 4
			
            fLocationKsi[0]  = 0.056715233422327;
            fLocationEta[0]  = 0.006358305388836;
            fLocationZeta[0] = 0.930568155800000;
            fWeight[0]       = 0.000046095400013;
			
            fLocationKsi[1]  = 0.269567441328218;
            fLocationEta[1]  = 0.030221018435891;
            fLocationZeta[1] = 0.669990521800000;
            fWeight[1]       = 0.001952267068526;
			
            fLocationKsi[2]  = 0.547280131652241;
            fLocationEta[2]  = 0.061355195073880;
            fLocationZeta[2] = 0.330009478200000;
            fWeight[2]       = 0.008046809490043;
			
            fLocationKsi[3]  = 0.760132339558132;
            fLocationEta[3]  = 0.085217908120935;
            fLocationZeta[3] = 0.069431844200000;
            fWeight[3]       = 0.008280118650457;
			
            fLocationKsi[4]  = 0.006358305388836;
            fLocationEta[4]  = 0.056715233422327;
            fLocationZeta[4] = 0.930568155800000;
            fWeight[4]       = 0.000046095400013;
			
            fLocationKsi[5]  = 0.030221018435891;
            fLocationEta[5]  = 0.269567441328218;
            fLocationZeta[5] = 0.669990521800000;
            fWeight[5]       = 0.001952267068526;
			
            fLocationKsi[6]  = 0.061355195073880;
            fLocationEta[6]  = 0.547280131652241;
            fLocationZeta[6] = 0.330009478200000;
            fWeight[6]       = 0.008046809490043;
			
            fLocationKsi[7]  = 0.085217908120935;
            fLocationEta[7]  = 0.760132339558132;
            fLocationZeta[7] = 0.069431844200000;
            fWeight[7]       = 0.008280118650457;
			
            fLocationKsi[8]  = 0.006358305388836;
            fLocationEta[8]  = 0.006358305388836;
            fLocationZeta[8] = 0.930568155800000;
            fWeight[8]       = 0.000046095400013;
			
            fLocationKsi[9]  = 0.030221018435891;
            fLocationEta[9]  = 0.030221018435891;
            fLocationZeta[9] = 0.669990521800000;
            fWeight[9]       = 0.001952267068526;
			
            fLocationKsi[10]  = 0.061355195073880;
            fLocationEta[10]  = 0.061355195073880;
            fLocationZeta[10] = 0.330009478200000;
            fWeight[10]       = 0.008046809490043;
			
            fLocationKsi[11]  = 0.085217908120935;
            fLocationEta[11]  = 0.085217908120935;
            fLocationZeta[11] = 0.069431844200000;
            fWeight[11]       = 0.008280118650457;
			
            fLocationKsi[12]  = 0.007505791914995;
            fLocationEta[12]  = 0.030963026142502;
            fLocationZeta[12] = 0.930568155800000;
            fWeight[12]       = 0.000093648935337;
			
            fLocationKsi[13]  = 0.035675020617490;
            fLocationEta[13]  = 0.147167228791255;
            fLocationZeta[13] = 0.669990521800000;
            fWeight[13]       = 0.003966290180995;
			
            fLocationKsi[14]  = 0.072427997550580;
            fLocationEta[14]  = 0.298781262124710;
            fLocationZeta[14] = 0.330009478200000;
            fWeight[14]       = 0.016348163621277;
			
            fLocationKsi[15]  = 0.100597226253075;
            fLocationEta[15]  = 0.414985464773463;
            fLocationZeta[15] = 0.069431844200000;
            fWeight[15]       = 0.016822162208359;
			
            fLocationKsi[16]  = 0.030963026142502;
            fLocationEta[16]  = 0.007505791914995;
            fLocationZeta[16] = 0.930568155800000;
            fWeight[16]       = 0.000093648935337;
			
            fLocationKsi[17]  = 0.147167228791255;
            fLocationEta[17]  = 0.035675020617490;
            fLocationZeta[17] = 0.669990521800000;
            fWeight[17]       = 0.003966290180995;
			
            fLocationKsi[18]  = 0.298781262124710;
            fLocationEta[18]  = 0.072427997550580;
            fLocationZeta[18] = 0.330009478200000;
            fWeight[18]       = 0.016348163621277;
			
            fLocationKsi[19]  = 0.414985464773463;
            fLocationEta[19]  = 0.100597226253075;
            fLocationZeta[19] = 0.069431844200000;
            fWeight[19]       = 0.016822162208359;
			
            fLocationKsi[20]  = 0.030963026142502;
            fLocationEta[20]  = 0.030963026142502;
            fLocationZeta[20] = 0.930568155800000;
            fWeight[20]       = 0.000093648935337;
			
            fLocationKsi[21]  = 0.147167228791255;
            fLocationEta[21]  = 0.147167228791255;
            fLocationZeta[21] = 0.669990521800000;
            fWeight[21]       = 0.003966290180995;
			
            fLocationKsi[22]  = 0.298781262124710;
            fLocationEta[22]  = 0.298781262124710;
            fLocationZeta[22] = 0.330009478200000;
            fWeight[22]       = 0.016348163621277;
			
            fLocationKsi[23]  = 0.414985464773463;
            fLocationEta[23]  = 0.414985464773463;
            fLocationZeta[23] = 0.069431844200000;
            fWeight[23]       = 0.016822162208359;
			
			
            break;
			
		case 65://integra ate grau 6
			
			fLocationKsi[0]  = 0.015636692350000;
			fLocationEta[0]  = 0.015636692350000;
			fLocationZeta[0] = 0.953089922950000;
			fWeight[0]       = -0.000019495360419;
			
			fLocationKsi[1]  = 0.076921781650000;
			fLocationEta[1]  = 0.076921781650000;
			fLocationZeta[1] = 0.769234655050000;
			fWeight[1]       = -0.000953069387246;
			
			fLocationKsi[2]  = 0.166666666666666;
			fLocationEta[2]  = 0.166666666666666;
			fLocationZeta[2] = 0.500000000000000;
			fWeight[2]       = -0.005318046024686;
			
			fLocationKsi[3]  = 0.256411551683333;
			fLocationEta[3]  = 0.256411551683333;
			fLocationZeta[3] = 0.230765344950000;
			fWeight[3]       = -0.010590123488919;
			
			fLocationKsi[4]  = 0.317696640983333;
			fLocationEta[4]  = 0.317696640983333;
			fLocationZeta[4] = 0.046910077050000;
			fWeight[4]       = -0.008047606477835;
			
			fLocationKsi[5]  = 0.022484378393151;
			fLocationEta[5]  = 0.012212849328424;
			fLocationZeta[5] = 0.953089922950000;
			fWeight[5]       = 0.000022890163273;
			
			fLocationKsi[6]  = 0.110607691612859;
			fLocationEta[6]  = 0.060078826668570;
			fLocationZeta[6] = 0.769234655050000;
			fWeight[6]       = 0.001119031062595;
			
			fLocationKsi[7]  = 0.239654033920961;
			fLocationEta[7]  = 0.130172983039519;
			fLocationZeta[7] = 0.500000000000000;
			fWeight[7]       = 0.006244098041094;
			
			fLocationKsi[8]  = 0.368700376229064;
			fLocationEta[8]  = 0.200267139410468;
			fLocationZeta[8] = 0.230765344950000;
			fWeight[8]       = 0.012434222837702;
			
			fLocationKsi[9]  = 0.456823689448772;
			fLocationEta[9]  = 0.248133116750614;
			fLocationZeta[9] = 0.046910077050000;
			fWeight[9]       = 0.009448967461072;
			
			fLocationKsi[10]  = 0.012212849328424;
			fLocationEta[10]  = 0.022484378393151;
			fLocationZeta[10] = 0.953089922950000;
			fWeight[10]       = 0.000022890163273;
			
			fLocationKsi[11]  = 0.060078826668570;
			fLocationEta[11]  = 0.110607691612859;
			fLocationZeta[11] = 0.769234655050000;
			fWeight[11]       = 0.001119031062595;
			
			fLocationKsi[12]  = 0.130172983039519;
			fLocationEta[12]  = 0.239654033920961;
			fLocationZeta[12] = 0.500000000000000;
			fWeight[12]       = 0.006244098041094;
			
			fLocationKsi[13]  = 0.200267139410468;
			fLocationEta[13]  = 0.368700376229064;
			fLocationZeta[13] = 0.230765344950000;
			fWeight[13]       = 0.012434222837702;
			
			fLocationKsi[14]  = 0.248133116750614;
			fLocationEta[14]  = 0.456823689448772;
			fLocationZeta[14] = 0.046910077050000;
			fWeight[14]       = 0.009448967461072;
			
			fLocationKsi[15]  = 0.012212849328424;
			fLocationEta[15]  = 0.012212849328424;
			fLocationZeta[15] = 0.953089922950000;
			fWeight[15]       = 0.000022890163273;
			
			fLocationKsi[16]  = 0.060078826668570;
			fLocationEta[16]  = 0.060078826668570;
			fLocationZeta[16] = 0.769234655050000;
			fWeight[16]       = 0.001119031062595;
			
			fLocationKsi[17]  = 0.130172983039519;
			fLocationEta[17]  = 0.130172983039519;
			fLocationZeta[17] = 0.500000000000000;
			fWeight[17]       = 0.006244098041094;
			
			fLocationKsi[18]  = 0.200267139410468;
			fLocationEta[18]  = 0.200267139410468;
			fLocationZeta[18] = 0.230765344950000;
			fWeight[18]       = 0.012434222837702;
			
			fLocationKsi[19]  = 0.248133116750614;
			fLocationEta[19]  = 0.248133116750614;
			fLocationZeta[19] = 0.046910077050000;
			fWeight[19]       = 0.009448967461072;
			
			fLocationKsi[20]  = 0.040799560759165;
			fLocationEta[20]  = 0.003055258145417;
			fLocationZeta[20] = 0.953089922950000;
			fWeight[20]       = 0.000006953421651;
			
			fLocationKsi[21]  = 0.200705803624282;
			fLocationEta[21]  = 0.015029770662859;
			fLocationZeta[21] = 0.769234655050000;
			fWeight[21]       = 0.000339931818126;
			
			fLocationKsi[22]  = 0.434869897097784;
			fLocationEta[22]  = 0.032565051451108;
			fLocationZeta[22] = 0.500000000000000;
			fWeight[22]       = 0.001896790599129;
			
			fLocationKsi[23]  = 0.669033990571286;
			fLocationEta[23]  = 0.050100332239357;
			fLocationZeta[23] = 0.230765344950000;
			fWeight[23]       = 0.003777185564802;
			
			fLocationKsi[24]  = 0.828940233436403;
			fLocationEta[24]  = 0.062074844756799;
			fLocationZeta[24] = 0.046910077050000;
			fWeight[24]       = 0.002870344529135;
			
			fLocationKsi[25]  = 0.003055258145417;
			fLocationEta[25]  = 0.040799560759165;
			fLocationZeta[25] = 0.953089922950000;
			fWeight[25]       = 0.000006953421651;
			
			fLocationKsi[26]  = 0.015029770662859;
			fLocationEta[26]  = 0.200705803624282;
			fLocationZeta[26] = 0.769234655050000;
			fWeight[26]       = 0.000339931818126;
			
			fLocationKsi[27]  = 0.032565051451108;
			fLocationEta[27]  = 0.434869897097784;
			fLocationZeta[27] = 0.500000000000000;
			fWeight[27]       = 0.001896790599129;
			
			fLocationKsi[28]  = 0.050100332239357;
			fLocationEta[28]  = 0.669033990571286;
			fLocationZeta[28] = 0.230765344950000;
			fWeight[28]       = 0.003777185564802;
			
			fLocationKsi[29]  = 0.062074844756799;
			fLocationEta[29]  = 0.828940233436403;
			fLocationZeta[29] = 0.046910077050000;
			fWeight[29]       = 0.002870344529135;
			
			fLocationKsi[30]  = 0.003055258145417;
			fLocationEta[30]  = 0.003055258145417;
			fLocationZeta[30] = 0.953089922950000;
			fWeight[30]       = 0.000006953421651;
			
			fLocationKsi[31]  = 0.015029770662859;
			fLocationEta[31]  = 0.015029770662859;
			fLocationZeta[31] = 0.769234655050000;
			fWeight[31]       = 0.000339931818126;
			
			fLocationKsi[32]  = 0.032565051451108;
			fLocationEta[32]  = 0.032565051451108;
			fLocationZeta[32] = 0.500000000000000;
			fWeight[32]       = 0.001896790599129;
			
			fLocationKsi[33]  = 0.050100332239357;
			fLocationEta[33]  = 0.050100332239357;
			fLocationZeta[33] = 0.230765344950000;
			fWeight[33]       = 0.003777185564802;
			
			fLocationKsi[34]  = 0.062074844756799;
			fLocationEta[34]  = 0.062074844756799;
			fLocationZeta[34] = 0.046910077050000;
			fWeight[34]       = 0.002870344529135;
			
			fLocationKsi[35]  = 0.029949466077934;
			fLocationEta[35]  = 0.014676544523875;
			fLocationZeta[35] = 0.953089922950000;
			fWeight[35]       = 0.000010051214246;
			
			fLocationKsi[36]  = 0.147330793406635;
			fLocationEta[36]  = 0.072198514108518;
			fLocationZeta[36] = 0.769234655050000;
			fWeight[36]       = 0.000491373557463;
			
			fLocationKsi[37]  = 0.319222094284905;
			fLocationEta[37]  = 0.156432748002438;
			fLocationZeta[37] = 0.500000000000000;
			fWeight[37]       = 0.002741822609003;
			
			fLocationKsi[38]  = 0.491113395163174;
			fLocationEta[38]  = 0.240666981896357;
			fLocationZeta[38] = 0.230765344950000;
			fWeight[38]       = 0.005459945227865;
			
			fLocationKsi[39]  = 0.608494722491875;
			fLocationEta[39]  = 0.298188951481000;
			fLocationZeta[39] = 0.046910077050000;
			fWeight[39]       = 0.004149100870294;
			
			fLocationKsi[40]  = 0.002284066448190;
			fLocationEta[40]  = 0.014676544523875;
			fLocationZeta[40] = 0.953089922950000;
			fWeight[40]       = 0.000010051214246;
			
			fLocationKsi[41]  = 0.011236037434847;
			fLocationEta[41]  = 0.072198514108518;
			fLocationZeta[41] = 0.769234655050000;
			fWeight[41]       = 0.000491373557463;
			
			fLocationKsi[42]  = 0.024345157712658;
			fLocationEta[42]  = 0.156432748002438;
			fLocationZeta[42] = 0.500000000000000;
			fWeight[42]       = 0.002741822609003;
			
			fLocationKsi[43]  = 0.037454277990469;
			fLocationEta[43]  = 0.240666981896357;
			fLocationZeta[43] = 0.230765344950000;
			fWeight[43]       = 0.005459945227865;
			
			fLocationKsi[44]  = 0.046406248977126;
			fLocationEta[44]  = 0.298188951481000;
			fLocationZeta[44] = 0.046910077050000;
			fWeight[44]       = 0.004149100870294;
			
			fLocationKsi[45]  = 0.014676544523875;
			fLocationEta[45]  = 0.029949466077934;
			fLocationZeta[45] = 0.953089922950000;
			fWeight[45]       = 0.000010051214246;
			
			fLocationKsi[46]  = 0.072198514108518;
			fLocationEta[46]  = 0.147330793406635;
			fLocationZeta[46] = 0.769234655050000;
			fWeight[46]       = 0.000491373557463;
			
			fLocationKsi[47]  = 0.156432748002438;
			fLocationEta[47]  = 0.319222094284905;
			fLocationZeta[47] = 0.500000000000000;
			fWeight[47]       = 0.002741822609003;
			
			fLocationKsi[48]  = 0.240666981896357;
			fLocationEta[48]  = 0.491113395163174;
			fLocationZeta[48] = 0.230765344950000;
			fWeight[48]       = 0.005459945227865;
			
			fLocationKsi[49]  = 0.298188951481000;
			fLocationEta[49]  = 0.608494722491875;
			fLocationZeta[49] = 0.046910077050000;
			fWeight[49]       = 0.004149100870294;
			
			fLocationKsi[50]  = 0.014676544523875;
			fLocationEta[50]  = 0.002284066448190;
			fLocationZeta[50] = 0.953089922950000;
			fWeight[50]       = 0.000010051214246;
			
			fLocationKsi[51]  = 0.072198514108518;
			fLocationEta[51]  = 0.011236037434847;
			fLocationZeta[51] = 0.769234655050000;
			fWeight[51]       = 0.000491373557463;
			
			fLocationKsi[52]  = 0.156432748002438;
			fLocationEta[52]  = 0.024345157712658;
			fLocationZeta[52] = 0.500000000000000;
			fWeight[52]       = 0.002741822609003;
			
			fLocationKsi[53]  = 0.240666981896357;
			fLocationEta[53]  = 0.037454277990469;
			fLocationZeta[53] = 0.230765344950000;
			fWeight[53]       = 0.005459945227865;
			
			fLocationKsi[54]  = 0.298188951481000;
			fLocationEta[54]  = 0.046406248977126;
			fLocationZeta[54] = 0.046910077050000;
			fWeight[54]       = 0.004149100870294;
			
			fLocationKsi[55]  = 0.029949466077934;
			fLocationEta[55]  = 0.002284066448190;
			fLocationZeta[55] = 0.953089922950000;
			fWeight[55]       = 0.000010051214246;
			
			fLocationKsi[56]  = 0.147330793406635;
			fLocationEta[56]  = 0.011236037434847;
			fLocationZeta[56] = 0.769234655050000;
			fWeight[56]       = 0.000491373557463;
			
			fLocationKsi[57]  = 0.319222094284905;
			fLocationEta[57]  = 0.024345157712658;
			fLocationZeta[57] = 0.500000000000000;
			fWeight[57]       = 0.002741822609003;
			
			fLocationKsi[58]  = 0.491113395163174;
			fLocationEta[58]  = 0.037454277990469;
			fLocationZeta[58] = 0.230765344950000;
			fWeight[58]       = 0.005459945227865;
			
			fLocationKsi[59]  = 0.608494722491875;
			fLocationEta[59]  = 0.046406248977126;
			fLocationZeta[59] = 0.046910077050000;
			fWeight[59]       = 0.004149100870294;
			
			fLocationKsi[60]  = 0.002284066448190;
			fLocationEta[60]  = 0.029949466077934;
			fLocationZeta[60] = 0.953089922950000;
			fWeight[60]       = 0.000010051214246;
			
			fLocationKsi[61]  = 0.011236037434847;
			fLocationEta[61]  = 0.147330793406635;
			fLocationZeta[61] = 0.769234655050000;
			fWeight[61]       = 0.000491373557463;
			
			fLocationKsi[62]  = 0.024345157712658;
			fLocationEta[62]  = 0.319222094284905;
			fLocationZeta[62] = 0.500000000000000;
			fWeight[62]       = 0.002741822609003;
			
			fLocationKsi[63]  = 0.037454277990469;
			fLocationEta[63]  = 0.491113395163174;
			fLocationZeta[63] = 0.230765344950000;
			fWeight[63]       = 0.005459945227865;
			
			fLocationKsi[64]  = 0.046406248977126;
			fLocationEta[64]  = 0.608494722491875;
			fLocationZeta[64] = 0.046910077050000;
			fWeight[64]       = 0.004149100870294;
			
			break;
			
		case 125://integra x7 , x2y2z3 etc.
			
			fLocationKsi[0]  = 0.000103228220029;
			fLocationEta[0]  = 0.044709521721163;
			fLocationZeta[0] = 0.953089922950000;
			fWeight[0]       = 0.000000171613515;
			
			fLocationKsi[1]  = 0.000507811909541;
			fLocationEta[1]  = 0.219940124837926;
			fLocationZeta[1] = 0.769234655050000;
			fWeight[1]       = 0.000008389667268;
			
			fLocationKsi[2]  = 0.001100277664419;
			fLocationEta[2]  = 0.476544961475000;
			fLocationZeta[2] = 0.500000000000000;
			fWeight[2]       = 0.000046813628953;
			
			fLocationKsi[3]  = 0.001692743419296;
			fLocationEta[3]  = 0.733149798112074;
			fLocationZeta[3] = 0.230765344950000;
			fWeight[3]       = 0.000093222606438;
			
			fLocationKsi[4]  = 0.002097327108808;
			fLocationEta[4]  = 0.908380401228837;
			fLocationZeta[4] = 0.046910077050000;
			fWeight[4]       = 0.000070841369531;
			
			fLocationKsi[5]  = 0.000507811909541;
			fLocationEta[5]  = 0.044709521721163;
			fLocationZeta[5] = 0.953089922950000;
			fWeight[5]       = 0.000000346685638;
			
			fLocationKsi[6]  = 0.002498085653322;
			fLocationEta[6]  = 0.219940124837926;
			fLocationZeta[6] = 0.769234655050000;
			fWeight[6]       = 0.000016948415497;
			
			fLocationKsi[7]  = 0.005412610056037;
			fLocationEta[7]  = 0.476544961475000;
			fLocationZeta[7] = 0.500000000000000;
			fWeight[7]       = 0.000094570715276;
			
			fLocationKsi[8]  = 0.008327134458752;
			fLocationEta[8]  = 0.733149798112074;
			fLocationZeta[8] = 0.230765344950000;
			fWeight[8]       = 0.000188323972481;
			
			fLocationKsi[9]  = 0.010317408202533;
			fLocationEta[9]  = 0.908380401228837;
			fLocationZeta[9] = 0.046910077050000;
			fWeight[9]       = 0.000143110438935;
			
			fLocationKsi[10]  = 0.001100277664418;
			fLocationEta[10]  = 0.044709521721163;
			fLocationZeta[10] = 0.953089922950000;
			fWeight[10]       = 0.000000412063923;
			
			fLocationKsi[11]  = 0.005412610056037;
			fLocationEta[11]  = 0.219940124837926;
			fLocationZeta[11] = 0.769234655050000;
			fWeight[11]       = 0.000020144562696;
			
			fLocationKsi[12]  = 0.011727519262500;
			fLocationEta[12]  = 0.476544961475000;
			fLocationZeta[12] = 0.500000000000000;
			fWeight[12]       = 0.000112404944487;
			
			fLocationKsi[13]  = 0.018042428468963;
			fLocationEta[13]  = 0.733149798112074;
			fLocationZeta[13] = 0.230765344950000;
			fWeight[13]       = 0.000223838273938;
			
			fLocationKsi[14]  = 0.022354760860582;
			fLocationEta[14]  = 0.908380401228837;
			fLocationZeta[14] = 0.046910077050000;
			fWeight[14]       = 0.000170098332207;
			
			fLocationKsi[15]  = 0.001692743419296;
			fLocationEta[15]  = 0.044709521721163;
			fLocationZeta[15] = 0.953089922950000;
			fWeight[15]       = 0.000000346685638;
			
			fLocationKsi[16]  = 0.008327134458752;
			fLocationEta[16]  = 0.219940124837926;
			fLocationZeta[16] = 0.769234655050000;
			fWeight[16]       = 0.000016948415497;
			
			fLocationKsi[17]  = 0.018042428468963;
			fLocationEta[17]  = 0.476544961475000;
			fLocationZeta[17] = 0.500000000000000;
			fWeight[17]       = 0.000094570715276;
			
			fLocationKsi[18]  = 0.027757722479174;
			fLocationEta[18]  = 0.733149798112074;
			fLocationZeta[18] = 0.230765344950000;
			fWeight[18]       = 0.000188323972481;
			
			fLocationKsi[19]  = 0.034392113518630;
			fLocationEta[19]  = 0.908380401228837;
			fLocationZeta[19] = 0.046910077050000;
			fWeight[19]       = 0.000143110438935;
			
			fLocationKsi[20]  = 0.002097327108808;
			fLocationEta[20]  = 0.044709521721163;
			fLocationZeta[20] = 0.953089922950000;
			fWeight[20]       = 0.000000171613515;
			
			fLocationKsi[21]  = 0.010317408202534;
			fLocationEta[21]  = 0.219940124837926;
			fLocationZeta[21] = 0.769234655050000;
			fWeight[21]       = 0.000008389667268;
			
			fLocationKsi[22]  = 0.022354760860581;
			fLocationEta[22]  = 0.476544961475000;
			fLocationZeta[22] = 0.500000000000000;
			fWeight[22]       = 0.000046813628953;
			
			fLocationKsi[23]  = 0.034392113518629;
			fLocationEta[23]  = 0.733149798112074;
			fLocationZeta[23] = 0.230765344950000;
			fWeight[23]       = 0.000093222606438;
			
			fLocationKsi[24]  = 0.042612194612355;
			fLocationEta[24]  = 0.908380401228837;
			fLocationZeta[24] = 0.046910077050000;
			fWeight[24]       = 0.000070841369531;
			
			fLocationKsi[25]  = 0.000507811909541;
			fLocationEta[25]  = 0.036084856937926;
			fLocationZeta[25] = 0.953089922950000;
			fWeight[25]       = 0.000001705455115;
			
			fLocationKsi[26]  = 0.002498085653322;
			fLocationEta[26]  = 0.177512700520108;
			fLocationZeta[26] = 0.769234655050000;
			fWeight[26]       = 0.000083374558187;
			
			fLocationKsi[27]  = 0.005412610056037;
			fLocationEta[27]  = 0.384617327525000;
			fLocationZeta[27] = 0.500000000000000;
			fWeight[27]       = 0.000465222935140;
			
			fLocationKsi[28]  = 0.008327134458752;
			fLocationEta[28]  = 0.591721954529893;
			fLocationZeta[28] = 0.230765344950000;
			fWeight[28]       = 0.000926424538286;
			
			fLocationKsi[29]  = 0.010317408202533;
			fLocationEta[29]  = 0.733149798112074;
			fLocationZeta[29] = 0.046910077050000;
			fWeight[29]       = 0.000704005021598;
			
			fLocationKsi[30]  = 0.002498085653322;
			fLocationEta[30]  = 0.036084856937926;
			fLocationZeta[30] = 0.953089922950000;
			fWeight[30]       = 0.000003445281080;
			
			fLocationKsi[31]  = 0.012288864861364;
			fLocationEta[31]  = 0.177512700520108;
			fLocationZeta[31] = 0.769234655050000;
			fWeight[31]       = 0.000168429403570;
			
			fLocationKsi[32]  = 0.026626322214946;
			fLocationEta[32]  = 0.384617327525000;
			fLocationZeta[32] = 0.500000000000000;
			fWeight[32]       = 0.000939821729752;
			
			fLocationKsi[33]  = 0.040963779568528;
			fLocationEta[33]  = 0.591721954529893;
			fLocationZeta[33] = 0.230765344950000;
			fWeight[33]       = 0.001871519751699;
			
			fLocationKsi[34]  = 0.050754558776570;
			fLocationEta[34]  = 0.733149798112074;
			fLocationZeta[34] = 0.046910077050000;
			fWeight[34]       = 0.001422198191827;
			
			fLocationKsi[35]  = 0.005412610056037;
			fLocationEta[35]  = 0.036084856937926;
			fLocationZeta[35] = 0.953089922950000;
			fWeight[35]       = 0.000004094995236;
			
			fLocationKsi[36]  = 0.026626322214946;
			fLocationEta[36]  = 0.177512700520108;
			fLocationZeta[36] = 0.769234655050000;
			fWeight[36]       = 0.000200191969608;
			
			fLocationKsi[37]  = 0.057691336237500;
			fLocationEta[37]  = 0.384617327525000;
			fLocationZeta[37] = 0.500000000000000;
			fWeight[37]       = 0.001117054143585;
			
			fLocationKsi[38]  = 0.088756350260054;
			fLocationEta[38]  = 0.591721954529893;
			fLocationZeta[38] = 0.230765344950000;
			fWeight[38]       = 0.002224452603354;
			
			fLocationKsi[39]  = 0.109970062418963;
			fLocationEta[39]  = 0.733149798112074;
			fLocationZeta[39] = 0.046910077050000;
			fWeight[39]       = 0.001690397585932;
			
			fLocationKsi[40]  = 0.008327134458752;
			fLocationEta[40]  = 0.036084856937926;
			fLocationZeta[40] = 0.953089922950000;
			fWeight[40]       = 0.000003445281080;
			
			fLocationKsi[41]  = 0.040963779568529;
			fLocationEta[41]  = 0.177512700520108;
			fLocationZeta[41] = 0.769234655050000;
			fWeight[41]       = 0.000168429403570;
			
			fLocationKsi[42]  = 0.088756350260054;
			fLocationEta[42]  = 0.384617327525000;
			fLocationZeta[42] = 0.500000000000000;
			fWeight[42]       = 0.000939821729752;
			
			fLocationKsi[43]  = 0.136548920951579;
			fLocationEta[43]  = 0.591721954529893;
			fLocationZeta[43] = 0.230765344950000;
			fWeight[43]       = 0.001871519751699;
			
			fLocationKsi[44]  = 0.169185566061356;
			fLocationEta[44]  = 0.733149798112074;
			fLocationZeta[44] = 0.046910077050000;
			fWeight[44]       = 0.001422198191827;
			
			fLocationKsi[45]  = 0.010317408202534;
			fLocationEta[45]  = 0.036084856937926;
			fLocationZeta[45] = 0.953089922950000;
			fWeight[45]       = 0.000001705455115;
			
			fLocationKsi[46]  = 0.050754558776570;
			fLocationEta[46]  = 0.177512700520108;
			fLocationZeta[46] = 0.769234655050000;
			fWeight[46]       = 0.000083374558187;
			
			fLocationKsi[47]  = 0.109970062418963;
			fLocationEta[47]  = 0.384617327525000;
			fLocationZeta[47] = 0.500000000000000;
			fWeight[47]       = 0.000465222935140;
			
			fLocationKsi[48]  = 0.169185566061356;
			fLocationEta[48]  = 0.591721954529893;
			fLocationZeta[48] = 0.230765344950000;
			fWeight[48]       = 0.000926424538286;
			
			fLocationKsi[49]  = 0.209622716635392;
			fLocationEta[49]  = 0.733149798112074;
			fLocationZeta[49] = 0.046910077050000;
			fWeight[49]       = 0.000704005021598;
			
			fLocationKsi[50]  = 0.001100277664418;
			fLocationEta[50]  = 0.023455038525000;
			fLocationZeta[50] = 0.953089922950000;
			fWeight[50]       = 0.000004392061882;
			
			fLocationKsi[51]  = 0.005412610056037;
			fLocationEta[51]  = 0.115382672475000;
			fLocationZeta[51] = 0.769234655050000;
			fWeight[51]       = 0.000214714662212;
			
			fLocationKsi[52]  = 0.011727519262500;
			fLocationEta[52]  = 0.250000000000000;
			fLocationZeta[52] = 0.500000000000000;
			fWeight[52]       = 0.001198089531675;
			
			fLocationKsi[53]  = 0.018042428468963;
			fLocationEta[53]  = 0.384617327525000;
			fLocationZeta[53] = 0.230765344950000;
			fWeight[53]       = 0.002385822919231;
			
			fLocationKsi[54]  = 0.022354760860582;
			fLocationEta[54]  = 0.476544961475000;
			fLocationZeta[54] = 0.046910077050000;
			fWeight[54]       = 0.001813025504363;
			
			fLocationKsi[55]  = 0.005412610056037;
			fLocationEta[55]  = 0.023455038525000;
			fLocationZeta[55] = 0.953089922950000;
			fWeight[55]       = 0.000008872639081;
			
			fLocationKsi[56]  = 0.026626322214946;
			fLocationEta[56]  = 0.115382672475000;
			fLocationZeta[56] = 0.769234655050000;
			fWeight[56]       = 0.000433756571315;
			
			fLocationKsi[57]  = 0.057691336237500;
			fLocationEta[57]  = 0.250000000000000;
			fLocationZeta[57] = 0.500000000000000;
			fWeight[57]       = 0.002420324732526;
			
			fLocationKsi[58]  = 0.088756350260054;
			fLocationEta[58]  = 0.384617327525000;
			fLocationZeta[58] = 0.230765344950000;
			fWeight[58]       = 0.004819728464506;
			
			fLocationKsi[59]  = 0.109970062418963;
			fLocationEta[59]  = 0.476544961475000;
			fLocationZeta[59] = 0.046910077050000;
			fWeight[59]       = 0.003662589775553;
			
			fLocationKsi[60]  = 0.011727519262500;
			fLocationEta[60]  = 0.023455038525000;
			fLocationZeta[60] = 0.953089922950000;
			fWeight[60]       = 0.000010545849213;
			
			fLocationKsi[61]  = 0.057691336237500;
			fLocationEta[61]  = 0.115382672475000;
			fLocationZeta[61] = 0.769234655050000;
			fWeight[61]       = 0.000515554769544;
			
			fLocationKsi[62]  = 0.125000000000000;
			fLocationEta[62]  = 0.250000000000000;
			fLocationZeta[62] = 0.500000000000000;
			fWeight[62]       = 0.002876751713329;
			
			fLocationKsi[63]  = 0.192308663762500;
			fLocationEta[63]  = 0.384617327525000;
			fLocationZeta[63] = 0.230765344950000;
			fWeight[63]       = 0.005728637125308;
			
			fLocationKsi[64]  = 0.238272480737500;
			fLocationEta[64]  = 0.476544961475000;
			fLocationZeta[64] = 0.046910077050000;
			fWeight[64]       = 0.004353284197964;
			
			fLocationKsi[65]  = 0.018042428468963;
			fLocationEta[65]  = 0.023455038525000;
			fLocationZeta[65] = 0.953089922950000;
			fWeight[65]       = 0.000008872639081;
			
			fLocationKsi[66]  = 0.088756350260054;
			fLocationEta[66]  = 0.115382672475000;
			fLocationZeta[66] = 0.769234655050000;
			fWeight[66]       = 0.000433756571315;
			
			fLocationKsi[67]  = 0.192308663762500;
			fLocationEta[67]  = 0.250000000000000;
			fLocationZeta[67] = 0.500000000000000;
			fWeight[67]       = 0.002420324732526;
			
			fLocationKsi[68]  = 0.295860977264946;
			fLocationEta[68]  = 0.384617327525000;
			fLocationZeta[68] = 0.230765344950000;
			fWeight[68]       = 0.004819728464506;
			
			fLocationKsi[69]  = 0.366574899056037;
			fLocationEta[69]  = 0.476544961475000;
			fLocationZeta[69] = 0.046910077050000;
			fWeight[69]       = 0.003662589775553;
			
			fLocationKsi[70]  = 0.022354760860582;
			fLocationEta[70]  = 0.023455038525000;
			fLocationZeta[70] = 0.953089922950000;
			fWeight[70]       = 0.000004392061882;
			
			fLocationKsi[71]  = 0.109970062418963;
			fLocationEta[71]  = 0.115382672475000;
			fLocationZeta[71] = 0.769234655050000;
			fWeight[71]       = 0.000214714662212;
			
			fLocationKsi[72]  = 0.238272480737500;
			fLocationEta[72]  = 0.250000000000000;
			fLocationZeta[72] = 0.500000000000000;
			fWeight[72]       = 0.001198089531675;
			
			fLocationKsi[73]  = 0.366574899056037;
			fLocationEta[73]  = 0.384617327525000;
			fLocationZeta[73] = 0.230765344950000;
			fWeight[73]       = 0.002385822919231;
			
			fLocationKsi[74]  = 0.454190200614418;
			fLocationEta[74]  = 0.476544961475000;
			fLocationZeta[74] = 0.046910077050000;
			fWeight[74]       = 0.001813025504363;
			
			fLocationKsi[75]  = 0.001692743419296;
			fLocationEta[75]  = 0.010825220112074;
			fLocationZeta[75] = 0.953089922950000;
			fWeight[75]       = 0.000005684974828;
			
			fLocationKsi[76]  = 0.008327134458752;
			fLocationEta[76]  = 0.053252644429892;
			fLocationZeta[76] = 0.769234655050000;
			fWeight[76]       = 0.000277921277655;
			
			fLocationKsi[77]  = 0.018042428468963;
			fLocationEta[77]  = 0.115382672475000;
			fLocationZeta[77] = 0.500000000000000;
			fWeight[77]       = 0.001550777063649;
			
			fLocationKsi[78]  = 0.027757722479174;
			fLocationEta[78]  = 0.177512700520108;
			fLocationZeta[78] = 0.230765344950000;
			fWeight[78]       = 0.003088149393892;
			
			fLocationKsi[79]  = 0.034392113518630;
			fLocationEta[79]  = 0.219940124837926;
			fLocationZeta[79] = 0.046910077050000;
			fWeight[79]       = 0.002346734775360;
			
			fLocationKsi[80]  = 0.008327134458752;
			fLocationEta[80]  = 0.010825220112074;
			fLocationZeta[80] = 0.953089922950000;
			fWeight[80]       = 0.000011484521662;
			
			fLocationKsi[81]  = 0.040963779568529;
			fLocationEta[81]  = 0.053252644429892;
			fLocationZeta[81] = 0.769234655050000;
			fWeight[81]       = 0.000561443635237;
			
			fLocationKsi[82]  = 0.088756350260054;
			fLocationEta[82]  = 0.115382672475000;
			fLocationZeta[82] = 0.500000000000000;
			fWeight[82]       = 0.003132807676347;
			
			fLocationKsi[83]  = 0.136548920951579;
			fLocationEta[83]  = 0.177512700520108;
			fLocationZeta[83] = 0.230765344950000;
			fWeight[83]       = 0.006238535733904;
			
			fLocationKsi[84]  = 0.169185566061356;
			fLocationEta[84]  = 0.219940124837926;
			fLocationZeta[84] = 0.046910077050000;
			fWeight[84]       = 0.004740764414778;
			
			fLocationKsi[85]  = 0.018042428468963;
			fLocationEta[85]  = 0.010825220112074;
			fLocationZeta[85] = 0.953089922950000;
			fWeight[85]       = 0.000013650282925;
			
			fLocationKsi[86]  = 0.088756350260054;
			fLocationEta[86]  = 0.053252644429892;
			fLocationZeta[86] = 0.769234655050000;
			fWeight[86]       = 0.000667321173023;
			
			fLocationKsi[87]  = 0.192308663762500;
			fLocationEta[87]  = 0.115382672475000;
			fLocationZeta[87] = 0.500000000000000;
			fWeight[87]       = 0.003723595321467;
			
			fLocationKsi[88]  = 0.295860977264946;
			fLocationEta[88]  = 0.177512700520108;
			fLocationZeta[88] = 0.230765344950000;
			fWeight[88]       = 0.007415004325658;
			
			fLocationKsi[89]  = 0.366574899056037;
			fLocationEta[89]  = 0.219940124837926;
			fLocationZeta[89] = 0.046910077050000;
			fWeight[89]       = 0.005634781965175;
			
			fLocationKsi[90]  = 0.027757722479174;
			fLocationEta[90]  = 0.010825220112074;
			fLocationZeta[90] = 0.953089922950000;
			fWeight[90]       = 0.000011484521662;
			
			fLocationKsi[91]  = 0.136548920951579;
			fLocationEta[91]  = 0.053252644429892;
			fLocationZeta[91] = 0.769234655050000;
			fWeight[91]       = 0.000561443635237;
			
			fLocationKsi[92]  = 0.295860977264946;
			fLocationEta[92]  = 0.115382672475000;
			fLocationZeta[92] = 0.500000000000000;
			fWeight[92]       = 0.003132807676347;
			
			fLocationKsi[93]  = 0.455173033578314;
			fLocationEta[93]  = 0.177512700520108;
			fLocationZeta[93] = 0.230765344950000;
			fWeight[93]       = 0.006238535733904;
			
			fLocationKsi[94]  = 0.563964232050719;
			fLocationEta[94]  = 0.219940124837926;
			fLocationZeta[94] = 0.046910077050000;
			fWeight[94]       = 0.004740764414778;
			
			fLocationKsi[95]  = 0.034392113518629;
			fLocationEta[95]  = 0.010825220112074;
			fLocationZeta[95] = 0.953089922950000;
			fWeight[95]       = 0.000005684974828;
			
			fLocationKsi[96]  = 0.169185566061356;
			fLocationEta[96]  = 0.053252644429892;
			fLocationZeta[96] = 0.769234655050000;
			fWeight[96]       = 0.000277921277655;
			
			fLocationKsi[97]  = 0.366574899056037;
			fLocationEta[97]  = 0.115382672475000;
			fLocationZeta[97] = 0.500000000000000;
			fWeight[97]       = 0.001550777063649;
			
			fLocationKsi[98]  = 0.563964232050718;
			fLocationEta[98]  = 0.177512700520108;
			fLocationZeta[98] = 0.230765344950000;
			fWeight[98]       = 0.003088149393892;
			
			fLocationKsi[99]  = 0.698757684593445;
			fLocationEta[99]  = 0.219940124837926;
			fLocationZeta[99] = 0.046910077050000;
			fWeight[99]       = 0.002346734775360;
			
			fLocationKsi[100]  = 0.002097327108808;
			fLocationEta[100]  = 0.002200555328837;
			fLocationZeta[100] = 0.953089922950000;
			fWeight[100]       = 0.000003486737214;
			
			fLocationKsi[101]  = 0.010317408202534;
			fLocationEta[101]  = 0.010825220112074;
			fLocationZeta[101] = 0.769234655050000;
			fWeight[101]       = 0.000170456068985;
			
			fLocationKsi[102]  = 0.022354760860581;
			fLocationEta[102]  = 0.023455038525000;
			fLocationZeta[102] = 0.500000000000000;
			fWeight[102]       = 0.000951130350194;
			
			fLocationKsi[103]  = 0.034392113518629;
			fLocationEta[103]  = 0.036084856937926;
			fLocationZeta[103] = 0.230765344950000;
			fWeight[103]       = 0.001894039242206;
			
			fLocationKsi[104]  = 0.042612194612355;
			fLocationEta[104]  = 0.044709521721163;
			fLocationZeta[104] = 0.046910077050000;
			fWeight[104]       = 0.001439311117646;
			
			fLocationKsi[105]  = 0.010317408202534;
			fLocationEta[105]  = 0.002200555328837;
			fLocationZeta[105] = 0.953089922950000;
			fWeight[105]       = 0.000007043744305;
			
			fLocationKsi[106]  = 0.050754558776570;
			fLocationEta[106]  = 0.010825220112074;
			fLocationZeta[106] = 0.769234655050000;
			fWeight[106]       = 0.000344347420345;
			
			fLocationKsi[107]  = 0.109970062418963;
			fLocationEta[107]  = 0.023455038525000;
			fLocationZeta[107] = 0.500000000000000;
			fWeight[107]       = 0.001921429283512;
			
			fLocationKsi[108]  = 0.169185566061356;
			fLocationEta[108]  = 0.036084856937926;
			fLocationZeta[108] = 0.230765344950000;
			fWeight[108]       = 0.003826249959697;
			
			fLocationKsi[109]  = 0.209622716635392;
			fLocationEta[109]  = 0.044709521721163;
			fLocationZeta[109] = 0.046910077050000;
			fWeight[109]       = 0.002907629358023;
			
			fLocationKsi[110]  = 0.022354760860582;
			fLocationEta[110]  = 0.002200555328837;
			fLocationZeta[110] = 0.953089922950000;
			fWeight[110]       = 0.000008372059842;
			
			fLocationKsi[111]  = 0.109970062418963;
			fLocationEta[111]  = 0.010825220112074;
			fLocationZeta[111] = 0.769234655050000;
			fWeight[111]       = 0.000409284761728;
			
			fLocationKsi[112]  = 0.238272480737500;
			fLocationEta[112]  = 0.023455038525000;
			fLocationZeta[112] = 0.500000000000000;
			fWeight[112]       = 0.002283774118863;
			
			fLocationKsi[113]  = 0.366574899056037;
			fLocationEta[113]  = 0.036084856937926;
			fLocationZeta[113] = 0.230765344950000;
			fWeight[113]       = 0.004547807564524;
			
			fLocationKsi[114]  = 0.454190200614418;
			fLocationEta[114]  = 0.044709521721163;
			fLocationZeta[114] = 0.046910077050000;
			fWeight[114]       = 0.003455952676520;
			
			fLocationKsi[115]  = 0.034392113518629;
			fLocationEta[115]  = 0.002200555328837;
			fLocationZeta[115] = 0.953089922950000;
			fWeight[115]       = 0.000007043744305;
			
			fLocationKsi[116]  = 0.169185566061356;
			fLocationEta[116]  = 0.010825220112074;
			fLocationZeta[116] = 0.769234655050000;
			fWeight[116]       = 0.000344347420345;
			
			fLocationKsi[117]  = 0.366574899056037;
			fLocationEta[117]  = 0.023455038525000;
			fLocationZeta[117] = 0.500000000000000;
			fWeight[117]       = 0.001921429283512;
			
			fLocationKsi[118]  = 0.563964232050718;
			fLocationEta[118]  = 0.036084856937926;
			fLocationZeta[118] = 0.230765344950000;
			fWeight[118]       = 0.003826249959697;
			
			fLocationKsi[119]  = 0.698757684593445;
			fLocationEta[119]  = 0.044709521721163;
			fLocationZeta[119] = 0.046910077050000;
			fWeight[119]       = 0.002907629358023;
			
			fLocationKsi[120]  = 0.042612194612355;
			fLocationEta[120]  = 0.002200555328837;
			fLocationZeta[120] = 0.953089922950000;
			fWeight[120]       = 0.000003486737214;
			
			fLocationKsi[121]  = 0.209622716635392;
			fLocationEta[121]  = 0.010825220112074;
			fLocationZeta[121] = 0.769234655050000;
			fWeight[121]       = 0.000170456068985;
			
			fLocationKsi[122]  = 0.454190200614419;
			fLocationEta[122]  = 0.023455038525000;
			fLocationZeta[122] = 0.500000000000000;
			fWeight[122]       = 0.000951130350194;
			
			fLocationKsi[123]  = 0.698757684593445;
			fLocationEta[123]  = 0.036084856937926;
			fLocationZeta[123] = 0.230765344950000;
			fWeight[123]       = 0.001894039242206;
			
			fLocationKsi[124]  = 0.865768206616482;
			fLocationEta[124]  = 0.044709521721163;
			fLocationZeta[124] = 0.046910077050000;
			fWeight[124]       = 0.001439311117646;
			
			
			break;
			
		case 150://integra x8 , x2y3z3 etc.
			
			fLocationKsi[0]  = 0.000074302285193;
			fLocationEta[0]  = 0.032181312753949;
			fLocationZeta[0] = 0.966234757100000;
			fWeight[0]       = 0.000000064293054;
			
			fLocationKsi[1]  = 0.000372763745059;
			fLocationEta[1]  = 0.161448959906104;
			fLocationZeta[1] = 0.830604693200000;
			fWeight[1]       = 0.000003407437258;
			
			fLocationKsi[2]  = 0.000837730303761;
			fLocationEta[2]  = 0.362832190675434;
			fLocationZeta[2] = 0.619309593000000;
			fWeight[2]       = 0.000022321053762;
			
			fLocationKsi[3]  = 0.001362825025076;
			fLocationEta[3]  = 0.590257732274566;
			fLocationZeta[3] = 0.380690407000000;
			fWeight[3]       = 0.000059072648340;
			
			fLocationKsi[4]  = 0.001827791583778;
			fLocationEta[4]  = 0.791640963043896;
			fLocationZeta[4] = 0.169395306800000;
			fWeight[4]       = 0.000081924515509;
			
			fLocationKsi[5]  = 0.002126253043644;
			fLocationEta[5]  = 0.920908610196051;
			fLocationZeta[5] = 0.033765242900000;
			fWeight[5]       = 0.000052648937772;
			
			fLocationKsi[6]  = 0.000365516186530;
			fLocationEta[6]  = 0.032181312753949;
			fLocationZeta[6] = 0.966234757100000;
			fWeight[6]       = 0.000000129881837;
			
			fLocationKsi[7]  = 0.001833741482062;
			fLocationEta[7]  = 0.161448959906104;
			fLocationZeta[7] = 0.830604693200000;
			fWeight[7]       = 0.000006883546223;
			
			fLocationKsi[8]  = 0.004121057450330;
			fLocationEta[8]  = 0.362832190675434;
			fLocationZeta[8] = 0.619309593000000;
			fWeight[8]       = 0.000045091954356;
			
			fLocationKsi[9]  = 0.006704162661744;
			fLocationEta[9]  = 0.590257732274566;
			fLocationZeta[9] = 0.380690407000000;
			fWeight[9]       = 0.000119335815909;
			
			fLocationKsi[10]  = 0.008991478630012;
			fLocationEta[10]  = 0.791640963043896;
			fLocationZeta[10] = 0.169395306800000;
			fWeight[10]       = 0.000165500094813;
			
			fLocationKsi[11]  = 0.010459703925544;
			fLocationEta[11]  = 0.920908610196051;
			fLocationZeta[11] = 0.033765242900000;
			fWeight[11]       = 0.000106358934673;
			
			fLocationKsi[12]  = 0.000791965073025;
			fLocationEta[12]  = 0.032181312753949;
			fLocationZeta[12] = 0.966234757100000;
			fWeight[12]       = 0.000000154375069;
			
			fLocationKsi[13]  = 0.003973173446948;
			fLocationEta[13]  = 0.161448959906104;
			fLocationZeta[13] = 0.830604693200000;
			fWeight[13]       = 0.000008181651464;
			
			fLocationKsi[14]  = 0.008929108162283;
			fLocationEta[14]  = 0.362832190675434;
			fLocationZeta[14] = 0.619309593000000;
			fWeight[14]       = 0.000053595435029;
			
			fLocationKsi[15]  = 0.014525930362717;
			fLocationEta[15]  = 0.590257732274566;
			fLocationZeta[15] = 0.380690407000000;
			fWeight[15]       = 0.000141840269722;
			
			fLocationKsi[16]  = 0.019481865078052;
			fLocationEta[16]  = 0.791640963043896;
			fLocationZeta[16] = 0.169395306800000;
			fWeight[16]       = 0.000196710249212;
			
			fLocationKsi[17]  = 0.022663073451975;
			fLocationEta[17]  = 0.920908610196051;
			fLocationZeta[17] = 0.033765242900000;
			fWeight[17]       = 0.000126416196734;
			
			fLocationKsi[18]  = 0.001218413959521;
			fLocationEta[18]  = 0.032181312753949;
			fLocationZeta[18] = 0.966234757100000;
			fWeight[18]       = 0.000000129881837;
			
			fLocationKsi[19]  = 0.006112605411834;
			fLocationEta[19]  = 0.161448959906104;
			fLocationZeta[19] = 0.830604693200000;
			fWeight[19]       = 0.000006883546223;
			
			fLocationKsi[20]  = 0.013737158874236;
			fLocationEta[20]  = 0.362832190675434;
			fLocationZeta[20] = 0.619309593000000;
			fWeight[20]       = 0.000045091954356;
			
			fLocationKsi[21]  = 0.022347698063690;
			fLocationEta[21]  = 0.590257732274566;
			fLocationZeta[21] = 0.380690407000000;
			fWeight[21]       = 0.000119335815909;
			
			fLocationKsi[22]  = 0.029972251526092;
			fLocationEta[22]  = 0.791640963043896;
			fLocationZeta[22] = 0.169395306800000;
			fWeight[22]       = 0.000165500094813;
			
			fLocationKsi[23]  = 0.034866442978405;
			fLocationEta[23]  = 0.920908610196051;
			fLocationZeta[23] = 0.033765242900000;
			fWeight[23]       = 0.000106358934673;
			
			fLocationKsi[24]  = 0.001509627860858;
			fLocationEta[24]  = 0.032181312753949;
			fLocationZeta[24] = 0.966234757100000;
			fWeight[24]       = 0.000000064293054;
			
			fLocationKsi[25]  = 0.007573583148838;
			fLocationEta[25]  = 0.161448959906104;
			fLocationZeta[25] = 0.830604693200000;
			fWeight[25]       = 0.000003407437258;
			
			fLocationKsi[26]  = 0.017020486020805;
			fLocationEta[26]  = 0.362832190675434;
			fLocationZeta[26] = 0.619309593000000;
			fWeight[26]       = 0.000022321053762;
			
			fLocationKsi[27]  = 0.027689035700358;
			fLocationEta[27]  = 0.590257732274566;
			fLocationZeta[27] = 0.380690407000000;
			fWeight[27]       = 0.000059072648340;
			
			fLocationKsi[28]  = 0.037135938572325;
			fLocationEta[28]  = 0.791640963043896;
			fLocationZeta[28] = 0.169395306800000;
			fWeight[28]       = 0.000081924515509;
			
			fLocationKsi[29]  = 0.043199893860305;
			fLocationEta[29]  = 0.920908610196051;
			fLocationZeta[29] = 0.033765242900000;
			fWeight[29]       = 0.000052648937772;
			
			fLocationKsi[30]  = 0.000365516186530;
			fLocationEta[30]  = 0.025973394974861;
			fLocationZeta[30] = 0.966234757100000;
			fWeight[30]       = 0.000000638929392;
			
			fLocationKsi[31]  = 0.001833741482062;
			fLocationEta[31]  = 0.130304740393387;
			fLocationZeta[31] = 0.830604693200000;
			fWeight[31]       = 0.000033862317408;
			
			fLocationKsi[32]  = 0.004121057450330;
			fLocationEta[32]  = 0.292840253909489;
			fLocationZeta[32] = 0.619309593000000;
			fWeight[32]       = 0.000221821430613;
			
			fLocationKsi[33]  = 0.006704162661744;
			fLocationEta[33]  = 0.476394401140511;
			fLocationZeta[33] = 0.380690407000000;
			fWeight[33]       = 0.000587050212979;
			
			fLocationKsi[34]  = 0.008991478630012;
			fLocationEta[34]  = 0.638929914656613;
			fLocationZeta[34] = 0.169395306800000;
			fWeight[34]       = 0.000814146743524;
			
			fLocationKsi[35]  = 0.010459703925544;
			fLocationEta[35]  = 0.743261260075139;
			fLocationZeta[35] = 0.033765242900000;
			fWeight[35]       = 0.000523212874327;
			
			fLocationKsi[36]  = 0.001798088474243;
			fLocationEta[36]  = 0.025973394974861;
			fLocationZeta[36] = 0.966234757100000;
			fWeight[36]       = 0.000001290735432;
			
			fLocationKsi[37]  = 0.009020748041113;
			fLocationEta[37]  = 0.130304740393387;
			fLocationZeta[37] = 0.830604693200000;
			fWeight[37]       = 0.000068407078233;
			
			fLocationKsi[38]  = 0.020272770881842;
			fLocationEta[38]  = 0.292840253909489;
			fLocationZeta[38] = 0.619309593000000;
			fWeight[38]       = 0.000448113334207;
			
			fLocationKsi[39]  = 0.032979873548050;
			fLocationEta[39]  = 0.476394401140511;
			fLocationZeta[39] = 0.380690407000000;
			fWeight[39]       = 0.001185931528607;
			
			fLocationKsi[40]  = 0.044231896388779;
			fLocationEta[40]  = 0.638929914656613;
			fLocationZeta[40] = 0.169395306800000;
			fWeight[40]       = 0.001644701374281;
			
			fLocationKsi[41]  = 0.051454555955649;
			fLocationEta[41]  = 0.743261260075139;
			fLocationZeta[41] = 0.033765242900000;
			fWeight[41]       = 0.001056970307001;
			
			fLocationKsi[42]  = 0.003895923962570;
			fLocationEta[42]  = 0.025973394974861;
			fLocationZeta[42] = 0.966234757100000;
			fWeight[42]       = 0.000001534143463;
			
			fLocationKsi[43]  = 0.019545283203307;
			fLocationEta[43]  = 0.130304740393387;
			fLocationZeta[43] = 0.830604693200000;
			fWeight[43]       = 0.000081307345607;
			
			fLocationKsi[44]  = 0.043925076545255;
			fLocationEta[44]  = 0.292840253909489;
			fLocationZeta[44] = 0.619309593000000;
			fWeight[44]       = 0.000532618943492;
			
			fLocationKsi[45]  = 0.071457595929745;
			fLocationEta[45]  = 0.476394401140511;
			fLocationZeta[45] = 0.380690407000000;
			fWeight[45]       = 0.001409575546192;
			
			fLocationKsi[46]  = 0.095837389271693;
			fLocationEta[46]  = 0.638929914656613;
			fLocationZeta[46] = 0.169395306800000;
			fWeight[46]       = 0.001954860615517;
			
			fLocationKsi[47]  = 0.111486748512430;
			fLocationEta[47]  = 0.743261260075139;
			fLocationZeta[47] = 0.033765242900000;
			fWeight[47]       = 0.001256294703244;
			
			fLocationKsi[48]  = 0.005993759450896;
			fLocationEta[48]  = 0.025973394974861;
			fLocationZeta[48] = 0.966234757100000;
			fWeight[48]       = 0.000001290735432;
			
			fLocationKsi[49]  = 0.030069818365500;
			fLocationEta[49]  = 0.130304740393387;
			fLocationZeta[49] = 0.830604693200000;
			fWeight[49]       = 0.000068407078233;
			
			fLocationKsi[50]  = 0.067577382208669;
			fLocationEta[50]  = 0.292840253909489;
			fLocationZeta[50] = 0.619309593000000;
			fWeight[50]       = 0.000448113334207;
			
			fLocationKsi[51]  = 0.109935318311439;
			fLocationEta[51]  = 0.476394401140511;
			fLocationZeta[51] = 0.380690407000000;
			fWeight[51]       = 0.001185931528607;
			
			fLocationKsi[52]  = 0.147442882154608;
			fLocationEta[52]  = 0.638929914656613;
			fLocationZeta[52] = 0.169395306800000;
			fWeight[52]       = 0.001644701374281;
			
			fLocationKsi[53]  = 0.171518941069212;
			fLocationEta[53]  = 0.743261260075139;
			fLocationZeta[53] = 0.033765242900000;
			fWeight[53]       = 0.001056970307001;
			
			fLocationKsi[54]  = 0.007426331738609;
			fLocationEta[54]  = 0.025973394974861;
			fLocationZeta[54] = 0.966234757100000;
			fWeight[54]       = 0.000000638929392;
			
			fLocationKsi[55]  = 0.037256824924551;
			fLocationEta[55]  = 0.130304740393387;
			fLocationZeta[55] = 0.830604693200000;
			fWeight[55]       = 0.000033862317408;
			
			fLocationKsi[56]  = 0.083729095640181;
			fLocationEta[56]  = 0.292840253909489;
			fLocationZeta[56] = 0.619309593000000;
			fWeight[56]       = 0.000221821430613;
			
			fLocationKsi[57]  = 0.136211029197745;
			fLocationEta[57]  = 0.476394401140511;
			fLocationZeta[57] = 0.380690407000000;
			fWeight[57]       = 0.000587050212979;
			
			fLocationKsi[58]  = 0.182683299913375;
			fLocationEta[58]  = 0.638929914656613;
			fLocationZeta[58] = 0.169395306800000;
			fWeight[58]       = 0.000814146743524;
			
			fLocationKsi[59]  = 0.212513793099317;
			fLocationEta[59]  = 0.743261260075139;
			fLocationZeta[59] = 0.033765242900000;
			fWeight[59]       = 0.000523212874327;
			
			fLocationKsi[60]  = 0.000791965073025;
			fLocationEta[60]  = 0.016882621450000;
			fLocationZeta[60] = 0.966234757100000;
			fWeight[60]       = 0.000001645436108;
			
			fLocationKsi[61]  = 0.003973173446948;
			fLocationEta[61]  = 0.084697653400000;
			fLocationZeta[61] = 0.830604693200000;
			fWeight[61]       = 0.000087205692025;
			
			fLocationKsi[62]  = 0.008929108162283;
			fLocationEta[62]  = 0.190345203500000;
			fLocationZeta[62] = 0.619309593000000;
			fWeight[62]       = 0.000571257162631;
			
			fLocationKsi[63]  = 0.014525930362717;
			fLocationEta[63]  = 0.309654796500000;
			fLocationZeta[63] = 0.380690407000000;
			fWeight[63]       = 0.001511831557760;
			
			fLocationKsi[64]  = 0.019481865078052;
			fLocationEta[64]  = 0.415302346600000;
			fLocationZeta[64] = 0.169395306800000;
			fWeight[64]       = 0.002096673695528;
			
			fLocationKsi[65]  = 0.022663073451975;
			fLocationEta[65]  = 0.483117378550000;
			fLocationZeta[65] = 0.033765242900000;
			fWeight[65]       = 0.001347431135097;
			
			fLocationKsi[66]  = 0.003895923962570;
			fLocationEta[66]  = 0.016882621450000;
			fLocationZeta[66] = 0.966234757100000;
			fWeight[66]       = 0.000003324033475;
			
			fLocationKsi[67]  = 0.019545283203307;
			fLocationEta[67]  = 0.084697653400000;
			fLocationZeta[67] = 0.830604693200000;
			fWeight[67]       = 0.000176168881911;
			
			fLocationKsi[68]  = 0.043925076545255;
			fLocationEta[68]  = 0.190345203500000;
			fLocationZeta[68] = 0.619309593000000;
			fWeight[68]       = 0.001154027143043;
			
			fLocationKsi[69]  = 0.071457595929745;
			fLocationEta[69]  = 0.309654796500000;
			fLocationZeta[69] = 0.380690407000000;
			fWeight[69]       = 0.003054131777233;
			
			fLocationKsi[70]  = 0.095837389271693;
			fLocationEta[70]  = 0.415302346600000;
			fLocationZeta[70] = 0.169395306800000;
			fWeight[70]       = 0.004235602655028;
			
			fLocationKsi[71]  = 0.111486748512430;
			fLocationEta[71]  = 0.483117378550000;
			fLocationZeta[71] = 0.033765242900000;
			fWeight[71]       = 0.002722017691860;
			
			fLocationKsi[72]  = 0.008441310725000;
			fLocationEta[72]  = 0.016882621450000;
			fLocationZeta[72] = 0.966234757100000;
			fWeight[72]       = 0.000003950882650;
			
			fLocationKsi[73]  = 0.042348826700000;
			fLocationEta[73]  = 0.084697653400000;
			fLocationZeta[73] = 0.830604693200000;
			fWeight[73]       = 0.000209390965626;
			
			fLocationKsi[74]  = 0.095172601750000;
			fLocationEta[74]  = 0.190345203500000;
			fLocationZeta[74] = 0.619309593000000;
			fWeight[74]       = 0.001371654603353;
			
			fLocationKsi[75]  = 0.154827398250000;
			fLocationEta[75]  = 0.309654796500000;
			fLocationZeta[75] = 0.380690407000000;
			fWeight[75]       = 0.003630082651645;
			
			fLocationKsi[76]  = 0.207651173300000;
			fLocationEta[76]  = 0.415302346600000;
			fLocationZeta[76] = 0.169395306800000;
			fWeight[76]       = 0.005034356353545;
			
			fLocationKsi[77]  = 0.241558689275000;
			fLocationEta[77]  = 0.483117378550000;
			fLocationZeta[77] = 0.033765242900000;
			fWeight[77]       = 0.003235338198027;
			
			fLocationKsi[78]  = 0.012986697487430;
			fLocationEta[78]  = 0.016882621450000;
			fLocationZeta[78] = 0.966234757100000;
			fWeight[78]       = 0.000003324033475;
			
			fLocationKsi[79]  = 0.065152370196693;
			fLocationEta[79]  = 0.084697653400000;
			fLocationZeta[79] = 0.830604693200000;
			fWeight[79]       = 0.000176168881911;
			
			fLocationKsi[80]  = 0.146420126954745;
			fLocationEta[80]  = 0.190345203500000;
			fLocationZeta[80] = 0.619309593000000;
			fWeight[80]       = 0.001154027143043;
			
			fLocationKsi[81]  = 0.238197200570255;
			fLocationEta[81]  = 0.309654796500000;
			fLocationZeta[81] = 0.380690407000000;
			fWeight[81]       = 0.003054131777233;
			
			fLocationKsi[82]  = 0.319464957328307;
			fLocationEta[82]  = 0.415302346600000;
			fLocationZeta[82] = 0.169395306800000;
			fWeight[82]       = 0.004235602655028;
			
			fLocationKsi[83]  = 0.371630630037570;
			fLocationEta[83]  = 0.483117378550000;
			fLocationZeta[83] = 0.033765242900000;
			fWeight[83]       = 0.002722017691860;
			
			fLocationKsi[84]  = 0.016090656376975;
			fLocationEta[84]  = 0.016882621450000;
			fLocationZeta[84] = 0.966234757100000;
			fWeight[84]       = 0.000001645436108;
			
			fLocationKsi[85]  = 0.080724479953052;
			fLocationEta[85]  = 0.084697653400000;
			fLocationZeta[85] = 0.830604693200000;
			fWeight[85]       = 0.000087205692025;
			
			fLocationKsi[86]  = 0.181416095337717;
			fLocationEta[86]  = 0.190345203500000;
			fLocationZeta[86] = 0.619309593000000;
			fWeight[86]       = 0.000571257162631;
			
			fLocationKsi[87]  = 0.295128866137283;
			fLocationEta[87]  = 0.309654796500000;
			fLocationZeta[87] = 0.380690407000000;
			fWeight[87]       = 0.001511831557760;
			
			fLocationKsi[88]  = 0.395820481521948;
			fLocationEta[88]  = 0.415302346600000;
			fLocationZeta[88] = 0.169395306800000;
			fWeight[88]       = 0.002096673695528;
			
			fLocationKsi[89]  = 0.460454305098025;
			fLocationEta[89]  = 0.483117378550000;
			fLocationZeta[89] = 0.033765242900000;
			fWeight[89]       = 0.001347431135097;
			
			fLocationKsi[90]  = 0.001218413959521;
			fLocationEta[90]  = 0.007791847925139;
			fLocationZeta[90] = 0.966234757100000;
			fWeight[90]       = 0.000002129811261;
			
			fLocationKsi[91]  = 0.006112605411834;
			fLocationEta[91]  = 0.039090566406613;
			fLocationZeta[91] = 0.830604693200000;
			fWeight[91]       = 0.000112876862235;
			
			fLocationKsi[92]  = 0.013737158874236;
			fLocationEta[92]  = 0.087850153090511;
			fLocationZeta[92] = 0.619309593000000;
			fWeight[92]       = 0.000739420954639;
			
			fLocationKsi[93]  = 0.022347698063690;
			fLocationEta[93]  = 0.142915191859489;
			fLocationZeta[93] = 0.380690407000000;
			fWeight[93]       = 0.001956876879307;
			
			fLocationKsi[94]  = 0.029972251526092;
			fLocationEta[94]  = 0.191674778543387;
			fLocationZeta[94] = 0.169395306800000;
			fWeight[94]       = 0.002713881885301;
			
			fLocationKsi[95]  = 0.034866442978405;
			fLocationEta[95]  = 0.222973497024861;
			fLocationZeta[95] = 0.033765242900000;
			fWeight[95]       = 0.001744081092365;
			
			fLocationKsi[96]  = 0.005993759450896;
			fLocationEta[96]  = 0.007791847925139;
			fLocationZeta[96] = 0.966234757100000;
			fWeight[96]       = 0.000004302545624;
			
			fLocationKsi[97]  = 0.030069818365500;
			fLocationEta[97]  = 0.039090566406613;
			fLocationZeta[97] = 0.830604693200000;
			fWeight[97]       = 0.000228028585655;
			
			fLocationKsi[98]  = 0.067577382208669;
			fLocationEta[98]  = 0.087850153090511;
			fLocationZeta[98] = 0.619309593000000;
			fWeight[98]       = 0.001493743812082;
			
			fLocationKsi[99]  = 0.109935318311439;
			fLocationEta[99]  = 0.142915191859489;
			fLocationZeta[99] = 0.380690407000000;
			fWeight[99]       = 0.003953191630739;
			
			fLocationKsi[100]  = 0.147442882154608;
			fLocationEta[100]  = 0.191674778543387;
			fLocationZeta[100] = 0.169395306800000;
			fWeight[100]       = 0.005482457925297;
			
			fLocationKsi[101]  = 0.171518941069212;
			fLocationEta[101]  = 0.222973497024861;
			fLocationZeta[101] = 0.033765242900000;
			fWeight[101]       = 0.003523311481973;
			
			fLocationKsi[102]  = 0.012986697487430;
			fLocationEta[102]  = 0.007791847925139;
			fLocationZeta[102] = 0.966234757100000;
			fWeight[102]       = 0.000005113923486;
			
			fLocationKsi[103]  = 0.065152370196693;
			fLocationEta[103]  = 0.039090566406613;
			fLocationZeta[103] = 0.830604693200000;
			fWeight[103]       = 0.000271030418215;
			
			fLocationKsi[104]  = 0.146420126954745;
			fLocationEta[104]  = 0.087850153090511;
			fLocationZeta[104] = 0.619309593000000;
			fWeight[104]       = 0.001775435342594;
			
			fLocationKsi[105]  = 0.238197200570255;
			fLocationEta[105]  = 0.142915191859489;
			fLocationZeta[105] = 0.380690407000000;
			fWeight[105]       = 0.004698688008274;
			
			fLocationKsi[106]  = 0.319464957328307;
			fLocationEta[106]  = 0.191674778543387;
			fLocationZeta[106] = 0.169395306800000;
			fWeight[106]       = 0.006516344694538;
			
			fLocationKsi[107]  = 0.371630630037570;
			fLocationEta[107]  = 0.222973497024861;
			fLocationZeta[107] = 0.033765242900000;
			fWeight[107]       = 0.004187740680475;
			
			fLocationKsi[108]  = 0.019979635523965;
			fLocationEta[108]  = 0.007791847925139;
			fLocationZeta[108] = 0.966234757100000;
			fWeight[108]       = 0.000004302545624;
			
			fLocationKsi[109]  = 0.100234922027887;
			fLocationEta[109]  = 0.039090566406613;
			fLocationZeta[109] = 0.830604693200000;
			fWeight[109]       = 0.000228028585655;
			
			fLocationKsi[110]  = 0.225262871700820;
			fLocationEta[110]  = 0.087850153090511;
			fLocationZeta[110] = 0.619309593000000;
			fWeight[110]       = 0.001493743812082;
			
			fLocationKsi[111]  = 0.366459082829073;
			fLocationEta[111]  = 0.142915191859489;
			fLocationZeta[111] = 0.380690407000000;
			fWeight[111]       = 0.003953191630739;
			
			fLocationKsi[112]  = 0.491487032502006;
			fLocationEta[112]  = 0.191674778543387;
			fLocationZeta[112] = 0.169395306800000;
			fWeight[112]       = 0.005482457925297;
			
			fLocationKsi[113]  = 0.571742319005928;
			fLocationEta[113]  = 0.222973497024861;
			fLocationZeta[113] = 0.033765242900000;
			fWeight[113]       = 0.003523311481973;
			
			fLocationKsi[114]  = 0.024754981015340;
			fLocationEta[114]  = 0.007791847925139;
			fLocationZeta[114] = 0.966234757100000;
			fWeight[114]       = 0.000002129811261;
			
			fLocationKsi[115]  = 0.124192134981553;
			fLocationEta[115]  = 0.039090566406613;
			fLocationZeta[115] = 0.830604693200000;
			fWeight[115]       = 0.000112876862235;
			
			fLocationKsi[116]  = 0.279103095035253;
			fLocationEta[116]  = 0.087850153090511;
			fLocationZeta[116] = 0.619309593000000;
			fWeight[116]       = 0.000739420954639;
			
			fLocationKsi[117]  = 0.454046703076821;
			fLocationEta[117]  = 0.142915191859489;
			fLocationZeta[117] = 0.380690407000000;
			fWeight[117]       = 0.001956876879307;
			
			fLocationKsi[118]  = 0.608957663130521;
			fLocationEta[118]  = 0.191674778543387;
			fLocationZeta[118] = 0.169395306800000;
			fWeight[118]       = 0.002713881885301;
			
			fLocationKsi[119]  = 0.708394817096734;
			fLocationEta[119]  = 0.222973497024861;
			fLocationZeta[119] = 0.033765242900000;
			fWeight[119]       = 0.001744081092365;
			
			fLocationKsi[120]  = 0.001509627860858;
			fLocationEta[120]  = 0.001583930146051;
			fLocationZeta[120] = 0.966234757100000;
			fWeight[120]       = 0.000001306266502;
			
			fLocationKsi[121]  = 0.007573583148838;
			fLocationEta[121]  = 0.007946346893896;
			fLocationZeta[121] = 0.830604693200000;
			fWeight[121]       = 0.000069230202086;
			
			fLocationKsi[122]  = 0.017020486020805;
			fLocationEta[122]  = 0.017858216324566;
			fLocationZeta[122] = 0.619309593000000;
			fWeight[122]       = 0.000453505360642;
			
			fLocationKsi[123]  = 0.027689035700358;
			fLocationEta[123]  = 0.029051860725434;
			fLocationZeta[123] = 0.380690407000000;
			fWeight[123]       = 0.001200201521631;
			
			fLocationKsi[124]  = 0.037135938572325;
			fLocationEta[124]  = 0.038963730156104;
			fLocationZeta[124] = 0.169395306800000;
			fWeight[124]       = 0.001664491620657;
			
			fLocationKsi[125]  = 0.043199893860305;
			fLocationEta[125]  = 0.045326146903949;
			fLocationZeta[125] = 0.033765242900000;
			fWeight[125]       = 0.001069688544558;
			
			fLocationKsi[126]  = 0.007426331738609;
			fLocationEta[126]  = 0.001583930146051;
			fLocationZeta[126] = 0.966234757100000;
			fWeight[126]       = 0.000002638858816;
			
			fLocationKsi[127]  = 0.037256824924551;
			fLocationEta[127]  = 0.007946346893896;
			fLocationZeta[127] = 0.830604693200000;
			fWeight[127]       = 0.000139855633419;
			
			fLocationKsi[128]  = 0.083729095640181;
			fLocationEta[128]  = 0.017858216324566;
			fLocationZeta[128] = 0.619309593000000;
			fWeight[128]       = 0.000916150430896;
			
			fLocationKsi[129]  = 0.136211029197745;
			fLocationEta[129]  = 0.029051860725434;
			fLocationZeta[129] = 0.380690407000000;
			fWeight[129]       = 0.002424591276377;
			
			fLocationKsi[130]  = 0.182683299913375;
			fLocationEta[130]  = 0.038963730156104;
			fLocationZeta[130] = 0.169395306800000;
			fWeight[130]       = 0.003362528534012;
			
			fLocationKsi[131]  = 0.212513793099317;
			fLocationEta[131]  = 0.045326146903949;
			fLocationZeta[131] = 0.033765242900000;
			fWeight[131]       = 0.002160935032020;
			
			fLocationKsi[132]  = 0.016090656376975;
			fLocationEta[132]  = 0.001583930146051;
			fLocationZeta[132] = 0.966234757100000;
			fWeight[132]       = 0.000003136497148;
			
			fLocationKsi[133]  = 0.080724479953052;
			fLocationEta[133]  = 0.007946346893896;
			fLocationZeta[133] = 0.830604693200000;
			fWeight[133]       = 0.000166229732585;
			
			fLocationKsi[134]  = 0.181416095337717;
			fLocationEta[134]  = 0.017858216324566;
			fLocationZeta[134] = 0.619309593000000;
			fWeight[134]       = 0.001088918890234;
			
			fLocationKsi[135]  = 0.295128866137283;
			fLocationEta[135]  = 0.029051860725434;
			fLocationZeta[135] = 0.380690407000000;
			fWeight[135]       = 0.002881822845797;
			
			fLocationKsi[136]  = 0.395820481521948;
			fLocationEta[136]  = 0.038963730156104;
			fLocationZeta[136] = 0.169395306800000;
			fWeight[136]       = 0.003996637141845;
			
			fLocationKsi[137]  = 0.460454305098025;
			fLocationEta[137]  = 0.045326146903949;
			fLocationZeta[137] = 0.033765242900000;
			fWeight[137]       = 0.002568446073461;
			
			fLocationKsi[138]  = 0.024754981015340;
			fLocationEta[138]  = 0.001583930146051;
			fLocationZeta[138] = 0.966234757100000;
			fWeight[138]       = 0.000002638858816;
			
			fLocationKsi[139]  = 0.124192134981553;
			fLocationEta[139]  = 0.007946346893896;
			fLocationZeta[139] = 0.830604693200000;
			fWeight[139]       = 0.000139855633419;
			
			fLocationKsi[140]  = 0.279103095035253;
			fLocationEta[140]  = 0.017858216324566;
			fLocationZeta[140] = 0.619309593000000;
			fWeight[140]       = 0.000916150430896;
			
			fLocationKsi[141]  = 0.454046703076821;
			fLocationEta[141]  = 0.029051860725434;
			fLocationZeta[141] = 0.380690407000000;
			fWeight[141]       = 0.002424591276377;
			
			fLocationKsi[142]  = 0.608957663130521;
			fLocationEta[142]  = 0.038963730156104;
			fLocationZeta[142] = 0.169395306800000;
			fWeight[142]       = 0.003362528534012;
			
			fLocationKsi[143]  = 0.708394817096734;
			fLocationEta[143]  = 0.045326146903949;
			fLocationZeta[143] = 0.033765242900000;
			fWeight[143]       = 0.002160935032020;
			
			fLocationKsi[144]  = 0.030671684893091;
			fLocationEta[144]  = 0.001583930146051;
			fLocationZeta[144] = 0.966234757100000;
			fWeight[144]       = 0.000001306266502;
			
			fLocationKsi[145]  = 0.153875376757266;
			fLocationEta[145]  = 0.007946346893896;
			fLocationZeta[145] = 0.830604693200000;
			fWeight[145]       = 0.000069230202086;
			
			fLocationKsi[146]  = 0.345811704654629;
			fLocationEta[146]  = 0.017858216324566;
			fLocationZeta[146] = 0.619309593000000;
			fWeight[146]       = 0.000453505360642;
			
			fLocationKsi[147]  = 0.562568696574208;
			fLocationEta[147]  = 0.029051860725434;
			fLocationZeta[147] = 0.380690407000000;
			fWeight[147]       = 0.001200201521631;
			
			fLocationKsi[148]  = 0.754505024471571;
			fLocationEta[148]  = 0.038963730156104;
			fLocationZeta[148] = 0.169395306800000;
			fWeight[148]       = 0.001664491620657;
			
			fLocationKsi[149]  = 0.877708716335746;
			fLocationEta[149]  = 0.045326146903949;
			fLocationZeta[149] = 0.033765242900000;
			fWeight[149]       = 0.001069688544558;
			
			
			break;
			
		default:
			PZError << "TPZIntRuleT3D creation : invalid number of integration points "
			" specified\n";
			fNumInt = 0;
	}
}


#include "pzquad.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include <math.h>
#include <stdlib.h>

//#pragma segment UTIL
//***** number of integration rules in PZINTVEC

#define NUMINT_RULES    23
#define NUMINT_RULEST   19
#define NUMINT_RULEST3D  8
#define NUMINT_RULESP3D  8

TPZIntRuleList  gIntRuleList;


int TPZInt1Point::GetMaxOrder() {return 1;}
int TPZInt1d::GetMaxOrder() {return NUMINT_RULES - 1;}
int TPZIntTriang::GetMaxOrder() {return NUMINT_RULEST - 1 ;}
int TPZIntQuad::GetMaxOrder() {return NUMINT_RULES - 2 ;}
int TPZIntCube3D::GetMaxOrder() {return NUMINT_RULES - 3;}
int TPZIntTetra3D::GetMaxOrder() {return NUMINT_RULEST3D - 1;}
int TPZIntPyram3D::GetMaxOrder() {return NUMINT_RULEST3D - 1;}
int TPZIntPrism3D::GetMaxOrder() {return NUMINT_RULESP3D - 1;}



//***************************************
//***************************************
TPZIntRule::TPZIntRule(int precision){

  if(precision < 1 && precision > NUMINT_RULES){
    PZError << "TPZIntRule creation precision = " << precision << " not allowable\n";
    fNumInt = 0;
    fLocation = NULL;
    fWeight = NULL;
    return;
  }
  int intpoints[] = {1,1,2,3,4,5,6,7,8,9,10,12,12};
  int numpoints = (precision+1) >> 1;
  if(!(precision%2)) numpoints++;
  if(numpoints > 12) exit(-1);
  fNumInt = (short) intpoints[numpoints];
  fLocation = new REAL[fNumInt];
  fWeight = new REAL[fNumInt];

  if(fLocation == NULL || fWeight == NULL){
    fNumInt = 0;
    return;
  }

  switch(fNumInt){

  case 1:
    fLocation[0] = 0.;
    fWeight[0] = 2.;
    break;

  case 2:
    fLocation[0] = -1./sqrt(3.);
    fWeight[0] = 1.;
    fLocation[1] = 1./sqrt(3.);
    fWeight[1] = 1.;
    break;

  case 3:
    fLocation[0] = -sqrt(3./5.);
    fWeight[0] = 5./9.;
    fLocation[1] = 0.;
    fWeight[1] = 8./9.;
    fLocation[2] = sqrt(3./5.);
    fWeight[2] = 5./9.;
    break;

  case 4:
    fLocation[0] = -0.8611363116;
    fWeight[0] =    0.3478548451;
    fLocation[1] = -0.3399810436;
    fWeight[1] =    0.6521451549;
    fLocation[2] =  0.3399810436;
    fWeight[2] =    0.6521451549;
    fLocation[3] =  0.8611363116;
    fWeight[3] =    0.3478548451;
    break;

  case 5:
    fLocation[0] = -0.9061798459;
    fWeight[0] =    0.2369268850;
    fLocation[1] = -0.5384693101;
    fWeight[1] =    0.4786286704;
    fLocation[2] =  0.0000000000;
    fWeight[2] =    0.5688888888;
    fLocation[3] =  0.5384693101;
    fWeight[3] =    0.4786286704;
    fLocation[4] =  0.9061798459;
    fWeight[4] =    0.2369268850;
    break;

  case 6:
    fLocation[0] = -0.9324695142;
    fWeight[0] =    0.1713244923;
    fLocation[1] = -0.6612093864;
    fWeight[1] =    0.3607615730;
    fLocation[2] = -0.2386191860;
    fWeight[2] =    0.4679139345;
    fLocation[3] =  0.2386191860;
    fWeight[3] =    0.4679139345;
    fLocation[4] =  0.6612093864;
    fWeight[4] =    0.3607615730;
    fLocation[5] =  0.9324695142;
    fWeight[5] =    0.1713244923;
    break;

  case 7:
    fLocation[0] = -0.9491079123;
    fWeight[0] =    0.1294849661;
    fLocation[1] = -0.7415311855;
    fWeight[1] =    0.2797053914;
    fLocation[2] = -0.4058451513;
    fWeight[2] =    0.3818300505;
    fLocation[3] =  0.0000000000;
    fWeight[3] =    0.4179591836;
    fLocation[4] =  0.4058451513;
    fWeight[4] =    0.3818300505;
    fLocation[5] =  0.7415311855;
    fWeight[5] =    0.2797053914;
    fLocation[6] =  0.9491079123;
    fWeight[6] =    0.1294849661;
    break;

  case 8:
    fLocation[0] = -0.9602898564;
    fWeight[0] =    0.1012285362;
    fLocation[1] = -0.7966664774;
    fWeight[1] =    0.2223810344;
    fLocation[2] = -0.5255324099;
    fWeight[2] =    0.3137066458;
    fLocation[3] = -0.1834346424;
    fWeight[3] =    0.3626837833;
    fLocation[4] =  0.1834346424;
    fWeight[4] =    0.3626837833;
    fLocation[5] =  0.5255324099;
    fWeight[5] =    0.3137066458;
    fLocation[6] =  0.7966664774;
    fWeight[6] =    0.2223810344;
    fLocation[7] =  0.9602898564;
    fWeight[7] =    0.1012285362;
    break;

  case 9:
    fLocation[0] = -0.9681602395;
    fWeight[0] =    0.0812743883;
    fLocation[1] = -0.8360311073;
    fWeight[1] =    0.1806481606;
    fLocation[2] = -0.6133714327;
    fWeight[2] =    0.2606106964;
    fLocation[3] = -0.3242534234;
    fWeight[3] =    0.3123470770;
    fLocation[4] =  0.0000000000;
    fWeight[4] =    0.3302393550;
    fLocation[5] =  0.3242534234;
    fWeight[5] =    0.3123470770;
    fLocation[6] =  0.6133714327;
    fWeight[6] =    0.2606106964;
    fLocation[7] =  0.8360311073;
    fWeight[7] =    0.1806481606;
    fLocation[8] =  0.9681602395;
    fWeight[8] =    0.0812743883;
    break;

  case 10:
    fLocation[0] = -0.9739065285;
    fWeight[0] =    0.0666713443;
    fLocation[1] = -0.8650633666;
    fWeight[1] =    0.1494513491;
    fLocation[2] = -0.6794095682;
    fWeight[2] =    0.2190863625;
    fLocation[3] = -0.4333953941;
    fWeight[3] =    0.2692667193;
    fLocation[4] = -0.1488743389;
    fWeight[4] =    0.2955242247;
    fLocation[5] =  0.1488743389;
    fWeight[5] =    0.2955242247;
    fLocation[6] =  0.4333953941;
    fWeight[6] =    0.2692667193;
    fLocation[7] =  0.6794095682;
    fWeight[7] =    0.2190863625;
    fLocation[8] =  0.8650633666;
    fWeight[8] =    0.1494513491;
    fLocation[9] =  0.9739065285;
    fWeight[9] =    0.0666713443;
    break;

  case 11:
  case 12:
    fLocation[0] = -0.9815606342;
    fWeight[0] =    0.0471753363;
    fLocation[1] = -0.9041172563;
    fWeight[1] =    0.1069393259;
    fLocation[2] = -0.7699026741;
    fWeight[2] =    0.1600783285;
    fLocation[3] = -0.5873179542;
    fWeight[3] =    0.2031674267;
    fLocation[4] = -0.3678314989;
    fWeight[4] =    0.2334925365;
    fLocation[5] = -0.1252334085;
    fWeight[5] =    0.2491470458;
    fLocation[6] =  0.1252334085;
    fWeight[6] =    0.2491470458;
    fLocation[7] =  0.3678314989;
    fWeight[7] =    0.2334925365;
    fLocation[8] =  0.5873179542;
    fWeight[8] =    0.2031674267;
    fLocation[9] =  0.7699026741;
    fWeight[9]  =    0.1600783285;
    fLocation[10] =  0.9041172563;
    fWeight[10] =    0.1069393259;
    fLocation[11] =  0.9815606342;
    fWeight[11] =    0.0471753363;
    fNumInt=12;
    break;

  default:
    PZError << "TPZIntRule1D creation : invalid number of integration points "
      " specified\n";
    //			PZError.show();

    fNumInt = 0;
  }
}

//***************************************
//***************************************
TPZIntRule::~TPZIntRule(){

  if (fLocation) delete []fLocation;
  if (fWeight)   delete []fWeight;

}

//***************************************
//***************************************
REAL TPZIntRule::Loc(int i) {

  if (fLocation && i>=0 && i<fNumInt)
    return fLocation[i];
  else {
    PZError << "ERROR(TPZIntRule::loc) Out of bounds!!\n";
    return 0.0;
  }
}

//***************************************
//***************************************
REAL TPZIntRule::W(int i) {

  if (fWeight && i>=0 && i<fNumInt)
    return fWeight[i];
  else {
    PZError << "ERROR(TPZIntRule::w) Out of bounds!!\n";
    return 0.0;
  }
}

//***************************************
//***************************************
TPZIntRuleT::TPZIntRuleT(int precision){

  if(precision < 0 && precision >= NUMINT_RULEST){
    PZError << "TPZIntRule creation precision = " << precision << " not available\n";
    //		PZError.show();
    precision = NUMINT_RULEST-1;
    PZError << "TPZIntRule creation precision gotten = " << precision << "\n";
  }

  fNumInt = (short) precision;
  int rule = 0;
  if (precision ==  0 || precision == 1) rule = 1;//fNumInt =  1;//integra constantes
  else if (precision <=  2) rule = 2;
  else if (precision <=  4) rule = 4;//fNumInt =  3;//integra até grau  2 : cte, x, y, xy, x2, y2
  else if (precision <=  6) rule = 6;//fNumInt =  4;//integra até grau  4 : anteriores + x3,x4,xy2,x2y,y3,y4
  else if (precision <=  8) rule = 8;//fNumInt = 13;//integra até grau  6 : anteriores + etc
  else if (precision <=  10) rule = 10;//fNumInt = 25;//integra até grau  8
  else if (precision <=  12) rule = 12;//fNumInt = 30;//integra até grau 10
  else if (precision <=  14) rule = 14;
  else if (precision <=  16) rule = 17;
  else if (precision <=  18) rule = 19;

  int NRGAUPO [20] = {1,1,3,4,6,7,12,13,16,19,25,27,33,37,42,48,52,61,70,73};
  fNumInt = NRGAUPO[rule];
  //  cout << "fNumInt " << fNumInt << "\t" << "precision" << precision << endl; 
  
  fLocationKsi = new REAL[fNumInt];
  fLocationEta = new REAL[fNumInt];
  fWeight      = new REAL[fNumInt];

  if(fLocationKsi == NULL || fLocationEta == NULL || fWeight == NULL){
    fNumInt = 0;
    return;
  }
  //rule = precision;
  //  switch(/*fNumInt*/NSELECT[precision]){
  switch (rule){
  case 0:
  case 1:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.5;
    break;

  case 2:
    fLocationKsi[0]= 0.166666666667;
    fLocationEta[0]= 0.166666666667;
    fWeight[0]= 0.166666666667;
    fLocationKsi[1]= 0.166666666667;
    fLocationEta[1]= 0.666666666667;
    fWeight[1]= 0.166666666667;
    fLocationKsi[2]= 0.666666666667;
    fLocationEta[2]= 0.166666666667;
    fWeight[2]= 0.166666666667;
    break;

  case 3:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= -0.28125;
    fLocationKsi[1]= 0.200000000000;
    fLocationEta[1]= 0.200000000000;
    fWeight[1]= 0.260417;
    fLocationKsi[2]= 0.200000000000;
    fLocationEta[2]= 0.600000000000;
    fWeight[2]= 0.260417;
    fLocationKsi[3]= 0.600000000000;
    fLocationEta[3]= 0.200000000000;
    fWeight[3]= 0.260417;
    break;

  case 4:
    fLocationKsi[0]= 0.445948490916;
    fLocationEta[0]= 0.445948490916;
    fWeight[0]= 0.111691;
    fLocationKsi[1]= 0.445948490916;
    fLocationEta[1]= 0.108103018168;
    fWeight[1]= 0.111691;
    fLocationKsi[2]= 0.108103018168;
    fLocationEta[2]= 0.445948490916;
    fWeight[2]= 0.111691;
    fLocationKsi[3]= 0.091576213510;
    fLocationEta[3]= 0.091576213510;
    fWeight[3]= 0.0549759;
    fLocationKsi[4]= 0.091576213510;
    fLocationEta[4]= 0.816847572980;
    fWeight[4]= 0.0549759;
    fLocationKsi[5]= 0.816847572980;
    fLocationEta[5]= 0.091576213510;
    fWeight[5]= 0.0549759;
    break;

  case 5:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.1125;
    fLocationKsi[1]= 0.470142064105;
    fLocationEta[1]= 0.470142064105;
    fWeight[1]= 0.0661971;
    fLocationKsi[2]= 0.470142064105;
    fLocationEta[2]= 0.059715871790;
    fWeight[2]= 0.0661971;
    fLocationKsi[3]= 0.059715871790;
    fLocationEta[3]= 0.470142064105;
    fWeight[3]= 0.0661971;
    fLocationKsi[4]= 0.101286507323;
    fLocationEta[4]= 0.101286507323;
    fWeight[4]= 0.0629696;
    fLocationKsi[5]= 0.101286507323;
    fLocationEta[5]= 0.797426985353;
    fWeight[5]= 0.0629696;
    fLocationKsi[6]= 0.797426985353;
    fLocationEta[6]= 0.101286507323;
    fWeight[6]= 0.0629696;
    break;

  case 6:
    fLocationKsi[0]= 0.249286745171;
    fLocationEta[0]= 0.249286745171;
    fWeight[0]= 0.0583931;
    fLocationKsi[1]= 0.249286745171;
    fLocationEta[1]= 0.501426509658;
    fWeight[1]= 0.0583931;
    fLocationKsi[2]= 0.501426509658;
    fLocationEta[2]= 0.249286745171;
    fWeight[2]= 0.0583931;
    fLocationKsi[3]= 0.063089014492;
    fLocationEta[3]= 0.063089014492;
    fWeight[3]= 0.0254225;
    fLocationKsi[4]= 0.063089014492;
    fLocationEta[4]= 0.873821971017;
    fWeight[4]= 0.0254225;
    fLocationKsi[5]= 0.873821971017;
    fLocationEta[5]= 0.063089014492;
    fWeight[5]= 0.0254225;
    fLocationKsi[6]= 0.310352451034;
    fLocationEta[6]= 0.636502499121;
    fWeight[6]= 0.0414255;
    fLocationKsi[7]= 0.636502499121;
    fLocationEta[7]= 0.053145049845;
    fWeight[7]= 0.0414255;
    fLocationKsi[8]= 0.053145049845;
    fLocationEta[8]= 0.310352451034;
    fWeight[8]= 0.0414255;
    fLocationKsi[9]= 0.310352451034;
    fLocationEta[9]= 0.053145049845;
    fWeight[9]= 0.0414255;
    fLocationKsi[10]= 0.636502499121;
    fLocationEta[10]= 0.310352451034;
    fWeight[10]= 0.0414255;
    fLocationKsi[11]= 0.053145049845;
    fLocationEta[11]= 0.636502499121;
    fWeight[11]= 0.0414255;
    break;

  case 7:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= -0.074785;
    fLocationKsi[1]= 0.260345966079;
    fLocationEta[1]= 0.260345966079;
    fWeight[1]= 0.0878076;
    fLocationKsi[2]= 0.260345966079;
    fLocationEta[2]= 0.479308067842;
    fWeight[2]= 0.0878076;
    fLocationKsi[3]= 0.479308067842;
    fLocationEta[3]= 0.260345966079;
    fWeight[3]= 0.0878076;
    fLocationKsi[4]= 0.065130102902;
    fLocationEta[4]= 0.065130102902;
    fWeight[4]= 0.0266736;
    fLocationKsi[5]= 0.065130102902;
    fLocationEta[5]= 0.869739794196;
    fWeight[5]= 0.0266736;
    fLocationKsi[6]= 0.869739794196;
    fLocationEta[6]= 0.065130102902;
    fWeight[6]= 0.0266736;
    fLocationKsi[7]= 0.312865496005;
    fLocationEta[7]= 0.638444188570;
    fWeight[7]= 0.0385569;
    fLocationKsi[8]= 0.638444188570;
    fLocationEta[8]= 0.048690315425;
    fWeight[8]= 0.0385569;
    fLocationKsi[9]= 0.048690315425;
    fLocationEta[9]= 0.312865496005;
    fWeight[9]= 0.0385569;
    fLocationKsi[10]= 0.312865496005;
    fLocationEta[10]= 0.048690315425;
    fWeight[10]= 0.0385569;
    fLocationKsi[11]= 0.638444188570;
    fLocationEta[11]= 0.312865496005;
    fWeight[11]= 0.0385569;
    fLocationKsi[12]= 0.048690315425;
    fLocationEta[12]= 0.638444188570;
    fWeight[12]= 0.0385569;
    break;

  case 8:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.0721578;
    fLocationKsi[1]= 0.459292588293;
    fLocationEta[1]= 0.459292588293;
    fWeight[1]= 0.0475458;
    fLocationKsi[2]= 0.459292588293;
    fLocationEta[2]= 0.081414823415;
    fWeight[2]= 0.0475458;
    fLocationKsi[3]= 0.081414823415;
    fLocationEta[3]= 0.459292588293;
    fWeight[3]= 0.0475458;
    fLocationKsi[4]= 0.170569307752;
    fLocationEta[4]= 0.170569307752;
    fWeight[4]= 0.0516087;
    fLocationKsi[5]= 0.170569307752;
    fLocationEta[5]= 0.658861384496;
    fWeight[5]= 0.0516087;
    fLocationKsi[6]= 0.658861384496;
    fLocationEta[6]= 0.170569307752;
    fWeight[6]= 0.0516087;
    fLocationKsi[7]= 0.050547228317;
    fLocationEta[7]= 0.050547228317;
    fWeight[7]= 0.0162292;
    fLocationKsi[8]= 0.050547228317;
    fLocationEta[8]= 0.898905543366;
    fWeight[8]= 0.0162292;
    fLocationKsi[9]= 0.898905543366;
    fLocationEta[9]= 0.050547228317;
    fWeight[9]= 0.0162292;
    fLocationKsi[10]= 0.263112829635;
    fLocationEta[10]= 0.728492392955;
    fWeight[10]= 0.0136152;
    fLocationKsi[11]= 0.728492392955;
    fLocationEta[11]= 0.008394777410;
    fWeight[11]= 0.0136152;
    fLocationKsi[12]= 0.008394777410;
    fLocationEta[12]= 0.263112829635;
    fWeight[12]= 0.0136152;
    fLocationKsi[13]= 0.263112829635;
    fLocationEta[13]= 0.008394777410;
    fWeight[13]= 0.0136152;
    fLocationKsi[14]= 0.728492392955;
    fLocationEta[14]= 0.263112829635;
    fWeight[14]= 0.0136152;
    fLocationKsi[15]= 0.008394777410;
    fLocationEta[15]= 0.728492392955;
    fWeight[15]= 0.0136152;
    break;

  case 9:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.0485679;
    fLocationKsi[1]= 0.489682519199;
    fLocationEta[1]= 0.489682519199;
    fWeight[1]= 0.0156674;
    fLocationKsi[2]= 0.489682519199;
    fLocationEta[2]= 0.020634961603;
    fWeight[2]= 0.0156674;
    fLocationKsi[3]= 0.020634961603;
    fLocationEta[3]= 0.489682519199;
    fWeight[3]= 0.0156674;
    fLocationKsi[4]= 0.437089591493;
    fLocationEta[4]= 0.437089591493;
    fWeight[4]= 0.0389138;
    fLocationKsi[5]= 0.437089591493;
    fLocationEta[5]= 0.125820817014;
    fWeight[5]= 0.0389138;
    fLocationKsi[6]= 0.125820817014;
    fLocationEta[6]= 0.437089591493;
    fWeight[6]= 0.0389138;
    fLocationKsi[7]= 0.188203535619;
    fLocationEta[7]= 0.188203535619;
    fWeight[7]= 0.0398239;
    fLocationKsi[8]= 0.188203535619;
    fLocationEta[8]= 0.623592928762;
    fWeight[8]= 0.0398239;
    fLocationKsi[9]= 0.623592928762;
    fLocationEta[9]= 0.188203535619;
    fWeight[9]= 0.0398239;
    fLocationKsi[10]= 0.044729513394;
    fLocationEta[10]= 0.044729513394;
    fWeight[10]= 0.0127888;
    fLocationKsi[11]= 0.044729513394;
    fLocationEta[11]= 0.910540973211;
    fWeight[11]= 0.0127888;
    fLocationKsi[12]= 0.910540973211;
    fLocationEta[12]= 0.044729513394;
    fWeight[12]= 0.0127888;
    fLocationKsi[13]= 0.221962989161;
    fLocationEta[13]= 0.741198598784;
    fWeight[13]= 0.0216418;
    fLocationKsi[14]= 0.741198598784;
    fLocationEta[14]= 0.036838412055;
    fWeight[14]= 0.0216418;
    fLocationKsi[15]= 0.036838412055;
    fLocationEta[15]= 0.221962989161;
    fWeight[15]= 0.0216418;
    fLocationKsi[16]= 0.221962989161;
    fLocationEta[16]= 0.036838412055;
    fWeight[16]= 0.0216418;
    fLocationKsi[17]= 0.741198598784;
    fLocationEta[17]= 0.221962989161;
    fWeight[17]= 0.0216418;
    fLocationKsi[18]= 0.036838412055;
    fLocationEta[18]= 0.741198598784;
    fWeight[18]= 0.0216418;
    break;

  case 10:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.045409;
    fLocationKsi[1]= 0.485577633384;
    fLocationEta[1]= 0.485577633384;
    fWeight[1]= 0.018363;
    fLocationKsi[2]= 0.485577633384;
    fLocationEta[2]= 0.028844733233;
    fWeight[2]= 0.018363;
    fLocationKsi[3]= 0.028844733233;
    fLocationEta[3]= 0.485577633384;
    fWeight[3]= 0.018363;
    fLocationKsi[4]= 0.109481575485;
    fLocationEta[4]= 0.109481575485;
    fWeight[4]= 0.0226605;
    fLocationKsi[5]= 0.109481575485;
    fLocationEta[5]= 0.781036849030;
    fWeight[5]= 0.0226605;
    fLocationKsi[6]= 0.781036849030;
    fLocationEta[6]= 0.109481575485;
    fWeight[6]= 0.0226605;
    fLocationKsi[7]= 0.307939838764;
    fLocationEta[7]= 0.550352941821;
    fWeight[7]= 0.036379;
    fLocationKsi[8]= 0.550352941821;
    fLocationEta[8]= 0.141707219415;
    fWeight[8]= 0.036379;
    fLocationKsi[9]= 0.141707219415;
    fLocationEta[9]= 0.307939838764;
    fWeight[9]= 0.036379;
    fLocationKsi[10]= 0.307939838764;
    fLocationEta[10]= 0.141707219415;
    fWeight[10]= 0.036379;
    fLocationKsi[11]= 0.550352941821;
    fLocationEta[11]= 0.307939838764;
    fWeight[11]= 0.036379;
    fLocationKsi[12]= 0.141707219415;
    fLocationEta[12]= 0.550352941821;
    fWeight[12]= 0.036379;
    fLocationKsi[13]= 0.246672560640;
    fLocationEta[13]= 0.728323904597;
    fWeight[13]= 0.0141636;
    fLocationKsi[14]= 0.728323904597;
    fLocationEta[14]= 0.025003534763;
    fWeight[14]= 0.0141636;
    fLocationKsi[15]= 0.025003534763;
    fLocationEta[15]= 0.246672560640;
    fWeight[15]= 0.0141636;
    fLocationKsi[16]= 0.246672560640;
    fLocationEta[16]= 0.025003534763;
    fWeight[16]= 0.0141636;
    fLocationKsi[17]= 0.728323904597;
    fLocationEta[17]= 0.246672560640;
    fWeight[17]= 0.0141636;
    fLocationKsi[18]= 0.025003534763;
    fLocationEta[18]= 0.728323904597;
    fWeight[18]= 0.0141636;
    fLocationKsi[19]= 0.066803251012;
    fLocationEta[19]= 0.923655933587;
    fWeight[19]= 0.00471083;
    fLocationKsi[20]= 0.923655933587;
    fLocationEta[20]= 0.009540815400;
    fWeight[20]= 0.00471083;
    fLocationKsi[21]= 0.009540815400;
    fLocationEta[21]= 0.066803251012;
    fWeight[21]= 0.00471083;
    fLocationKsi[22]= 0.066803251012;
    fLocationEta[22]= 0.009540815400;
    fWeight[22]= 0.00471083;
    fLocationKsi[23]= 0.923655933587;
    fLocationEta[23]= 0.066803251012;
    fWeight[23]= 0.00471083;
    fLocationKsi[24]= 0.009540815400;
    fLocationEta[24]= 0.923655933587;
    fWeight[24]= 0.00471083;
    break;

  case 11:
    fLocationKsi[0]= 0.534611048271;
    fLocationEta[0]= 0.534611048271;
    fWeight[0]= 0.000463503;
    fLocationKsi[1]= 0.534611048271;
    fLocationEta[1]= -0.069222096542;
    fWeight[1]= 0.000463503;
    fLocationKsi[2]= -0.069222096542;
    fLocationEta[2]= 0.534611048271;
    fWeight[2]= 0.000463503;
    fLocationKsi[3]= 0.398969302966;
    fLocationEta[3]= 0.398969302966;
    fWeight[3]= 0.0385748;
    fLocationKsi[4]= 0.398969302966;
    fLocationEta[4]= 0.202061394068;
    fWeight[4]= 0.0385748;
    fLocationKsi[5]= 0.202061394068;
    fLocationEta[5]= 0.398969302966;
    fWeight[5]= 0.0385748;
    fLocationKsi[6]= 0.203309900431;
    fLocationEta[6]= 0.203309900431;
    fWeight[6]= 0.0296615;
    fLocationKsi[7]= 0.203309900431;
    fLocationEta[7]= 0.593380199137;
    fWeight[7]= 0.0296615;
    fLocationKsi[8]= 0.593380199137;
    fLocationEta[8]= 0.203309900431;
    fWeight[8]= 0.0296615;
    fLocationKsi[9]= 0.119350912283;
    fLocationEta[9]= 0.119350912283;
    fWeight[9]= 0.0180923;
    fLocationKsi[10]= 0.119350912283;
    fLocationEta[10]= 0.761298175435;
    fWeight[10]= 0.0180923;
    fLocationKsi[11]= 0.761298175435;
    fLocationEta[11]= 0.119350912283;
    fWeight[11]= 0.0180923;
    fLocationKsi[12]= 0.032364948111;
    fLocationEta[12]= 0.032364948111;
    fWeight[12]= 0.00682987;
    fLocationKsi[13]= 0.032364948111;
    fLocationEta[13]= 0.935270103777;
    fWeight[13]= 0.00682987;
    fLocationKsi[14]= 0.935270103777;
    fLocationEta[14]= 0.032364948111;
    fWeight[14]= 0.00682987;
    fLocationKsi[15]= 0.356620648261;
    fLocationEta[15]= 0.593201213428;
    fWeight[15]= 0.0261686;
    fLocationKsi[16]= 0.593201213428;
    fLocationEta[16]= 0.050178138310;
    fWeight[16]= 0.0261686;
    fLocationKsi[17]= 0.050178138310;
    fLocationEta[17]= 0.356620648261;
    fWeight[17]= 0.0261686;
    fLocationKsi[18]= 0.356620648261;
    fLocationEta[18]= 0.050178138310;
    fWeight[18]= 0.0261686;
    fLocationKsi[19]= 0.593201213428;
    fLocationEta[19]= 0.356620648261;
    fWeight[19]= 0.0261686;
    fLocationKsi[20]= 0.050178138310;
    fLocationEta[20]= 0.593201213428;
    fWeight[20]= 0.0261686;
    fLocationKsi[21]= 0.171488980304;
    fLocationEta[21]= 0.807489003160;
    fWeight[21]= 0.0103538;
    fLocationKsi[22]= 0.807489003160;
    fLocationEta[22]= 0.021022016536;
    fWeight[22]= 0.0103538;
    fLocationKsi[23]= 0.021022016536;
    fLocationEta[23]= 0.171488980304;
    fWeight[23]= 0.0103538;
    fLocationKsi[24]= 0.171488980304;
    fLocationEta[24]= 0.021022016536;
    fWeight[24]= 0.0103538;
    fLocationKsi[25]= 0.807489003160;
    fLocationEta[25]= 0.171488980304;
    fWeight[25]= 0.0103538;
    fLocationKsi[26]= 0.021022016536;
    fLocationEta[26]= 0.807489003160;
    fWeight[26]= 0.0103538;
    break;

  case 12:
    fLocationKsi[0]= 0.488217389774;
    fLocationEta[0]= 0.488217389774;
    fWeight[0]= 0.0128655;
    fLocationKsi[1]= 0.488217389774;
    fLocationEta[1]= 0.023565220452;
    fWeight[1]= 0.0128655;
    fLocationKsi[2]= 0.023565220452;
    fLocationEta[2]= 0.488217389774;
    fWeight[2]= 0.0128655;
    fLocationKsi[3]= 0.439724392294;
    fLocationEta[3]= 0.439724392294;
    fWeight[3]= 0.0218463;
    fLocationKsi[4]= 0.439724392294;
    fLocationEta[4]= 0.120551215411;
    fWeight[4]= 0.0218463;
    fLocationKsi[5]= 0.120551215411;
    fLocationEta[5]= 0.439724392294;
    fWeight[5]= 0.0218463;
    fLocationKsi[6]= 0.271210385012;
    fLocationEta[6]= 0.271210385012;
    fWeight[6]= 0.0314291;
    fLocationKsi[7]= 0.271210385012;
    fLocationEta[7]= 0.457579229976;
    fWeight[7]= 0.0314291;
    fLocationKsi[8]= 0.457579229976;
    fLocationEta[8]= 0.271210385012;
    fWeight[8]= 0.0314291;
    fLocationKsi[9]= 0.127576145542;
    fLocationEta[9]= 0.127576145542;
    fWeight[9]= 0.0173981;
    fLocationKsi[10]= 0.127576145542;
    fLocationEta[10]= 0.744847708917;
    fWeight[10]= 0.0173981;
    fLocationKsi[11]= 0.744847708917;
    fLocationEta[11]= 0.127576145542;
    fWeight[11]= 0.0173981;
    fLocationKsi[12]= 0.021317350453;
    fLocationEta[12]= 0.021317350453;
    fWeight[12]= 0.00308313;
    fLocationKsi[13]= 0.021317350453;
    fLocationEta[13]= 0.957365299094;
    fWeight[13]= 0.00308313;
    fLocationKsi[14]= 0.957365299094;
    fLocationEta[14]= 0.021317350453;
    fWeight[14]= 0.00308313;
    fLocationKsi[15]= 0.275713269686;
    fLocationEta[15]= 0.608943235780;
    fWeight[15]= 0.0201858;
    fLocationKsi[16]= 0.608943235780;
    fLocationEta[16]= 0.115343494535;
    fWeight[16]= 0.0201858;
    fLocationKsi[17]= 0.115343494535;
    fLocationEta[17]= 0.275713269686;
    fWeight[17]= 0.0201858;
    fLocationKsi[18]= 0.275713269686;
    fLocationEta[18]= 0.115343494535;
    fWeight[18]= 0.0201858;
    fLocationKsi[19]= 0.608943235780;
    fLocationEta[19]= 0.275713269686;
    fWeight[19]= 0.0201858;
    fLocationKsi[20]= 0.115343494535;
    fLocationEta[20]= 0.608943235780;
    fWeight[20]= 0.0201858;
    fLocationKsi[21]= 0.281325580990;
    fLocationEta[21]= 0.695836086788;
    fWeight[21]= 0.0111784;
    fLocationKsi[22]= 0.695836086788;
    fLocationEta[22]= 0.022838332222;
    fWeight[22]= 0.0111784;
    fLocationKsi[23]= 0.022838332222;
    fLocationEta[23]= 0.281325580990;
    fWeight[23]= 0.0111784;
    fLocationKsi[24]= 0.281325580990;
    fLocationEta[24]= 0.022838332222;
    fWeight[24]= 0.0111784;
    fLocationKsi[25]= 0.695836086788;
    fLocationEta[25]= 0.281325580990;
    fWeight[25]= 0.0111784;
    fLocationKsi[26]= 0.022838332222;
    fLocationEta[26]= 0.695836086788;
    fWeight[26]= 0.0111784;
    fLocationKsi[27]= 0.116251915908;
    fLocationEta[27]= 0.858014033544;
    fWeight[27]= 0.00865812;
    fLocationKsi[28]= 0.858014033544;
    fLocationEta[28]= 0.025734050548;
    fWeight[28]= 0.00865812;
    fLocationKsi[29]= 0.025734050548;
    fLocationEta[29]= 0.116251915908;
    fWeight[29]= 0.00865812;
    fLocationKsi[30]= 0.116251915908;
    fLocationEta[30]= 0.025734050548;
    fWeight[30]= 0.00865812;
    fLocationKsi[31]= 0.858014033544;
    fLocationEta[31]= 0.116251915908;
    fWeight[31]= 0.00865812;
    fLocationKsi[32]= 0.025734050548;
    fLocationEta[32]= 0.858014033544;
    fWeight[32]= 0.00865812;
    break;

  case 13:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.0262605;
    fLocationKsi[1]= 0.495048184940;
    fLocationEta[1]= 0.495048184940;
    fWeight[1]= 0.00564007;
    fLocationKsi[2]= 0.495048184940;
    fLocationEta[2]= 0.009903630121;
    fWeight[2]= 0.00564007;
    fLocationKsi[3]= 0.009903630121;
    fLocationEta[3]= 0.495048184940;
    fWeight[3]= 0.00564007;
    fLocationKsi[4]= 0.468716635110;
    fLocationEta[4]= 0.468716635110;
    fWeight[4]= 0.0157118;
    fLocationKsi[5]= 0.468716635110;
    fLocationEta[5]= 0.062566729781;
    fWeight[5]= 0.0157118;
    fLocationKsi[6]= 0.062566729781;
    fLocationEta[6]= 0.468716635110;
    fWeight[6]= 0.0157118;
    fLocationKsi[7]= 0.414521336801;
    fLocationEta[7]= 0.414521336801;
    fWeight[7]= 0.0235363;
    fLocationKsi[8]= 0.414521336801;
    fLocationEta[8]= 0.170957326397;
    fWeight[8]= 0.0235363;
    fLocationKsi[9]= 0.170957326397;
    fLocationEta[9]= 0.414521336801;
    fWeight[9]= 0.0235363;
    fLocationKsi[10]= 0.229399572043;
    fLocationEta[10]= 0.229399572043;
    fWeight[10]= 0.0236818;
    fLocationKsi[11]= 0.229399572043;
    fLocationEta[11]= 0.541200855914;
    fWeight[11]= 0.0236818;
    fLocationKsi[12]= 0.541200855914;
    fLocationEta[12]= 0.229399572043;
    fWeight[12]= 0.0236818;
    fLocationKsi[13]= 0.114424495196;
    fLocationEta[13]= 0.114424495196;
    fWeight[13]= 0.0155838;
    fLocationKsi[14]= 0.114424495196;
    fLocationEta[14]= 0.771151009607;
    fWeight[14]= 0.0155838;
    fLocationKsi[15]= 0.771151009607;
    fLocationEta[15]= 0.114424495196;
    fWeight[15]= 0.0155838;
    fLocationKsi[16]= 0.024811391363;
    fLocationEta[16]= 0.024811391363;
    fWeight[16]= 0.00398789;
    fLocationKsi[17]= 0.024811391363;
    fLocationEta[17]= 0.950377217273;
    fWeight[17]= 0.00398789;
    fLocationKsi[18]= 0.950377217273;
    fLocationEta[18]= 0.024811391363;
    fWeight[18]= 0.00398789;
    fLocationKsi[19]= 0.268794997059;
    fLocationEta[19]= 0.636351174562;
    fWeight[19]= 0.0184242;
    fLocationKsi[20]= 0.636351174562;
    fLocationEta[20]= 0.094853828380;
    fWeight[20]= 0.0184242;
    fLocationKsi[21]= 0.094853828380;
    fLocationEta[21]= 0.268794997059;
    fWeight[21]= 0.0184242;
    fLocationKsi[22]= 0.268794997059;
    fLocationEta[22]= 0.094853828380;
    fWeight[22]= 0.0184242;
    fLocationKsi[23]= 0.636351174562;
    fLocationEta[23]= 0.268794997059;
    fWeight[23]= 0.0184242;
    fLocationKsi[24]= 0.094853828380;
    fLocationEta[24]= 0.636351174562;
    fWeight[24]= 0.0184242;
    fLocationKsi[25]= 0.291730066734;
    fLocationEta[25]= 0.690169159987;
    fWeight[25]= 0.00870073;
    fLocationKsi[26]= 0.690169159987;
    fLocationEta[26]= 0.018100773279;
    fWeight[26]= 0.00870073;
    fLocationKsi[27]= 0.018100773279;
    fLocationEta[27]= 0.291730066734;
    fWeight[27]= 0.00870073;
    fLocationKsi[28]= 0.291730066734;
    fLocationEta[28]= 0.018100773279;
    fWeight[28]= 0.00870073;
    fLocationKsi[29]= 0.690169159987;
    fLocationEta[29]= 0.291730066734;
    fWeight[29]= 0.00870073;
    fLocationKsi[30]= 0.018100773279;
    fLocationEta[30]= 0.690169159987;
    fWeight[30]= 0.00870073;
    fLocationKsi[31]= 0.126357385492;
    fLocationEta[31]= 0.851409537834;
    fWeight[31]= 0.00776089;
    fLocationKsi[32]= 0.851409537834;
    fLocationEta[32]= 0.022233076674;
    fWeight[32]= 0.00776089;
    fLocationKsi[33]= 0.022233076674;
    fLocationEta[33]= 0.126357385492;
    fWeight[33]= 0.00776089;
    fLocationKsi[34]= 0.126357385492;
    fLocationEta[34]= 0.022233076674;
    fWeight[34]= 0.00776089;
    fLocationKsi[35]= 0.851409537834;
    fLocationEta[35]= 0.126357385492;
    fWeight[35]= 0.00776089;
    fLocationKsi[36]= 0.022233076674;
    fLocationEta[36]= 0.851409537834;
    fWeight[36]= 0.00776089;
    break;

  case 14:
    fLocationKsi[0]= 0.488963910362;
    fLocationEta[0]= 0.488963910362;
    fWeight[0]= 0.0109418;
    fLocationKsi[1]= 0.488963910362;
    fLocationEta[1]= 0.022072179276;
    fWeight[1]= 0.0109418;
    fLocationKsi[2]= 0.022072179276;
    fLocationEta[2]= 0.488963910362;
    fWeight[2]= 0.0109418;
    fLocationKsi[3]= 0.417644719340;
    fLocationEta[3]= 0.417644719340;
    fWeight[3]= 0.0163942;
    fLocationKsi[4]= 0.417644719340;
    fLocationEta[4]= 0.164710561319;
    fWeight[4]= 0.0163942;
    fLocationKsi[5]= 0.164710561319;
    fLocationEta[5]= 0.417644719340;
    fWeight[5]= 0.0163942;
    fLocationKsi[6]= 0.273477528309;
    fLocationEta[6]= 0.273477528309;
    fWeight[6]= 0.0258871;
    fLocationKsi[7]= 0.273477528309;
    fLocationEta[7]= 0.453044943382;
    fWeight[7]= 0.0258871;
    fLocationKsi[8]= 0.453044943382;
    fLocationEta[8]= 0.273477528309;
    fWeight[8]= 0.0258871;
    fLocationKsi[9]= 0.177205532413;
    fLocationEta[9]= 0.177205532413;
    fWeight[9]= 0.0210813;
    fLocationKsi[10]= 0.177205532413;
    fLocationEta[10]= 0.645588935175;
    fWeight[10]= 0.0210813;
    fLocationKsi[11]= 0.645588935175;
    fLocationEta[11]= 0.177205532413;
    fWeight[11]= 0.0210813;
    fLocationKsi[12]= 0.061799883091;
    fLocationEta[12]= 0.061799883091;
    fWeight[12]= 0.00721685;
    fLocationKsi[13]= 0.061799883091;
    fLocationEta[13]= 0.876400233818;
    fWeight[13]= 0.00721685;
    fLocationKsi[14]= 0.876400233818;
    fLocationEta[14]= 0.061799883091;
    fWeight[14]= 0.00721685;
    fLocationKsi[15]= 0.019390961249;
    fLocationEta[15]= 0.019390961249;
    fWeight[15]= 0.0024617;
    fLocationKsi[16]= 0.019390961249;
    fLocationEta[16]= 0.961218077503;
    fWeight[16]= 0.0024617;
    fLocationKsi[17]= 0.961218077503;
    fLocationEta[17]= 0.019390961249;
    fWeight[17]= 0.0024617;
    fLocationKsi[18]= 0.172266687821;
    fLocationEta[18]= 0.770608554775;
    fWeight[18]= 0.0123329;
    fLocationKsi[19]= 0.770608554775;
    fLocationEta[19]= 0.057124757404;
    fWeight[19]= 0.0123329;
    fLocationKsi[20]= 0.057124757404;
    fLocationEta[20]= 0.172266687821;
    fWeight[20]= 0.0123329;
    fLocationKsi[21]= 0.172266687821;
    fLocationEta[21]= 0.057124757404;
    fWeight[21]= 0.0123329;
    fLocationKsi[22]= 0.770608554775;
    fLocationEta[22]= 0.172266687821;
    fWeight[22]= 0.0123329;
    fLocationKsi[23]= 0.057124757404;
    fLocationEta[23]= 0.770608554775;
    fWeight[23]= 0.0123329;
    fLocationKsi[24]= 0.336861459796;
    fLocationEta[24]= 0.570222290847;
    fWeight[24]= 0.0192858;
    fLocationKsi[25]= 0.570222290847;
    fLocationEta[25]= 0.092916249357;
    fWeight[25]= 0.0192858;
    fLocationKsi[26]= 0.092916249357;
    fLocationEta[26]= 0.336861459796;
    fWeight[26]= 0.0192858;
    fLocationKsi[27]= 0.336861459796;
    fLocationEta[27]= 0.092916249357;
    fWeight[27]= 0.0192858;
    fLocationKsi[28]= 0.570222290847;
    fLocationEta[28]= 0.336861459796;
    fWeight[28]= 0.0192858;
    fLocationKsi[29]= 0.092916249357;
    fLocationEta[29]= 0.570222290847;
    fWeight[29]= 0.0192858;
    fLocationKsi[30]= 0.298372882136;
    fLocationEta[30]= 0.686980167808;
    fWeight[30]= 0.00721815;
    fLocationKsi[31]= 0.686980167808;
    fLocationEta[31]= 0.014646950056;
    fWeight[31]= 0.00721815;
    fLocationKsi[32]= 0.014646950056;
    fLocationEta[32]= 0.298372882136;
    fWeight[32]= 0.00721815;
    fLocationKsi[33]= 0.298372882136;
    fLocationEta[33]= 0.014646950056;
    fWeight[33]= 0.00721815;
    fLocationKsi[34]= 0.686980167808;
    fLocationEta[34]= 0.298372882136;
    fWeight[34]= 0.00721815;
    fLocationKsi[35]= 0.014646950056;
    fLocationEta[35]= 0.686980167808;
    fWeight[35]= 0.00721815;
    fLocationKsi[36]= 0.118974497697;
    fLocationEta[36]= 0.879757171370;
    fWeight[36]= 0.00250511;
    fLocationKsi[37]= 0.879757171370;
    fLocationEta[37]= 0.001268330933;
    fWeight[37]= 0.00250511;
    fLocationKsi[38]= 0.001268330933;
    fLocationEta[38]= 0.118974497697;
    fWeight[38]= 0.00250511;
    fLocationKsi[39]= 0.118974497697;
    fLocationEta[39]= 0.001268330933;
    fWeight[39]= 0.00250511;
    fLocationKsi[40]= 0.879757171370;
    fLocationEta[40]= 0.118974497697;
    fWeight[40]= 0.00250511;
    fLocationKsi[41]= 0.001268330933;
    fLocationEta[41]= 0.879757171370;
    fWeight[41]= 0.00250511;
    break;

  case 15:
    fLocationKsi[0]= 0.506972916858;
    fLocationEta[0]= 0.506972916858;
    fWeight[0]= 0.000958438;
    fLocationKsi[1]= 0.506972916858;
    fLocationEta[1]= -0.013945833716;
    fWeight[1]= 0.000958438;
    fLocationKsi[2]= -0.013945833716;
    fLocationEta[2]= 0.506972916858;
    fWeight[2]= 0.000958438;
    fLocationKsi[3]= 0.431406354283;
    fLocationEta[3]= 0.431406354283;
    fWeight[3]= 0.0221245;
    fLocationKsi[4]= 0.431406354283;
    fLocationEta[4]= 0.137187291434;
    fWeight[4]= 0.0221245;
    fLocationKsi[5]= 0.137187291434;
    fLocationEta[5]= 0.431406354283;
    fWeight[5]= 0.0221245;
    fLocationKsi[6]= 0.277693644847;
    fLocationEta[6]= 0.277693644847;
    fWeight[6]= 0.0255933;
    fLocationKsi[7]= 0.277693644847;
    fLocationEta[7]= 0.444612710306;
    fWeight[7]= 0.0255933;
    fLocationKsi[8]= 0.444612710306;
    fLocationEta[8]= 0.277693644847;
    fWeight[8]= 0.0255933;
    fLocationKsi[9]= 0.126464891041;
    fLocationEta[9]= 0.126464891041;
    fWeight[9]= 0.0118439;
    fLocationKsi[10]= 0.126464891041;
    fLocationEta[10]= 0.747070217917;
    fWeight[10]= 0.0118439;
    fLocationKsi[11]= 0.747070217917;
    fLocationEta[11]= 0.126464891041;
    fWeight[11]= 0.0118439;
    fLocationKsi[12]= 0.070808385975;
    fLocationEta[12]= 0.070808385975;
    fWeight[12]= 0.00664489;
    fLocationKsi[13]= 0.070808385975;
    fLocationEta[13]= 0.858383228051;
    fWeight[13]= 0.00664489;
    fLocationKsi[14]= 0.858383228051;
    fLocationEta[14]= 0.070808385975;
    fWeight[14]= 0.00664489;
    fLocationKsi[15]= 0.018965170241;
    fLocationEta[15]= 0.018965170241;
    fWeight[15]= 0.00237446;
    fLocationKsi[16]= 0.018965170241;
    fLocationEta[16]= 0.962069659518;
    fWeight[16]= 0.00237446;
    fLocationKsi[17]= 0.962069659518;
    fLocationEta[17]= 0.018965170241;
    fWeight[17]= 0.00237446;
    fLocationKsi[18]= 0.261311371140;
    fLocationEta[18]= 0.604954466893;
    fWeight[18]= 0.019275;
    fLocationKsi[19]= 0.604954466893;
    fLocationEta[19]= 0.133734161967;
    fWeight[19]= 0.019275;
    fLocationKsi[20]= 0.133734161967;
    fLocationEta[20]= 0.261311371140;
    fWeight[20]= 0.019275;
    fLocationKsi[21]= 0.261311371140;
    fLocationEta[21]= 0.133734161967;
    fWeight[21]= 0.019275;
    fLocationKsi[22]= 0.604954466893;
    fLocationEta[22]= 0.261311371140;
    fWeight[22]= 0.019275;
    fLocationKsi[23]= 0.133734161967;
    fLocationEta[23]= 0.604954466893;
    fWeight[23]= 0.019275;
    fLocationKsi[24]= 0.388046767090;
    fLocationEta[24]= 0.575586555513;
    fWeight[24]= 0.0136079;
    fLocationKsi[25]= 0.575586555513;
    fLocationEta[25]= 0.036366677397;
    fWeight[25]= 0.0136079;
    fLocationKsi[26]= 0.036366677397;
    fLocationEta[26]= 0.388046767090;
    fWeight[26]= 0.0136079;
    fLocationKsi[27]= 0.388046767090;
    fLocationEta[27]= 0.036366677397;
    fWeight[27]= 0.0136079;
    fLocationKsi[28]= 0.575586555513;
    fLocationEta[28]= 0.388046767090;
    fWeight[28]= 0.0136079;
    fLocationKsi[29]= 0.036366677397;
    fLocationEta[29]= 0.575586555513;
    fWeight[29]= 0.0136079;
    fLocationKsi[30]= 0.285712220050;
    fLocationEta[30]= 0.724462663077;
    fWeight[30]= 0.00109104;
    fLocationKsi[31]= 0.724462663077;
    fLocationEta[31]= -0.010174883127;
    fWeight[31]= 0.00109104;
    fLocationKsi[32]= -0.010174883127;
    fLocationEta[32]= 0.285712220050;
    fWeight[32]= 0.00109104;
    fLocationKsi[33]= 0.285712220050;
    fLocationEta[33]= -0.010174883127;
    fWeight[33]= 0.00109104;
    fLocationKsi[34]= 0.724462663077;
    fLocationEta[34]= 0.285712220050;
    fWeight[34]= 0.00109104;
    fLocationKsi[35]= -0.010174883127;
    fLocationEta[35]= 0.724462663077;
    fWeight[35]= 0.00109104;
    fLocationKsi[36]= 0.215599664072;
    fLocationEta[36]= 0.747556466052;
    fWeight[36]= 0.0107527;
    fLocationKsi[37]= 0.747556466052;
    fLocationEta[37]= 0.036843869876;
    fWeight[37]= 0.0107527;
    fLocationKsi[38]= 0.036843869876;
    fLocationEta[38]= 0.215599664072;
    fWeight[38]= 0.0107527;
    fLocationKsi[39]= 0.215599664072;
    fLocationEta[39]= 0.036843869876;
    fWeight[39]= 0.0107527;
    fLocationKsi[40]= 0.747556466052;
    fLocationEta[40]= 0.215599664072;
    fWeight[40]= 0.0107527;
    fLocationKsi[41]= 0.036843869876;
    fLocationEta[41]= 0.747556466052;
    fWeight[41]= 0.0107527;
    fLocationKsi[42]= 0.103575616576;
    fLocationEta[42]= 0.883964574092;
    fWeight[42]= 0.00383697;
    fLocationKsi[43]= 0.883964574092;
    fLocationEta[43]= 0.012459809331;
    fWeight[43]= 0.00383697;
    fLocationKsi[44]= 0.012459809331;
    fLocationEta[44]= 0.103575616576;
    fWeight[44]= 0.00383697;
    fLocationKsi[45]= 0.103575616576;
    fLocationEta[45]= 0.012459809331;
    fWeight[45]= 0.00383697;
    fLocationKsi[46]= 0.883964574092;
    fLocationEta[46]= 0.103575616576;
    fWeight[46]= 0.00383697;
    fLocationKsi[47]= 0.012459809331;
    fLocationEta[47]= 0.883964574092;
    fWeight[47]= 0.00383697;
    break;

  case 16:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.0234378;
    fLocationKsi[1]= 0.497380541948;
    fLocationEta[1]= 0.497380541948;
    fWeight[1]= 0.00320294;
    fLocationKsi[2]= 0.497380541948;
    fLocationEta[2]= 0.005238916103;
    fWeight[2]= 0.00320294;
    fLocationKsi[3]= 0.005238916103;
    fLocationEta[3]= 0.497380541948;
    fWeight[3]= 0.00320294;
    fLocationKsi[4]= 0.413469438549;
    fLocationEta[4]= 0.413469438549;
    fWeight[4]= 0.0208551;
    fLocationKsi[5]= 0.413469438549;
    fLocationEta[5]= 0.173061122901;
    fWeight[5]= 0.0208551;
    fLocationKsi[6]= 0.173061122901;
    fLocationEta[6]= 0.413469438549;
    fWeight[6]= 0.0208551;
    fLocationKsi[7]= 0.470458599067;
    fLocationEta[7]= 0.470458599067;
    fWeight[7]= 0.0134457;
    fLocationKsi[8]= 0.470458599067;
    fLocationEta[8]= 0.059082801866;
    fWeight[8]= 0.0134457;
    fLocationKsi[9]= 0.059082801866;
    fLocationEta[9]= 0.470458599067;
    fWeight[9]= 0.0134457;
    fLocationKsi[10]= 0.240553749970;
    fLocationEta[10]= 0.240553749970;
    fWeight[10]= 0.0210663;
    fLocationKsi[11]= 0.240553749970;
    fLocationEta[11]= 0.518892500061;
    fWeight[11]= 0.0210663;
    fLocationKsi[12]= 0.518892500061;
    fLocationEta[12]= 0.240553749970;
    fWeight[12]= 0.0210663;
    fLocationKsi[13]= 0.147965794223;
    fLocationEta[13]= 0.147965794223;
    fWeight[13]= 0.0150001;
    fLocationKsi[14]= 0.147965794223;
    fLocationEta[14]= 0.704068411555;
    fWeight[14]= 0.0150001;
    fLocationKsi[15]= 0.704068411555;
    fLocationEta[15]= 0.147965794223;
    fWeight[15]= 0.0150001;
    fLocationKsi[16]= 0.075465187657;
    fLocationEta[16]= 0.075465187657;
    fWeight[16]= 0.00710005;
    fLocationKsi[17]= 0.075465187657;
    fLocationEta[17]= 0.849069624685;
    fWeight[17]= 0.00710005;
    fLocationKsi[18]= 0.849069624685;
    fLocationEta[18]= 0.075465187657;
    fWeight[18]= 0.00710005;
    fLocationKsi[19]= 0.016596402623;
    fLocationEta[19]= 0.016596402623;
    fWeight[19]= 0.00179123;
    fLocationKsi[20]= 0.016596402623;
    fLocationEta[20]= 0.966807194754;
    fWeight[20]= 0.00179123;
    fLocationKsi[21]= 0.966807194754;
    fLocationEta[21]= 0.016596402623;
    fWeight[21]= 0.00179123;
    fLocationKsi[22]= 0.296555596580;
    fLocationEta[22]= 0.599868711175;
    fWeight[22]= 0.0163866;
    fLocationKsi[23]= 0.599868711175;
    fLocationEta[23]= 0.103575692245;
    fWeight[23]= 0.0163866;
    fLocationKsi[24]= 0.103575692245;
    fLocationEta[24]= 0.296555596580;
    fWeight[24]= 0.0163866;
    fLocationKsi[25]= 0.296555596580;
    fLocationEta[25]= 0.103575692245;
    fWeight[25]= 0.0163866;
    fLocationKsi[26]= 0.599868711175;
    fLocationEta[26]= 0.296555596580;
    fWeight[26]= 0.0163866;
    fLocationKsi[27]= 0.103575692245;
    fLocationEta[27]= 0.599868711175;
    fWeight[27]= 0.0163866;
    fLocationKsi[28]= 0.337723063403;
    fLocationEta[28]= 0.642193524942;
    fWeight[28]= 0.00764915;
    fLocationKsi[29]= 0.642193524942;
    fLocationEta[29]= 0.020083411655;
    fWeight[29]= 0.00764915;
    fLocationKsi[30]= 0.020083411655;
    fLocationEta[30]= 0.337723063403;
    fWeight[30]= 0.00764915;
    fLocationKsi[31]= 0.337723063403;
    fLocationEta[31]= 0.020083411655;
    fWeight[31]= 0.00764915;
    fLocationKsi[32]= 0.642193524942;
    fLocationEta[32]= 0.337723063403;
    fWeight[32]= 0.00764915;
    fLocationKsi[33]= 0.020083411655;
    fLocationEta[33]= 0.642193524942;
    fWeight[33]= 0.00764915;
    fLocationKsi[34]= 0.204748281643;
    fLocationEta[34]= 0.799592720971;
    fWeight[34]= 0.00119312;
    fLocationKsi[35]= 0.799592720971;
    fLocationEta[35]= -0.004341002614;
    fWeight[35]= 0.00119312;
    fLocationKsi[36]= -0.004341002614;
    fLocationEta[36]= 0.204748281643;
    fWeight[36]= 0.00119312;
    fLocationKsi[37]= 0.204748281643;
    fLocationEta[37]= -0.004341002614;
    fWeight[37]= 0.00119312;
    fLocationKsi[38]= 0.799592720971;
    fLocationEta[38]= 0.204748281643;
    fWeight[38]= 0.00119312;
    fLocationKsi[39]= -0.004341002614;
    fLocationEta[39]= 0.799592720971;
    fWeight[39]= 0.00119312;
    fLocationKsi[40]= 0.189358492131;
    fLocationEta[40]= 0.768699721401;
    fWeight[40]= 0.0095424;
    fLocationKsi[41]= 0.768699721401;
    fLocationEta[41]= 0.041941786468;
    fWeight[41]= 0.0095424;
    fLocationKsi[42]= 0.041941786468;
    fLocationEta[42]= 0.189358492131;
    fWeight[42]= 0.0095424;
    fLocationKsi[43]= 0.189358492131;
    fLocationEta[43]= 0.041941786468;
    fWeight[43]= 0.0095424;
    fLocationKsi[44]= 0.768699721401;
    fLocationEta[44]= 0.189358492131;
    fWeight[44]= 0.0095424;
    fLocationKsi[45]= 0.041941786468;
    fLocationEta[45]= 0.768699721401;
    fWeight[45]= 0.0095424;
    fLocationKsi[46]= 0.085283615683;
    fLocationEta[46]= 0.900399064087;
    fWeight[46]= 0.00342503;
    fLocationKsi[47]= 0.900399064087;
    fLocationEta[47]= 0.014317320231;
    fWeight[47]= 0.00342503;
    fLocationKsi[48]= 0.014317320231;
    fLocationEta[48]= 0.085283615683;
    fWeight[48]= 0.00342503;
    fLocationKsi[49]= 0.085283615683;
    fLocationEta[49]= 0.014317320231;
    fWeight[49]= 0.00342503;
    fLocationKsi[50]= 0.900399064087;
    fLocationEta[50]= 0.085283615683;
    fWeight[50]= 0.00342503;
    fLocationKsi[51]= 0.014317320231;
    fLocationEta[51]= 0.900399064087;
    fWeight[51]= 0.00342503;
    break;

  case 17:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.0167186;
    fLocationKsi[1]= 0.497170540557;
    fLocationEta[1]= 0.497170540557;
    fWeight[1]= 0.00254671;
    fLocationKsi[2]= 0.497170540557;
    fLocationEta[2]= 0.005658918886;
    fWeight[2]= 0.00254671;
    fLocationKsi[3]= 0.005658918886;
    fLocationEta[3]= 0.497170540557;
    fWeight[3]= 0.00254671;
    fLocationKsi[4]= 0.482176322625;
    fLocationEta[4]= 0.482176322625;
    fWeight[4]= 0.00733543;
    fLocationKsi[5]= 0.482176322625;
    fLocationEta[5]= 0.035647354751;
    fWeight[5]= 0.00733543;
    fLocationKsi[6]= 0.035647354751;
    fLocationEta[6]= 0.482176322625;
    fWeight[6]= 0.00733543;
    fLocationKsi[7]= 0.450239969021;
    fLocationEta[7]= 0.450239969021;
    fWeight[7]= 0.0121754;
    fLocationKsi[8]= 0.450239969021;
    fLocationEta[8]= 0.099520061958;
    fWeight[8]= 0.0121754;
    fLocationKsi[9]= 0.099520061958;
    fLocationEta[9]= 0.450239969021;
    fWeight[9]= 0.0121754;
    fLocationKsi[10]= 0.400266239377;
    fLocationEta[10]= 0.400266239377;
    fWeight[10]= 0.0155538;
    fLocationKsi[11]= 0.400266239377;
    fLocationEta[11]= 0.199467521245;
    fWeight[11]= 0.0155538;
    fLocationKsi[12]= 0.199467521245;
    fLocationEta[12]= 0.400266239377;
    fWeight[12]= 0.0155538;
    fLocationKsi[13]= 0.252141267971;
    fLocationEta[13]= 0.252141267971;
    fWeight[13]= 0.0156286;
    fLocationKsi[14]= 0.252141267971;
    fLocationEta[14]= 0.495717464058;
    fWeight[14]= 0.0156286;
    fLocationKsi[15]= 0.495717464058;
    fLocationEta[15]= 0.252141267971;
    fWeight[15]= 0.0156286;
    fLocationKsi[16]= 0.162047004658;
    fLocationEta[16]= 0.162047004658;
    fWeight[16]= 0.0124078;
    fLocationKsi[17]= 0.162047004658;
    fLocationEta[17]= 0.675905990683;
    fWeight[17]= 0.0124078;
    fLocationKsi[18]= 0.675905990683;
    fLocationEta[18]= 0.162047004658;
    fWeight[18]= 0.0124078;
    fLocationKsi[19]= 0.075875882261;
    fLocationEta[19]= 0.075875882261;
    fWeight[19]= 0.00702804;
    fLocationKsi[20]= 0.075875882261;
    fLocationEta[20]= 0.848248235479;
    fWeight[20]= 0.00702804;
    fLocationKsi[21]= 0.848248235479;
    fLocationEta[21]= 0.075875882261;
    fWeight[21]= 0.00702804;
    fLocationKsi[22]= 0.015654726968;
    fLocationEta[22]= 0.015654726968;
    fWeight[22]= 0.00159734;
    fLocationKsi[23]= 0.015654726968;
    fLocationEta[23]= 0.968690546064;
    fWeight[23]= 0.00159734;
    fLocationKsi[24]= 0.968690546064;
    fLocationEta[24]= 0.015654726968;
    fWeight[24]= 0.00159734;
    fLocationKsi[25]= 0.334319867364;
    fLocationEta[25]= 0.655493203809;
    fWeight[25]= 0.00405983;
    fLocationKsi[26]= 0.655493203809;
    fLocationEta[26]= 0.010186928827;
    fWeight[26]= 0.00405983;
    fLocationKsi[27]= 0.010186928827;
    fLocationEta[27]= 0.334319867364;
    fWeight[27]= 0.00405983;
    fLocationKsi[28]= 0.334319867364;
    fLocationEta[28]= 0.010186928827;
    fWeight[28]= 0.00405983;
    fLocationKsi[29]= 0.655493203809;
    fLocationEta[29]= 0.334319867364;
    fWeight[29]= 0.00405983;
    fLocationKsi[30]= 0.010186928827;
    fLocationEta[30]= 0.655493203809;
    fWeight[30]= 0.00405983;
    fLocationKsi[31]= 0.292221537797;
    fLocationEta[31]= 0.572337590532;
    fWeight[31]= 0.0134029;
    fLocationKsi[32]= 0.572337590532;
    fLocationEta[32]= 0.135440871671;
    fWeight[32]= 0.0134029;
    fLocationKsi[33]= 0.135440871671;
    fLocationEta[33]= 0.292221537797;
    fWeight[33]= 0.0134029;
    fLocationKsi[34]= 0.292221537797;
    fLocationEta[34]= 0.135440871671;
    fWeight[34]= 0.0134029;
    fLocationKsi[35]= 0.572337590532;
    fLocationEta[35]= 0.292221537797;
    fWeight[35]= 0.0134029;
    fLocationKsi[36]= 0.135440871671;
    fLocationEta[36]= 0.572337590532;
    fWeight[36]= 0.0134029;
    fLocationKsi[37]= 0.319574885423;
    fLocationEta[37]= 0.626001190286;
    fWeight[37]= 0.00923;
    fLocationKsi[38]= 0.626001190286;
    fLocationEta[38]= 0.054423924291;
    fWeight[38]= 0.00923;
    fLocationKsi[39]= 0.054423924291;
    fLocationEta[39]= 0.319574885423;
    fWeight[39]= 0.00923;
    fLocationKsi[40]= 0.319574885423;
    fLocationEta[40]= 0.054423924291;
    fWeight[40]= 0.00923;
    fLocationKsi[41]= 0.626001190286;
    fLocationEta[41]= 0.319574885423;
    fWeight[41]= 0.00923;
    fLocationKsi[42]= 0.054423924291;
    fLocationEta[42]= 0.626001190286;
    fWeight[42]= 0.00923;
    fLocationKsi[43]= 0.190704224192;
    fLocationEta[43]= 0.796427214974;
    fWeight[43]= 0.00423843;
    fLocationKsi[44]= 0.796427214974;
    fLocationEta[44]= 0.012868560834;
    fWeight[44]= 0.00423843;
    fLocationKsi[45]= 0.012868560834;
    fLocationEta[45]= 0.190704224192;
    fWeight[45]= 0.00423843;
    fLocationKsi[46]= 0.190704224192;
    fLocationEta[46]= 0.012868560834;
    fWeight[46]= 0.00423843;
    fLocationKsi[47]= 0.796427214974;
    fLocationEta[47]= 0.190704224192;
    fWeight[47]= 0.00423843;
    fLocationKsi[48]= 0.012868560834;
    fLocationEta[48]= 0.796427214974;
    fWeight[48]= 0.00423843;
    fLocationKsi[49]= 0.180483211649;
    fLocationEta[49]= 0.752351005938;
    fWeight[49]= 0.0091464;
    fLocationKsi[50]= 0.752351005938;
    fLocationEta[50]= 0.067165782414;
    fWeight[50]= 0.0091464;
    fLocationKsi[51]= 0.067165782414;
    fLocationEta[51]= 0.180483211649;
    fWeight[51]= 0.0091464;
    fLocationKsi[52]= 0.180483211649;
    fLocationEta[52]= 0.067165782414;
    fWeight[52]= 0.0091464;
    fLocationKsi[53]= 0.752351005938;
    fLocationEta[53]= 0.180483211649;
    fWeight[53]= 0.0091464;
    fLocationKsi[54]= 0.067165782414;
    fLocationEta[54]= 0.752351005938;
    fWeight[54]= 0.0091464;
    fLocationKsi[55]= 0.080711313680;
    fLocationEta[55]= 0.904625504096;
    fWeight[55]= 0.00333282;
    fLocationKsi[56]= 0.904625504096;
    fLocationEta[56]= 0.014663182225;
    fWeight[56]= 0.00333282;
    fLocationKsi[57]= 0.014663182225;
    fLocationEta[57]= 0.080711313680;
    fWeight[57]= 0.00333282;
    fLocationKsi[58]= 0.080711313680;
    fLocationEta[58]= 0.014663182225;
    fWeight[58]= 0.00333282;
    fLocationKsi[59]= 0.904625504096;
    fLocationEta[59]= 0.080711313680;
    fWeight[59]= 0.00333282;
    fLocationKsi[60]= 0.014663182225;
    fLocationEta[60]= 0.904625504096;
    fWeight[60]= 0.00333282;
    break;

  case 18:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.015405;
    fLocationKsi[1]= 0.493344808631;
    fLocationEta[1]= 0.493344808631;
    fWeight[1]= 0.00453622;
    fLocationKsi[2]= 0.493344808631;
    fLocationEta[2]= 0.013310382738;
    fWeight[2]= 0.00453622;
    fLocationKsi[3]= 0.013310382738;
    fLocationEta[3]= 0.493344808631;
    fWeight[3]= 0.00453622;
    fLocationKsi[4]= 0.469210594242;
    fLocationEta[4]= 0.469210594242;
    fWeight[4]= 0.00938066;
    fLocationKsi[5]= 0.469210594242;
    fLocationEta[5]= 0.061578811516;
    fWeight[5]= 0.00938066;
    fLocationKsi[6]= 0.061578811516;
    fLocationEta[6]= 0.469210594242;
    fWeight[6]= 0.00938066;
    fLocationKsi[7]= 0.436281395887;
    fLocationEta[7]= 0.436281395887;
    fWeight[7]= 0.00972055;
    fLocationKsi[8]= 0.436281395887;
    fLocationEta[8]= 0.127437208226;
    fWeight[8]= 0.00972055;
    fLocationKsi[9]= 0.127437208226;
    fLocationEta[9]= 0.436281395887;
    fWeight[9]= 0.00972055;
    fLocationKsi[10]= 0.394846170673;
    fLocationEta[10]= 0.394846170673;
    fWeight[10]= 0.013877;
    fLocationKsi[11]= 0.394846170673;
    fLocationEta[11]= 0.210307658653;
    fWeight[11]= 0.013877;
    fLocationKsi[12]= 0.210307658653;
    fLocationEta[12]= 0.394846170673;
    fWeight[12]= 0.013877;
    fLocationKsi[13]= 0.249794568803;
    fLocationEta[13]= 0.249794568803;
    fWeight[13]= 0.0161281;
    fLocationKsi[14]= 0.249794568803;
    fLocationEta[14]= 0.500410862394;
    fWeight[14]= 0.0161281;
    fLocationKsi[15]= 0.500410862394;
    fLocationEta[15]= 0.249794568803;
    fWeight[15]= 0.0161281;
    fLocationKsi[16]= 0.161432193744;
    fLocationEta[16]= 0.161432193744;
    fWeight[16]= 0.012537;
    fLocationKsi[17]= 0.161432193744;
    fLocationEta[17]= 0.677135612512;
    fWeight[17]= 0.012537;
    fLocationKsi[18]= 0.677135612512;
    fLocationEta[18]= 0.161432193744;
    fWeight[18]= 0.012537;
    fLocationKsi[19]= 0.076598227485;
    fLocationEta[19]= 0.076598227485;
    fWeight[19]= 0.00763596;
    fLocationKsi[20]= 0.076598227485;
    fLocationEta[20]= 0.846803545029;
    fWeight[20]= 0.00763596;
    fLocationKsi[21]= 0.846803545029;
    fLocationEta[21]= 0.076598227485;
    fWeight[21]= 0.00763596;
    fLocationKsi[22]= 0.024252439353;
    fLocationEta[22]= 0.024252439353;
    fWeight[22]= 0.00339696;
    fLocationKsi[23]= 0.024252439353;
    fLocationEta[23]= 0.951495121293;
    fWeight[23]= 0.00339696;
    fLocationKsi[24]= 0.951495121293;
    fLocationEta[24]= 0.024252439353;
    fWeight[24]= 0.00339696;
    fLocationKsi[25]= 0.043146367217;
    fLocationEta[25]= 0.043146367217;
    fWeight[25]= -0.00111155;
    fLocationKsi[26]= 0.043146367217;
    fLocationEta[26]= 0.913707265566;
    fWeight[26]= -0.00111155;
    fLocationKsi[27]= 0.913707265566;
    fLocationEta[27]= 0.043146367217;
    fWeight[27]= -0.00111155;
    fLocationKsi[28]= 0.358911494941;
    fLocationEta[28]= 0.632657968857;
    fWeight[28]= 0.00316596;
    fLocationKsi[29]= 0.632657968857;
    fLocationEta[29]= 0.008430536202;
    fWeight[29]= 0.00316596;
    fLocationKsi[30]= 0.008430536202;
    fLocationEta[30]= 0.358911494941;
    fWeight[30]= 0.00316596;
    fLocationKsi[31]= 0.358911494941;
    fLocationEta[31]= 0.008430536202;
    fWeight[31]= 0.00316596;
    fLocationKsi[32]= 0.632657968857;
    fLocationEta[32]= 0.358911494941;
    fWeight[32]= 0.00316596;
    fLocationKsi[33]= 0.008430536202;
    fLocationEta[33]= 0.632657968857;
    fWeight[33]= 0.00316596;
    fLocationKsi[34]= 0.294402476752;
    fLocationEta[34]= 0.574410971511;
    fWeight[34]= 0.0136288;
    fLocationKsi[35]= 0.574410971511;
    fLocationEta[35]= 0.131186551737;
    fWeight[35]= 0.0136288;
    fLocationKsi[36]= 0.131186551737;
    fLocationEta[36]= 0.294402476752;
    fWeight[36]= 0.0136288;
    fLocationKsi[37]= 0.294402476752;
    fLocationEta[37]= 0.131186551737;
    fWeight[37]= 0.0136288;
    fLocationKsi[38]= 0.574410971511;
    fLocationEta[38]= 0.294402476752;
    fWeight[38]= 0.0136288;
    fLocationKsi[39]= 0.131186551737;
    fLocationEta[39]= 0.574410971511;
    fWeight[39]= 0.0136288;
    fLocationKsi[40]= 0.325017801642;
    fLocationEta[40]= 0.624779046793;
    fWeight[40]= 0.00883839;
    fLocationKsi[41]= 0.624779046793;
    fLocationEta[41]= 0.050203151566;
    fWeight[41]= 0.00883839;
    fLocationKsi[42]= 0.050203151566;
    fLocationEta[42]= 0.325017801642;
    fWeight[42]= 0.00883839;
    fLocationKsi[43]= 0.325017801642;
    fLocationEta[43]= 0.050203151566;
    fWeight[43]= 0.00883839;
    fLocationKsi[44]= 0.624779046793;
    fLocationEta[44]= 0.325017801642;
    fWeight[44]= 0.00883839;
    fLocationKsi[45]= 0.050203151566;
    fLocationEta[45]= 0.624779046793;
    fWeight[45]= 0.00883839;
    fLocationKsi[46]= 0.184737559666;
    fLocationEta[46]= 0.748933176523;
    fWeight[46]= 0.00918974;
    fLocationKsi[47]= 0.748933176523;
    fLocationEta[47]= 0.066329263811;
    fWeight[47]= 0.00918974;
    fLocationKsi[48]= 0.066329263811;
    fLocationEta[48]= 0.184737559666;
    fWeight[48]= 0.00918974;
    fLocationKsi[49]= 0.184737559666;
    fLocationEta[49]= 0.066329263811;
    fWeight[49]= 0.00918974;
    fLocationKsi[50]= 0.748933176523;
    fLocationEta[50]= 0.184737559666;
    fWeight[50]= 0.00918974;
    fLocationKsi[51]= 0.066329263811;
    fLocationEta[51]= 0.748933176523;
    fWeight[51]= 0.00918974;
    fLocationKsi[52]= 0.218796800013;
    fLocationEta[52]= 0.769207005420;
    fWeight[52]= 0.00405237;
    fLocationKsi[53]= 0.769207005420;
    fLocationEta[53]= 0.011996194566;
    fWeight[53]= 0.00405237;
    fLocationKsi[54]= 0.011996194566;
    fLocationEta[54]= 0.218796800013;
    fWeight[54]= 0.00405237;
    fLocationKsi[55]= 0.218796800013;
    fLocationEta[55]= 0.011996194566;
    fWeight[55]= 0.00405237;
    fLocationKsi[56]= 0.769207005420;
    fLocationEta[56]= 0.218796800013;
    fWeight[56]= 0.00405237;
    fLocationKsi[57]= 0.011996194566;
    fLocationEta[57]= 0.769207005420;
    fWeight[57]= 0.00405237;
    fLocationKsi[58]= 0.101179597136;
    fLocationEta[58]= 0.883962302273;
    fWeight[58]= 0.00381706;
    fLocationKsi[59]= 0.883962302273;
    fLocationEta[59]= 0.014858100590;
    fWeight[59]= 0.00381706;
    fLocationKsi[60]= 0.014858100590;
    fLocationEta[60]= 0.101179597136;
    fWeight[60]= 0.00381706;
    fLocationKsi[61]= 0.101179597136;
    fLocationEta[61]= 0.014858100590;
    fWeight[61]= 0.00381706;
    fLocationKsi[62]= 0.883962302273;
    fLocationEta[62]= 0.101179597136;
    fWeight[62]= 0.00381706;
    fLocationKsi[63]= 0.014858100590;
    fLocationEta[63]= 0.883962302273;
    fWeight[63]= 0.00381706;
    fLocationKsi[64]= 0.020874755283;
    fLocationEta[64]= 1.014347260005;
    fWeight[64]= 2.30938e-05;
    fLocationKsi[65]= 1.014347260005;
    fLocationEta[65]= -0.035222015288;
    fWeight[65]= 2.30938e-05;
    fLocationKsi[66]= -0.035222015288;
    fLocationEta[66]= 0.020874755283;
    fWeight[66]= 2.30938e-05;
    fLocationKsi[67]= 0.020874755283;
    fLocationEta[67]= -0.035222015288;
    fWeight[67]= 2.30938e-05;
    fLocationKsi[68]= 1.014347260005;
    fLocationEta[68]= 0.020874755283;
    fWeight[68]= 2.30938e-05;
    fLocationKsi[69]= -0.035222015288;
    fLocationEta[69]= 1.014347260005;
    fWeight[69]= 2.30938e-05;
    break;

  case 19:
    fLocationKsi[0]= 0.333333333333;
    fLocationEta[0]= 0.333333333333;
    fWeight[0]= 0.0164532;
    fLocationKsi[1]= 0.489609987073;
    fLocationEta[1]= 0.489609987073;
    fWeight[1]= 0.00516537;
    fLocationKsi[2]= 0.489609987073;
    fLocationEta[2]= 0.020780025854;
    fWeight[2]= 0.00516537;
    fLocationKsi[3]= 0.020780025854;
    fLocationEta[3]= 0.489609987073;
    fWeight[3]= 0.00516537;
    fLocationKsi[4]= 0.454536892698;
    fLocationEta[4]= 0.454536892698;
    fWeight[4]= 0.0111936;
    fLocationKsi[5]= 0.454536892698;
    fLocationEta[5]= 0.090926214604;
    fWeight[5]= 0.0111936;
    fLocationKsi[6]= 0.090926214604;
    fLocationEta[6]= 0.454536892698;
    fWeight[6]= 0.0111936;
    fLocationKsi[7]= 0.401416680649;
    fLocationEta[7]= 0.401416680649;
    fWeight[7]= 0.0151331;
    fLocationKsi[8]= 0.401416680649;
    fLocationEta[8]= 0.197166638701;
    fWeight[8]= 0.0151331;
    fLocationKsi[9]= 0.197166638701;
    fLocationEta[9]= 0.401416680649;
    fWeight[9]= 0.0151331;
    fLocationKsi[10]= 0.255551654403;
    fLocationEta[10]= 0.255551654403;
    fWeight[10]= 0.0152455;
    fLocationKsi[11]= 0.255551654403;
    fLocationEta[11]= 0.488896691194;
    fWeight[11]= 0.0152455;
    fLocationKsi[12]= 0.488896691194;
    fLocationEta[12]= 0.255551654403;
    fWeight[12]= 0.0152455;
    fLocationKsi[13]= 0.177077942152;
    fLocationEta[13]= 0.177077942152;
    fWeight[13]= 0.0120796;
    fLocationKsi[14]= 0.177077942152;
    fLocationEta[14]= 0.645844115696;
    fWeight[14]= 0.0120796;
    fLocationKsi[15]= 0.645844115696;
    fLocationEta[15]= 0.177077942152;
    fWeight[15]= 0.0120796;
    fLocationKsi[16]= 0.110061053228;
    fLocationEta[16]= 0.110061053228;
    fWeight[16]= 0.0080254;
    fLocationKsi[17]= 0.110061053228;
    fLocationEta[17]= 0.779877893544;
    fWeight[17]= 0.0080254;
    fLocationKsi[18]= 0.779877893544;
    fLocationEta[18]= 0.110061053228;
    fWeight[18]= 0.0080254;
    fLocationKsi[19]= 0.055528624252;
    fLocationEta[19]= 0.055528624252;
    fWeight[19]= 0.00404229;
    fLocationKsi[20]= 0.055528624252;
    fLocationEta[20]= 0.888942751496;
    fWeight[20]= 0.00404229;
    fLocationKsi[21]= 0.888942751496;
    fLocationEta[21]= 0.055528624252;
    fWeight[21]= 0.00404229;
    fLocationKsi[22]= 0.012621863777;
    fLocationEta[22]= 0.012621863777;
    fWeight[22]= 0.00103968;
    fLocationKsi[23]= 0.012621863777;
    fLocationEta[23]= 0.974756272446;
    fWeight[23]= 0.00103968;
    fLocationKsi[24]= 0.974756272446;
    fLocationEta[24]= 0.012621863777;
    fWeight[24]= 0.00103968;
    fLocationKsi[25]= 0.395754787357;
    fLocationEta[25]= 0.600633794795;
    fWeight[25]= 0.00194244;
    fLocationKsi[26]= 0.600633794795;
    fLocationEta[26]= 0.003611417848;
    fWeight[26]= 0.00194244;
    fLocationKsi[27]= 0.003611417848;
    fLocationEta[27]= 0.395754787357;
    fWeight[27]= 0.00194244;
    fLocationKsi[28]= 0.395754787357;
    fLocationEta[28]= 0.003611417848;
    fWeight[28]= 0.00194244;
    fLocationKsi[29]= 0.600633794795;
    fLocationEta[29]= 0.395754787357;
    fWeight[29]= 0.00194244;
    fLocationKsi[30]= 0.003611417848;
    fLocationEta[30]= 0.600633794795;
    fWeight[30]= 0.00194244;
    fLocationKsi[31]= 0.307929983880;
    fLocationEta[31]= 0.557603261589;
    fWeight[31]= 0.0127871;
    fLocationKsi[32]= 0.557603261589;
    fLocationEta[32]= 0.134466754531;
    fWeight[32]= 0.0127871;
    fLocationKsi[33]= 0.134466754531;
    fLocationEta[33]= 0.307929983880;
    fWeight[33]= 0.0127871;
    fLocationKsi[34]= 0.307929983880;
    fLocationEta[34]= 0.134466754531;
    fWeight[34]= 0.0127871;
    fLocationKsi[35]= 0.557603261589;
    fLocationEta[35]= 0.307929983880;
    fWeight[35]= 0.0127871;
    fLocationKsi[36]= 0.134466754531;
    fLocationEta[36]= 0.557603261589;
    fWeight[36]= 0.0127871;
    fLocationKsi[37]= 0.264566948407;
    fLocationEta[37]= 0.720987025817;
    fWeight[37]= 0.00444045;
    fLocationKsi[38]= 0.720987025817;
    fLocationEta[38]= 0.014446025776;
    fWeight[38]= 0.00444045;
    fLocationKsi[39]= 0.014446025776;
    fLocationEta[39]= 0.264566948407;
    fWeight[39]= 0.00444045;
    fLocationKsi[40]= 0.264566948407;
    fLocationEta[40]= 0.014446025776;
    fWeight[40]= 0.00444045;
    fLocationKsi[41]= 0.720987025817;
    fLocationEta[41]= 0.264566948407;
    fWeight[41]= 0.00444045;
    fLocationKsi[42]= 0.014446025776;
    fLocationEta[42]= 0.720987025817;
    fWeight[42]= 0.00444045;
    fLocationKsi[43]= 0.358539352206;
    fLocationEta[43]= 0.594527068956;
    fWeight[43]= 0.00806227;
    fLocationKsi[44]= 0.594527068956;
    fLocationEta[44]= 0.046933578838;
    fWeight[44]= 0.00806227;
    fLocationKsi[45]= 0.046933578838;
    fLocationEta[45]= 0.358539352206;
    fWeight[45]= 0.00806227;
    fLocationKsi[46]= 0.358539352206;
    fLocationEta[46]= 0.046933578838;
    fWeight[46]= 0.00806227;
    fLocationKsi[47]= 0.594527068956;
    fLocationEta[47]= 0.358539352206;
    fWeight[47]= 0.00806227;
    fLocationKsi[48]= 0.046933578838;
    fLocationEta[48]= 0.594527068956;
    fWeight[48]= 0.00806227;
    fLocationKsi[49]= 0.157807405969;
    fLocationEta[49]= 0.839331473681;
    fWeight[49]= 0.00124597;
    fLocationKsi[50]= 0.839331473681;
    fLocationEta[50]= 0.002861120351;
    fWeight[50]= 0.00124597;
    fLocationKsi[51]= 0.002861120351;
    fLocationEta[51]= 0.157807405969;
    fWeight[51]= 0.00124597;
    fLocationKsi[52]= 0.157807405969;
    fLocationEta[52]= 0.002861120351;
    fWeight[52]= 0.00124597;
    fLocationKsi[53]= 0.839331473681;
    fLocationEta[53]= 0.157807405969;
    fWeight[53]= 0.00124597;
    fLocationKsi[54]= 0.002861120351;
    fLocationEta[54]= 0.839331473681;
    fWeight[54]= 0.00124597;
    fLocationKsi[55]= 0.075050596976;
    fLocationEta[55]= 0.701087978926;
    fWeight[55]= 0.00912142;
    fLocationKsi[56]= 0.701087978926;
    fLocationEta[56]= 0.223861424098;
    fWeight[56]= 0.00912142;
    fLocationKsi[57]= 0.223861424098;
    fLocationEta[57]= 0.075050596976;
    fWeight[57]= 0.00912142;
    fLocationKsi[58]= 0.075050596976;
    fLocationEta[58]= 0.223861424098;
    fWeight[58]= 0.00912142;
    fLocationKsi[59]= 0.701087978926;
    fLocationEta[59]= 0.075050596976;
    fWeight[59]= 0.00912142;
    fLocationKsi[60]= 0.223861424098;
    fLocationEta[60]= 0.701087978926;
    fWeight[60]= 0.00912142;
    fLocationKsi[61]= 0.142421601113;
    fLocationEta[61]= 0.822931324070;
    fWeight[61]= 0.00512928;
    fLocationKsi[62]= 0.822931324070;
    fLocationEta[62]= 0.034647074817;
    fWeight[62]= 0.00512928;
    fLocationKsi[63]= 0.034647074817;
    fLocationEta[63]= 0.142421601113;
    fWeight[63]= 0.00512928;
    fLocationKsi[64]= 0.142421601113;
    fLocationEta[64]= 0.034647074817;
    fWeight[64]= 0.00512928;
    fLocationKsi[65]= 0.822931324070;
    fLocationEta[65]= 0.142421601113;
    fWeight[65]= 0.00512928;
    fLocationKsi[66]= 0.034647074817;
    fLocationEta[66]= 0.822931324070;
    fWeight[66]= 0.00512928;
    fLocationKsi[67]= 0.065494628083;
    fLocationEta[67]= 0.924344252621;
    fWeight[67]= 0.00189996;
    fLocationKsi[68]= 0.924344252621;
    fLocationEta[68]= 0.010161119296;
    fWeight[68]= 0.00189996;
    fLocationKsi[69]= 0.010161119296;
    fLocationEta[69]= 0.065494628083;
    fWeight[69]= 0.00189996;
    fLocationKsi[70]= 0.065494628083;
    fLocationEta[70]= 0.010161119296;
    fWeight[70]= 0.00189996;
    fLocationKsi[71]= 0.924344252621;
    fLocationEta[71]= 0.065494628083;
    fWeight[71]= 0.00189996;
    fLocationKsi[72]= 0.010161119296;
    fLocationEta[72]= 0.924344252621;
    fWeight[72]= 0.00189996;
    break;

  default:
    PZError << "TPZIntRuleT creation : invalid number of integration points "
      " specified\n";
    //			PZError.show();

    fNumInt = 0;
  }
}

// ***************************************
// ***************************************

TPZIntRuleT::~TPZIntRuleT(){

  if (fLocationKsi) delete []fLocationKsi;
  if (fLocationEta) delete []fLocationEta;
  if (fWeight)      delete []fWeight;

}

//***************************************
//***************************************
void TPZIntRuleT::Loc(int i, TPZVec<REAL> &Points) {

  if (fLocationKsi && fLocationEta && i>=0 && i<fNumInt){
    Points[0] = fLocationKsi[i];
    Points[1] = fLocationEta[i];
    return;
  }
  else {
    PZError << "ERROR(TPZIntRuleT::loc) Out of bounds!!\n";
  }
}

//***************************************
//***************************************
REAL TPZIntRuleT::W(int i) {

  // Quando eu encontrar o IDIOTA que introduziu este bug!!!
  // Philippe
  if (fWeight && i>=0 && i<fNumInt)
    return fWeight[i];
    //    return fWeight[i]*0.5;
  else {
    PZError << "ERROR(TPZIntRule::w) Out of bounds!!\n";
    return 0.0;
  }
}

//***************************************
//***************************************
TPZIntRuleList::TPZIntRuleList(){

  static int	first = 1;

  if(first != 1) {
    PZError << "second initialization of the integration rule list\n"
      " something fishy is going on!\n";
    //		PZError.show();
  }

  first++;

  intavail    = NUMINT_RULES;//reta, quadrilatero, cubo
  intavailT   = NUMINT_RULEST;//triangulo
  intavailT3D = NUMINT_RULEST3D;//tetraedro
  intavailP3D = NUMINT_RULESP3D;//piramide

  intlist    = new TPZIntRule*[intavail];
  intlistT   = new TPZIntRuleT*[intavailT];
  intlistT3D = new TPZIntRuleT3D*[intavailT3D];
  intlistP3D = new TPZIntRuleP3D*[intavailP3D];


  if(intlist == NULL || intlistT == NULL || intlistT3D == NULL || intlistP3D == NULL){
    PZError << "TPZIntRuleList unable to initialize a list of integration rules\n";
    //		PZError.show();
    intavail = 0;

  } else {

    int i;
    for(i = 1; i<=intavail; ++i) {
      intlist[i-1] = new TPZIntRule(i-1);
      if(intlist[i-1] == NULL) {
	PZError << "TPZIntRuleList error: some integration rules"
	  "could not be initialized\n";
	//				PZError.show();
      }
    }

    for(i = 0; i<intavailT; ++i) {
      intlistT[i] = new TPZIntRuleT(i);
      if(intlistT[i] == NULL) {
	PZError << "TPZIntRuleList error: some integration rules"
	  "for triangles could not be initialized\n";
	//				PZError.show();
      }
    }
    for(i = 1; i<=intavailT3D; ++i) {
      intlistT3D[i-1] = new TPZIntRuleT3D(i);
      if(intlistT3D[i-1] == NULL) {
	PZError << "TPZIntRuleList error: some integration rules"
	  "for tetrahedros could not be initialized\n";
	//				PZError.show();
      }
    }
    for(i = 1; i<=intavailP3D; ++i) {
      intlistP3D[i-1] = new TPZIntRuleP3D(i);
      if(intlistP3D[i-1] == NULL) {
	PZError << "TPZIntRuleList error: some integration rules"
	  "for pyramides could not be initialized\n";
	//				PZError.show();
      }
    }


  }
}

//***************************************
//***************************************
TPZIntRuleList::~TPZIntRuleList(){

  if (!intlist) return;
  int i;
  for(i=0 ; i<intavail ; ++i)    if (intlist[i])    delete intlist[i];
  for(i=0 ; i<intavailT ; ++i)   if (intlistT[i])   delete intlistT[i];
  for(i=0 ; i<intavailT3D ; ++i) if (intlistT3D[i]) delete intlistT3D[i];
  for(i=0 ; i<intavailP3D ; ++i) if (intlistP3D[i]) delete intlistP3D[i];
  delete []intlistT3D;
  delete []intlistP3D;
  delete []intlist;
  delete []intlistT;
}

//***************************************
//***************************************
TPZIntRule* TPZIntRuleList::GetRule(int fNumInt) {

  if (fNumInt < 0 || fNumInt > intavail) {
    PZError << "\nERROR(TPZIntRuleList::getrule)-> Numint = " << fNumInt;
    //		PZError.show();
    fNumInt = intavailT-1;
    PZError << "\n                     precision obtained = " << fNumInt;
    //return NULL;
  }
  if(fNumInt == 0) fNumInt = 1;
  //return intlist[fNumInt-1];
  return intlist[fNumInt];
}

//**************************************
//**************************************
TPZIntRuleT* TPZIntRuleList::GetRuleT(int precision) {

  // <<<<<>>>>>
  if (precision < 0 || precision >= intavailT) {
    PZError << "\nERROR(TPZIntRuleList::getrule)-> precision required = " << precision;
    precision = intavailT-1;
    PZError << "\n                                 precision obtained = " << precision;
    //		PZError.show();
  }

  return intlistT[precision];
}

//**************************************
//**************************************
TPZInt1d::TPZInt1d(int OrdX){
  fOrdKsi 	= OrdX;
  fIntP 	= gIntRuleList.GetRule(OrdX);
}

void TPZInt1d::SetOrder(TPZVec<int> &ord){
  if(ord.NElements() < 1) {
    cout << "TPZINt1d::SetOrder: NULL number of integration points specified\n";
    return;
  }
  fOrdKsi = ord[0];
  fIntP   = gIntRuleList.GetRule(ord[0]);
}

void TPZInt1d::GetOrder(TPZVec<int> &ord) {
  ord[0] = (short) fOrdKsi;
}

int TPZInt1d::NPoints(){
  if(fIntP) return fIntP->NInt();
  PZError << "Null Pointer passed to method TPZInt1d::TPZInt1d(TPZIntRule *)\n";
  return 0;
}

void TPZInt1d::Point(int ip, TPZVec<REAL> &pos, REAL &w){
  if((fIntP) && ((ip >= 0) && (ip < NPoints()))){
    pos[0] 	= fIntP->Loc(ip);
    w        = fIntP->W(ip);
    return;
  }
  if(!fIntP)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZInt1d::TPZInt1d(TPZIntRule *)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << endl;
}

//**************************************
//**************************************
TPZIntQuad::TPZIntQuad(int OrdK, int OrdE){
  fOrdKsi = OrdK;
  fOrdEta = OrdE;
  fIntKsi = gIntRuleList.GetRule(OrdK);
  fIntEta = gIntRuleList.GetRule(OrdE);
}


int TPZIntQuad::NPoints(){
  if (!fIntKsi || !fIntEta){
    PZError << "Null Pointer passed to method TPZInt1d::TPZInt1d(TPZIntRule *)\n";
    return 0;
  }
  return fIntKsi->NInt() * fIntEta->NInt();
}

void TPZIntQuad::Point(int ip, TPZVec<REAL> &pos, REAL &w){
  if((fIntEta) && (fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    int ik, ie;
    ik = ip/fIntEta->NInt();
    ie = ip - (ip/fIntEta->NInt())*(fIntEta->NInt());
    pos[0] 	= fIntKsi->Loc(ik);
    pos[1]	= fIntEta->Loc(ie);
    w        = fIntKsi->W(ik)*fIntEta->W(ie);
    return;
  }
  if(!fIntKsi || !fIntEta)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZInt1d::TPZInt1d(TPZIntRule *)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << endl;
}

void TPZIntQuad::SetOrder(TPZVec<int> &ord){
  fOrdKsi = ord[0];
  fOrdEta = ord[1];
  fIntKsi = gIntRuleList.GetRule(ord[0]);
  fIntEta = gIntRuleList.GetRule(ord[1]);
}

void TPZIntQuad::GetOrder(TPZVec<int> &ord) {
  ord[0] = (short) fOrdKsi;
  ord[1] = (short) fOrdEta;
}

//**************************************
//**************************************
TPZIntTriang::TPZIntTriang(	int OrdK){
  fOrdKsi = OrdK;
  fIntKsi = gIntRuleList.GetRuleT(OrdK);
}

int  TPZIntTriang::NPoints(){
  if (!fIntKsi){
    PZError << "Null Pointer passed to method TPZIntTriang::NPoints()\n";
    return 0;
  }
  return fIntKsi->NInt();
}

void TPZIntTriang::Point(int ip, TPZVec<REAL> &pos, REAL &w){
  if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    fIntKsi->Loc(ip, pos);
    w = fIntKsi->W(ip);
    return;
  }
  if(!fIntKsi)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntTriang::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << endl;
}

void TPZIntTriang::SetOrder(TPZVec<int> &ord){
  fOrdKsi = ord[0];
  if(ord[1]>ord[0]) fOrdKsi = ord[1];
  fIntKsi = gIntRuleList.GetRuleT(fOrdKsi);
  if(ord[0] <0 || ord[0]>11 || ord[1] <0 || ord[1]>11) fOrdKsi = 10;//havendo erro assume a maxima ordem
}

void TPZIntTriang::GetOrder(TPZVec<int> &ord) {
  ord[0] = (short) fOrdKsi;
  ord[1] = ord[0];
}
//##############################################################################
//##############################################################################
TPZIntCube3D::TPZIntCube3D(int OrdK, int OrdE, int OrdZ){
  fOrdKsi = OrdK;
  fOrdEta = OrdE;
  fOrdZeta = OrdZ;
  fIntKsi = gIntRuleList.GetRule(OrdK);
  fIntEta = gIntRuleList.GetRule(OrdE);
  fIntZeta = gIntRuleList.GetRule(OrdZ);
}
//------------------------------------------------------------------------------
void TPZIntCube3D::SetOrder(TPZVec<int> &ord){
  fOrdKsi = ord[0];
  fOrdEta = ord[1];
  fOrdZeta = ord[2];
  fIntKsi = gIntRuleList.GetRule(ord[0]);
  fIntEta = gIntRuleList.GetRule(ord[1]);
  fIntZeta = gIntRuleList.GetRule(ord[2]);
}
//------------------------------------------------------------------------------
int TPZIntCube3D::NPoints(){
  if (!fIntKsi || !fIntEta|| !fIntZeta){
    PZError << "Null Pointer passed to method TPZIntCube3D::NPoints()\n";
    return 0;
  }
  return fIntKsi->NInt() * fIntEta->NInt() * fIntZeta->NInt();
}
//------------------------------------------------------------------------------
void TPZIntCube3D::Point(int ip, TPZVec<REAL> &pos, REAL &w){

  if((fIntZeta) && (fIntEta) && (fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    int ik, ie , iz;
    int IK=0;

    if (ip > fIntKsi->NInt() * fIntEta->NInt() * fIntZeta->NInt()) {
      cout<<"\nPonto de Integracao numero "<<ip<<" fora do intervalo"<<endl;
      return;
    }

    while ( ip >= IK*fIntZeta->NInt()*fIntEta->NInt() ) IK++;
    IK--;
    int IP = ip - IK*fIntZeta->NInt()*fIntEta->NInt();


    if (IP <  fIntZeta->NInt()) {ik = IK; ie = 0; iz = IP;} else
      if (IP >=  fIntZeta->NInt() && IP < fIntZeta->NInt()*fIntEta->NInt()) {
	for (int i=0;i<fIntEta->NInt();i++) {
	  if (IP >= fIntZeta->NInt()+i*fIntEta->NInt() &&
	      IP < fIntZeta->NInt()+(i+1)*fIntEta->NInt())
	    {ik = IK; ie = (i+1); iz = IP - (i+1)*fIntZeta->NInt();}
	}
      }
    pos[0] 	= fIntKsi->Loc(ik);
    pos[1]	= fIntEta->Loc(ie);
    pos[2]	= fIntZeta->Loc(iz);
    w        = fIntKsi->W(ik)*fIntEta->W(ie)*fIntZeta->W(iz);
    return;
  }
  if(!fIntKsi || !fIntEta || !fIntZeta)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntCube3D::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << endl;
}
//------------------------------------------------------------------------------
void TPZIntCube3D::GetOrder(TPZVec<int> &ord) {
  ord[0] = (short) fOrdKsi;
  ord[1] = (short) fOrdEta;
  ord[2] = (short) fOrdZeta;
}
//##############################################################################
//##############################################################################
TPZIntTetra3D::TPZIntTetra3D(int OrdK){
  fOrdKsi = OrdK;
  fIntKsi = gIntRuleList.GetRuleT3D(OrdK);
}
//**************************************
TPZIntRuleT3D* TPZIntRuleList::GetRuleT3D(int precision) {

  if (precision >= intavailT3D) {
    PZError << "\nERROR(TPZIntRuleList::getrule)-> precision required = " << precision << endl;
    precision = intavailT3D-1;
    PZError << "\nERROR(TPZIntRuleList::getrule)-> precision gotten = " << precision << endl;
  }
  if (precision < 0) {
    PZError << "\nERROR(TPZIntRuleList::getrule)-> precision required = " << precision << endl;
    precision = 0;
    PZError << "\nERROR(TPZIntRuleList::getrule)-> precision gotten = " << precision << endl;
  }

  return intlistT3D[precision];
}
//------------------------------------------------------------------------------
void TPZIntTetra3D::SetOrder(TPZVec<int> &ord){

  fOrdKsi = (ord[1] > ord[0]) ? ord[1] : ord[0];
  fOrdKsi = (fOrdKsi > ord[2]) ? fOrdKsi : ord[2];
  fIntKsi = gIntRuleList.GetRuleT3D(fOrdKsi);
}
//------------------------------------------------------------------------------
int TPZIntTetra3D::NPoints(){
  if (!fIntKsi){
    PZError << "Null Pointer passed to method TPZIntTetra3D::NPoints()\n";
    return 0;
  }
  return fIntKsi->NInt();
}
//------------------------------------------------------------------------------
void TPZIntTetra3D::Point(int ip, TPZVec<REAL> &pos, REAL &w){

  if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    fIntKsi->Loc(ip, pos);
    w = fIntKsi->W(ip);
    return;
  }
  if(!fIntKsi)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntTetra3D::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << endl;
}
//------------------------------------------------------------------------------
REAL TPZIntRuleT3D::W(int i) {

  if (fWeight && i>=0 && i<fNumInt)
    return fWeight[i];
  else {
    PZError << "ERROR(TPZIntRuleT3D::w) Out of bounds!!\n";
    return 0.0;
  }
}
//------------------------------------------------------------------------------
void TPZIntTetra3D::GetOrder(TPZVec<int> &ord) {
  ord[0] = (short) fOrdKsi;
  ord[1] = ord[0];
  ord[2] = ord[0];
}
//------------------------------------------------------------------------------
TPZIntRuleT3D::~TPZIntRuleT3D(){

  if (fLocationKsi) delete []fLocationKsi;
  if (fLocationEta) delete []fLocationEta;
  if (fLocationZeta) delete []fLocationZeta;
  if (fWeight)   delete []fWeight;

}
//------------------------------------------------------------------------------
void TPZIntRuleT3D::Loc(int i, TPZVec<REAL> &Points) {

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

  fNumInt = (short) precision;

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

  case 24://integra até grau 4

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
//##############################################################################
//##############################################################################
TPZIntPrism3D::TPZIntPrism3D(int OrdK,int OrdL) : fIntRule1D(OrdK), fIntTriang(OrdL) {
  fOrdKsi = OrdK;
  fOrdKti = OrdL;
}
TPZIntPrism3D::~TPZIntPrism3D() {}
//------------------------------------------------------------------------------
void TPZIntPrism3D::SetOrder(TPZVec<int> &ord){

  fOrdKsi = ord[0];//ordem na reta : zeta
  fOrdKti = (ord[1] > ord[2]) ? ord[1] : ord[2];//ordem no plano XY
  TPZVec<int> prc1(1),prc2(2);
  prc1[0] = ord[0];
  prc2[0] = ord[1];
  prc2[1] = ord[2];
  fIntRule1D.SetOrder(prc1);
  fIntTriang.SetOrder(prc2);
}
//------------------------------------------------------------------------------
int TPZIntPrism3D::NPoints(){
  if (!&fIntRule1D || !&fIntTriang){
    PZError << "Null Pointer passed to method TPZIntPrism3D::NPoints()\n";
    return 0;
  }
  return fIntRule1D.NPoints()*fIntTriang.NPoints();
}
//------------------------------------------------------------------------------
void TPZIntPrism3D::Point(int ip, TPZVec<REAL> &pos, REAL &w){

  if((&fIntRule1D && &fIntTriang) && ((ip >= 0) && (ip < NPoints()))){

    REAL v;
    TPZVec<REAL> ps(2);
    fIntTriang.Point(ip % fIntTriang.NPoints(), ps, v);
    pos[0] = ps[0]; pos[1] = ps[1];
    fIntRule1D.Point(ip / fIntTriang.NPoints(), ps, w);
    pos[2] = ps[0];
    w *= v;
    return;
  }
  if(!&fIntRule1D || !&fIntTriang)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntPrism3D::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << endl;
}
//------------------------------------------------------------------------------
void TPZIntPrism3D::GetOrder(TPZVec<int> &ord) {
  ord[0] = fOrdKsi;
  ord[1] = fOrdKti;
  ord[2] = ord[1];
}
//##############################################################################
//##############################################################################
TPZIntPyram3D::TPZIntPyram3D(int OrdK){
  fOrdKsi = OrdK;
  fIntKsi = gIntRuleList.GetRuleP3D(OrdK);
}
//**************************************
TPZIntRuleP3D* TPZIntRuleList::GetRuleP3D(int precision) {

  // <<<<<>>>>>
  if (precision < 1 || precision > intavailP3D) {
    PZError << "\nERROR(TPZIntRuleList::getrule)-> precision = " << precision;
    //		PZError.show();
    precision = intavailP3D;
  }

  return intlistP3D[precision-1];
}
//------------------------------------------------------------------------------
void TPZIntPyram3D::SetOrder(TPZVec<int> &ord){

  fOrdKsi = (ord[1] > ord[0]) ? ord[1] : ord[0];
  fOrdKsi = (fOrdKsi > ord[2]) ? fOrdKsi : ord[2];
  fIntKsi = gIntRuleList.GetRuleP3D(fOrdKsi);
}
//------------------------------------------------------------------------------
int TPZIntPyram3D::NPoints(){
  if (!fIntKsi){
    PZError << "Null Pointer passed to method TPZIntPyram3D::NPoints()\n";
    return 0;
  }
  return fIntKsi->NInt();
}
//------------------------------------------------------------------------------
void TPZIntPyram3D::Point(int ip, TPZVec<REAL> &pos, REAL &w){

  if((fIntKsi) && ((ip >= 0) && (ip < NPoints()))){
    fIntKsi->Loc(ip, pos);
    w = fIntKsi->W(ip);
    return;
  }
  if(!fIntKsi)
    PZError 	<< "Null Pointer passed to method "
		<< "TPZIntPyram3D::Point(..)\n";
  if((ip < 0) || (ip >= NPoints()))
    PZError 	<< "ip = " << ip << ", Out of Range: 0 -> "
		<< NPoints() << endl;
}
//------------------------------------------------------------------------------
REAL TPZIntRuleP3D::W(int i) {

  if (fWeight && i>=0 && i<fNumInt)
    return fWeight[i];
  else {
    PZError << "ERROR(TPZIntRuleP3D::w) Out of bounds!!\n";
    return 0.0;
  }
}
//------------------------------------------------------------------------------
void TPZIntPyram3D::GetOrder(TPZVec<int> &ord) {
  ord[0] = (short) fOrdKsi;
  ord[1] = ord[0];
  ord[2] = ord[0];
}
//------------------------------------------------------------------------------
TPZIntRuleP3D::~TPZIntRuleP3D(){

  if (fLocationKsi) delete []fLocationKsi;
  if (fLocationEta) delete []fLocationEta;
  if (fLocationZeta) delete []fLocationZeta;
  if (fWeight)   delete []fWeight;

}
//------------------------------------------------------------------------------
void TPZIntRuleP3D::Loc(int i, TPZVec<REAL> &Points) {

  if (fLocationKsi && fLocationEta && fLocationZeta && i>=0 && i<fNumInt){
    Points[0] = fLocationKsi[i];
    Points[1] = fLocationEta[i];
    Points[2] = fLocationZeta[i];
    return;
  }
  else {
    PZError << "ERROR(TPZIntRuleP3D::loc) Out of bounds!!\n";
  }
}
//------------------------------------------------------------------------------
TPZIntRuleP3D::TPZIntRuleP3D(int precision){

  if(precision < 1 && precision > NUMINT_RULEST3D){
    PZError << "TPZIntRule creation precision = " << precision << " not available\n";
    //		PZError.show();
    precision = NUMINT_RULEST3D;
  }

  fNumInt = (short) precision;
  //2(ord-1)
  if (precision >= 2)  fNumInt =  12;//2 integra lineares e bilineares (xy) e quadraticas
  if (precision >  2)  fNumInt =  36;//3 e 4 integra de grau 4
  if (precision >  4)  fNumInt =  96;//5 e 6 integra de grau 6
  if (precision >  6)  fNumInt = 200;//7 e 8 integra de grau 8

  fLocationKsi  = new REAL[fNumInt];
  fLocationEta  = new REAL[fNumInt];
  fLocationZeta = new REAL[fNumInt];
  fWeight = new REAL[fNumInt];

  if(fLocationKsi == NULL || fLocationEta == NULL || fLocationZeta == NULL || fWeight == NULL){
    fNumInt = 0;
    return;
  }

  switch(fNumInt){

  case 1:
    break;

  case 12://integra quadraticas  xy , xz , x2

    fLocationKsi[0]  = -0.065068336845;
    fLocationEta[0]  = -0.065068336845;
    fLocationZeta[0] = 0.887298334621;
    fWeight[0]       = 0.003528240383;

    fLocationKsi[1]  = -0.288675134595;
    fLocationEta[1]  = -0.288675134595;
    fLocationZeta[1] = 0.500000000000;
    fWeight[1]       = 0.111111111111;

    fLocationKsi[2]  = -0.512281932345;
    fLocationEta[2]  = -0.512281932345;
    fLocationZeta[2] = 0.112701665379;
    fWeight[2]       = 0.218693981839;

    fLocationKsi[3]  = -0.065068336845;
    fLocationEta[3]  = 0.065068336845;
    fLocationZeta[3] = 0.887298334621;
    fWeight[3]       = 0.003528240383;

    fLocationKsi[4]  = -0.288675134595;
    fLocationEta[4]  = 0.288675134595;
    fLocationZeta[4] = 0.500000000000;
    fWeight[4]       = 0.111111111111;

    fLocationKsi[5]  = -0.512281932345;
    fLocationEta[5]  = 0.512281932345;
    fLocationZeta[5] = 0.112701665379;
    fWeight[5]       = 0.218693981839;

    fLocationKsi[6]  = 0.065068336845;
    fLocationEta[6]  = -0.065068336845;
    fLocationZeta[6] = 0.887298334621;
    fWeight[6]       = 0.003528240383;

    fLocationKsi[7]  = 0.288675134595;
    fLocationEta[7]  = -0.288675134595;
    fLocationZeta[7] = 0.500000000000;
    fWeight[7]       = 0.111111111111;

    fLocationKsi[8]  = 0.512281932345;
    fLocationEta[8]  = -0.512281932345;
    fLocationZeta[8] = 0.112701665379;
    fWeight[8]       = 0.218693981839;

    fLocationKsi[9]  = 0.065068336845;
    fLocationEta[9]  = 0.065068336845;
    fLocationZeta[9] = 0.887298334621;
    fWeight[9]       = 0.003528240383;

    fLocationKsi[10]  = 0.288675134595;
    fLocationEta[10]  = 0.288675134595;
    fLocationZeta[10] = 0.500000000000;
    fWeight[10]       = 0.111111111111;

    fLocationKsi[11]  = 0.512281932345;
    fLocationEta[11]  = 0.512281932345;
    fLocationZeta[11] = 0.112701665379;
    fWeight[11]       = 0.218693981839;
    break;


  case 36://integra x4 , y4 , x3y , z4 , y2z2
    fLocationKsi[0]  = -0.053781675257;
    fLocationEta[0]  = -0.053781675257;
    fLocationZeta[0] = 0.930568155800;
    fWeight[0]       = 0.000258785806;

    fLocationKsi[1]  = -0.255624242632;
    fLocationEta[1]  = -0.255624242632;
    fLocationZeta[1] = 0.669990521800;
    fWeight[1]       = 0.010960291203;

    fLocationKsi[2]  = -0.518972426610;
    fLocationEta[2]  = -0.518972426610;
    fLocationZeta[2] = 0.330009478200;
    fWeight[2]       = 0.045175876132;

    fLocationKsi[3]  = -0.720814993985;
    fLocationEta[3]  = -0.720814993985;
    fLocationZeta[3] = 0.069431844200;
    fWeight[3]       = 0.046485705294;

    fLocationKsi[4]  = -0.053781675257;
    fLocationEta[4]  = 0.000000000000;
    fLocationZeta[4] = 0.930568155800;
    fWeight[4]       = 0.000414057290;

    fLocationKsi[5]  = -0.255624242632;
    fLocationEta[5]  = 0.000000000000;
    fLocationZeta[5] = 0.669990521800;
    fWeight[5]       = 0.017536465925;

    fLocationKsi[6]  = -0.518972426610;
    fLocationEta[6]  = 0.000000000000;
    fLocationZeta[6] = 0.330009478200;
    fWeight[6]       = 0.072281401811;

    fLocationKsi[7]  = -0.720814993985;
    fLocationEta[7]  = 0.000000000000;
    fLocationZeta[7] = 0.069431844200;
    fWeight[7]       = 0.074377128471;

    fLocationKsi[8]  = -0.053781675257;
    fLocationEta[8]  = 0.053781675257;
    fLocationZeta[8] = 0.930568155800;
    fWeight[8]       = 0.000258785806;

    fLocationKsi[9]  = -0.255624242632;
    fLocationEta[9]  = 0.255624242632;
    fLocationZeta[9] = 0.669990521800;
    fWeight[9]       = 0.010960291203;

    fLocationKsi[10]  = -0.518972426610;
    fLocationEta[10]  = 0.518972426610;
    fLocationZeta[10] = 0.330009478200;
    fWeight[10]       = 0.045175876132;

    fLocationKsi[11]  = -0.720814993985;
    fLocationEta[11]  = 0.720814993985;
    fLocationZeta[11] = 0.069431844200;
    fWeight[11]       = 0.046485705294;

    fLocationKsi[12]  = 0.000000000000;
    fLocationEta[12]  = -0.053781675257;
    fLocationZeta[12] = 0.930568155800;
    fWeight[12]       = 0.000414057290;

    fLocationKsi[13]  = 0.000000000000;
    fLocationEta[13]  = -0.255624242632;
    fLocationZeta[13] = 0.669990521800;
    fWeight[13]       = 0.017536465925;

    fLocationKsi[14]  = 0.000000000000;
    fLocationEta[14]  = -0.518972426610;
    fLocationZeta[14] = 0.330009478200;
    fWeight[14]       = 0.072281401811;

    fLocationKsi[15]  = 0.000000000000;
    fLocationEta[15]  = -0.720814993985;
    fLocationZeta[15] = 0.069431844200;
    fWeight[15]       = 0.074377128471;

    fLocationKsi[16]  = 0.000000000000;
    fLocationEta[16]  = 0.000000000000;
    fLocationZeta[16] = 0.930568155800;
    fWeight[16]       = 0.000662491664;

    fLocationKsi[17]  = 0.000000000000;
    fLocationEta[17]  = 0.000000000000;
    fLocationZeta[17] = 0.669990521800;
    fWeight[17]       = 0.028058345479;

    fLocationKsi[18]  = 0.000000000000;
    fLocationEta[18]  = 0.000000000000;
    fLocationZeta[18] = 0.330009478200;
    fWeight[18]       = 0.115650242898;

    fLocationKsi[19]  = 0.000000000000;
    fLocationEta[19]  = 0.000000000000;
    fLocationZeta[19] = 0.069431844200;
    fWeight[19]       = 0.119003405553;

    fLocationKsi[20]  = 0.000000000000;
    fLocationEta[20]  = 0.053781675257;
    fLocationZeta[20] = 0.930568155800;
    fWeight[20]       = 0.000414057290;

    fLocationKsi[21]  = 0.000000000000;
    fLocationEta[21]  = 0.255624242632;
    fLocationZeta[21] = 0.669990521800;
    fWeight[21]       = 0.017536465925;

    fLocationKsi[22]  = 0.000000000000;
    fLocationEta[22]  = 0.518972426610;
    fLocationZeta[22] = 0.330009478200;
    fWeight[22]       = 0.072281401811;

    fLocationKsi[23]  = 0.000000000000;
    fLocationEta[23]  = 0.720814993985;
    fLocationZeta[23] = 0.069431844200;
    fWeight[23]       = 0.074377128471;

    fLocationKsi[24]  = 0.053781675257;
    fLocationEta[24]  = -0.053781675257;
    fLocationZeta[24] = 0.930568155800;
    fWeight[24]       = 0.000258785806;

    fLocationKsi[25]  = 0.255624242632;
    fLocationEta[25]  = -0.255624242632;
    fLocationZeta[25] = 0.669990521800;
    fWeight[25]       = 0.010960291203;

    fLocationKsi[26]  = 0.518972426610;
    fLocationEta[26]  = -0.518972426610;
    fLocationZeta[26] = 0.330009478200;
    fWeight[26]       = 0.045175876132;

    fLocationKsi[27]  = 0.720814993985;
    fLocationEta[27]  = -0.720814993985;
    fLocationZeta[27] = 0.069431844200;
    fWeight[27]       = 0.046485705294;

    fLocationKsi[28]  = 0.053781675257;
    fLocationEta[28]  = 0.000000000000;
    fLocationZeta[28] = 0.930568155800;
    fWeight[28]       = 0.000414057290;

    fLocationKsi[29]  = 0.255624242632;
    fLocationEta[29]  = 0.000000000000;
    fLocationZeta[29] = 0.669990521800;
    fWeight[29]       = 0.017536465925;

    fLocationKsi[30]  = 0.518972426610;
    fLocationEta[30]  = 0.000000000000;
    fLocationZeta[30] = 0.330009478200;
    fWeight[30]       = 0.072281401811;

    fLocationKsi[31]  = 0.720814993985;
    fLocationEta[31]  = 0.000000000000;
    fLocationZeta[31] = 0.069431844200;
    fWeight[31]       = 0.074377128471;

    fLocationKsi[32]  = 0.053781675257;
    fLocationEta[32]  = 0.053781675257;
    fLocationZeta[32] = 0.930568155800;
    fWeight[32]       = 0.000258785806;

    fLocationKsi[33]  = 0.255624242632;
    fLocationEta[33]  = 0.255624242632;
    fLocationZeta[33] = 0.669990521800;
    fWeight[33]       = 0.010960291203;

    fLocationKsi[34]  = 0.518972426610;
    fLocationEta[34]  = 0.518972426610;
    fLocationZeta[34] = 0.330009478200;
    fWeight[34]       = 0.045175876132;

    fLocationKsi[35]  = 0.720814993985;
    fLocationEta[35]  = 0.720814993985;
    fLocationZeta[35] = 0.069431844200;
    fWeight[35]       = 0.046485705294;

    break;

/*  case 80://x6y2 , x4y4 , x2y2z4
    
    fLocationKsi[0]  = -0.040395970728;
    fLocationEta[0]  = -0.040395970728;
    fLocationZeta[0] = 0.953089922950;
    fWeight[0]       = 0.000031543709;
    
    fLocationKsi[1]  = -0.198720417995;
    fLocationEta[1]  = -0.198720417995;
    fLocationZeta[1] = 0.769234655050;
    fWeight[1]       = 0.001542076812;
    
    fLocationKsi[2]  = -0.430568155800;
    fLocationEta[2]  = -0.430568155800;
    fLocationZeta[2] = 0.500000000000;
    fWeight[2]       = 0.008604657297;

    fLocationKsi[3]  = -0.662415893605;
    fLocationEta[3]  = -0.662415893605;
    fLocationZeta[3] = 0.230765344950;
    fWeight[3]       = 0.017134936955;
    
    fLocationKsi[4]  = -0.820740340872;
    fLocationEta[4]  = -0.820740340872;
    fLocationZeta[4] = 0.046910077050;
    fWeight[4]       = 0.013021116305;
    
    fLocationKsi[5]  = -0.040395970728;
    fLocationEta[5]  = -0.015948536951;
    fLocationZeta[5] = 0.953089922950;
    fWeight[5]       = 0.000059136957;
    
    fLocationKsi[6]  = -0.198720417995;
    fLocationEta[6]  = -0.078455842803;
    fLocationZeta[6] = 0.769234655050;
    fWeight[6]       = 0.002891027494;
    
    fLocationKsi[7]  = -0.430568155800;
    fLocationEta[7]  = -0.169990521800;
    fLocationZeta[7] = 0.500000000000;
    fWeight[7]       = 0.016131687239;
    
    fLocationKsi[8]  = -0.662415893605;
    fLocationEta[8]  = -0.261525200797;
    fLocationZeta[8] = 0.230765344950;
    fWeight[8]       = 0.032123934084;
    
    fLocationKsi[9]  = -0.820740340872;
    fLocationEta[9]  = -0.324032506649;
    fLocationZeta[9] = 0.046910077050;
    fWeight[9]       = 0.024411498156;
    
    fLocationKsi[10]  = -0.040395970728;
    fLocationEta[10]  = 0.015948536951;
    fLocationZeta[10] = 0.953089922950;
    fWeight[10]       = 0.000059136957;
    
    fLocationKsi[11]  = -0.198720417995;
    fLocationEta[11]  = 0.078455842803;
    fLocationZeta[11] = 0.769234655050;
    fWeight[11]       = 0.002891027494;
    
    fLocationKsi[12]  = -0.430568155800;
    fLocationEta[12]  = 0.169990521800;
    fLocationZeta[12] = 0.500000000000;
    fWeight[12]       = 0.016131687239;
    
    fLocationKsi[13]  = -0.662415893605;
    fLocationEta[13]  = 0.261525200797;
    fLocationZeta[13] = 0.230765344950;
    fWeight[13]       = 0.032123934084;
    
    fLocationKsi[14]  = -0.820740340872;
    fLocationEta[14]  = 0.324032506649;
    fLocationZeta[14] = 0.046910077050;
    fWeight[14]       = 0.024411498156;
    
    fLocationKsi[15]  = -0.040395970728;
    fLocationEta[15]  = 0.040395970728;
    fLocationZeta[15] = 0.953089922950;
    fWeight[15]       = 0.000031543709;
    
    fLocationKsi[16]  = -0.198720417995;
    fLocationEta[16]  = 0.198720417995;
    fLocationZeta[16] = 0.769234655050;
    fWeight[16]       = 0.001542076812;
    
    fLocationKsi[17]  = -0.430568155800;
    fLocationEta[17]  = 0.430568155800;
    fLocationZeta[17] = 0.500000000000;
    fWeight[17]       = 0.008604657297;
    
    fLocationKsi[18]  = -0.662415893605;
    fLocationEta[18]  = 0.662415893605;
    fLocationZeta[18] = 0.230765344950;
    fWeight[18]       = 0.017134936955;
    
    fLocationKsi[19]  = -0.820740340872;
    fLocationEta[19]  = 0.820740340872;
    fLocationZeta[19] = 0.046910077050;
    fWeight[19]       = 0.013021116305;
    
    fLocationKsi[20]  = -0.015948536951;
    fLocationEta[20]  = -0.040395970728;
    fLocationZeta[20] = 0.953089922950;
    fWeight[20]       = 0.000059136957;
    
    fLocationKsi[21]  = -0.078455842803;
    fLocationEta[21]  = -0.198720417995;
    fLocationZeta[21] = 0.769234655050;
    fWeight[21]       = 0.002891027494;
    
    fLocationKsi[22]  = -0.169990521800;
    fLocationEta[22]  = -0.430568155800;
    fLocationZeta[22] = 0.500000000000;
    fWeight[22]       = 0.016131687239;
    
    fLocationKsi[23]  = -0.261525200797;
    fLocationEta[23]  = -0.662415893605;
    fLocationZeta[23] = 0.230765344950;
    fWeight[23]       = 0.032123934084;
    
    fLocationKsi[24]  = -0.324032506649;
    fLocationEta[24]  = -0.820740340872;
    fLocationZeta[24] = 0.046910077050;
    fWeight[24]       = 0.024411498156;
    
    fLocationKsi[25]  = -0.015948536951;
    fLocationEta[25]  = -0.015948536951;
    fLocationZeta[25] = 0.953089922950;
    fWeight[25]       = 0.000110867738;
    
    fLocationKsi[26]  = -0.078455842803;
    fLocationEta[26]  = -0.078455842803;
    fLocationZeta[26] = 0.769234655050;
    fWeight[26]       = 0.005419989399;
    
    fLocationKsi[27]  = -0.169990521800;
    fLocationEta[27]  = -0.169990521800;
    fLocationZeta[27] = 0.500000000000;
    fWeight[27]       = 0.030243079324;
    
    fLocationKsi[28]  = -0.261525200797;
    fLocationEta[28]  = -0.261525200797;
    fLocationZeta[28] = 0.230765344950;
    fWeight[28]       = 0.060224741049;
    
    fLocationKsi[29]  = -0.324032506649;
    fLocationEta[29]  = -0.324032506649;
    fLocationZeta[29] = 0.046910077050;
    fWeight[29]       = 0.045765756811;
    
    fLocationKsi[30]  = -0.015948536951;
    fLocationEta[30]  = 0.015948536951;
    fLocationZeta[30] = 0.953089922950;
    fWeight[30]       = 0.000110867738;
    
    fLocationKsi[31]  = -0.078455842803;
    fLocationEta[31]  = 0.078455842803;
    fLocationZeta[31] = 0.769234655050;
    fWeight[31]       = 0.005419989399;
    
    fLocationKsi[32]  = -0.169990521800;
    fLocationEta[32]  = 0.169990521800;
    fLocationZeta[32] = 0.500000000000;
    fWeight[32]       = 0.030243079324;
    
    fLocationKsi[33]  = -0.261525200797;
    fLocationEta[33]  = 0.261525200797;
    fLocationZeta[33] = 0.230765344950;
    fWeight[33]       = 0.060224741049;
    
    fLocationKsi[34]  = -0.324032506649;
    fLocationEta[34]  = 0.324032506649;
    fLocationZeta[34] = 0.046910077050;
    fWeight[34]       = 0.045765756811;
    
    fLocationKsi[35]  = -0.015948536951;
    fLocationEta[35]  = 0.040395970728;
    fLocationZeta[35] = 0.953089922950;
    fWeight[35]       = 0.000059136957;

    fLocationKsi[36]  = -0.078455842803;
    fLocationEta[36]  = 0.198720417995;
    fLocationZeta[36] = 0.769234655050;
    fWeight[36]       = 0.002891027494;
    
    fLocationKsi[37]  = -0.169990521800;
    fLocationEta[37]  = 0.430568155800;
    fLocationZeta[37] = 0.500000000000;
    fWeight[37]       = 0.016131687239;
    
    fLocationKsi[38]  = -0.261525200797;
    fLocationEta[38]  = 0.662415893605;
    fLocationZeta[38] = 0.230765344950;
    fWeight[38]       = 0.032123934084;
    
    fLocationKsi[39]  = -0.324032506649;
    fLocationEta[39]  = 0.820740340872;
    fLocationZeta[39] = 0.046910077050;
    fWeight[39]       = 0.024411498156;
    
    fLocationKsi[40]  = 0.015948536951;
    fLocationEta[40]  = -0.040395970728;
    fLocationZeta[40] = 0.953089922950;
    fWeight[40]       = 0.000059136957;
    
    fLocationKsi[41]  = 0.078455842803;
    fLocationEta[41]  = -0.198720417995;
    fLocationZeta[41] = 0.769234655050;
    fWeight[41]       = 0.002891027494;
    
    fLocationKsi[42]  = 0.169990521800;
    fLocationEta[42]  = -0.430568155800;
    fLocationZeta[42] = 0.500000000000;
    fWeight[42]       = 0.016131687239;
    
    fLocationKsi[43]  = 0.261525200797;
    fLocationEta[43]  = -0.662415893605;
    fLocationZeta[43] = 0.230765344950;
    fWeight[43]       = 0.032123934084;
    
    fLocationKsi[44]  = 0.324032506649;
    fLocationEta[44]  = -0.820740340872;
    fLocationZeta[44] = 0.046910077050;
    fWeight[44]       = 0.024411498156;
    
    fLocationKsi[45]  = 0.015948536951;
    fLocationEta[45]  = -0.015948536951;
    fLocationZeta[45] = 0.953089922950;
    fWeight[45]       = 0.000110867738;
    
    fLocationKsi[46]  = 0.078455842803;
    fLocationEta[46]  = -0.078455842803;
    fLocationZeta[46] = 0.769234655050;
    fWeight[46]       = 0.005419989399;
    
    fLocationKsi[47]  = 0.169990521800;
    fLocationEta[47]  = -0.169990521800;
    fLocationZeta[47] = 0.500000000000;
    fWeight[47]       = 0.030243079324;
    
    fLocationKsi[48]  = 0.261525200797;
    fLocationEta[48]  = -0.261525200797;
    fLocationZeta[48] = 0.230765344950;
    fWeight[48]       = 0.060224741049;
    
    fLocationKsi[49]  = 0.324032506649;
    fLocationEta[49]  = -0.324032506649;
    fLocationZeta[49] = 0.046910077050;
    fWeight[49]       = 0.045765756811;
    
    fLocationKsi[50]  = 0.015948536951;
    fLocationEta[50]  = 0.015948536951;
    fLocationZeta[50] = 0.953089922950;
    fWeight[50]       = 0.000110867738;
    
    fLocationKsi[51]  = 0.078455842803;
    fLocationEta[51]  = 0.078455842803;
    fLocationZeta[51] = 0.769234655050;
    fWeight[51]       = 0.005419989399;
    
    fLocationKsi[52]  = 0.169990521800;
    fLocationEta[52]  = 0.169990521800;
    fLocationZeta[52] = 0.500000000000;
    fWeight[52]       = 0.030243079324;
    
    fLocationKsi[53]  = 0.261525200797;
    fLocationEta[53]  = 0.261525200797;
    fLocationZeta[53] = 0.230765344950;
    fWeight[53]       = 0.060224741049;
    
    fLocationKsi[54]  = 0.324032506649;
    fLocationEta[54]  = 0.324032506649;
    fLocationZeta[54] = 0.046910077050;
    fWeight[54]       = 0.045765756811;
    
    fLocationKsi[55]  = 0.015948536951;
    fLocationEta[55]  = 0.040395970728;
    fLocationZeta[55] = 0.953089922950;
    fWeight[55]       = 0.000059136957;
    
    fLocationKsi[56]  = 0.078455842803;
    fLocationEta[56]  = 0.198720417995;
    fLocationZeta[56] = 0.769234655050;
    fWeight[56]       = 0.002891027494;
    
    fLocationKsi[57]  = 0.169990521800;
    fLocationEta[57]  = 0.430568155800;
    fLocationZeta[57] = 0.500000000000;
    fWeight[57]       = 0.016131687239;
    
    fLocationKsi[58]  = 0.261525200797;
    fLocationEta[58]  = 0.662415893605;
    fLocationZeta[58] = 0.230765344950;
    fWeight[58]       = 0.032123934084;
    
    fLocationKsi[59]  = 0.324032506649;
    fLocationEta[59]  = 0.820740340872;
    fLocationZeta[59] = 0.046910077050;
    fWeight[59]       = 0.024411498156;
    
    fLocationKsi[60]  = 0.040395970728;
    fLocationEta[60]  = -0.040395970728;
    fLocationZeta[60] = 0.953089922950;
    fWeight[60]       = 0.000031543709;
    
    fLocationKsi[61]  = 0.198720417995;
    fLocationEta[61]  = -0.198720417995;
    fLocationZeta[61] = 0.769234655050;
    fWeight[61]       = 0.001542076812;
    
    fLocationKsi[62]  = 0.430568155800;
    fLocationEta[62]  = -0.430568155800;
    fLocationZeta[62] = 0.500000000000;
    fWeight[62]       = 0.008604657297;
    
    fLocationKsi[63]  = 0.662415893605;
    fLocationEta[63]  = -0.662415893605;
    fLocationZeta[63] = 0.230765344950;
    fWeight[63]       = 0.017134936955;
    
    fLocationKsi[64]  = 0.820740340872;
    fLocationEta[64]  = -0.820740340872;
    fLocationZeta[64] = 0.046910077050;
    fWeight[64]       = 0.013021116305;
    
    fLocationKsi[65]  = 0.040395970728;
    fLocationEta[65]  = -0.015948536951;
    fLocationZeta[65] = 0.953089922950;
    fWeight[65]       = 0.000059136957;
    
    fLocationKsi[66]  = 0.198720417995;
    fLocationEta[66]  = -0.078455842803;
    fLocationZeta[66] = 0.769234655050;
    fWeight[66]       = 0.002891027494;
    
    fLocationKsi[67]  = 0.430568155800;
    fLocationEta[67]  = -0.169990521800;
    fLocationZeta[67] = 0.500000000000;
    fWeight[67]       = 0.016131687239;
    
    fLocationKsi[68]  = 0.662415893605;
    fLocationEta[68]  = -0.261525200797;
    fLocationZeta[68] = 0.230765344950;
    fWeight[68]       = 0.032123934084;
    
    fLocationKsi[69]  = 0.820740340872;
    fLocationEta[69]  = -0.324032506649;
    fLocationZeta[69] = 0.046910077050;
    fWeight[69]       = 0.024411498156;
    
    fLocationKsi[70]  = 0.040395970728;
    fLocationEta[70]  = 0.015948536951;
    fLocationZeta[70] = 0.953089922950;
    fWeight[70]       = 0.000059136957;
    
    fLocationKsi[71]  = 0.198720417995;
    fLocationEta[71]  = 0.078455842803;
    fLocationZeta[71] = 0.769234655050;
    fWeight[71]       = 0.002891027494;
    
    fLocationKsi[72]  = 0.430568155800;
    fLocationEta[72]  = 0.169990521800;
    fLocationZeta[72] = 0.500000000000;
    fWeight[72]       = 0.016131687239;
    
    fLocationKsi[73]  = 0.662415893605;
    fLocationEta[73]  = 0.261525200797;
    fLocationZeta[73] = 0.230765344950;
    fWeight[73]       = 0.032123934084;
    
    fLocationKsi[74]  = 0.820740340872;
    fLocationEta[74]  = 0.324032506649;
    fLocationZeta[74] = 0.046910077050;
    fWeight[74]       = 0.024411498156;
    
    fLocationKsi[75]  = 0.040395970728;
    fLocationEta[75]  = 0.040395970728;
    fLocationZeta[75] = 0.953089922950;
    fWeight[75]       = 0.000031543709;
    
    fLocationKsi[76]  = 0.198720417995;
    fLocationEta[76]  = 0.198720417995;
    fLocationZeta[76] = 0.769234655050;
    fWeight[76]       = 0.001542076812;

    fLocationKsi[77]  = 0.430568155800;
    fLocationEta[77]  = 0.430568155800;
    fLocationZeta[77] = 0.500000000000;
    fWeight[77]       = 0.008604657297;
    
    fLocationKsi[78]  = 0.662415893605;
    fLocationEta[78]  = 0.662415893605;
    fLocationZeta[78] = 0.230765344950;
    fWeight[78]       = 0.017134936955;
    
    fLocationKsi[79]  = 0.820740340872;
    fLocationEta[79]  = 0.820740340872;
    fLocationZeta[79] = 0.046910077050;
    fWeight[79]       = 0.013021116305;
    
    break;*/
    
/*  case 150://x8 , x8y4 , x6z6

    break;      */

  /*    case 150://x8 , x8y4 , x6z6

  fLocationKsi[0]  = -0.030597382608;
  fLocationEta[0]  = -0.030597382608;
  fLocationZeta[0] = 0.966234757100;
  fWeight[0]       = 0.000005482238;

  fLocationKsi[1]  = -0.153502613012;
  fLocationEta[1]  = -0.153502613012;
  fLocationZeta[1] = 0.830604693200;
  fWeight[1]       = 0.000290550557;

  fLocationKsi[2]  = -0.344973974351;
  fLocationEta[2]  = -0.344973974351;
  fLocationZeta[2] = 0.619309593000;
  fWeight[2]       = 0.001903305658;

  fLocationKsi[3]  = -0.561205871549;
  fLocationEta[3]  = -0.561205871549;
  fLocationZeta[3] = 0.380690407000;
  fWeight[3]       = 0.005037096680;

  fLocationKsi[4]  = -0.752677232888;
  fLocationEta[4]  = -0.752677232888;
  fLocationZeta[4] = 0.169395306800;
  fWeight[4]       = 0.006985664545;

  fLocationKsi[5]  = -0.875582463292;
  fLocationEta[5]  = -0.875582463292;
  fLocationZeta[5] = 0.033765242900;
  fWeight[5]       = 0.004489349929;

  fLocationKsi[6]  = -0.030597382608;
  fLocationEta[6]  = -0.018181547050;
  fLocationZeta[6] = 0.966234757100;
  fWeight[6]       = 0.000011074963;

  fLocationKsi[7]  = -0.153502613012;
  fLocationEta[7]  = -0.091214173987;
  fLocationZeta[7] = 0.830604693200;
  fWeight[7]       = 0.000586956719;

  fLocationKsi[8]  = -0.344973974351;
  fLocationEta[8]  = -0.204990100819;
  fLocationZeta[8] = 0.619309593000;
  fWeight[8]       = 0.003844969541;

  fLocationKsi[9]  = -0.561205871549;
  fLocationEta[9]  = -0.333479209281;
  fLocationZeta[9] = 0.380690407000;
  fWeight[9]       = 0.010175708369;

  fLocationKsi[10]  = -0.752677232888;
  fLocationEta[10]  = -0.447255136113;
  fLocationZeta[10] = 0.169395306800;
  fWeight[10]       = 0.014112114515;

  fLocationKsi[11]  = -0.875582463292;
  fLocationEta[11]  = -0.520287763050;
  fLocationZeta[11] = 0.033765242900;
  fWeight[11]       = 0.009069175867;

  fLocationKsi[12]  = -0.030597382608;
  fLocationEta[12]  = 0.000000000000;
  fLocationZeta[12] = 0.966234757100;
  fWeight[12]       = 0.000013163489;

  fLocationKsi[13]  = -0.153502613012;
  fLocationEta[13]  = 0.000000000000;
  fLocationZeta[13] = 0.830604693200;
  fWeight[13]       = 0.000697645536;

  fLocationKsi[14]  = -0.344973974351;
  fLocationEta[14]  = 0.000000000000;
  fLocationZeta[14] = 0.619309593000;
  fWeight[14]       = 0.004570057301;

  fLocationKsi[15]  = -0.561205871549;
  fLocationEta[15]  = 0.000000000000;
  fLocationZeta[15] = 0.380690407000;
  fWeight[15]       = 0.012094652462;

  fLocationKsi[16]  = -0.752677232888;
  fLocationEta[16]  = 0.000000000000;
  fLocationZeta[16] = 0.169395306800;
  fWeight[16]       = 0.016773389564;

  fLocationKsi[17]  = -0.875582463292;
  fLocationEta[17]  = 0.000000000000;
  fLocationZeta[17] = 0.033765242900;
  fWeight[17]       = 0.010779449081;

  fLocationKsi[18]  = -0.030597382608;
  fLocationEta[18]  = 0.018181547050;
  fLocationZeta[18] = 0.966234757100;
  fWeight[18]       = 0.000011074963;

  fLocationKsi[19]  = -0.153502613012;
  fLocationEta[19]  = 0.091214173987;
  fLocationZeta[19] = 0.830604693200;
  fWeight[19]       = 0.000586956719;

  fLocationKsi[20]  = -0.344973974351;
  fLocationEta[20]  = 0.204990100819;
  fLocationZeta[20] = 0.619309593000;
  fWeight[20]       = 0.003844969541;

  fLocationKsi[21]  = -0.561205871549;
  fLocationEta[21]  = 0.333479209281;
  fLocationZeta[21] = 0.380690407000;
  fWeight[21]       = 0.010175708369;

  fLocationKsi[22]  = -0.752677232888;
  fLocationEta[22]  = 0.447255136113;
  fLocationZeta[22] = 0.169395306800;
  fWeight[22]       = 0.014112114515;

  fLocationKsi[23]  = -0.875582463292;
  fLocationEta[23]  = 0.520287763050;
  fLocationZeta[23] = 0.033765242900;
  fWeight[23]       = 0.009069175867;

  fLocationKsi[24]  = -0.030597382608;
  fLocationEta[24]  = 0.030597382608;
  fLocationZeta[24] = 0.966234757100;
  fWeight[24]       = 0.000005482238;

  fLocationKsi[25]  = -0.153502613012;
  fLocationEta[25]  = 0.153502613012;
  fLocationZeta[25] = 0.830604693200;
  fWeight[25]       = 0.000290550557;

  fLocationKsi[26]  = -0.344973974351;
  fLocationEta[26]  = 0.344973974351;
  fLocationZeta[26] = 0.619309593000;
  fWeight[26]       = 0.001903305658;

  fLocationKsi[27]  = -0.561205871549;
  fLocationEta[27]  = 0.561205871549;
  fLocationZeta[27] = 0.380690407000;
  fWeight[27]       = 0.005037096680;

  fLocationKsi[28]  = -0.752677232888;
  fLocationEta[28]  = 0.752677232888;
  fLocationZeta[28] = 0.169395306800;
  fWeight[28]       = 0.006985664545;

  fLocationKsi[29]  = -0.875582463292;
  fLocationEta[29]  = 0.875582463292;
  fLocationZeta[29] = 0.033765242900;
  fWeight[29]       = 0.004489349929;

  fLocationKsi[30]  = -0.018181547050;
  fLocationEta[30]  = -0.030597382608;
  fLocationZeta[30] = 0.966234757100;
  fWeight[30]       = 0.000011074963;

  fLocationKsi[31]  = -0.091214173987;
  fLocationEta[31]  = -0.153502613012;
  fLocationZeta[31] = 0.830604693200;
  fWeight[31]       = 0.000586956719;

  fLocationKsi[32]  = -0.204990100819;
  fLocationEta[32]  = -0.344973974351;
  fLocationZeta[32] = 0.619309593000;
  fWeight[32]       = 0.003844969541;

  fLocationKsi[33]  = -0.333479209281;
  fLocationEta[33]  = -0.561205871549;
  fLocationZeta[33] = 0.380690407000;
  fWeight[33]       = 0.010175708369;

  fLocationKsi[34]  = -0.447255136113;
  fLocationEta[34]  = -0.752677232888;
  fLocationZeta[34] = 0.169395306800;
  fWeight[34]       = 0.014112114515;

  fLocationKsi[35]  = -0.520287763050;
  fLocationEta[35]  = -0.875582463292;
  fLocationZeta[35] = 0.033765242900;
  fWeight[35]       = 0.009069175867;

  fLocationKsi[36]  = -0.018181547050;
  fLocationEta[36]  = -0.018181547050;
  fLocationZeta[36] = 0.966234757100;
  fWeight[36]       = 0.000022373124;

  fLocationKsi[37]  = -0.091214173987;
  fLocationEta[37]  = -0.091214173987;
  fLocationZeta[37] = 0.830604693200;
  fWeight[37]       = 0.001185742656;

  fLocationKsi[38]  = -0.204990100819;
  fLocationEta[38]  = -0.204990100819;
  fLocationZeta[38] = 0.619309593000;
  fWeight[38]       = 0.007767428585;

  fLocationKsi[39]  = -0.333479209281;
  fLocationEta[39]  = -0.333479209281;
  fLocationZeta[39] = 0.380690407000;
  fWeight[39]       = 0.020556492637;

  fLocationKsi[40]  = -0.447255136113;
  fLocationEta[40]  = -0.447255136113;
  fLocationZeta[40] = 0.169395306800;
  fWeight[40]       = 0.028508637198;

  fLocationKsi[41]  = -0.520287763050;
  fLocationEta[41]  = -0.520287763050;
  fLocationZeta[41] = 0.033765242900;
  fWeight[41]       = 0.018321127156;

  fLocationKsi[42]  = -0.018181547050;
  fLocationEta[42]  = 0.000000000000;
  fLocationZeta[42] = 0.966234757100;
  fWeight[42]       = 0.000026592268;

  fLocationKsi[43]  = -0.091214173987;
  fLocationEta[43]  = 0.000000000000;
  fLocationZeta[43] = 0.830604693200;
  fWeight[43]       = 0.001409351055;

  fLocationKsi[44]  = -0.204990100819;
  fLocationEta[44]  = 0.000000000000;
  fLocationZeta[44] = 0.619309593000;
  fWeight[44]       = 0.009232217144;

  fLocationKsi[45]  = -0.333479209281;
  fLocationEta[45]  = 0.000000000000;
  fLocationZeta[45] = 0.380690407000;
  fWeight[45]       = 0.024433054218;

  fLocationKsi[46]  = -0.447255136113;
  fLocationEta[46]  = 0.000000000000;
  fLocationZeta[46] = 0.169395306800;
  fWeight[46]       = 0.033884821240;

  fLocationKsi[47]  = -0.520287763050;
  fLocationEta[47]  = 0.000000000000;
  fLocationZeta[47] = 0.033765242900;
  fWeight[47]       = 0.021776141535;

  fLocationKsi[48]  = -0.018181547050;
  fLocationEta[48]  = 0.018181547050;
  fLocationZeta[48] = 0.966234757100;
  fWeight[48]       = 0.000022373124;

  fLocationKsi[49]  = -0.091214173987;
  fLocationEta[49]  = 0.091214173987;
  fLocationZeta[49] = 0.830604693200;
  fWeight[49]       = 0.001185742656;

  fLocationKsi[50]  = -0.204990100819;
  fLocationEta[50]  = 0.204990100819;
  fLocationZeta[50] = 0.619309593000;
  fWeight[50]       = 0.007767428585;

  fLocationKsi[51]  = -0.333479209281;
  fLocationEta[51]  = 0.333479209281;
  fLocationZeta[51] = 0.380690407000;
  fWeight[51]       = 0.020556492637;

  fLocationKsi[52]  = -0.447255136113;
  fLocationEta[52]  = 0.447255136113;
  fLocationZeta[52] = 0.169395306800;
  fWeight[52]       = 0.028508637198;

  fLocationKsi[53]  = -0.520287763050;
  fLocationEta[53]  = 0.520287763050;
  fLocationZeta[53] = 0.033765242900;
  fWeight[53]       = 0.018321127156;

  fLocationKsi[54]  = -0.018181547050;
  fLocationEta[54]  = 0.030597382608;
  fLocationZeta[54] = 0.966234757100;
  fWeight[54]       = 0.000011074963;

  fLocationKsi[55]  = -0.091214173987;
  fLocationEta[55]  = 0.153502613012;
  fLocationZeta[55] = 0.830604693200;
  fWeight[55]       = 0.000586956719;

  fLocationKsi[56]  = -0.204990100819;
  fLocationEta[56]  = 0.344973974351;
  fLocationZeta[56] = 0.619309593000;
  fWeight[56]       = 0.003844969541;

  fLocationKsi[57]  = -0.333479209281;
  fLocationEta[57]  = 0.561205871549;
  fLocationZeta[57] = 0.380690407000;
  fWeight[57]       = 0.010175708369;

  fLocationKsi[58]  = -0.447255136113;
  fLocationEta[58]  = 0.752677232888;
  fLocationZeta[58] = 0.169395306800;
  fWeight[58]       = 0.014112114515;

  fLocationKsi[59]  = -0.520287763050;
  fLocationEta[59]  = 0.875582463292;
  fLocationZeta[59] = 0.033765242900;
  fWeight[59]       = 0.009069175867;

  fLocationKsi[60]  = 0.000000000000;
  fLocationEta[60]  = -0.030597382608;
  fLocationZeta[60] = 0.966234757100;
  fWeight[60]       = 0.000013163489;

  fLocationKsi[61]  = 0.000000000000;
  fLocationEta[61]  = -0.153502613012;
  fLocationZeta[61] = 0.830604693200;
  fWeight[61]       = 0.000697645536;

  fLocationKsi[62]  = 0.000000000000;
  fLocationEta[62]  = -0.344973974351;
  fLocationZeta[62] = 0.619309593000;
  fWeight[62]       = 0.004570057301;

  fLocationKsi[63]  = 0.000000000000;
  fLocationEta[63]  = -0.561205871549;
  fLocationZeta[63] = 0.380690407000;
  fWeight[63]       = 0.012094652462;

  fLocationKsi[64]  = 0.000000000000;
  fLocationEta[64]  = -0.752677232888;
  fLocationZeta[64] = 0.169395306800;
  fWeight[64]       = 0.016773389564;

  fLocationKsi[65]  = 0.000000000000;
  fLocationEta[65]  = -0.875582463292;
  fLocationZeta[65] = 0.033765242900;
  fWeight[65]       = 0.010779449081;

  fLocationKsi[66]  = 0.000000000000;
  fLocationEta[66]  = -0.018181547050;
  fLocationZeta[66] = 0.966234757100;
  fWeight[66]       = 0.000026592268;

  fLocationKsi[67]  = 0.000000000000;
  fLocationEta[67]  = -0.091214173987;
  fLocationZeta[67] = 0.830604693200;
  fWeight[67]       = 0.001409351055;

  fLocationKsi[68]  = 0.000000000000;
  fLocationEta[68]  = -0.204990100819;
  fLocationZeta[68] = 0.619309593000;
  fWeight[68]       = 0.009232217144;

  fLocationKsi[69]  = 0.000000000000;
  fLocationEta[69]  = -0.333479209281;
  fLocationZeta[69] = 0.380690407000;
  fWeight[69]       = 0.024433054218;

  fLocationKsi[70]  = 0.000000000000;
  fLocationEta[70]  = -0.447255136113;
  fLocationZeta[70] = 0.169395306800;
  fWeight[70]       = 0.033884821240;

  fLocationKsi[71]  = 0.000000000000;
  fLocationEta[71]  = -0.520287763050;
  fLocationZeta[71] = 0.033765242900;
  fWeight[71]       = 0.021776141535;

  fLocationKsi[72]  = 0.000000000000;
  fLocationEta[72]  = 0.000000000000;
  fLocationZeta[72] = 0.966234757100;
  fWeight[72]       = 0.000031607061;

  fLocationKsi[73]  = 0.000000000000;
  fLocationEta[73]  = 0.000000000000;
  fLocationZeta[73] = 0.830604693200;
  fWeight[73]       = 0.001675127725;

  fLocationKsi[74]  = 0.000000000000;
  fLocationEta[74]  = 0.000000000000;
  fLocationZeta[74] = 0.619309593000;
  fWeight[74]       = 0.010973236827;

  fLocationKsi[75]  = 0.000000000000;
  fLocationEta[75]  = 0.000000000000;
  fLocationZeta[75] = 0.380690407000;
  fWeight[75]       = 0.029040661213;

  fLocationKsi[76]  = 0.000000000000;
  fLocationEta[76]  = 0.000000000000;
  fLocationZeta[76] = 0.169395306800;
  fWeight[76]       = 0.040274850828;

  fLocationKsi[77]  = 0.000000000000;
  fLocationEta[77]  = 0.000000000000;
  fLocationZeta[77] = 0.033765242900;
  fWeight[77]       = 0.025882705584;

  fLocationKsi[78]  = 0.000000000000;
  fLocationEta[78]  = 0.018181547050;
  fLocationZeta[78] = 0.966234757100;
  fWeight[78]       = 0.000026592268;

  fLocationKsi[79]  = 0.000000000000;
  fLocationEta[79]  = 0.091214173987;
  fLocationZeta[79] = 0.830604693200;
  fWeight[79]       = 0.001409351055;

  fLocationKsi[80]  = 0.000000000000;
  fLocationEta[80]  = 0.204990100819;
  fLocationZeta[80] = 0.619309593000;
  fWeight[80]       = 0.009232217144;

  fLocationKsi[81]  = 0.000000000000;
  fLocationEta[81]  = 0.333479209281;
  fLocationZeta[81] = 0.380690407000;
  fWeight[81]       = 0.024433054218;

  fLocationKsi[82]  = 0.000000000000;
  fLocationEta[82]  = 0.447255136113;
  fLocationZeta[82] = 0.169395306800;
  fWeight[82]       = 0.033884821240;

  fLocationKsi[83]  = 0.000000000000;
  fLocationEta[83]  = 0.520287763050;
  fLocationZeta[83] = 0.033765242900;
  fWeight[83]       = 0.021776141535;

  fLocationKsi[84]  = 0.000000000000;
  fLocationEta[84]  = 0.030597382608;
  fLocationZeta[84] = 0.966234757100;
  fWeight[84]       = 0.000013163489;

  fLocationKsi[85]  = 0.000000000000;
  fLocationEta[85]  = 0.153502613012;
  fLocationZeta[85] = 0.830604693200;
  fWeight[85]       = 0.000697645536;

  fLocationKsi[86]  = 0.000000000000;
  fLocationEta[86]  = 0.344973974351;
  fLocationZeta[86] = 0.619309593000;
  fWeight[86]       = 0.004570057301;

  fLocationKsi[87]  = 0.000000000000;
  fLocationEta[87]  = 0.561205871549;
  fLocationZeta[87] = 0.380690407000;
  fWeight[87]       = 0.012094652462;

  fLocationKsi[88]  = 0.000000000000;
  fLocationEta[88]  = 0.752677232888;
  fLocationZeta[88] = 0.169395306800;
  fWeight[88]       = 0.016773389564;

  fLocationKsi[89]  = 0.000000000000;
  fLocationEta[89]  = 0.875582463292;
  fLocationZeta[89] = 0.033765242900;
  fWeight[89]       = 0.010779449081;

  fLocationKsi[90]  = 0.018181547050;
  fLocationEta[90]  = -0.030597382608;
  fLocationZeta[90] = 0.966234757100;
  fWeight[90]       = 0.000011074963;

  fLocationKsi[91]  = 0.091214173987;
  fLocationEta[91]  = -0.153502613012;
  fLocationZeta[91] = 0.830604693200;
  fWeight[91]       = 0.000586956719;

  fLocationKsi[92]  = 0.204990100819;
  fLocationEta[92]  = -0.344973974351;
  fLocationZeta[92] = 0.619309593000;
  fWeight[92]       = 0.003844969541;

  fLocationKsi[93]  = 0.333479209281;
  fLocationEta[93]  = -0.561205871549;
  fLocationZeta[93] = 0.380690407000;
  fWeight[93]       = 0.010175708369;

  fLocationKsi[94]  = 0.447255136113;
  fLocationEta[94]  = -0.752677232888;
  fLocationZeta[94] = 0.169395306800;
  fWeight[94]       = 0.014112114515;

  fLocationKsi[95]  = 0.520287763050;
  fLocationEta[95]  = -0.875582463292;
  fLocationZeta[95] = 0.033765242900;
  fWeight[95]       = 0.009069175867;

  fLocationKsi[96]  = 0.018181547050;
  fLocationEta[96]  = -0.018181547050;
  fLocationZeta[96] = 0.966234757100;
  fWeight[96]       = 0.000022373124;

  fLocationKsi[97]  = 0.091214173987;
  fLocationEta[97]  = -0.091214173987;
  fLocationZeta[97] = 0.830604693200;
  fWeight[97]       = 0.001185742656;

  fLocationKsi[98]  = 0.204990100819;
  fLocationEta[98]  = -0.204990100819;
  fLocationZeta[98] = 0.619309593000;
  fWeight[98]       = 0.007767428585;

  fLocationKsi[99]  = 0.333479209281;
  fLocationEta[99]  = -0.333479209281;
  fLocationZeta[99] = 0.380690407000;
  fWeight[99]       = 0.020556492637;

  fLocationKsi[100]  = 0.447255136113;
  fLocationEta[100]  = -0.447255136113;
  fLocationZeta[100] = 0.169395306800;
  fWeight[100]       = 0.028508637198;

  fLocationKsi[101]  = 0.520287763050;
  fLocationEta[101]  = -0.520287763050;
  fLocationZeta[101] = 0.033765242900;
  fWeight[101]       = 0.018321127156;

  fLocationKsi[102]  = 0.018181547050;
  fLocationEta[102]  = 0.000000000000;
  fLocationZeta[102] = 0.966234757100;
  fWeight[102]       = 0.000026592268;

  fLocationKsi[103]  = 0.091214173987;
  fLocationEta[103]  = 0.000000000000;
  fLocationZeta[103] = 0.830604693200;
  fWeight[103]       = 0.001409351055;

  fLocationKsi[104]  = 0.204990100819;
  fLocationEta[104]  = 0.000000000000;
  fLocationZeta[104] = 0.619309593000;
  fWeight[104]       = 0.009232217144;

  fLocationKsi[105]  = 0.333479209281;
  fLocationEta[105]  = 0.000000000000;
  fLocationZeta[105] = 0.380690407000;
  fWeight[105]       = 0.024433054218;

  fLocationKsi[106]  = 0.447255136113;
  fLocationEta[106]  = 0.000000000000;
  fLocationZeta[106] = 0.169395306800;
  fWeight[106]       = 0.033884821240;

  fLocationKsi[107]  = 0.520287763050;
  fLocationEta[107]  = 0.000000000000;
  fLocationZeta[107] = 0.033765242900;
  fWeight[107]       = 0.021776141535;

  fLocationKsi[108]  = 0.018181547050;
  fLocationEta[108]  = 0.018181547050;
  fLocationZeta[108] = 0.966234757100;
  fWeight[108]       = 0.000022373124;

  fLocationKsi[109]  = 0.091214173987;
  fLocationEta[109]  = 0.091214173987;
  fLocationZeta[109] = 0.830604693200;
  fWeight[109]       = 0.001185742656;

  fLocationKsi[110]  = 0.204990100819;
  fLocationEta[110]  = 0.204990100819;
  fLocationZeta[110] = 0.619309593000;
  fWeight[110]       = 0.007767428585;

  fLocationKsi[111]  = 0.333479209281;
  fLocationEta[111]  = 0.333479209281;
  fLocationZeta[111] = 0.380690407000;
  fWeight[111]       = 0.020556492637;

  fLocationKsi[112]  = 0.447255136113;
  fLocationEta[112]  = 0.447255136113;
  fLocationZeta[112] = 0.169395306800;
  fWeight[112]       = 0.028508637198;

  fLocationKsi[113]  = 0.520287763050;
  fLocationEta[113]  = 0.520287763050;
  fLocationZeta[113] = 0.033765242900;
  fWeight[113]       = 0.018321127156;

  fLocationKsi[114]  = 0.018181547050;
  fLocationEta[114]  = 0.030597382608;
  fLocationZeta[114] = 0.966234757100;
  fWeight[114]       = 0.000011074963;

  fLocationKsi[115]  = 0.091214173987;
  fLocationEta[115]  = 0.153502613012;
  fLocationZeta[115] = 0.830604693200;
  fWeight[115]       = 0.000586956719;

  fLocationKsi[116]  = 0.204990100819;
  fLocationEta[116]  = 0.344973974351;
  fLocationZeta[116] = 0.619309593000;
  fWeight[116]       = 0.003844969541;

  fLocationKsi[117]  = 0.333479209281;
  fLocationEta[117]  = 0.561205871549;
  fLocationZeta[117] = 0.380690407000;
  fWeight[117]       = 0.010175708369;

  fLocationKsi[118]  = 0.447255136113;
  fLocationEta[118]  = 0.752677232888;
  fLocationZeta[118] = 0.169395306800;
  fWeight[118]       = 0.014112114515;

  fLocationKsi[119]  = 0.520287763050;
  fLocationEta[119]  = 0.875582463292;
  fLocationZeta[119] = 0.033765242900;
  fWeight[119]       = 0.009069175867;

  fLocationKsi[120]  = 0.030597382608;
  fLocationEta[120]  = -0.030597382608;
  fLocationZeta[120] = 0.966234757100;
  fWeight[120]       = 0.000005482238;

  fLocationKsi[121]  = 0.153502613012;
  fLocationEta[121]  = -0.153502613012;
  fLocationZeta[121] = 0.830604693200;
  fWeight[121]       = 0.000290550557;

  fLocationKsi[122]  = 0.344973974351;
  fLocationEta[122]  = -0.344973974351;
  fLocationZeta[122] = 0.619309593000;
  fWeight[122]       = 0.001903305658;

  fLocationKsi[123]  = 0.561205871549;
  fLocationEta[123]  = -0.561205871549;
  fLocationZeta[123] = 0.380690407000;
  fWeight[123]       = 0.005037096680;

  fLocationKsi[124]  = 0.752677232888;
  fLocationEta[124]  = -0.752677232888;
  fLocationZeta[124] = 0.169395306800;
  fWeight[124]       = 0.006985664545;

  fLocationKsi[125]  = 0.875582463292;
  fLocationEta[125]  = -0.875582463292;
  fLocationZeta[125] = 0.033765242900;
  fWeight[125]       = 0.004489349929;

  fLocationKsi[126]  = 0.030597382608;
  fLocationEta[126]  = -0.018181547050;
  fLocationZeta[126] = 0.966234757100;
  fWeight[126]       = 0.000011074963;

  fLocationKsi[127]  = 0.153502613012;
  fLocationEta[127]  = -0.091214173987;
  fLocationZeta[127] = 0.830604693200;
  fWeight[127]       = 0.000586956719;

  fLocationKsi[128]  = 0.344973974351;
  fLocationEta[128]  = -0.204990100819;
  fLocationZeta[128] = 0.619309593000;
  fWeight[128]       = 0.003844969541;

  fLocationKsi[129]  = 0.561205871549;
  fLocationEta[129]  = -0.333479209281;
  fLocationZeta[129] = 0.380690407000;
  fWeight[129]       = 0.010175708369;

  fLocationKsi[130]  = 0.752677232888;
  fLocationEta[130]  = -0.447255136113;
  fLocationZeta[130] = 0.169395306800;
  fWeight[130]       = 0.014112114515;

  fLocationKsi[131]  = 0.875582463292;
  fLocationEta[131]  = -0.520287763050;
  fLocationZeta[131] = 0.033765242900;
  fWeight[131]       = 0.009069175867;

  fLocationKsi[132]  = 0.030597382608;
  fLocationEta[132]  = 0.000000000000;
  fLocationZeta[132] = 0.966234757100;
  fWeight[132]       = 0.000013163489;

  fLocationKsi[133]  = 0.153502613012;
  fLocationEta[133]  = 0.000000000000;
  fLocationZeta[133] = 0.830604693200;
  fWeight[133]       = 0.000697645536;

  fLocationKsi[134]  = 0.344973974351;
  fLocationEta[134]  = 0.000000000000;
  fLocationZeta[134] = 0.619309593000;
  fWeight[134]       = 0.004570057301;

  fLocationKsi[135]  = 0.561205871549;
  fLocationEta[135]  = 0.000000000000;
  fLocationZeta[135] = 0.380690407000;
  fWeight[135]       = 0.012094652462;

  fLocationKsi[136]  = 0.752677232888;
  fLocationEta[136]  = 0.000000000000;
  fLocationZeta[136] = 0.169395306800;
  fWeight[136]       = 0.016773389564;

  fLocationKsi[137]  = 0.875582463292;
  fLocationEta[137]  = 0.000000000000;
  fLocationZeta[137] = 0.033765242900;
  fWeight[137]       = 0.010779449081;

  fLocationKsi[138]  = 0.030597382608;
  fLocationEta[138]  = 0.018181547050;
  fLocationZeta[138] = 0.966234757100;
  fWeight[138]       = 0.000011074963;

  fLocationKsi[139]  = 0.153502613012;
  fLocationEta[139]  = 0.091214173987;
  fLocationZeta[139] = 0.830604693200;
  fWeight[139]       = 0.000586956719;

  fLocationKsi[140]  = 0.344973974351;
  fLocationEta[140]  = 0.204990100819;
  fLocationZeta[140] = 0.619309593000;
  fWeight[140]       = 0.003844969541;

  fLocationKsi[141]  = 0.561205871549;
  fLocationEta[141]  = 0.333479209281;
  fLocationZeta[141] = 0.380690407000;
  fWeight[141]       = 0.010175708369;

  fLocationKsi[142]  = 0.752677232888;
  fLocationEta[142]  = 0.447255136113;
  fLocationZeta[142] = 0.169395306800;
  fWeight[142]       = 0.014112114515;

  fLocationKsi[143]  = 0.875582463292;
  fLocationEta[143]  = 0.520287763050;
  fLocationZeta[143] = 0.033765242900;
  fWeight[143]       = 0.009069175867;

  fLocationKsi[144]  = 0.030597382608;
  fLocationEta[144]  = 0.030597382608;
  fLocationZeta[144] = 0.966234757100;
  fWeight[144]       = 0.000005482238;

  fLocationKsi[145]  = 0.153502613012;
  fLocationEta[145]  = 0.153502613012;
  fLocationZeta[145] = 0.830604693200;
  fWeight[145]       = 0.000290550557;

  fLocationKsi[146]  = 0.344973974351;
  fLocationEta[146]  = 0.344973974351;
  fLocationZeta[146] = 0.619309593000;
  fWeight[146]       = 0.001903305658;

  fLocationKsi[147]  = 0.561205871549;
  fLocationEta[147]  = 0.561205871549;
  fLocationZeta[147] = 0.380690407000;
  fWeight[147]       = 0.005037096680;

  fLocationKsi[148]  = 0.752677232888;
  fLocationEta[148]  = 0.752677232888;
  fLocationZeta[148] = 0.169395306800;
  fWeight[148]       = 0.006985664545;

  fLocationKsi[149]  = 0.875582463292;
  fLocationEta[149]  = 0.875582463292;
  fLocationZeta[149] = 0.033765242900;
  fWeight[149]       = 0.004489349929;

  break;        */

  case 96:
    fLocationKsi[0]  = -0.029076476731;
    fLocationEta[0]  = -0.029076476731;
    fLocationZeta[0] = 0.966234757100;
    fWeight[0]       = 0.000011817492;
    
    fLocationKsi[1]  = -0.145872449700;
    fLocationEta[1]  = -0.145872449700;
    fLocationZeta[1] = 0.830604693200;
    fWeight[1]       = 0.000626309699;
    
    fLocationKsi[2]  = -0.327826332945;
    fLocationEta[2]  = -0.327826332945;
    fLocationZeta[2] = 0.619309593000;
    fWeight[2]       = 0.004102758586;
    
    fLocationKsi[3]  = -0.533309978655;
    fLocationEta[3]  = -0.533309978655;
    fLocationZeta[3] = 0.380690407000;
    fWeight[3]       = 0.010857946841;
    
    fLocationKsi[4]  = -0.715263861900;
    fLocationEta[4]  = -0.715263861900;
    fLocationZeta[4] = 0.169395306800;
    fWeight[4]       = 0.015058272472;
    
    fLocationKsi[5]  = -0.832059834869;
    fLocationEta[5]  = -0.832059834869;
    fLocationZeta[5] = 0.033765242900;
    fWeight[5]       = 0.009677225986;
    
    fLocationKsi[6]  = -0.029076476731;
    fLocationEta[6]  = -0.011479542519;
    fLocationZeta[6] = 0.966234757100;
    fWeight[6]       = 0.000022154989;
    
    fLocationKsi[7]  = -0.145872449700;
    fLocationEta[7]  = -0.057591193187;
    fLocationZeta[7] = 0.830604693200;
    fWeight[7]       = 0.001174181822;
    
    fLocationKsi[8]  = -0.327826332945;
    fLocationEta[8]  = -0.129427521860;
    fLocationZeta[8] = 0.619309593000;
    fWeight[8]       = 0.007691697187;
    
    fLocationKsi[9]  = -0.533309978655;
    fLocationEta[9]  = -0.210553521740;
    fLocationZeta[9] = 0.380690407000;
    fWeight[9]       = 0.020356069562;
    
    fLocationKsi[10]  = -0.715263861900;
    fLocationEta[10]  = -0.282389850413;
    fLocationZeta[10] = 0.169395306800;
    fWeight[10]       = 0.028230681769;
    
    fLocationKsi[11]  = -0.832059834869;
    fLocationEta[11]  = -0.328501501081;
    fLocationZeta[11] = 0.033765242900;
    fWeight[11]       = 0.018142498598;
    
    fLocationKsi[12]  = -0.029076476731;
    fLocationEta[12]  = 0.011479542519;
    fLocationZeta[12] = 0.966234757100;
    fWeight[12]       = 0.000022154989;
    
    fLocationKsi[13]  = -0.145872449700;
    fLocationEta[13]  = 0.057591193187;
    fLocationZeta[13] = 0.830604693200;
    fWeight[13]       = 0.001174181822;
    
    fLocationKsi[14]  = -0.327826332945;
    fLocationEta[14]  = 0.129427521860;
    fLocationZeta[14] = 0.619309593000;
    fWeight[14]       = 0.007691697187;
    
    fLocationKsi[15]  = -0.533309978655;
    fLocationEta[15]  = 0.210553521740;
    fLocationZeta[15] = 0.380690407000;
    fWeight[15]       = 0.020356069562;
    
    fLocationKsi[16]  = -0.715263861900;
    fLocationEta[16]  = 0.282389850413;
    fLocationZeta[16] = 0.169395306800;
    fWeight[16]       = 0.028230681769;
    
    fLocationKsi[17]  = -0.832059834869;
    fLocationEta[17]  = 0.328501501081;
    fLocationZeta[17] = 0.033765242900;
    fWeight[17]       = 0.018142498598;
    
    fLocationKsi[18]  = -0.029076476731;
    fLocationEta[18]  = 0.029076476731;
    fLocationZeta[18] = 0.966234757100;
    fWeight[18]       = 0.000011817492;
    
    fLocationKsi[19]  = -0.145872449700;
    fLocationEta[19]  = 0.145872449700;
    fLocationZeta[19] = 0.830604693200;
    fWeight[19]       = 0.000626309699;
    
    fLocationKsi[20]  = -0.327826332945;
    fLocationEta[20]  = 0.327826332945;
    fLocationZeta[20] = 0.619309593000;
    fWeight[20]       = 0.004102758586;
    
    fLocationKsi[21]  = -0.533309978655;
    fLocationEta[21]  = 0.533309978655;
    fLocationZeta[21] = 0.380690407000;
    fWeight[21]       = 0.010857946841;
    
    fLocationKsi[22]  = -0.715263861900;
    fLocationEta[22]  = 0.715263861900;
    fLocationZeta[22] = 0.169395306800;
    fWeight[22]       = 0.015058272472;
    
    fLocationKsi[23]  = -0.832059834869;
    fLocationEta[23]  = 0.832059834869;
    fLocationZeta[23] = 0.033765242900;
    fWeight[23]       = 0.009677225986;
    
    fLocationKsi[24]  = -0.011479542519;
    fLocationEta[24]  = -0.029076476731;
    fLocationZeta[24] = 0.966234757100;
    fWeight[24]       = 0.000022154989;
    
    fLocationKsi[25]  = -0.057591193187;
    fLocationEta[25]  = -0.145872449700;
    fLocationZeta[25] = 0.830604693200;
    fWeight[25]       = 0.001174181822;
    
    fLocationKsi[26]  = -0.129427521860;
    fLocationEta[26]  = -0.327826332945;
    fLocationZeta[26] = 0.619309593000;
    fWeight[26]       = 0.007691697187;
    
    fLocationKsi[27]  = -0.210553521740;
    fLocationEta[27]  = -0.533309978655;
    fLocationZeta[27] = 0.380690407000;
    fWeight[27]       = 0.020356069562;
    
    fLocationKsi[28]  = -0.282389850413;
    fLocationEta[28]  = -0.715263861900;
    fLocationZeta[28] = 0.169395306800;
    fWeight[28]       = 0.028230681769;
    
    fLocationKsi[29]  = -0.328501501081;
    fLocationEta[29]  = -0.832059834869;
    fLocationZeta[29] = 0.033765242900;
    fWeight[29]       = 0.018142498598;
    
    fLocationKsi[30]  = -0.011479542519;
    fLocationEta[30]  = -0.011479542519;
    fLocationZeta[30] = 0.966234757100;
    fWeight[30]       = 0.000041535339;
    
    fLocationKsi[31]  = -0.057591193187;
    fLocationEta[31]  = -0.057591193187;
    fLocationZeta[31] = 0.830604693200;
    fWeight[31]       = 0.002201311832;
    
    fLocationKsi[32]  = -0.129427521860;
    fLocationEta[32]  = -0.129427521860;
    fLocationZeta[32] = 0.619309593000;
    fWeight[32]       = 0.014420104030;
    
    fLocationKsi[33]  = -0.210553521740;
    fLocationEta[33]  = -0.210553521740;
    fLocationZeta[33] = 0.380690407000;
    fWeight[33]       = 0.038162792109;
    
    fLocationKsi[34]  = -0.282389850413;
    fLocationEta[34]  = -0.282389850413;
    fLocationZeta[34] = 0.169395306800;
    fWeight[34]       = 0.052925818325;
    
    fLocationKsi[35]  = -0.328501501081;
    fLocationEta[35]  = -0.328501501081;
    fLocationZeta[35] = 0.033765242900;
    fWeight[35]       = 0.034012872682;
    
    fLocationKsi[36]  = -0.011479542519;
    fLocationEta[36]  = 0.011479542519;
    fLocationZeta[36] = 0.966234757100;
    fWeight[36]       = 0.000041535339;
    
    fLocationKsi[37]  = -0.057591193187;
    fLocationEta[37]  = 0.057591193187;
    fLocationZeta[37] = 0.830604693200;
    fWeight[37]       = 0.002201311832;
    
    fLocationKsi[38]  = -0.129427521860;
    fLocationEta[38]  = 0.129427521860;
    fLocationZeta[38] = 0.619309593000;
    fWeight[38]       = 0.014420104030;
    
    fLocationKsi[39]  = -0.210553521740;
    fLocationEta[39]  = 0.210553521740;
    fLocationZeta[39] = 0.380690407000;
    fWeight[39]       = 0.038162792109;
    
    fLocationKsi[40]  = -0.282389850413;
    fLocationEta[40]  = 0.282389850413;
    fLocationZeta[40] = 0.169395306800;
    fWeight[40]       = 0.052925818325;
    
    fLocationKsi[41]  = -0.328501501081;
    fLocationEta[41]  = 0.328501501081;
    fLocationZeta[41] = 0.033765242900;
    fWeight[41]       = 0.034012872682;
    
    fLocationKsi[42]  = -0.011479542519;
    fLocationEta[42]  = 0.029076476731;
    fLocationZeta[42] = 0.966234757100;
    fWeight[42]       = 0.000022154989;
    
    fLocationKsi[43]  = -0.057591193187;
    fLocationEta[43]  = 0.145872449700;
    fLocationZeta[43] = 0.830604693200;
    fWeight[43]       = 0.001174181822;
    
    fLocationKsi[44]  = -0.129427521860;
    fLocationEta[44]  = 0.327826332945;
    fLocationZeta[44] = 0.619309593000;
    fWeight[44]       = 0.007691697187;
    
    fLocationKsi[45]  = -0.210553521740;
    fLocationEta[45]  = 0.533309978655;
    fLocationZeta[45] = 0.380690407000;
    fWeight[45]       = 0.020356069562;
    
    fLocationKsi[46]  = -0.282389850413;
    fLocationEta[46]  = 0.715263861900;
    fLocationZeta[46] = 0.169395306800;
    fWeight[46]       = 0.028230681769;
    
    fLocationKsi[47]  = -0.328501501081;
    fLocationEta[47]  = 0.832059834869;
    fLocationZeta[47] = 0.033765242900;
    fWeight[47]       = 0.018142498598;
    
    fLocationKsi[48]  = 0.011479542519;
    fLocationEta[48]  = -0.029076476731;
    fLocationZeta[48] = 0.966234757100;
    fWeight[48]       = 0.000022154989;
    
    fLocationKsi[49]  = 0.057591193187;
    fLocationEta[49]  = -0.145872449700;
    fLocationZeta[49] = 0.830604693200;
    fWeight[49]       = 0.001174181822;
    
    fLocationKsi[50]  = 0.129427521860;
    fLocationEta[50]  = -0.327826332945;
    fLocationZeta[50] = 0.619309593000;
    fWeight[50]       = 0.007691697187;
    
    fLocationKsi[51]  = 0.210553521740;
    fLocationEta[51]  = -0.533309978655;
    fLocationZeta[51] = 0.380690407000;
    fWeight[51]       = 0.020356069562;

    fLocationKsi[52]  = 0.282389850413;
    fLocationEta[52]  = -0.715263861900;
    fLocationZeta[52] = 0.169395306800;
    fWeight[52]       = 0.028230681769;
					  
    fLocationKsi[53]  = 0.328501501081;
    fLocationEta[53]  = -0.832059834869;
    fLocationZeta[53] = 0.033765242900;
    fWeight[53]       = 0.018142498598;
    
    fLocationKsi[54]  = 0.011479542519;
    fLocationEta[54]  = -0.011479542519;
    fLocationZeta[54] = 0.966234757100;
    fWeight[54]       = 0.000041535339;
    
    fLocationKsi[55]  = 0.057591193187;
    fLocationEta[55]  = -0.057591193187;
    fLocationZeta[55] = 0.830604693200;
    fWeight[55]       = 0.002201311832;
    
    fLocationKsi[56]  = 0.129427521860;
    fLocationEta[56]  = -0.129427521860;
    fLocationZeta[56] = 0.619309593000;
    fWeight[56]       = 0.014420104030;
    
    fLocationKsi[57]  = 0.210553521740;
    fLocationEta[57]  = -0.210553521740;
    fLocationZeta[57] = 0.380690407000;
    fWeight[57]       = 0.038162792109;
    
    fLocationKsi[58]  = 0.282389850413;
    fLocationEta[58]  = -0.282389850413;
    fLocationZeta[58] = 0.169395306800;
    fWeight[58]       = 0.052925818325;
    
    fLocationKsi[59]  = 0.328501501081;
    fLocationEta[59]  = -0.328501501081;
    fLocationZeta[59] = 0.033765242900;
    fWeight[59]       = 0.034012872682;
    
    fLocationKsi[60]  = 0.011479542519;
    fLocationEta[60]  = 0.011479542519;
    fLocationZeta[60] = 0.966234757100;
    fWeight[60]       = 0.000041535339;
    
    fLocationKsi[61]  = 0.057591193187;
    fLocationEta[61]  = 0.057591193187;
    fLocationZeta[61] = 0.830604693200;
    fWeight[61]       = 0.002201311832;
    
    fLocationKsi[62]  = 0.129427521860;
    fLocationEta[62]  = 0.129427521860;
    fLocationZeta[62] = 0.619309593000;
    fWeight[62]       = 0.014420104030;
    
    fLocationKsi[63]  = 0.210553521740;
    fLocationEta[63]  = 0.210553521740;
    fLocationZeta[63] = 0.380690407000;
    fWeight[63]       = 0.038162792109;
    
    fLocationKsi[64]  = 0.282389850413;
    fLocationEta[64]  = 0.282389850413;
    fLocationZeta[64] = 0.169395306800;
    fWeight[64]       = 0.052925818325;
    
    fLocationKsi[65]  = 0.328501501081;
    fLocationEta[65]  = 0.328501501081;
    fLocationZeta[65] = 0.033765242900;
    fWeight[65]       = 0.034012872682;
    
    fLocationKsi[66]  = 0.011479542519;
    fLocationEta[66]  = 0.029076476731;
    fLocationZeta[66] = 0.966234757100;
    fWeight[66]       = 0.000022154989;
    
    fLocationKsi[67]  = 0.057591193187;
    fLocationEta[67]  = 0.145872449700;
    fLocationZeta[67] = 0.830604693200;
    fWeight[67]       = 0.001174181822;
    
    fLocationKsi[68]  = 0.129427521860;
    fLocationEta[68]  = 0.327826332945;
    fLocationZeta[68] = 0.619309593000;
    fWeight[68]       = 0.007691697187;
    
    fLocationKsi[69]  = 0.210553521740;
    fLocationEta[69]  = 0.533309978655;
    fLocationZeta[69] = 0.380690407000;
    fWeight[69]       = 0.020356069562;
    
    fLocationKsi[70]  = 0.282389850413;
    fLocationEta[70]  = 0.715263861900;
    fLocationZeta[70] = 0.169395306800;
    fWeight[70]       = 0.028230681769;

    fLocationKsi[71]  = 0.328501501081;
    fLocationEta[71]  = 0.832059834869;
    fLocationZeta[71] = 0.033765242900;
    fWeight[71]       = 0.018142498598;
    
    fLocationKsi[72]  = 0.029076476731;
    fLocationEta[72]  = -0.029076476731;
    fLocationZeta[72] = 0.966234757100;
    fWeight[72]       = 0.000011817492;
    
    fLocationKsi[73]  = 0.145872449700;
    fLocationEta[73]  = -0.145872449700;
    fLocationZeta[73] = 0.830604693200;
    fWeight[73]       = 0.000626309699;
    
    fLocationKsi[74]  = 0.327826332945;
    fLocationEta[74]  = -0.327826332945;
    fLocationZeta[74] = 0.619309593000;
    fWeight[74]       = 0.004102758586;
    
    fLocationKsi[75]  = 0.533309978655;
    fLocationEta[75]  = -0.533309978655;
    fLocationZeta[75] = 0.380690407000;
    fWeight[75]       = 0.010857946841;
    
    fLocationKsi[76]  = 0.715263861900;
    fLocationEta[76]  = -0.715263861900;
    fLocationZeta[76] = 0.169395306800;
    fWeight[76]       = 0.015058272472;
    
    fLocationKsi[77]  = 0.832059834869;
    fLocationEta[77]  = -0.832059834869;
    fLocationZeta[77] = 0.033765242900;
    fWeight[77]       = 0.009677225986;
    
    fLocationKsi[78]  = 0.029076476731;
    fLocationEta[78]  = -0.011479542519;
    fLocationZeta[78] = 0.966234757100;
    fWeight[78]       = 0.000022154989;
    
    fLocationKsi[79]  = 0.145872449700;
    fLocationEta[79]  = -0.057591193187;
    fLocationZeta[79] = 0.830604693200;
    fWeight[79]       = 0.001174181822;
    
    fLocationKsi[80]  = 0.327826332945;
    fLocationEta[80]  = -0.129427521860;
    fLocationZeta[80] = 0.619309593000;
    fWeight[80]       = 0.007691697187;
    
    fLocationKsi[81]  = 0.533309978655;
    fLocationEta[81]  = -0.210553521740;
    fLocationZeta[81] = 0.380690407000;
    fWeight[81]       = 0.020356069562;
    
    fLocationKsi[82]  = 0.715263861900;
    fLocationEta[82]  = -0.282389850413;
    fLocationZeta[82] = 0.169395306800;
    fWeight[82]       = 0.028230681769;
    
    fLocationKsi[83]  = 0.832059834869;
    fLocationEta[83]  = -0.328501501081;
    fLocationZeta[83] = 0.033765242900;
    fWeight[83]       = 0.018142498598;
    
    fLocationKsi[84]  = 0.029076476731;
    fLocationEta[84]  = 0.011479542519;
    fLocationZeta[84] = 0.966234757100;
    fWeight[84]       = 0.000022154989;
    
    fLocationKsi[85]  = 0.145872449700;
    fLocationEta[85]  = 0.057591193187;
    fLocationZeta[85] = 0.830604693200;
    fWeight[85]       = 0.001174181822;
    
    fLocationKsi[86]  = 0.327826332945;
    fLocationEta[86]  = 0.129427521860;
    fLocationZeta[86] = 0.619309593000;
    fWeight[86]       = 0.007691697187;
    
    fLocationKsi[87]  = 0.533309978655;
    fLocationEta[87]  = 0.210553521740;
    fLocationZeta[87] = 0.380690407000;
    fWeight[87]       = 0.020356069562;
    
    fLocationKsi[88]  = 0.715263861900;
    fLocationEta[88]  = 0.282389850413;
    fLocationZeta[88] = 0.169395306800;
    fWeight[88]       = 0.028230681769;
    
    fLocationKsi[89]  = 0.832059834869;
    fLocationEta[89]  = 0.328501501081;
    fLocationZeta[89] = 0.033765242900;
    fWeight[89]       = 0.018142498598;
    
    fLocationKsi[90]  = 0.029076476731;
    fLocationEta[90]  = 0.029076476731;
    fLocationZeta[90] = 0.966234757100;
    fWeight[90]       = 0.000011817492;
    
    fLocationKsi[91]  = 0.145872449700;
    fLocationEta[91]  = 0.145872449700;
    fLocationZeta[91] = 0.830604693200;
    fWeight[91]       = 0.000626309699;
    
    fLocationKsi[92]  = 0.327826332945;
    fLocationEta[92]  = 0.327826332945;
    fLocationZeta[92] = 0.619309593000;
    fWeight[92]       = 0.004102758586;
    
    fLocationKsi[93]  = 0.533309978655;
    fLocationEta[93]  = 0.533309978655;
    fLocationZeta[93] = 0.380690407000;
    fWeight[93]       = 0.010857946841;
    
    fLocationKsi[94]  = 0.715263861900;
    fLocationEta[94]  = 0.715263861900;
    fLocationZeta[94] = 0.169395306800;
    fWeight[94]       = 0.015058272472;
    
    fLocationKsi[95]  = 0.832059834869;
    fLocationEta[95]  = 0.832059834869;
    fLocationZeta[95] = 0.033765242900;
    fWeight[95]       = 0.009677225986;
    
    
    break;
    
  case 200:
    fLocationKsi[0]  = -0.017992265904;
    fLocationEta[0]  = -0.017992265904;
    fLocationZeta[0] = 0.980144928200;
    fWeight[0]       = 0.000001120068;
    
    fLocationKsi[1]  = -0.092128370088;
    fLocationEta[1]  = -0.092128370088;
    fLocationZeta[1] = 0.898333238700;
    fWeight[1]       = 0.000064514066;

    fLocationKsi[2]  = -0.214976483841;
    fLocationEta[2]  = -0.214976483841;
    fLocationZeta[2] = 0.762766204950;
    fWeight[2]       = 0.000495536359;

    fLocationKsi[3]  = -0.369977534959;
    fLocationEta[3]  = -0.369977534959;
    fLocationZeta[3] = 0.591717321200;
    fWeight[3]       = 0.001696870666;

    fLocationKsi[4]  = -0.536202310941;
    fLocationEta[4]  = -0.536202310941;
    fLocationZeta[4] = 0.408282678800;
    fWeight[4]       = 0.003564145260;

    fLocationKsi[5]  = -0.691203362059;
    fLocationEta[5]  = -0.691203362059;
    fLocationZeta[5] = 0.237233795050;
    fWeight[5]       = 0.005122775204;

    fLocationKsi[6]  = -0.814051475812;
    fLocationEta[6]  = -0.814051475812;
    fLocationZeta[6] = 0.101666761300;
    fWeight[6]       = 0.005036993351;

    fLocationKsi[7]  = -0.888187579996;
    fLocationEta[7]  = -0.888187579996;
    fLocationZeta[7] = 0.019855071800;
    fWeight[7]       = 0.002729494631;

    fLocationKsi[8]  = -0.017992265904;
    fLocationEta[8]  = -0.010691346814;
    fLocationZeta[8] = 0.980144928200;
    fWeight[8]       = 0.000002262710;

    fLocationKsi[9]  = -0.092128370088;
    fLocationEta[9]  = -0.054744430817;
    fLocationZeta[9] = 0.898333238700;
    fWeight[9]       = 0.000130328316;

    fLocationKsi[10]  = -0.214976483841;
    fLocationEta[10]  = -0.127743117953;
    fLocationZeta[10] = 0.762766204950;
    fWeight[10]       = 0.001001059498;

    fLocationKsi[11]  = -0.369977534959;
    fLocationEta[11]  = -0.219847692379;
    fLocationZeta[11] = 0.591717321200;
    fWeight[11]       = 0.003427939175;

    fLocationKsi[12]  = -0.536202310941;
    fLocationEta[12]  = -0.318621617721;
    fLocationZeta[12] = 0.408282678800;
    fWeight[12]       = 0.007200120438;

    fLocationKsi[13]  = -0.691203362059;
    fLocationEta[13]  = -0.410726192147;
    fLocationZeta[13] = 0.237233795050;
    fWeight[13]       = 0.010348792137;

    fLocationKsi[14]  = -0.814051475812;
    fLocationEta[14]  = -0.483724879283;
    fLocationZeta[14] = 0.101666761300;
    fWeight[14]       = 0.010175499630;

    fLocationKsi[15]  = -0.888187579996;
    fLocationEta[15]  = -0.527777963286;
    fLocationZeta[15] = 0.019855071800;
    fWeight[15]       = 0.005513998067;

    fLocationKsi[16]  = -0.017992265904;
    fLocationEta[16]  = 0.000000000000;
    fLocationZeta[16] = 0.980144928200;
    fWeight[16]       = 0.000002689414;

    fLocationKsi[17]  = -0.092128370088;
    fLocationEta[17]  = 0.000000000000;
    fLocationZeta[17] = 0.898333238700;
    fWeight[17]       = 0.000154905745;

    fLocationKsi[18]  = -0.214976483841;
    fLocationEta[18]  = 0.000000000000;
    fLocationZeta[18] = 0.762766204950;
    fWeight[18]       = 0.001189840185;

    fLocationKsi[19]  = -0.369977534959;
    fLocationEta[19]  = 0.000000000000;
    fLocationZeta[19] = 0.591717321200;
    fWeight[19]       = 0.004074382979;

    fLocationKsi[20]  = -0.536202310941;
    fLocationEta[20]  = 0.000000000000;
    fLocationZeta[20] = 0.408282678800;
    fWeight[20]       = 0.008557925524;

    fLocationKsi[21]  = -0.691203362059;
    fLocationEta[21]  = 0.000000000000;
    fLocationZeta[21] = 0.237233795050;
    fWeight[21]       = 0.012300376520;

    fLocationKsi[22]  = -0.814051475812;
    fLocationEta[22]  = 0.000000000000;
    fLocationZeta[22] = 0.101666761300;
    fWeight[22]       = 0.012094404358;

    fLocationKsi[23]  = -0.888187579996;
    fLocationEta[23]  = 0.000000000000;
    fLocationZeta[23] = 0.019855071800;
    fWeight[23]       = 0.006553832704;

    fLocationKsi[24]  = -0.017992265904;
    fLocationEta[24]  = 0.010691346814;
    fLocationZeta[24] = 0.980144928200;
    fWeight[24]       = 0.000002262710;

    fLocationKsi[25]  = -0.092128370088;
    fLocationEta[25]  = 0.054744430817;
    fLocationZeta[25] = 0.898333238700;
    fWeight[25]       = 0.000130328316;

    fLocationKsi[26]  = -0.214976483841;
    fLocationEta[26]  = 0.127743117953;
    fLocationZeta[26] = 0.762766204950;
    fWeight[26]       = 0.001001059498;

    fLocationKsi[27]  = -0.369977534959;
    fLocationEta[27]  = 0.219847692379;
    fLocationZeta[27] = 0.591717321200;
    fWeight[27]       = 0.003427939175;

    fLocationKsi[28]  = -0.536202310941;
    fLocationEta[28]  = 0.318621617721;
    fLocationZeta[28] = 0.408282678800;
    fWeight[28]       = 0.007200120438;

    fLocationKsi[29]  = -0.691203362059;
    fLocationEta[29]  = 0.410726192147;
    fLocationZeta[29] = 0.237233795050;
    fWeight[29]       = 0.010348792137;

    fLocationKsi[30]  = -0.814051475812;
    fLocationEta[30]  = 0.483724879283;
    fLocationZeta[30] = 0.101666761300;
    fWeight[30]       = 0.010175499630;

    fLocationKsi[31]  = -0.888187579996;
    fLocationEta[31]  = 0.527777963286;
    fLocationZeta[31] = 0.019855071800;
    fWeight[31]       = 0.005513998067;

    fLocationKsi[32]  = -0.017992265904;
    fLocationEta[32]  = 0.017992265904;
    fLocationZeta[32] = 0.980144928200;
    fWeight[32]       = 0.000001120068;

    fLocationKsi[33]  = -0.092128370088;
    fLocationEta[33]  = 0.092128370088;
    fLocationZeta[33] = 0.898333238700;
    fWeight[33]       = 0.000064514066;

    fLocationKsi[34]  = -0.214976483841;
    fLocationEta[34]  = 0.214976483841;
    fLocationZeta[34] = 0.762766204950;
    fWeight[34]       = 0.000495536359;

    fLocationKsi[35]  = -0.369977534959;
    fLocationEta[35]  = 0.369977534959;
    fLocationZeta[35] = 0.591717321200;
    fWeight[35]       = 0.001696870666;

    fLocationKsi[36]  = -0.536202310941;
    fLocationEta[36]  = 0.536202310941;
    fLocationZeta[36] = 0.408282678800;
    fWeight[36]       = 0.003564145260;

    fLocationKsi[37]  = -0.691203362059;
    fLocationEta[37]  = 0.691203362059;
    fLocationZeta[37] = 0.237233795050;
    fWeight[37]       = 0.005122775204;

    fLocationKsi[38]  = -0.814051475812;
    fLocationEta[38]  = 0.814051475812;
    fLocationZeta[38] = 0.101666761300;
    fWeight[38]       = 0.005036993351;

    fLocationKsi[39]  = -0.888187579996;
    fLocationEta[39]  = 0.888187579996;
    fLocationZeta[39] = 0.019855071800;
    fWeight[39]       = 0.002729494631;

    fLocationKsi[40]  = -0.010691346814;
    fLocationEta[40]  = -0.017992265904;
    fLocationZeta[40] = 0.980144928200;
    fWeight[40]       = 0.000002262710;

    fLocationKsi[41]  = -0.054744430817;
    fLocationEta[41]  = -0.092128370088;
    fLocationZeta[41] = 0.898333238700;
    fWeight[41]       = 0.000130328316;

    fLocationKsi[42]  = -0.127743117953;
    fLocationEta[42]  = -0.214976483841;
    fLocationZeta[42] = 0.762766204950;
    fWeight[42]       = 0.001001059498;

    fLocationKsi[43]  = -0.219847692379;
    fLocationEta[43]  = -0.369977534959;
    fLocationZeta[43] = 0.591717321200;
    fWeight[43]       = 0.003427939175;

    fLocationKsi[44]  = -0.318621617721;
    fLocationEta[44]  = -0.536202310941;
    fLocationZeta[44] = 0.408282678800;
    fWeight[44]       = 0.007200120438;

    fLocationKsi[45]  = -0.410726192147;
    fLocationEta[45]  = -0.691203362059;
    fLocationZeta[45] = 0.237233795050;
    fWeight[45]       = 0.010348792137;

    fLocationKsi[46]  = -0.483724879283;
    fLocationEta[46]  = -0.814051475812;
    fLocationZeta[46] = 0.101666761300;
    fWeight[46]       = 0.010175499630;

    fLocationKsi[47]  = -0.527777963286;
    fLocationEta[47]  = -0.888187579996;
    fLocationZeta[47] = 0.019855071800;
    fWeight[47]       = 0.005513998067;

    fLocationKsi[48]  = -0.010691346814;
    fLocationEta[48]  = -0.010691346814;
    fLocationZeta[48] = 0.980144928200;
    fWeight[48]       = 0.000004571022;

    fLocationKsi[49]  = -0.054744430817;
    fLocationEta[49]  = -0.054744430817;
    fLocationZeta[49] = 0.898333238700;
    fWeight[49]       = 0.000263283200;

    fLocationKsi[50]  = -0.127743117953;
    fLocationEta[50]  = -0.127743117953;
    fLocationZeta[50] = 0.762766204950;
    fWeight[50]       = 0.002022293825;

    fLocationKsi[51]  = -0.219847692379;
    fLocationEta[51]  = -0.219847692379;
    fLocationZeta[51] = 0.591717321200;
    fWeight[51]       = 0.006924963241;

    fLocationKsi[52]  = -0.318621617721;
    fLocationEta[52]  = -0.318621617721;
    fLocationZeta[52] = 0.408282678800;
    fWeight[52]       = 0.014545348333;

    fLocationKsi[53]  = -0.410726192147;
    fLocationEta[53]  = -0.410726192147;
    fLocationZeta[53] = 0.237233795050;
    fWeight[53]       = 0.020906148412;

    fLocationKsi[54]  = -0.483724879283;
    fLocationEta[54]  = -0.483724879283;
    fLocationZeta[54] = 0.101666761300;
    fWeight[54]       = 0.020556070951;

    fLocationKsi[55]  = -0.527777963286;
    fLocationEta[55]  = -0.527777963286;
    fLocationZeta[55] = 0.019855071800;
    fWeight[55]       = 0.011139122364;

    fLocationKsi[56]  = -0.010691346814;
    fLocationEta[56]  = 0.000000000000;
    fLocationZeta[56] = 0.980144928200;
    fWeight[56]       = 0.000005433029;

    fLocationKsi[57]  = -0.054744430817;
    fLocationEta[57]  = 0.000000000000;
    fLocationZeta[57] = 0.898333238700;
    fWeight[57]       = 0.000312933379;

    fLocationKsi[58]  = -0.127743117953;
    fLocationEta[58]  = 0.000000000000;
    fLocationZeta[58] = 0.762766204950;
    fWeight[58]       = 0.002403659785;

    fLocationKsi[59]  = -0.219847692379;
    fLocationEta[59]  = 0.000000000000;
    fLocationZeta[59] = 0.591717321200;
    fWeight[59]       = 0.008230878940;

    fLocationKsi[60]  = -0.318621617721;
    fLocationEta[60]  = 0.000000000000;
    fLocationZeta[60] = 0.408282678800;
    fWeight[60]       = 0.017288322999;

    fLocationKsi[61]  = -0.410726192147;
    fLocationEta[61]  = 0.000000000000;
    fLocationZeta[61] = 0.237233795050;
    fWeight[61]       = 0.024848648388;

    fLocationKsi[62]  = -0.483724879283;
    fLocationEta[62]  = 0.000000000000;
    fLocationZeta[62] = 0.101666761300;
    fWeight[62]       = 0.024432553010;

    fLocationKsi[63]  = -0.527777963286;
    fLocationEta[63]  = 0.000000000000;
    fLocationZeta[63] = 0.019855071800;
    fWeight[63]       = 0.013239747921;

    fLocationKsi[64]  = -0.010691346814;
    fLocationEta[64]  = 0.010691346814;
    fLocationZeta[64] = 0.980144928200;
    fWeight[64]       = 0.000004571022;

    fLocationKsi[65]  = -0.054744430817;
    fLocationEta[65]  = 0.054744430817;
    fLocationZeta[65] = 0.898333238700;
    fWeight[65]       = 0.000263283200;

    fLocationKsi[66]  = -0.127743117953;
    fLocationEta[66]  = 0.127743117953;
    fLocationZeta[66] = 0.762766204950;
    fWeight[66]       = 0.002022293825;

    fLocationKsi[67]  = -0.219847692379;
    fLocationEta[67]  = 0.219847692379;
    fLocationZeta[67] = 0.591717321200;
    fWeight[67]       = 0.006924963241;

    fLocationKsi[68]  = -0.318621617721;
    fLocationEta[68]  = 0.318621617721;
    fLocationZeta[68] = 0.408282678800;
    fWeight[68]       = 0.014545348333;

    fLocationKsi[69]  = -0.410726192147;
    fLocationEta[69]  = 0.410726192147;
    fLocationZeta[69] = 0.237233795050;
    fWeight[69]       = 0.020906148412;

    fLocationKsi[70]  = -0.483724879283;
    fLocationEta[70]  = 0.483724879283;
    fLocationZeta[70] = 0.101666761300;
    fWeight[70]       = 0.020556070951;

    fLocationKsi[71]  = -0.527777963286;
    fLocationEta[71]  = 0.527777963286;
    fLocationZeta[71] = 0.019855071800;
    fWeight[71]       = 0.011139122364;

    fLocationKsi[72]  = -0.010691346814;
    fLocationEta[72]  = 0.017992265904;
    fLocationZeta[72] = 0.980144928200;
    fWeight[72]       = 0.000002262710;

    fLocationKsi[73]  = -0.054744430817;
    fLocationEta[73]  = 0.092128370088;
    fLocationZeta[73] = 0.898333238700;
    fWeight[73]       = 0.000130328316;

    fLocationKsi[74]  = -0.127743117953;
    fLocationEta[74]  = 0.214976483841;
    fLocationZeta[74] = 0.762766204950;
    fWeight[74]       = 0.001001059498;

    fLocationKsi[75]  = -0.219847692379;
    fLocationEta[75]  = 0.369977534959;
    fLocationZeta[75] = 0.591717321200;
    fWeight[75]       = 0.003427939175;

    fLocationKsi[76]  = -0.318621617721;
    fLocationEta[76]  = 0.536202310941;
    fLocationZeta[76] = 0.408282678800;
    fWeight[76]       = 0.007200120438;

    fLocationKsi[77]  = -0.410726192147;
    fLocationEta[77]  = 0.691203362059;
    fLocationZeta[77] = 0.237233795050;
    fWeight[77]       = 0.010348792137;

    fLocationKsi[78]  = -0.483724879283;
    fLocationEta[78]  = 0.814051475812;
    fLocationZeta[78] = 0.101666761300;
    fWeight[78]       = 0.010175499630;

    fLocationKsi[79]  = -0.527777963286;
    fLocationEta[79]  = 0.888187579996;
    fLocationZeta[79] = 0.019855071800;
    fWeight[79]       = 0.005513998067;

    fLocationKsi[80]  = 0.000000000000;
    fLocationEta[80]  = -0.017992265904;
    fLocationZeta[80] = 0.980144928200;
    fWeight[80]       = 0.000002689414;

    fLocationKsi[81]  = 0.000000000000;
    fLocationEta[81]  = -0.092128370088;
    fLocationZeta[81] = 0.898333238700;
    fWeight[81]       = 0.000154905745;

    fLocationKsi[82]  = 0.000000000000;
    fLocationEta[82]  = -0.214976483841;
    fLocationZeta[82] = 0.762766204950;
    fWeight[82]       = 0.001189840185;

    fLocationKsi[83]  = 0.000000000000;
    fLocationEta[83]  = -0.369977534959;
    fLocationZeta[83] = 0.591717321200;
    fWeight[83]       = 0.004074382979;

    fLocationKsi[84]  = 0.000000000000;
    fLocationEta[84]  = -0.536202310941;
    fLocationZeta[84] = 0.408282678800;
    fWeight[84]       = 0.008557925524;

    fLocationKsi[85]  = 0.000000000000;
    fLocationEta[85]  = -0.691203362059;
    fLocationZeta[85] = 0.237233795050;
    fWeight[85]       = 0.012300376520;

    fLocationKsi[86]  = 0.000000000000;
    fLocationEta[86]  = -0.814051475812;
    fLocationZeta[86] = 0.101666761300;
    fWeight[86]       = 0.012094404358;

    fLocationKsi[87]  = 0.000000000000;
    fLocationEta[87]  = -0.888187579996;
    fLocationZeta[87] = 0.019855071800;
    fWeight[87]       = 0.006553832704;

    fLocationKsi[88]  = 0.000000000000;
    fLocationEta[88]  = -0.010691346814;
    fLocationZeta[88] = 0.980144928200;
    fWeight[88]       = 0.000005433029;

    fLocationKsi[89]  = 0.000000000000;
    fLocationEta[89]  = -0.054744430817;
    fLocationZeta[89] = 0.898333238700;
    fWeight[89]       = 0.000312933379;

    fLocationKsi[90]  = 0.000000000000;
    fLocationEta[90]  = -0.127743117953;
    fLocationZeta[90] = 0.762766204950;
    fWeight[90]       = 0.002403659785;

    fLocationKsi[91]  = 0.000000000000;
    fLocationEta[91]  = -0.219847692379;
    fLocationZeta[91] = 0.591717321200;
    fWeight[91]       = 0.008230878940;

    fLocationKsi[92]  = 0.000000000000;
    fLocationEta[92]  = -0.318621617721;
    fLocationZeta[92] = 0.408282678800;
    fWeight[92]       = 0.017288322999;

    fLocationKsi[93]  = 0.000000000000;
    fLocationEta[93]  = -0.410726192147;
    fLocationZeta[93] = 0.237233795050;
    fWeight[93]       = 0.024848648388;

    fLocationKsi[94]  = 0.000000000000;
    fLocationEta[94]  = -0.483724879283;
    fLocationZeta[94] = 0.101666761300;
    fWeight[94]       = 0.024432553010;

    fLocationKsi[95]  = 0.000000000000;
    fLocationEta[95]  = -0.527777963286;
    fLocationZeta[95] = 0.019855071800;
    fWeight[95]       = 0.013239747921;

    fLocationKsi[96]  = 0.000000000000;
    fLocationEta[96]  = 0.000000000000;
    fLocationZeta[96] = 0.980144928200;
    fWeight[96]       = 0.000006457595;

    fLocationKsi[97]  = 0.000000000000;
    fLocationEta[97]  = 0.000000000000;
    fLocationZeta[97] = 0.898333238700;
    fWeight[97]       = 0.000371946632;

    fLocationKsi[98]  = 0.000000000000;
    fLocationEta[98]  = 0.000000000000;
    fLocationZeta[98] = 0.762766204950;
    fWeight[98]       = 0.002856944076;

    fLocationKsi[99]  = 0.000000000000;
    fLocationEta[99]  = 0.000000000000;
    fLocationZeta[99] = 0.591717321200;
    fWeight[99]       = 0.009783065378;

    fLocationKsi[100]  = 0.000000000000;
    fLocationEta[100]  = 0.000000000000;
    fLocationZeta[100] = 0.408282678800;
    fWeight[100]       = 0.020548570255;

    fLocationKsi[101]  = 0.000000000000;
    fLocationEta[101]  = 0.000000000000;
    fLocationZeta[101] = 0.237233795050;
    fWeight[101]       = 0.029534628500;

    fLocationKsi[102]  = 0.000000000000;
    fLocationEta[102]  = 0.000000000000;
    fLocationZeta[102] = 0.101666761300;
    fWeight[102]       = 0.029040065487;

    fLocationKsi[103]  = 0.000000000000;
    fLocationEta[103]  = 0.000000000000;
    fLocationZeta[103] = 0.019855071800;
    fWeight[103]       = 0.015736511305;

    fLocationKsi[104]  = 0.000000000000;
    fLocationEta[104]  = 0.010691346814;
    fLocationZeta[104] = 0.980144928200;
    fWeight[104]       = 0.000005433029;

    fLocationKsi[105]  = 0.000000000000;
    fLocationEta[105]  = 0.054744430817;
    fLocationZeta[105] = 0.898333238700;
    fWeight[105]       = 0.000312933379;

    fLocationKsi[106]  = 0.000000000000;
    fLocationEta[106]  = 0.127743117953;
    fLocationZeta[106] = 0.762766204950;
    fWeight[106]       = 0.002403659785;

    fLocationKsi[107]  = 0.000000000000;
    fLocationEta[107]  = 0.219847692379;
    fLocationZeta[107] = 0.591717321200;
    fWeight[107]       = 0.008230878940;

    fLocationKsi[108]  = 0.000000000000;
    fLocationEta[108]  = 0.318621617721;
    fLocationZeta[108] = 0.408282678800;
    fWeight[108]       = 0.017288322999;

    fLocationKsi[109]  = 0.000000000000;
    fLocationEta[109]  = 0.410726192147;
    fLocationZeta[109] = 0.237233795050;
    fWeight[109]       = 0.024848648388;

    fLocationKsi[110]  = 0.000000000000;
    fLocationEta[110]  = 0.483724879283;
    fLocationZeta[110] = 0.101666761300;
    fWeight[110]       = 0.024432553010;

    fLocationKsi[111]  = 0.000000000000;
    fLocationEta[111]  = 0.527777963286;
    fLocationZeta[111] = 0.019855071800;
    fWeight[111]       = 0.013239747921;

    fLocationKsi[112]  = 0.000000000000;
    fLocationEta[112]  = 0.017992265904;
    fLocationZeta[112] = 0.980144928200;
    fWeight[112]       = 0.000002689414;

    fLocationKsi[113]  = 0.000000000000;
    fLocationEta[113]  = 0.092128370088;
    fLocationZeta[113] = 0.898333238700;
    fWeight[113]       = 0.000154905745;

    fLocationKsi[114]  = 0.000000000000;
    fLocationEta[114]  = 0.214976483841;
    fLocationZeta[114] = 0.762766204950;
    fWeight[114]       = 0.001189840185;

    fLocationKsi[115]  = 0.000000000000;
    fLocationEta[115]  = 0.369977534959;
    fLocationZeta[115] = 0.591717321200;
    fWeight[115]       = 0.004074382979;

    fLocationKsi[116]  = 0.000000000000;
    fLocationEta[116]  = 0.536202310941;
    fLocationZeta[116] = 0.408282678800;
    fWeight[116]       = 0.008557925524;

    fLocationKsi[117]  = 0.000000000000;
    fLocationEta[117]  = 0.691203362059;
    fLocationZeta[117] = 0.237233795050;
    fWeight[117]       = 0.012300376520;

    fLocationKsi[118]  = 0.000000000000;
    fLocationEta[118]  = 0.814051475812;
    fLocationZeta[118] = 0.101666761300;
    fWeight[118]       = 0.012094404358;

    fLocationKsi[119]  = 0.000000000000;
    fLocationEta[119]  = 0.888187579996;
    fLocationZeta[119] = 0.019855071800;
    fWeight[119]       = 0.006553832704;

    fLocationKsi[120]  = 0.010691346814;
    fLocationEta[120]  = -0.017992265904;
    fLocationZeta[120] = 0.980144928200;
    fWeight[120]       = 0.000002262710;

    fLocationKsi[121]  = 0.054744430817;
    fLocationEta[121]  = -0.092128370088;
    fLocationZeta[121] = 0.898333238700;
    fWeight[121]       = 0.000130328316;

    fLocationKsi[122]  = 0.127743117953;
    fLocationEta[122]  = -0.214976483841;
    fLocationZeta[122] = 0.762766204950;
    fWeight[122]       = 0.001001059498;

    fLocationKsi[123]  = 0.219847692379;
    fLocationEta[123]  = -0.369977534959;
    fLocationZeta[123] = 0.591717321200;
    fWeight[123]       = 0.003427939175;

    fLocationKsi[124]  = 0.318621617721;
    fLocationEta[124]  = -0.536202310941;
    fLocationZeta[124] = 0.408282678800;
    fWeight[124]       = 0.007200120438;

    fLocationKsi[125]  = 0.410726192147;
    fLocationEta[125]  = -0.691203362059;
    fLocationZeta[125] = 0.237233795050;
    fWeight[125]       = 0.010348792137;

    fLocationKsi[126]  = 0.483724879283;
    fLocationEta[126]  = -0.814051475812;
    fLocationZeta[126] = 0.101666761300;
    fWeight[126]       = 0.010175499630;

    fLocationKsi[127]  = 0.527777963286;
    fLocationEta[127]  = -0.888187579996;
    fLocationZeta[127] = 0.019855071800;
    fWeight[127]       = 0.005513998067;

    fLocationKsi[128]  = 0.010691346814;
    fLocationEta[128]  = -0.010691346814;
    fLocationZeta[128] = 0.980144928200;
    fWeight[128]       = 0.000004571022;

    fLocationKsi[129]  = 0.054744430817;
    fLocationEta[129]  = -0.054744430817;
    fLocationZeta[129] = 0.898333238700;
    fWeight[129]       = 0.000263283200;

    fLocationKsi[130]  = 0.127743117953;
    fLocationEta[130]  = -0.127743117953;
    fLocationZeta[130] = 0.762766204950;
    fWeight[130]       = 0.002022293825;

    fLocationKsi[131]  = 0.219847692379;
    fLocationEta[131]  = -0.219847692379;
    fLocationZeta[131] = 0.591717321200;
    fWeight[131]       = 0.006924963241;

    fLocationKsi[132]  = 0.318621617721;
    fLocationEta[132]  = -0.318621617721;
    fLocationZeta[132] = 0.408282678800;
    fWeight[132]       = 0.014545348333;

    fLocationKsi[133]  = 0.410726192147;
    fLocationEta[133]  = -0.410726192147;
    fLocationZeta[133] = 0.237233795050;
    fWeight[133]       = 0.020906148412;

    fLocationKsi[134]  = 0.483724879283;
    fLocationEta[134]  = -0.483724879283;
    fLocationZeta[134] = 0.101666761300;
    fWeight[134]       = 0.020556070951;

    fLocationKsi[135]  = 0.527777963286;
    fLocationEta[135]  = -0.527777963286;
    fLocationZeta[135] = 0.019855071800;
    fWeight[135]       = 0.011139122364;

    fLocationKsi[136]  = 0.010691346814;
    fLocationEta[136]  = 0.000000000000;
    fLocationZeta[136] = 0.980144928200;
    fWeight[136]       = 0.000005433029;

    fLocationKsi[137]  = 0.054744430817;
    fLocationEta[137]  = 0.000000000000;
    fLocationZeta[137] = 0.898333238700;
    fWeight[137]       = 0.000312933379;

    fLocationKsi[138]  = 0.127743117953;
    fLocationEta[138]  = 0.000000000000;
    fLocationZeta[138] = 0.762766204950;
    fWeight[138]       = 0.002403659785;

    fLocationKsi[139]  = 0.219847692379;
    fLocationEta[139]  = 0.000000000000;
    fLocationZeta[139] = 0.591717321200;
    fWeight[139]       = 0.008230878940;

    fLocationKsi[140]  = 0.318621617721;
    fLocationEta[140]  = 0.000000000000;
    fLocationZeta[140] = 0.408282678800;
    fWeight[140]       = 0.017288322999;

    fLocationKsi[141]  = 0.410726192147;
    fLocationEta[141]  = 0.000000000000;
    fLocationZeta[141] = 0.237233795050;
    fWeight[141]       = 0.024848648388;

    fLocationKsi[142]  = 0.483724879283;
    fLocationEta[142]  = 0.000000000000;
    fLocationZeta[142] = 0.101666761300;
    fWeight[142]       = 0.024432553010;

    fLocationKsi[143]  = 0.527777963286;
    fLocationEta[143]  = 0.000000000000;
    fLocationZeta[143] = 0.019855071800;
    fWeight[143]       = 0.013239747921;

    fLocationKsi[144]  = 0.010691346814;
    fLocationEta[144]  = 0.010691346814;
    fLocationZeta[144] = 0.980144928200;
    fWeight[144]       = 0.000004571022;

    fLocationKsi[145]  = 0.054744430817;
    fLocationEta[145]  = 0.054744430817;
    fLocationZeta[145] = 0.898333238700;
    fWeight[145]       = 0.000263283200;

    fLocationKsi[146]  = 0.127743117953;
    fLocationEta[146]  = 0.127743117953;
    fLocationZeta[146] = 0.762766204950;
    fWeight[146]       = 0.002022293825;

    fLocationKsi[147]  = 0.219847692379;
    fLocationEta[147]  = 0.219847692379;
    fLocationZeta[147] = 0.591717321200;
    fWeight[147]       = 0.006924963241;

    fLocationKsi[148]  = 0.318621617721;
    fLocationEta[148]  = 0.318621617721;
    fLocationZeta[148] = 0.408282678800;
    fWeight[148]       = 0.014545348333;

    fLocationKsi[149]  = 0.410726192147;
    fLocationEta[149]  = 0.410726192147;
    fLocationZeta[149] = 0.237233795050;
    fWeight[149]       = 0.020906148412;

    fLocationKsi[150]  = 0.483724879283;
    fLocationEta[150]  = 0.483724879283;
    fLocationZeta[150] = 0.101666761300;
    fWeight[150]       = 0.020556070951;

    fLocationKsi[151]  = 0.527777963286;
    fLocationEta[151]  = 0.527777963286;
    fLocationZeta[151] = 0.019855071800;
    fWeight[151]       = 0.011139122364;

    fLocationKsi[152]  = 0.010691346814;
    fLocationEta[152]  = 0.017992265904;
    fLocationZeta[152] = 0.980144928200;
    fWeight[152]       = 0.000002262710;

    fLocationKsi[153]  = 0.054744430817;
    fLocationEta[153]  = 0.092128370088;
    fLocationZeta[153] = 0.898333238700;
    fWeight[153]       = 0.000130328316;

    fLocationKsi[154]  = 0.127743117953;
    fLocationEta[154]  = 0.214976483841;
    fLocationZeta[154] = 0.762766204950;
    fWeight[154]       = 0.001001059498;

    fLocationKsi[155]  = 0.219847692379;
    fLocationEta[155]  = 0.369977534959;
    fLocationZeta[155] = 0.591717321200;
    fWeight[155]       = 0.003427939175;

    fLocationKsi[156]  = 0.318621617721;
    fLocationEta[156]  = 0.536202310941;
    fLocationZeta[156] = 0.408282678800;
    fWeight[156]       = 0.007200120438;

    fLocationKsi[157]  = 0.410726192147;
    fLocationEta[157]  = 0.691203362059;
    fLocationZeta[157] = 0.237233795050;
    fWeight[157]       = 0.010348792137;

    fLocationKsi[158]  = 0.483724879283;
    fLocationEta[158]  = 0.814051475812;
    fLocationZeta[158] = 0.101666761300;
    fWeight[158]       = 0.010175499630;

    fLocationKsi[159]  = 0.527777963286;
    fLocationEta[159]  = 0.888187579996;
    fLocationZeta[159] = 0.019855071800;
    fWeight[159]       = 0.005513998067;

    fLocationKsi[160]  = 0.017992265904;
    fLocationEta[160]  = -0.017992265904;
    fLocationZeta[160] = 0.980144928200;
    fWeight[160]       = 0.000001120068;

    fLocationKsi[161]  = 0.092128370088;
    fLocationEta[161]  = -0.092128370088;
    fLocationZeta[161] = 0.898333238700;
    fWeight[161]       = 0.000064514066;

    fLocationKsi[162]  = 0.214976483841;
    fLocationEta[162]  = -0.214976483841;
    fLocationZeta[162] = 0.762766204950;
    fWeight[162]       = 0.000495536359;

    fLocationKsi[163]  = 0.369977534959;
    fLocationEta[163]  = -0.369977534959;
    fLocationZeta[163] = 0.591717321200;
    fWeight[163]       = 0.001696870666;

    fLocationKsi[164]  = 0.536202310941;
    fLocationEta[164]  = -0.536202310941;
    fLocationZeta[164] = 0.408282678800;
    fWeight[164]       = 0.003564145260;

    fLocationKsi[165]  = 0.691203362059;
    fLocationEta[165]  = -0.691203362059;
    fLocationZeta[165] = 0.237233795050;
    fWeight[165]       = 0.005122775204;

    fLocationKsi[166]  = 0.814051475812;
    fLocationEta[166]  = -0.814051475812;
    fLocationZeta[166] = 0.101666761300;
    fWeight[166]       = 0.005036993351;

    fLocationKsi[167]  = 0.888187579996;
    fLocationEta[167]  = -0.888187579996;
    fLocationZeta[167] = 0.019855071800;
    fWeight[167]       = 0.002729494631;

    fLocationKsi[168]  = 0.017992265904;
    fLocationEta[168]  = -0.010691346814;
    fLocationZeta[168] = 0.980144928200;
    fWeight[168]       = 0.000002262710;

    fLocationKsi[169]  = 0.092128370088;
    fLocationEta[169]  = -0.054744430817;
    fLocationZeta[169] = 0.898333238700;
    fWeight[169]       = 0.000130328316;

    fLocationKsi[170]  = 0.214976483841;
    fLocationEta[170]  = -0.127743117953;
    fLocationZeta[170] = 0.762766204950;
    fWeight[170]       = 0.001001059498;

    fLocationKsi[171]  = 0.369977534959;
    fLocationEta[171]  = -0.219847692379;
    fLocationZeta[171] = 0.591717321200;
    fWeight[171]       = 0.003427939175;

    fLocationKsi[172]  = 0.536202310941;
    fLocationEta[172]  = -0.318621617721;
    fLocationZeta[172] = 0.408282678800;
    fWeight[172]       = 0.007200120438;

    fLocationKsi[173]  = 0.691203362059;
    fLocationEta[173]  = -0.410726192147;
    fLocationZeta[173] = 0.237233795050;
    fWeight[173]       = 0.010348792137;

    fLocationKsi[174]  = 0.814051475812;
    fLocationEta[174]  = -0.483724879283;
    fLocationZeta[174] = 0.101666761300;
    fWeight[174]       = 0.010175499630;

    fLocationKsi[175]  = 0.888187579996;
    fLocationEta[175]  = -0.527777963286;
    fLocationZeta[175] = 0.019855071800;
    fWeight[175]       = 0.005513998067;

    fLocationKsi[176]  = 0.017992265904;
    fLocationEta[176]  = 0.000000000000;
    fLocationZeta[176] = 0.980144928200;
    fWeight[176]       = 0.000002689414;

    fLocationKsi[177]  = 0.092128370088;
    fLocationEta[177]  = 0.000000000000;
    fLocationZeta[177] = 0.898333238700;
    fWeight[177]       = 0.000154905745;

    fLocationKsi[178]  = 0.214976483841;
    fLocationEta[178]  = 0.000000000000;
    fLocationZeta[178] = 0.762766204950;
    fWeight[178]       = 0.001189840185;

    fLocationKsi[179]  = 0.369977534959;
    fLocationEta[179]  = 0.000000000000;
    fLocationZeta[179] = 0.591717321200;
    fWeight[179]       = 0.004074382979;

    fLocationKsi[180]  = 0.536202310941;
    fLocationEta[180]  = 0.000000000000;
    fLocationZeta[180] = 0.408282678800;
    fWeight[180]       = 0.008557925524;

    fLocationKsi[181]  = 0.691203362059;
    fLocationEta[181]  = 0.000000000000;
    fLocationZeta[181] = 0.237233795050;
    fWeight[181]       = 0.012300376520;

    fLocationKsi[182]  = 0.814051475812;
    fLocationEta[182]  = 0.000000000000;
    fLocationZeta[182] = 0.101666761300;
    fWeight[182]       = 0.012094404358;

    fLocationKsi[183]  = 0.888187579996;
    fLocationEta[183]  = 0.000000000000;
    fLocationZeta[183] = 0.019855071800;
    fWeight[183]       = 0.006553832704;

    fLocationKsi[184]  = 0.017992265904;
    fLocationEta[184]  = 0.010691346814;
    fLocationZeta[184] = 0.980144928200;
    fWeight[184]       = 0.000002262710;

    fLocationKsi[185]  = 0.092128370088;
    fLocationEta[185]  = 0.054744430817;
    fLocationZeta[185] = 0.898333238700;
    fWeight[185]       = 0.000130328316;

    fLocationKsi[186]  = 0.214976483841;
    fLocationEta[186]  = 0.127743117953;
    fLocationZeta[186] = 0.762766204950;
    fWeight[186]       = 0.001001059498;

    fLocationKsi[187]  = 0.369977534959;
    fLocationEta[187]  = 0.219847692379;
    fLocationZeta[187] = 0.591717321200;
    fWeight[187]       = 0.003427939175;

    fLocationKsi[188]  = 0.536202310941;
    fLocationEta[188]  = 0.318621617721;
    fLocationZeta[188] = 0.408282678800;
    fWeight[188]       = 0.007200120438;

    fLocationKsi[189]  = 0.691203362059;
    fLocationEta[189]  = 0.410726192147;
    fLocationZeta[189] = 0.237233795050;
    fWeight[189]       = 0.010348792137;

    fLocationKsi[190]  = 0.814051475812;
    fLocationEta[190]  = 0.483724879283;
    fLocationZeta[190] = 0.101666761300;
    fWeight[190]       = 0.010175499630;

    fLocationKsi[191]  = 0.888187579996;
    fLocationEta[191]  = 0.527777963286;
    fLocationZeta[191] = 0.019855071800;
    fWeight[191]       = 0.005513998067;

    fLocationKsi[192]  = 0.017992265904;
    fLocationEta[192]  = 0.017992265904;
    fLocationZeta[192] = 0.980144928200;
    fWeight[192]       = 0.000001120068;

    fLocationKsi[193]  = 0.092128370088;
    fLocationEta[193]  = 0.092128370088;
    fLocationZeta[193] = 0.898333238700;
    fWeight[193]       = 0.000064514066;

    fLocationKsi[194]  = 0.214976483841;
    fLocationEta[194]  = 0.214976483841;
    fLocationZeta[194] = 0.762766204950;
    fWeight[194]       = 0.000495536359;

    fLocationKsi[195]  = 0.369977534959;
    fLocationEta[195]  = 0.369977534959;
    fLocationZeta[195] = 0.591717321200;
    fWeight[195]       = 0.001696870666;

    fLocationKsi[196]  = 0.536202310941;
    fLocationEta[196]  = 0.536202310941;
    fLocationZeta[196] = 0.408282678800;
    fWeight[196]       = 0.003564145260;

    fLocationKsi[197]  = 0.691203362059;
    fLocationEta[197]  = 0.691203362059;
    fLocationZeta[197] = 0.237233795050;
    fWeight[197]       = 0.005122775204;

    fLocationKsi[198]  = 0.814051475812;
    fLocationEta[198]  = 0.814051475812;
    fLocationZeta[198] = 0.101666761300;
    fWeight[198]       = 0.005036993351;

    fLocationKsi[199]  = 0.888187579996;
    fLocationEta[199]  = 0.888187579996;
    fLocationZeta[199] = 0.019855071800;
    fWeight[199]       = 0.002729494631;

    break;

  default:
    PZError << "TPZIntRuleP3D creation : invalid number of integration points "
      " specified\n";
    //			PZError.show();

    fNumInt = 0;
  }
}
////////////////////////////////////////////////////////////////////////////////
TPZInt1Point::~TPZInt1Point() {}
TPZInt1Point::TPZInt1Point() {
}

void  TPZInt1Point::SetOrder(TPZVec<int> &ord) {
  //   fOrdKsi = ord[0];
}

int TPZInt1Point::NPoints() {
   return 1;
}

void TPZInt1Point::Point(int ip, TPZVec<REAL> &pos, REAL &w) {
	if(ip!=0) {
   	cout << "TPZInt1Point:: Bad number point " << ip << endl;
   	return;
   }
   w = 1.;
}

void TPZInt1Point::GetOrder(TPZVec<int> &/* ord */) {
}



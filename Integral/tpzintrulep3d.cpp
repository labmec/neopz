/**
 * @file
 * @brief Contains the implementation of the TPZIntRuleP3D methods. Cubature rule for pyramids.
 */

#include <cmath>

#include "tpzintrulep3d.h"
#include "pzerror.h"
#include "pzvec.h"
#include "tpzgaussrule.h"

#include <fstream>

using namespace std;

TPZIntRuleP3D::TPZIntRuleP3D(int order) {
	if(order < 0)
		order = 1;
	if (order > NRULESPYRAMID_ORDER){
		PZError << "TPZGaussRule creation precision = " << order << " not available\n";
		order = NRULESPYRAMID_ORDER;
	}
	ComputingCubatureRuleForPyramid(order);
}

TPZIntRuleP3D::~TPZIntRuleP3D(){
	fLocationKsi.Resize(0);
	fLocationEta.Resize(0);
	fLocationZeta.Resize(0);
	fWeight.Resize(0);
	fNumInt = 0;
}

void TPZIntRuleP3D::Loc(int i, TPZVec<REAL> &Points) const {
	if (i>=0 && i<fNumInt){
		Points[0] = fLocationKsi[i];
		Points[1] = fLocationEta[i];
		Points[2] = fLocationZeta[i];
		return;
	}
	else {
		PZError << "ERROR(TPZIntRuleP3D::loc) Out of bounds!!\n";
	}
}

REAL TPZIntRuleP3D::W(int i) const {
	if(i>=0 && i<fNumInt)
		return fWeight[i];
	else {
		PZError << "ERROR(TPZIntRuleP3D::w) Out of bounds!!\n";
		return 0.0;
	}
}

/** @brief Auxiliar function to compute the linear transformation [-1,1] into [0,1] : T(ksi) = 1/2 * ksi + 1/2*/
long double TransformM1AndP1ToZeroP1(long double ksi) {
	return 0.5L * (ksi + 1.0L);
}
/** @brief Auxiliar function to compute the linear transformation from [-1,1] to [-(1-z),(1-z)] ; Tz(ksi) = (1-z) * ksi */
long double TransformM1AndP1ToM1PZAndP1MZ(long double zeta,long double ksi) {
	return ((1.0L - zeta)*ksi);
}
/** 
 * Computes the points and weights for Pyramid quadrature rule 
 * order is the degree of the polinomial will be integrated exactly with this cubature rule
 */
void TPZIntRuleP3D::ComputingCubatureRuleForPyramid(int order) {
	int i;
	int j;
	int k;
	int l;

	// Cleaning the vectors
	fWeight.Resize(0);
	fLocationKsi.Resize(0);
	fLocationEta.Resize(0);
	fLocationZeta.Resize(0);
	
	long double wi, wj, wk, xi, xj, xk;
	long double volume;

	TPZVec<long double> leg_w;
	TPZVec<long double> leg_x;
	TPZVec<long double> jacobi_w;
	TPZVec<long double> jacobi_x;

	// Determining the number of points needed for order
	int plane_order = (int)(0.51*(order+2));
	int zeta_order = plane_order + 1;
	fNumInt = plane_order * plane_order * zeta_order;

	fWeight.Resize(fNumInt,0.0L);
	fLocationKsi.Resize(fNumInt,0.0L);
	fLocationEta.Resize(fNumInt,0.0L);
	fLocationZeta.Resize(fNumInt,0.0L);

	/** Computing the rule for Legender format to use for find the integration points at plane z = cte
	 * measure of the plane (1-z)^2 */
	leg_w.Resize(plane_order,0.0L);
	leg_x.Resize(plane_order,0.0L);
	TPZGaussRule intrule(2);
	intrule.ComputingGaussLegendreQuadrature(plane_order,leg_x,leg_w);

	/** Computing the rule to jacobi format */
	jacobi_w.Resize(zeta_order,0.0L);
	jacobi_x.Resize(zeta_order,0.0L);
	intrule.ComputingGaussJacobiQuadrature(zeta_order,2.0L,0.0L,jacobi_x,jacobi_w);
	
	volume = 4.0L / 3.0L;

	l = 0;
	for ( k = 0; k < zeta_order; k++ )
	{
		xk = ( jacobi_x[k] + 1.0L ) / 2.0L;
		wk = jacobi_w[k] / 2.0L;
		for ( j = 0; j < plane_order; j++ )
		{
			xj = leg_x[j];
			wj = leg_w[j];
			for ( i = 0; i < plane_order; i++ )
			{
				xi = leg_x[i];
				wi = leg_w[i];
				fWeight[l] = wi * wj * wk / 4.0L / volume;
				fLocationKsi[l] = xi * ( 1.0L - xk );
				fLocationEta[l] = xj * ( 1.0L - xk );
				fLocationZeta[l] = xk;
				l = l + 1;
			}
		}
	}
	/** Redimension of the weights to take a volume of the pyramid */
	for(k=0;k<fNumInt;k++) 
		fWeight[k] *= volume;
}

void TPZIntRuleP3D::Print(std::ostream &out) {
	int i;
	out << "\nQuadrature For Pyramid - " << " \tNumber of points : " << fNumInt  << std::endl;
	for(i=0;i<fNumInt;i++) {
		out << i << " : \t" << fLocationKsi[i] << " \t" << fLocationEta[i] << " \t" << fLocationZeta[i] << " \t" << fWeight[i] << std::endl;
	}
	out << std::endl;
}

// OBSOLETE
void TPZIntRuleP3D::ComputingCubatureRuleForPyramidOldVersion(int order) {
	if(order > 8) order = 8;
	int NRGAUPO[11] = { 1, 1, 12, 12, 36, 36, 96, 96, 96, 200, 200};
	fNumInt = NRGAUPO[order];
	fLocationKsi.Resize(fNumInt,0.0L);
	fLocationEta.Resize(fNumInt,0.0L);
	fLocationZeta.Resize(fNumInt,0.0L);
	fWeight.Resize(fNumInt,0.0L);
	
	switch(fNumInt){
			
		case 1:
		{
			fLocationKsi[0]  = 0.0L;
			fLocationEta[0]  = 0.0L;
			fLocationZeta[0] = 1.L/4.L;
			fWeight[0]       = (4.0L/3.0L);			
		}
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

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

TPZIntRuleP3D::TPZIntRuleP3D(int &order) {
	if(order < 0)
    {
		order = 1;
    }
//	if (order > NRULESPYRAMID_ORDER){
//		PZError << "TPZGaussRule creation precision = " << order << " not available\n";
//		order = NRULESPYRAMID_ORDER;
//	}
	// Computing integration points and weights to cubature rule for pyramid
	fOrder = ComputingCubatureRuleForPyramid(order);
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
int TPZIntRuleP3D::ComputingCubatureRuleForPyramid(int order) {
	int i, j, k, l;

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
	intrule.ComputingGaussLegendreQuadrature(&plane_order,leg_x,leg_w);

	/** Computing the rule to jacobi format */
	jacobi_w.Resize(zeta_order,0.0L);
	jacobi_x.Resize(zeta_order,0.0L);
	intrule.ComputingGaussJacobiQuadrature(&zeta_order,2.0L,0.0L,jacobi_x,jacobi_w);
	
	//volume = pztopology::TPZPyramid::RefElVolume();      
	volume = 4.0L/3.0L;

	l = 0;
	for(k=0;k<zeta_order;k++) {
		xk = (jacobi_x[k] + 1.0L) / 2.0L;
		wk = jacobi_w[k] / 2.0L;
		for(j=0;j<plane_order;j++) {
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
    {
		fWeight[k] *= volume;
    }
    return order;
}

void TPZIntRuleP3D::Print(std::ostream &out) {
	int i;
	out << "\nQuadrature For Pyramid - " << " \tNumber of points : " << fNumInt  << std::endl;
	for(i=0;i<fNumInt;i++) {
		out << i << " : \t" << fLocationKsi[i] << " \t" << fLocationEta[i] << " \t" << fLocationZeta[i] << " \t" << fWeight[i] << std::endl;
	}
	out << std::endl;
}

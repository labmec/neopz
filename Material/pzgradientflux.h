/**
 * \file
 * @brief Contains the TPZGradientFlux class.
 */
//$Id: pzgradientflux.h,v 1.1 2009-08-28 19:43:44 fortiago Exp $

#ifndef TPZGRADIENTFLUX_H
#define TPZGRADIENTFLUX_H

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzreal.h"
#include "pzmaterialdata.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZGradientFlux{
	
public:
	
	TPZGradientFlux();
	
	TPZGradientFlux(const TPZGradientFlux &cp);
	
	~TPZGradientFlux();
	
	/** @brief Computes numerical flux */
	void ComputeFlux(TPZVec<REAL> &solL, TPZVec<REAL> &solR, const TPZVec<REAL> &normal, TPZVec<REAL> & F);
	
	/** @brief Apply limiter */
	void ApplyLimiter(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright);
	
private:
	
	void ApplyVanAlbadaLimiter(REAL &soll, REAL &solr,
							   const TPZVec<REAL>& gradL, const TPZVec<REAL> &gradR, 
							   const TPZVec<REAL> &normal, 
							   const TPZVec<REAL> &dL, const TPZVec<REAL> & dR);
	
	/** @brief It corrects \f$ soll \f$ and \f$ solr \f$ values */
	void ApplyMinModLimiter(REAL &soll, REAL &solr,
							const TPZVec<REAL>& gradL, const TPZVec<REAL> &gradR, 
							const TPZVec<REAL> &normal, 
							const TPZVec<REAL> &dL, const TPZVec<REAL> & dR);
	
	REAL Dot(const TPZVec<REAL> &A, const TPZVec<REAL> &B){
		double result = 0.;
		int n = A.NElements();
		for(int i = 0; i < n; i++) result += A[i]*B[i];
		return result;
	}
	
};

#endif

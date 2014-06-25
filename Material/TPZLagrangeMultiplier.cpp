//
//  TPZLagrangeMultiplier.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#include "TPZLagrangeMultiplier.h"
#include "pzaxestools.h"



/** @brief Unique identifier for serialization purposes */
int TPZLagrangeMultiplier::ClassId() const
{
    return TPZLagrangeMultiplierID;
}

/** @brief Saves the element data to a stream */
void TPZLagrangeMultiplier::Write(TPZStream &buf, int withclassid)
{
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
}

/** @brief Reads the element data from a stream */
void TPZLagrangeMultiplier::Read(TPZStream &buf, void *context)
{
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fNStateVariables);
    
}

/**
 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZLagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
	TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	
	TPZFNMatrix<660> dphiL, dphiR;
	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
	
	
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr;
    
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl; il++) {
		for(jr=0; jr<nrowr; jr++) {
			ek(il,jr+nrowl) += weight * fMultiplier * (phiL(il) * phiR(jr));
		}
	}
	
    //	// 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr; ir++) {
		for(jl=0; jl<nrowl; jl++) {
			ek(ir+nrowl,jl) += weight * fMultiplier * (phiR(ir) * phiL(jl));
		}
	}
    
}

/**
 * @brief It computes a contribution to residual vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZLagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
    
}



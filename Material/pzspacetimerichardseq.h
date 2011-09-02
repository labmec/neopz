/**
 * \file
 * @brief Contains the TPZSpaceTimeRichardsEq class which implements a 1D space-time Richards' equation.
 */
//$Id: pzspacetimerichardseq.h,v 1.2 2008-10-20 11:56:22 longhin Exp $

#ifndef PZSPACETIMERICHARDSEQ_H
#define PZSPACETIMERICHARDSEQ_H

#include "pzmaterial.h"

class TPZCompMesh;

/**
 * @ingroup material
 * @brief Implemenents a 1D space-time Richards' equation
 * @author Tiago Forti <forti@simworx.com.br>
 */
class TPZSpaceTimeRichardsEq : public TPZMaterial
{
	
protected:
	
	/** @brief Soil parameters */
	REAL fAlpha, fN, fThetaS, fThetaR, fKs;
	
	/** @brief Computes Se coeficient which allows the computation of K and Theta coefficients */
	REAL Se(REAL sol);
	
	/** @brief Computes Theta coefficient from Se coefficient */
	REAL Theta(REAL Se);
	
	REAL DKDsol(REAL sol);
	
	REAL DCDsol(REAL sol);
	
public:
	
    TPZSpaceTimeRichardsEq();
    
    TPZSpaceTimeRichardsEq(int id);
	
    TPZSpaceTimeRichardsEq(int matid, REAL Alpha, REAL N, REAL ThetaS, REAL ThetaR, REAL Ks);
	
    ~TPZSpaceTimeRichardsEq();
	
    void Set(REAL Alpha, REAL N, REAL ThetaS, REAL ThetaR, REAL Ks);
	
	/** @brief It returns the integrable dimension of the material */
	virtual int Dimension();
	
	/** @brief It returns the number of state variables associated with the material */
	virtual int NStateVariables();
	
	//   virtual int HasForcingFunction();
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since April 16, 2007
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);
	
	/** @brief Computes coeficient C(sol) based on current solution sol */
	REAL C_Coef(REAL sol);
	
	/** @brief Computes coeficient K(sol) based on current solution sol */  
	REAL K_Coef(REAL sol);
	
	void AnalysisOfParameters(REAL sol0, REAL solL, char* filename);
	
	static int main();
	
	static TPZCompMesh * CreateMesh(REAL L, REAL Time, int p, int ndiv);
	
	static void DirichletT0(TPZVec<REAL> &x, TPZVec<REAL> &f);
	
};

#endif

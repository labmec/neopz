/**
 * @file
 * @brief Contains the TPZSpaceTimeRichardsEq class which implements a 1D space-time Richards' equation.
 */

#ifndef PZSPACETIMERICHARDSEQ_H
#define PZSPACETIMERICHARDSEQ_H

#include "TPZMaterial.h"

class TPZCompMesh;

/**
 * @ingroup material
 * @brief Implemenents a 1D space-time Richards' equation
 * @author Tiago Forti <forti@simworx.com.br>
 */
class TPZSpaceTimeRichardsEq : public TPZMaterial
{
    static STATE TCoeff, LCoeff, deltaDerivada;
protected:
	
	/** @brief Soil parameters */
	STATE fAlpha, fN, fThetaS, fThetaR, fKs;
	
	/** @brief Computes Se coeficient which allows the computation of K and Theta coefficients */
	STATE Se(STATE sol);
	
	/** @brief Computes Theta coefficient from Se coefficient */
	STATE Theta(STATE Se);
	
	STATE DKDsol(STATE sol);
	
	STATE DCDsol(STATE sol);
	
public:
	
    TPZSpaceTimeRichardsEq();
    
    TPZSpaceTimeRichardsEq(int id);
	
    TPZSpaceTimeRichardsEq(int matid, STATE Alpha, STATE N, STATE ThetaS, STATE ThetaR, STATE Ks);
	
    ~TPZSpaceTimeRichardsEq();
	
    void Set(STATE Alpha, STATE N, STATE ThetaS, STATE ThetaR, STATE Ks);
	
	/** @brief It returns the integrable dimension of the material */
	virtual int Dimension() const override;
	
	/** @brief It returns the number of state variables associated with the material */
	virtual int NStateVariables() const override;
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since April 16, 2007
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
	
	/** @brief Computes coeficient C(sol) based on current solution sol */
	STATE C_Coef(STATE sol);
	
	/** @brief Computes coeficient K(sol) based on current solution sol */  
	STATE K_Coef(STATE sol);
	
	void AnalysisOfParameters(STATE sol0, STATE solL, char* filename);
	
	static int main();
	
	static TPZCompMesh * CreateMesh(REAL L, REAL Time, int p, int ndiv);
	
	static void DirichletT0(TPZVec<REAL> &x, TPZVec<STATE> &f);
    public:
virtual int ClassId() const override;

	
};

#endif

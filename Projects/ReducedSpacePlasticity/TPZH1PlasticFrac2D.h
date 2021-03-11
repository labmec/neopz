//
//  TPZH1PlasticFrac2D.h
//  PZ
//
//  Edited from Agnaldo Farias by Nathan Shauer on 04/08/14.
//
//

#ifndef TPZH1PlasticFrac2D_H
#define TPZH1PlasticFrac2D_H

#include <iostream>

#include "TPZMaterial.h"

#include "pzvec.h"
#include "pzfmatrix.h"
#include "TPZMatElastoPlastic2D.h"

#include "pznlfluidstructureData.h"

#include <iostream>

/**
 * @ingroup material
 * @brief Material  to validate the reduce space.
 */
/**
 **@ingroup elasticity equation.
 * \f$  div(T(u)) + fxy = 0 \f$ (Eq. 1)
 *
 *@ingroup pressure equation (1d)
 * \f$ * div(Q) + dw/dt + qL= 0 (Eq. 2)  \f$
 *
 */


template <class T, class TMEM = TPZElastoPlasticMem>
class TPZH1PlasticFrac2D : public TPZMatElastoPlastic2D<T,TMEM>
{
	
protected:
	
	/** @brief Forcing vector */
	TPZVec<REAL>  ff;
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief term that multiplies the Laplacian operator, outflow to the poros medio and right side
	 * @note \f$fXf f$ => vetor de carga
	 */
	REAL fXf;
	
	/** @brief Uses plain stress
	 * @note \f$fPlaneStress = 1\f$ => Plain stress state
	 * @note \f$fPlaneStress != 1\f$ => Plain Strain state
	 */
	int fPlaneStress;
	
	/** @brief If set,the material will calculate contribution of a plastic material
   */
  bool fSetRunPlasticity;
  
	REAL fmatId;
	REAL fE;//young
	REAL fPoiss;//poisson
	REAL fVisc;//viscosity
	

	/** @brief State: one ou one+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	EState gState;	
public:

	
	TPZH1PlasticFrac2D();
	
	TPZH1PlasticFrac2D(int matid, int dim, REAL young, REAL poiss, REAL visc);
	
	virtual ~TPZH1PlasticFrac2D();
                int ClassId() const override;

	
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name() override { return "TPZH1PlasticFrac2D"; }
	
	//virtual int Dimension() const {return 2;}
	
	virtual int NStateVariables();
	
	/** @brief Set plane problem
	 * planestress = 1 => Plain stress state
	 * planestress != 1 => Plain Strain state
	 */
	void SetfPlaneProblem(int planestress)
	{
		fPlaneStress = planestress;
	}
  
  /** @brief if IsPlasticity is true, it will calculate the contribution of a plastic material
   */
  void SetRunPlasticity(bool IsPlasticity = true);
  
	int MatId()
	{
		return fmatId;
	}
	
  void SetLastState(){ gState = ELastState; }
	void SetCurrentState(){ gState = ECurrentState; }
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
	 * @param datavec [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) override;

	virtual void ContributePlastic(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc) override;
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
	//public:
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) override;
	
	virtual void ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
	virtual void ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef,TPZBndCond &bc);
	virtual void ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
	virtual void ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef,TPZBndCond &bc);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
//	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
//	{
//		DebugStop();
//	}
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
//	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
//																		 REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc)
//	{
//		DebugStop();
//	}
	
	/** @name Contribute methods
	 * @{
	 */
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
//	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
//	{
//		DebugStop();
//	}
	
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since October 07, 2011
	 */
//	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc)
//	{
//		DebugStop();
//	}
	
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
	
	virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override;
	
	void ContributePressure(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
        
        void Read(TPZStream& buf, void* context) override {
            DebugStop();
        }

        void Write(TPZStream &buf, int withclassid) const override{
            DebugStop();
        }

};


template <class T, class TMEM>
int TPZH1PlasticFrac2D<T,TMEM>::ClassId() const{
    return Hash("TPZH1PlasticFrac2D") ^ TPZMatElastoPlastic2D<T,TMEM>::ClassId() << 1;
}

#endif

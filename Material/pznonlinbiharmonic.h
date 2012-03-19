/**
 * \file
 * @brief Contains the TPZNonLinBiharmonic class which implements a discontinuous Galerkin formulation for the non-linear bi-harmonic equation.
 */
// -*- c++ -*-
//$Id: pznonlinbiharmonic.h,v 1.8 2009-09-01 19:44:47 phil Exp $

#ifndef TPZNONLINBIHARMONICHPP
#define TPZNONLINBIHARMONICHPP

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief This class implements discontinuous Galerkin formulation for the non-linear bi-harmonic equation.
 * @since Jan 31, 2005
 * @author Igor Mozolevski e Paulo Bosing
 */
class TPZNonLinBiharmonic : public TPZDiscontinuousGalerkin {
	
private:
	REAL  fXf;
	
	public :
	
	static REAL gLambda1, gLambda2, gSigmaA,gSigmaB, gL_alpha, gM_alpha, gL_betta,
	gM_betta, g_teta, Re;
	static int NorP;
	
	/** @brief Inicialisation of biharmonic material */
	TPZNonLinBiharmonic(int nummat, REAL f);
	
	virtual ~TPZNonLinBiharmonic();
	
	/** @brief Returns the number of norm errors. Default is 3: energy, L2,  H1, semi-norm H2 and H2. */
	virtual int NEvalErrors() {return 8;}
	
	void SetMaterial(REAL &xfin){
		fXf = xfin;
	}
	
	int Dimension() { return 2;}
	
	/** @brief Returns one because of scalar problem */
	int NStateVariables(){
		return 1;
	};
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZBiharmonic"; }
	
	/** @brief Implements integral over  element's volume */
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix<REAL> &ek,
                            TPZFMatrix<REAL> &ef);
	/** @brief Implements integral over  element's volume */
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
							TPZFMatrix<REAL> &ef)
	{
		TPZDiscontinuousGalerkin::Contribute(data,weight,ef);
	}
	/** @brief Implements boundary conditions for continuous Galerkin */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ek,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc);
	
	/** @brief Implements boundary conditions for continuous Galerkin */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc)
	{
		TPZDiscontinuousGalerkin::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 0;}
	
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout);
public:

	virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<REAL> &Solout)
	{
		TPZDiscontinuousGalerkin::SolutionDisc(data,dataleft,dataright,var,Solout);
	}

	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux);
	
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values);//Cedric
	
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<REAL> &ek,
									 TPZFMatrix<REAL> &ef);
	
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<REAL> &ek,
									   TPZFMatrix<REAL> &ef,
									   TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<REAL> &ef)
	{
		TPZDiscontinuousGalerkin::ContributeInterface(data,dataleft,dataright,weight,ef);
	}
	
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<REAL> &ef,
									   TPZBndCond &bc)
	{
		TPZDiscontinuousGalerkin::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
	
};

#endif

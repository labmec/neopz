/**
 * \file
 * @brief Contains the TPZCoupledTransportDarcyBC class.
 */

//$Id: pzcoupledtransportdarcyBC.h,v 1.6 2009-09-01 19:44:47 phil Exp $

#ifndef MATCOUPLEDTRANSPDARCYBC
#define MATCOUPLEDTRANSPDARCYBC

#include <iostream>

#include "pzreal.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzmaterialid.h"
#include "pzcoupledtransportdarcy.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZCoupledTransportDarcyBC : public TPZBndCond{
	
#warning THIS CLASS IS NOT THREADSAFE!!!
protected:
	
	TPZBndCond * fMaterials[2];
	
	TPZBndCond * GetNonNullMaterial(){
		if (this->fMaterials[0]) return this->fMaterials[0];
		if (this->fMaterials[1]) return this->fMaterials[1];
		PZError << "Error! - "  << __PRETTY_FUNCTION__ << std::endl;
		exit (-1);
		// the code will never reach this point
		return 0;
	}
	
	void UpdateConvectionDir(TPZFMatrix &dsol);
	void UpdateConvectionDirInterface(TPZFMatrix &dsolL, TPZFMatrix &dsolR, TPZFMatrix &phiL, TPZFMatrix &phiR);
	
	public :
	
	TPZCoupledTransportDarcyBC(TPZCoupledTransportDarcy * material, int id);
	
	~TPZCoupledTransportDarcyBC();
	
	TPZBndCond * GetCurrentMaterial(){
		const int eq = TPZCoupledTransportDarcy::CurrentEquation();
		if (eq == 0 || eq == 1) return this->fMaterials[eq];
		else {
			PZError << "Error! - " << __PRETTY_FUNCTION__ << std::endl;
			exit (-1);
		}
		// the code will never reach this point
		return 0;
	}
	
	virtual int HasForcingFunction() {
		TPZBndCond * bc = this->GetCurrentMaterial();
		if (bc) return bc->HasForcingFunction();
		return 0;
	}
	
	void SetMaterial(int eq, TPZBndCond * mat){
		if (eq == 0 || eq == 1) this->fMaterials[eq] = mat;
		else {
			PZError << "Error! - " << __PRETTY_FUNCTION__ << std::endl;
			exit (-1);
		}
	}
	
	/** @brief Returns the integrable dimension of the material */
	int Dimension() { 
		return this->GetNonNullMaterial()->Dimension();
	}
	
	virtual int NFluxes(){ return this->GetNonNullMaterial()->NFluxes(); }
	
	int NStateVariables() { return this->GetNonNullMaterial()->NStateVariables(); }
	
	/** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
	virtual int NEvalErrors() {return this->GetNonNullMaterial()->NEvalErrors();}
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux){
		flux.Fill(0.);
	}
	
	void Print(std::ostream & out = std::cout) {
		out << " Boundary condition number = " << Id() << "\n";
	}
	
	void Contribute(TPZMaterialData &data,
					REAL weight,
					TPZFMatrix &ek,
					TPZFMatrix &ef);
	
	void Contribute(TPZMaterialData &data,
					REAL weight,
					TPZFMatrix &ef)
	{
		TPZBndCond::Contribute(data,weight,ef);
	}
	
	void ContributeBC(TPZMaterialData &data,
					  REAL weight,
					  TPZFMatrix &ek,
					  TPZFMatrix &ef,
					  TPZBndCond &bc) {  }
	
    void ContributeBC(TPZMaterialData &data,
					  REAL weight,
					  TPZFMatrix &ef,
					  TPZBndCond &bc)
	{
		TPZBndCond::ContributeBC(data,weight,ef,bc);
	}
	
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol, TPZFMatrix &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &uexact,TPZFMatrix &duexact,TPZVec<REAL> &val){
		val.Fill(0.);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix &ef);
	
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix &ek,
									   TPZFMatrix &ef,
									   TPZBndCond &bc) {
		//NOTHING TO BE DONE HERE
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix &ef,
									   TPZBndCond &bc)
	{
		TPZBndCond::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
	
};

#endif

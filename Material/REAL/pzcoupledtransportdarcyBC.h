/**
 * \file
 * @brief Contains the TPZCoupledTransportDarcyBC class.
 */

#ifndef MATCOUPLEDTRANSPDARCYBC
#define MATCOUPLEDTRANSPDARCYBC

#include <iostream>

#include "pzreal.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzcoupledtransportdarcy.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZCoupledTransportDarcyBC : public TPZBndCond{
	
// #warning THIS CLASS IS NOT THREADSAFE!!!
protected:
	
	TPZBndCond * fMaterials[2];
	
	TPZBndCond * GetNonNullMaterial() const {
		if (this->fMaterials[0]) return this->fMaterials[0];
		if (this->fMaterials[1]) return this->fMaterials[1];
		PZError << "Error! - "  << __PRETTY_FUNCTION__ << std::endl;
		exit (-1);
		// the code will never reach this point
		return 0;
	}
	
	void UpdateConvectionDir(TPZFMatrix<STATE> &dsol);
	void UpdateConvectionDirInterface(TPZFMatrix<STATE> &dsolL, TPZFMatrix<STATE> &dsolR, TPZFMatrix<REAL> &phiL, TPZFMatrix<REAL> &phiR);
	
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
	
	virtual int HasForcingFunction()  override {
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
	int Dimension() const  override {
		return this->GetNonNullMaterial()->Dimension();
	}
	
	virtual int NStateVariables() const override { return this->GetNonNullMaterial()->NStateVariables(); }
	
	/** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
	virtual int NEvalErrors()  override {return this->GetNonNullMaterial()->NEvalErrors();}
	
	void Print(std::ostream & out = std::cout)  override {
		out << " Boundary condition number = " << Id() << "\n";
	}
	
	void Contribute(TPZMaterialData &data,
					REAL weight,
					TPZFMatrix<STATE> &ek,
					TPZFMatrix<STATE> &ef) override;
	
	void Contribute(TPZMaterialData &data,
					REAL weight,
					TPZFMatrix<STATE> &ef) override
	{
		TPZBndCond::Contribute(data,weight,ef);
	}
	
	void ContributeBC(TPZMaterialData &data,
					  REAL weight,
					  TPZFMatrix<STATE> &ek,
					  TPZFMatrix<STATE> &ef,
					  TPZBndCond &bc)  override {  }
	
    void ContributeBC(TPZMaterialData &data,
					  REAL weight,
					  TPZFMatrix<STATE> &ef,
					  TPZBndCond &bc) override
	{
		TPZBndCond::ContributeBC(data,weight,ef,bc);
	}
	
	
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val) override {
		val.Fill(0.);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix<STATE> &ef) override;
	
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ek,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override {
		//NOTHING TO BE DONE HERE
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override
	{
		TPZBndCond::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
public:
virtual int ClassId() const override;

};

#endif

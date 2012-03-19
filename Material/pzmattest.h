/**
 * \file
 * @brief Contains the TPZMaterialTest class. Test.
 */
#ifndef MATTESTHPP
#define MATTESTHPP


#include "pzmaterial.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZMaterialTest : public TPZMaterial {
	
	REAL fNumMat;//material id
	REAL fAlfa;//intensity of singularity
	REAL fX0;//singularity
	TPZFMatrix<REAL> fXf;//fonte
	
	public :
	
	TPZMaterialTest(int nummat, REAL alfa, REAL x0);
	
	virtual ~TPZMaterialTest();
	
	REAL Alfa() {return fAlfa;}
	
	REAL   X0() {return fX0;}
	
	void SetMaterial(TPZFMatrix<REAL> &xkin){
		fXf = xkin;
	}
	
	int Dimension() { return 2;}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMaterialTest"; }
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<REAL> &ek,
							TPZFMatrix<REAL> &ef);
	
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ek,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc);
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<REAL> &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 2;}
	
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	/** @brief Compute the value of the flux function to be used by ZZ error estimator */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
		        TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values);//Cedric
};

#endif

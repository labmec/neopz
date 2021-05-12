/**
 * @file
 * @brief Contains the TPZMaterialTest class. Test.
 */

#ifndef MATTESTHPP
#define MATTESTHPP


#include "TPZMaterial.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZMaterialTest : public TPZMaterial {
	
	int fNumMat;//material id
	STATE fAlfa;//intensity of singularity
	STATE fX0;//singularity
	TPZFMatrix<STATE> fXf;//fonte
	
	public :
	
	TPZMaterialTest(int nummat, STATE alfa, STATE x0);
	
	virtual ~TPZMaterialTest();
	
	STATE Alfa() {return fAlfa;}
	
	STATE   X0() {return fX0;}
	
	void SetMaterial(TPZFMatrix<STATE> &xkin){
		fXf = xkin;
	}
	
//	int Dimension() { return 2;}
	
	virtual int NStateVariables() const override;
	
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name() override { return "TPZMaterialTest"; }
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef) override;
	
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override;
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override
	{
		TPZMaterial::Solution(data,var,Solout);
	}
protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
		        TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
};

#endif

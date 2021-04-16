/**
 * \file
 * @brief Contains the TPZMatHybrid class.
 */
#ifndef HIBRIDAH
#define HIBRIDAH


#include "TPZMaterial.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZMatHybrid : public TPZMaterial {
	
	int fNumMat; /**< material id */
	TPZFMatrix<STATE> fXf; /**< source */
	
	public :
	
	TPZMatHybrid(int nummat);
	
	virtual ~TPZMatHybrid();
	
	void SetMaterial(TPZFMatrix<STATE> &xkin){
		fXf = xkin;
	}
	
	virtual int NStateVariables() const override;
	
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name() override { return "TPZMatHybrid"; }
	
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> & ek,
							TPZFMatrix<STATE> & ef) override;
	
    virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> & ek,
							  TPZFMatrix<STATE> & ef,
							  TPZBndCond & bc) override;
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> & ef) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> & ef,
							  TPZBndCond & bc) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<STATE> &Solout) override
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
};

#endif

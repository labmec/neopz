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
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatHybrid"; }
	
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> & ek,
							TPZFMatrix<STATE> & ef);
	
    virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> & ek,
							  TPZFMatrix<STATE> & ef,
							  TPZBndCond & bc);
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> & ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> & ef,
							  TPZBndCond & bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 3;}
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<STATE> &Solout)
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
};

#endif

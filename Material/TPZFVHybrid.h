/**
 * \file
 * @brief Contains the TPZMatHybrid class.
 */
#ifndef HIBRIDAH
#define HIBRIDAH


#include "pzmaterial.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZMatHybrid : public TPZMaterial {
	
	REAL fNumMat; /**< material id */
	TPZFMatrix fXf; /**< source */
	
	public :
	
	TPZMatHybrid(int nummat);
	
	virtual ~TPZMatHybrid();
	
	void SetMaterial(TPZFMatrix &xkin){
		fXf = xkin;
	}
	
	int Dimension() { return fXf.Rows();}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatHybrid"; }
	
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix & ek,
							TPZFMatrix & ef);
	
    virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix & ek,
							  TPZFMatrix & ef,
							  TPZBndCond & bc);
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix & ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix & ef,
							  TPZBndCond & bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 3;}
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<REAL> &Solout)
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
};

#endif

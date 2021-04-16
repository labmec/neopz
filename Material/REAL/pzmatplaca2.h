/**
 * @file
 * @brief Contains the TPZMatPlaca2 class.
 */

#ifndef TPZMATPLACA2HPP
#define TPZMATPLACA2HPP

#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"

class TPZBndCond;


/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZMatPlaca2 : public TPZMaterial{
	
protected:
	int fIdfMax;
	STATE fE1, fE2, fG12, fG13, fG23, fh,ff,fmi,fni1,fni2;
	TPZFMatrix<STATE> fnaxes;
	TPZFMatrix<STATE> fRmat,fRmatT;
	TPZFMatrix<STATE> fKxxR,fKyyR,fKxyR,fKyxR,fBx0R,fB0xR,fBy0R,fB0yR,fB00R;
	TPZFMatrix<STATE> fKxx,fKyy,fKxy,fKyx,fBx0,fB0x,fBy0,fB0y,fB00;
	TPZVec<STATE> fXF;
	
	/** @brief Modify the direction of the fibres for the plate */
	void SetNAxes(TPZFMatrix<STATE> &n);
	
	public :
	
	TPZMatPlaca2(int num, STATE h, STATE f, STATE E1 , STATE E2 ,
				 STATE ni1 , STATE ni2 , STATE G12 , STATE G13 ,
				 STATE G23 , TPZFMatrix<STATE> &naxes, TPZVec<STATE> &xf);
	
	virtual int NStateVariables() const override { return fIdfMax; }
	
	virtual int Dimension() const  override { return 2; }
	
	void Print(std::ostream & out) override;
	
	virtual std::string Name() override { return "TPZMatPlaca2"; }
	
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

protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
						TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
public:
	
	/** @brief Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name) override;
	
	/**
	 * @brief Returns the number of variables associated with the variable indexed by var.
	 * @param var is obtained by calling VariableIndex
	 */
	virtual int NSolutionVariables(int var) override;
	
protected:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	/** @brief Returns the solution associated with the var index based on
	 * the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override
	{
		TPZMaterial::Solution(data,var,Solout);
	}
    public:
virtual int ClassId() const override;
 
};

#endif

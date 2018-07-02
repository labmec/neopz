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
	
	virtual int NStateVariables() { return fIdfMax; }
	
	virtual int Dimension() const { return 2; }
	
	void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatPlaca2"; }
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCond &bc);
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int NFluxes();
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &fl);
	
	virtual void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
						TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
	
	/** @brief Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	/**
	 * @brief Returns the number of variables associated with the variable indexed by var.
	 * @param var is obtained by calling VariableIndex
	 */
	virtual int NSolutionVariables(int var);
	
protected:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on
	 * the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
    public:
virtual int ClassId() const;
 
};

#endif

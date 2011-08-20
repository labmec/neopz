/**
 * \file
 * @brief Contains the TPZMatPlaca2 class.
 */
#ifndef TPZMATPLACA2HPP
#define TPZMATPLACA2HPP

#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
//#include "pzreal.h"
class TPZBndCond;

//const Float BIGNUMBER = 1.e9;

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZMatPlaca2 : public TPZMaterial{
	
protected:
	int fIdfMax;
	REAL fE1, fE2, fG12, fG13, fG23, fh,ff,fmi,fni1,fni2;
	TPZFMatrix fnaxes;
	TPZFMatrix fRmat,fRmatT;
	TPZFMatrix fKxxR,fKyyR,fKxyR,fKyxR,fBx0R,fB0xR,fBy0R,fB0yR,fB00R;
	TPZFMatrix fKxx,fKyy,fKxy,fKyx,fBx0,fB0x,fBy0,fB0y,fB00;
	TPZVec<REAL> fXF;
	
	/** @brief Modify the direction of the fibres for the plate */
	void SetNAxes(TPZFMatrix &n);
	
	public :
	
	TPZMatPlaca2(int num, REAL h, REAL f, REAL E1 , REAL E2 ,
				 REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
				 REAL G23 , TPZFMatrix &naxes, TPZVec<REAL> &xf);
	
	virtual int NStateVariables() { return fIdfMax; }
	
	//  int NFluxes() { return NStateVariables(); }
	
	int Dimension() { return 2; }
	
	void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatPlaca2"; }
	
	//  virtual TPZBndCond *CreateBC(int num,int typ,TPZFMatrix &val1,TPZFMatrix &val2);
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef,
                              TPZBndCond &bc);
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int NFluxes();
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &fl);
	
	virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
						TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
	
	/** @brief Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	/**
	 * @brief Returns the number of variables associated with the variable indexed by var.
	 * @param var is obtained by calling VariableIndex
	 */
	virtual int NSolutionVariables(int var);
	
protected:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on
	 * the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	
};

#endif


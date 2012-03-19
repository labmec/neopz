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
	TPZFMatrix<REAL> fnaxes;
	TPZFMatrix<REAL> fRmat,fRmatT;
	TPZFMatrix<REAL> fKxxR,fKyyR,fKxyR,fKyxR,fBx0R,fB0xR,fBy0R,fB0yR,fB00R;
	TPZFMatrix<REAL> fKxx,fKyy,fKxy,fKyx,fBx0,fB0x,fBy0,fB0y,fB00;
	TPZVec<REAL> fXF;
	
	/** @brief Modify the direction of the fibres for the plate */
	void SetNAxes(TPZFMatrix<REAL> &n);
	
	public :
	
	TPZMatPlaca2(int num, REAL h, REAL f, REAL E1 , REAL E2 ,
				 REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
				 REAL G23 , TPZFMatrix<REAL> &naxes, TPZVec<REAL> &xf);
	
	virtual int NStateVariables() { return fIdfMax; }
	
	//  int NFluxes() { return NStateVariables(); }
	
	int Dimension() { return 2; }
	
	void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatPlaca2"; }
	
	//  virtual TPZBndCond *CreateBC(int num,int typ,TPZFMatrix<REAL> &val1,TPZFMatrix<REAL> &val2);
	
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
	
	virtual int NFluxes();
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &axes, TPZVec<REAL> &fl);
	
	virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
						TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values);
	
	/** @brief Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	/**
	 * @brief Returns the number of variables associated with the variable indexed by var.
	 * @param var is obtained by calling VariableIndex
	 */
	virtual int NSolutionVariables(int var);
	
protected:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on
	 * the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	
};

#endif


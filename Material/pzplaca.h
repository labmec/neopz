/**
 * \file
 * @brief Contains the TPZPlaca class.
 */
#ifndef PLACAHPP
#define PLACAHPP

#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
//#include "pzreal.h"
class TPZBndCond;
template<class T>
class TPZVec;

//const Float BIGNUMBER = 1.e9;

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZPlaca : public TPZMaterial{
	
	TPZFMatrix fnaxes;
	REAL fE1,fE2,fG12,fG13,fG23,fh,ff,fmi,fni1,fni2;
	TPZFMatrix fRmat, fRmatT;
	TPZFMatrix fKxxR,fKyyR, fKxyR, fKyxR, fBx0R, fB0xR,
	fBy0R,fB0yR, fB00R;
	TPZVec<REAL> fXF;
	public :
	
	TPZPlaca(int num, REAL h, REAL f, REAL E1 , REAL E2 ,
			 REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
			 REAL G23 , TPZFMatrix &naxes, TPZVec<REAL> &xf);
	
	virtual int NStateVariables() { return 6; }
	
	//  int NFluxes() { return NStateVariables(); }
	
	int Dimension() { return 2; }
	
	void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZPlaca"; }
	
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
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
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
	
	/** @brief Exact solution for tests */
	void SetExactFunction( void (*fp)(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &uexact,TPZFMatrix &duexact) )
    {
		fExactFunction = fp;
    }
	
protected:
	
	void (*fExactFunction)(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &uexact,TPZFMatrix &duexact);
	
};

#endif

/**
 * @file
 * @brief Contains the TPZPlaca class.
 */

#ifndef PLACAHPP
#define PLACAHPP

#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"

class TPZBndCond;
template<class T>
class TPZVec;

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZPlaca : public TPZMaterial{
	
	TPZFMatrix<STATE> fnaxes;
	STATE fE1,fE2,fG12,fG13,fG23,fh,ff,fmi,fni1,fni2;
	TPZFMatrix<STATE> fRmat, fRmatT;
	TPZFMatrix<STATE> fKxxR,fKyyR, fKxyR, fKyxR, fBx0R, fB0xR,
	fBy0R,fB0yR, fB00R;
	TPZVec<STATE> fXF;
	public :
	
	TPZPlaca(int num, STATE h, STATE f, STATE E1 , STATE E2 ,
			 STATE ni1 , STATE ni2 , STATE G12 , STATE G13 ,
			 STATE G23 , TPZFMatrix<STATE> &naxes, TPZVec<STATE> &xf);
	
	virtual int NStateVariables() { return 6; }
	
	int Dimension() const { return 2; }
	
	void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZPlaca"; }
	
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
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
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
	
	/** @brief Exact solution for tests */
	void SetExactFunction( void (*fp)(TPZFMatrix<REAL> &axes,TPZVec<REAL> &x,TPZFMatrix<STATE> &uexact,TPZFMatrix<STATE> &duexact) )
    {
		fExactFunction = fp;
    }
	
protected:
	
	void (*fExactFunction)(TPZFMatrix<REAL> &axes,TPZVec<REAL> &x,TPZFMatrix<STATE> &uexact,TPZFMatrix<STATE> &duexact);
    public:
virtual int ClassId() const;

};

#endif

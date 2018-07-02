/**
 * @file
 * @brief Contains the TPZMatOrthotropic class.
 */

#ifndef ORTHOTROPICHPP
#define ORTHOTROPICHPP

#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"

/**
 * @ingroup material
 * @brief Implements a orthotropic material.
 */
class TPZMatOrthotropic : public TPZMaterial {
	
	
	TPZFMatrix<STATE> fKXX,fKYY,fKZZ,fKXY,fKYX,fKXZ,fKZX,fKYZ,fKZY;
 	TPZFMatrix<STATE> fLocAxs;
	STATE fEppx,fEppy,fEppz,fVxy,fVyx,fVyz,fVzy,fVzx,fVxz;
	STATE fNumNom,fGxy,fGzx,fGyz;
	TPZFMatrix<STATE> fXf;
	
	public :
	
	TPZMatOrthotropic(int nummat,TPZFMatrix<STATE> naxes,STATE eppx,STATE eppy,
					  STATE eppz,STATE vxy,STATE vyz,STATE vzx,
		              STATE gxy,STATE gyz,STATE gzx);
	
	virtual ~TPZMatOrthotropic();
	
	void SetMaterial(TPZFMatrix<STATE> &xkin){
		fXf = xkin;
	}
	
//	int Dimension() { return 3;}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatOrthotropic"; }
	
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
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
		        TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);//Cedric
	
	/** @brief Verifies the consistency of the axles */
	/** In case that they are not normal or they are not linearly independent calculates new axles based on the input data */
	void Normalize(TPZFMatrix<STATE> &naxes);
};

#endif


/**
 * \file
 * @brief Contains the TPZMatOrthotropic class.
 */
#ifndef ORTHOTROPICHPP
#define ORTHOTROPICHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
//#include "checkconv.h"
#include "pzvec.h"

/**
 * @ingroup material
 * @brief Implements a orthotropic material.
 */
class TPZMatOrthotropic : public TPZMaterial {
	
	
	TPZFMatrix<REAL> fKXX,fKYY,fKZZ,fKXY,fKYX,fKXZ,fKZX,fKYZ,fKZY;
 	TPZFMatrix<REAL> fLocAxs;
	REAL fEppx,fEppy,fEppz,fVxy,fVyx,fVyz,fVzy,fVzx,fVxz;
	REAL fNumNom,fGxy,fGzx,fGyz;
	TPZFMatrix<REAL> fXf;//fonte
	
	public :
	
	TPZMatOrthotropic(int nummat,TPZFMatrix<REAL> naxes,REAL eppx,REAL eppy,
					  REAL eppz,REAL vxy,REAL vyz,REAL vzx,
		              REAL gxy,REAL gyz,REAL gzx);
	
	virtual ~TPZMatOrthotropic();
	
	void SetMaterial(TPZFMatrix<REAL> &xkin){
		fXf = xkin;
	}
	
	int Dimension() { return 3;}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatOrthotropic"; }
	
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
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 3;}
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout);
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<REAL> &Solout)
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }
		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
		        TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values);//Cedric
	
	/** @brief Verifies the consistency of the axles */
	/** In case that they are not normal or they are not linearly independent calculates new axles based on the input data */
	void Normalize(TPZFMatrix<REAL> &naxes);
};

#endif


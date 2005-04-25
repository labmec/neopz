#ifndef ORTHOTROPICHPP
#define ORTHOTROPICHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
//#include "checkconv.h"
#include "pzvec.h"

class TPZMatOrthotropic : public TPZMaterial {


        TPZFMatrix fKXX,fKYY,fKZZ,fKXY,fKYX,fKXZ,fKZX,fKYZ,fKZY;
 	TPZFMatrix fLocAxs;
	REAL fEppx,fEppy,fEppz,fVxy,fVyx,fVyz,fVzy,fVzx,fVxz;
	REAL fNumNom,fGxy,fGzx,fGyz;
        TPZFMatrix fXf;//fonte         

public :

TPZMatOrthotropic(int nummat,TPZFMatrix naxes,REAL eppx,REAL eppy,
                              REAL eppz,REAL vxy,REAL vyz,REAL vzx,
		              REAL gxy,REAL gyz,REAL gzx);

virtual ~TPZMatOrthotropic();

void SetMaterial(TPZFMatrix &xkin){
   fXf = xkin;
}

int Dimension() { return 3;}

int NStateVariables();

virtual void Print(std::ostream & out);

char *Name() { return "TPZMatOrthotropic"; }

virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);


virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

virtual int VariableIndex(char *name);

virtual int NSolutionVariables(int var);

virtual int NFluxes(){ return 3;}

virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

/**compute the value of the flux function to be used by ZZ error estimator*/
virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);

void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				  TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
		        TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric

/**
 * It verifies the consistency of the axles.
 * In case that they are not normal or they are not linearly
 * independent calculates new axles based on the input data
 */
void TPZMatOrthotropic::Normalize(TPZFMatrix &naxes);
};

#endif


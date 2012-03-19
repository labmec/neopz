/**
 * \file
 * @brief Contains the TPZMultPlaca class.
 */
#ifndef MULTPLACAHPP
#define MULTPLACAHPP

#include "pzmatplaca2.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZMultPlaca :public TPZMatPlaca2 {
	
public:
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout);
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<REAL> &Solout)
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }
		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
	
	TPZMultPlaca(int num, REAL h, TPZVec<REAL> &esp, REAL f, REAL E1 , REAL E2 ,
				 REAL ni1 , REAL ni2 , REAL G12 , REAL G13 , REAL G23 ,
				 TPZFMatrix<REAL> &naxes, TPZVec<REAL> &xf,
				 int camadaref, int camadaatual);
	
private:
	TPZFMatrix<REAL> fT;
};


#endif

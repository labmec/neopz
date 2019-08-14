/**
 * @file
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
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<STATE> &Solout) override
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }
		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
	
	TPZMultPlaca(int num, STATE h, TPZVec<STATE> &esp, STATE f, STATE E1 , STATE E2 ,
				 STATE ni1 , STATE ni2 , STATE G12 , STATE G13 , STATE G23 ,
				 TPZFMatrix<STATE> &naxes, TPZVec<STATE> &xf,
				 int camadaref, int camadaatual);
    public:
int ClassId() const override;

private:
	TPZFMatrix<STATE> fT;
};


#endif

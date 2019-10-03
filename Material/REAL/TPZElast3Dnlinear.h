//
//  TPZElast3Dnlinear.h
//  PZ
//
//  Created by Cesar Lucci (Caju) on 22/10/13.
//  *** Esta classe implementa o material elastico linear 3D no formato de matriz tangente e residuo.
//      A motivacao foi o acoplamento com formulacao nao linar de fluxo de fluido no interior da fratura,
//      em que o sistema global eh resolvido pelo metodo de Newton-Raphson.
//

#ifndef __PZ__TPZElast3Dnlinear__
#define __PZ__TPZElast3Dnlinear__

#include <iostream>

#include "pzelast3d.h"

class TPZElast3Dnlinear : public TPZElasticity3D
{
public:
    TPZElast3Dnlinear();
    
    /**
	 * @brief Class constructor
	 * @param nummat - material ID.
	 * @param E - Young's modulus.
	 * @param poisson - poisson's ratio
	 * @param force - external forces
	 */
	TPZElast3Dnlinear(int nummat, STATE E, STATE poisson, TPZVec<STATE> &force,
                      STATE preStressXX = 0., STATE preStressYY = 0., STATE preStressZZ = 0.);
    
    ~TPZElast3Dnlinear();

    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef) override ;
    
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override;
    
    public:
int ClassId() const override;

    
    virtual void FillDataRequirements(TPZMaterialData &data) override;
    
    /** @brief Creates a new material from the current object */
	virtual TPZMaterial * NewMaterial()  override { return new TPZElast3Dnlinear(*this);}
    
protected:
    
    void ContributeVecShapeAux(TPZMaterialData &data,
                               REAL weight,
                               TPZFMatrix<STATE> &ek,
                               TPZFMatrix<STATE> &ef);
    
    void ContributeVecShapeBCAux(TPZMaterialData &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,
                                 TPZFMatrix<STATE> &ef,
                                 TPZBndCond &bc);
};

#endif /* defined(__PZ__TPZElast3Dnlinear__) */

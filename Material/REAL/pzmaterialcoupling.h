/*
 * @file
 * @brief Contains the declaration of the TPZMaterialCoupling class.
 */

#ifndef PZMATERIALCOUPLING_H
#define PZMATERIALCOUPLING_H

#include "pzpoisson3d.h"
#include "pzmaterialdata.h"
#include "pzfmatrix.h"

/**
 * @brief Implemented a Poisson Problem coupling the interpolation spaces H(div) and H1
 * @author Denise de Siqueira
 * @since 8/1/11.
 */
class TPZMaterialCoupling : public TPZMatPoisson3d {
	
public:
	/** @brief Default constructor */
	TPZMaterialCoupling();
	/** @brief Simple constructor */
	TPZMaterialCoupling(int nummat, int dim);

	/** @brief Method to possibilite the coupling between H(div) and H1 */
	virtual void ContributeInterface(TPZMaterialData &data,TPZMaterialData &dataleft,TPZMaterialData &dataright, 
									 REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override ;
	
	virtual void ContributeInterface2(TPZMaterialData &data, TPZMaterialData &dataleft,TPZMaterialData &dataright, REAL weight,TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
    virtual void InitMaterialData(TPZMaterialData &data);		
	
	/** @brief Destructor */
    virtual ~TPZMaterialCoupling();
    public:
int ClassId() const override;
};

#endif

/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCompElHDivConstant_H
#define TPZCompElHDivConstant_H

#include "pzelctemp.h"
#include "TPZBndCond.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"
#include "pzelchdiv.h"

/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivConstant : public TPZCompElHDiv<TSHAPE>  {

protected:
  ///! Indexes of the connects associated with the elements
  TPZManVector<int64_t,TSHAPE::NSides> fConnectIndexes =
    TPZManVector<int64_t,TSHAPE::NSides>(TSHAPE::NSides,-1);
private:

    int fSideOrient = 1;

public:
	    
	TPZCompElHDivConstant();
    
    TPZCompElHDivConstant(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
    virtual void InitMaterialData(TPZMaterialData &data) override;

    // void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override;
    void ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data) override;
    //  void ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data) override;
    void SetSideOrient(int orient);

    int GetSideOrient();

	/** @brief Compute the solution for a given variable */
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;

    void AdjustConnects();
    /**
   * @brief Number of shapefunctions of the connect associated
   * @param connect connect number
   * @return number of shape functions
   */
    int NConnectShapeF(int connect, int order) const override;

protected:
	//  //@{
    // /** @brief Compute the solution using Hdiv structure */
	// void ReallyComputeSolution(TPZMaterialDataT<STATE> &data) override{
    //     ComputeSolutionHDivConstantT(data);
    // }
    // void ReallyComputeSolution(TPZMaterialDataT<CSTATE> &data) override{
    //     ComputeSolutionHDivConstantT(data);
    // }

	// template<class TVar>
    // void ComputeSolutionHDivConstantT(TPZMaterialDataT<TVar> &data);
    // template<class TVar>
    // void ComputeRequiredDataT(TPZMaterialDataT<TVar> &data, TPZVec<REAL>&qsi);
	
};


#endif

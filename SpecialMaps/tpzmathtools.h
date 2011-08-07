/**
 * @file
 * @brief Contains the TPZMathTools class.
 */
#ifndef TPZMATHTOOLS_H
#define TPZMATHTOOLS_H

#include "pzvec.h"
#include "pzgeoel.h"

/**
 * @author Paulo Cesar de Alvarenga Lucci (Caju)
 * @since 2007
 * @brief ??
 */
class TPZMathTools
{
	
public:
	
    TPZMathTools();
	~TPZMathTools();
	
	void JacobianConv(TPZGeoEl &Object, TPZVec< REAL > StartPoint);
	void JacobianConv(TPZGeoElSide &Object, TPZVec< REAL > StartPoint);
	
	/**
     * @brief Type in Function method (.cpp arquive) the function to be integrated
	 */
	static void Function(const TPZVec<REAL> &x, REAL &fx);
	
	REAL IntegrateFunction(void (func)(const TPZVec<REAL> &coord, REAL &result), TPZGeoEl *Geo);
	//      REAL IntegrateFunction(TPZGeoEl *Geo);
	
};

#endif

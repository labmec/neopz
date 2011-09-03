/**
 * @file
 * @brief Contains ConvTest class which implements methods to evaluate jacobians by obtained convergence order to geometric element.
 */
#ifndef CONVTEST_H
#define CONVTEST_H


#include "pzvec.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"

/**
 * @brief Implements methods to evaluate jacobians by obtained convergence order to geometric element. \ref geometry "Geometry"
 * @author Paulo Cesar de Alvarenga Lucci
 * @since 2007
 */
class ConvTest {

public:

	/** @brief Default constructor */
    ConvTest();
	/** @brief Default destructor */
   ~ConvTest();

	/** @brief Evaluates the Jacobian by obtained Convergence Order */
    void JacobianConv(TPZGeoEl &Object, TPZVec< REAL > QsiEta);
	/** @brief Evaluates the Jacobian by obtained Convergence Order to computational element and its side */
    void JacobianConv(TPZGeoElSide &Object, TPZVec< REAL > QsiEta);

};

#endif

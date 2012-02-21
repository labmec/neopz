/**
 * @file
 * @brief Contains the TPZIntPoints class which defines integration rules.
 */

#ifndef TPZINTPOINTS_H
#define TPZINTPOINTS_H

#include <string>

#include "pzreal.h"
#include "pzvec.h"

/**
 * @ingroup integral
 * @brief Abstract class defining integration rules. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntPoints {
public:
	/** @brief Default destructor */
    virtual ~TPZIntPoints() {
    }
	
	/** @brief Returns the name of the cubature rule */
	virtual void Name(std::string &name) {
	}
	
	/** @brief Returns number of points for the cubature rule related. */
    virtual int NPoints() const = 0;
	/**
	 * @brief Returns i-th point at master element and related weight.
	 * @param i Index of the integration point at cubature rule (the sequence is not important)
	 * @param pos Vector (3d) to get the coordinates of the point at master element.
	 * @param w It gets the weight related with integration point.
	 */
	virtual void Point(int i, TPZVec<REAL> &pos, REAL &w) const = 0;
	/**
	 * @brief Sets the order of the cubature rule. 
	 * @param ord Vector of orders for each dimension, to reach a exactly integration for polinomial corresponding
	 * @param type Type of the integration rule, mainly important for 1d rules.
	 */
	virtual void SetOrder(TPZVec<int> &ord,int type = 0) = 0;
	/**
	 * @brief Gets the order of the integration rule for each dimension of the master element.
	 * @param ord Vector (3d) to get the orders of the polinomial integrated exactly.
	 */
	virtual void GetOrder(TPZVec<int> &ord) const = 0;
	/**
	 * @brief Returns the minimum order to integrate polinomials exactly for all implemented cubature rules
	 * @note To get real maxime order for especific rule call GetRealMaxOrder(), but for Gauss rule 1D the implementation
	 * can to compute a cubature rule up to on hundred order.
	 */
    virtual int GetMaxOrder() const;
	/** @brief Returns the dimension of the master element related for the cubature rule */
    virtual int Dimension() const = 0;
	/** @brief Make a clone of the related cubature rule */
	virtual TPZIntPoints *Clone() const = 0;
    virtual TPZIntPoints *PrismExtend(int order) = 0;
	
	/**
	 * @brief Sets the type of gaussian quadrature as Lobatto, Raud or Legendre rule.
	 * @note It is util for one-dimensional, quadrangular and hexahedra master element.
	 */
	virtual void SetType(int type,int order) {
	}

	/**
	 * @brief Prints information of the cubature rule.
	 * @param out Ostream to print.
	 */
    virtual void Print(std::ostream &out);
};

#endif

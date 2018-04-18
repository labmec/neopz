/**
 * @file
 * @brief Contains the TPZIntRuleList class which creates instances of all integration rules for rapid selection.
 */

#ifndef TPZINTRULELIST_H
#define TPZINTRULELIST_H

#include "pzvec.h"

class TPZGaussRule;
class TPZGaussLegendreRule;
class TPZGaussLobattoRule;
class TPZIntRuleT;
class TPZIntRuleT3D;
class TPZIntRuleP3D;

/**
 * @ingroup integral
 * @brief Creates instances of all integration rules for rapid selection. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRuleList {
	   
    
	/** @brief Pointer to an array of integration rules (Gauss Legendre) for line, quad and cube elements */
    TPZVec<TPZGaussRule* >	fintlist;

	/** @brief Pointer to an array of integration rules for triangle */
    TPZVec<TPZIntRuleT* >   fintlistT;
	/** @brief Pointer to an array of integration rules for tetrahedra */
    TPZVec<TPZIntRuleT3D* > fintlistT3D;
	/** @brief Pointer to an array of integration rules for pyramid */
    TPZVec<TPZIntRuleP3D* > fintlistP3D;
	
    public :
	
    /**
     * @ingroup integral
     * @brief Static variable with list of all integration rules
     */
    static TPZIntRuleList gIntRuleList;
    
	/** 
	 * @brief Method which initializes all integration rule vectors 
	 * @note Should be called only once!
	 */
	TPZIntRuleList();

	/** @brief Destructor of all integration rule vectors */
	~TPZIntRuleList();
	
	/** 
	 * @brief Returns a pointer to an gaussian integration rule with numint points. This method computes the number of points to right Gaussian quadrature 
	 * @param order Degree of the polinomial for which the integration is exact. 
	 * @param type Type values: 0 - Gauss Legendre (default), 1 - Gauss Lobatto, 2 - Gauss Jacobi
	 */
	TPZGaussRule *GetRule(int order,int type = 0);
	
	/** 
	 * @brief Returns a pointer to an integration rule for a triangle 
	 * @param order Degree of the polinomial for which the integration is exact. 
	 */
	TPZIntRuleT *GetRuleT(int order);
	/** 
	 * @brief Returns a pointer to an integration rule for a tetrahedra 
	 * @param order Degree of the polinomial for which the integration is exact. 
	 */
	TPZIntRuleT3D *GetRuleT3D(int order);
	/** 
	 * @brief Returns a pointer to an integration rule for a  pyramid
	 * @param order Degree of the polinomial for which the integration is exact. 
	 */
	TPZIntRuleP3D *GetRuleP3D(int order);
};



#endif

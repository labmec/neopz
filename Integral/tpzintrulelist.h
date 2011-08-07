/**
 * @file
 * @brief Contains the TPZIntRuleList class which creates instances of all integration rules for rapid selection.
 */
//
// C++ Interface: tpzintrulelist
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef TPZINTRULELIST_H
#define TPZINTRULELIST_H

class TPZIntRule;
class TPZIntRuleT;
class TPZIntRuleT3D;
class TPZIntRuleP3D;

/**
 * @ingroup integral
 * @brief Creates instances of all integration rules for rapid selection. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRuleList {
	
	/** @brief Number of integration rules available */
    int		intavail;
	/** @brief Number of integration rules available for triangles */
    int        	intavailT;
    int        	intavailT3D;
    int        	intavailP3D;
	/** @brief Pointer to an array of integration rules */
    TPZIntRule	**intlist;
	/** @brief Pointer to an array of integration rules */
    TPZIntRuleT   **intlistT;
    TPZIntRuleT3D **intlistT3D;
    TPZIntRuleP3D **intlistP3D;
	
    public :
	
	/** @brief Method which initializes all integration rules 
	 * @note should be called only once!
	 */
	TPZIntRuleList();
	
	~TPZIntRuleList();
	
	/** @brief Returns a pointer to an integration
	 * rule with numint integration points
	 */
	TPZIntRule *GetRule(int numint);
	
	/** @brief Returns a pointer to an integration
	 * rule for a triangle
	 */
	TPZIntRuleT *GetRuleT(int numint);
	TPZIntRuleT3D *GetRuleT3D(int numint);
	TPZIntRuleP3D *GetRuleP3D(int numint);
};

/** @ingroup integral
 * @brief Extern variable with list of all integration rules
 */
extern  TPZIntRuleList  gIntRuleList;

#endif

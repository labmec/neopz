/**
 * @file
 * @brief Contains the TPZIntRuleP3D class which defines the integration rule for pyramid.
 */
//
// C++ Interface: tpzintrulep3d
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef TPZINTRULEP3D_H
#define TPZINTRULEP3D_H

#include "pzreal.h"
template<class T>
class TPZVec;

/**
 * @ingroup integral
 * @brief Integration rule for pyramid. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRuleP3D {
	
	friend class TPZIntRuleList;
	
	/** @brief Number of integration points for this object */
	short	 fNumInt;
	/** @brief Location of the integration point Ksi */
	REAL	*fLocationKsi;
	/** @brief Location of the integration point Eta */
	REAL	*fLocationEta;
	/** @brief Location of the integration point ZEta */
	REAL	*fLocationZeta;
	/** @brief Weight of the integration point */
	REAL	*fWeight;
	
	TPZIntRuleP3D(int i = 2);
	~TPZIntRuleP3D();
	
public:
	
    enum {NUMINT_RULESP3D = 8};
	
	/** @brief Returns number of integration points */
    int NInt() const{ return fNumInt; }
	
	/** @brief Returns location of the ith point */
	void Loc(int i, TPZVec<REAL> &pos) const;
	
	/** @brief Returns weight for the ith point */
	REAL W(int i) const;
	
};

#endif

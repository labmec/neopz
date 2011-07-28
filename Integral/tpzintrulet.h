//
// C++ Interface: tpzintrulet
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZINTRULET_H
#define TPZINTRULET_H

#include "pzreal.h"
template<class T>
class TPZVec;

/**
 * @ingroup integral
 * @brief Integration rule (points and weights) for triangles
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRuleT {
	friend class TPZIntRuleList;
	
	/// number of integration points for this object
	short	   fNumInt;
	/// location of the integration point Ksi
	REAL	*fLocationKsi;
	/// location of the integration point Eta
	REAL	*fLocationEta;
	/// weight of the integration point
	REAL	*fWeight;
	
	TPZIntRuleT(int i);
	~TPZIntRuleT();
	
public:
    enum {NUMINT_RULEST = 19};
	
	///return number of integration points
    short NInt() const { return fNumInt;}
	
	/// return location of the ith point
    void Loc(int i, TPZVec<REAL> &pos) const;
	
	///return weight for the ith point
    REAL W(int i) const;
	
};

#endif

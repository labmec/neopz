//
// C++ Interface: tpzintrulet3d
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZINTRULET3D_H
#define TPZINTRULET3D_H

#include "pzreal.h"
template<class T>
class TPZVec;

/**
 * @ingroup integral
 * @brief Integration rule for tetrahedra
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRuleT3D {
	
    friend class TPZIntRuleList;
    
	/// number of integration points for this object
    short	 fNumInt;
	/// location of the integration point Ksi
    REAL	*fLocationKsi;
	/// location of the integration point Eta
    REAL	*fLocationEta;
	/// location of the integration point ZEta
    REAL	*fLocationZeta;
	/// weight of the integration point
    REAL	*fWeight;
	
    TPZIntRuleT3D(int i = 2);
    ~TPZIntRuleT3D();
	
public:
	enum {NUMINT_RULEST3D = 8};
	
 	///return number of integration points
	int NInt() const { return fNumInt;}
	
	/// return location of the ith point
	void Loc(int i, TPZVec<REAL> &pos) const;
	
	///return weight for the ith point
	REAL W(int i) const;
	
};

#endif

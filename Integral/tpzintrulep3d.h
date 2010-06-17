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
//
#ifndef TPZINTRULEP3D_H
#define TPZINTRULEP3D_H

#include "pzreal.h"
template<class T>
class TPZVec;

/**
Integration rule for pyramid

	@author Philippe R. B. Devloo <phil@fec.unicamp.br>
*/
class TPZIntRuleP3D{

  friend class TPZIntRuleList;
  
  
  short	 fNumInt;		// number of integration points for this object
  REAL	*fLocationKsi;	// location of the integration point Ksi
  REAL	*fLocationEta;	// location of the integration point Eta
  REAL	*fLocationZeta;	// location of the integration point Eta
  REAL	*fWeight;		// weight of the integration point

  TPZIntRuleP3D(int i = 2);
  ~TPZIntRuleP3D();

  public:

    enum {NUMINT_RULESP3D = 8};

    int NInt() const{ return fNumInt;}	//return number of integration points

    void Loc(int i, TPZVec<REAL> &pos) const;			   //return location of the ith pot

    REAL W(int i) const;						//return weight for the ith point

};

#endif

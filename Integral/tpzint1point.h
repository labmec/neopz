//
// C++ Interface: tpzint1point
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef TPZINT1POINT_H
#define TPZINT1POINT_H

#include "tpzintpoints.h"
#include "tpzprinteg.h"

/**
 * @ingroup integral
 * @brief Integration rule for one point
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZInt1Point : public TPZIntPoints {
	
public:
	
    enum {Dim = 0};
    TPZInt1Point(int order=0);
    TPZInt1Point(TPZVec<int> &ord);
	TPZInt1Point(const TPZInt1Point &copy ) : TPZIntPoints(copy)
	{
	}
    virtual ~TPZInt1Point();
    void SetOrder(TPZVec<int> &ord);
    int NPoints() const;
    void Point(int ip, TPZVec<REAL> &pos, REAL &w) const;
    void GetOrder(TPZVec<int> &ord) const;
    int GetMaxOrder() const;  
    int Dimension() const
    {
		return Dim;
    }
    TPZIntPoints *PrismExtend(int order);
	TPZIntPoints *Clone() const
	{
		return new TPZInt1Point(*this);
	}
	
};

inline TPZInt1Point::~TPZInt1Point() 
{
}

inline TPZInt1Point::TPZInt1Point(int order) {
}

inline void  TPZInt1Point::SetOrder(TPZVec<int> &ord) {
}

inline int TPZInt1Point::NPoints() const{
	return 1;
}

inline void TPZInt1Point::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
#ifndef NODEBUG
	if(ip!=0) {
		std::cout << "TPZInt1Point:: Bad number point " << ip << std::endl;
		return;
	}
#endif
	w = 1.;
}

inline void TPZInt1Point::GetOrder(TPZVec<int> &/* ord */) const{
}


inline int TPZInt1Point::GetMaxOrder() const {return 1;}

inline TPZIntPoints *TPZInt1Point::PrismExtend(int order) 
{
	return new TPZPrInteg<TPZInt1Point> (order);
}

#endif

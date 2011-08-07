/**
 * @file
 * @brief Contains the TPZPrInteg class which defines prismatic extension of an integration rule.
 */
//
// C++ Interface: tpzprinteg
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef TPZPRINTEG_H
#define TPZPRINTEG_H

#include "pzreal.h"
#include "pzvec.h"

class TPZIntPoints;

#include "tpzintrulelist.h"
#include "tpzintrule.h"

/**
 * @ingroup integral
 * @brief Prismatic extension of an integration rule. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
template< class TFather>
class TPZPrInteg : public TFather
{
public:
	enum {Dim = TFather::Dim+1};
	int fOrdKsi;
	TPZIntRule *fIntP;
public:
    TPZPrInteg(int order = 0) : TFather(order)
    {
		if(order>0)
		{
			fOrdKsi = order;
			fIntP   = gIntRuleList.GetRule(order);
		}
		else
		{
			fOrdKsi = 0;
			fIntP = 0;
		}
    }
    TPZPrInteg(TPZVec<int> &order) : TFather()
    {
		SetOrder(order);
    }
	
	TPZPrInteg(const TPZPrInteg &copy ) : TFather(copy), fOrdKsi(copy.fOrdKsi), fIntP(copy.fIntP)
	{
	}
    virtual ~TPZPrInteg();
    
    int NPoints() const
    {
		return TFather::NPoints()*fIntP->NInt();
    }
    void Point(int ip, TPZVec<REAL> &pos, REAL &w) const
    {
		int ipf = ip/fIntP->NInt();
		int iploc = ip%(fIntP->NInt());
		TFather::Point(ipf,pos,w);
		pos[Dim-1] = fIntP->Loc(iploc);
		w *= fIntP->W(iploc);
    }
    void SetOrder(TPZVec<int> &ord)
    {
#ifndef NODEBUG
		if(ord.NElements() < Dim) {
			std::cout << "TPZPrInteg::SetOrder: number of integration points specified smaller than dimension\n";
			return;
		}
#endif
		TFather::SetOrder(ord);
		fOrdKsi = ord[Dim-1];
		fIntP   = gIntRuleList.GetRule(ord[Dim-1]);
    }
    void GetOrder(TPZVec<int> &ord) const
    {
#ifndef NODEBUG
		if(ord.NElements() < Dim) {
			std::cout << "TPZPrInteg::GetOrder: number of integration points specified smaller than dimension\n";
			return;
		}
#endif
		TFather::GetOrder(ord);
		ord[Dim-1] = fOrdKsi;
    }
    
    int Dimension() const
    {
		return Dim;
    }
    
    virtual TPZIntPoints *PrismExtend(int order);
	
    int GetMaxOrder() const
    {
		int fatmax = TFather::GetMaxOrder();
		return (fatmax > fOrdKsi) ? fatmax : fOrdKsi;
		
    }
	virtual TPZIntPoints *Clone() const
	{
		return new TPZPrInteg<TFather>(*this);
	}
	
};

#endif

//
// C++ Interface: tpzintpoints
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZINTPOINTS_H
#define TPZINTPOINTS_H

#include "pzreal.h"
template<class T>
class TPZVec;
#include "pzvec.h"
             
/**
Abstract class defining integration rules

	@author Philippe R. B. Devloo <phil@fec.unicamp.br>
*/
class TPZIntPoints{
public:
    virtual ~TPZIntPoints()
    {
    }
    virtual int NPoints() = 0;
    virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) = 0;
    virtual void SetOrder(TPZVec<int> &ord) = 0;
    virtual void GetOrder(TPZVec<int> &ord) = 0;
    virtual int GetMaxOrder() = 0;
    virtual int Dimension() = 0;
    virtual TPZIntPoints *PrismExtend(int order) = 0;
    virtual void Print(std::ostream &out)
    {
      int np = NPoints();
      out << "integration rule numpoints " << np << std::endl;
      int ip;
      TPZVec<REAL> pos(Dimension());
      REAL w;
      for(ip=0; ip<np; ip++)
      {
        Point(ip,pos,w);
        out << "ip " << ip << " pos " << pos << " w " << w << std::endl;
      }
    }
};

#endif

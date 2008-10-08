//
// C++ Implementation: tpzpoint
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzpoint.h"
#include "pzquad.h"
#include "tpzint1point.h"
#include "pzeltype.h"

namespace pztopology {

TPZPoint::TPZPoint()
{
}


TPZPoint::~TPZPoint()
{
}

TPZIntPoints *TPZPoint::CreateSideIntegrationRule(int side, int order)
{
  return new IntruleType(order);
}


MElementType TPZPoint::Type()
{
  return EPoint;
}

MElementType TPZPoint::Type(int side)
{
  switch(side) {
    case 0:
      return EPoint;
    default:
      return ENoType;
  }
}

}

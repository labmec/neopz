//
// C++ Implementation: pzmaterialid
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@corona>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pzmaterialid.h"

#include "pzmat2dlin.h"

void RegisterMaterialClasses() {

TPZSaveable::Register(TPZMAT2DLINID,Restore<TPZMat2dLin>);

}
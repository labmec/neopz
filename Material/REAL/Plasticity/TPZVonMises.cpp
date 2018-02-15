/*
 *  TPZVonMises.cpp
 *  ElastoPlasticModels
 *
 *  Created by Diogo Cecilio on 12/17/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZVonMises.h"

int TPZVonMises::ClassId() const{
    return Hash("TPZVonMises") ^ VONMISESPARENT::ClassId() << 1;
}
/*
 *  TPZWillamWarnke.cpp
 *  ElastoPlasticModels
 *
 *  Created by Diogo Cecilio on 12/14/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZWillamWarnke.h"

int TPZWillamWarnke::ClassId() const{
    return Hash("TPZWillamWarnke") ^ WILLAMWARNKEPARENT::ClassId() << 1;
}
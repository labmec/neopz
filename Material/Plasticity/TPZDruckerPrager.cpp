/*
 *  TPZDruckerPrager.cpp
 *  ElastoPlasticModels
 *
 *  Created by Diogo Cecilio on 6/11/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZDruckerPrager.h"

int TPZDruckerPrager::ClassId() const{
    return Hash("TPZDruckerPrager") ^ DRUCKERPARENT::ClassId() << 1;
}
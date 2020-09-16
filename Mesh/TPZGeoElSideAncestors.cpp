//
// Created by Gustavo A. Batistela on 15/09/2020.
//

#include "TPZGeoElSideAncestors.h"

void TPZGeoElSideAncestors::BuildAncestors() {

    int nlowerlevels = 0;
    if (!fCurrent) {
        fAncestors.resize(nlowerlevels);
        return;
    }

    TPZGeoElSide lower = fCurrent.StrictFather();

    TPZGeoElSide neighbour = fCurrent.Neighbour();
    TPZGeoElSide neighLower;
    while (neighbour != fCurrent) {
        neighLower = neighbour.StrictFather();
        if (neighLower) break;
        neighbour = neighbour.Neighbour();
    }

    while(lower || neighLower) {
        nlowerlevels++;
        fAncestors.resize(nlowerlevels);
        fAncestors[nlowerlevels - 1] = std::make_pair(TPZGeoElSideAncestors(lower), TPZGeoElSideAncestors(neighLower));

        lower = lower.StrictFather();
        neighLower = neighLower.StrictFather();
    }
}

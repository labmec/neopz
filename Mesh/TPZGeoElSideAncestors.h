//
// Created by Gustavo A. Batistela on 15/09/2020.
//

#ifndef TPZGEOELSIDEANCESTORS_H
#define TPZGEOELSIDEANCESTORS_H

#include "pzgeoelside.h"
#include "pzgeoel.h"

class TPZGeoElSideAncestors {
    // The TPZGeoElSide from which we derive the ancestors tree
    TPZGeoElSide fCurrent;

    // Data structure with a pair o TPZGeoElSideAncestors for the lower side and its neighbour
    TPZVec<std::pair<TPZGeoElSideAncestors, TPZGeoElSideAncestors>> fAncestors;

    // Compute the ancestor data structure
    void BuildAncestors();

public:

    TPZGeoElSideAncestors() : fCurrent(), fAncestors() {}

    TPZGeoElSideAncestors(const TPZGeoElSideAncestors &cp) = default;

    TPZGeoElSideAncestors &operator=(const TPZGeoElSideAncestors &cp) = default;

    TPZGeoElSideAncestors(const TPZGeoElSide &current) : fCurrent(current) {
        BuildAncestors();
    }

    void SetCurrent(const TPZGeoElSide &current) {
        fCurrent = current;
        BuildAncestors();
    }
};

#endif

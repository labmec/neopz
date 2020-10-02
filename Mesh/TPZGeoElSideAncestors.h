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

    // Data structure containing the ancestor tree. It is filled by AddAncestor method.
    TPZStack<std::pair<TPZGeoElSide, TPZGeoElSide>> fAncestors;

    // Compute the ancestor data structure
    void BuildAncestors();

    // Fills fAncestor data structure as follows:
    // Insert a pair containing the given side and its father side, if it exists.
    // For every neighbour, insert a pair of its side and its father side, if it exists, and call AddAncestor for the
    // father.
    void AddAncestor(TPZGeoElSide gelside);

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
    
    TPZTransform<REAL> BuildTransform(TPZGeoElSide larger);
    
    /// return true is a (strict) larger element with matid exists
    TPZGeoElSide HasLarger(int matid);
    
    /// return the element/side of the larger element
    TPZGeoElSide LargeSide(TPZGeoEl *large);
};

#endif

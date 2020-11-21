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
    TPZStack<std::pair<TPZGeoElSide, TPZGeoElSide>> fAncestors;

    // Compute the ancestor data structure
    void BuildAncestors();
    
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
    
    /// return true is a (strict) larger element with matid exists
    TPZGeoElSide HasLarger(const std::set<int> &matid);
    
    /// return true is a (strict) larger element with matid exists
    TPZGeoElSide HasLargerorEqual(int matid);

    /// return true is a (strict) larger element with matid exists
    TPZGeoElSide HasLargerorEqual(const std::set<int> &matid);

/// return the element/side of the larger element
    TPZGeoElSide LargeSide(TPZGeoEl *large);
};

#endif

//
// Created by Gustavo A. Batistela on 15/09/2020.
//

#ifndef TPZGEOELSIDEPARTITION_H
#define TPZGEOELSIDEPARTITION_H

#include "pzgeoelside.h"
#include "pzgeoel.h"

class TPZGeoElSidePartition {
    // the TPZGeoElSide that defines the partition
    TPZGeoElSide fCurrent;
    // the same dimension objects whose closure form the fCurrent set
    TPZVec<TPZGeoElSidePartition> fPartition;

    // compute the partition data structure
    void BuildPartition();

public:

    TPZGeoElSidePartition() : fCurrent(), fPartition() {}

    TPZGeoElSidePartition(const TPZGeoElSidePartition &cp) = default;

    TPZGeoElSidePartition &operator=(const TPZGeoElSidePartition &cp) = default;

    TPZGeoElSidePartition(const TPZGeoElSide &current) : fCurrent(current) {
        BuildPartition();
    }

    void SetCurrent(const TPZGeoElSide &current) {
        fCurrent = current;
        BuildPartition();
    }

    /// checks whether an element with MaterialID matid is neighbour of a partition of fCurrent
    TPZGeoElSide HasHigherLevelNeighbour(int matid) const;

    // fills a stack of neighbours which matid is neighbour of a partition of fCurrent
    // returns true if a neighbour of the given material id is found for every subelement of the partition
    bool HigherLevelNeighbours(TPZStack<TPZGeoElSide> &neighbours, int matid) const;
};

#endif

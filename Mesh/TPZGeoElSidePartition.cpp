//
// Created by Gustavo A. Batistela on 15/09/2020.
//

#include "TPZGeoElSidePartition.h"

// compute the partition data structure
void TPZGeoElSidePartition::BuildPartition() {
    if (!fCurrent) {
        fPartition.resize(0);
        return;
    }
    TPZStack<TPZGeoElSide> subels, subelfit;
    if (fCurrent.HasSubElement() && fCurrent.NSubElements() > 1) {
        // this method returns all subelements
        fCurrent.GetSubElements2(subels);
    } else {
        TPZGeoElSide neighbour = fCurrent.Neighbour();
        while (neighbour != fCurrent) {
            if (neighbour.HasSubElement() && neighbour.NSubElements() > 1) {
                neighbour.GetSubElements2(subels);
                break;
            }
            neighbour = neighbour.Neighbour();
        }
    }
    // filter out the subels of the same dimension
    int dim = fCurrent.Dimension();
    int nel = subels.size();
    for (int el = 0; el < nel; el++) {
        if (subels[el].Dimension() == dim) {
            subelfit.Push(subels[el]);
        }
    }
    fPartition.resize(subelfit.size());
    for (int el = 0; el < subelfit.size(); el++) {
        TPZGeoElSidePartition sidepart(subelfit[el]);
        fPartition[el] = TPZGeoElSidePartition(sidepart);
    }
}

// checks whether an element with MaterialID matid is neighbour of a partition of fCurrent
TPZGeoElSide TPZGeoElSidePartition::HasHigherLevelNeighbour(int matid) const {
    int nel = fPartition.size();
    for (int el = 0; el < nel; el++) {
        // As hasneighbour does not verify the element itself, we do it here
        if (fPartition[el].fCurrent.Element()->MaterialId() == matid) {
            return fPartition[el].fCurrent;
        }
        TPZGeoElSide nextlev = fPartition[el].fCurrent.HasNeighbour(matid);
        if (nextlev) {
            return nextlev;
        }
    }
    for (int el = 0; el < nel; el++) {
        TPZGeoElSide higher = fPartition[el].HasHigherLevelNeighbour(matid);
        if (higher) {
            return higher;
        }
    }
    return TPZGeoElSide();
}

bool TPZGeoElSidePartition::HigherLevelNeighbours(TPZStack<TPZGeoElSide> &neighbours, int matid) const {
    bool foundNeighbourForEverySubElement = true;
    for (auto &partition : fPartition) {
        // First check the partition side itself
        TPZGeoElSide neighbour = partition.fCurrent.HasNeighbour(matid);
        if (neighbour) {
            neighbours.Push(neighbour);
        } else if (partition.fPartition.size() == 0) {
            return false;
        } else {
            foundNeighbourForEverySubElement = partition.HigherLevelNeighbours(neighbours, matid);
        }
    }
    return foundNeighbourForEverySubElement;
}
